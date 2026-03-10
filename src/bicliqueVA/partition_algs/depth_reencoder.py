
from ..utils.graph import Graph, erdos_renyi
from ..utils.formula_to_graph import form2graph
from ..utils.bicliques_to_formula import biclique_partition_to_formula


from typing import Optional
from eznf import modeler
import math
import itertools
import time


def pre_formula(input_formula_filename):
    enc = modeler.Modeler()
    

    with open(input_formula_filename, 'r') as f:
        for line in f:
            if line.startswith('p cnf'):
                tokens = line.split(' ')
                n_vars = int(tokens[2])
                for v in range(1, n_vars + 1):
                    enc.add_var(v)
            elif line.startswith('c') or line.strip() == '':
                continue
            else:
                tokens = line.split(' ')[:-1]
                clause = [int(lit) for lit in tokens]
                assert 0 not in clause, f"Clause should not contain 0; {clause}"
                if len(clause) != 2:
                    enc.add_clause(clause)
                elif clause[0] > 0 or clause[1] > 0:
                    enc.add_clause(clause)
    return enc


def _estimate_default_r(n):
    best_r = 3
    best_est = 1e18
    for r in range(3, min(12, n) + 1):
        est_s = n * (1 << (r - 1))
        est_a = n * math.ceil(n / (2 * r))
        est = est_s + est_a
        if est < best_est:
            best_est = est
            best_r = r
    return best_r


def _partition_vertices(vertices, block_size):
    n = len(vertices)
    num_blocks = math.ceil(n / block_size)
    blocks = []
    vertex_to_block = {}
    for i in range(1, num_blocks + 1):
        lo = (i - 1) * block_size
        hi = min(i * block_size, n)
        block = vertices[lo:hi]
        blocks.append(block)
        for v in block:
            vertex_to_block[v] = i - 1
    return blocks, vertex_to_block


def _collect_depth_groups(graph, r):
    vertices = graph.vertices
    n = graph.n_vertices()
    if n == 0:
        return 0, []

    parts, vertex_to_part = _partition_vertices(vertices, r)
    num_parts = len(parts)
    is_edge = graph.is_edge

    def R(i, j):
        dij = min(abs(i - j), num_parts - abs(i - j))
        return (i < j and dij % 2 == 0) or (i > j and dij % 2 == 1)

    within_part_clauses = 0
    for part in parts:
        for (u, v) in itertools.combinations(part, 2):
            if is_edge(u, v):
                within_part_clauses += 1

    groups_data = []

    for i, part_i in enumerate(parts):
        groups = {}
        for v in vertices:
            j = vertex_to_part[v]
            if R(j, i):
                mask = 0
                for idx, u in enumerate(part_i):
                    if is_edge(u, v):
                        mask |= (1 << idx)
                if mask != 0:
                    if mask not in groups:
                        groups[mask] = set()
                    groups[mask].add(v)

        for mask, right_set in groups.items():
            s_size = mask.bit_count()
            if s_size > 0 and len(right_set) > 0:
                groups_data.append((s_size, tuple(sorted(right_set))))

    return within_part_clauses, groups_data


def _estimate_depth_clause_count_from_groups(vertices, within_part_clauses, groups_data, p):
    _, vertex_to_block = _partition_vertices(vertices, p)

    s_weight = 0
    a_weight = 0
    direct_clauses = 0
    used_z_pairs = set()

    for s_size, right_tuple in groups_data:
        a_size = len(right_tuple)
        if s_size + a_size >= s_size * a_size:
            direct_clauses += s_size * a_size
            continue

        s_weight += s_size
        right_by_block = {}
        for v in right_tuple:
            b_v = vertex_to_block[v]
            if b_v not in right_by_block:
                right_by_block[b_v] = []
            right_by_block[b_v].append(v)

        for u_list in right_by_block.values():
            k = len(u_list)
            a_weight += (k // 2) + (k % 2)
            for t in range(k // 2):
                used_z_pairs.add((u_list[2 * t], u_list[2 * t + 1]))

    b_clauses = 2 * len(used_z_pairs)
    return within_part_clauses + s_weight + a_weight + direct_clauses + b_clauses


def _adaptive_depth_parameters(graph):
    n = graph.n_vertices()
    if n <= 4:
        return max(2, n), max(2, n)

    lg = lambda x: math.log2(x) if x > 0 else 0.0
    p_base = max(4, 2 * int(n // (lg(n) * lg(n))) + 10)
    r_base = _estimate_default_r(n)

    r_candidates = sorted({
        max(2, r_base - 1),
        max(2, r_base),
        min(n, r_base + 1),
    })

    p_candidates = sorted({
        max(4, p_base // 4),
        max(4, p_base // 3),
        max(4, p_base // 2),
        max(4, p_base),
        max(4, int(round(math.sqrt(n)))),
        8,
        12,
    })
    p_candidates = [p for p in p_candidates if p <= n]
    if len(p_candidates) == 0:
        p_candidates = [max(2, n)]

    best_cost = None
    best_r = r_base
    best_p = p_base
    vertices = graph.vertices
    for r in r_candidates:
        within_part_clauses, groups_data = _collect_depth_groups(graph, r=r)
        for p in p_candidates:
            est_cost = _estimate_depth_clause_count_from_groups(
                vertices=vertices,
                within_part_clauses=within_part_clauses,
                groups_data=groups_data,
                p=p,
            )
            if best_cost is None or est_cost < best_cost:
                best_cost = est_cost
                best_r = r
                best_p = p

    return best_r, best_p


def depth_reencode(input_formula_filename, forced_r: Optional[int] = None, forced_p: Optional[int] = None):
    g = form2graph(input_formula_filename)
    enc = pre_formula(input_formula_filename)
    n = g.n_vertices()
    print(f"Number of clauses in original formula = {enc.n_clauses()}")

    adaptive_r, adaptive_p = _adaptive_depth_parameters(g)
    r = forced_r if forced_r is not None else adaptive_r
    p = forced_p if forced_p is not None else adaptive_p
    print(f"Adaptive depth parameters: r={adaptive_r}, p={adaptive_p}")
    if forced_r is not None or forced_p is not None:
        print(f"Using overridden depth parameters: r={r}, p={p}")

    vc = enc.n_vars()

    print(f"p = {p}")
    num_blocks = math.ceil(n / p)
    blocks = []
    vertex_to_block = {}
    for i in range(1, num_blocks + 1):
        lo = (i - 1) * p
        hi = min(i * p, n)
        # print(f"lo = {lo}, hi = {hi}")
        blocks.append(g.vertices[lo:hi])  # B_i
        for v in blocks[-1]:
            vertex_to_block[v] = i-1

    b_clauses = 0
    # Track which z_{u,v} have actually been instantiated.
    # This is equivalent to "create all, then delete unused", but cheaper.
    b_marked = {}
    print(f"b_clauses = {b_clauses}")

    print(f"depth using r = {r}")
   
    num_parts = math.ceil(n / r)
    def R(i, j):
        dij = min(abs(i-j), num_parts - abs(i-j))
        return (i < j and dij % 2 == 0) or (i > j and dij % 2 == 1)
        
   
    parts = []
    vertex_to_part = {}
    for i in range(1, num_parts + 1):
        lo = (i - 1) * r
        hi = min(i * r, n)
        parts.append(g.vertices[lo:hi])  # P_i
        for v in parts[-1]:
            vertex_to_part[v] = i-1

    
    # edges within parts
    within_part_clauses = 0
    for P_i in parts:
        for (u, v) in itertools.combinations(P_i, 2):
            if g.is_edge(u, v):
                enc.add_clause([-u, -v])
                within_part_clauses += 1
    
    print(f"Number of within-part clauses = {within_part_clauses}, enc.n_clauses() = {enc.n_clauses()}")

    vc = enc.n_vars()
    # accross parts
    S_weight = 0
    A_weight = 0
    direct_clauses = 0
    sum_right_sets = 0
    n_right_sets = 0
    for i, P_i in enumerate(parts):  
        A = {}  
        for v in g.vertices:
            j = vertex_to_part[v]
            if R(j, i):
                mask = 0
                for idx, u in enumerate(P_i):
                    if g.is_edge(u, v):
                        mask |= (1 << idx)   # A(S) ← A(S) ∪ {v}
                if mask != 0:  
                    if mask not in A:
                        A[mask] = set()
                    A[mask].add(v)

        new_vertices = set()
        used_right = set()
        new_edges = []

        for mask, right_set in A.items():
            S = tuple(u for bit, u in enumerate(P_i) if (mask >> bit) & 1)
            if len(S) > 0 and len(right_set) > 0:
                if len(S) + len(right_set) >= len(S) * len(right_set):
                    for v in S:
                        for u in right_set:
                            enc.add_clause([-v, -u])
                            direct_clauses += 1
                    continue
                sum_right_sets += len(right_set)
                n_right_sets += 1
                enc.add_var(vc + 1)
                new_vertices.add(vc + 1)
                for v in right_set:
                    used_right.add(v)

                new_edges.extend([(v, vc + 1) for v in right_set])
      
                for v in S:
                    enc.add_clause([-v, vc+1])
                S_weight += len(S)
                vc += 1

        g2 = Graph(vertices=list(new_vertices) + list(used_right))
        for (u, v) in new_edges:
            g2.add_edge(u, v)


        right_by_block = {}
        for v in used_right:
            b_v = vertex_to_block[v]
            if b_v not in right_by_block:
                right_by_block[b_v] = set()
            right_by_block[b_v].add(v)
        
        for v in new_vertices:
            adj = g2.neighbors(v)
            adj_per_block = {}
            for u in adj:
                b_u = vertex_to_block[u]
                if b_u not in adj_per_block:
                    adj_per_block[b_u] = set()
                adj_per_block[b_u].add(u)

            for b_u, u_set in adj_per_block.items():
                u_sorted = tuple(sorted(u_set))
                for k in range(len(u_sorted)//2):
                    pair = (u_sorted[2*k], u_sorted[2*k+1])
                    if pair not in b_marked:
                        z_name = f"z_{pair[0], pair[1]}"
                        enc.add_var(z_name)
                        enc.add_clause([-pair[0], z_name])
                        enc.add_clause([-pair[1], z_name])
                        b_clauses += 2
                        b_marked[pair] = True
                    enc.add_clause([-v, f"-z_{pair[0], pair[1]}"])
                    # enc.add_clause([-v, f"-z_{u_sorted[2*k], u_sorted[2*k+1]}"])
                    # b_clauses += 2
                    A_weight += 1
                if len(u_sorted) % 2 == 1:
                    enc.add_clause([-v, -u_sorted[-1]])
                    A_weight += 1

        
        

        # originally = len(right_set)
        # reencoded = 0
        # # for right 
        # rB = tuple(sorted(right_set))
        # for k in range(len(rB)//2):
        #     if (rB[2*k], rB[2*k+1]) not in b_marked:
        #         enc.add_var(f"z_{rB[2*k], rB[2*k+1]}")
        #         enc.add_clause([-rB[2*k], f"z_{rB[2*k], rB[2*k+1]}"])
        #         enc.add_clause([-rB[2*k+1], f"z_{rB[2*k], rB[2*k+1]}"])
        #         b_clauses += 2
        #         b_marked[(rB[2*k], rB[2*k+1])] = True
        #     enc.add_clause([-(vc+1), f"-z_{rB[2*k], rB[2*k+1]}"])
        #     A_weight += 1
        #     reencoded += 1
            
        # if len(rB) % 2 == 1:
        #     enc.add_clause([-(vc+1), -rB[-1]])
        #     A_weight += 1
        #     reencoded += 1
                
                # print(f"Reencoded {originally} clauses into {reencoded}")
              
        # print(f"used right = {len(used_right)}")
    avg_right_size = (sum_right_sets / n_right_sets) if n_right_sets > 0 else 0.0
    print(f"Average right set size = {avg_right_size}")
    print(f"Within part clauses = {within_part_clauses}")
    print(f"Number of b_clauses = {b_clauses}")
    print(f"Number of direct clauses = {direct_clauses}")
    print(f"Total weight of across-part bicliques: |A|={A_weight}, |S|={S_weight}, total={A_weight + S_weight}")
    print(f"Total clauses = {within_part_clauses + b_clauses + direct_clauses + A_weight + S_weight}?")
    print(f"total clauses in enc = {enc.n_clauses()}")
    print(f"Theoretical prediction for A-weight {n* math.ceil(n / (2*r))}, S-weight {n* (2**r) // 2}")
    # return bicliques
    return enc


def run_depth_reencode(
    input_formula_filename,
    timeout=None,
    output_formula_filename='tmp/depth_reencoded.cnf',
    forced_r: Optional[int] = None,
    forced_p: Optional[int] = None,
):
    start = time.perf_counter_ns()
    encoding = depth_reencode(
        input_formula_filename=input_formula_filename,
        forced_r=forced_r,
        forced_p=forced_p,
    )
    end = time.perf_counter_ns()
    encoding.serialize(output_formula_filename)
    vars = encoding.n_vars()
    cls = encoding.n_clauses()
    runtime_ms = (end - start) / 1e6
    return vars, cls, runtime_ms


if __name__ == "__main__":
    # g = erdos_renyi(1000, 0.5)
    # print("Generated graph,  # edges =", g.n_edges)
    # B = biclique_partition_components(g, delta_r=0)
    formula = "tmp/random_2cnf_1000_0.5.cnf"
    vars, cls, runtime = run_biclique_partition_on_formula(formula)
    print(f"Reencoded formula with {vars} variables and {cls} clauses.")
    print(f"Runtime: {runtime} ms")
    # edges = g.n_edges
    # print(f"Graph with {len(g.vertices)} vertices and {edges} edges")
    # total_weight = 0
    # for i, (L, R) in enumerate(B):
    #     total_weight += min(len(L)*len(R), len(L)+len(R))
    # print(f"Found {len(B)} bicliques with total weight {total_weight}")
