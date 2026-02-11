
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

def depth_reencode(input_formula_filename):
    g = form2graph(input_formula_filename)
    enc = pre_formula(input_formula_filename)
    n = g.n_vertices()
    print(f"Number of clauses in original formula = {enc.n_clauses()}")

    min_est = 1e9
    min_pr = -1
    for pr in range(3, 12):
        est_S = n * (2**pr // 2)
        est_A = n * math.ceil(n / (2*pr))
        print(f"pr = {pr}, est_S = {est_S}, est_A = {est_A}, est_total = {est_S + est_A}")
        if est_S + est_A < min_est:
            min_est = est_S + est_A
            min_pr = pr
    
    vc = enc.n_vars()


   

    lg = lambda x: math.log2(x) if x > 0 else 0.0
    p = 2*int(n // (lg(n)*lg(n))) + 10
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
    b_marked = {}

    for block in blocks:
        for (u, v) in itertools.combinations(block, 2):
            enc.add_var(f"z_{u, v}")
            b_marked[(u, v)] = False
            enc.add_clause([-u, f"z_{u, v}"])
            enc.add_clause([-v, f"z_{u, v}"])
            b_clauses += 2  
    print(f"b_clauses = {b_clauses}")

    # r = max(1, int(math.floor(lg(n) - 1.1 * lg(lg(n))))) -1
    # r = min_pr 
    r = 5
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
                # if len(S) + len(right_set) >= len(S) * len(right_set):
                #     for v in S:
                #         for u in right_set:
                #             enc.add_clause([-v, -u])
                #             direct_clauses+=1
                #     continue
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
                    # if (u_sorted[2*k], u_sorted[2*k+1]) not in b_marked:
                        # enc.add_var(f"z_{u_sorted[2*k], u_sorted[2*k+1]}")
                    enc.add_clause([-v, f"-z_{u_sorted[2*k], u_sorted[2*k+1]}"])
                    # enc.add_clause([-v, f"-z_{u_sorted[2*k], u_sorted[2*k+1]}"])
                    # b_clauses += 2
                    A_weight += 1
                    b_marked[(u_sorted[2*k], u_sorted[2*k+1])] = True
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
    print(f"Average right set size = {sum_right_sets / n_right_sets}")
    print(f"Within part clauses = {within_part_clauses}")
    print(f"Number of b_clauses = {b_clauses}")
    print(f"Number of direct clauses = {direct_clauses}")
    print(f"Total weight of across-part bicliques: |A|={A_weight}, |S|={S_weight}, total={A_weight + S_weight}")
    print(f"Total clauses = {within_part_clauses + b_clauses + direct_clauses + A_weight + S_weight}?")
    print(f"total clauses in enc = {enc.n_clauses()}")
    print(f"Theoretical prediction for A-weight {n* math.ceil(n / (2*r))}, S-weight {n* (2**r) // 2}")
    # return bicliques
    return enc


def run_depth_reencode(input_formula_filename, timeout=None, output_formula_filename='tmp/depth_reencoded.cnf'):
    start = time.perf_counter_ns()
    encoding = depth_reencode(input_formula_filename=input_formula_filename)
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