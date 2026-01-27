
try:
    from utils.graph import Graph, erdos_renyi
    from utils.formula_to_graph import form2graph
    from utils.bicliques_to_formula import biclique_partition_to_formula
except ImportError:
    from graph import Graph, erdos_renyi
    from bicliques_to_formula import biclique_partition_to_formula
    from formula_to_graph import form2graph


from typing import Optional
from eznf import modeler
import math
import itertools
import time

def biclique_partition(graph: Graph, forced_r: Optional[int] = None, delta_r = 0):
    n = graph.n_vertices()

    lg = lambda x: math.log2(x) if x > 0 else 0.0
    # print(f"lg(n) = {lg(n)}, lg lg n = {lg(lg(n))}")
    r = max(1, int(math.floor(lg(n) - 1.1 * lg(lg(n)))))    # guard for tiny n
    # print(f"lg(n) - 1.1 lg(lg(n)) = {lg(n) - 1.1 * lg(lg(n))}")
    if forced_r is not None:
        r = forced_r
    else:
        r += delta_r
    # print(f"Using r = {r} for n = {n}")

    # 2–3: R(i, j)
    def R(i, j):
        return (i < j and abs(i - j) % 2 == 0) or (i > j and abs(i - j) % 2 == 1)

    num_blocks = math.ceil(n / r)
    blocks = []
    vertex_to_block = {}
    for i in range(1, num_blocks + 1):
        lo = (i - 1) * r
        hi = min(i * r, n)
        block = graph.vertices[lo:hi]
        blocks.append(block)  # P_i
        for v in block:
            vertex_to_block[v] = i-1

    bicliques = []  # list of (L, R) bicliques
    internal_edges = 0
    for P_i in blocks:
        for (u, v) in itertools.combinations(P_i, 2):
            if graph.is_edge(u, v):
                bicliques.append(((u,), (v,)))
                internal_edges += 1

    # print(f"There are {internal_edges} internal edges in blocks")
    A_weight = 0
    S_weight = 0
    for i, P_i in enumerate(blocks):  
        A = {}  
        for v in graph.vertices:
            j = vertex_to_block[v]
            if R(j, i):
                mask = 0
                for idx, u in enumerate(P_i):
                    if graph.is_edge(u, v):
                        mask |= (1 << idx)   # A(S) ← A(S) ∪ {v}
                if mask != 0:  
                    if mask not in A:
                        A[mask] = set()
                    A[mask].add(v)

        for mask, right_set in A.items():
            S = tuple(u for bit, u in enumerate(P_i) if (mask >> bit) & 1)
            if len(S) > 0 and len(right_set) > 0:
                A_weight += len(right_set)
                S_weight += len(S)
                bicliques.append((tuple(sorted(S)), tuple(sorted(right_set))))

    # print(f"Total weight of across-part bicliques: |A|={A_weight}, |S|={S_weight}, total={A_weight + S_weight}")
    # print(f"Theoretical prediction for A-weight {n* math.ceil(n / (2*r))}, S-weight {n* (2**r) // 2}")
    return bicliques

def biclique_partition_to_maximal_covering(graph, bicliques):
    final_bicliques = []

    def quality(l_size, r_size):
        return l_size * r_size - l_size - r_size - 1
    
    print(f"There are {len(bicliques)} bicliques")

    marked_graph = graph.copy()
    

    for ind, (L, R) in enumerate(bicliques):
        cur_L = set(L)
        cur_R = set(R)
        # print(f"Processing biclique {ind}")
        cnt = 0
        while True:
            cnt += 1
            # print(f"Iteration {cnt}")
            improved = False
            current_q = quality(len(cur_L), len(cur_R))
            # Try to find a vertex to add to L that improves the overall quality
            # note: adding to L restricts cur_R to the neighbors of the new vertex
            for v in graph.vertices:
                if v in cur_L:
                    continue
                
                new_R = cur_R.intersection(graph.neighbors(v))
                if quality(len(cur_L) + 1, len(new_R)) > current_q:
                    cur_L.add(v)
                    cur_R = new_R
                    improved = True
                    break
            
            if not improved:
                break

        final_bicliques.append((tuple(sorted(cur_L)), tuple(sorted(cur_R))))
    print(f"There are {len(final_bicliques)} final bicliques")  
    return final_bicliques

def biclique_partition_components(graph, forced_r: Optional[int] = None, delta_r = 0):
    visited = set()
    subgraphs = []
    for start in graph.vertices:
        if start in visited:
            continue
        stack = [start]
        visited.add(start)
        comp = []
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in graph.vertices:
                if v not in visited and graph.is_edge(u, v):
                    visited.add(v)
                    stack.append(v)
        subgraphs.append(graph.get_induced_subgraph(comp))

    all_bicliques = []
    for g in subgraphs:
        all_bicliques.extend(biclique_partition(g, forced_r=forced_r, delta_r=delta_r))

    # Deduplicate (L,R) pairs
    seen = set()
    unique = []
    for L, R in all_bicliques:
        key = (tuple(L), tuple(R))
        if key not in seen:
            seen.add(key)
            unique.append((L, R))
    return unique



def run_biclique_partition_on_formula(formula_filename, forced_r: Optional[int] = None, delta_r = 0):
    start = time.perf_counter_ns()
    graph = form2graph(formula_filename)
    
    start_inner = time.perf_counter_ns()
    bicliques = biclique_partition(graph, forced_r, delta_r)
    end_inner = time.perf_counter_ns()
    
    resulting_formula = biclique_partition_to_formula(bicliques, original_formula_filename=formula_filename)
    vars = resulting_formula.n_vars()
    cls = resulting_formula.n_clauses()
    
    resulting_formula.serialize(f"tmp/result_bp.cnf")
    end = time.perf_counter_ns()
    runtime_ms = (end - start) / 1e6 
    # runtime_ms = (end_inner - start_inner) / 1e6
    return vars, cls, runtime_ms

def run_biclique_partition_on_formula_plus_covering(formula_filename, forced_r: Optional[int] = None, delta_r = 0):
    start = time.perf_counter_ns()
    graph = form2graph(formula_filename)
    
    start_inner = time.perf_counter_ns()
    bicliques = biclique_partition(graph, forced_r, delta_r)
    bicliques = biclique_partition_to_maximal_covering(graph, bicliques)
    end_inner = time.perf_counter_ns()

    
    resulting_formula = biclique_partition_to_formula(bicliques, original_formula_filename=formula_filename)
    vars = resulting_formula.n_vars()
    cls = resulting_formula.n_clauses()
    
    resulting_formula.serialize(f"../tmp/result_bp.cnf")
    end = time.perf_counter_ns()
    runtime_ms = (end - start) / 1e6 
    # runtime_ms = (end_inner - start_inner) / 1e6
    return vars, cls, runtime_ms

def run_biclique_partition_on_formula_plus_finder(formula_filename, forced_r: Optional[int] = None, delta_r = 0):
    start = time.perf_counter_ns()
    graph = form2graph(formula_filename)
    
    start_inner = time.perf_counter_ns()
    bicliques = biclique_partition(graph, forced_r, delta_r)
    end_inner = time.perf_counter_ns()
    
    resulting_formula = biclique_partition_to_formula(bicliques, original_formula_filename=formula_filename)
    vars = resulting_formula.n_vars()
    cls = resulting_formula.n_clauses()
    
    resulting_formula.serialize(f"tmp/result.cnf")


    end = time.perf_counter_ns()
    runtime_ms = (end - start) / 1e6 
    # runtime_ms = (end_inner - start_inner) / 1e6
    return vars, cls, runtime_ms





if __name__ == "__main__":
    # g = erdos_renyi(1000, 0.5)
    # print("Generated graph,  # edges =", g.n_edges)
    # B = biclique_partition_components(g, delta_r=0)
    formula = "../tmp/random_2cnf_500_0.9.cnf"
    vars, cls, runtime = run_biclique_partition_on_formula(formula)
    print(f"Reencoded formula with {vars} variables and {cls} clauses.")
    print(f"Runtime: {runtime} ms")
    
    vars, cls, runtime = run_biclique_partition_on_formula_plus_covering(formula)
    print(f"Reencoded formula with {vars} variables and {cls} clauses.")
    print(f"Runtime: {runtime} ms")
    # edges = g.n_edges
    # print(f"Graphy with {len(g.vertices)} vertices and {edges} edges")
    # total_weight = 0
    # for i, (L, R) in enumerate(B):
    #     total_weight += min(len(L)*len(R), len(L)+len(R))
    # print(f"Found {len(B)} bicliques with total weight {total_weight}")