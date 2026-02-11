import itertools
import random
from typing import Optional
from math import log2
import time

try:
    from ..utils.graph import Graph, erdos_renyi
    from ..utils.formula_to_graph import form2graph
    from ..utils.bicliques_to_formula import biclique_partition_to_formula
except ImportError:
    from graph import Graph, erdos_renyi
    # from bicliques_to_formula import biclique_partition_to_formula
    from formula_to_graph import form2graph


# from ..utils.graph import Graph  
# from ..utils.bicliques_to_formula import biclique_partition_to_formula
# from ..utils.formula_to_graph import form2graph

def lg(x):
    return log2(x)

def find_balanced_biclique(graph: Graph, r_multiplier: float = 10):
    vertices = graph.vertices
    n = graph.n_vertices()

    r = int(r_multiplier * int(lg(n) - 2 * lg(lg(n))))
    n_edges = sum(graph.degree(v) for v in vertices) // 2
    gamma = 2* n_edges / (n * (n - 1))

    vertices_by_degree = sorted(vertices, key=lambda v: graph.degree(v), reverse=True)
    D = vertices_by_degree[:r]
    V_minus_D = [v for v in vertices if v not in D]

    V_star = []
    neighborhoods = {}
    eps = 1/ lg(n)
   #print(f"Parameters: n={n}, r={r}, gamma={gamma:.4f}, eps={eps:.4f}, (1-eps)^2 * gamma * r = {(1-eps)**2 * gamma * r:.4f}")
    for v in V_minus_D:
        neighborhood = []
        for u in D:
            if graph.is_edge(u, v):
                neighborhood.append(u)
            if len(neighborhood) >= (1-eps)**2 * gamma * r:
                V_star.append(v)
                neighborhoods[v] = neighborhood
                break
        if len(V_star) >= gamma * n / lg(n):
            break
    # pigeonhole time 
    neighborhood_buckets = {}
    for v in V_star:
        key = tuple(sorted(neighborhoods[v]))
        if key not in neighborhood_buckets:
            neighborhood_buckets[key] = []
        neighborhood_buckets[key].append(v)
    if len(neighborhood_buckets) == 0:
        return 0, [], []
    max_bucket_index = max(neighborhood_buckets.keys(), key=lambda k: len(neighborhood_buckets[k]))
    U = list(max_bucket_index)
    W = neighborhood_buckets[max_bucket_index]

    # validate biclique
    for u in U:
        for w in W:
            assert graph.is_edge(u, w), "Not a biclique!"
    t = min(len(U), len(W))
    #assert t >= (1-eps)**2 * gamma * r, "Biclique size too small!"
    return min(len(U), len(W)), U, W


# def experiments():
#     ns = [1000, 2000, 4000, 8000, 16000]
#     ps = [0.1, 0.3, 0.5, 0.7, 0.9]
#     r_multipliers = [1,2,5,10,20, 30]
#     data = {}
#     for n in ns:
#         for p in ps:
#             g = random_graph(n, p)
#             for r in r_multipliers:
#                 start_time = time.perf_counter_ns()
#                 size, U, W = find_balanced_biclique(g, r_multiplier=r)
#                 runtime = (time.perf_counter_ns() - start_time)/ 1e6
#                 data[(n, p, r)] = (size, runtime)
#                 print(f"n={n}, p={p:.1f}, r={r} => balanced biclique size: {size}, runtime: {runtime:.2f} ms")

def partition_to_bicliques(graph: Graph, r_multiplier: int = 10):
    partition = []
    total_processed_edges = 0
    while True:
        # print(f"Remaining edges: {graph.n_edges()}")
        if graph.n_edges() < 0:
            print("Negtive edges, something went wrong!")
            break
        size, U, W = find_balanced_biclique(graph, r_multiplier=r_multiplier)
        
        # if size <= 1:
        #     print(f"Found biclique {U}, {W}")
        #     break
        print(f"Reduced {len(U)*len(W)} edges into {len(U)+len(W)} weight biclique ")
        total_processed_edges += len(U)*len(W)
        # print(f"Total processed edges so far: {total_processed_edges}")
        partition.append((U, W))
        # remove edges of the biclique
        for u in U:
            for w in W:
                assert graph.is_edge(u, w), "Not a biclique during removal!"
                if graph.is_edge(u, w):
                    graph.remove_edge(u, w)
                
    for u in graph.vertices:
        for v in graph.adjacency.get(u, []):
            if u < v:
                partition.append(([u], [v]))
    return partition


def run_finder_partition_on_formula(formula_filename):
    start = time.perf_counter_ns()
    graph = form2graph(formula_filename)
    end = time.perf_counter_ns()
    runtime_ms = (end - start) / 1e6 

    bicliques = partition_to_bicliques(graph)
    
    resulting_formula = biclique_partition_to_formula(bicliques, original_formula_filename=formula_filename)
    vars = resulting_formula.n_vars()
    cls = resulting_formula.n_clauses()
    
    end = time.perf_counter_ns()
    runtime_ms = (end - start) / 1e6 
    
    return vars, cls, runtime_ms
    
    
    
    

if __name__ == "__main__":
    # experiments()
    g = form2graph(f"tmp/result_bp.cnf")
    print(f"|V| = {g.n_vertices()}, |E| = {g.n_edges()}")
    t, U, W = find_balanced_biclique(g)
    print(f"t = {t}, |U| = {len(U)}, |W| = {len(W)}")
    # n = 1000
    # p = 0.85
    # g = random_graph(n, p)
    # partition = partition_to_bicliques(g, r_multiplier=3)
    # total_weight = 0
    # for L, R in partition:
    #     if len(L) == 1 :
    #         total_weight += len(R)
    #     elif len(R) == 1:
    #         total_weight += len(L)
    #     else:
    #         total_weight += len(L) + len(R)
    # print(f"Biclique partition size: {len(partition)}, total weight: {total_weight}")
    # size, U, W = find_balanced_biclique(g)
    # print(f"Found balanced biclique of size {size} with U: {U} and W: {W}")