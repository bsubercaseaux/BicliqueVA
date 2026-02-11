from eznf import modeler
from .graph import Graph, erdos_renyi

def graph_to_formula(graph: Graph):
    n = graph.n_vertices()
    enc = modeler.Modeler()
    for i in range(n):
        enc.add_var(f"x{i}")

    for u in range(n):
        for v in range(u + 1, n):
            if graph.is_edge(u, v):
                enc.add_clause([f"-x{u}", f"-x{v}"])
    return enc

def graph_to_formula_is(graph: Graph, is_bound):
    n = graph.n_vertices()
    enc = modeler.Modeler()
    for i in range(n):
        enc.add_var(f"x{i}")

    for u in range(n):
        for v in range(u + 1, n):
            if graph.is_edge(u, v):
                enc.add_clause([f"-x{u}", f"-x{v}"])
        
    enc.at_least_k([f"x{i}" for i in range(n)], is_bound)
    
    return enc

def random_2cnf(n, p):
    g = erdos_renyi(n, p)
    return graph_to_formula(g)

def clique_2cnf(n):
    g = Graph(list(range(n)))
    for i in range(n):
        for j in range(i+1, n):
            g.add_edge(i, j)
    return graph_to_formula(g)

if __name__ == "__main__":
    n = 20
    p = 0.3
    f = clique_2cnf(n)
    f.serialize(f"clique_{n}.cnf")