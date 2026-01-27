import random

class Graph:
    def __init__(self, vertices=None):
        self.vertices = vertices if vertices is not None else []
        self._n_vertices = len(self.vertices)
        self.adjacency = {}
        self._n_edges = 0

    def add_edge(self, u, v):
        if self.is_edge(u, v):
            return  # edge already exists
        if u not in self.adjacency:
            self.adjacency[u] = set()
        if v not in self.adjacency:
            self.adjacency[v] = set()
        
        self.adjacency[u].add(v)
        self.adjacency[v].add(u)
        self._n_edges += 1

    def remove_edge(self, u, v):
        if not self.is_edge(u, v):
            return
        if u in self.adjacency and v in self.adjacency[u]:
            self.adjacency[u].remove(v)
        if v in self.adjacency and u in self.adjacency[v]:
            self.adjacency[v].remove(u)
        
        self._n_edges -= 1

    def is_edge(self, u, v):
        return u in self.adjacency and v in self.adjacency[u]

    def degree(self, v):
        return len(self.adjacency.get(v, []))

    def n_vertices(self):
        return len(self.vertices)

    def n_edges(self):
        return self._n_edges

    def neighbors(self, v):
        return self.adjacency.get(v, [])
    
    def to_adj_list(self):
        adj = []
        for v in self.vertices:
            adj.append(list(self.adjacency.get(v, [])))
        return adj
    
    def get_induced_subgraph(self, vertex_subset):
        subgraph = Graph(vertices=vertex_subset)
        vertex_set = set(vertex_subset)
        for u in vertex_subset:
            for v in self.adjacency.get(u, []):
                if v in vertex_set:
                    subgraph.add_edge(u, v)
        return subgraph
    
    @staticmethod
    def copy(graph):
        new_graph = Graph(vertices=list(graph.vertices))
        for u in graph.vertices:
            for v in graph.adjacency.get(u, []):
                new_graph.add_edge(u, v)
        return new_graph

def erdos_renyi(n, p):
    """Generates a random graph with n vertices and edge probability p."""
    g = Graph(vertices=list(range(n)))
    for u in range(n):
        for v in range(u + 1, n):
            if random.random() < p:
                g.add_edge(u, v)
    return g