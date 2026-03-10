from .graph import Graph

# for now only considers clauses of the form -x -y.
def form2graph(formula_filename):
    with open(formula_filename, 'r') as f:
        header = f.readline()
        tokens = header.split(' ')
        n_vars = int(tokens[2])
        
        edges = []
        
        for line in f:
            tokens = line.split(' ')
            if line.startswith('c'):
                continue
            if len(tokens) == 3:
                if int(tokens[0]) < 0 and int(tokens[1]) < 0:
                    u = abs(int(tokens[0]))
                    v = abs(int(tokens[1]))
                    edges.append((u, v))
    vertices = list(range(1, n_vars + 1))
    g = Graph(vertices=vertices)
    for u, v in edges:
        g.add_edge(u, v)
    return g


