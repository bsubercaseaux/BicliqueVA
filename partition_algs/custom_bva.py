import itertools
from pqdict import pqdict

try:
    from utils.graph import Graph, erdos_renyi
    from utils.formula_to_graph import form2graph
    from utils.bicliques_to_formula import biclique_partition_to_formula
except ImportError:
    from graph import Graph, erdos_renyi
    # from bicliques_to_formula import biclique_partition_to_formula
    from formula_to_graph import form2graph
   
import random
import time
import matplotlib.pyplot as plt
random.seed(42) # ensure reproducibility


CLASSICAL_STRATEGY = 'C'
BEN_STRATEGY = 'B'

class BVA:
    def __init__(self, graph: Graph, verbosity: int = 0, strategy: str = CLASSICAL_STRATEGY):
        self.graph = Graph.copy(graph)
        self.queue = pqdict.maxpq()
        self.vertices = set(graph.vertices)
        self.bicliques = []
        for v in graph.vertices:
            self.queue[v] = (graph.degree(v), v)
        self.history = []
        self.verbosity = verbosity
        self.n_edges = graph.n_edges()
        self.original_vertices = list(graph.vertices)
        self.next_v = max(self.vertices) + 1 if self.vertices else 0

        self.depth_history = []
        self.depth = {}
        for v in graph.vertices:
            self.depth[v] = 0


        assert strategy in [CLASSICAL_STRATEGY, BEN_STRATEGY], f"Unknown strategy {strategy}"
        self.strategy = strategy
        # print(f"BVA initialized with strategy {self.strategy} on graph with {self.graph.n_vertices()} vertices and {self.graph.n_edges()} edges.")

    def run(self):
        it = 0
        while len(self.queue):
            it += 1
            if it % 100 == 0:
                print(f"  BVA iteration {it}, remaining edges: {self.n_edges}")
            should_continue =self.dumb_step()
            if not should_continue:
                break

        
        for u in range(len(self.original_vertices)):
            for v in range(u + 1, len(self.original_vertices)):
                if self.graph.is_edge(self.original_vertices[u], self.original_vertices[v]):
                    self.bicliques.append(((self.original_vertices[u],), (self.original_vertices[v],)))
                    # print(f"Edge between original vertices {self.original_vertices[u]} and {self.original_vertices[v]} still exists!")
        return self.graph

    def augment_biclique(self, L, R, u, new_R):
        L_new = L | {u}
        new_quality = len(L_new) * len(new_R) - (len(L_new) + len(new_R)) - 1
        return L_new, new_R, new_quality
    
    def replace_biclique(self, L, R):
        self.bicliques.append((L, R))
        new_vertex = self.next_v
        self.next_v += 1
        self.graph.vertices.append(new_vertex)
        self.graph.vertices.append(-new_vertex)
        self.vertices.add(new_vertex)
        self.vertices.add(-new_vertex)
        max_depth = -1
        
        for v in L:
            for u in R:
                self.graph.remove_edge(v, u)
                self.n_edges -= 1
            self.graph.add_edge(v, new_vertex)
            self.n_edges += 1
            if v in self.queue:
                self.queue[v] = (self.queue[v][0] - len(R) + 1, v)

        for u in R:
            self.graph.add_edge(-new_vertex, u)
            self.n_edges += 1
            if u in self.queue:
                self.queue[u] = (self.queue[u][0] - len(L) + 1, u)

        
        self.queue[new_vertex] = (self.graph.degree(new_vertex), new_vertex) # self.graph.degree(new_vertex)
        self.queue[-new_vertex] = (self.graph.degree(-new_vertex), -new_vertex) # self.graph.degree(-new_vertex)
        self.history.append((L, R, new_vertex))


        for v in L.union(R):
            if self.depth[v] > max_depth:
                max_depth = self.depth[v]
        self.depth[new_vertex] = max_depth + 1
        self.depth[-new_vertex] = max_depth + 1
        self.depth_history.append(max_depth + 1)

        return new_vertex

    def print_history(self):
        step_num = 1
        for (L, R, new_vertex) in self.history:
            print(f"[Step {step_num}] Biclique L: {sorted(L)}, R: {sorted(R)} replaced with new vertex {new_vertex}")
            step_num += 1
        

    def print_history_summary(self):
        step_num = 1
        biclique_size = {}
        for (L, R, new_vertex) in self.history:
            bic_size = (min(len(L), len(R)), max(len(L), len(R)))
            if bic_size not in biclique_size:
                biclique_size[bic_size] = 1
            else:
                biclique_size[bic_size] += 1
        print("Biclique size summary (min_size, max_size): count")
        for bic_size, count in sorted(biclique_size.items()):
            print(f"  {bic_size}: {count}")
            # print(f"[Step {step_num}] Biclique L: {L}, R: {R} replaced with new vertex {new_vertex}")
            # step_num += 1

    def print_depth_stats(self):
        D = {}
        for d in self.depth_history:
            if d not in D:
                D[d] = 1
            else:
                D[d] += 1
        print("Depth statistics (depth: count):")
        for d, count in sorted(D.items()):
            print(f"  {d}: {count}")

    def max_codegree(self, L, R):
        max_u = None
        max_codeg = -1
        best_new_R = None
        
        if not R:
            return None, 0, set()
            
        candidates = set()
        for r in R:
            candidates.update(self.graph.neighbors(r))
        candidates -= L

        # l = L.pop()
        # # depth_L = self.depth[l]
        # L.add(l)
        # final_candidates = set()
        # for u in candidates:
        #     # if self.depth[u] == depth_L:
        #     if True:
        #         final_candidates.add(u)
        
        old_quality = len(L) * len(R) - (len(L) + len(R)) 
        def new_quality(codeg):
            return (len(L) + 1) * codeg - (len(L) + 1 + codeg) 
        for u in candidates:
            neighbors_u = set(self.graph.neighbors(u))
            new_R = R & neighbors_u
            codeg = len(new_R)
            # test for direct quality improvement
            if new_quality(codeg) > old_quality:
                return u, codeg, new_R
            if codeg > max_codeg:
                max_codeg = codeg
                max_u = u
                best_new_R = new_R
        return max_u, max_codeg, best_new_R
    
    def dumb_step(self):
        k_24s = []
        # pairs = list(itertools.combinations(self.graph.vertices, 2))
        # random.shuffle(pairs)
        for u, v in itertools.combinations(self.graph.vertices, 2):
            neighbors_u = set(self.graph.neighbors(u))
            neighbors_v = set(self.graph.neighbors(v))
            common_neighbors = neighbors_u & neighbors_v
            if len(common_neighbors) >= 4:
                random_4 = random.sample(list(common_neighbors), 4)
                self.replace_biclique(set([u, v]), set(random_4))
                return True
        # print(f"  Found {len(k_24s)} K_{2,4} bicliques.")
        if not k_24s:
            return False
        # L, R = random.choice(k_24s)
        # self.replace_biclique(set(L), set(R))
        # return True

    def step(self):
        v, (deg_v, _) = self.queue.popitem()
        L = {v}
        R = set(self.graph.neighbors(v)) 
        quality = 0
        best_quality = quality   
        best_L, best_R = L, R

        while True:
            u, codeg_u, new_R = self.max_codegree(L, R)
            # if u is None or (self.strategy == BEN_STRATEGY and codeg_u <= 1):
            #     break
            if self.verbosity > 0:
                print(f"  Considering adding vertex {u} with codegree {codeg_u} to biclique L: {L}, R: {R}")
            
            new_L, new_R, new_quality = self.augment_biclique(L, R, u, new_R)
           
            if new_quality > best_quality:
                best_quality = new_quality
                best_L, best_R = new_L, new_R
                L = new_L
                R = new_R
            elif self.strategy == CLASSICAL_STRATEGY: # break if no improvement
                break
            elif self.strategy == BEN_STRATEGY: # augment anyways 
                L = new_L
                R = new_R
    
        if len(L) == 1:
            return
     
        self.queue[v] = (deg_v, v)

        if best_quality > 0:
            self.replace_biclique(best_L, best_R)



def experiment(n = 200, p = 0.5, verbosity: int = 0, strategy: str = CLASSICAL_STRATEGY):
    G = erdos_renyi(n, p)
    if verbosity > 0:
        print(f"Generated random graph with {G.n_vertices()} vertices and {G.n_edges()} edges.")
    B = BVA(G, verbosity=verbosity, strategy=strategy)
    resulting_graph = B.run()
    
    final_size = B.n_edges
    if verbosity > 0:
        print(f"Original size: {G.n_edges()} , Final size: {final_size}")
    # B.print_history_summary()
    # B.print_depth_stats()
    # B.print_history()
    return G.n_edges(), final_size


def experiments():
    ns = list(range(100, 400, 25))
    ps = [0.1, 0.3, 0.5]
    strategies = [CLASSICAL_STRATEGY]

    DATA = []
    for n in ns:
        for p in ps:
            for strategy in strategies:
                original_size, final_size = experiment(n=n, p=p, verbosity=0, strategy=strategy)
                reduction = original_size / final_size
                print(f"n={n}, p={p}, strategy={strategy}: original size={original_size}, final size={final_size}, reduction={reduction:.2f}%")
                DATA.append((n, p, strategy, original_size, final_size, reduction))
    # Plotting
    plt.figure(figsize=(12, 8))
    for strategy in strategies:
        for p in ps:
            xs = [data[0] for data in DATA if data[2] == strategy and data[1] == p]
            ys = [data[5] for data in DATA if data[2] == strategy and data[1] == p]
            marker = 'o' if strategy == CLASSICAL_STRATEGY else 's'
            plt.plot(xs, ys, marker=marker, label=f"p={p}, strategy={strategy}")
    plt.title("BVA Reduction Percentage vs Number of Vertices")
    plt.xlabel("Number of Vertices (n)")
    plt.ylabel("Reduction Percentage (%)")
    plt.legend()
    plt.grid(True)
    plt.savefig("bva_experiment_results.png")
    plt.show()



def run_custom_bva_on_formula(formula_filename: str, output_file_path: str):
   
    g = form2graph(formula_filename)
    original_edges = g.n_edges()
    start_time = time.perf_counter_ns()
    B = BVA(g)
    B.run()
    final_edges = B.n_edges
    # print(f"Custom BVA reduced edges from {original_edges} to {final_edges}, reduction: {100.0 * (original_edges - final_edges) / original_edges:.2f}%")
    # factored_edges = sum([len(L) * len(R) for (L, R) in B.bicliques])
    final_vertices = len(B.graph.vertices)

    # aux_vars = (final_vertices - g.n_vertices())//2
    # n_steps = len(B.history)
    # print(f"aux_vars = {aux_vars}, n_steps = {n_steps}")

    # print(f"Final formula has {final_vertices} variables ({aux_vars} auxiliary) and {final_edges} clauses (original had {g.n_vertices()} variables and {original_edges} clauses)."  )

    biclique_formula = biclique_partition_to_formula(B.bicliques, formula_filename)
    biclique_formula.serialize(output_file_path)
    return biclique_formula.n_vars(), final_edges, (time.perf_counter_ns() - start_time) / 1e6
    

if __name__ == "__main__":
    # experiments()
    start = time.perf_counter_ns()
    ratios = []
    for n in range(10, 301, 20):
        original_res, res = experiment(n=n, p=0.5, verbosity=1, strategy=CLASSICAL_STRATEGY)
        print(f"n={n}: original size={original_res}, final size={res}, reduction={original_res / res:.2f}%")
        ratios.append(original_res / res)

    ns = list(range(10, 301, 20))
    plt.figure(figsize=(10, 6))
    plt.plot(ns, ratios, marker='o', linestyle='-')
    plt.title("Reduction Ratio (Original Size / Final Size) vs n")
    plt.xlabel("Number of Vertices (n)")
    plt.ylabel("Reduction Ratio")
    plt.grid(True)
    plt.show()
    # original_res, res = experiment(n=80, p=0.5, verbosity=0, strategy=CLASSICAL_STRATEGY)
    end = time.perf_counter_ns()

    print(f"Time: {(end - start) / 1e9} s, result: {res}, original: {original_res}")
