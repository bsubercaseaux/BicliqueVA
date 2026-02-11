import math
import random
from collections import defaultdict

# ---------- utilities ----------

def h2(gamma: float) -> float:
    """Binary entropy in base 2: h2(x) = -x log2 x - (1-x) log2(1-x)."""
    if gamma <= 0 or gamma >= 1:
        return 0.0
    return -(gamma * math.log2(gamma) + (1 - gamma) * math.log2(1 - gamma))

def choose_leq_threshold(L: int, y: int, T: int) -> bool:
    """Return True iff C(L, y) <= T (y >= 0, L >= y). Uses fast incremental product."""
    if y < 0 or y > L: 
        return y == 0  # C(L,0)=1 <= T for T>=1 in our calls
    # small shortcuts
    if y == 0 or y == L:
        return 1 <= T
    if y == 1 or y == L - 1:
        return L <= T
    # compute C(L, y) exactly but stop if exceeds T
    y = min(y, L - y)
    num = 1
    for i in range(1, y + 1):
        num = (num * (L - y + i)) // i
        if num > T:
            return False
    return num <= T

def precompute_max_len_per_y(r: int, T: int):
    """
    For each y in [1..r], compute X[y] = max L in [y..r] such that C(L,y) < T.
    (We never use y=0 because empty S are ignored.)
    """
    X = [0] * (r + 1)
    for y in range(1, r + 1):
        # binary search on L in [y .. r]
        lo, hi, ans = y, r, y - 1
        while lo <= hi:
            mid = (lo + hi) // 2
            if choose_leq_threshold(mid, y, T - 1):  # strictly < T
                ans = mid
                lo = mid + 1
            else:
                hi = mid - 1
        X[y] = ans
    return X

# ---------- core algorithm (Theorem 7) ----------

class Graph:
    """Undirected simple graph given by adjacency sets; vertices are 0..n-1."""
    def __init__(self, n):
        self.n = n
        self.adj = [set() for _ in range(n)]
        self.m = 0

    def add_edge(self, u, v):
        if u == v or v in self.adj[u]:
            return
        self.adj[u].add(v)
        self.adj[v].add(u)
        self.m += 1

def biclique_partition_density(G: Graph, gamma: float = None):
    """
    Density-aware biclique partition (Theorem 7).
    Returns a list of bicliques (L, R) with L, R as tuples of vertex ids.

    Time: O(m).  Weight: (1/2 + o(1)) * h2(gamma) * n^2 / log2 n (as in Theorem 7).
    """
    n = G.n
    # estimate density if not provided
    if gamma is None:
        denom = n * (n - 1) // 2
        gamma = (G.m / denom) if denom else 0.0

    # 1) set r = ceil( lg^2 n / h2(gamma) ), guard small cases
    lg = lambda x: math.log2(max(x, 2))
    h = max(h2(gamma), 1e-9)
    print(f"lg(n) = {lg(n):.4f}, h2(gamma) = {h:.4f}, lg(n)^2 = {lg(n)**2:.4f}")
    # r = max(1, math.ceil((lg(n) ** 2) / h))
    r = 8
    print(f"Using r = {r} for n = {n}, gamma = {gamma:.4f}, h2(gamma) = {h:.4f}")

    # 2) partition vertices into blocks P_i of size at most r, and record positions in each block
    parts = []
    pos_in_part = []  # length n; for v -> (i, j) i=block index, j=1..|P_i|
    pos_in_part = [(-1, -1)] * n
    order = list(range(n))  # any fixed labeling works (paper labels (i,j)) :contentReference[oaicite:2]{index=2}
    for i in range(0, n, r):
        block = order[i:i + r]
        parts.append(block)
        bidx = len(parts) - 1
        for j, v in enumerate(block, start=1):  # positions 1..|P_i|
            pos_in_part[v] = (bidx, j)

    t = len(parts)

    # 3) orient parts by a tournament R on [t]; the parity/circular tournament is fine here
    def circ_dist(i, j):
        d = abs(i - j)
        return min(d, t - d)

    def R_rel(i, j):  # True if edge oriented i -> j
        if i == j:
            return False
        if i < j:
            return circ_dist(i, j) % 2 == 0
        else:
            return circ_dist(i, j) % 2 == 1
    # (For Theorem 7 we only minimize weight; almost-regularity is not required.) :contentReference[oaicite:3]{index=3}

    # 4) Precompute threshold and X[y]
    #    T = n / r^4 (paper’s Eq. (4) threshold), and X[y] = max allowed slice length for |S|=y
    T = max(10, n // (r ** 4))
    X = precompute_max_len_per_y(r, T)  # uses strictly "< T" as in Eq. (4) :contentReference[oaicite:4]{index=4}

    # 5) Edge-driven accumulation of neighbor positions per (v, target part i)
    #    For each edge {u,v}, push the position of the endpoint in part P_i to key (other, i) iff R(g(other),i).
    neigh_pos = defaultdict(list)  # key=(v, i) -> sorted list of positions of neighbors of v inside P_i
    within_part_edges = []         # list of (u,v) with g(u)==g(v)
    for u in range(n):
        for v in G.adj[u]:
            if u < v:
                iu, ju = pos_in_part[u]
                iv, jv = pos_in_part[v]
                if iu == iv:
                    within_part_edges.append((u, v))  # emit later as ({u},{v})
                else:
                    if R_rel(iu, iv):
                        neigh_pos[(u, iv)].append(jv)
                    if R_rel(iv, iu):
                        neigh_pos[(v, iu)].append(ju)

    # sort positions for each (v,i)
    for key in neigh_pos:
        neigh_pos[key].sort()

    # 6) Build A_i(a,b,S) buckets using the “slicing” rule and then emit bicliques (S, A)
    #    We key buckets by (i, a, b, mask) where mask encodes S within [a..b].
    Ai = defaultdict(list)

    def pack_mask(a, b, positions):
        """Build a bitmask for S = positions ∩ [a,b]. a,b are 1-based positions; positions sorted."""
        mask = 0
        for p in positions:
            if a <= p <= b:
                mask |= 1 << (p - a)  # bit 0 corresponds to position a
        return mask

    # per-block fast access to vertex ids at positions
    block_arrays = [list(block) for block in parts]  # P_i as list of vertex ids

    for (v, i) in neigh_pos.keys():
        P = neigh_pos[(v, i)]         # sorted neighbor positions in P_i
        mlen = len(parts[i])
        idx = 0
        while idx < len(P):
            a = P[idx]  # next slice starts at the next neighbor position (we skip empty-S slices)
            # find maximal t >= 1 such that (P[idx+t-1] - a + 1) <= X[t]
            # (i.e., we can include t neighbors within allowed length)
            tmax = 1
            # try to extend t while feasible and within bounds
            while idx + tmax < len(P):
                cand_t = tmax + 1
                if cand_t > r:
                    break
                needed_len = P[idx + cand_t - 1] - a + 1
                if needed_len <= X[cand_t]:
                    tmax = cand_t
                else:
                    break
            # length cap also cannot cross the next neighbor if we want exactly tmax neighbors
            next_break = P[idx + tmax] - 1 if (idx + tmax < len(P)) else mlen
            b_cap = a + X[tmax] - 1
            b = min(b_cap, next_break)
            # build mask for the tmax neighbors we’re packing here
            positions = P[idx: idx + tmax]
            mask = pack_mask(a, b, positions)
            if mask != 0:  # S non-empty by construction
                Ai[(i, a, b, mask)].append(v)
            idx += tmax

    bicliques = []

    # Emit across-parts bicliques (S, A_i(a,b,S))
    for (i, a, b, mask), Averts in Ai.items():
        # reconstruct S as actual vertex ids from mask bits over P_i[a..b]
        P_i = block_arrays[i]
        Lverts = []
        for off in range(b - a + 1):
            if (mask >> off) & 1:
                Lverts.append(P_i[a - 1 + off])
        if Lverts and Averts:
            print(f"Emitting biclique with |L|={len(Lverts)}, |R|={len(Averts)} from part {i}, slice [{a},{b}]")
            bicliques.append((tuple(Lverts), tuple(Averts)))

    # Emit within-part single-edge bicliques ({u},{v})
    for (u, v) in within_part_edges:
        bicliques.append(((u,), (v,)))

    return bicliques


if __name__ == "__main__":
    # example usage
    n = 1000
    p = 0.9

    random_graph = Graph(n)
    for u in range(n):
        for v in range(u + 1, n):
            if random.random() < p:
                random_graph.add_edge(u, v)
    print(f"Generated random graph with {n} vertices and {random_graph.m} edges.")

    bicliques = biclique_partition_density(random_graph)
    total_weight = 0
    for L, R in bicliques:
       # print("len(L) =", len(L), "len(R) =", len(R))
        total_weight += len(L) + len(R)
    print(f"Number of bicliques: {len(bicliques)}")
    print(f"Total weight of bicliques: {total_weight}")