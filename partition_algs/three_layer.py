"""
three_layer_hybrid_reencoding.py

Hybrid reencoding of the independent-set constraints of a CNF formula:

  - Parse original.cnf, keep all clauses EXCEPT binary negative ones (-u -v 0).
  - Use form2graph(original.cnf) to get the conflict graph G on variables.
  - Partition vertices into blocks and build "patterns" à la Nechiporuk.
  - For each biclique S × R coming from a pattern:
        * pairwise cost  = |S| * |R|
        * 3-layer cost   = |S| + |R| + 1
    Use 3-layer iff 3-layer cost < pairwise cost.
  - Add:
        * pairwise clauses (-u ∨ -v) for edges assigned to pairwise;
        * 3-layer clauses

              (¬x_u  ∨ A)
              (¬B     ∨ ¬x_v)
              (¬A     ∨ ¬B)

          for the bicliques assigned to 3-layer.

Guarantee: on the part of the formula coming from binary (-u -v) clauses,
the hybrid encoding never uses more clauses than the original pairwise encoding
(and often fewer).
"""

import math
from typing import Dict, Set, Tuple, List, Optional

try:
    # Your project layout
    from utils.graph import Graph
    from utils.formula_to_graph import form2graph
except ImportError:
    # Flat layout (as in the uploaded files)
    from graph import Graph
    from formula_to_graph import form2graph

from eznf import modeler


# --------------------------------------------------------------------
#  Block size choice
# --------------------------------------------------------------------

def _choose_block_size(n: int) -> int:
    """
    Block size r ≈ log n − 1.1 log log n, with guards for small n.
    This is a heuristic choice; you can override via forced_r.
    """
    if n <= 1:
        return 1
    lg = lambda x: math.log2(x) if x > 0 else 0.0
    r = int(math.floor(lg(n) - 1.1 * lg(lg(n))))
    return max(1, r)


# --------------------------------------------------------------------
#  Nechiporuk-style partition: edges -> bicliques or pairwise
# --------------------------------------------------------------------

def three_layer_partition(
    graph: Graph,
    forced_r: Optional[int] = None,
    delta_r: int = 0,
) -> Tuple[Dict[str, Set[int]],
           Dict[str, Set[int]],
           List[Tuple[str, str]],
           Set[Tuple[int, int]]]:
    """
    Partition the edges of 'graph' into:

      * a family of 3-layer bicliques:
            left_groups:  dict[L_id -> set of vertices]
            right_groups: dict[R_id -> set of vertices]
            index_edges:  list of (L_id, R_id)

        where each biclique S×R corresponding to some (L_id,R_id) is encoded via:
            (¬x_u ∨ L_id)    for u in S
            (¬R_id ∨ ¬x_v)   for v in R
            (¬L_id ∨ ¬R_id)

      * a set of edges that remain pairwise:
            pairwise_edges: set of (u,v) with u < v

    Construction:

      - Vertices are partitioned into contiguous blocks P_0, …, P_{t-1}
        of size ≈ r.
      - For each vertex v and each block P_i we build a pattern

            S_i(v) = { u in P_i : {u,v} ∈ E and u < v and block[u] == i }

        i.e., edges are oriented from the smaller endpoint into the pattern.
      - For each (i, S) with S_i(v)=S ≠ ∅, we get a biclique S × R_{i,S}
        where R_{i,S} = { v : S_i(v) = S }.
      - For that biclique, we compare pairwise cost |S||R_{i,S}| with
        3-layer cost |S| + |R_{i,S}| + 1 and choose the cheaper.

    By construction each undirected edge {u,v} is covered by exactly one such
    pattern (the block of min(u,v)), so sums of costs over patterns are exact.
    """

    vertices = sorted(graph.vertices)
    n = len(vertices)
    if n == 0:
        return {}, {}, [], set()

    # Block size and partition
    if forced_r is not None:
        r = forced_r
    else:
        r = _choose_block_size(n)
    r += delta_r
    if r <= 0:
        r = 1

    parts: List[List[int]] = []
    for i in range(0, n, r):
        parts.append(vertices[i:i + r])

    block_of: Dict[int, int] = {v: i for i, part in enumerate(parts) for v in part}

    # (part_idx, S) -> set of "right" vertices v
    pattern_to_right: Dict[Tuple[int, frozenset], Set[int]] = {}

    # Assign each edge {u,v} to exactly one pattern:
    # for edge {u,v} with u<v and block_of[u]=i, it contributes u to S_i(v).
    for v in vertices:
        v_id = v
        for i, part in enumerate(parts):
            S: Set[int] = set()
            for u in part:
                if graph.is_edge(u, v_id) and u < v_id and block_of[u] == i:
                    S.add(u)
            if not S:
                continue
            key = (i, frozenset(S))
            if key not in pattern_to_right:
                pattern_to_right[key] = set()
            pattern_to_right[key].add(v_id)

    left_groups: Dict[str, Set[int]] = {}
    right_groups: Dict[str, Set[int]] = {}
    index_edges: List[Tuple[str, str]] = []
    pairwise_edges: Set[Tuple[int, int]] = set()

    left_key_to_id: Dict[Tuple[int, frozenset], str] = {}
    right_key_to_id: Dict[Tuple[int, frozenset], str] = {}
    left_ctr = 0
    right_ctr = 0

    for key, R in pattern_to_right.items():
        _, S_frozen = key
        S = set(S_frozen)

        s = len(S)
        t = len(R)

        pair_cost = s * t
        three_cost = s + t + 1

        if three_cost < pair_cost:
            # Use 3-layer encoding for this biclique S × R
            if key not in left_key_to_id:
                left_ctr += 1
                L_id = f"A_{left_ctr}"
                left_key_to_id[key] = L_id
                left_groups[L_id] = S
            else:
                L_id = left_key_to_id[key]

            if key not in right_key_to_id:
                right_ctr += 1
                R_id = f"B_{right_ctr}"
                right_key_to_id[key] = R_id
                right_groups[R_id] = set(R)
            else:
                R_id = right_key_to_id[key]

            index_edges.append((L_id, R_id))
        else:
            # Keep pairwise encoding for these edges
            for u in S:
                for v in R:
                    if u == v:
                        continue
                    e = (u, v) if u < v else (v, u)
                    pairwise_edges.add(e)

    return left_groups, right_groups, index_edges, pairwise_edges


# --------------------------------------------------------------------
#  CNF builder: hybrid 2-layer + 3-layer ISP encoding
# --------------------------------------------------------------------

def hybrid_three_layer_isp_encoding_from_graph(
    graph: Graph,
    original_formula_filename: Optional[str] = None,
    forced_r: Optional[int] = None,
    delta_r: int = 0,
) -> modeler.Modeler:
    """
    Build a hybrid CNF encoding of the independent-set constraints
    of 'graph':

      - If 'original_formula_filename' is given:
          * read that DIMACS CNF;
          * copy *all* non-binary clauses;
          * drop all binary (-u -v) clauses;
          * then add our hybrid ISP encoding for the resulting conflict graph.

      - If no original formula is given:
          * create one variable per vertex in graph.vertices;
          * only add the ISP part (pairwise + 3-layer).

    The ISP part is never asymptotically worse (in #clauses)
    than the original pairwise encoding and is often smaller.
    """

    enc = modeler.Modeler()

    # 1) Copy original formula if provided, dropping (-u -v) clauses
    max_var_in_original = 0
    if original_formula_filename is not None:
        with open(original_formula_filename, "r") as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("c"):
                    continue
                if stripped.startswith("p"):
                    toks = stripped.split()
                    if len(toks) >= 4 and toks[1] == "cnf":
                        n_vars = int(toks[2])
                        max_var_in_original = max(max_var_in_original, n_vars)
                        # Declare base variables 1..n_vars
                        for v in range(1, n_vars + 1):
                            enc.add_var(v)
                    continue

                # Clause line: ints ending with 0
                toks = stripped.split()
                lits = [int(t) for t in toks if t not in ("", "0")]
                if not lits:
                    continue

                # Drop binary negative clauses: they'll be reencoded
                if len(lits) == 2 and lits[0] < 0 and lits[1] < 0:
                    continue

                enc.add_clause(lits)
    else:
        # No original CNF: declare variables from the graph
        for v in graph.vertices:
            enc.add_var(v)
        max_var_in_original = max(graph.vertices) if graph.vertices else 0

    # 2) Nechiporuk-style partition of edges into 3-layer vs pairwise
    left_groups, right_groups, index_edges, pairwise_edges = three_layer_partition(
        graph, forced_r=forced_r, delta_r=delta_r
    )

    # 3) Add auxiliary vars for each left/right type
    for L_id in left_groups.keys():
        enc.add_var(L_id)
    for R_id in right_groups.keys():
        enc.add_var(R_id)

    # 4) Base → left-type: (¬x_u ∨ A_L) for all u in each left group
    for L_id, S in left_groups.items():
        for u in S:
            # pattern: [-u, L_id]
            enc.add_clause([-u, L_id])

    # 5) Right-type → base: (¬B_R ∨ ¬x_v) for all v in each right group
    for R_id, T in right_groups.items():
        for v in T:
            # negative aux literal as string "-R_id"
            enc.add_clause([f"-{R_id}", -v])

    # 6) Type–type clauses: (¬A_L ∨ ¬B_R) for each index edge
    for (L_id, R_id) in index_edges:
        enc.add_clause([f"-{L_id}", f"-{R_id}"])

    # 7) Pairwise edges that stayed 2-layer: (¬x_u ∨ ¬x_v)
    for (u, v) in pairwise_edges:
        enc.add_clause([-u, -v])

    return enc


def hybrid_three_layer_reencode_formula(
    original_formula_filename: str,
    forced_r: Optional[int] = None,
    delta_r: int = 0,
) -> modeler.Modeler:
    """
    Convenience wrapper:

      - Build conflict graph via form2graph(original_formula_filename).
      - Build the hybrid encoding.
      - Return the Modeler object.
    """
    g = form2graph(original_formula_filename)
    enc = hybrid_three_layer_isp_encoding_from_graph(
        g,
        original_formula_filename=original_formula_filename,
        forced_r=forced_r,
        delta_r=delta_r,
    )
    return enc


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 3:
        print(
            "Usage: python three_layer_hybrid_reencoding.py "
            "<input.cnf> <output.cnf> [forced_r] [delta_r]"
        )
        sys.exit(1)

    in_cnf = sys.argv[1]
    out_cnf = sys.argv[2]
    forced_r = int(sys.argv[3]) if len(sys.argv) >= 4 else None
    delta_r = int(sys.argv[4]) if len(sys.argv) >= 5 else 0

    # Build hybrid encoding
    g = form2graph(in_cnf)
    enc = hybrid_three_layer_isp_encoding_from_graph(
        g, original_formula_filename=in_cnf, forced_r=forced_r, delta_r=delta_r
    )

    # Adjust this to your Modeler API
    try:
        enc.serialize(out_cnf)
    except AttributeError:
        try:
            enc.serialize(out_cnf)
        except AttributeError:
            raise RuntimeError(
                "Don't know how to write DIMACS for this Modeler. "
                "Replace the last lines with the appropriate call."
            )
