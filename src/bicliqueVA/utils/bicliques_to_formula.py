from eznf import modeler
from typing import Optional 

def biclique_partition_to_formula(bicliques, original_formula_filename: Optional[str] = None):
    enc = modeler.Modeler()
    if original_formula_filename:
        with open(original_formula_filename, 'r') as f:
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
    # print(f"Clauses thus far: {enc.n_clauses()}")

    b_ind = 0


    for L, R in bicliques:
        if len(L)*len(R) - len(L) - len(R) < 1:
            for u in L:
                for v in R:
                    enc.add_clause([-u, -v])
        else:
            enc.add_var(f"y_{b_ind}")
            for u in L:
                enc.add_clause([-u, f"y_{b_ind}"])
            for v in R:
                enc.add_clause([-v, f"-y_{b_ind}"])
            b_ind += 1
       
    return enc