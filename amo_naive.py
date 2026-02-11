from eznf import modeler
import itertools

enc = modeler.Modeler()
N = 10
for i in range(N):
    enc.add_var(f"x_{i}")

for i, j in itertools.combinations(range(N), 2):
    enc.add_clause([f"-x_{i}", f"-x_{j}"])

enc.serialize("amo_naive.cnf")
