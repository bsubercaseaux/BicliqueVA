
import math
from .partition_algs.cpp.density_aware import precompute_max_len_per_y, h2

n = 1000
gamma = 0.6
lg = lambda x: math.log2(max(x, 2))
h = max(h2(gamma), 1e-9)
r = max(1, math.ceil((lg(n) ** 2) / h))
print(f"r = {r}")

T = max(1, n // (r ** 4))
print(f"T = {T}")

X = precompute_max_len_per_y(r, T)
print(f"X[1] = {X[1]}")
print(f"X[2] = {X[2]}")
print(f"X[10] = {X[10]}")
