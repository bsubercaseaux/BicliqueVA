import math
import numpy as np
import matplotlib.pyplot as plt

def S(n, r):
    # S(n, r) := ceil(n / r) * (r * 2**r // 2)
    # return math.ceil(n / r) * (r * (2 ** r) // 2)
    return n * (2 ** r) // 2

def A(n, r):
    # A(n, r) := ceil(n / r) * ceil(n / 2)
    return n * n / (4 * r)
    # return math.ceil(n / r) * math.ceil(n / 2)

def best_ratio(n):
    """
    Compute (min_r S(n,r)+A(n,r)) / (n^2 / (2 log2 n))
    We search r from 1 up to min(n, ceil(log2 n) + r_extra).
    """
    if n <= 1:
        raise ValueError("n must be >= 2")

    r_max = min(n, int(math.ceil(math.log2(n))) + 4)
    best_val = None

    for r in range(1, r_max + 1):
        val = S(n, r) + A(n, r)
        if best_val is None or val < best_val:
            best_val = val

    denom = (n ** 2) / (2.0 * math.log2(n))
    return best_val / denom

# ---- parameters you can tweak ----
N_MIN = 4
N_MAX = int(1e60)
STEP = N_MAX // 10000  # you can increase this for fewer points (e.g., 5 or 10)

# compute values
n_values = list(range(N_MIN, N_MAX + 1, STEP))
ratios = [best_ratio(n) for n in n_values]

# plot
plt.figure()
plt.plot(n_values, ratios)
plt.xscale('log')          # x axis on log scale
plt.xlabel('n (log scale)')
plt.ylabel('min_r (S(n,r) + A(n,r)) / (n^2 / (2 log2 n))')
plt.title('Ratio vs n')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

plt.tight_layout()
plt.show()
