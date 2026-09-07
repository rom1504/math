"""Independent integer replay of an actual width optimum with no row repair."""
import itertools
import json

n = 9
code = 134624811
A = [[0 if i == j else 1 for j in range(n)] for i in range(n)]
for e, (i, j) in enumerate(itertools.combinations(range(1, n), 2)):
    A[i][j] = A[j][i] = 1 - 2 * ((code >> e) & 1)

def extrema(B):
    m = len(B)
    energies = [sum(B[i][j] * x[i] * x[j]
                    for i, j in itertools.combinations(range(m), 2))
                for tail in itertools.product((-1, 1), repeat=m-1)
                for x in [(1,) + tail]]
    return min(energies), max(energies)

result = {"n": n, "code": code, "matrix": A, "extrema": extrema(A),
          "row_replacements": []}
for deleted in range(n):
    vv = [i for i in range(n) if i != deleted]
    B = [[A[i][j] for j in vv] for i in vv]
    child_spins = [(1,) + t for t in itertools.product((-1, 1), repeat=n-2)]
    energies = [sum(B[i][j]*x[i]*x[j]
                    for i, j in itertools.combinations(range(n-1), 2))
                for x in child_spins]
    best_cap, best_twice_width, best_width_cap = n*n, n*n, n*n
    hist = {}
    # Fixing the first row sign loses no extension: global row negation
    # is switching the added vertex, which preserves the entire spectrum.
    for tail in itertools.product((-1, 1), repeat=n-2):
        a = (1,) + tail
        fields = [abs(sum(ai*xi for ai, xi in zip(a, x))) for x in child_spins]
        lo = min(q-h for q, h in zip(energies, fields))
        hi = max(q+h for q, h in zip(energies, fields))
        cap, width = max(hi, -lo), hi-lo
        best_cap = min(best_cap, cap)
        if width < best_twice_width:
            best_twice_width, best_width_cap = width, cap
        elif width == best_twice_width:
            best_width_cap = min(best_width_cap, cap)
        key = str((lo, hi))
        hist[key] = hist.get(key, 0) + 1
    result["row_replacements"].append({"deleted_vertex": deleted,
        "best_cap": best_cap, "best_twice_width": best_twice_width,
        "best_cap_at_minimum_width": best_width_cap, "endpoint_histogram": hist})
print(json.dumps(result, indent=2))
