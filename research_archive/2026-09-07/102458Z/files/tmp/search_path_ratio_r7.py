from functools import lru_cache
from itertools import product
import random
import sys


def evaluate(A):
    n = len(A)

    def energy(vertices, spin):
        return sum(
            A[i][j] * spin[r] * spin[s]
            for r, i in enumerate(vertices)
            for s, j in enumerate(vertices)
            if r != s
        )

    @lru_cache(None)
    def endpoints(vertices):
        spins = tuple(product((-1, 1), repeat=len(vertices)))
        vals = tuple(energy(vertices, x) for x in spins)
        hi, lo = max(vals), min(vals)
        return (
            hi,
            -lo,
            tuple(x for x, v in zip(spins, vals) if v == hi),
            tuple(x for x, v in zip(spins, vals) if v == lo),
        )

    def cross_l1(vertices, p, child_positions):
        child = set(child_positions)
        other = [j for j in range(len(vertices)) if j not in child]
        return sum(
            abs(sum(A[vertices[i]][vertices[j]] * p[i] for i in child))
            for j in other
        )

    @lru_cache(None)
    def best_path(vertices):
        if len(vertices) <= 1:
            return 0
        _, _, ps, ns = endpoints(vertices)
        best = -10**9
        for p in ps:
            for neg in ns:
                sides = tuple(
                    tuple(i for i in range(len(vertices)) if p[i] * neg[i] == bit)
                    for bit in (-1, 1)
                )
                if any(not side for side in sides):
                    continue
                for positions in sides:
                    child_vertices = tuple(vertices[i] for i in positions)
                    h = energy(child_vertices, tuple(p[i] for i in positions))
                    fields = [
                        sum(A[vertices[i]][vertices[j]] * p[i] for i in positions)
                        for j in range(len(vertices))
                        if j not in positions
                    ]
                    L = sum(abs(v) for v in fields)
                    CP, CN, _, _ = endpoints(child_vertices)
                    Q = max(CP, CN)
                    reward = max(0, 2 * L - Q + abs(h))
                    best = max(best, reward + best_path(child_vertices))
        return best

    root = tuple(range(n))
    P, N, _, _ = endpoints(root)
    return P, N, best_path(root)


def random_matrix(n, rng):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            A[i][j] = A[j][i] = rng.choice((-1, 1))
    return tuple(tuple(row) for row in A)


def main():
    n = int(sys.argv[1])
    trials = int(sys.argv[2])
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 13579
    rng = random.Random(seed)
    worst = (10,)
    for trial in range(trials):
        A = random_matrix(n, rng)
        P, N, V = evaluate(A)
        row = (V / (P + N), P, N, V, trial, A)
        if row[0] < worst[0]:
            worst = row
            print("new", worst[:5], flush=True)
    print("FINAL", worst[:5])
    for row in worst[5]:
        print(row)


if __name__ == "__main__":
    main()
