import numpy as np
from ortools.sat.python import cp_model


def walsh(m):
    q = 1 << m
    return np.array(
        [[1 if (bin(i & j).count("1") & 1) == 0 else -1 for j in range(q)] for i in range(q)],
        dtype=np.int64,
    )


def mats(m, a):
    R = walsh(m)
    q = len(R)
    d = np.array(
        [1 if (bin(a & v).count("1") & 1) == 0 else -1 for v in range(q)],
        dtype=np.int64,
    )
    C = np.kron(R, (d[:, None] * R) * d[None, :])
    W = np.kron(R, R)
    return C, W


def solve(m, labels, edges, roots=(), limit=30):
    cs = []
    w = None
    for a in labels:
        c, w = mats(m, a)
        cs.append(c)
    n = len(w)
    total = n * len(labels)
    coef = {}
    constant = 0

    def add(i, j, c):
        nonlocal constant
        if i == j:
            constant += c
        else:
            if i > j:
                i, j = j, i
            coef[(i, j)] = coef.get((i, j), 0) + int(c)

    for block, c in enumerate(cs):
        for i in range(n):
            for j in range(i, n):
                add(block * n + i, block * n + j, int(c[i, j]) if i < j else int(c[i, i]) // 2)
    for u, v in edges:
        for i in range(n):
            for j in range(n):
                add(u * n + i, v * n + j, int(w[i, j]))

    model = cp_model.CpModel()
    bits = [model.NewBoolVar(f"b{i}") for i in range(total)]
    terms = []
    c0 = constant
    for (i, j), c in coef.items():
        if not c:
            continue
        z = model.NewBoolVar(f"z{i}_{j}")
        model.Add(z >= bits[i] - bits[j])
        model.Add(z >= bits[j] - bits[i])
        model.Add(z <= bits[i] + bits[j])
        model.Add(z <= 2 - bits[i] - bits[j])
        c0 += c
        terms.append(-2 * c * z)

    omega = (1 << m) - 1
    y = np.array(
        [
            1
            if ((bin(u & v).count("1") + bin(omega & u).count("1")) & 1) == 0
            else -1
            for u in range(1 << m)
            for v in range(1 << m)
        ],
        dtype=np.int64,
    )
    for block, mult in roots:
        h = w @ y
        for i, c in enumerate(h):
            c = int(c) * int(mult)
            c0 += c
            terms.append(-2 * c * bits[block * n + i])
    model.Maximize(sum(terms) + c0)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = limit
    solver.parameters.num_search_workers = 8
    status = solver.Solve(model)
    return solver.StatusName(status), solver.ObjectiveValue(), solver.BestObjectiveBound(), solver.WallTime()


if __name__ == "__main__":
    for labels in ([4, 2], [4, 7]):
        for edges, roots in [([(0, 1)], ()), ([(0, 1)], ((0, 1),)), ([], ((0, 1), (1, 1)))]:
            print(labels, edges, roots, solve(3, labels, edges, roots, 60))
