#!/usr/bin/env python3
"""Exact verifier for the Wave 20 puncture-cycle memo."""

import itertools


A6 = [
    [0, 1, 1, 1, 1, 1],
    [1, 0, -1, -1, 1, 1],
    [1, -1, 0, 1, -1, 1],
    [1, -1, 1, 0, 1, -1],
    [1, 1, -1, 1, 0, -1],
    [1, 1, 1, -1, -1, 0],
]

A8 = [
    [0, 1, 1, 1, 1, 1, 1, 1],
    [1, 0, 1, -1, 1, 1, -1, -1],
    [1, 1, 0, 1, -1, 1, -1, -1],
    [1, -1, 1, 0, -1, -1, -1, 1],
    [1, 1, -1, -1, 0, -1, 1, -1],
    [1, 1, 1, -1, -1, 0, 1, 1],
    [1, -1, -1, -1, 1, 1, 0, 1],
    [1, -1, -1, 1, -1, 1, 1, 0],
]

A9 = [
    [0, 1, -1, 1, 1, 1, 1, -1, -1],
    [1, 0, 1, -1, -1, -1, 1, -1, -1],
    [-1, 1, 0, 1, 1, 1, 1, 1, -1],
    [1, -1, 1, 0, 1, 1, 1, -1, 1],
    [1, -1, 1, 1, 0, 1, -1, 1, -1],
    [1, -1, 1, 1, 1, 0, -1, -1, -1],
    [1, 1, 1, 1, -1, -1, 0, 1, 1],
    [-1, -1, 1, -1, 1, -1, 1, 0, -1],
    [-1, -1, -1, 1, -1, -1, 1, -1, 0],
]


def edge_score(A, q, vertices=None):
    if vertices is None:
        vertices = range(len(A))
    vertices = list(vertices)
    return sum(
        A[u][v] * q[u][v]
        for pos, u in enumerate(vertices)
        for v in vertices[pos + 1 :]
    )


def matrix_norm(A):
    n = len(A)
    best = -1
    grounds = []
    for tail in itertools.product((-1, 1), repeat=n - 1):
        x = (1,) + tail
        h = sum(A[i][j] * x[i] * x[j] for i in range(n) for j in range(i + 1, n))
        if abs(h) > best:
            best, grounds = abs(h), []
        if abs(h) == best:
            grounds.append((1 if h > 0 else -1, x))
    return best, grounds


def principal(A, omitted):
    keep = [v for v in range(len(A)) if v not in set(omitted)]
    return [[A[u][v] for v in keep] for u in keep], keep


def child_extensions(A, omitted):
    n = len(A)
    keep = [v for v in range(n) if v != omitted]
    B, _ = principal(A, (omitted,))
    m, grounds = matrix_norm(B)
    out = []
    for t, local_x in grounds:
        x0 = {v: local_x[k] for k, v in enumerate(keep)}
        h = t * sum(A[omitted][j] * x0[j] for j in keep)
        missing_signs = (-1, 1) if h == 0 else ((1,) if h > 0 else (-1,))
        for xi in missing_signs:
            x = dict(x0)
            x[omitted] = xi
            q = [[0] * n for _ in range(n)]
            for u in range(n):
                for v in range(u + 1, n):
                    q[u][v] = q[v][u] = t * x[u] * x[v]
            score = edge_score(A, q)
            fields = [sum(A[j][k] * q[j][k] for k in range(n) if k != j) for j in range(n)]
            out.append({"t": t, "x": tuple(x[v] for v in range(n)), "q": q,
                        "score": score, "fields": fields})
    # Deduplicate a full oriented word that may have two projective descriptions.
    unique = {}
    for z in out:
        key = tuple(z["q"][u][v] for u in range(n) for v in range(u + 1, n))
        unique[key] = z
    return m, list(unique.values())


def build(A):
    n = len(A)
    M, _ = matrix_norm(A)
    rows = []
    child_M = []
    for i in range(n):
        mi, ext = child_extensions(A, i)
        child_M.append(mi)
        rows.append(ext)
    d = [M - z for z in child_M]
    for i in range(n):
        for z in rows[i]:
            z["s"] = M - z["score"]
            z["u"] = [z["s"] + z["fields"][j] - d[j] for j in range(n)]
            assert 0 <= z["s"] <= d[i]
            assert z["fields"][i] == d[i] - z["s"]
            assert z["u"][i] == 0 and min(z["u"]) >= 0
    return M, child_M, d, rows


def check_cut_residual(A, M, child_M, rows):
    n = len(A)
    for i in range(n):
        for zi in rows[i]:
            for j in range(n):
                if i == j:
                    continue
                # Any q_j ground gauges B_j.  The result is independent of its choice.
                zj = rows[j][0]
                support = []
                for a in range(n):
                    if a == j:
                        continue
                    for b in range(a + 1, n):
                        if b == j:
                            continue
                        if zi["q"][a][b] != zj["q"][a][b]:
                            support.append((a, b))
                residual = 2 * sum(A[a][b] * zj["q"][a][b] for a, b in support)
                direct = child_M[j] - edge_score(A, zi["q"], [v for v in range(n) if v != j])
                assert residual == direct == zi["u"][j]
                # The full ratio has augmented-cut form, hence cut/complement support.
                ratio = [[0] * n for _ in range(n)]
                for a in range(n):
                    for b in range(a + 1, n):
                        ratio[a][b] = ratio[b][a] = zi["q"][a][b] * zj["q"][a][b]
                orient = ratio[0][1] * ratio[0][2] * ratio[1][2]
                x = [1] + [ratio[0][v] * orient for v in range(1, n)]
                assert all(ratio[a][b] == orient * x[a] * x[b]
                           for a in range(n) for b in range(a + 1, n))


def check_assignments(rows, d):
    n = len(rows)
    chosen = [row[0] for row in rows]
    diagonal = sum(z["fields"][i] for i, z in enumerate(chosen))
    for perm in itertools.permutations(range(n)):
        lhs = sum(chosen[i]["u"][perm[i]] for i in range(n))
        rhs = sum(chosen[i]["fields"][perm[i]] for i in range(n)) - diagonal
        assert lhs == rhs >= 0


def terminal_data(A, deleted, q):
    C, keep = principal(A, deleted)
    terminal_M, _ = matrix_norm(C)
    carried = edge_score(A, q, keep)
    return terminal_M, carried, terminal_M - carried


def check_codimension_two(A, M, d, rows):
    n = len(A)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            for z in rows[i]:
                terminal_M, carried, refresh = terminal_data(A, (i, j), z["q"])
                edge = A[i][j] * z["q"][i][j]
                predicted_carried = M - d[i] - d[j] + z["s"] - z["u"][j] + edge
                assert carried == predicted_carried
                assert terminal_M == predicted_carried + refresh


def cycle_boundary(A, cycle, selected):
    total_u = sum(selected[a]["u"][b] for a, b in zip(cycle, cycle[1:] + cycle[:1]))
    boundary = 0
    for a, b in zip(cycle, cycle[1:] + cycle[:1]):
        qa, qb = selected[a]["q"], selected[b]["q"]
        boundary -= sum(A[b][v] * (qb[b][v] - qa[b][v]) for v in range(len(A)) if v != b)
    assert total_u == boundary
    assert 0 <= total_u <= 2 * len(cycle) * (len(A) - 1)
    return total_u


def find_zero_cycle(A, rows, length, terminal_minimum):
    n = len(A)
    for cycle in itertools.permutations(range(n), length):
        # Remove cyclic rotations from the search.
        if cycle[0] != min(cycle):
            continue
        candidates = []
        ok = True
        for a, b in zip(cycle, cycle[1:] + cycle[:1]):
            choices = [z for z in rows[a] if z["u"][b] == 0]
            if not choices:
                ok = False
                break
            candidates.append(choices)
        if not ok:
            continue
        selected = {a: choices[0] for a, choices in zip(cycle, candidates)}
        assert cycle_boundary(A, cycle, selected) == 0
        terminal_M, _, _ = terminal_data(A, cycle, selected[cycle[0]]["q"])
        if terminal_M > terminal_minimum:
            refreshes = [terminal_data(A, cycle, selected[a]["q"])[2] for a in cycle]
            return cycle, terminal_M, terminal_M - terminal_minimum, refreshes, selected
    raise AssertionError("no requested zero cycle")


def gauge_minimum(n):
    edges = [(i, j) for i in range(1, n) for j in range(i + 1, n)]
    ans = None
    for bits in itertools.product((-1, 1), repeat=len(edges)):
        A = [[0] * n for _ in range(n)]
        for j in range(1, n):
            A[0][j] = A[j][0] = 1
        for (i, j), bit in zip(edges, bits):
            A[i][j] = A[j][i] = bit
        value, _ = matrix_norm(A)
        ans = value if ans is None else min(ans, value)
    return ans


def find_refresh_wall(A, rows):
    n = len(A)
    for i in range(n):
        for j in range(i + 1, n):
            for zi in rows[i]:
                if zi["u"][j] != 0:
                    continue
                for zj in rows[j]:
                    if zj["u"][i] != 0:
                        continue
                    vi = terminal_data(A, (i, j), zi["q"])[2]
                    vj = terminal_data(A, (i, j), zj["q"])[2]
                    if vi > 0 and vj > 0:
                        return (i, j), vi, vj
    raise AssertionError("no refresh wall")


def main():
    known = {
        "A6": (A6, 5, {2: 1, 3: 3, 4: 4, 5: 4, 6: 5}),
        "A8": (A8, 10, {6: 5, 5: 4}),
        "A9": (A9, 12, {7: 9, 6: 5}),
    }
    built = {}
    for name, (A, expected_M, _) in known.items():
        M, child_M, d, rows = build(A)
        assert M == expected_M
        check_cut_residual(A, M, child_M, rows)
        check_codimension_two(A, M, d, rows)
        if name in ("A8", "A9"):
            check_assignments(rows, d)
        built[name] = (M, child_M, d, rows)
        print(name, "M", M, "child_M", child_M, "extensions", [len(r) for r in rows])

    # Independently certify the small global minima used below.
    assert gauge_minimum(4) == 4
    assert gauge_minimum(5) == 4
    assert gauge_minimum(6) == 5

    _, _, _, rows6 = built["A6"]
    print("A6 zero-2-cycle refresh wall", find_refresh_wall(A6, rows6))

    for name in ("A8", "A9"):
        A, _, minima = known[name]
        _, _, _, rows = built[name]
        z2 = find_zero_cycle(A, rows, 2, minima[len(A) - 2])
        z3 = find_zero_cycle(A, rows, 3, minima[len(A) - 3])
        print(name, "zero 2-cycle: cycle, terminal M, excess, refreshes", z2[:4])
        print(name, "zero 3-cycle: cycle, terminal M, excess, refreshes", z3[:4])


if __name__ == "__main__":
    main()
