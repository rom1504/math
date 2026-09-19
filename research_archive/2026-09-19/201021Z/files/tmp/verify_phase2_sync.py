from itertools import product


def ultra(vals, eta):
    # vals are the three edges xy, xz, yz of one triangle.
    a, b, c = vals
    return (
        c + eta >= min(a, b)
        and b + eta >= min(a, c)
        and a + eta >= min(b, c)
    )


def cancellation(r1, r2, i, j):
    q_i = (r1[i] + r2[i]) / 2
    q_j = (r1[j] + r2[j]) / 2
    return (
        (abs(r1[i] - r1[j]) + abs(r2[i] - r2[j])) / 2
        - abs(q_i - q_j)
    )


grid = (0.0, 0.5, 1.0)
exact_count = 0
for r1 in product(grid, repeat=3):
    if not ultra(r1, 0.0):
        continue
    for r2 in product(grid, repeat=3):
        if not ultra(r2, 0.0):
            continue
        if not ultra(tuple(a + b for a, b in zip(r1, r2)), 0.0):
            continue
        exact_count += 1
        q = tuple((a + b) / 2 for a, b in zip(r1, r2))
        # On three states all pair labels are adjacent. Exact no-crossing
        # therefore says q-order must be coordinatewise order.
        for i in range(3):
            for j in range(3):
                if q[i] <= q[j]:
                    assert r1[i] <= r1[j] and r2[i] <= r2[j]
                if q[i] == q[j]:
                    assert cancellation(r1, r2, i, j) == 0

print({"exact_profiles_checked": exact_count})
