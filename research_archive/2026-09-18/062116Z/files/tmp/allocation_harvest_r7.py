"""Exact small-order audit of two-channel allocation path harvesting.

All energies use the ledger's doubled normalization x^T A x.  For every
signed complete graph and every positive/negative endpoint pair, we form the
split and the Q-benchmarked payoffs mu_X^sigma from (10.407).

Two recursively defined quantities are audited:

  F_sum(U): largest obligation that can be routed down ONE nested child path
            if both orientation buckets on the selected edge may be added;
  F_one(U): same, but an actual successor orientation must be selected, so the
            selected edge contributes max_sigma (mu_X^sigma)_+.

At a selected child X, at most |P(X)-N(X)| may be passed onward.  Endpoint
pairs and all subsequent decompositions are optimized.  Thus F_one/R is a
generous test of any proposed single-tower purification: it allows the tower
to choose the best endpoint pair afresh at every node.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
from itertools import combinations


def edge_list(n: int):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


@lru_cache(None)
def spins(n: int):
    # Quotient the irrelevant global sign by fixing the last coordinate +1.
    if n == 0:
        return ((),)
    out = []
    for bits in range(1 << (n - 1)):
        out.append(tuple(-1 if bits >> i & 1 else 1 for i in range(n - 1)) + (1,))
    return tuple(out)


def matrix_from_code(n: int, code: int):
    a = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(edge_list(n)):
        a[i][j] = a[j][i] = 1 if code >> k & 1 else -1
    return tuple(tuple(row) for row in a)


def code_from_matrix(a):
    code = 0
    for k, (i, j) in enumerate(edge_list(len(a))):
        if a[i][j] == 1:
            code |= 1 << k
    return code


def induced(a, ids):
    return tuple(tuple(a[i][j] for j in ids) for i in ids)


def energy(a, x):
    return sum(a[i][j] * x[i] * x[j] for i in range(len(a)) for j in range(len(a)))


@lru_cache(None)
def endpoint_data(n: int, code: int):
    a = matrix_from_code(n, code)
    rows = tuple((energy(a, x), x) for x in spins(n))
    p = max((v for v, _ in rows), default=0)
    mn = min((v for v, _ in rows), default=0)
    ps = tuple(x for v, x in rows if v == p)
    ns = tuple(x for v, x in rows if v == mn)
    return p, -mn, ps, ns


def split_record(a, pspin, nspin):
    n = len(a)
    ratio = tuple(pspin[i] * nspin[i] for i in range(n))
    s = tuple(i for i in range(n) if ratio[i] == 1)
    t = tuple(i for i in range(n) if ratio[i] == -1)
    if not s or not t:
        return None
    out = []
    for x, y in ((s, t), (t, s)):
        ax = induced(a, x)
        px, nx, _, _ = endpoint_data(len(x), code_from_matrix(ax))
        qx = max(px, nx)
        hx = energy(ax, tuple(pspin[i] for i in x))
        # Fix p on retained X and expose every sign on peeled Y.
        fields = [sum(a[j][i] * pspin[i] for i in x) for j in y]
        lx = sum(abs(v) for v in fields)
        mu_plus = 2 * lx - (qx - hx)
        mu_minus = 2 * lx - (qx + hx)
        out.append(
            {
                "x": x,
                "y": y,
                "P": px,
                "N": nx,
                "Q": qx,
                "h": hx,
                "L": lx,
                "mu": (mu_plus, mu_minus),
                "d": abs(px - nx),
                "code": code_from_matrix(ax),
            }
        )
    return tuple(out)


@lru_cache(None)
def path_values(n: int, code: int):
    """Return optimal (F_sum, F_one) and exact witnesses."""
    p, nn, ps, ns = endpoint_data(n, code)
    r = p + nn
    if n <= 1 or r == 0:
        return 0, 0, None, None
    a = matrix_from_code(n, code)
    best_sum = (-1, None)
    best_one = (-1, None)
    for xp in ps:
        for xn in ns:
            rec = split_record(a, xp, xn)
            if rec is None:
                continue
            for shore_index, shore in enumerate(rec):
                child_sum, child_one, _, _ = path_values(len(shore["x"]), shore["code"])
                positive = tuple(max(0, z) for z in shore["mu"])
                here_sum = sum(positive)
                here_one = max(positive)
                value_sum = here_sum + min(shore["d"], child_sum)
                value_one = here_one + min(shore["d"], child_one)
                common = {
                    "p": xp,
                    "n": xn,
                    "split": tuple(z["x"] for z in rec),
                    "shore": shore_index,
                    "data": shore,
                    "edge_positive": positive,
                }
                if value_sum > best_sum[0]:
                    best_sum = (value_sum, common)
                if value_one > best_one[0]:
                    best_one = (value_one, common)
    if best_sum[0] < 0:
        raise AssertionError((n, code, p, nn))
    return best_sum[0], best_one[0], best_sum[1], best_one[1]


@lru_cache(None)
def liberal_one_path(n: int, code: int):
    """Best actual-orientation payoff sum on one path.

    This deliberately drops the imbalance-obligation restriction, so it is
    an upper bound on what any one-path purification of (10.408) could use.
    """
    p, nn, ps, ns = endpoint_data(n, code)
    if n <= 1 or p + nn == 0:
        return 0, ()
    a = matrix_from_code(n, code)
    best = (-1, None)
    for xp in ps:
        for xn in ns:
            rec = split_record(a, xp, xn)
            if rec is None:
                continue
            for shore_index, shore in enumerate(rec):
                child_value, child_trace = liberal_one_path(len(shore["x"]), shore["code"])
                positive = tuple(max(0, z) for z in shore["mu"])
                here = max(positive)
                value = here + child_value
                trace = (
                    {
                        "p": xp,
                        "n": xn,
                        "split": tuple(z["x"] for z in rec),
                        "shore": shore_index,
                        "data": shore,
                        "edge_positive": positive,
                        "chosen": here,
                    },
                ) + child_trace
                if value > best[0]:
                    best = (value, trace)
    return best


def fixed_pair_values(a, xp, xn):
    """One-step branch capacities, without optimizing the root endpoint pair."""
    rec = split_record(a, xp, xn)
    rows = []
    for shore in rec:
        child_sum, child_one, _, _ = path_values(len(shore["x"]), shore["code"])
        positive = tuple(max(0, z) for z in shore["mu"])
        rows.append(
            {
                **shore,
                "positive": positive,
                "F_sum_branch": sum(positive) + min(shore["d"], child_sum),
                "F_one_branch": max(positive) + min(shore["d"], child_one),
            }
        )
    return rows


def show_witness(n, code):
    a = matrix_from_code(n, code)
    p, nn, ps, ns = endpoint_data(n, code)
    fs, fo, ws, wo = path_values(n, code)
    print("matrix n/code", n, code)
    for row in a:
        print(list(row))
    print("P N R F_sum F_one ratios", p, nn, p + nn, fs, fo, fs / (p + nn), fo / (p + nn))
    print("best-sum-root", ws)
    print("best-one-root", wo)
    print("endpoint multiplicities", len(ps), len(ns))


def audit_all(max_n: int):
    for n in range(2, max_n + 1):
        worst_sum = None
        worst_one = None
        count = 1 << (n * (n - 1) // 2)
        for code in range(count):
            p, nn, _, _ = endpoint_data(n, code)
            r = p + nn
            fs, fo, _, _ = path_values(n, code)
            zsum = (fs / r, fs, r, code)
            zone = (fo / r, fo, r, code)
            if worst_sum is None or zsum < worst_sum:
                worst_sum = zsum
            if worst_one is None or zone < worst_one:
                worst_one = zone
        print("AUDIT", n, "matrices", count, "worst F_sum/R", worst_sum, "worst F_one/R", worst_one)
        show_witness(n, worst_one[3])


def ledger_examples():
    a6 = (
        (0, -1, 1, -1, -1, -1),
        (-1, 0, 1, -1, 1, 1),
        (1, 1, 0, -1, 1, -1),
        (-1, -1, -1, 0, 1, -1),
        (-1, 1, 1, 1, 0, -1),
        (-1, 1, -1, -1, -1, 0),
    )
    p6 = (-1, -1, -1, 1, -1, 1)
    n6 = (-1, -1, 1, -1, 1, 1)
    print("LEDGER n=6 fixed endpoint pair")
    print("P/N", endpoint_data(6, code_from_matrix(a6))[:2])
    for z in fixed_pair_values(a6, p6, n6):
        print(z)
    print("LEDGER n=6 optimized tree")
    show_witness(6, code_from_matrix(a6))

    a12 = (
        (0,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1),
        (-1,0,-1,-1,-1,-1,1,-1,1,-1,-1,-1),
        (-1,-1,0,1,-1,1,1,-1,1,1,1,-1),
        (-1,-1,1,0,-1,1,1,1,1,1,1,1),
        (-1,-1,-1,-1,0,1,1,1,1,-1,1,1),
        (-1,-1,1,1,1,0,1,-1,-1,-1,-1,-1),
        (-1,1,1,1,1,1,0,-1,1,-1,-1,-1),
        (1,-1,-1,1,1,-1,-1,0,-1,1,1,1),
        (-1,1,1,1,1,-1,1,-1,0,1,1,1),
        (1,-1,1,1,-1,-1,-1,1,1,0,1,-1),
        (-1,-1,1,1,1,-1,-1,1,1,1,0,-1),
        (-1,-1,-1,1,1,-1,-1,1,1,-1,-1,0),
    )
    print("LEDGER n=12 optimized tree")
    show_witness(12, code_from_matrix(a12))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", type=int, default=0)
    parser.add_argument("--ledger", action="store_true")
    args = parser.parse_args()
    if args.audit:
        audit_all(args.audit)
    if args.ledger:
        ledger_examples()


if __name__ == "__main__":
    main()
