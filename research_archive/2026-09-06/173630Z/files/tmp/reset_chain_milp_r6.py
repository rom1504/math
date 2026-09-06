"""MILP search for long singleton-peel alternating reset chains.

The inherited spin is switched to all ones.  On suffix R_t={t,...,n-1},
the prescribed endpoint flips only t and has orientation tau_t alternating
with t.  Endpoint domination is imposed against every Boolean state on every
suffix.  The objective maximizes reset cost minus C times the exact full
energy range P+N.
"""

from __future__ import annotations

import argparse
import itertools

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


def solve(
    n: int,
    constant: float,
    tau0: int = 1,
    time_limit: float = 300,
    strict: bool = False,
    block_sizes: tuple[int, ...] | None = None,
):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    edge_index = {e: k for k, e in enumerate(edges)}
    ne = len(edges)
    # Binary b_e represents signing e_e=2b_e-1; final variables are U,L.
    nv = ne + 2
    U, L = ne, ne + 1

    rows: list[dict[int, float]] = []
    lows: list[float] = []
    highs: list[float] = []

    def affine_energy(vertices, signs):
        """Return coeff dict and constant for doubled Boolean energy."""
        coeff = {}
        const = 0.0
        for ai, i in enumerate(vertices):
            for aj in range(ai + 1, len(vertices)):
                j = vertices[aj]
                z = signs[ai] * signs[aj]
                k = edge_index[(i, j)]
                # 2 e_ij z = 4 b_ij z - 2z.
                coeff[k] = 4.0 * z
                const -= 2.0 * z
        return coeff, const

    def add_ge(lhs_coeff, lhs_const, rhs_coeff, rhs_const):
        # lhs-rhs >= 0.
        row = dict(lhs_coeff)
        for k, v in rhs_coeff.items():
            row[k] = row.get(k, 0.0) - v
        rows.append(row)
        lows.append(rhs_const - lhs_const)
        highs.append(np.inf)

    if block_sizes is None:
        transitions = n - 2 if strict else n - 1
        block_sizes = (1,) * transitions + (n - transitions,)
    assert sum(block_sizes) == n and all(d > 0 for d in block_sizes)
    transitions = len(block_sizes) - 1
    starts = [0]
    for d in block_sizes:
        starts.append(starts[-1] + d)

    # Prescribed absolute endpoint on every suffix of blocks.
    for t in range(transitions):
        start, cut = starts[t], starts[t + 1]
        vertices = list(range(start, n))
        tau = tau0 * (-1 if t % 2 else 1)
        x = [-1] * (cut - start) + [1] * (n - cut)
        ex, cx = affine_energy(vertices, x)
        qcoeff = {k: tau * v for k, v in ex.items()}
        qconst = tau * cx
        # Representatives modulo global negation.
        for head in itertools.product((-1, 1), repeat=len(vertices) - 1):
            z = list(head) + [1]
            ez, cz = affine_energy(vertices, z)
            add_ge(qcoeff, qconst, ez, cz)
            add_ge(qcoeff, qconst, {k: -v for k, v in ez.items()}, -cz)
            if strict:
                # The opposite one-sided extremum is smaller.  All doubled
                # energies are congruent mod 4, so strictness means a gap >=4:
                # q >= (-tau) E(z)+4.
                opp = {k: -tau * v for k, v in ez.items()}
                add_ge(qcoeff, qconst, opp, -tau * cz + 4.0)

    # U and L are exact at optimum because their difference is penalized.
    vertices = list(range(n))
    for head in itertools.product((-1, 1), repeat=n - 1):
        z = list(head) + [1]
        ez, cz = affine_energy(vertices, z)
        # U >= E(z).
        add_ge({U: 1.0}, 0.0, ez, cz)
        # E(z) >= L.
        add_ge(ez, cz, {L: 1.0}, 0.0)

    rr, cc, vv = [], [], []
    for i, row in enumerate(rows):
        for j, v in row.items():
            rr.append(i)
            cc.append(j)
            vv.append(v)
    matrix = coo_matrix((vv, (rr, cc)), shape=(len(rows), nv)).tocsr()

    # scipy minimizes.  S=-4 sum_t tau_t sum_{j>t} e_tj.
    # In b variables, -S has coefficient +8 tau_t and constant omitted.
    objective = np.zeros(nv)
    s_const = 0.0
    for t in range(transitions):
        tau = tau0 * (-1 if t % 2 else 1)
        for i in range(starts[t], starts[t + 1]):
            for j in range(starts[t + 1], n):
                k = edge_index[(i, j)]
                objective[k] += 8.0 * tau  # coefficient of -S
                s_const += 4.0 * tau       # S constant after e=2b-1
    objective[U] = constant
    objective[L] = -constant

    max_energy = n * (n - 1)
    lower = np.r_[np.zeros(ne), -max_energy, -max_energy]
    upper = np.r_[np.ones(ne), max_energy, max_energy]
    result = milp(
        objective,
        integrality=np.r_[np.ones(ne), 0, 0],
        bounds=Bounds(lower, upper),
        constraints=LinearConstraint(matrix, np.array(lows), np.array(highs)),
        options={"time_limit": time_limit},
    )
    if result.x is None:
        return result, None
    signing = np.zeros((n, n), dtype=int)
    for k, (i, j) in enumerate(edges):
        signing[i, j] = signing[j, i] = 1 if result.x[k] > 0.5 else -1
    # Recompute exact certificate quantities.
    vals = []
    for head in itertools.product((-1, 1), repeat=n - 1):
        z = np.array(list(head) + [1], dtype=int)
        vals.append(int(z @ signing @ z))
    p, neg = max(vals), -min(vals)
    costs = []
    for t in range(transitions):
        tau = tau0 * (-1 if t % 2 else 1)
        cross = int(
            signing[starts[t] : starts[t + 1], starts[t + 1] :].sum()
        )
        costs.append(-4 * tau * cross)
    cert = {
        "matrix": signing.tolist(),
        "P": p,
        "N": neg,
        "R": p + neg,
        "costs": costs,
        "S": sum(costs),
        "ratio": sum(costs) / (p + neg),
        "objective_S_minus_C_R": sum(costs) - constant * (p + neg),
        "block_sizes": block_sizes,
    }
    return result, cert


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("--constant", type=float, default=2.0)
    ap.add_argument("--tau0", type=int, choices=(-1, 1), default=1)
    ap.add_argument("--time-limit", type=float, default=300)
    ap.add_argument("--strict", action="store_true")
    ap.add_argument("--blocks", help="comma-separated ordered block sizes")
    args = ap.parse_args()
    blocks = None if args.blocks is None else tuple(int(x) for x in args.blocks.split(","))
    result, cert = solve(
        args.n, args.constant, args.tau0, args.time_limit, args.strict, blocks
    )
    print("status", result.status, result.message, "mip_gap", getattr(result, "mip_gap", None))
    print(cert)


if __name__ == "__main__":
    main()
