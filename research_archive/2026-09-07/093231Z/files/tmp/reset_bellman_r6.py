"""Exact reset-path dynamic program for the quadratic-signing descent.

All energies use the ledger's doubled normalization x^T A x.  By switching
the fixed inherited spin to 1, a state is (principal mask, inherited
orientation rho).  A transition chooses an absolute endpoint (tau, x), and
the next mask is the nonempty proper agreement shore {i: x_i=1}.
"""

from __future__ import annotations

import argparse
import itertools
import random
from dataclasses import dataclass
from functools import lru_cache


@dataclass(frozen=True)
class Metric:
    p: int
    n: int
    h: int
    # (endpoint orientation tau, agreement mask)
    transitions: tuple[tuple[int, int], ...]

    @property
    def q(self) -> int:
        return max(self.p, self.n)

    @property
    def r(self) -> int:
        return self.p + self.n


def edge_list(n: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


class Instance:
    def __init__(self, n: int, signs: tuple[int, ...]):
        self.n = n
        self.edges = edge_list(n)
        assert len(signs) == len(self.edges)
        self.signs = signs

    @lru_cache(None)
    def energy(self, mask: int, plus: int) -> int:
        """Energy of + on plus and - on mask\plus."""
        ans = 0
        for a, (i, j) in zip(self.signs, self.edges):
            bi, bj = 1 << i, 1 << j
            if mask & bi and mask & bj:
                ans += 2 * a * (1 if bool(plus & bi) == bool(plus & bj) else -1)
        return ans

    @lru_cache(None)
    def metric(self, mask: int) -> Metric:
        # Enumerating all plus masks includes the global-sign duplicate, which
        # is useful because each representative exposes the other shore.
        energies = [(self.energy(mask, plus), plus) for plus in submasks(mask)]
        p = max(e for e, _ in energies)
        n = -min(e for e, _ in energies)
        q = max(p, n)
        h = self.energy(mask, mask)
        tr = []
        for e, plus in energies:
            if plus == 0 or plus == mask:
                continue
            if e == q:
                tr.append((1, plus))
            if e == -q:
                tr.append((-1, plus))
        return Metric(p, n, h, tuple(sorted(set(tr))))

    @lru_cache(None)
    def max_reset(self, mask: int, rho: int):
        """Maximum cumulative ordinary reset deficit; return cost, path."""
        m = self.metric(mask)
        best = (0, ())
        for tau, child in m.transitions:
            reset = tau == -rho
            # a=q-tau*H=q+rho*H at an opposite-orientation step.
            a = m.q - tau * m.h if reset else 0
            tail, path = self.max_reset(child, tau)
            row = (a + tail, ((mask, rho, tau, child, a),) + path)
            if row[0] > best[0]:
                best = row
        return best

    @lru_cache(None)
    def max_stat(self, mask: int, rho: int, stat: str):
        """Maximum cumulative selected statistic over reset transitions."""
        m = self.metric(mask)
        best = (0, ())
        for tau, child in m.transitions:
            is_reset = tau == -rho
            charge = 0
            if is_reset:
                old_cap = m.p if rho == 1 else m.n
                h_old = rho * m.h
                vals = {
                    "a": m.q + h_old,
                    "b": m.q - old_cap,
                    "u": old_cap - h_old,
                    "g": m.q - h_old,
                    "q": m.q,
                    "one": 1,
                }
                charge = vals[stat]
            tail, path = self.max_stat(child, tau, stat)
            row = (charge + tail, ((mask, rho, tau, child, charge),) + path)
            if row[0] > best[0]:
                best = row
        return best

    @lru_cache(None)
    def max_stat_strict(self, mask: int, rho: int, stat: str):
        """As max_stat, but never change orientation across a P=N tie."""
        m = self.metric(mask)
        best = (0, ())
        old_cap = m.p if rho == 1 else m.n
        opposite_cap = m.n if rho == 1 else m.p
        for tau, child in m.transitions:
            is_reset = tau == -rho
            if is_reset and opposite_cap <= old_cap:
                continue
            charge = 0
            if is_reset:
                h_old = rho * m.h
                vals = {
                    "a": m.q + h_old,
                    "b": m.q - old_cap,
                    "u": old_cap - h_old,
                    "g": m.q - h_old,
                    "q": m.q,
                    "one": 1,
                }
                charge = vals[stat]
            tail, path = self.max_stat_strict(child, tau, stat)
            row = (charge + tail, ((mask, rho, tau, child, charge),) + path)
            if row[0] > best[0]:
                best = row
        return best

    @lru_cache(None)
    def max_reset_with_count(self, mask: int, rho: int, count: int):
        """Maximum reset cost among paths with exactly count resets."""
        if count == 0:
            return 0, ()
        best = None
        m = self.metric(mask)
        for tau, child in m.transitions:
            is_reset = int(tau == -rho)
            if is_reset > count:
                continue
            tail = self.max_reset_with_count(child, tau, count - is_reset)
            if tail is None:
                continue
            a = m.q - tau * m.h if is_reset else 0
            row = (a + tail[0], ((mask, rho, tau, child, a),) + tail[1])
            if best is None or row[0] > best[0]:
                best = row
        return best

    def matrix_rows(self) -> list[list[int]]:
        rows = [[0] * self.n for _ in range(self.n)]
        for a, (i, j) in zip(self.signs, self.edges):
            rows[i][j] = rows[j][i] = a
        return rows

    def describe_path(self, path) -> list[dict]:
        rows = []
        for mask, rho, tau, child, a in path:
            m, mc = self.metric(mask), self.metric(child)
            rows.append(
                dict(
                    mask=mask,
                    vertices=bits(mask),
                    rho=rho,
                    tau=tau,
                    reset=tau == -rho,
                    child=child,
                    child_vertices=bits(child),
                    a=a,
                    P=m.p,
                    N=m.n,
                    Q=m.q,
                    H=m.h,
                    R=m.r,
                    childP=mc.p,
                    childN=mc.n,
                    childQ=mc.q,
                    childH=mc.h,
                    childR=mc.r,
                )
            )
        return rows


def submasks(mask: int):
    sub = mask
    while True:
        yield sub
        if sub == 0:
            break
        sub = (sub - 1) & mask


def bits(mask: int) -> list[int]:
    return [i for i in range(mask.bit_length()) if mask & (1 << i)]


def signing_from_code(n: int, code: int) -> tuple[int, ...]:
    return tuple(1 if code & (1 << k) else -1 for k in range(n * (n - 1) // 2))


def print_certificate(tag: str, ratio: float, inst: Instance, rho: int, cost: int, path):
    full = (1 << inst.n) - 1
    m = inst.metric(full)
    print(tag, "ratio", ratio, "cost", cost, "rho", rho,
          "P,N,Q,H,R", (m.p, m.n, m.q, m.h, m.r))
    print("matrix", inst.matrix_rows())
    for row in inst.describe_path(path):
        print(" ", row)


def exhaustive(n: int):
    full = (1 << n) - 1
    max_by_den = {"R": None, "Q": None}
    max_resets = None
    local_range_failure = None
    local_q_failure = None
    total = 1 << (n * (n - 1) // 2)
    for code in range(total):
        inst = Instance(n, signing_from_code(n, code))
        for rho in (1, -1):
            cost, path = inst.max_reset(full, rho)
            m = inst.metric(full)
            for key, den in (("R", m.r), ("Q", m.q)):
                ratio = cost / den if den else 0.0
                old = max_by_den[key]
                if old is None or ratio > old[0]:
                    max_by_den[key] = (ratio, inst, rho, cost, path)
            resets = sum(a > 0 for *_, a in path)
            if max_resets is None or resets > max_resets[0]:
                max_resets = (resets, inst, rho, cost, path)
            for row in inst.describe_path(path):
                if not row["reset"]:
                    continue
                dr = row["R"] - row["childR"]
                dq = row["Q"] - row["childQ"]
                if row["a"] > dr:
                    gap = row["a"] - dr
                    if local_range_failure is None or gap > local_range_failure[0]:
                        local_range_failure = (gap, inst, rho, cost, path, row)
                if row["a"] > 2 * dq:
                    gap = row["a"] - 2 * dq
                    if local_q_failure is None or gap > local_q_failure[0]:
                        local_q_failure = (gap, inst, rho, cost, path, row)
    print("exhaustive n", n, "signings", total)
    for key, row in max_by_den.items():
        print_certificate("max cost/" + key, row[0], row[1], row[2], row[3], row[4])
    print("max positive resets", max_resets[0])
    print_certificate("max-reset-count", max_resets[0], max_resets[1], max_resets[2], max_resets[3], max_resets[4])
    print("worst local a-(R-Rchild)", None if local_range_failure is None else local_range_failure[0])
    if local_range_failure:
        print_certificate("local-range-failure", local_range_failure[0], local_range_failure[1], local_range_failure[2], local_range_failure[3], local_range_failure[4])
    print("worst local a-2(Q-Qchild)", None if local_q_failure is None else local_q_failure[0])
    if local_q_failure:
        print_certificate("local-q-failure", local_q_failure[0], local_q_failure[1], local_q_failure[2], local_q_failure[3], local_q_failure[4])


def random_scan(n: int, trials: int, seed: int):
    rng = random.Random(seed)
    full = (1 << n) - 1
    max_row = None
    for _ in range(trials):
        signs = tuple(rng.choice((-1, 1)) for _ in edge_list(n))
        inst = Instance(n, signs)
        for rho in (1, -1):
            cost, path = inst.max_reset(full, rho)
            m = inst.metric(full)
            ratio = cost / m.r if m.r else 0.0
            if max_row is None or ratio > max_row[0]:
                max_row = (ratio, inst, rho, cost, path)
    print("random n,trials", n, trials)
    print_certificate("max cost/R", max_row[0], max_row[1], max_row[2], max_row[3], max_row[4])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustive", type=int)
    ap.add_argument("--random", type=int)
    ap.add_argument("--trials", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=6102026)
    args = ap.parse_args()
    if args.exhaustive:
        exhaustive(args.exhaustive)
    if args.random:
        random_scan(args.random, args.trials, args.seed)


if __name__ == "__main__":
    main()
