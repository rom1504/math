#!/usr/bin/env python3
"""Exact-data MILP for minimum path-cover/Carleson mass K.

The combinatorial data (endpoints, capacities, imbalances, and endpoint cuts)
are computed with Python integers.  The optimization is the finite MILP
described in ledger (10.407)--(10.412): choose one endpoint-pair action at
every active block, choose a conserved allocation of the full root range,
and minimize theta_root.  Once the endpoint actions are fixed, the problem is
an LP.  HiGHS solves both models; printed rational values are reconstructed
from the LP/MILP output and should be independently certified before being
called a theorem.

No file outside /home/math/quadra/tmp is written.
"""

from __future__ import annotations

import argparse
from collections import defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import isqrt

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


def popcount(x: int) -> int:
    # The project venv currently uses Python 3.9.
    return bin(x).count("1")


def frac(x: float, denominator: int = 1000000) -> Fraction:
    return Fraction(float(x)).limit_denominator(denominator)


@dataclass(frozen=True)
class Shore:
    child: int
    caps: tuple[int, int]
    imbalance: int
    P: int
    N: int
    h: int
    L: int


@dataclass(frozen=True)
class Action:
    shores: tuple[Shore, Shore]
    p: tuple[int, ...]
    n: tuple[int, ...]

    @property
    def key(self):
        return tuple(
            (s.child, s.caps, s.imbalance, s.P, s.N, s.h, s.L)
            for s in self.shores
        )


class EndpointSystem:
    def __init__(self, A):
        self.A = tuple(tuple(int(v) for v in row) for row in A)
        self.n = len(A)
        assert all(len(row) == self.n for row in self.A)
        assert all(self.A[i][i] == 0 for i in range(self.n))
        assert all(self.A[i][j] == self.A[j][i] for i in range(self.n) for j in range(self.n))
        self._endpoint_cache = {}
        self._action_cache = {}

    def ids(self, mask):
        return tuple(i for i in range(self.n) if mask >> i & 1)

    def energy(self, ids, spin):
        return sum(
            self.A[i][j] * spin[ii] * spin[jj]
            for ii, i in enumerate(ids)
            for jj, j in enumerate(ids)
        )

    def endpoints(self, mask):
        if mask in self._endpoint_cache:
            return self._endpoint_cache[mask]
        ids = self.ids(mask)
        m = len(ids)
        if m <= 1:
            ans = (0, 0, ((1,) if m else (),), ((1,) if m else (),))
            self._endpoint_cache[mask] = ans
            return ans
        # Quotient global sign by fixing the first local coordinate to +1.
        rows = []
        for tail in product((-1, 1), repeat=m - 1):
            x = (1,) + tail
            rows.append((self.energy(ids, x), x))
        P = max(e for e, _ in rows)
        mn = min(e for e, _ in rows)
        ps = tuple(x for e, x in rows if e == P)
        ns = tuple(x for e, x in rows if e == mn)
        ans = P, -mn, ps, ns
        self._endpoint_cache[mask] = ans
        return ans

    def actions(self, mask):
        if mask in self._action_cache:
            return self._action_cache[mask]
        ids = self.ids(mask)
        if len(ids) <= 1:
            self._action_cache[mask] = ()
            return ()
        P, N, ps, ns = self.endpoints(mask)
        unique = {}
        for p in ps:
            for n in ns:
                smask = sum(1 << ids[k] for k in range(len(ids)) if p[k] * n[k] == 1)
                tmask = mask ^ smask
                if not smask or not tmask:
                    continue
                shores = []
                for xmask, ymask in ((smask, tmask), (tmask, smask)):
                    xids = self.ids(xmask)
                    yids = self.ids(ymask)
                    pos = {i: k for k, i in enumerate(ids)}
                    px = tuple(p[pos[i]] for i in xids)
                    PX, NX, _, _ = self.endpoints(xmask)
                    QX = max(PX, NX)
                    h = self.energy(xids, px)
                    fields = [sum(self.A[j][i] * p[pos[i]] for i in xids) for j in yids]
                    L = sum(abs(v) for v in fields)
                    caps = (max(0, 2 * L - (QX - h)), max(0, 2 * L - (QX + h)))
                    shores.append(Shore(xmask, caps, abs(PX - NX), PX, NX, h, L))
                shores.sort(key=lambda z: z.child)
                act = Action(tuple(shores), p, n)
                # Endpoint labels with identical exact optimization data are equivalent.
                unique.setdefault(act.key, act)
        acts = tuple(unique.values())
        if not acts:
            raise AssertionError(("no endpoint split", mask, P, N, len(ps), len(ns)))
        self._action_cache[mask] = acts
        return acts

    def reachable(self):
        root = (1 << self.n) - 1
        seen = {root}
        q = deque([root])
        incoming = defaultdict(list)
        while q:
            U = q.popleft()
            if popcount(U) <= 1:
                continue
            for k, act in enumerate(self.actions(U)):
                for r, shore in enumerate(act.shores):
                    X = shore.child
                    incoming[X].append((U, k, r))
                    if X not in seen:
                        seen.add(X)
                        q.append(X)
        return tuple(sorted(seen, key=lambda x: (-popcount(x), x))), incoming


class Model:
    def __init__(self, sys: EndpointSystem, fixed=None, root_obligation=None):
        self.sys = sys
        self.root = (1 << sys.n) - 1
        self.subsets, self.incoming = sys.reachable()
        self.fixed = fixed or {}
        self.root_obligation = root_obligation
        self.names = []
        self.index = {}
        self.lb = []
        self.ub = []
        self.integrality = []
        self.rows = []
        self.row_lb = []
        self.row_ub = []
        self.z = {}
        self.w = {}
        self.theta = {}
        self.a = {}
        self.flow = {}
        self.q = {}
        self._build()

    def var(self, name, lb=0.0, ub=np.inf, integer=0):
        j = len(self.names)
        self.names.append(name)
        self.index[name] = j
        self.lb.append(lb)
        self.ub.append(ub)
        self.integrality.append(integer)
        return j

    def constraint(self, terms, lb=-np.inf, ub=np.inf):
        self.rows.append(dict(terms))
        self.row_lb.append(lb)
        self.row_ub.append(ub)

    def _build(self):
        S = self.sys
        # Node obligation and congestion variables.
        for U in self.subsets:
            P, N, _, _ = S.endpoints(U)
            upper_w = 0 if popcount(U) <= 1 else P + N
            if U == self.root:
                root_w = P + N if self.root_obligation is None else self.root_obligation
                if not (0 <= root_w <= P + N):
                    raise ValueError(("invalid root obligation", root_w, P + N))
                self.w[U] = self.var(("w", U), root_w, root_w)
            else:
                self.w[U] = self.var(("w", U), 0, upper_w)
            upper_theta = 4 * popcount(U)
            if popcount(U) <= 1:
                self.theta[U] = self.var(("theta", U), 0, 0)
            else:
                self.theta[U] = self.var(("theta", U), 0, upper_theta)

        # Action, allocation, outgoing obligation, and edge-max variables.
        for U in self.subsets:
            if popcount(U) <= 1:
                continue
            acts = S.actions(U)
            for k, act in enumerate(acts):
                if U in self.fixed:
                    lo = hi = 1 if k == self.fixed[U] else 0
                    integer = 0
                else:
                    lo, hi, integer = 0, 1, 1
                z = self.z[U, k] = self.var(("z", U, k), lo, hi, integer)
                for r, shore in enumerate(act.shores):
                    X = shore.child
                    for sig, cap in enumerate(shore.caps):
                        av = self.a[U, k, r, sig] = self.var(("a", U, k, r, sig), 0, cap)
                        self.constraint({av: 1, z: -cap}, ub=0)
                    fv = self.flow[U, k, r] = self.var(("flow", U, k, r), 0, shore.imbalance)
                    self.constraint({fv: 1, z: -shore.imbalance}, ub=0)
                    Mx = 4 * popcount(X)
                    qv = self.q[U, k, r] = self.var(("q", U, k, r), 0, Mx)
                    self.constraint({qv: 1, z: -Mx}, ub=0)
                    # q >= sum_sigma a_sigma/c_sigma.
                    load = {qv: -1}
                    for sig, cap in enumerate(shore.caps):
                        if cap:
                            load[self.a[U, k, r, sig]] = Fraction(1, cap)
                    self.constraint(load, ub=0)
                    # If this action is selected, q >= theta_child.
                    # theta_X - q + Mx*z <= Mx.
                    self.constraint({self.theta[X]: 1, qv: -1, z: Mx}, ub=Mx)

        # One selected action per active node, with activation generated by its parent.
        for U in self.subsets:
            incoming_z = defaultdict(int)
            for V, k, r in self.incoming.get(U, ()):
                incoming_z[self.z[V, k]] += 1
            if U == self.root:
                terms = {self.z[U, k]: 1 for k in range(len(S.actions(U)))}
                self.constraint(terms, lb=1, ub=1)
            elif popcount(U) > 1:
                terms = {self.z[U, k]: 1 for k in range(len(S.actions(U)))}
                for j, coef in incoming_z.items():
                    terms[j] = terms.get(j, 0) - coef
                self.constraint(terms, lb=0, ub=0)
                # A subset cannot be activated by two alternative parents.
                self.constraint(incoming_z, ub=1)
            else:
                self.constraint(incoming_z, ub=1)

            # w_U equals the sum of incoming obligation flows (except at root).
            if U != self.root:
                terms = {self.w[U]: 1}
                for V, k, r in self.incoming.get(U, ()):
                    terms[self.flow[V, k, r]] = terms.get(self.flow[V, k, r], 0) - 1
                self.constraint(terms, lb=0, ub=0)

            if popcount(U) > 1:
                # Local conservation.
                terms = {self.w[U]: 1}
                for k, act in enumerate(S.actions(U)):
                    for r, shore in enumerate(act.shores):
                        terms[self.flow[U, k, r]] = -1
                        for sig in range(2):
                            terms[self.a[U, k, r, sig]] = -1
                self.constraint(terms, lb=0, ub=0)
                # theta is the sum of the two edge maxima for the selected action.
                terms = {self.theta[U]: 1}
                for k, act in enumerate(S.actions(U)):
                    for r in range(2):
                        terms[self.q[U, k, r]] = -1
                self.constraint(terms, lb=0, ub=0)
                # Inactive nodes have theta zero.
                terms = {self.theta[U]: 1}
                for k in range(len(S.actions(U))):
                    terms[self.z[U, k]] = terms.get(self.z[U, k], 0) - 4 * popcount(U)
                self.constraint(terms, ub=0)

    def solve(self, time_limit=300, mip_gap=1e-10, presolve=False):
        rr, cc, vv = [], [], []
        for i, row in enumerate(self.rows):
            for j, v in row.items():
                rr.append(i)
                cc.append(j)
                vv.append(float(v))
        mat = coo_matrix((vv, (rr, cc)), shape=(len(self.rows), len(self.names))).tocsr()
        c = np.zeros(len(self.names))
        c[self.theta[self.root]] = 1
        options = {"time_limit": time_limit, "mip_rel_gap": mip_gap, "presolve": presolve}
        res = milp(
            c,
            integrality=np.array(self.integrality),
            bounds=Bounds(np.array(self.lb), np.array(self.ub)),
            constraints=LinearConstraint(mat, np.array(self.row_lb), np.array(self.row_ub)),
            options=options,
        )
        return res

    def selected_actions(self, x, threshold=0.5):
        return {U: k for (U, k), j in self.z.items() if x[j] > threshold}

    def describe(self, res):
        print("status", res.status, res.message)
        print("nodes/actions/vars/constraints", len(self.subsets), len(self.z), len(self.names), len(self.rows))
        print("objective", res.fun, "~", frac(res.fun))
        if hasattr(res, "mip_gap"):
            print("mip_gap", res.mip_gap, "mip_node_count", res.mip_node_count, "dual_bound", res.mip_dual_bound)
        selected = self.selected_actions(res.x)
        for U in sorted(selected, key=lambda q: (-popcount(q), q)):
            k = selected[U]
            w = res.x[self.w[U]]
            th = res.x[self.theta[U]]
            if w < 1e-8 and U != self.root:
                continue
            print("block", self.sys.ids(U), "w", frac(w), "theta", frac(th), "action", k)
            act = self.sys.actions(U)[k]
            for r, shore in enumerate(act.shores):
                av = tuple(frac(res.x[self.a[U, k, r, s]]) for s in range(2))
                flow = frac(res.x[self.flow[U, k, r]])
                qv = frac(res.x[self.q[U, k, r]])
                print("  ->", self.sys.ids(shore.child), "caps", shore.caps, "d", shore.imbalance,
                      "a", av, "pass", flow, "q", qv)
        return selected


def matrix_10409():
    return (
        (0,-1,1,-1,-1,-1),(-1,0,1,-1,1,1),(1,1,0,-1,1,-1),
        (-1,-1,-1,0,1,-1),(-1,1,1,1,0,-1),(-1,1,-1,-1,-1,0),
    )


def matrix_10420():
    return (
        (0,1,-1,1,-1,-1),(1,0,1,-1,-1,-1),(-1,1,0,-1,-1,-1),
        (1,-1,-1,0,-1,-1),(-1,-1,-1,-1,0,-1),(-1,-1,-1,-1,-1,0),
    )


def strict12():
    return (
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


def family_10416(m):
    k = isqrt(m)
    if k * k != m:
        raise ValueError("m must be a square")
    s = (1,) * ((m + k) // 2) + (-1,) * ((m - k) // 2)
    C = tuple(tuple(0 if i == j else s[i] * s[j] for j in range(m)) for i in range(m))
    A = []
    for i in range(2 * m):
        row = []
        for j in range(2 * m):
            if i < m and j < m:
                v = C[i][j]
            elif i >= m and j >= m:
                v = -C[i-m][j-m]
            else:
                v = 1
            row.append(v)
        A.append(tuple(row))
    return tuple(A)


def positive_clique(n):
    return tuple(tuple(0 if i == j else 1 for j in range(n)) for i in range(n))


class ConcavePL:
    """Increasing concave piecewise-linear function, constant after last knot."""

    def __init__(self, knots):
        clean = []
        for x, y in knots:
            x, y = Fraction(x), Fraction(y)
            if clean and x == clean[-1][0]:
                clean[-1] = (x, max(y, clean[-1][1]))
            elif not clean or (x > clean[-1][0] and y >= clean[-1][1]):
                clean.append((x, y))
            else:
                raise ValueError((clean[-1] if clean else None, (x, y)))
        self.knots = tuple(clean)

    def __call__(self, x):
        x = Fraction(x)
        if x <= self.knots[0][0]:
            return self.knots[0][1]
        for (x0, y0), (x1, y1) in zip(self.knots, self.knots[1:]):
            if x <= x1:
                return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
        return self.knots[-1][1]

    def inverse(self, target):
        target = Fraction(target)
        if target <= self.knots[0][1]:
            return self.knots[0][0]
        for (x0, y0), (x1, y1) in zip(self.knots, self.knots[1:]):
            if target <= y1:
                if y1 == y0:
                    return x0
                return x0 + (x1 - x0) * (target - y0) / (y1 - y0)
        if target == self.knots[-1][1]:
            return self.knots[-1][0]
        raise ValueError(("target exceeds range", target, self.knots[-1]))


def local_clique_service(m, q):
    """Maximum local allocation on one K_m endpoint shore at load q."""
    q = Fraction(q)
    c = Fraction(m * m, 2)
    if q <= 1:
        return c * q
    if q <= 2:
        return c + m * (q - 1)
    return c + m


_clique_service_cache = {}


def clique_service(m):
    """H_m(K): maximum K_m obligation serviceable with congestion at most K.

    This uses only the exact balanced endpoint geometry of the positive
    clique.  It is also valid for every switching/sign reversal of K_m.
    """
    if m in _clique_service_cache:
        return _clique_service_cache[m]
    if m == 1:
        ans = ConcavePL(((0, 0),))
    elif m == 2:
        ans = ConcavePL(((0, 0), (2, 4)))
    elif m % 2:
        raise ValueError("structured clique recursion currently requires a power of two")
    else:
        h = m // 2
        child = clique_service(h)
        d = max(0, h * (h - 2))
        qd = child.inverse(d) if d else Fraction(0)
        qmax = max(Fraction(2), qd)
        qs = {Fraction(0), Fraction(1), Fraction(2), qd, qmax}
        qs.update(x for x, _ in child.knots if x <= qd)
        qs = sorted(q for q in qs if 0 <= q <= qmax)

        def raw(q):
            return 2 * (local_clique_service(m, q) + min(Fraction(d), child(q)))

        R = Fraction(m * m)
        out = []
        for q0, q1 in zip(qs, qs[1:]):
            y0, y1 = raw(q0), raw(q1)
            K0, K1 = 2 * q0, 2 * q1
            if not out:
                out.append((K0, min(R, y0)))
            if y0 < R < y1:
                qstar = q0 + (q1 - q0) * (R - y0) / (y1 - y0)
                out.append((2 * qstar, R))
                break
            out.append((K1, min(R, y1)))
            if y1 >= R:
                break
        if len(qs) == 1:
            out = [(0, min(R, raw(0)))]
        ans = ConcavePL(out)
    _clique_service_cache[m] = ans
    return ans


def structured_family_k(m):
    """Exact K_min for (10.416), conditional only on its forced root cut.

    At the root the endpoints are unique, and both children are switched
    copies of K_m.  Endpoint choices inside a clique are all isomorphic, so
    the scalar recursion is the full optimization, not a heuristic policy.
    """
    if m < 2 or m & (m - 1):
        raise ValueError("m must be a power of two")
    child = clique_service(m)
    d = m * (m - 2)
    c = m * m + m

    def edge_service(q):
        q = Fraction(q)
        local = c * min(q, Fraction(2))
        return local + min(Fraction(d), child(q))

    target_edge = Fraction(2 * m * m)
    # Assemble every possible slope-change of edge_service and invert exactly.
    qd = child.inverse(d) if d else Fraction(0)
    qs = sorted({Fraction(0), Fraction(2), qd, *(x for x, _ in child.knots if x <= qd)})
    qstar = None
    for q0, q1 in zip(qs, qs[1:]):
        y0, y1 = edge_service(q0), edge_service(q1)
        if y0 <= target_edge <= y1:
            qstar = q0 if y1 == y0 else q0 + (q1 - q0) * (target_edge - y0) / (y1 - y0)
            break
    if qstar is None:
        if edge_service(qs[-1]) == target_edge:
            qstar = qs[-1]
        else:
            raise AssertionError((m, target_edge, [(q, edge_service(q)) for q in qs]))
    return 2 * qstar, child


def run(name, A, time_limit, fixed=None, root_obligation=None):
    print("===", name, "n", len(A), "===")
    sys = EndpointSystem(A)
    root = (1 << len(A)) - 1
    P, N, ps, ns = sys.endpoints(root)
    print("P N R endpoint multiplicities root actions", P, N, P + N, len(ps), len(ns), len(sys.actions(root)))
    model = Model(sys, fixed=fixed, root_obligation=root_obligation)
    res = model.solve(time_limit=time_limit)
    if res.x is not None:
        chosen = model.describe(res)
        # Re-solve the chosen endpoint tree as a continuous LP consistency check.
        fixed_model = Model(sys, fixed=chosen, root_obligation=root_obligation)
        fixed_res = fixed_model.solve(time_limit=time_limit)
        print("fixed-tree LP", fixed_res.status, fixed_res.fun, "~", frac(fixed_res.fun))
    return sys, model, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", choices=("10409", "10420", "strict12", "family", "clique", "all"), default="all")
    ap.add_argument("--m", type=int, default=4)
    ap.add_argument("--time", type=float, default=300)
    ap.add_argument("--obligation", type=float)
    args = ap.parse_args()
    cases = []
    if args.case in ("10409", "all"):
        cases.append(("10.409", matrix_10409()))
    if args.case in ("10420", "all"):
        cases.append(("10.420", matrix_10420()))
    if args.case in ("strict12", "all"):
        cases.append(("strict reset n=12", strict12()))
    if args.case == "family":
        cases.append((f"10.416 m={args.m}", family_10416(args.m)))
    if args.case == "clique":
        cases.append((f"K_{args.m}", positive_clique(args.m)))
    for name, A in cases:
        run(name, A, args.time, root_obligation=args.obligation)


if __name__ == "__main__":
    main()
