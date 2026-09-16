#!/usr/bin/env python3
"""Exact checks for the Wave 31 harmonic-composition identities."""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np
from scipy.special import logsumexp

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import A4, A6, oriented_cuts, restriction_key


def flip_vertex(d: np.ndarray, i: int) -> np.ndarray:
    out = d.copy()
    out[i, :] *= -1
    out[:, i] *= -1
    out[i, i] = 0
    return out


class HarmonicData:
    def __init__(self, A: np.ndarray, beta: float):
        self.A = np.asarray(A, dtype=int)
        self.beta = float(beta)
        self.n = len(A)
        self.cuts = oriented_cuts(self.n)
        self.energy = np.einsum("bij,ij->b", self.cuts, self.A).astype(float)
        self.logZA = float(logsumexp(self.beta * self.energy))
        self.nu = np.exp(self.beta * self.energy - self.logZA)
        self.cut_index = {
            restriction_key(d, tuple(range(self.n))): j
            for j, d in enumerate(self.cuts)
        }
        self.F: dict[tuple[int, ...], dict[tuple[int, ...], float]] = {}
        self.groups: dict[tuple[int, ...], dict[tuple[int, ...], list[int]]] = {}
        # The matrix representation of an oriented cut has 2^m states only
        # for m >= 3; on <=2 vertices it forgets the global orientation bit.
        for m in range(3, self.n + 1):
            for S in itertools.combinations(range(self.n), m):
                groups: dict[tuple[int, ...], list[int]] = {}
                for j, d in enumerate(self.cuts):
                    groups.setdefault(restriction_key(d, S), []).append(j)
                assert len(groups) == 2**m
                assert {len(v) for v in groups.values()} == {2 ** (self.n - m)}
                vals: dict[tuple[int, ...], float] = {}
                for key, inds in groups.items():
                    child = np.asarray(
                        [np.sum(self.A[np.ix_(S, S)] * self.cuts[j][np.ix_(S, S)])
                         for j in inds], dtype=float
                    )
                    assert np.ptp(child) < 1e-12
                    external = self.energy[inds] - child[0]
                    assert abs(float(np.mean(external))) < 1e-12
                    vals[key] = float(logsumexp(self.beta * external)
                                      - (self.n - m) * math.log(2))
                self.F[S] = vals
                self.groups[S] = groups

    def f(self, S: tuple[int, ...], d: np.ndarray) -> float:
        return self.F[S][restriction_key(d, S)]

    def recurrence_checks(self) -> None:
        V = tuple(range(self.n))
        assert max(abs(v) for v in self.F[V].values()) < 1e-12
        for m in range(3, self.n):
            for S in itertools.combinations(range(self.n), m):
                for i in set(range(self.n)) - set(S):
                    R = tuple(sorted(S + (i,)))
                    for key, inds in self.groups[S].items():
                        dplus = self.cuts[inds[0]]
                        dminus = flip_vertex(dplus, i)
                        a = float(np.sum(self.A[np.ix_(R, R)] * dplus[np.ix_(R, R)])
                                  - np.sum(self.A[np.ix_(S, S)] * dplus[np.ix_(S, S)]))
                        assert abs(a - 2 * sum(self.A[i, j] * dplus[i, j] for j in S)) < 1e-12
                        fp, fm = self.f(R, dplus), self.f(R, dminus)
                        lhs = self.F[S][key]
                        direct = float(np.logaddexp(self.beta * a + fp,
                                                    -self.beta * a + fm) - math.log(2))
                        b = self.beta * a + (fp - fm) / 2
                        centered = (fp + fm) / 2 + math.log(math.cosh(b))
                        assert abs(lhs - direct) < 2e-11
                        assert abs(lhs - centered) < 2e-11

    def reverse_kl_checks(self) -> None:
        for S, groups in self.groups.items():
            k = self.n - len(S)
            for key, inds in groups.items():
                child = float(np.sum(self.A[np.ix_(S, S)]
                                     * self.cuts[inds[0]][np.ix_(S, S)]))
                logw = self.beta * (self.energy[inds] - child)
                logp = logw - logsumexp(logw)
                reverse_kl = -k * math.log(2) - float(np.mean(logp))
                assert abs(reverse_kl - self.F[S][key]) < 2e-11

    def block_chain_checks(self) -> None:
        # F_S = D(U_{R\S} || parent marginal on R\S | y)
        #       + E_{U_{R\S}} F_R.
        for m in range(3, self.n + 1):
            for S in itertools.combinations(range(self.n), m):
                outside = [i for i in range(self.n) if i not in S]
                for r in range(m, self.n + 1):
                    for add in itertools.combinations(outside, r - m):
                        R = tuple(sorted(S + add))
                        q = 2 ** (r - m)
                        for keyS, inds in self.groups[S].items():
                            byR: dict[tuple[int, ...], list[int]] = {}
                            for j in inds:
                                byR.setdefault(restriction_key(self.cuts[j], R), []).append(j)
                            assert len(byR) == q
                            logmass = []
                            fR = []
                            for keyR, js in byR.items():
                                logmass.append(float(logsumexp(self.beta * self.energy[js])))
                                fR.append(self.F[R][keyR])
                            logmass = np.asarray(logmass)
                            logp = logmass - logsumexp(logmass)
                            d_rev = -math.log(q) - float(np.mean(logp))
                            rhs = d_rev + float(np.mean(fR))
                            assert abs(rhs - self.F[S][keyS]) < 3e-11

    def endpoint_stats(self, m: int):
        selectors = list(itertools.combinations(range(self.n), m))
        Fall = np.asarray([[self.f(S, d) for d in self.cuts] for S in selectors])
        U = np.mean(np.exp(-Fall), axis=0)
        Phi = -np.log(U)
        posterior = np.exp(-Fall) / np.sum(np.exp(-Fall), axis=0)[None, :]
        variational = np.sum(posterior * Fall, axis=0) + np.sum(
            posterior * np.log(posterior * len(selectors)), axis=0
        )
        assert np.allclose(Phi, variational, rtol=2e-12, atol=2e-12)

        # Exact soft endpoint response under a vertex flip.  The cost is
        # unchanged for selectors not containing the flipped vertex.
        k = self.n - m
        for j, d in enumerate(self.cuts):
            for i in range(self.n):
                df = flip_vertex(d, i)
                jf = self.cut_index[restriction_key(df, tuple(range(self.n)))]
                delta = Fall[:, jf] - Fall[:, j]
                for si, S in enumerate(selectors):
                    if i not in S:
                        assert abs(delta[si]) < 2e-12
                    else:
                        assert abs(delta[si]) <= 4 * self.beta * k + 2e-11
                        if m >= 4:
                            S0 = tuple(x for x in S if x != i)
                            a = float(np.sum(self.A[np.ix_(S, S)] * d[np.ix_(S, S)])
                                      - np.sum(self.A[np.ix_(S0, S0)] * d[np.ix_(S0, S0)]))
                            residual = (Fall[si, j] - Fall[si, jf]) / 2
                            b = self.beta * a + residual
                            assert abs(delta[si] + 2 * residual) < 2e-12
                            # b is the renormalized parent-cavity half log-odds.
                            assert math.isfinite(b)
                response = -math.log(float(posterior[:, j] @ np.exp(-delta)))
                assert abs(response - (Phi[jf] - Phi[j])) < 3e-11

        f = U / float(self.nu @ U)
        ent = float((self.nu * f) @ np.log(f))
        meanF = float(np.mean(Fall @ self.nu))
        freegap = -math.log(float(self.nu @ U))
        varPhi = float(self.nu @ (Phi - self.nu @ Phi) ** 2)
        assert ent <= freegap + 2e-12
        assert freegap <= meanF + 2e-12
        return ent, meanF, freegap, varPhi, float(np.ptp(Phi))


def run(A: np.ndarray, beta: float) -> None:
    H = HarmonicData(A, beta)
    H.recurrence_checks()
    H.reverse_kl_checks()
    H.block_chain_checks()
    print(f"n={H.n} beta={beta}")
    print("m entropy mean_reverse_KL free_gap var_Phi osc_Phi")
    for m in range(3, H.n + 1):
        print(m, *H.endpoint_stats(m))


def main() -> None:
    run(A4, 0.5)
    run(A6, 0.5)
    print("PASS harmonic_composition_r31")


if __name__ == "__main__":
    main()
