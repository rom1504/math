#!/usr/bin/env python3
"""Wave 26 audit of the endpoint parent-cut likelihood ratio.

Everything is enumerated on the n-bit oriented-cut cube.  The cube
coordinates are the global orientation and switches of vertices 1,...,n-1.
No claim about an asymptotic signing family is inferred from these finite
calculations.
"""

from __future__ import annotations

from itertools import product
from fractions import Fraction
import math

import numpy as np
from scipy.special import logsumexp

from endpoint_transport_r25 import A4, A6, A8, A9, oriented_cuts, restriction_key


class Endpoint:
    def __init__(self, a_mat: np.ndarray):
        self.a_mat = a_mat
        self.n = len(a_mat)
        self.m = self.n - 1
        self.selectors = [tuple(j for j in range(self.n) if j != i) for i in range(self.n)]
        self.cuts = oriented_cuts(self.n)
        self.nd = len(self.cuts)
        self.ns = len(self.selectors)
        self.energy = np.einsum("bij,ij->b", self.cuts, a_mat, optimize=True).astype(float)
        self.c = np.empty((self.ns, self.nd), dtype=float)
        self.groups: list[list[list[int]]] = []
        for si, ids in enumerate(self.selectors):
            self.c[si] = np.asarray([
                np.sum(a_mat[np.ix_(ids, ids)] * d[np.ix_(ids, ids)]) for d in self.cuts
            ])
            by_key: dict[tuple[int, ...], list[int]] = {}
            for di, d in enumerate(self.cuts):
                by_key.setdefault(restriction_key(d, ids), []).append(di)
            assert len(by_key) == 2**self.m
            self.groups.append(list(by_key.values()))

        # Exact n-cube flip maps: global orientation, then vertex switches.
        index = {tuple(d.ravel()): i for i, d in enumerate(self.cuts)}
        maps = []
        maps.append(np.asarray([index[tuple((-d).ravel())] for d in self.cuts], dtype=int))
        for v in range(1, self.n):
            move = []
            for d in self.cuts:
                e = d.copy()
                e[v, :] *= -1
                e[:, v] *= -1
                e[v, v] = 0
                move.append(index[tuple(e.ravel())])
            maps.append(np.asarray(move, dtype=int))
        self.flip_maps = maps
        for flip in maps:
            assert np.array_equal(flip[flip], np.arange(self.nd))

    def evaluate(self, beta: float, gamma: float) -> dict[str, object]:
        log_za = float(logsumexp(beta * self.energy))
        nu = np.exp(beta * self.energy - log_za)
        ell = np.empty((self.ns, self.nd), dtype=float)
        log_zs = np.empty(self.ns, dtype=float)
        for si, groups in enumerate(self.groups):
            child_c = np.asarray([self.c[si, ds[0]] for ds in groups])
            log_zs[si] = float(logsumexp(gamma * child_c))
            for ds in groups:
                log_k = float(logsumexp(beta * (self.energy[ds] - self.c[si, ds])))
                ell[si, ds] = log_k + (beta - gamma) * self.c[si, ds]

        log_z0 = float(logsumexp(log_zs) - math.log(self.ns))
        log_c = log_za - log_z0
        log_l = log_c - ell                         # joint density dM0/d(U x nu)
        log_f = logsumexp(log_l, axis=0) - math.log(self.ns)  # E[L | D]
        f = np.exp(log_f)
        l = np.exp(log_l)
        w = np.exp(log_zs - log_z0)                 # E[L | S]
        pi = w / self.ns
        mu = nu * f

        assert abs(float(np.sum(mu)) - 1.0) < 3e-10
        assert np.allclose(np.mean(l, axis=0), f, rtol=2e-11, atol=2e-11)
        assert np.allclose(np.dot(l, nu), w, rtol=2e-11, atol=2e-11)

        ent_cut = float(np.dot(mu, log_f))
        total = float(np.sum((nu[None, :] / self.ns) * l * log_l))
        rho = l / (self.ns * f[None, :])
        selector_cond = float(np.sum(mu[None, :] * np.where(rho > 0, rho * np.log(rho * self.ns), 0.0)))
        d_pi = float(np.dot(pi, np.log(w)))
        mutual_info = selector_cond - d_pi

        directed_selectors = [(s, t) for s in range(self.ns) for t in range(self.ns) if s != t]
        selector_h_mean = float(np.mean([
            0.5 * np.dot(nu, (np.sqrt(l[s]) - np.sqrt(l[t])) ** 2)
            for s, t in directed_selectors
        ]))

        # q_S/nu=L_S/w_S.  Convexity/data processing says Ent(f) is at
        # most the pi-weighted channel KL, with exact gap I(S;D).
        r = l / w[:, None]
        avg_channel_kl = float(np.sum(pi[:, None] * nu[None, :] * r * np.log(r)))
        assert abs(total - (ent_cut + selector_cond)) < 2e-9
        assert abs(total - (d_pi + avg_channel_kl)) < 2e-9
        assert abs(avg_channel_kl - (ent_cut + mutual_info)) < 2e-9

        # Sum of one-coordinate conditional entropies for the Gibbs heat-bath
        # chain, plus its Jeffreys, Hellinger, and Metropolis Dirichlet forms.
        cond_entropy = 0.0
        heatbath_jeff = 0.0
        heatbath_hell = 0.0
        metropolis_jeff = 0.0
        coordinate_cond = []
        coordinate_jeff = []
        for flip in self.flip_maps:
            ce = 0.0
            hj = 0.0
            hh = 0.0
            mj = 0.0
            for i in range(self.nd):
                j = int(flip[i])
                if i >= j:
                    continue
                mass = nu[i] + nu[j]
                af = (nu[i] * f[i] + nu[j] * f[j]) / mass
                ce += nu[i] * f[i] * math.log(f[i] / af) + nu[j] * f[j] * math.log(f[j] / af)
                conductance = nu[i] * nu[j] / mass
                jeff = (f[i] - f[j]) * (log_f[i] - log_f[j])
                hj += conductance * jeff
                hh += conductance * (math.sqrt(f[i]) - math.sqrt(f[j])) ** 2
                mj += min(nu[i], nu[j]) * jeff
            cond_entropy += ce
            heatbath_jeff += hj
            heatbath_hell += hh
            metropolis_jeff += mj
            coordinate_cond.append(ce)
            coordinate_jeff.append(hj)

        # Ground graph for diagnosing the beta -> infinity bottleneck.
        qmax = float(np.max(self.energy))
        ground = set(int(i) for i in np.where(self.energy == qmax)[0])
        ground_edges = []
        for k, flip in enumerate(self.flip_maps):
            for i in ground:
                j = int(flip[i])
                if j in ground and i < j:
                    ground_edges.append((i, j, k))
        ground_adj = {i: set() for i in ground}
        for i, j, _ in ground_edges:
            ground_adj[i].add(j)
            ground_adj[j].add(i)
        remaining = set(ground)
        ground_components = []
        while remaining:
            stack = [remaining.pop()]
            component = set(stack)
            while stack:
                i = stack.pop()
                for j in ground_adj[i]:
                    if j not in component:
                        component.add(j)
                        remaining.discard(j)
                        stack.append(j)
            ground_components.append(component)

        assert cond_entropy <= heatbath_jeff + 2e-9

        return {
            "ent_cut": ent_cut,
            "total": total,
            "selector_cond": selector_cond,
            "d_pi": d_pi,
            "mutual_info": mutual_info,
            "avg_channel_kl": avg_channel_kl,
            "selector_h_mean": selector_h_mean,
            "cond_entropy": cond_entropy,
            "heatbath_jeff": heatbath_jeff,
            "heatbath_hell": heatbath_hell,
            "metropolis_jeff": metropolis_jeff,
            "coordinate_cond": coordinate_cond,
            "coordinate_jeff": coordinate_jeff,
            "f": f,
            "nu": nu,
            "mu": mu,
            "energy": self.energy,
            "qmax": qmax,
            "ground": ground,
            "ground_edges": ground_edges,
            "ground_components": ground_components,
        }

    def zero_gamma_limiting_energy(self) -> Fraction:
        """Mean best extension energy for uniform S and uniform child cut."""
        total = 0
        count = 0
        for groups in self.groups:
            for ds in groups:
                total += int(np.max(self.energy[ds]))
                count += 1
        return Fraction(total, count)


def abstract_two_mode(eps: float, delta: float = 0.5) -> dict[str, float]:
    """Positive endpoint-shaped two-bit Gibbs bottleneck.

    nu puts mass (1-2 eps)/2 on 00 and 11 and eps on 01 and 10.  It is an
    exact ferromagnetic two-spin Gibbs law.  L_1 and L_2 are normalized,
    positive, and each depends on only one retained coordinate.
    """
    assert 0 < eps < 0.25 and 0 < delta < 1
    heavy = (1.0 - 2.0 * eps) / 2.0
    nu = np.asarray([heavy, eps, eps, heavy])  # 00, 01, 10, 11
    hi, lo = 2.0 - delta, delta
    likelihood = np.asarray([
        [hi, hi, lo, lo],  # retained coordinate 1
        [hi, lo, hi, lo],  # retained coordinate 2
    ])
    assert np.all(likelihood > 0)
    assert np.allclose(likelihood @ nu, np.ones(2))
    f = np.mean(likelihood, axis=0)
    assert abs(float(np.dot(nu, f)) - 1.0) < 1e-14
    ent = float(np.dot(nu * f, np.log(f)))
    selector = float(np.sum((nu[None, :] / 2.0) * likelihood * np.log(likelihood / f[None, :])))
    selector_h = 0.5 * float(np.sum(nu * (np.sqrt(likelihood[0]) - np.sqrt(likelihood[1])) ** 2))

    flips = (np.asarray([2, 3, 0, 1]), np.asarray([1, 0, 3, 2]))
    cond_entropy = 0.0
    hb_jeff = 0.0
    hb_hell = 0.0
    for flip in flips:
        for i in range(4):
            j = int(flip[i])
            if i >= j:
                continue
            mass = nu[i] + nu[j]
            af = (nu[i] * f[i] + nu[j] * f[j]) / mass
            cond_entropy += nu[i] * f[i] * math.log(f[i] / af) + nu[j] * f[j] * math.log(f[j] / af)
            conductance = nu[i] * nu[j] / mass
            hb_jeff += conductance * (f[i] - f[j]) * (math.log(f[i]) - math.log(f[j]))
            hb_hell += conductance * (math.sqrt(f[i]) - math.sqrt(f[j])) ** 2

    k_delta = hi * math.log(hi) + lo * math.log(lo)
    assert abs(ent - heavy * k_delta) < 2e-14
    assert abs(selector - eps * k_delta) < 2e-14
    assert abs(selector_h - eps * (math.sqrt(hi) - math.sqrt(lo)) ** 2) < 2e-14
    hb_exact = 2.0 * eps * (1.0 - 2.0 * eps) * (1.0 - delta) * math.log(hi / lo)
    assert abs(hb_jeff - hb_exact) < 2e-14
    return {
        "ent_cut": ent,
        "selector_entropy": selector,
        "selector_h": selector_h,
        "cond_entropy": cond_entropy,
        "heatbath_jeff": hb_jeff,
        "heatbath_hell": hb_hell,
    }


def ratio(x: float, y: float) -> float:
    return x / y if y > 0 else math.inf


def main() -> None:
    endpoints = {name: Endpoint(a) for name, a in (("A4", A4), ("A6", A6), ("A8", A8), ("A9", A9))}
    expected_q = {"A4": 8, "A6": 10, "A8": 20, "A9": 24}
    for name, ep in endpoints.items():
        assert int(np.max(ep.energy)) == expected_q[name]

    print("Likelihood/data-processing identities and matched low-temperature Gibbs-flip audit")
    print("name beta Entcut sum_cond HB_Jeff Metro_Jeff Ent/sum_cond Ent/HB_Jeff")
    records: dict[str, list[tuple[float, dict[str, object]]]] = {}
    for name, ep in endpoints.items():
        rows = []
        for beta in (0.1, 0.5, 1.0, 2.0, 4.0, 8.0):
            z = ep.evaluate(beta, beta)
            rows.append((beta, z))
            print(
                f"{name} {beta:4.1f} {z['ent_cut']:.10g} {z['cond_entropy']:.10g} "
                f"{z['heatbath_jeff']:.10g} {z['metropolis_jeff']:.10g} "
                f"{ratio(z['ent_cut'], z['cond_entropy']):.8g} "
                f"{ratio(z['ent_cut'], z['heatbath_jeff']):.8g}"
            )
        records[name] = rows

    # Reproduce the A4 coefficient audit and likelihood identities at the
    # intentionally mismatched temperatures from Wave 25.
    a4 = endpoints["A4"].evaluate(0.5, 0.1)
    assert abs(a4["ent_cut"] - 0.6276970227309644) < 2e-11
    print("A4 beta=.5 gamma=.1:", {k: a4[k] for k in (
        "ent_cut", "selector_cond", "d_pi", "mutual_info", "avg_channel_kl",
        "cond_entropy", "heatbath_jeff", "heatbath_hell", "metropolis_jeff"
    )})

    # A9 has a nonconstant endpoint likelihood on its ground face at matched
    # low temperature.  Its ground graph has multiple cube-flip components;
    # record rather than assume the resulting bottleneck.
    ep9 = endpoints["A9"]
    z9 = ep9.evaluate(8.0, 8.0)
    print("A9 ground count/edges:", len(z9["ground"]), len(z9["ground_edges"]))
    print("A9 ground component sizes:", sorted((len(x) for x in z9["ground_components"]), reverse=True))
    assert sorted((len(x) for x in z9["ground_components"]), reverse=True) == [15, 8, 2]
    ground_f = z9["f"][sorted(z9["ground"])]
    print("A9 beta=8 ground f range:", float(np.min(ground_f)), float(np.max(ground_f)))

    # Exact-minimizer signing endpoints disprove temperature-uniform
    # absorption into selector entropy/Hellinger.  At gamma=0, the limiting
    # cut law chooses a uniform child cut and then a best parent extension.
    print("gamma=0 low-temperature slopes")
    for name, ep in endpoints.items():
        ebar = ep.zero_gamma_limiting_energy()
        deficit = Fraction(int(max(ep.energy)), 1) - ebar
        z4 = ep.evaluate(4.0, 0.0)
        z8 = ep.evaluate(8.0, 0.0)
        slope = (z8["ent_cut"] - z4["ent_cut"]) / 4.0
        assert abs(slope - float(deficit)) < 4e-4
        assert z8["selector_cond"] <= math.log(ep.ns) + 1e-10
        assert z8["selector_h_mean"] <= 1.0 + 1e-10
        print(
            name, "Ebar", ebar, "deficit", deficit, "observed slope", slope,
            "Ent/selector", ratio(z8["ent_cut"], z8["selector_cond"]),
            "mean H", z8["selector_h_mean"]
        )

    # Exact positive two-mode obstruction to any generic dimension-free
    # parent approximate-tensorization or Gibbs-flip inequality based only on
    # the endpoint conditional-expectation structure.
    print("abstract positive two-mode obstruction, delta=1/2")
    previous = None
    for eps in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5):
        z = abstract_two_mode(eps)
        required = ratio(z["ent_cut"], z["heatbath_jeff"])
        print(eps, z, "Ent/HB_Jeff", required, "Ent/sum_cond", ratio(z["ent_cut"], z["cond_entropy"]))
        if previous is not None:
            assert required > 8.0 * previous
        previous = required

    print("all parent-cut identities and finite audits verified")


if __name__ == "__main__":
    main()
