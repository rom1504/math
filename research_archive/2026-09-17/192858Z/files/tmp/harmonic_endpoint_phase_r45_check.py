#!/usr/bin/env python3
"""Wave 45B checks for affine endpoint phases and the tropical-star wall."""

from __future__ import annotations

import math
import sys
from fractions import Fraction

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import A6, A8, A9
from harmonic_load_increments_r44_check import a8_tropical_wall, tropical_data
from harmonic_resonance_r34_check import Endpoint


def phase_certificate(A: np.ndarray, m: int) -> dict[str, object]:
    """Exact tropical affine/eta/degree certificate for a project endpoint."""
    H = Endpoint(A, 1.0, m)
    exponent, coefficient = tropical_data(H)
    energy = H.energies.astype(int)
    endpoint_line = energy - exponent
    active_mask = endpoint_line == int(np.max(endpoint_line))
    active = np.flatnonzero(active_mask)
    base_top = int(np.max(energy))
    ground = active[energy[active] == base_top]
    affine = len(ground) > 0
    z0 = sum((coefficient[d] for d in active), Fraction())
    zinf = sum((coefficient[d] for d in ground), Fraction())

    degree = np.zeros(H.N, dtype=int)
    load = [Fraction() for _ in range(H.N)]
    for bit in range(1, H.n):
        for x0, y0 in H.edges[bit]:
            x, y = int(x0), int(y0)
            ax, ay = bool(active_mask[x]), bool(active_mask[y])
            if ax and ay and energy[x] != energy[y]:
                high, low = (x, y) if energy[x] > energy[y] else (y, x)
                h = int(energy[high] - energy[low])
                kappa = Fraction(h) * coefficient[low] / (
                    coefficient[high] + coefficient[low]
                )
                for d in (x, y):
                    degree[d] += 1
                    load[d] += kappa
            elif ax != ay:
                xactive, other = (x, y) if ax else (y, x)
                if energy[other] > energy[xactive]:
                    h = int(energy[other] - energy[xactive])
                    degree[xactive] += 1
                    load[xactive] += h

    endpoint_mean = sum(
        (coefficient[d] * load[d] for d in active), Fraction()
    ) / z0
    delta = max((int(degree[d]) for d in active), default=0)
    return {
        "n": H.n,
        "m": m,
        "active": len(active),
        "ground": len(ground),
        "affine": affine,
        "eta": str(zinf / z0),
        "Delta": delta,
        "endpoint_mean": str(endpoint_mean),
    }


def envelope_facets(A: np.ndarray, m: int) -> dict[str, object]:
    """Find exact strict facets of max_d(E_d-s a_d) on [0,1]."""
    H = Endpoint(A, 1.0, m)
    exponent, _ = tropical_data(H)
    energy = H.energies.astype(int)
    points = {Fraction(), Fraction(1)}
    for i in range(H.N):
        for j in range(i):
            den = int(exponent[i] - exponent[j])
            if den:
                s = Fraction(int(energy[i] - energy[j]), den)
                if 0 < s < 1:
                    points.add(s)
    pts = sorted(points)
    facets: list[tuple[str, str, int, int]] = []
    for left, right in zip(pts, pts[1:]):
        mid = (left + right) / 2
        vals = [Fraction(int(energy[d])) - mid * int(exponent[d]) for d in range(H.N)]
        top = max(vals)
        winners = [d for d, value in enumerate(vals) if value == top]
        if winners:
            d = winners[0]
            facet = (str(left), str(right), int(energy[d]), int(exponent[d]))
            if not facets or facets[-1][2:] != facet[2:]:
                facets.append(facet)
            else:
                facets[-1] = (facets[-1][0], str(right), facet[2], facet[3])
    return {"m": m, "facets": facets}


def active_active_edge_integral(A: float, B: float, h: float, r: float, zinf: float) -> float:
    """Worst-case single-edge integral with Z replaced by its lower bound."""
    f = lambda u: h * h * B * B * math.exp(-(r + 2.0 * h) * u) / (
        zinf * (A + B * math.exp(-h * u))
    )
    return quad(f, 0.0, math.inf, epsabs=2e-12)[0]


def star_stats(k: int, h: float = 1.0) -> dict[str, float]:
    """Integrated variance and endpoint cost of the affine tropical star."""
    def variance(u: float) -> float:
        x = math.exp(-h * u)
        z = k + 2.0 + k * x
        p = x / (1.0 + x)
        mean = (k * h * p + k * x * h * p) / z
        second = ((k * h * p) ** 2 + k * x * (h * p) ** 2) / z
        return second - mean * mean

    integrated = quad(variance, 0.0, math.inf, epsabs=2e-11)[0]
    endpoint = k * h / (2.0 * (k + 1.0))
    eta = (k + 2.0) / (2.0 * k + 2.0)
    return {
        "k": float(k),
        "eta": eta,
        "I": integrated,
        "endpoint": endpoint,
        "ratio": integrated / endpoint,
        "ratio_over_k": integrated / endpoint / k,
        "theorem_rhs": k / eta * endpoint,
    }


def main() -> None:
    certs = []
    for name, A, sizes in (
        ("A6", A6, range(3, 7)),
        ("A8", A8, range(3, 9)),
        ("A9", A9, range(3, 10)),
    ):
        for m in sizes:
            z = phase_certificate(A, m)
            assert z["affine"]
            certs.append((name, z))
    a8 = phase_certificate(A8, 4)
    assert a8["eta"] == "4/5" and a8["Delta"] == 1
    wall = a8_tropical_wall()
    assert abs(float(wall["J2_over_endpoint_limit"]) - 0.0133226596003) < 2e-12
    print("project certificates", certs)
    print("A8 wall", a8, wall)

    generic = np.load("/home/math/quadra/tmp/r45_interior_matrix.npy")
    generic_cert = phase_certificate(generic, 4)
    generic_facets = envelope_facets(generic, 4)
    assert not generic_cert["affine"]
    assert len(generic_facets["facets"]) >= 2
    print("generic interior-facet witness", generic_cert, generic_facets)

    for pars in ((2.0, 3.0, 1.5, 0.0, 4.0), (7.0, 0.4, 2.0, 3.0, 1.2)):
        A, B, h, r, zinf = pars
        val = active_active_edge_integral(A, B, h, r, zinf)
        bound = h * B / zinf
        assert val <= bound * (1.0 + 2e-11)
        print("single-edge", pars, {"integral": val, "bound": bound})

    stars = [star_stats(k) for k in (10, 100, 1000, 10000)]
    for z in stars:
        assert z["I"] <= z["theorem_rhs"] * (1.0 + 1e-10)
    assert abs(stars[-1]["ratio_over_k"] - 0.25) < 0.002
    print("affine stars", stars)
    print("PASS harmonic_endpoint_phase_r45_check")


if __name__ == "__main__":
    main()
