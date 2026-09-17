"""Replay exact signed-covariance identities and numerical covariance diagnostics.

The Gaussian-sign CLT is proved analytically in the linked artifact;
this script does not promote finite floating output to a CLT proof.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np


def check(a: np.ndarray, words: np.ndarray, weights: np.ndarray) -> dict:
    n = len(a)
    denominator = int(weights.sum())
    energy = np.einsum("bi,ij,bj->b", words, a, words)//2
    mark = np.where(energy >= 0, 1, -1)
    sigma_num = words.T @ (weights[:, None]*words)
    signed_num = words.T @ ((weights*mark)[:, None]*words)
    mean_mark_num = int(weights @ mark)
    hollow_num = signed_num-mean_mark_num*np.eye(n, dtype=np.int64)
    mean_energy_num = int(weights @ np.abs(energy))
    assert np.trace(hollow_num) == 0
    assert int(np.sum(a*hollow_num)) == 2*mean_energy_num
    frob_num = int(np.sum(hollow_num*hollow_num))
    assert frob_num*n*(n-1) >= 4*mean_energy_num**2
    # Exact trace comparison, before any eigenvalue calculation.
    assert frob_num <= int(np.sum(sigma_num*sigma_num))
    sigma = sigma_num/denominator
    k0 = hollow_num/denominator
    k = signed_num/denominator
    l = float(np.linalg.eigvalsh(sigma)[-1])
    t = 1/(2*(l+1))
    identity = np.eye(n)
    for polarity in (-1, 1):
        eig = np.linalg.eigvalsh(identity+polarity*t*k0)
        assert eig[0] >= .5-1e-12 and eig[-1] <= 1.5+1e-12
    kappa2 = 2/np.pi
    b = kappa2*np.arcsin(t*k0)
    # Gaussian orthant covariance formula is exact analytically.
    for polarity in (-1, 1):
        assert np.linalg.eigvalsh(identity+polarity*b)[0] >= -1e-12
    v = np.einsum("bi,ij,bj->b", words, b, words)/n
    signed_v = float((weights*mark) @ v/denominator)
    trace_v = float(np.sum(k*b)/n)
    lower = kappa2*t*frob_num/(denominator**2*n)
    assert abs(signed_v-trace_v) < 1e-11
    assert signed_v >= lower-1e-12
    assert np.max(np.abs(v)) <= 1+1e-12
    f = lambda z: (np.sqrt(1+z)+np.sqrt(1-z))/2
    average_f = float(weights @ f(v)/denominator)
    assert average_f <= f(lower)+1e-12
    assert average_f <= 1-lower**2/8+1e-12
    return {"n": n, "support": len(words), "covariance_norm": l,
            "average_energy_normalized": mean_energy_num/(denominator*n**1.5),
            "signed_variance_statistic": signed_v,
            "guaranteed_statistic": lower,
            "gaussian_mixture_mean_factor": average_f}


def main() -> None:
    rng = np.random.default_rng(20260917)
    reports = []
    for n in (4, 6, 8, 10, 12):
        words = np.asarray([(1,)+w for w in itertools.product((-1, 1), repeat=n-1)], dtype=np.int64)
        for _ in range(12):
            raw = rng.choice((-1, 1), size=(n,n))
            a = np.triu(raw, 1)
            a += a.T
            energy = np.einsum("bi,ij,bj->b", words, a, words)//2
            ground = words[np.abs(energy) == np.abs(energy).max()]
            for weighted in (False, True):
                weights = rng.integers(1, 8, size=len(ground), dtype=np.int64) if weighted else np.ones(len(ground), dtype=np.int64)
                reports.append(check(a, ground, weights))
    # An exactly isotropic high-energy query law, in both eigensectors.
    h4 = np.ones((4,4), dtype=np.int64)-2*np.eye(4, dtype=np.int64)
    walsh4 = np.asarray([[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]], dtype=np.int64)
    h, basis = h4, walsh4
    for exponent in (1, 2, 3):
        if exponent > 1:
            h = np.kron(h,h4)
            basis = np.kron(basis,walsh4)
        a = h-np.diag(np.diag(h))
        report = check(a, basis, np.ones(len(basis), dtype=np.int64))
        assert abs(report["covariance_norm"]-1) < 1e-12
        reports.append(report)
    result = {"status": "PASS", "exact_signed_identity_cases": len(reports),
              "float_diagnostics": "PSD, arcsine trace, and concavity; no finite CLT claim",
              "isotropic_high_energy_cases": reports[-3:],
              "largest_signed_variance_statistic": max(r["signed_variance_statistic"] for r in reports)}
    target = Path(__file__).resolve().parents[1]/"tmp/paper_portfolio_2026_09_17/bernoulli/signed_covariance_response.json"
    target.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
