"""Exact-integer physical checks and floating scalar replay of shared phases.

The b=8 toy uses exact zero/ferromagnetic slices; it tests the composition
identities, not the separate quantitative epsilon-band size condition.
"""
import itertools
import json
import math
from pathlib import Path

import numpy as np


def entropy(counts, denominator):
    values = np.asarray(counts, dtype=float) / denominator
    positive = values > 0
    return float(-np.sum(values[positive] * np.log(values[positive])))


def response_counts(words, counts, denominator):
    queries = np.asarray([
        [a] * 8 + [b] * 8 for a, b in itertools.product((-1, 1), repeat=2)
    ], dtype=np.int64)
    numerators = np.abs(words @ queries.T).T @ counts
    # sqrt(n)=4, so these are exact rational normalized responses.
    assert np.all(numerators == numerators[0])
    return {"numerator": int(numerators[0]), "denominator": int(4 * denominator),
            "value": float(numerators[0] / (4 * denominator))}


def check_physical():
    block = np.asarray(list(itertools.product((-1, 1), repeat=8)), dtype=np.int64)
    sums = np.sum(block, axis=1)
    cold = sums == 0
    hot = np.abs(sums) == 8
    q_counts = cold.astype(np.int64) + 5 * hot.astype(np.int64)
    assert int(q_counts.sum()) == 80
    # All 65,536 actual physical words at n=16.
    words = np.hstack((np.repeat(block, 256, axis=0), np.tile(block, (256, 1))))
    rho_counts = (np.outer(cold, cold).astype(np.int64)
                  + 175 * np.outer(hot, hot).astype(np.int64)).reshape(-1)
    eta_counts = np.outer(q_counts, q_counts).reshape(-1)
    rho_den, eta_den = 5600, 6400
    assert int(rho_counts.sum()) == rho_den
    assert int(eta_counts.sum()) == eta_den
    for counts, den in ((rho_counts, rho_den), (eta_counts, eta_den)):
        assert np.all(words.T @ counts == 0)
        gram = words.T @ (counts[:, None] * words)
        assert np.array_equal(gram, den * np.eye(16, dtype=np.int64))
        marg = counts.reshape(256, 256).sum(axis=1)
        assert np.array_equal(marg * 80, q_counts * den)

    phase_entropy = -(7 / 8) * math.log(7 / 8) - (1 / 8) * math.log(1 / 8)
    total_correlation = 2 * entropy(q_counts, 80) - entropy(rho_counts, rho_den)
    assert abs(total_correlation - phase_entropy) < 1e-12

    # Exact full-support repair delta=1/64, keeping identical repaired marginals.
    qp_counts = 63 * 256 * q_counts + 80
    qp_den = 64 * 80 * 256
    rp_counts = 63 * 65536 * rho_counts + 5600
    rp_den = 64 * 5600 * 65536
    ep_counts = np.outer(qp_counts, qp_counts).reshape(-1)
    ep_den = qp_den * qp_den
    for counts, den in ((rp_counts, rp_den), (ep_counts, ep_den)):
        assert np.min(counts) > 0
        assert int(counts.sum()) == den
        assert np.all(words.T @ counts == 0)
        gram = words.T @ (counts[:, None] * words)
        assert np.array_equal(gram, den * np.eye(16, dtype=np.int64))
        marg = counts.reshape(256, 256).sum(axis=1)
        assert np.array_equal(marg * qp_den, qp_counts * den)

    return {
        "physical_order": 16, "complete_words": 65536,
        "all_covariance_and_marginal_checks": "exact integer PASS",
        "shared_response": response_counts(words, rho_counts, rho_den),
        "product_response": response_counts(words, eta_counts, eta_den),
        "shared_fullsupport_response": response_counts(words, rp_counts, rp_den),
        "product_fullsupport_response": response_counts(words, ep_counts, ep_den),
        "total_correlation": total_correlation,
        "exact_formula_phase_entropy": phase_entropy,
    }


def scalar_sequence(r):
    length = 1
    while length <= 2 * r:
        length *= 2
    base = np.zeros(length)
    base[0], base[1], base[-1] = 7 / 8, 1 / 16, 1 / 16
    mass = np.fft.ifft(np.fft.fft(base) ** r).real
    assert float(np.min(mass)) > -1e-12
    mass = np.maximum(mass, 0)
    mass /= mass.sum()
    j = np.arange(length)
    abs_sum = np.minimum(j, length - j)
    product = math.sqrt(8 / r) * float(mass @ abs_sum)
    m = r // 2
    central = math.exp(math.lgamma(2 * m + 1) - 2 * math.lgamma(m + 1)
                       - 2 * m * math.log(2))
    mu_r = r * central
    shared = mu_r / math.sqrt(8 * r)
    return {"blocks": r, "product_response": product, "shared_response": shared,
            "product_limit": math.sqrt(2 / math.pi),
            "shared_limit": math.sqrt(2 / math.pi) / math.sqrt(8)}


def main():
    result = {"status": "PASS", "finite_physical": check_physical(),
              "scalar_fft_sequence": [scalar_sequence(r) for r in
                                      (1, 2, 4, 16, 64, 256, 1024, 4096)],
              "scope": "Exact integer finite physical identities; FFT limits are numerical replay."}
    path = Path("tmp/paper_portfolio_2026_09_17/bernoulli/shared_phase_composition.json")
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
