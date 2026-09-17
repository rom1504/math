"""Exact finite flip/anchor and martingale replay, not a concentration proof."""

from itertools import product
from pathlib import Path
import json
import math

import numpy as np


def main():
    rng = np.random.default_rng(2026091721)
    counts = {"flip_worlds": 0, "escape_implications": 0,
              "martingale_nodes": 0, "gaussian_half_deficit": 0,
              "clone_parameter_checks": 0, "binomial_mode_checks": 0}
    for n in range(3, 6):
        words = np.array(list(product((-1, 1), repeat=n)), dtype=float)
        ii, jj = np.triu_indices(n, 1)
        queries = words[:, ii] * words[:, jj]
        queries = np.concatenate((queries, -queries))
        bits = np.array(list(product((0, 1), repeat=len(ii))), dtype=float)
        for trial in range(4):
            signing = rng.choice((-1, 1), size=len(ii))
            energy = queries @ signing
            cap = energy.max()
            ground = int(energy.argmax())
            deficit = cap-energy
            increments = queries-queries[ground]
            positive = np.unique(deficit[deficit > 0])
            threshold = float(positive[min(trial, len(positive)-1)])
            outside = deficit >= threshold
            rho = (0.07, 0.13, 0.21, 0.24)[trial]
            s = 1-2*rho
            centered = -2*(bits-rho)*signing
            actual = (1-2*bits)*signing
            assert np.all(np.isin(actual, (-1, 1)))
            assert np.max(np.abs(actual-s*signing-centered)) < 1e-12
            new_energy = actual @ queries.T
            f = np.maximum(0, np.max(centered @ increments[outside].T
                                    -s*deficit[outside]/2, axis=1))
            best = new_energy.max(axis=1)
            escape = np.any((new_energy[:, outside] >= best[:, None]-1e-12), axis=1)
            assert np.all(f[escape] >= s*threshold/2-1e-10)
            counts["escape_implications"] += int(escape.sum())
            counts["flip_worlds"] += len(bits)

            # At each prefix, integrate the remaining independent bits.
            # Each conditional difference is bounded by the literal edge
            # Lipschitz constant4, so each martingale increment has the
            # exact Bernoulli conditional variance rho(1-rho)d^2.
            future = f.copy()
            for level in range(len(ii)-1, -1, -1):
                children = future.reshape(-1, 2)
                diff = children[:, 1]-children[:, 0]
                assert np.max(np.abs(diff)) <= 4+1e-10
                means = (1-rho)*children[:, 0]+rho*children[:, 1]
                variance = ((1-rho)*(children[:, 0]-means)**2
                            +rho*(children[:, 1]-means)**2)
                assert np.max(np.abs(variance-rho*(1-rho)*diff**2)) < 1e-10
                assert np.all(variance <= 16*rho+1e-10)
                counts["martingale_nodes"] += len(diff)
                future = means
            probabilities = np.prod(np.where(bits == 1, rho, 1-rho), axis=1)
            assert abs(future[0]-probabilities @ f) < 1e-10

            for _ in range(10):
                gaussian = rng.normal(size=len(ii))
                response = np.abs(increments[outside] @ gaussian)
                nonzero = response > 1e-12
                amplitude = 0.49*s*np.min(deficit[outside][nonzero]/response[nonzero])
                f_gaussian = max(0, np.max(amplitude*(increments[outside] @ gaussian)
                                          -s*deficit[outside]/2))
                assert f_gaussian < 1e-9
                counts["gaussian_half_deficit"] += 1

    for n in (10**6, 10**9, 10**12, 10**15):
        for gamma in (0.05, 0.15, 0.25, 0.32, 0.45, 0.6, 2/3):
            eta = n**(-gamma)
            p = int(math.floor(eta**(-0.5)+1e-10))
            r = max(1, int(math.floor(math.sqrt(n)*eta**0.75+1e-10)))
            q = p*r
            b = n/math.sqrt(p)
            assert r >= math.sqrt(n)*eta**0.75/2-1e-8
            assert p >= eta**(-0.5)/2
            assert p*r*r <= n*(1+1e-10)
            assert q <= math.sqrt(n)*eta**0.25*(1+1e-10)
            assert eta*n**1.5/r <= 2*n*eta**0.25*(1+1e-10)
            assert b <= math.sqrt(2)*n*eta**0.25*(1+1e-10)
            counts["clone_parameter_checks"] += 1

    for d in range(2, 129):
        for k in range(1,d):
            rho=k/d
            logprob=(math.lgamma(d+1)-math.lgamma(k+1)-math.lgamma(d-k+1)
                     +k*math.log(rho)+(d-k)*math.log1p(-rho))
            assert math.exp(logprob)>=1/(d+1)-1e-12
            assert (d-k+1)/(d-k)>1
            assert k/(k+1)<1
            counts['binomial_mode_checks']+=1

    result = {"status": "PASS", "checks": counts,
              "scope": "Exact support, signed escape implication, Doob variance, and algebra; no asymptotic theorem inferred from sampling."}
    output = Path("tmp/paper_portfolio_2026_09_17/localization/actual_sign_noise_stability_audit.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
