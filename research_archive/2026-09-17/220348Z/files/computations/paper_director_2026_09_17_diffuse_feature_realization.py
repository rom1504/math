"""Finite-binomial diagnostics for the physical diffuse-feature theorem.

No Monte Carlo and no solver. Probabilities/analytic integrals are floating
point diagnostics; exact target normalization/covariance follow the formulas
in the proof. Full covariance repair is checked using exchangeable blocks.
"""

import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import gammaln, logsumexp


ROOT = Path(__file__).resolve().parents[1]
V = 1 / 16
K = 64
P = 15 / 128
assert (1 - P) * V + P * (2 * K + 1) * V == 1
LOG_MOMENT = gammaln(2 * K + 1) - K * math.log(2) - gammaln(K + 1)


def log_weight(y):
    y = np.asarray(y, dtype=float)
    with np.errstate(divide="ignore"):
        log_polynomial = np.logaddexp(
            math.log1p(-P),
            math.log(P) + 2 * K * np.log(np.abs(y))
            - LOG_MOMENT - K * math.log(V),
        )
    return -.5 * math.log(V) - .5 * (1 / V - 1) * y * y + log_polynomial


def main():
    kappa = math.sqrt(2 / math.pi)
    active_mean = math.exp(
        (2 * K + .5) * math.log(2) + 2 * gammaln(K + 1)
        - .5 * math.log(math.pi) - gammaln(2 * K + 1)
    )
    target_mean = math.sqrt(V) * ((1 - P) * kappa + P * active_mean)

    def target_kl_integrand(y):
        lw = float(log_weight(y))
        return math.exp(-y * y / 2 - .5 * math.log(2 * math.pi) + lw) * lw

    target_kl, quadrature_error = quad(target_kl_integrand, -12, 12,
                                     epsabs=1e-11, points=[-4, -3, -2, 0, 2, 3, 4])
    rows = []
    for n in (16, 32, 64, 128, 256, 512, 1024, 4096, 16384):
        j = np.arange(n + 1)
        y = (n - 2 * j) / math.sqrt(n)
        log_uniform = (gammaln(n + 1) - gammaln(j + 1)
                       - gammaln(n - j + 1) - n * math.log(2))
        lw = log_weight(y)
        log_z = float(logsumexp(log_uniform + lw))
        weights = np.exp(log_uniform + lw - log_z)
        var = float(weights @ (y * y))
        response = float(weights @ np.abs(y))
        kl_tilt = float(weights @ lw - log_z)
        off = (var - 1) / (n - 1)
        frob = abs(off) * math.sqrt(n * (n - 1))
        theta = 4 * frob / (1 + 4 * frob)
        if theta:
            d = -(1 - theta) * off / theta
            latent_off = math.sin(math.pi * d / 2)
            eig_small = 1 - latent_off
            eig_large = 1 + (n - 1) * latent_off
            assert min(eig_small, eig_large) > .5
            assert max(eig_small, eig_large) < 1.5
            repair_var = 1 + (n - 1) * d
            repaired_off = (1 - theta) * off + theta * d
            assert abs(repaired_off) < 1e-14
            repair_kl_upper = -.5 * ((n - 1) * math.log1p(-latent_off)
                                       + math.log1p((n - 1) * latent_off))
        else:
            repair_var, repair_kl_upper = 1, 0
        rows.append({
            "n": n, "normalization": math.exp(log_z),
            "tilted_feature_variance": var,
            "covariance_frobenius_error": frob,
            "repair_mass": theta,
            "tilted_query_response": response,
            "exact_isotropy_repaired_response_upper":
                (1 - theta) * response + theta * math.sqrt(repair_var),
            "tilt_KL": kl_tilt,
            "repaired_KL_upper": (1 - theta) * kl_tilt + theta * repair_kl_upper,
            "target_response": target_mean,
            "target_KL_quadrature": target_kl,
        })
    result = {"status": "PASS finite binomial floating diagnostics; analytic proof separate",
              "target": {"v_exact": "1/16", "k": K, "p_exact": "15/128",
                         "exact_variance_identity": True,
                         "absolute_mean": target_mean,
                         "KL_quadrature": target_kl,
                         "quadrature_error_reported": quadrature_error},
              "cases": rows,
              "scope": "One diffuse feature, no original-signing optimizer assertion."}
    output = ROOT / "tmp/paper_portfolio_2026_09_17/director/diffuse_feature_realization.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
