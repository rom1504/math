"""One full-gradient innovation update; quadrature exploration, not certification."""

import argparse
import json
import math
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import ndtr, ndtri, roots_hermitenorm, roots_legendre


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nodes", type=int, default=512)
    ap.add_argument("--rectangle-nodes", type=int, default=128)
    args = ap.parse_args()
    x, wx = roots_hermitenorm(args.nodes)
    wx /= math.sqrt(2 * math.pi)
    h2 = (x * x - 1) / math.sqrt(2)
    phi = lambda z: np.exp(-z * z / 2) / math.sqrt(2 * math.pi)
    psi = lambda k, t: k * (2 * ndtr(k / math.sqrt(t)) - 1) + 2 * math.sqrt(t) * phi(k / math.sqrt(t))
    # These rounded constants DEFINE the starting response; its moments
    # are recomputed below, not assumed equal to the constants.
    old_c = np.array([0.479929033657, -0.249506953601])
    old_a = np.array([0.518240271543, -0.193133028234])
    old_t, old_b = 0.086387266103, 0.099365026252
    old_threshold = np.maximum(psi(old_a[0] + old_a[1] * h2, old_t) - old_b, 0)
    old_mean, old_sd = old_c[0] * x, abs(old_c[1])
    old_zp = (old_threshold - old_mean) / old_sd
    old_zm = (-old_threshold - old_mean) / old_sd
    q = ndtr(-old_zp) + ndtr(old_zm)
    Hbar = 1 - q
    a = np.array([wx @ (x * (ndtr(-old_zp) - ndtr(old_zm))),
                  np.sign(old_c[1]) * (wx @ (phi(old_zp) + phi(old_zm)))])
    p = float(wx @ q)
    t = p - a @ a
    K = a[0] + a[1] * h2
    s = 2 * ndtr(K / math.sqrt(t)) - 1
    B = float(wx @ (Hbar * phi(K / math.sqrt(t)))) / math.sqrt(t)
    A = np.array([wx @ (Hbar * s), wx @ (Hbar * s * h2)]) - 2 * B * a
    wnorm2 = float(wx @ (Hbar * s * s - 4 * B * Hbar * K * s + 4 * B * B * K * K))
    sigma2 = wnorm2 - A @ A
    sigma = math.sqrt(sigma2)
    hnew_H = (s - 2 * B * K - A[0] - A[1] * h2) / sigma
    hnew_O = (-2 * B * K - A[0] - A[1] * h2) / sigma
    new_threshold_raw = psi(K, t) - B
    new_threshold = np.maximum(new_threshold_raw, 0)
    new_sd = math.sqrt(A[1] ** 2 + sigma2)
    zp = (new_threshold - A[0] * x) / new_sd
    zm = (-new_threshold - A[0] * x) / new_sd
    qnew = ndtr(-zp) + ndtr(zm)
    densities = phi(zp) + phi(zm)
    an = np.array([wx @ (x * (ndtr(-zp) - ndtr(zm))),
                   A[1] / new_sd * (wx @ densities),
                   sigma / new_sd * (wx @ densities)])
    pn = float(wx @ qnew)
    # Integrate old-center/new-support joint probability. Transform the
    # old G2 center interval into a bounded normal-CDF interval.
    y0 = (-old_threshold - old_mean) / old_c[1]
    y1 = (old_threshold - old_mean) / old_c[1]
    ulo, uhi = ndtr(np.minimum(y0, y1)), ndtr(np.maximum(y0, y1))
    gl, gw = roots_legendre(args.rectangle_nodes)
    uniforms = ulo[:, None] + (uhi - ulo)[:, None] * (gl[None, :] + 1) / 2
    y = ndtri(np.clip(uniforms, 1e-300, 1 - 1e-16))
    means = A[0] * x[:, None] + A[1] * y
    support_cond = ndtr((means - new_threshold[:, None]) / sigma) + ndtr((-means - new_threshold[:, None]) / sigma)
    Hqnew = (uhi - ulo) * (support_cond @ gw) / 2

    def value(theta):
        at = (1 - theta) * np.r_[a, 0.] + theta * an
        pt = (1 - theta) * p + theta * pn
        tt = pt - at @ at
        base = at[0] + at[1] * h2
        center_mass = Hbar - theta * Hqnew
        outer_mass = theta * (q - qnew + Hqnew)
        return float(wx @ (center_mass * psi(base + at[2] * hnew_H, tt)
                          + outer_mass * psi(base + at[2] * hnew_O, tt)))

    opt = minimize_scalar(lambda v: -value(v), bounds=(0, 1), method="bounded", options={"xatol": 1e-13})
    C = value(0)
    linear_max = wx @ (A[0] * x * (ndtr(-zp) - ndtr(zm)) + new_sd * densities - new_threshold_raw * qnew)
    gap = float(linear_max - A @ a - wx @ ((B - psi(K, t)) * q))
    print(json.dumps({"status": "quadrature heuristic only", "nodes": args.nodes,
                      "rectangle_nodes": args.rectangle_nodes, "initial": C,
                      "full_step": value(1), "best_theta": float(opt.x), "best_value": -float(opt.fun),
                      "gain": -float(opt.fun) - C, "full_gradient_gap": gap,
                      "innovation_variance": sigma2, "a": a.tolist(), "A": A.tolist(),
                      "B": B, "p": p, "tau2": float(t), "new_a": an.tolist(), "new_p": pn}, indent=2))


if __name__ == "__main__":
    main()
