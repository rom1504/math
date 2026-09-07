"""Exhaustive small tensor-rotation variance-profile diagnostic."""

import json
import math

import numpy as np
from scipy.linalg import hadamard


reports = []
for h_order in [1, 2, 4, 8]:
    n = 2 * h_order
    words = ((np.arange(1 << (n - 1), dtype=np.uint64)[:, None]
              >> np.arange(n, dtype=np.uint64)[None, :]) & 1)
    signs = (1 - 2 * words.astype(np.int8)).astype(float)
    h = hadamard(h_order)
    angle_reports = []
    for theta in [0, math.pi / 16, math.pi / 8, 3 * math.pi / 16, math.pi / 4]:
        c, s = math.cos(theta), math.sin(theta)
        r = np.array([[c, s], [-s, c]])
        matrix = math.sqrt(2) * np.kron(r, h)
        response = np.abs(signs @ matrix).sum(axis=1)
        argmax = int(np.argmax(response))
        beta = float(response[argmax])
        assert np.allclose((matrix * matrix).sum(axis=1), n)
        angle_reports.append({
            "theta_over_pi": theta / math.pi,
            "bilinear_constant": beta / n ** 1.5,
            "witness_x": signs[argmax].astype(int).tolist(),
        })
    reports.append({"n": n, "angles": angle_reports})
print(json.dumps({"status": "exhaustive finite diagnostic", "reports": reports}, indent=2))

# A reproducible witness falsifies the tempting all-order cos(pi/8) upper
# suggested by the first four exhaustive orders.  This is only a lower
# witness; the alternating optimization does not certify a global maximum.
rng = np.random.default_rng(20260907)
witness_reports = []
for h_order in [16, 32, 64, 128]:
    n = 2 * h_order
    c, s = math.cos(math.pi / 8), math.sin(math.pi / 8)
    h = hadamard(h_order)
    matrix = math.sqrt(2) * np.kron(np.array([[c, s], [-s, c]]), h)
    x = rng.choice([-1, 1], size=(256, n))
    best = -1.0
    best_x = None
    for iteration in range(100):
        y = np.where(x @ matrix >= 0, 1, -1)
        x = np.where(y @ matrix.T >= 0, 1, -1)
        values = np.abs(x @ matrix).sum(axis=1)
        index = int(values.argmax())
        if values[index] > best:
            best = float(values[index])
            best_x = x[index].copy()
    best_y = np.where(best_x @ matrix >= 0, 1, -1)
    x1, x2 = np.split(best_x, 2)
    y1, y2 = np.split(best_y, 2)
    coeff_a = int(x1 @ h @ y1 + x2 @ h @ y2)
    coeff_b = int(x1 @ h @ y2 - x2 @ h @ y1)
    report = {
        "n": n,
        "certified_witness_coeff_a": coeff_a,
        "certified_witness_coeff_b": coeff_b,
        "witness_value_over_n_1_5": best / n ** 1.5,
    }
    if n == 64:
        report["witness_x"] = best_x.tolist()
        report["witness_y"] = best_y.tolist()
    witness_reports.append(report)
print(json.dumps({"status": "explicit lower witnesses", "reports": witness_reports}, indent=2))
