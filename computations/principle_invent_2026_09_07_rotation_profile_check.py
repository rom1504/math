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
