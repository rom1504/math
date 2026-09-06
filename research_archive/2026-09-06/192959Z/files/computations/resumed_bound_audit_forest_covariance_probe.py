"""Monte Carlo falsification probe, not a covariance-operator certificate.

Tests higher local Hermite responses on actual hollow Sylvester signings.
All responses use the exact own-free degree-three old tree field.
"""
import argparse
import json
import numpy as np
from scipy.linalg import hadamard


def probe(n, samples, batch, seed):
    rng = np.random.default_rng(seed)
    a = hadamard(n).astype(float)
    np.fill_diagonal(a, 0.0)
    b = a / np.sqrt(n - 1)
    q = b @ b
    dirs = [np.ones(n) / np.sqrt(n)]
    dirs.extend(b[j] for j in [0, 1, min(7, n - 1)])
    for _ in range(4):
        w = rng.choice([-1.0, 1.0], size=n)
        dirs.append(w / np.linalg.norm(w))
    directions = np.stack(dirs, axis=1)
    target = np.einsum("ia,ij,ja->a", directions, q ** 3, directions)
    sum_y = np.zeros((3, len(dirs)))
    sum_y2 = np.zeros_like(sum_y)
    sum_y4 = np.zeros_like(sum_y)
    sum_local2 = np.zeros(3)
    sum_local_linear = np.zeros((3, 2))
    seen = 0
    while seen < samples:
        count = min(batch, samples - seen)
        spins = rng.choice([-1.0, 1.0], size=(count, n))
        g = spins @ b
        raw_star = (spins * (g * g - 1)) @ b / np.sqrt(2)
        common = np.einsum("bi,bi->b", spins, g)
        own_coeff = np.sqrt(2) / (n - 1) * (
            common[:, None] - 2 * spins * g
        )
        star = raw_star - spins * own_coeff
        responses = [
            (g ** 3 - 3 * g) / np.sqrt(6),
            g * (star ** 2 - 1) / np.sqrt(2),
            (star ** 3 - 3 * star) / np.sqrt(6),
        ]
        for k, response in enumerate(responses):
            projected = response @ directions
            sum_y[k] += projected.sum(axis=0)
            sum_y2[k] += (projected ** 2).sum(axis=0)
            sum_y4[k] += (projected ** 4).sum(axis=0)
            sum_local2[k] += np.mean(response ** 2, axis=1).sum()
            sum_local_linear[k, 0] += np.mean(response * g, axis=1).sum()
            sum_local_linear[k, 1] += np.mean(response * star, axis=1).sum()
        seen += count
    second = sum_y2 / samples
    stderr = np.sqrt(np.maximum(0, sum_y4 / samples - second ** 2) / samples)
    return {
        "n": n,
        "samples": samples,
        "method": "Monte Carlo diagnostic, not proof",
        "response_names": ["h3_edge", "edge_h2_star", "h3_star"],
        "direction_names": ["constant", "row0", "row1", "row7"]
        + [f"random{i}" for i in range(4)],
        "predicted_schur_cubic": target.tolist(),
        "empirical_projected_second_moments": second.tolist(),
        "second_moment_standard_errors": stderr.tolist(),
        "local_second_moments": (sum_local2 / samples).tolist(),
        "local_edge_star_first_coefficients": (sum_local_linear / samples).tolist(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orders", nargs="+", type=int, default=[64, 128, 256])
    parser.add_argument("--samples", type=int, default=65536)
    parser.add_argument("--batch", type=int, default=1024)
    args = parser.parse_args()
    for n in args.orders:
        print(json.dumps(probe(n, args.samples, args.batch, 20260906 + n)), flush=True)


if __name__ == "__main__":
    main()
