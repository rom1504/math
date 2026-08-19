#!/usr/bin/env python3
"""Audit labelled Fourier switch certificates on exact small child minimizers.

The bridge is enumerated modulo row/column switching by fixing its first row
and column positive.  Every remaining switch is recovered by minimizing the
convolution over the rank-one switching group.  This is intended only as a
finite falsifier for the one-character certificate, not as asymptotic
evidence.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np


def spins(n):
    z = np.arange(1 << n, dtype=np.uint64)[:, None]
    bits = ((z >> np.arange(n, dtype=np.uint64)) & 1).astype(np.int8)
    return 1 - 2 * bits


def energies(A, X):
    return np.einsum("bi,ij,bj->b", X, A, X, optimize=True) / 2.0


def logmeanexp(v, axis=-1):
    vmax = np.max(v, axis=axis, keepdims=True)
    return np.squeeze(vmax, axis=axis) + np.log(
        np.mean(np.exp(v - vmax), axis=axis)
    )


def pressure(A, u, X):
    h = energies(A, X)
    return float(logmeanexp(np.logaddexp(u * h, -u * h) - math.log(2.0)))


def enumerate_child_minimizer(n, beta):
    X = spins(n)
    interior = [(i, j) for i in range(1, n) for j in range(i + 1, n)]
    best = float("inf")
    best_A = None
    u = beta / math.sqrt(n)
    for mask in range(1 << len(interior)):
        A = np.ones((n, n), dtype=np.int8)
        np.fill_diagonal(A, 0)
        for bit, (i, j) in enumerate(interior):
            if (mask >> bit) & 1:
                A[i, j] = A[j, i] = -1
        val = pressure(A, u, X)
        if val < best - 1e-13:
            best, best_A = val, A.copy()
    return best_A, best


def fwht_rows(a, normalized):
    """Walsh transform along the final axis, on a fresh float array."""
    out = np.array(a, dtype=np.float64, copy=True)
    n = out.shape[-1]
    h = 1
    while h < n:
        view = out.reshape(out.shape[:-1] + (-1, 2 * h))
        left = view[..., :h].copy()
        right = view[..., h : 2 * h].copy()
        view[..., :h] = left + right
        view[..., h : 2 * h] = left - right
        h *= 2
    if normalized:
        out /= n
    return out


def audit(n, beta, batch_size):
    X = spins(n)
    A, P = enumerate_child_minimizer(n, beta)
    t = beta / math.sqrt(2 * n)

    # Unique group representatives: x is arbitrary and y_0=+1.
    Ytail = spins(n - 1)
    Y = np.ones((1 << (n - 1), n), dtype=np.int8)
    Y[:, 1:] = Ytail
    Xg = np.repeat(X, len(Y), axis=0)
    Yg = np.tile(Y, (len(X), 1))
    group_size = len(Xg)

    h_x = energies(A, Xg)
    h_y = energies(A, Yg)

    variable_edges = [(i, j) for i in range(1, n) for j in range(1, n)]
    features = np.stack(
        [Xg[:, i] * Yg[:, j] for i, j in variable_edges], axis=1
    ).astype(np.float64)
    base_energy = Xg.sum(axis=1).astype(np.float64) * Yg.sum(axis=1)

    total_bridges = 1 << len(variable_edges)
    records = []
    for epsilon in (1, -1):
        w = np.cosh(t * (h_x + epsilon * h_y))
        mean_w = float(np.mean(w))
        a = w / mean_w
        ahat = fwht_rows(a[None, :], normalized=True)[0]

        best_one = {"certificate": float("inf")}
        best_full = {"certificate": float("inf")}
        for start in range(0, total_bridges, batch_size):
            stop = min(total_bridges, start + batch_size)
            masks = np.arange(start, stop, dtype=np.uint64)[:, None]
            bit_positions = np.arange(len(variable_edges), dtype=np.uint64)
            bits = ((masks >> bit_positions) & 1).astype(np.float64)
            bridge_energy = base_energy[None, :] - 2.0 * bits.dot(features.T)
            scaled = t * bridge_energy
            psi = logmeanexp(scaled, axis=1)
            row_max = np.max(scaled, axis=1, keepdims=True)
            kval = np.exp(scaled - row_max)
            kval /= np.mean(kval, axis=1, keepdims=True)
            bhat = fwht_rows(kval, normalized=True)

            products = bhat * ahat[None, :]
            max_mode_index = np.argmax(np.abs(products[:, 1:]), axis=1) + 1
            max_mode = np.abs(products[np.arange(stop - start), max_mode_index])
            one_log_gain = np.log(np.maximum(1.0 - max_mode, 1e-300))
            base = math.log(mean_w) + psi - 2.0 * P
            one_cert = base + one_log_gain

            # Fourier inversion of the exact normalized convolution.
            convolution = fwht_rows(products, normalized=False)
            min_conv = np.maximum(np.min(convolution, axis=1), 1e-300)
            full_cert = base + np.log(min_conv)

            j = int(np.argmin(one_cert))
            if one_cert[j] < best_one["certificate"]:
                best_one = {
                    "certificate": float(one_cert[j]),
                    "bridge_gauge_mask": int(start + j),
                    "mode_index": int(max_mode_index[j]),
                    "mode_product_abs": float(max_mode[j]),
                    "base_before_switch_gain": float(base[j]),
                    "log_gain": float(one_log_gain[j]),
                }
            j = int(np.argmin(full_cert))
            if full_cert[j] < best_full["certificate"]:
                best_full = {
                    "certificate": float(full_cert[j]),
                    "bridge_gauge_mask": int(start + j),
                    "base_before_switch_gain": float(base[j]),
                    "log_gain": float(math.log(min_conv[j])),
                }

        records.append(
            {
                "epsilon": epsilon,
                "mean_internal_kernel": mean_w,
                "best_one_character": best_one,
                "best_full_switch_convolution": best_full,
            }
        )

    return {
        "child_order": n,
        "parent_order": 2 * n,
        "beta": beta,
        "contracted_temperature": t,
        "child_pressure": P,
        "child_matrix": A.tolist(),
        "switch_group_size": group_size,
        "bridge_switching_classes": total_bridges,
        "orientations": records,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=4)
    parser.add_argument("--beta", type=float, default=4.0)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.n > 5:
        raise SystemExit("Exact bridge enumeration is intentionally capped at n=5")
    result = audit(args.n, args.beta, args.batch_size)
    payload = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n")
    print(payload)


if __name__ == "__main__":
    main()
