"""Exhaustive finite cube audit of second-layer linear corrections."""

import json
import math
import numpy as np


def audit(A):
    n = len(A)
    m = n - 1
    B = A / math.sqrt(m)
    Q = B @ B
    B3 = Q @ B
    B4 = Q @ Q
    masks = np.arange(1 << n, dtype=np.uint64)
    S = (2 * ((masks[:, None] >> np.arange(n, dtype=np.uint64)) & 1)
         .astype(np.int64) - 1)
    G = S @ B
    Y = (S * (G * G - 1)) @ B / math.sqrt(2)
    W = (S * (Y * Y - 1)) @ B / math.sqrt(2)
    covY = Y.T @ Y / len(S)
    predCov = ((m - 3) / m) * Q + (2 / m) * B4
    observed_pairs = (S * (Y * Y)).T @ S / len(S)
    predicted_pairs = (4 / m) * B * np.diag(B3)[None, :]
    predicted_pairs += (8 * (m - 3) / (m * m)) * Q
    # Formula is E[S_j S_l Y_j²]; its diagonal is Var(Y_j), not the
    # off-diagonal expression.  B3's diagonal depends on l, hence columns.
    np.fill_diagonal(predicted_pairs, np.diag(predCov))
    predL = (B * (2 * np.diag(B4)[None, :] - 3) / m
             + (4 / m) * Q * np.diag(B3)[None, :]
             + (8 * (m - 3) / (m * m)) * (B3 - B)) / math.sqrt(2)
    observedL = W.T @ S / len(S)
    errors = {
        "covariance_Y": float(np.max(np.abs(covY - predCov))),
        "quadratic_endpoint": float(np.max(np.abs(observed_pairs
                                                   - predicted_pairs))),
        "linear_W": float(np.max(np.abs(observedL - predL))),
    }
    assert max(errors.values()) < 1e-10, errors
    return errors


def main():
    seed = 20260905
    rng = np.random.default_rng(seed)
    records = []
    for n in range(3, 11):
        for sample in range(5):
            if sample == 0:
                A = np.ones((n, n), dtype=np.int64) - np.eye(n, dtype=np.int64)
            else:
                A = rng.choice([-1, 1], size=(n, n))
                A = np.triu(A, 1)
                A += A.T
            records.append({"n": n, "sample": sample, "errors": audit(A)})
    print(json.dumps({"seed": seed, "signings": len(records),
                      "verified": True, "records": records}, indent=2))


if __name__ == "__main__":
    main()
