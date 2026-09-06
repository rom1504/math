"""Exact field checks for the fixed-basis, seed-only tail counterexample."""

import json

import numpy as np


H_PAIR = (np.array([[1, 1], [-1, 1]], dtype=np.int64),
          np.array([[1, 1], [1, -1]], dtype=np.int64))


def hadamard(order):
    out = np.ones((1, 1), dtype=np.int64)
    while len(out) < order:
        out = np.block([[out, out], [out, -out]])
    assert len(out) == order
    return out


def base_basis(index):
    out = np.ones((1, 1), dtype=np.int64)
    for position in range(5, -1, -1):
        out = np.kron(out, H_PAIR[(index >> position) & 1])
    assert np.array_equal(out@out.T, 64*np.eye(64, dtype=np.int64))
    return out


def check(t, rng):
    base = [base_basis(i) for i in range(64)]
    carrier = hadamard(t)
    base_phase = np.array([1-2*(bin(i >> 1).count("1") % 2) for i in range(64)])
    seed_phase = np.array([1-2*(bin(i).count("1") % 2) for i in range(64)])
    base_seed = np.outer(seed_phase, seed_phase)
    x_base = np.repeat(base_phase[:, None], 62, axis=1)
    spectra_base = np.array([base[i][:62].T@x_base[i] for i in range(64)])
    fields_base = np.array([base[i][:62]@(base_seed[i]*spectra_base[:, i]) for i in range(64)])
    assert np.array_equal(fields_base, -64*x_base)
    m, kept = 64*t, 62*t
    bases = [np.kron(matrix, carrier)[:kept] for matrix in base]
    seed = rng.choice([-1, 1], (m, m))
    seed = np.triu(seed, 1)+np.triu(seed, 1).T+np.eye(m, dtype=np.int64)
    for u in range(t):
        seed[u::t, u::t] = base_seed
    x = np.array([np.kron(x_base[i], carrier[:, u]) for i in range(64) for u in range(t)])
    assert np.all(np.abs(x) == 1)
    spectra = np.array([bases[i//t].T@x[i] for i in range(m)])
    expected = np.zeros_like(spectra)
    for i in range(64):
        for u in range(t):
            expected[i*t+u, u::t] = t*spectra_base[i]
    assert np.array_equal(spectra, expected)
    fields = np.array([bases[i//t]@(seed[i]*spectra[:, i]) for i in range(m)])
    assert np.array_equal(fields, -m*x)
    assert np.array_equal(fields-x, -(m+1)*x)
    n = m*kept
    cap = (m+1)*n//2
    assert -int(np.sum(x*(fields-x))) == 2*cap
    return {"t": t, "fibre_order": m, "retained_order": n,
            "negative_hollow_eigenvalue": -(m+1), "exact_global_cap": cap,
            "event_probability_negative_log2": 2016*t,
            "aligned_field_dispersion": 0}


def main():
    rng = np.random.default_rng(260906)
    print(json.dumps({"seed": 260906, "cases": [check(t, rng) for t in (1, 2, 4, 8)],
                      "status": "integer spectrum, field and global-cap witness identities pass"}, indent=2))


if __name__ == "__main__":
    main()
