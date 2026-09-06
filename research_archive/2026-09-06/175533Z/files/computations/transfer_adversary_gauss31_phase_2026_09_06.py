"""Exact three-phase independence certificate and bounded lower-witness search."""

import argparse
import json

import numpy as np


def multiply(a, b, modulus):
    product = [0] * 9
    for i in range(5):
        for j in range(5):
            product[i + j] += a[i] * b[j]
    # Work in (Z/modulus)[x]/(x^5+x^2+1).
    for k in range(8, 4, -1):
        product[k - 5] -= product[k]
        product[k - 3] -= product[k]
    return tuple(v % modulus for v in product[:5])


def power(a, exponent, modulus):
    result = (1, 0, 0, 0, 0)
    while exponent:
        if exponent & 1:
            result = multiply(result, a, modulus)
        a = multiply(a, a, modulus)
        exponent >>= 1
    return result


def exact_certificate():
    x = (0, 1, 0, 0, 0)
    field_orbit = [power(x, exponent, 2) for exponent in range(31)]
    assert len(set(field_orbit)) == 31
    assert power(x, 31, 2) == (1, 0, 0, 0, 0)
    signs = []
    for element in field_orbit:
        trace = [0] * 5
        term = element
        for _ in range(5):
            trace = [(a + b) % 2 for a, b in zip(trace, term)]
            term = multiply(term, term, 2)
        assert trace[1:] == [0] * 4
        signs.append(1 - 2 * trace[0])
    root = power(x, 32, 32)
    assert power(root, 31, 32) == (1, 0, 0, 0, 0)
    assert tuple(v % 2 for v in root) == x
    representatives = [1, 3, 5]
    valuation_matrix = []
    residues = []
    for embedding in representatives:
        row, residue_row = [], []
        for character in representatives:
            terms = [power(root, embedding * character * t, 32) for t in range(31)]
            gauss = tuple(sum(signs[t] * terms[t][k] for t in range(31)) % 32 for k in range(5))
            assert any(gauss)
            valuation = min((v & -v).bit_length() - 1 for v in gauss if v)
            row.append(2 * valuation - 5)
            residue_row.append(gauss)
        valuation_matrix.append(row)
        residues.append(residue_row)
    a = valuation_matrix
    determinant = (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
                   - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                   + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))
    assert determinant != 0
    orbits = []
    for representative in representatives:
        orbit = sorted({representative * 2 ** j % 31 for j in range(5)})
        orbits.extend([orbit, sorted({-j % 31 for j in orbit})])
    assert sorted(sum(orbits, [])) == list(range(1, 31))
    return {"field_polynomial": "x^5+x^2+1", "primitive_orbit_length": 31,
            "trace_signs": signs, "root_mod32": root,
            "squared_phase_valuation_matrix": valuation_matrix,
            "determinant": determinant, "gauss_residues_mod32": residues,
            "paired_frobenius_orbits": orbits}


def phase_components():
    indices = np.arange(31)
    angle = 2 * np.pi * (indices[:, None] + indices[None, :]) / 31
    real, imaginary = [], []
    for representative in [1, 3, 5]:
        orbit = {representative * 2 ** j % 31 for j in range(5)}
        real.append(sum((2 * np.cos(j * angle) / 31 for j in orbit)))
        imaginary.append(sum((-2 * np.sin(j * angle) / 31 for j in orbit)))
    return np.array(real), np.array(imaginary)


def phase_search(restarts, iterations, seed):
    rng = np.random.default_rng(seed)
    real, imaginary = phase_components()
    matrix = np.array([[-1.0, 2.0], [2.0, 4.0]])
    best, witness = 0.0, None
    for restart in range(restarts):
        f = rng.choice([-1.0, 1.0], size=(31, 2))
        g = rng.choice([-1.0, 1.0], size=(31, 2))
        for step in range(iterations):
            alpha = np.array([np.sum(g * (r @ f @ matrix)) for r in real])
            beta = np.array([np.sum(g * (r @ f @ matrix)) for r in imaginary])
            theta = np.arctan2(beta, alpha)
            kernel = np.einsum("j,jab->ab", np.cos(theta), real) + np.einsum("j,jab->ab", np.sin(theta), imaginary)
            g = np.where(kernel @ f @ matrix >= 0, 1.0, -1.0)
            f = np.where(kernel @ g @ matrix >= 0, 1.0, -1.0)
            value = float(np.sum(g * (kernel @ f @ matrix)) / 62)
            if value > best:
                best = value
                witness = {"theta": theta.tolist(), "f": f.astype(int).tolist(), "g": g.astype(int).tolist()}
            if step % 10 == 9:
                f *= 1 - 2 * (rng.random(f.shape) < rng.choice([1 / 62, 2 / 62, 4 / 62, 0.25]))
                g *= 1 - 2 * (rng.random(g.shape) < rng.choice([1 / 62, 2 / 62, 4 / 62, 0.25]))
    return {"status": "numerical lower witness only", "restarts": restarts,
            "iterations": iterations, "seed": seed, "lower_witness": best,
            "unamplified_value": 3.5, "witness": witness}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--restarts", type=int, default=0)
    parser.add_argument("--iterations", type=int, default=300)
    parser.add_argument("--seed", type=int, default=260906)
    args = parser.parse_args()
    print(json.dumps({"exact": exact_certificate()}), flush=True)
    if args.restarts:
        print(json.dumps({"numerical": phase_search(args.restarts, args.iterations, args.seed)}), flush=True)
