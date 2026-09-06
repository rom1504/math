"""Checks for the rooted-walk / proportional-thinning theorem.

Run from the repository with .venv/bin/python -B and OPENBLAS_NUM_THREADS=1.
The graph checks use exact integer sums. Gaussian energies are numerical
expectations, not certificates of exact Boolean maxima. No file is written
unless --output is explicitly supplied.
"""

import argparse
import itertools
import json
import math

import numpy as np


def conference_submatrix(prime, labels):
    """Paley order prime+1, infinity labelled prime; prime must be 1 mod4."""
    assert prime % 4 == 1
    residues = np.zeros(prime, dtype=np.int8)
    for k in range(1, prime):
        residues[(k * k) % prime] = 1
    character = 2 * residues - 1
    character[0] = 0
    labels = np.asarray(labels, dtype=np.int64)
    finite = labels != prime
    a = np.ones((len(labels), len(labels)), dtype=np.int64)
    idx = np.flatnonzero(finite)
    difference = (labels[idx, None] - labels[None, idx]) % prime
    a[np.ix_(idx, idx)] = character[difference]
    np.fill_diagonal(a, 0)
    return a


def injective_sum(a, vertices, edges, root=0, root_label=0):
    others = [v for v in range(vertices) if v != root]
    available = [v for v in range(len(a)) if v != root_label]
    value = 0
    for lab in itertools.permutations(available, len(others)):
        mapping = dict(zip(others, lab))
        mapping[root] = root_label
        term = 1
        for u, v in edges:
            term *= int(a[mapping[u], mapping[v]])
        value += term
    return value


def degree_two_identity(a, vertices, edges, vertex, root=0):
    adjacent = []
    remaining = []
    for u, v in edges:
        if u == vertex:
            adjacent.append(v)
        elif v == vertex:
            adjacent.append(u)
        else:
            remaining.append((u, v))
    assert len(adjacent) == 2
    old_vertices = [v for v in range(vertices) if v != vertex]
    renumber = {v: j for j, v in enumerate(old_vertices)}

    def summed(new_edges):
        if any(u == v for u, v in new_edges):
            return 0
        new_edges = [(renumber[u], renumber[v]) for u, v in new_edges]
        return injective_sum(a, vertices - 1, new_edges, renumber[root])

    if adjacent[0] == adjacent[1]:
        rhs = (len(a) - vertices + 1) * summed(remaining)
    else:
        u, v = adjacent
        rhs = -sum(summed(remaining + [(u, w), (w, v)]) for w in old_vertices)
    lhs = injective_sum(a, vertices, edges, root)
    assert lhs == rhs, (vertices, edges, vertex, lhs, rhs)
    return {"vertices": vertices, "edges": edges, "eliminated": vertex, "sum": lhs}


def graph_checks():
    a = conference_submatrix(5, range(6))
    assert np.array_equal(a @ a, 5 * np.eye(6, dtype=np.int64))
    records = []
    # A doubled path, simple four-cycle, four-cycle plus doubled chord,
    # and a triangle with a doubled leaf. The last two sums vanish or
    # lose a full N power by the equality-case argument.
    cases = [
        (3, [(0, 1), (0, 1), (1, 2), (1, 2)], 2),
        (4, [(0, 1), (1, 2), (2, 3), (3, 0)], 1),
        (4, [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (0, 2)], 1),
        (4, [(0, 1), (1, 2), (2, 0), (2, 3), (2, 3)], 3),
    ]
    for case in cases:
        records.append(degree_two_identity(a, *case))
    return records


def semicircle_moment(k):
    if k % 2:
        return 0
    j = k // 2
    return math.comb(2 * j, j) // (j + 1)


def polynomial(degree):
    coeff = np.zeros(degree + 1)
    p0 = np.array([1.0])
    p1 = np.array([0.0, 1.0])
    ps = [p0]
    if degree:
        ps.append(p1)
    for j in range(2, degree + 1):
        nxt = np.r_[0.0, ps[-1]]
        nxt[: len(ps[-2])] -= ps[-2]
        ps.append(nxt)
    for j, pj in enumerate(ps):
        weight = math.sqrt(2 / (degree + 2)) * math.sin((j + 1) * math.pi / (degree + 2))
        coeff[: len(pj)] += weight * pj
    square = np.convolve(coeff, coeff)
    nu = sum(c * semicircle_moment(j) for j, c in enumerate(square))
    lam = sum(c * semicircle_moment(j + 1) for j, c in enumerate(square))
    expected = 2 * math.cos(math.pi / (degree + 2))
    assert abs(nu - 1) < 1e-10
    assert abs(lam - expected) < 1e-10
    return coeff, nu, lam


def gaussian_test(a, retention, parent_n, degree, delta):
    m = len(a)
    b = a / math.sqrt(m)
    coeff, nu, lam = polynomial(degree)
    ident = np.eye(m)
    g = coeff[-1] * ident
    for coefficient in coeff[-2::-1]:
        g = g @ b + coefficient * ident
    k = g @ g
    diag = np.diag(k)
    good = diag <= nu + delta
    covariance = k * good[:, None] * good[None, :] / (nu + delta)
    np.fill_diagonal(covariance, 1.0)
    covariance = np.clip(covariance, -1.0, 1.0)
    rounded = float(np.sum(a * np.arcsin(covariance)) / math.pi / m**1.5)
    normalized = k / np.sqrt(diag[:, None] * diag[None, :])
    np.fill_diagonal(normalized, 1.0)
    standard_rounding = float(np.sum(a * np.arcsin(np.clip(normalized, -1, 1))) / math.pi / m**1.5)
    return {
        "m": m,
        "retained_fraction": m / parent_n,
        "degree": degree,
        "target_iterated_limit": lam / math.pi,
        "diag_mean": float(np.mean(diag)),
        "diag_mean_squared_error_from_1": float(np.mean((diag - 1) ** 2)),
        "discarded_fraction": float(np.mean(~good)),
        "trimmed_gaussian_energy": rounded,
        "standard_gaussian_energy": standard_rounding,
        "linear_trace_energy": float(np.trace(b @ k) / m / math.pi),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=1009)
    parser.add_argument("--retention", type=float, default=0.1)
    parser.add_argument("--samples", type=int, default=4)
    parser.add_argument("--degree", type=int, default=3)
    parser.add_argument("--delta", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=9060737)
    parser.add_argument("--bernoulli", action="store_true")
    parser.add_argument("--output")
    args = parser.parse_args()
    # Trial division is enough at the modest prime orders used here.
    assert args.prime >= 5 and all(args.prime % d for d in range(2, math.isqrt(args.prime) + 1))
    rng = np.random.default_rng(args.seed)
    records = []
    for _ in range(args.samples):
        if args.bernoulli:
            labels = np.flatnonzero(rng.random(args.prime + 1) < args.retention)
        else:
            labels = np.sort(rng.choice(args.prime + 1, int(args.retention * (args.prime + 1)), replace=False))
        if len(labels) < 2:
            continue
        a = conference_submatrix(args.prime, labels)
        records.append(gaussian_test(a, args.retention, args.prime + 1, args.degree, args.delta))
    result = {
        "parameters": vars(args),
        "exact_degree_two_identity_checks": graph_checks(),
        "polynomial_coefficients": polynomial(args.degree)[0].tolist(),
        "samples": records,
        "scope": "Exact graph identities; numerical Gaussian expected-energy lower bounds, not exact caps.",
    }
    rendered = json.dumps(result, indent=2)
    print(rendered)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as stream:
            stream.write(rendered + "\n")


if __name__ == "__main__":
    main()
