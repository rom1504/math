#!/usr/bin/env python3
"""Wave 26 audits for cancellation-preserving vertex deletion identities."""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import quad

from check_finite_bridge_r16 import A5
from check_response_dual_r16 import A8, A9, spins
from verify_compatible_replacement_r12 import A6


def oriented_states(A: np.ndarray):
    X0 = spins(len(A))
    X = np.tile(X0, (2, 1))
    sigma = np.repeat(np.asarray((-1, 1), dtype=np.int64), len(X0))
    energy = sigma * np.einsum("bi,ij,bj->b", X, A, X)
    return X, sigma, energy.astype(float)


def qnorm(A: np.ndarray) -> int:
    return int(np.max(oriented_states(A)[2]))


def cycle_polynomial_from_primal(A: np.ndarray, beta: float) -> float:
    _, _, energy = oriented_states(A)
    n = len(A)
    edge_count = n * (n - 1) // 2
    return float(
        np.exp(beta * energy).mean() / math.cosh(2 * beta) ** edge_count
    )


def multivariate_cycle_polynomial(A: np.ndarray, edge_parameters: np.ndarray) -> float:
    """Positive primal evaluation of the signed multivariate cycle sum."""
    n = len(A)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    edge_parameters = np.asarray(edge_parameters, dtype=float)
    assert len(edge_parameters) == len(edges)
    assert np.all(np.abs(edge_parameters) < 1)
    couplings = np.arctanh(edge_parameters)
    X, sigma, _ = oriented_states(A)
    exponent = np.zeros(len(X))
    for coupling, (i, j) in zip(couplings, edges):
        exponent += coupling * sigma * A[i, j] * X[:, i] * X[:, j]
    return float(np.exp(exponent).mean() / np.prod(np.cosh(couplings)))


def child_cavity(A: np.ndarray, i: int, beta: float):
    keep = [j for j in range(len(A)) if j != i]
    C = A[np.ix_(keep, keep)]
    row = A[i, keep]
    Y, sigma, energy = oriented_states(C)
    shift = float(np.max(beta * energy))
    weights = np.exp(beta * energy - shift)
    weights /= weights.sum()
    fields = Y @ row

    def moments(K: float):
        cosh = np.cosh(K * fields)
        sinh = np.sinh(K * fields)
        F = float(weights @ cosh)
        mean_h = float(weights @ (fields * sinh) / F)
        second_h = float(weights @ (fields * fields * cosh) / F)
        variance_h = second_h - mean_h * mean_h
        return F, mean_h, variance_h

    F, _, _ = moments(2 * beta)
    kappa = math.log(F) / beta
    v0 = float(weights @ (fields * fields))
    jensen = math.log(math.cosh(2 * beta * math.sqrt(v0))) / beta

    integral_mean = quad(lambda K: moments(K)[1], 0, 2 * beta, epsabs=1e-11)[0]
    integral_var = quad(
        lambda K: (2 * beta - K) * moments(K)[2],
        0,
        2 * beta,
        epsabs=1e-11,
    )[0]
    assert math.isclose(beta * kappa, integral_mean, rel_tol=2e-10, abs_tol=2e-10)
    assert math.isclose(beta * kappa, integral_var, rel_tol=2e-10, abs_tol=2e-10)
    assert kappa + 1e-11 >= jensen
    return {
        "kappa": kappa,
        "v0": v0,
        "jensen": jensen,
        "fields": fields,
        "weights": weights,
        "moments": moments,
    }


def audit_matrix(A: np.ndarray, name: str, beta: float) -> None:
    A = np.asarray(A, dtype=np.int64)
    n = len(A)
    q = qnorm(A)
    P = cycle_polynomial_from_primal(A, beta)
    kappas = []
    v0s = []
    jensen = []
    curvatures = []
    decrements = []
    deficit_log_ratios = []
    for i in range(n):
        data = child_cavity(A, i, beta)
        kappas.append(data["kappa"])
        v0s.append(data["v0"])
        jensen.append(data["jensen"])

        keep = [j for j in range(n) if j != i]
        child = A[np.ix_(keep, keep)]
        P_child = cycle_polynomial_from_primal(child, beta)
        ratio = P_child / P
        bc = ratio / math.cosh(2 * beta) ** (n - 1)
        assert math.isclose(bc, math.exp(-beta * data["kappa"]), rel_tol=2e-11)

        # Exact positive deficit-partition recurrence.
        _, _, parent_energy = oriented_states(A)
        _, _, child_energy = oriented_states(child)
        q_child = int(np.max(child_energy))
        D_parent = float(np.exp(-beta * (q - parent_energy)).sum())
        D_child = float(np.exp(-beta * (q_child - child_energy)).sum())
        decrement = q - q_child
        deficit_log_ratio = math.log(D_parent / (2 * D_child))
        assert math.isclose(
            data["kappa"],
            decrement + deficit_log_ratio / beta,
            rel_tol=2e-11,
            abs_tol=2e-11,
        )
        decrements.append(decrement)
        deficit_log_ratios.append(deficit_log_ratio)

        # The star high-temperature interpolation has exact deletion
        # curvature rho^2(v_i-(n-1)) at t=0.
        rho = math.tanh(2 * beta)
        curvature = rho * rho * (data["v0"] - (n - 1))
        curvatures.append(curvature)

        eps = 1e-3
        curvature_values = []
        edges = [(u, v) for u in range(n) for v in range(u + 1, n)]
        for t_small in (-eps, 0.0, eps):
            parameters = np.asarray(
                [rho * t_small if i in (u, v) else rho for u, v in edges]
            )
            curvature_values.append(
                math.log(multivariate_cycle_polynomial(A, parameters))
            )
        finite_curvature = (
            curvature_values[0] - 2 * curvature_values[1] + curvature_values[2]
        ) / (eps * eps)
        assert math.isclose(finite_curvature, curvature, rel_tol=2e-5, abs_tol=2e-5)

        alpha = q / n**1.5
        moment_threshold = 9 * alpha * alpha * n / 16
        if data["v0"] + 1e-12 >= moment_threshold:
            target = alpha * (n**1.5 - (n - 1) ** 1.5)
            assert data["kappa"] + 1e-11 >= target - math.log(2) / beta

        # Check the full interpolation identity at several interior points:
        # P_i(t)/P_i(0)=F_i(K(t))/cosh(K(t))^(n-1).
        for t in (0.0, 0.2, 0.6, 1.0):
            K = math.atanh(rho * t)
            F_t = data["moments"](K)[0]
            normalized = F_t / math.cosh(K) ** (n - 1)
            assert normalized > 0
            edges = [(u, v) for u in range(n) for v in range(u + 1, n)]
            parameters = np.asarray(
                [rho * t if i in (u, v) else rho for u, v in edges]
            )
            P_direct = multivariate_cycle_polynomial(A, parameters)
            assert math.isclose(P_direct / P_child, normalized, rel_tol=2e-11)
            if t == 1.0:
                assert math.isclose(P / P_child, normalized, rel_tol=2e-11)

        # The same multivariate interpolation is exactly the Boolean-noise
        # average of all star-competitor cycle polynomials.  Check this
        # independently at one interior point.
        t_noise = 0.37
        flip_probability = (1 - t_noise) / 2
        noise_average = 0.0
        for mask in range(1 << (n - 1)):
            twisted = A.copy()
            chosen = 0
            for bit, j in enumerate(keep):
                if (mask >> bit) & 1:
                    twisted[i, j] *= -1
                    twisted[j, i] *= -1
                    chosen += 1
            probability = (
                flip_probability**chosen
                * (1 - flip_probability) ** (n - 1 - chosen)
            )
            noise_average += probability * cycle_polynomial_from_primal(
                twisted, beta
            )
        edges = [(u, v) for u in range(n) for v in range(u + 1, n)]
        parameters = np.asarray(
            [rho * t_noise if i in (u, v) else rho for u, v in edges]
        )
        P_noise_direct = multivariate_cycle_polynomial(A, parameters)
        assert math.isclose(noise_average, P_noise_direct, rel_tol=3e-11)

    print(
        name,
        f"beta={beta}",
        {
            "n": n,
            "Q": q,
            "kappa_range": (min(kappas), max(kappas)),
            "child_field_second_moment_range": (min(v0s), max(v0s)),
            "jensen_lower_bound_range": (min(jensen), max(jensen)),
            "star_logP_curvature_range": (min(curvatures), max(curvatures)),
            "decrement_range": (min(decrements), max(decrements)),
            "log_parent_over_two_child_D_range": (
                min(deficit_log_ratios),
                max(deficit_log_ratios),
            ),
        },
    )


def main() -> None:
    for beta in (0.25, 0.5, 1.0):
        for A, name in ((A5, "A5"), (A6, "A6"), (A8, "A8"), (A9, "A9")):
            audit_matrix(A, name, beta)
    print("PASS: deletion polynomial, star interpolation, susceptibility, and Jensen bound")


if __name__ == "__main__":
    main()
