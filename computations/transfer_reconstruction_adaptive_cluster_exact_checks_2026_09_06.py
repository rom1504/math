"""Bounded exact checks for adaptive product clusters and their obstructions."""

from fractions import Fraction as F
from itertools import product
import json

from continued_feedback_ternary_latent_interval_certificate_2026_09_06 import (
    entropy_interval, log_interval,
)


def sign(value):
    return (value > 0)-(value < 0)


def adaptive_checks():
    n, root_n = 4, F(2)
    beta, c_bound = F(1), F(3, 2)  # Every hollow order-four signing.
    b = beta/root_n
    delta = F(1, 8)
    q = 1-2*delta
    k_bound = 1/(delta*(1-delta/2))+2*beta*c_bound
    spins = list(product((-1, 1), repeat=n))
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    matrices, nonconstant = 0, 0
    for signs in product((-1, 1), repeat=len(edges)):
        energies = {x: sum(a*x[i]*x[j] for a, (i, j) in zip(signs, edges))
                    for x in spins}
        x = max(spins, key=lambda z: abs(energies[z]))
        sigma = sign(energies[x])
        cap = abs(energies[x])
        matrix = [[F(0) for _ in range(n)] for _ in range(n)]
        for a, (i, j) in zip(signs, edges):
            matrix[i][j] = matrix[j][i] = F(sigma*a*x[i]*x[j])
        fields = [sum(row)/root_n for row in matrix]
        assert min(fields) >= 0
        mean = sum(fields)/n
        assert mean == F(2*cap, n)/root_n
        d1 = sum(abs(v-mean) for v in fields)/n
        signs_field = [F(sign(v-mean)) for v in fields]
        signs_mean = sum(signs_field)/n
        epsilon = min(delta/4, beta*q*d1/k_bound)
        u = [-epsilon*(v-signs_mean) for v in signs_field]
        assert sum(u) == 0
        assert sum(v*v for v in u) <= n*epsilon**2
        probabilities = [delta+v for v in u]
        assert min(probabilities) >= delta/2
        assert max(probabilities) <= 3*delta/2
        quadratic_u = sum(matrix[i][j]*u[i]*u[j]
                          for i in range(n) for j in range(n))
        expected_energy = q*q*cap-2*q*root_n*sum(
            fields[i]*u[i] for i in range(n))+2*quadratic_u
        # Replay the product law itself, not just its quadratic expansion.
        direct_expectation = F(0)
        for y in spins:
            weight = F(1)
            for i in range(n):
                weight *= probabilities[i] if y[i] == -1 else 1-probabilities[i]
            direct_expectation += weight*sum(matrix[i][j]*y[i]*y[j]
                                            for i, j in edges)
        assert expected_energy == direct_expectation
        if d1:
            h_new_low = sum(entropy_interval(v)[0] for v in probabilities)
            h_base_high = n*entropy_interval(delta)[1]
            actual_gain_low = h_new_low-h_base_high+b*(expected_energy-q*q*cap)
            proposed_gain = n*min((beta*q*d1)**2/k_bound,
                                  beta*q*delta*d1/4)
            assert actual_gain_low >= proposed_gain > 0
            nonconstant += 1
        matrices += 1
    return {"all_order_four_signings": matrices,
            "nonconstant_field_cases_with_positive_exact_gain": nonconstant}


def hub_checks():
    outputs = []
    for r in (1, 3, 5, 7, 9):
        core_n, root_core = 4**r, 2**r
        total_n = core_n+1
        field_sum = core_n*(root_core+3)
        second = F(core_n**2+core_n*(root_core+2)**2, total_n**2)
        mean_squared = F(field_sum**2, total_n**3)
        variance = second-mean_squared
        assert variance >= 0
        # d1^2 is rational even though the normalized fields need not be.
        raw_mean = F(field_sum, total_n)
        raw_d1 = (abs(core_n-raw_mean)
                  +core_n*abs(root_core+2-raw_mean))/total_n
        d1_squared = raw_d1**2/total_n
        outputs.append({"r": r, "order": total_n,
                        "variance_exact": str(variance),
                        "variance_display": float(variance),
                        "d1_squared_exact": str(d1_squared),
                        "d1_squared_display": float(d1_squared)})
    assert outputs[-1]["variance_display"] > 0.99
    assert outputs[-1]["d1_squared_display"] < 0.00002
    return outputs


def compression_check():
    k = [[1-2*(i == j) for j in range(3)] for i in range(3)]
    square = [[sum(k[i][a]*k[a][j] for a in range(3))
               for j in range(3)] for i in range(3)]
    assert square == [[4*(i == j)-1 for j in range(3)] for i in range(3)]
    assert [sum(row) for row in k] == [1, 1, 1]
    cap = max(abs(sum(x[i]*x[j] for i in range(3) for j in range(i+1, 3)))
              for x in product((-1, 1), repeat=3))
    assert cap == 3
    # Exact strict entropy separation in the illustrative spectral net bound.
    log17_upper = log_interval(17)[1]
    h_upper = entropy_interval(F(1, 64))[1]
    log2_lower = log_interval(2)[0]
    assert log17_upper/31+h_upper < log2_lower
    return {"three_coordinate_hadamard_compression": "passed",
            "middle_boolean_eigenvalue": 1, "hollow_cap": cap,
            "spectral_net_entropy_upper": str(log17_upper/31+h_upper)}


def main():
    print(json.dumps({"adaptive": adaptive_checks(), "hub": hub_checks(),
                      "compression": compression_check()}, indent=2))


if __name__ == "__main__":
    main()
