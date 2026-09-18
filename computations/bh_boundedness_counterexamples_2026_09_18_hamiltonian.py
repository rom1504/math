#!/usr/bin/env python3
"""Finite diagnostic replay of exact identities in the phase-entropy proof.

The theorem is analytic. Numerical entropy values are diagnostics only.
"""
import json
import math

import numpy as np


def fwht(values):
    values = np.array(values, copy=True)
    width = 1
    while width < len(values):
        for start in range(0, len(values), 2*width):
            left = values[start:start+width].copy()
            right = values[start+width:start+2*width].copy()
            values[start:start+width] = left+right
            values[start+width:start+2*width] = left-right
        width *= 2
    return values


def main():
    rng = np.random.RandomState(20260918)
    cases = 0
    max_ratio = 0.0
    operator_cases = 0
    for n in range(2, 9):
        words = np.arange(2**n)
        signs = 1-2*((words[:, None] >> np.arange(n)) & 1)
        for trial in range(12):
            matrix = rng.randint(-3, 4, size=(n, n))
            matrix = np.triu(matrix, 1)
            matrix = matrix+matrix.T
            energy = np.sum((signs @ matrix)*signs, axis=1)/2
            cap = np.max(np.abs(energy))
            if not cap:
                continue
            if n <= 5 and trial < 3:
                J = np.zeros((2**n, 2**n))
                for i in range(n):
                    J[words, words ^ (1 << i)] = 0.5
                C = 1j*J*(energy[None, :]-energy[:, None])
                Cnorm = max(abs(np.linalg.eigvalsh(C)))
                assert Cnorm+1e-10 >= cap
                for amplitude in [0.1, 0.4]:
                    phase_values = np.exp(1j*amplitude*energy/cap)
                    velocity_matrix = J*(phase_values.conj()[:, None]*phase_values[None, :]-1)
                    velocity = max(abs(np.linalg.eigvalsh(velocity_matrix)))
                    assert velocity+1e-10 >= amplitude*(1-amplitude)
                    operator_cases += 1
            for phase in [0.01, 0.1, 1, 4]:
                t = phase/cap
                values = np.exp(1j*t*energy)
                coeff = fwht(values)/len(words)
                probabilities = np.abs(coeff)**2
                entropy = -sum(p*math.log(p) for p in probabilities if p > 0)
                assert abs(sum(probabilities)-1) < 1e-12
                for i in range(n):
                    marginal = probabilities[((words >> i) & 1).astype(bool)].sum()
                    identity = np.mean(np.sin(t*(signs @ matrix[:, i]))**2)
                    assert abs(marginal-identity) < 1e-12
                row_sum = np.linalg.norm(matrix, axis=1).sum()
                assert row_sum <= 4*math.sqrt(3)*cap+1e-12
                assert entropy <= 2*abs(t)*row_sum+1e-12
                assert entropy <= 8*math.sqrt(3)*phase+1e-12
                max_ratio = max(max_ratio, entropy/phase)
                cases += 1
    print(json.dumps({"status": "PASS", "phase_cases": cases,
                      "small_amplitude_operator_cases": operator_cases,
                      "maximum_entropy_over_abs_t_cap_diagnostic": max_ratio,
                      "proved_constant": "8 sqrt(3)"}, indent=2))


if __name__ == "__main__":
    main()
