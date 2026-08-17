#!/usr/bin/env python3
"""Toy checks for the dynamic broadcast-incidence decomposition."""

import numpy as np


def main():
    rho = 0.4
    P = np.array([[(1 + rho) / 2, (1 - rho) / 2],
                  [(1 - rho) / 2, (1 + rho) / 2]])
    one = np.ones(2)
    mode = np.array([1.0, -1.0])

    # Two visible self-loops: scalar means form two independently pumpable
    # cycle coordinates, while centred modes obey the geometric shell.
    means = [0.3, -0.2]
    centred = [0.7 * mode, -0.5 * mode]
    terminal = 0.9 * mode
    for letter in range(2):
        for depth in range(1, 20):
            d = np.linalg.matrix_power(P, depth) @ terminal
            for s in range(depth):
                d += np.linalg.matrix_power(P, s) @ centred[letter]
            d += depth * means[letter] * one
            scalar = depth * means[letter]
            shell = rho ** depth * 0.9 + abs(centred[letter][0]) * (
                1 - rho ** depth) / (1 - rho)
            assert np.linalg.norm(d - scalar * one) / np.sqrt(2) <= shell + 1e-12
            # Orthogonality makes the scalar cycle rate visible.
            assert np.linalg.norm(d) / np.sqrt(2) + 1e-12 >= abs(scalar)
    print("dynamic broadcast-incidence checks passed")


if __name__ == "__main__":
    main()
