"""Independent exact-integer replay of the finite clique inequalities.

This does not search for minimizers or numerically test the spectral-core
lemma. The latter has a reconstructed proof in the companion audit.
"""

import itertools

import numpy as np


rng = np.random.default_rng(2609061953)
cases = 0
shell_cases = 0
for n in range(3, 10):
    tails = np.array(list(itertools.product((-1, 1), repeat=n - 1)), dtype=np.int64)
    spins = np.column_stack((np.ones(len(tails), dtype=np.int64), tails))
    for trial in range(20):
        upper = np.triu(rng.choice((-1, 1), size=(n, n)), 1)
        matrix = upper + upper.T
        energy = np.einsum("bi,ij,bj->b", spins, matrix, spins) // 2
        witness = spins[np.argmax(np.abs(energy))]
        orientation = 1 if energy[np.argmax(np.abs(energy))] > 0 else -1
        matrix = orientation * matrix * np.outer(witness, witness)
        old = np.einsum("bi,ij,bj->b", spins, matrix, spins) // 2
        maximum = int(np.max(np.abs(old)))
        assert int(np.sum(matrix) // 2) == maximum
        for r in range(1, n // 2 + 1):
            subset = rng.choice(n, size=r, replace=False)
            block = np.zeros_like(matrix)
            block[np.ix_(subset, subset)] = matrix[np.ix_(subset, subset)]
            block_energy = np.einsum("bi,ij,bj->b", spins, block, spins) // 2
            u = int(np.max(np.abs(block_energy)))
            modified = matrix.copy()
            modified[np.ix_(subset, subset)] = 1
            np.fill_diagonal(modified, 0)
            values = np.einsum("bi,ij,bj->b", spins, modified, spins) // 2
            positive = int(np.max(values))
            negative = int(-np.min(values))
            cap = max(positive, negative)
            clique_cap = r * (r - 1) // 2
            assert positive >= maximum + clique_cap - u
            assert 2 * negative <= 2 * maximum + r + 2 * u
            assert 2 * (positive - negative) >= r * r - 2 * r - 4 * u
            assert maximum + clique_cap - u <= cap <= maximum + clique_cap + u
            assert not np.sum(values)  # Uniform signs are isotropic.
            assert np.all(
                (positive + negative) * np.abs(values)
                <= 2 * positive * negative + (positive - negative) * values
            )
            if positive > negative:
                shell = spins[values == positive]
                denominator = len(shell)
                gram_difference = shell.T @ shell - denominator * np.eye(n, dtype=np.int64)
                frobenius_squared_numerator = int(np.sum(gram_difference**2))
                assert (
                    frobenius_squared_numerator * n * (n - 1)
                    >= 4 * positive**2 * denominator**2
                )
                shell_cases += 1
            cases += 1

print(f"PASS: {cases} integer replacement/chord checks; {shell_cases} shell covariance checks.")
