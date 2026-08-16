#!/usr/bin/env python3
"""Exact checks for the weighted-automaton response benchmark.

No optimization solver or floating-point arithmetic is used.
"""

from itertools import product


BLOCKS = ((0, 1), (2, 3))
BLOCK_OF = {i: a for a, block in enumerate(BLOCKS) for i in block}

S = {
    "A": ((0, -2), (-1, 1)),
    "B": ((-1, 0), (2, -2)),
    # For u_empty=(0,0), row maxima are respectively (K,0) and (0,K).
    "Q0": ((4, 4), (0, 0)),
    "Q1": ((0, 0), (4, 4)),
}


def eta(i: int, j: int) -> int:
    """Microscopic row perturbation from Example WA.2."""

    position = BLOCKS[BLOCK_OF[j]].index(j)
    return ((0, -1) if i in (0, 2) else (-2, 0))[position]


def lift(small):
    return tuple(
        tuple(small[BLOCK_OF[i]][BLOCK_OF[j]] + eta(i, j) for j in range(4))
        for i in range(4)
    )


T = {letter: lift(small) for letter, small in S.items()}


def row_times_matrix(row, matrix):
    return tuple(max(row[i] + matrix[i][j] for i in range(len(row))) for j in range(len(matrix[0])))


def matrix_times_col(matrix, col):
    return tuple(max(matrix[i][j] + col[j] for j in range(len(col))) for i in range(len(matrix)))


def aggregate(row):
    return tuple(max(row[i] for i in block) for block in BLOCKS)


def suffix_vector(word, matrices, dimension):
    col = (0,) * dimension
    for letter in reversed(word):
        col = matrix_times_col(matrices[letter], col)
    return col


def check_lumpability():
    checks = 0
    for letter, matrix in T.items():
        for a, source_block in enumerate(BLOCKS):
            for b, target_block in enumerate(BLOCKS):
                for i in source_block:
                    value = max(matrix[i][j] for j in target_block)
                    assert value == S[letter][a][b]
                    checks += 1
    return checks


def check_forward_updates():
    checks = 0
    for row in product(range(-2, 3), repeat=4):
        coarse = aggregate(row)
        for letter in S:
            raw_next = aggregate(row_times_matrix(row, T[letter]))
            coarse_next = row_times_matrix(coarse, S[letter])
            assert raw_next == coarse_next
            checks += 1
    return checks


def words(alphabet, max_depth):
    yield ()
    for depth in range(1, max_depth + 1):
        yield from product(alphabet, repeat=depth)


def check_suffix_factorization_and_responses():
    factor_checks = 0
    response_checks = 0
    all_words = tuple(words(tuple(S), 3))
    for word in all_words:
        raw_suffix = suffix_vector(word, T, 4)
        coarse_suffix = suffix_vector(word, S, 2)
        assert raw_suffix == tuple(coarse_suffix[BLOCK_OF[i]] for i in range(4))
        factor_checks += 1

        for row in product(range(-1, 2), repeat=4):
            raw_value = max(row[i] + raw_suffix[i] for i in range(4))
            coarse = aggregate(row)
            coarse_value = max(coarse[a] + coarse_suffix[a] for a in range(2))
            assert raw_value == coarse_value
            response_checks += 1
    return factor_checks, response_checks


def check_pin_queries():
    checks = 0
    u0 = suffix_vector(("Q0",), S, 2)
    u1 = suffix_vector(("Q1",), S, 2)
    assert u0 == (4, 0)
    assert u1 == (0, 4)
    for p0, p1 in product(range(-2, 3), repeat=2):
        assert max(p0 + u0[0], p1 + u0[1]) == p0 + 4
        assert max(p0 + u1[0], p1 + u1[1]) == p1 + 4
        checks += 2
    return checks


def check_affine_line_falsifier():
    p = 8
    bound = 3
    constant = 2 * bound
    vectors = {
        t: tuple(constant * (2 * i * t - i * i) for i in range(1, p + 1))
        for t in range(1, p + 1)
    }

    # Collinearity: h(t)-h(1)=(t-1)(h(2)-h(1)).
    direction = tuple(vectors[2][i] - vectors[1][i] for i in range(p))
    checks = 0
    for t in range(1, p + 1):
        assert tuple(vectors[t][i] - vectors[1][i] for i in range(p)) == tuple(
            (t - 1) * direction[i] for i in range(p)
        )
        for i in range(1, p + 1):
            if i == t:
                continue
            assert vectors[t][t - 1] - vectors[t][i - 1] == constant * (t - i) ** 2
            assert vectors[t][t - 1] - vectors[t][i - 1] >= 2 * bound
            checks += 1
    return checks


def main():
    lump = check_lumpability()
    updates = check_forward_updates()
    factors, responses = check_suffix_factorization_and_responses()
    pins = check_pin_queries()
    affine = check_affine_line_falsifier()
    total = lump + updates + factors + responses + pins + affine
    print(f"lumpability identities: {lump}")
    print(f"forward aggregation updates: {updates}")
    print(f"suffix factorizations: {factors}")
    print(f"raw/quotient responses: {responses}")
    print(f"pin exposures: {pins}")
    print(f"affine-line identities/exposures: {affine}")
    print(f"all {total} exact weighted-automaton checks passed")


if __name__ == "__main__":
    main()
