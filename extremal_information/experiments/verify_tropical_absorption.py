#!/usr/bin/env python3
"""Exact finite checks for clamp absorption and selector-reset converse."""

from __future__ import annotations

import itertools
from collections import Counter
from fractions import Fraction


def maxplus_projective(matrix, z):
    """Column-output convention, with projective coordinate u_2-u_1=z."""
    u = (Fraction(0), Fraction(z))
    out = tuple(max(u[a] + matrix[a][b] for a in range(2)) for b in range(2))
    return out[1] - out[0]


def clip(z, lo, hi):
    return min(max(z, lo), hi)


def check_clamps():
    checks = 0
    for delta in (Fraction(1, 9), Fraction(1, 3), Fraction(3, 4)):
        s0 = ((0, 0), (-1, 0))
        sd = ((0, delta), (-1, 0))
        shat = ((0, 0), (-1 + delta, delta))
        for z in (Fraction(k, 12) for k in range(-36, 49)):
            p0 = maxplus_projective(s0, z)
            pd = maxplus_projective(sd, z)
            drift = maxplus_projective(shat, z)
            assert p0 == clip(z, 0, 1)
            assert pd == clip(z, delta, 1)
            assert drift == clip(z + delta, 0, 1)
            assert maxplus_projective(s0, p0) == p0
            assert maxplus_projective(sd, pd) == pd
            checks += 5

        # Half-oscillation Hilbert convention in two coordinates.
        assert abs(maxplus_projective(s0, 0) - maxplus_projective(sd, 0)) / 2 == delta / 2
        # The kernel difference has nonzero rectangular circulation.
        error = tuple(tuple(Fraction(sd[i][j]) - Fraction(s0[i][j]) for j in range(2)) for i in range(2))
        assert error[0][0] + error[1][1] - error[0][1] - error[1][0] == -delta
        # Coherent translation reaches macroscopic drift after ceil(1/delta) uses.
        z = Fraction(0)
        for _ in range((delta.denominator + delta.numerator - 1) // delta.numerator):
            z = maxplus_projective(shat, z)
        assert z > 0
    return checks


def apply_selector(vector, sigma):
    return tuple(vector[sigma[j]] for j in range(len(sigma)))


def product_map(word, r):
    out = tuple(range(r))
    for sigma in word:
        out = apply_selector(out, sigma)
    return out


def reset_free(word, r):
    for left in range(len(word)):
        for right in range(left + 1, len(word) + 1):
            if len(set(product_map(word[left:right], r))) == 1:
                return False
    return True


def adversarial_selector_errors(word, r):
    """Construct the errors in the proof of the selector lower bound."""
    suffix_maps = []
    for s in range(len(word)):
        mapping = product_map(word[s + 1 :], r)
        pair = next((j, k) for j in range(r) for k in range(r) if mapping[j] != mapping[k])
        suffix_maps.append((mapping, pair))
    chosen_pair, count = Counter(pair for _, pair in suffix_maps).most_common(1)[0]
    errors = []
    for mapping, pair in suffix_maps:
        eta = [0] * r
        if pair == chosen_pair:
            j, k = pair
            eta[mapping[j]] = 1
            eta[mapping[k]] = -1
        errors.append(tuple(eta))
    return errors, count


def check_selector_converse():
    checked = 0
    for r, max_depth in ((2, 5), (3, 3)):
        selectors = list(itertools.product(range(r), repeat=r))
        for depth in range(1, max_depth + 1):
            # Exhaust r=2.  For r=3 use a deterministic prefix of the much
            # larger word set, plus every permutation word.
            words = itertools.product(selectors, repeat=depth)
            limit = None if r == 2 else 3000
            for index, word in enumerate(words):
                if limit is not None and index >= limit:
                    break
                if not reset_free(word, r):
                    continue
                errors, repeated = adversarial_selector_errors(word, r)
                e = (0,) * r
                for sigma, eta in zip(word, errors):
                    e = tuple(x + y for x, y in zip(apply_selector(e, sigma), eta))
                hilbert = Fraction(max(e) - min(e), 2)
                assert hilbert >= depth // (r * (r - 1))
                assert repeated >= depth // (r * (r - 1))
                checked += 1

        for word in itertools.product(itertools.permutations(range(r)), repeat=max_depth):
            assert reset_free(word, r)
            errors, repeated = adversarial_selector_errors(word, r)
            assert repeated >= max_depth // (r * (r - 1))
            checked += 1
    return checked


def main():
    clamp_checks = check_clamps()
    selector_checks = check_selector_converse()
    print(f"exact clamp/idempotence checks: {clamp_checks}")
    print(f"reset-free selector lower-bound checks: {selector_checks}")


if __name__ == "__main__":
    main()
