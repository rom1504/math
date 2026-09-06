"""Exact-integer finite Hadamard replay of fixed mixed-charge feedback.

All matrix-vector products, inertia comparisons, and final energies are exact
integer operations. Random seeds choose examples, not interval certification
of an asymptotic Gaussian expectation or a universal original-class bound.
"""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np


def hadamard(vector):
    out = vector.copy()
    width = 1
    while width < len(out):
        rows = out.reshape(-1, 2*width)
        left = rows[:, :width].copy()
        right = rows[:, width:].copy()
        rows[:, :width] = left+right
        rows[:, width:] = left-right
        width *= 2
    return out


def run(power, steps, seed):
    if power < 2 or power > 30 or power % 2:
        raise ValueError('use a positive even power <=30 for rational, int64-safe arithmetic')
    if steps < 0:
        raise ValueError('steps must be nonnegative')
    n = 2**power
    sqrt_n = math.isqrt(n)
    spin = np.random.default_rng(seed).choice(np.array([-1, 1], dtype=np.int64), size=n)
    schedule = [Fraction(1, 2)]*4+[Fraction(3, 10)]*4+[Fraction(1, 5)]*8
    schedule += [Fraction(1, 10)]*max(0, steps-len(schedule))
    trace = []
    for step in range(steps+1):
        field = hadamard(spin)
        numerator = int(np.dot(spin, field))
        energy = Fraction(numerator, 2*n*sqrt_n)
        trace.append(dict(step=step, energy_exact=str(energy), energy=float(energy),
                          quadratic_numerator=numerator))
        if step == steps:
            break
        alpha = schedule[step]
        score = alpha.denominator*field+alpha.numerator*sqrt_n*spin
        # Equivariant Boolean convention at an exact zero score.
        spin = np.where(score == 0, spin, np.sign(score))
    return dict(status='exact finite matrix calculation, not an asymptotic certificate',
                matrix='normalized symmetric Sylvester Hadamard',
                power=power, n=n, steps=steps, seed=seed,
                schedule=[str(x) for x in schedule], trace=trace)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--power', type=int, default=18)
    parser.add_argument('--steps', type=int, default=40)
    parser.add_argument('--seed', type=int, default=617)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run(args.power, args.steps, args.seed)
    encoded = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(encoded)
    else:
        print(encoded, end='')


if __name__ == '__main__':
    main()
