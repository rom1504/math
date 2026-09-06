#!/usr/bin/env python3
"""Finite 8-block recovery entropy, with an optional 4-wise marginal LP.

Floating LP values are discovery only. A dual must be enclosed rigorously
before it is used as a mathematical upper bound.
"""

import argparse
import itertools
import json
import math
from pathlib import Path
from fractions import Fraction
from functools import lru_cache

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix


# A fixed discovery dual, rounded outward and checked using rational arithmetic.
# Row order is exactly marginal_matrix's degree/coordinate/symbol order.
CERTIFIED_DUAL = {
    0:-39810963604,1:-53047564697,3:1373970815674,7:902270215988,
    17:-955132548179,18:-745336112683,45:-1619996279256,46:-1645081806875,
    49:-4174886125730,50:-3668107923281,93:-1690833592983,129:1128267212379,
    130:2204434303748,132:930535412337,137:883737608001,138:788925177859,
    139:518926783928,140:758272401849,225:103232190676,297:1326153313766,
    298:1128421513724,299:2550035451282,305:275639366051,306:877089153392,
    313:3296315662224,315:1689123455341,316:1416270959298,345:2490310908423,
    346:4894089781372,348:2186805636370,577:-452495662678,578:-429475107560,
    579:-858950215119,580:-897540908565,583:-338388484570,584:-429475107560,
    593:-203681327392,594:-695211799135,596:-667771494383,600:-171328460485,
    609:-540906049872,611:-510253273861,612:-388266092037,615:-357613316027,
    673:10802731814,674:19518279229,676:27440304752,677:107861332536,
    678:74398888733,1025:-169032721641,1026:-173938966552,1137:-27440304752,
    1138:6556456967,1140:103615057688,1142:-54880609504,1153:-2283355738192,
    1155:-2100063004348,1157:-2191397410162,1158:-2192021332378,
    1233:-169032721641,1234:-43484741638,1235:-130454224914,
    1249:-774929804525,1250:-736339111080,1253:-736339111080,
    1254:-391190110658,1255:-277082932550,1345:7698035469,
    1346:36819698438,1347:35138340221,1351:69172565345,
}


def data():
    patterns = np.asarray(list(itertools.product((-1, 0, 1), repeat=8)), dtype=np.int8)
    hadamard = np.asarray([[1 if bin(i & j).count('1') % 2 == 0 else -1
                           for j in range(8)] for i in range(8)], dtype=np.int8)
    transformed = patterns @ hadamard
    masses = np.prod(np.where(patterns == 0, 2, 15).astype(np.int64), axis=1)
    probabilities = masses / float(32**8)
    return patterns, transformed, masses, probabilities


def marginal_matrix(patterns):
    rows, columns, values = [], [], []
    rhs = []
    row = 0
    for degree in range(5):
        for coordinates in itertools.combinations(range(8), degree):
            for symbols in itertools.product((-1, 1), repeat=degree):
                kept = np.ones(len(patterns), dtype=bool)
                for coordinate, symbol in zip(coordinates, symbols):
                    kept &= patterns[:, coordinate] == symbol
                indices = np.flatnonzero(kept)
                rows.extend([row] * len(indices))
                columns.extend(indices.tolist())
                values.extend([1.0] * len(indices))
                rhs.append((15/32)**degree)
                row += 1
    return csr_matrix((values, (rows, columns)), shape=(row, len(patterns))), np.asarray(rhs)


def symmetry_compression(patterns, selected):
    selected = set(selected)
    parent = np.arange(len(patterns))
    powers = 3**np.arange(7, -1, -1)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    actions = 0
    for c1 in range(1, 8):
        for c2 in range(1, 8):
            if c2 == c1:
                continue
            for c4 in range(1, 8):
                if c4 in (c1, c2, c1 ^ c2):
                    continue
                linear = [((c1 if f & 1 else 0) ^ (c2 if f & 2 else 0)
                           ^ (c4 if f & 4 else 0)) for f in range(8)]
                permutation = [sum((bin(c & x).count('1') % 2) * bit
                                   for c, bit in ((c1,1),(c2,2),(c4,4)))
                               for x in range(8)]
                for offset in range(8):
                    if {linear[f] ^ offset for f in selected} != selected:
                        continue
                    signs = np.asarray([(-1)**(bin(offset & x).count('1'))
                                        for x in range(8)], dtype=np.int8)
                    for global_sign in (-1, 1):
                        transformed = patterns[:, permutation] * signs * global_sign
                        indices = (transformed.astype(int)+1) @ powers
                        for a, b in enumerate(indices):
                            ra, rb = find(a), find(int(b))
                            if ra != rb:
                                parent[rb] = ra
                        actions += 1
    roots = np.asarray([find(a) for a in range(len(patterns))])
    _, inverse = np.unique(roots, return_inverse=True)
    sizes = np.bincount(inverse)
    compression = csr_matrix((1/sizes[inverse], (np.arange(len(patterns)), inverse)),
                             shape=(len(patterns), len(sizes)))
    return compression, inverse, actions


@lru_cache(None)
def logarithm_bounds(numerator, denominator=1, terms=35):
    """Exact rational log enclosure via 2*atanh and range reduction."""
    ratio = Fraction(numerator, denominator)
    power = numerator.bit_length() - denominator.bit_length()
    scaled = ratio / (Fraction(2)**power)
    if scaled < 1:
        power -= 1
        scaled *= 2
    if scaled >= 2:
        power += 1
        scaled /= 2

    def series(value):
        z = (value-1)/(value+1)
        total, term = Fraction(0), z
        for j in range(terms):
            total += 2*term/(2*j+1)
            term *= z*z
        remainder = 2*term/((2*terms+1)*(1-z*z))
        return total, total+remainder

    lower, upper = series(scaled)
    two_lower, two_upper = series(Fraction(2))
    if power >= 0:
        return lower+power*two_lower, upper+power*two_upper
    return lower+power*two_upper, upper+power*two_lower


def certify_dual(matrix, rhs, orbit, masses, inverse_output, dual=None, scale=10**10):
    sizes = np.bincount(orbit)
    indicator = csr_matrix((np.ones(len(orbit), dtype=np.int64),
                            (np.arange(len(orbit)), orbit)),
                           shape=(len(orbit), len(sizes)))
    counts = (matrix.astype(np.int64) @ indicator).tocsr()
    if dual is None:
        integers = np.zeros(matrix.shape[0], dtype=np.int64)
        for j, value in CERTIFIED_DUAL.items():
            integers[j] = value
    else:
        integers = np.rint(-dual*scale).astype(np.int64)
    representatives = np.asarray([np.flatnonzero(orbit == j)[0]
                                  for j in range(len(sizes))])
    output_mass = np.zeros(int(inverse_output.max())+1, dtype=np.int64)
    np.add.at(output_mass, inverse_output, masses)
    assert np.array_equal(output_mass[inverse_output],
                          output_mass[inverse_output[representatives]][orbit])
    maximum_sum = int(np.max(np.asarray(counts.sum(axis=0))))
    assert max(abs(int(v)) for v in integers)*maximum_sum < 2**62
    numerators = counts.T @ integers
    required = 0
    for j, representative in enumerate(representatives):
        _, upper = logarithm_bounds(int(output_mass[inverse_output[representative]]), 32**8)
        polynomial = Fraction(int(numerators[j]), scale*int(sizes[j]))
        gap = (upper-polynomial)*scale
        required = max(required, -((-gap.numerator)//gap.denominator))
    integers[0] += required
    numerators = counts.T @ integers
    for j, representative in enumerate(representatives):
        _, upper = logarithm_bounds(int(output_mass[inverse_output[representative]]), 32**8)
        assert Fraction(int(numerators[j]), scale*int(sizes[j])) >= upper
    objective = sum(Fraction(int(value))*Fraction(float(target))
                    for value, target in zip(integers, rhs))/scale
    two_lower, two_upper = logarithm_bounds(2)
    count_upper = Fraction(23,16)*two_upper+objective/8
    sqrt15_lower = Fraction(3872983346207416, 10**15)
    assert sqrt15_lower**2 < 15
    exponential_partial = sum(Fraction(16)**j/Fraction(math.factorial(j)) for j in range(50))
    assert exponential_partial > 8000000
    tilted_upper = Fraction(15,16)*two_upper+objective/8+Fraction(1,16000000)+4-sqrt15_lower
    return {'certificate':'exact rational averaged dual, rational atanh log upper bounds',
            'dual_scale':scale, 'constant_correction_units':required,
            'objective_upper_fraction':[str(objective.numerator),str(objective.denominator)],
            'objective_upper_decimal':float(objective), 'count_rate_upper':float(count_upper),
            't4_tilted_flat_upper':float(tilted_upper),
            'nonzero_dual_entries':[[int(j),int(value)] for j,value in enumerate(integers) if value]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--lp', action='store_true')
    parser.add_argument('--symmetry', action='store_true')
    parser.add_argument('--certify', action='store_true')
    parser.add_argument('--verify-only', action='store_true')
    parser.add_argument('--output', type=str)
    parser.add_argument('--only', type=str, nargs='*')
    args = parser.parse_args()
    patterns, transformed, masses, probabilities = data()
    choices = {
        'l1':[0], 'l2':[0,1], 'l3':[0,1,2],
        'l4_plane':[0,1,2,3], 'l4_nonplane':[0,1,2,4],
        'l5':[0,1,2,3,4], 'l6':[0,1,2,3,4,5],
        'l7':[0,1,2,3,4,5,6], 'l8':list(range(8)),
    }
    matrix, rhs = marginal_matrix(patterns) if args.lp or args.verify_only else (None, None)
    p = 15/16
    records = []
    for name, selected in choices.items():
        if args.verify_only and name != 'l4_nonplane':
            continue
        if args.only and name not in args.only:
            continue
        output, inverse = np.unique(transformed[:, selected], axis=0, return_inverse=True)
        output_mass = np.zeros(len(output), dtype=np.int64)
        np.add.at(output_mass, inverse, masses)
        log_probability = np.log(output_mass[inverse].astype(float)) - 8*np.log(32)
        iid_objective = float(probabilities @ log_probability)
        constant = len(selected)/8*np.log(2)+p*np.log(2)
        result = {'status':'finite entropy / floating LP discovery; not a rigorous dual',
                  'projection':name, 'selected':selected, 'distinct_outputs':len(output),
                  'iid_E_log_PY':iid_objective, 'iid_count_rate':constant+iid_objective/8}
        if args.verify_only:
            compression, orbit, actions = symmetry_compression(patterns, selected)
            result['status'] = 'fixed rational dual verified without calling an optimizer'
            result['rigorous_dual'] = certify_dual(matrix, rhs, orbit, masses, inverse)
        if args.lp:
            if args.symmetry:
                compression, orbit, actions = symmetry_compression(patterns, selected)
                solve_matrix = matrix @ compression
                solve_cost = np.asarray(compression.T @ log_probability)
                result.update({'symmetry_actions':actions, 'symmetry_orbits':compression.shape[1],
                               'cost_symmetry_error':float(np.max(np.abs(log_probability-solve_cost[orbit])))})
            else:
                solve_matrix, solve_cost = matrix, log_probability
            solved = linprog(-solve_cost, A_eq=solve_matrix, b_eq=rhs,
                             bounds=(0,None), method='highs')
            result.update({'lp_success':solved.success, 'lp_message':solved.message})
            if solved.success:
                maximum = float(solve_cost @ solved.x)
                result.update({'lp_E_log_PY':maximum, 'lp_count_rate':constant+maximum/8,
                               'marginal_residual':float(np.max(np.abs(solve_matrix@solved.x-rhs))),
                               'nonzero_variables':int(np.count_nonzero(solved.x>1e-10)),
                               'dual_residual':float(np.max(solve_matrix.T @ solved.eqlin.marginals
                                                            + solve_cost))})
                if args.certify:
                    assert args.symmetry and name == 'l4_nonplane'
                    result['rigorous_dual'] = certify_dual(matrix, rhs, orbit, masses, inverse,
                                                           solved.eqlin.marginals)
        records.append(result)
        print(json.dumps(result), flush=True)
    if args.output:
        destination = Path(args.output)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(records, indent=2)+"\n")


if __name__ == '__main__':
    main()
