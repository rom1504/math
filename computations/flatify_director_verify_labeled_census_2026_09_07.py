"""Remove reliance on graph-atlas completeness for the minimizing packing census.

The C++ helper independently exhausts labeled switched signings. This checker
covers its complete minimizer list by explicit residual vertex permutations of
the stored rational-LP representatives. No floating-point conclusion is used.
"""
import hashlib
import itertools
import json
from fractions import Fraction as F
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'flatify_director_labeled_cap_census_2026_09_07'
source = ROOT / 'computations' / (PREFIX + '.cpp')
binary = ROOT / 'tmp' / PREFIX
subprocess.run(['g++', '-O3', '-std=c++17', str(source), '-o', str(binary)], check=True)
data = json.loads((ROOT / 'computations/results/flatify_director_small_packing_census_2026_09_07.json').read_text())
reports = []
for record in data['records']:
    n = record['n']
    brute = json.loads(subprocess.check_output([str(binary), str(n)]))
    assert brute['M'] == record['M']
    covered = set()
    edges = list(itertools.combinations(range(1, n), 2))
    for row in record['rows']:
        if row['Q'] != brute['M']:
            continue
        a = row['matrix']
        for perm0 in itertools.permutations(range(1, n)):
            perm = (0,) + perm0
            covered.add(sum((a[perm[i]][perm[j]] < 0) << e
                            for e, (i, j) in enumerate(edges)))
        for sign, cert in zip([1, -1], row['packing']):
            d = list(map(F, cert['dual']))
            assert all(v >= 0 for v in d)
            energies = {}
            for z in itertools.product([-1, 0, 1], repeat=n):
                support = sum((z[i] != 0) << i for i in range(n))
                h = sign * sum(a[i][j]*z[i]*z[j]
                               for i in range(n) for j in range(i+1, n))
                assert h <= sum(d[i]*z[i]**2 for i in range(n))
                energies[support] = max(energies.get(support, -n*n), h)
            load = [F(0)]*n
            value = F(0)
            for atom in cert['primal']:
                s, w = atom['support'], F(atom['weight'])
                assert w >= 0 and energies[s] == atom['energy']
                value += w*energies[s]
                for i in range(n):
                    load[i] += w*((s >> i) & 1)
            assert max(load) <= 1 and value == sum(d) == F(cert['value'])
    assert covered == set(brute['minimizers'])
    lower = F(record['minimum_packing_among_cap_minimizers'])
    next_caps = [int(k) for k in brute['histogram'] if int(k) > brute['M']]
    # For every signing, ternary packing dominates its full-spin cap.
    # Consequently all nonminimizers are excluded without their LPs.
    assert not next_caps or lower <= min(next_caps)
    reports.append(dict(n=n, M=brute['M'], minimum_packing=str(lower),
                        labeled_signings=brute['signings'],
                        labeled_minimizers=len(covered),
                        cap_histogram=brute['histogram'],
                        complete_minimizer_cover=True,
                        rational_certificates_rechecked=True))
    print(json.dumps(reports[-1]), flush=True)
out = dict(status='EXACT FINITE CERTIFICATION; no asymptotic inference',
           cpp_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
           checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           reports=reports)
(ROOT / 'computations/results/flatify_director_labeled_cap_census_2026_09_07.json').write_text(json.dumps(out, indent=2)+'\n')
