"""Exact integer regression for the full one-hole clipping tradeoff."""
import json
from fractions import Fraction
from pathlib import Path
import numpy as np
from principle_director_stratified_selector_check_2026_09_07 import paley12, sylvester

records = []
for h in [sylvester(4), sylvester(8), paley12(), sylvester(16)]:
    n = len(h)
    labels = np.arange(1, (1 << n)-1, dtype=np.int64)
    words = ((labels[:,None] >> np.arange(n)) & 1).astype(np.int64)
    sizes = words.sum(axis=1)
    squares = (words @ h.T)**2
    for numerator in range(2, 2*n+1):
        c = Fraction(numerator, 2)
        a,b = c.numerator,c.denominator
        expected = min(1/c, Fraction(2,n), Fraction(n,1)/c-1)
        raw = np.minimum(b*squares, a).sum(axis=1)-a
        slack = raw*expected.denominator-sizes*a*expected.numerator
        assert np.all(slack >= 0)
        assert np.any(slack == 0)
        records.append(dict(order=n, c=str(c), minimum=str(expected),
                            packets=len(words), equality_packets=int((slack==0).sum())))
result = dict(status='PASS', checked_packet_temperature_pairs=sum(r['packets'] for r in records),
              records=records)
Path('computations/results/principle_construct_2026_09_07_packet_curve_check.json').write_text(
    json.dumps(result, indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k != 'records'}))
