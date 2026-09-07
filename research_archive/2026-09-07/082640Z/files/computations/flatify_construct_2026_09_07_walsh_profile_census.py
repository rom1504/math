"""Exact small Walsh magnitude-profile counts for the rank-two diagnostic."""
import collections
import json
from pathlib import Path
import numpy as np
from flatify_construct_2026_09_07_rank_two_weave import sylvester


records = []
for k in [4, 8, 16]:
    h = sylvester(k)
    x = 1-2*((np.arange(1 << k, dtype=np.int64)[:, None] >> np.arange(k)) & 1)
    spectra = np.abs(x @ h)
    profiles = collections.Counter(tuple(np.bincount(row, minlength=k+1)) for row in spectra)
    rows = []
    for profile, count in sorted(profiles.items(), key=lambda pair: -pair[1]):
        rows.append(dict(histogram=list(map(int, profile)), count=count,
                         entropy_per_coordinate=float(np.log(count)/k),
                         support=int(k-profile[0])))
    record = dict(k=k, profiles=len(profiles), total=sum(profiles.values()), rows=rows)
    records.append(record)
    print(json.dumps(dict(k=k, profiles=len(profiles), top=rows[:3])), flush=True)
Path('computations/results/flatify_construct_2026_09_07_walsh_profile_census.json').write_text(
    json.dumps(dict(status='exact finite counts, no asymptotic entropy assertion', records=records), indent=2)+'\n')
