#!/usr/bin/env python3
"""Create balanced search groups from the complete order-nine census."""
import json
from pathlib import Path
from twisted_chiral_driver_2026_09_18 import ensure_binaries

ROOT=Path(__file__).resolve().parents[1]
census=ROOT/'computations/results/twisted_chiral_2026_09_18_classify9.json'
classes=json.loads(census.read_text())['classes']
for name,indices in [('all',range(9)),('group1',[1,3,5,7]),('group2',[2,4,6,8]),('quick',[8,7,5,3,6,4])]:
    rows=[dict(label=f'm9_census{c["class"]}',matrix=c['matrix'],cap=12,
        source=str(census.relative_to(ROOT)),root_orbit_size=c['root_orbit_size']) for index in indices for c in classes if c['class']==index]
    (ROOT/f'computations/results/twisted_chiral_2026_09_18_census9_{name}_seeds.json').write_text(json.dumps(rows,indent=2)+'\n')
source=ROOT/'computations/results/twisted_chiral_2026_09_18_sampled_m10_seeds.json'
rows=json.loads(source.read_text())[1:]
(ROOT/'computations/results/twisted_chiral_2026_09_18_sampled_m10_new_seeds.json').write_text(json.dumps(rows,indent=2)+'\n')
ensure_binaries()
