#!/usr/bin/env python3
"""Concrete nonoptimal children retained for the twisted-double campaign."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
source=ROOT/'computations/results/conference_order10_gf9.json'
data=json.loads(source.read_text())
rows=[dict(label='m10_conference',matrix=data['conference_matrix'],source=str(source.relative_to(ROOT)),
    evidence='nonoptimal conference child; exact child cap independently checked by driver')]
(ROOT/'computations/results/twisted_chiral_2026_09_18_alternative_conference_seed.json').write_text(json.dumps(rows,indent=2)+'\n')
