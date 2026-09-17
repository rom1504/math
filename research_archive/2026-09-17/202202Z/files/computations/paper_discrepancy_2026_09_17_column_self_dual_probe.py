"""Try a rational primal law as a code-side dual; no optimization."""
from fractions import Fraction as F
from pathlib import Path
import itertools
import json
import math
import numpy as np
import argparse

ROOT=Path(__file__).resolve().parents[1]
data=json.loads((ROOT/"tmp/paper_portfolio_2026_09_17/discrepancy/full_column_game_certificates.json").read_text())
parser=argparse.ArgumentParser()
parser.add_argument("--full",action="store_true")
args=parser.parse_args()
certificates={}
for row in data["cases"]:
    if not (row["n"]==14 and not row["isotropy_required"] and row["deficit_window"] in (2,4)):
        continue
    n=row["n"]
    words=np.asarray([(1,)+h for h in itertools.product((-1,1),repeat=n-1)],dtype=np.int64)
    atoms=row["primal_law"]
    support=np.asarray([a["word"] for a in atoms],dtype=np.int64)
    assert set(map(tuple,support)) <= {tuple(a["word"]) for a in row["dual_code_law"]}
    weights=[F(a["weight"]) for a in atoms]
    denominator=math.lcm(*(w.denominator for w in weights))
    numerators=np.asarray([w.numerator*(denominator//w.denominator) for w in weights],dtype=object)
    response=np.abs(words@support.T).astype(object)@numerators
    minimum=F(int(min(response)),denominator)
    if not args.full:
        print(row["case"],row["deficit_window"],"minimum",minimum,"target",row["upper_exact"])
    if minimum==F(row["upper_exact"]):
        row["lower_exact"]=str(minimum)
        row["dual_equality_coefficients"]=[str(minimum)]
        row["dual_code_law"]=atoms
        row["dual_code_support_size"]=len(atoms)
        row["endpoint_recovery_before_self_dual"]=row["status"]
        row["status"]="exact rational primal is also an all-column code-side dual"
        certificates[(row["case"],row["deficit_window"])]=row
for row in data["cases"]:
    key=(row["case"],row["deficit_window"])
    if row["isotropy_required"] and key in certificates:
        certificate=certificates[key]
        assert row["upper_exact"]==certificate["upper_exact"]
        row["lower_exact"]=certificate["lower_exact"]
        row["dual_equality_coefficients"]=certificate["dual_equality_coefficients"]+["0"]*(row["n"]*(row["n"]-1)//2)
        row["dual_code_law"]=certificate["dual_code_law"]
        row["dual_code_support_size"]=certificate["dual_code_support_size"]
        row["endpoint_recovery_before_transfer"]=row["status"]
        row["status"]="exact rational primal and transferred unrestricted all-column dual agree"
if args.full:
    print(json.dumps(data,separators=(",",":")))
