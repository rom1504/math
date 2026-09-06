#!/usr/bin/env python3
"""Freeze the last pure finite-cell oracle as a compact rational policy.

Only monotonicity and rational grid membership are asserted here. The
separate exact certificate reconstructs every mathematical quantity and
does not trust optimizer scores or moments.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ascent",type=Path,default=Path("computations/results/resumed_response_rich_core_cell_ascent_2026_09_06.json"))
    parser.add_argument("--old-policy",type=Path,default=Path("computations/results/resumed_response_rich_core_rectangle_birth_policy_2026_09_06.json"))
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    data=json.loads(args.ascent.read_text())
    old=json.loads(args.old_policy.read_text())["first"]["rectangle_policy"]
    grid=[F(float(x)).limit_denominator(1000000) for x in data["z_endpoints"]]
    assert all(abs(float(x)-y)<1e-13 for x,y in zip(grid,data["z_endpoints"]))
    assert grid[0]==-8 and grid[-1]==8
    final=data["policies"][-1]
    if "fstar" in final:
        policy=final["fstar"]
    else:
        policy=[]
        for runs in final["row_runs"]:
            row=[]
            for end,value in runs:
                row.extend([value]*(end-len(row)))
            policy.append(row)
    assert len(policy)==len(old)
    rows=[]
    for row,prior in zip(policy,old):
        assert len(row)==len(grid)-1 and set(row)<={-1,0,1}
        assert all(a<=b for a,b in zip(row,row[1:]))
        assert row[-1]==1
        nneg=sum(x==-1 for x in row)
        nnonpos=sum(x<=0 for x in row)
        rows.append({"left":prior["left"],"right":prior["right"],
                     "lower_Z":str(grid[nneg]),"upper_Z":str(grid[nnonpos])})
    result={"status":"frozen_rational_policy_not_a_numeric_certificate",
            "source_ascent":str(args.ascent),"first":{"rectangle_policy":rows},
            "theta":"1","V_tail_policy":"sign(V) outside |V|>2",
            "negative_V_extension":"odd simultaneous (V,Z) reflection"}
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"rectangles":len(rows),"theta":"1","output":str(args.output)}))


if __name__=="__main__":
    main()
