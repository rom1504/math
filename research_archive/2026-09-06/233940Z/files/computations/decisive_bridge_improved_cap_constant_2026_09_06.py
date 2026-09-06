"""Outward rational conversion of the replayed E certificate to cap."""
import json
from fractions import Fraction as F
from pathlib import Path
from continued_feedback_ternary_latent_interval_certificate_2026_09_06 import sqrt_interval


def main():
    p=F(31,32)
    a=F(19678127864847,800000000000000)
    root_lower,root_upper=sqrt_interval(p)
    lower=F(1,2)-a/(8*root_lower)
    upper=F(1,2)-a/(8*root_upper)
    assert upper<F(496876095,10**9)
    result={"status":"exact outward rational arithmetic",
            "p":str(p),"a_E":str(a),"sqrt_lower":str(root_lower),"sqrt_upper":str(root_upper),
            "cap_lower":str(lower),"cap_upper":str(upper),
            "cap_lower_display":float(lower),"cap_upper_display":float(upper),
            "upper_decimal":"0.496876095",
            "scope":"constant evaluation only; all-order proof in companion artifact"}
    Path("computations/decisive_bridge_improved_cap_constant_2026_09_06.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":main()
