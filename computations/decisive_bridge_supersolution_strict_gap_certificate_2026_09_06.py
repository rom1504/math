"""Exact rational-log certificate for a finite source with E<T."""
import json
from fractions import Fraction as F
from pathlib import Path
from continued_feedback_ternary_latent_interval_certificate_2026_09_06 import log_interval, entropy_interval


def main():
    radius = 2048
    m = 2 * radius + 1
    variance = F(radius * (radius + 1), 3)
    t_lower = -entropy_interval(F(1,4))[1] - F(1,8)
    t_lower -= log_interval(variance)[1] / 16 + log_interval(F(4,3))[1] / 16
    # Classical Archimedean bound pi<22/7 is sufficient here.
    e_upper = -log_interval(F(m))[0] / 4 + log_interval(F(2))[1] / 4
    e_upper += log_interval(F(22,7))[1] / 8 + F(1,48)
    gap = t_lower - e_upper
    assert gap > F(3,20)
    result = {"status":"exact rational interval certificate",
              "source":"(3/4)delta_0+(1/4)Uniform{-2048,...,2048}",
              "temperature":1,"uniform_cardinality":m,"uniform_variance":str(variance),
              "T_lower":str(t_lower),"T_lower_display":float(t_lower),
              "E_upper":str(e_upper),"E_upper_display":float(e_upper),
              "gap_lower":str(gap),"gap_lower_display":float(gap),
              "gap_exceeds":"3/20",
              "scope":"E<T certificate; strict BT<T additionally uses analytic equality theorem"}
    Path("computations/decisive_bridge_supersolution_strict_gap_certificate_2026_09_06.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
