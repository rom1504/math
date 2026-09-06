"""Exact 128-bin certificate for a uniform purified full-center gain."""

import argparse
import json
from fractions import Fraction as F
from pathlib import Path

from fresh_limit_mask_ascent_certificate import cdf_point, phi_large
from fresh_limit_rooted_lower_certificate import DIGITS, I, Phi, phi


def point(alpha):
    mass = 2*Phi(alpha)-1
    projection_cap = 2*phi(alpha)
    old_value_cap = projection_cap*mass.sqrt()
    derivative_sign = projection_cap-2*alpha*mass
    return mass, projection_cap, old_value_cap, derivative_sign


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--bins", type=int, default=128)
    args = parser.parse_args()
    assert args.bins >= 1
    left, right = F(47, 100), F(22, 25)
    threshold, target = F(43, 100), F(1, 40000)
    left_data, right_data = point(left), point(right)
    assert left_data[2].hi < threshold and right_data[2].hi < threshold
    assert left_data[3].lo > 0 and right_data[3].hi < 0
    width = (right-left)/args.bins
    records = []
    minimum_lower = None
    minimum_index = None
    for index in range(args.bins):
        a = left+index*width
        b = a+width
        mass_lower, projection_upper, _, _ = point(a)
        mass_upper, _, _, _ = point(b)
        support_lower = 1-mass_upper
        variance_lower = support_lower-projection_upper*projection_upper
        assert variance_lower.lo > 0
        sd_lower = variance_lower.sqrt()
        center_upper = projection_upper/mass_lower.sqrt()
        argument = center_upper/sd_lower
        assert 0 < argument.lo <= argument.hi < 8
        gain_lower = 2*mass_lower*(sd_lower*phi_large(argument)
                                  -center_upper*(1-cdf_point(argument)))
        assert gain_lower.lo > target
        if minimum_lower is None or gain_lower.lo < minimum_lower:
            minimum_lower, minimum_index = gain_lower.lo, index
        records.append({"left": str(a), "right": str(b),
                        "gain_lower_bound": gain_lower.json()})
    result = {
        "method": "exact_fraction_outward_intervals_and_monotone_Darboux_bounds",
        "digits": DIGITS,
        "old_value_threshold": str(threshold),
        "threshold_parameter_interval": [str(left), str(right)],
        "left_old_value_cap": left_data[2].json(),
        "right_old_value_cap": right_data[2].json(),
        "left_derivative_sign": left_data[3].json(),
        "right_derivative_sign": right_data[3].json(),
        "bins": args.bins,
        "minimum_bin_index": minimum_index,
        "minimum_bin_lower_bound": I(minimum_lower).json(),
        "certified_uniform_gain": str(target),
        "bin_certificates": records,
        "verified": True,
    }
    rendered = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
        print(json.dumps({key: value for key, value in result.items()
                          if key != "bin_certificates"}, indent=2))
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
