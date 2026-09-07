"""Directed-interval certificate for the paired-row anti-weave phase.

All accepted boxes are outward-rounded. Density interpolation is justified
by convexity of E at fixed source amplitude; no floating optimization is
used to accept a box. Exploratory targets are rational inputs, not bounds
until every corresponding rectangle has passed.
"""
import argparse
from fractions import Fraction as F
import json
import time

import mpmath as mp

from flatify_independent_2026_09_07_ternary_interval import (
    Certificate, endpoints, rational_iv, range_iv)


P = F(24, 25)
T = F(7, 2)
AMP2 = F(25, 12)
LOW_PRECISION = F(6, 25)
# Numeric pilots selected these rational targets; the verifier must prove
# each one afresh. The extra 1/2000 keeps a safe strict certificate margin.
TARGET_MILLI = [0, -224, -378, -465, -529, -581, -623, -660, -691,
                -720, -745, -768, -789, -808, -826, -843, -859,
                -874, -888, -901, -832]
TARGETS = [F(0)] + [F(x, 1000) + F(1, 2000) for x in TARGET_MILLI[1:]]


class FixedAmplitudeCertificate(Certificate):
    def __init__(self, index, precision=40):
        self.index = index
        density = P * F(index, 20)
        super().__init__(density, T, TARGETS[index], precision)
        self.amplitude_sq = rational_iv(AMP2)
        self.variance = density * AMP2
        self.vv = rational_iv(self.variance)

    def gain_arguments(self, lam, z):
        a2 = self.amplitude_sq
        A = mp.iv.exp(a2 * lam * z * z) - 1
        h = a2 * lam * z
        B = (mp.iv.exp(h) - mp.iv.exp(-h)) ** 2 / 2
        return A, B

    def point_upper(self, lam, z):
        ll, zz = rational_iv(lam), rational_iv(z)
        A, B = self.gain_arguments(ll, zz)
        cminus = mp.iv.ln(ll * (2 * self.tt - ll) / self.tt ** 2) / 4 - ll * self.vv
        return cminus + self.optimized_gain(A, B)

    def monotone_upper(self, rect):
        l, u, a, b = rect
        ll, uu, aa, bb = map(rational_iv, rect)
        A = mp.iv.exp(self.amplitude_sq * ll * aa * aa) - 1
        h = self.amplitude_sq * uu * bb
        B = (mp.iv.exp(h) - mp.iv.exp(-h)) ** 2 / 2
        # c_t increases on (0,t], while -lambda*v decreases.
        cminus = mp.iv.ln(uu * (2 * self.tt - uu) / self.tt ** 2) / 4 - ll * self.vv
        return endpoints(cminus + self.optimized_gain(A, B))[1]

    def meanvalue_upper(self, rect):
        l, u, a, b = rect
        if a == 0:
            return None
        lam, z = range_iv(l, u), range_iv(a, b)
        A, B = self.gain_arguments(lam, z)
        if endpoints(A)[0] <= 0 or endpoints(B)[0] <= 0:
            return None
        unclipped = (self.pp * B - A) / ((1 - self.pp) * A * B)
        vlo, vhi = endpoints(unclipped)
        vlo = max(F(0), min(F(1), vlo))
        vhi = max(F(0), min(F(1), vhi))
        weight = range_iv(vlo, vhi)
        expA = mp.iv.exp(self.amplitude_sq * lam * z * z)
        q = 2 * self.amplitude_sq * lam * z
        sinhq = (mp.iv.exp(q) - mp.iv.exp(-q)) / 2
        dl = (1 / lam - 1 / (2 * self.tt - lam)) / 4 - self.vv
        dl += (2 * self.pp * weight * self.amplitude_sq * z * sinhq / (1 + weight * B)
               - weight * self.amplitude_sq * z * z * expA / (1 + weight * A))
        dz = (2 * self.pp * weight * self.amplitude_sq * lam * sinhq / (1 + weight * B)
              - 2 * weight * self.amplitude_sq * lam * z * expA / (1 + weight * A))
        dlmax = max(map(abs, endpoints(dl)))
        dzmax = max(map(abs, endpoints(dz)))
        center = endpoints(self.point_upper((l + u) / 2, (a + b) / 2))[1]
        return center + dlmax * (u - l) / 2 + dzmax * (b - a) / 2

    def verify(self, max_boxes):
        # Source X is sqrt(AMP2)-bounded and centered. Its AMP2-subgaussian
        # bound gives J_lambda=lambda*variance for lambda<=1/(2*AMP2).
        assert LOW_PRECISION == 1 / (2 * AMP2)
        tv = self.tt * self.vv
        rho = (mp.iv.sqrt(1 + 16 * tv * tv) - 1) / (4 * tv)
        gaussian = -tv * (1 - rho) + mp.iv.ln(1 - rho * rho) / 4
        assert endpoints(gaussian)[1] < self.target

        stack = [(LOW_PRECISION, T, F(0), F(1))]
        checked = accepted = meanvalue_accepted = 0
        worst = None
        start = time.monotonic()
        while stack:
            rect = stack.pop()
            checked += 1
            upper = self.monotone_upper(rect)
            used_meanvalue = False
            if upper >= self.target:
                bound = self.meanvalue_upper(rect)
                if bound is not None and bound < upper:
                    upper = bound
                    used_meanvalue = True
            if upper < self.target:
                accepted += 1
                meanvalue_accepted += used_meanvalue
                worst = upper if worst is None else max(worst, upper)
            else:
                l, u, a, b = rect
                if u - l >= T * (b - a):
                    mid = (l + u) / 2
                    stack.extend(((l, mid, a, b), (mid, u, a, b)))
                else:
                    mid = (a + b) / 2
                    stack.extend(((l, u, a, mid), (l, u, mid, b)))
            if checked % 20000 == 0:
                print(json.dumps({"progress_index": self.index, "checked": checked,
                                  "accepted": accepted, "pending": len(stack),
                                  "seconds": time.monotonic() - start}), flush=True)
            if checked >= max_boxes:
                return {"status": "BUDGET_STOP_NOT_CERTIFIED", "index": self.index,
                        "checked": checked, "accepted": accepted, "pending": len(stack)}
        return {"status": "DIRECTED_INTERVAL_CERTIFICATE", "index": self.index,
                "density": str(self.p), "amplitude_squared": str(AMP2),
                "t": str(T), "E_upper": str(self.target),
                "gaussian_interval": list(map(str, endpoints(gaussian))),
                "checked": checked, "accepted": accepted,
                "meanvalue_accepted": meanvalue_accepted,
                "worst_accepted": str(worst), "seconds": time.monotonic() - start}


def entropy_iv(theta):
    if theta == 0 or theta == 1:
        return mp.iv.mpf(0)
    th = rational_iv(theta)
    return -th * mp.iv.ln(th) - (1 - th) * mp.iv.ln(1 - th)


def chord_verification(precision=40):
    mp.iv.dps = precision
    alpha = P / 2
    rows = []
    for j in range(10):
        lo, hi = F(j, 20), F(j + 1, 20)
        left = (TARGETS[j] + TARGETS[20 - j]) / 2
        right = (TARGETS[j + 1] + TARGETS[19 - j]) / 2
        slope = (right - left) / (hi - lo)
        intercept = left - slope * lo
        aa, bb, dd = map(rational_iv, (alpha, slope, intercept))

        def derivative(theta):
            th = rational_iv(theta)
            return aa * mp.iv.ln((1 - th) / th) + bb

        if lo > 0 and endpoints(derivative(lo))[1] <= 0:
            bound = aa * entropy_iv(lo) + bb * rational_iv(lo) + dd
            location = "left endpoint"
        elif endpoints(derivative(hi))[0] >= 0:
            bound = aa * entropy_iv(hi) + bb * rational_iv(hi) + dd
            location = "right endpoint"
        else:
            # Global entropy-plus-affine maximum is always an upper bound,
            # including when interval comparison cannot locate its maximizer.
            bound = aa * mp.iv.ln(1 + mp.iv.exp(bb / aa)) + dd
            location = "global entropy-affine maximum"
        row_exponent = bound + aa * mp.iv.ln(2)
        cap = (rational_iv(T) + row_exponent) / (2 * rational_iv(T) * mp.iv.sqrt(rational_iv(P)))
        upper = endpoints(cap)[1]
        assert upper < F(499, 1000)
        rows.append({"theta_interval": [str(lo), str(hi)], "bound_location": location,
                     "row_upper": str(endpoints(row_exponent)[1]),
                     "cap_upper": str(upper),
                     "cap_float_DISPLAY_ONLY": float(upper)})
    return {"status": "DIRECTED_DENSITY_CHORD_CERTIFICATE_CONDITIONAL_ON_POINT_BOUNDS",
            "target_cap": "499/1000", "intervals": rows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--indices", default=",".join(str(j) for j in range(1, 21)))
    parser.add_argument("--max-boxes", type=int, default=1000000)
    parser.add_argument("--precision", type=int, default=40)
    parser.add_argument("--chords-only", action="store_true")
    args = parser.parse_args()
    if args.chords_only:
        print(json.dumps(chord_verification(args.precision), indent=2))
        return
    results = []
    for index in map(int, args.indices.split(",")):
        if not 1 <= index <= 20:
            raise ValueError("index must be between 1 and 20")
        result = FixedAmplitudeCertificate(index, args.precision).verify(args.max_boxes)
        results.append(result)
        print("POINT_RESULT " + json.dumps(result), flush=True)
        if result["status"] != "DIRECTED_INTERVAL_CERTIFICATE":
            break
    print("FINAL_RESULT " + json.dumps({"point_results": results,
                                        "chords": chord_verification(args.precision)}), flush=True)


if __name__ == "__main__":
    main()
