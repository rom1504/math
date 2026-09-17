"""Exhaust all shift extrema of the binomial-vs-normal absolute-mean error.

This is a floating-point stress test, not a proof of the uniform theorem.
For every q, all atoms and all interior stationary points are tested.
"""

from __future__ import annotations

import argparse
import json
import math

import numpy as np
from scipy.special import ndtr, ndtri
from scipy.stats import binom


SQRT_TWO_PI = math.sqrt(2*math.pi)


def normal_density(value):
    return np.exp(-0.5*value*value)/SQRT_TWO_PI


def extrema(q):
    indices = np.arange(q+1)
    atoms = (2*indices-q)/math.sqrt(q)
    probabilities = binom.pmf(indices,q,0.5)
    probabilities /= probabilities.sum()
    tail_prob = np.concatenate((np.cumsum(probabilities[:0:-1])[::-1],[0.0]))
    tail_moment = np.concatenate((np.cumsum(
        (probabilities*atoms)[:0:-1])[::-1],[0.0]))
    # Both arrays at index k refer strictly to atoms with index > k.
    nonnegative = atoms >= 0
    locations = atoms[nonnegative]
    errors = 2*(tail_moment[nonnegative]-locations*tail_prob[nonnegative]
                -normal_density(locations)+locations*ndtr(-locations))
    kinds = ["atom"]*len(locations)
    valid = (tail_prob[:-1] > 0) & (tail_prob[:-1] < 1)
    roots = np.full(q,np.nan)
    # Phi(root)=CDF on the gap, equivalently Phi(-root)=strict upper tail.
    roots[valid] = -ndtri(tail_prob[:-1][valid])
    inside = valid & (roots >= atoms[:-1]-1e-12) & (roots <= atoms[1:]+1e-12) & (roots >= -1e-12)
    gap_roots = roots[inside]
    # At a stationary point the t*tail terms cancel exactly.
    gap_errors = 2*(tail_moment[:-1][inside]-normal_density(gap_roots))
    locations = np.concatenate((locations,gap_roots))
    errors = np.concatenate((errors,gap_errors))
    kinds.extend(["stationary"]*len(gap_roots))
    maximum = int(np.argmax(np.abs(errors)))
    if q <= 20:
        for test in (0.0,0.123,0.8,1.3,2.2):
            direct = float(np.dot(probabilities,np.abs(atoms-test)))
            gaussian = 2*float(normal_density(np.asarray(test)))+test*(2*float(ndtr(test))-1)
            index = int(np.searchsorted(atoms,test,side="right"))-1
            if index < 0:
                stop = -test
            else:
                stop = tail_moment[index]-test*tail_prob[index]
            cancellation_safe = 2*(stop-float(normal_density(np.asarray(test)))+test*float(ndtr(-test)))
            assert abs((direct-gaussian)-cancellation_safe) < 1e-12
    return {"q":q,"max_absolute_error":float(abs(errors[maximum])),
            "q_times_max_error":float(q*abs(errors[maximum])),
            "signed_error_at_max":float(errors[maximum]),
            "maximizing_nonnegative_shift":float(locations[maximum]),
            "extremum_type":kinds[maximum],
            "tested_extrema":len(locations),
            "stationary_points":len(gap_roots)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-q",type=int,default=1000)
    parser.add_argument("--full",action="store_true")
    args = parser.parse_args()
    cases = [extrema(q) for q in range(1,args.max_q+1)]
    assert all(case["q_times_max_error"] < 3 for case in cases)
    worst = max(cases,key=lambda case:case["q_times_max_error"])
    report = {"status":"PASS numerical stress test at every shift extremum",
              "not_a_uniform_analytic_proof":True,"max_q":args.max_q,
              "total_extrema":sum(case["tested_extrema"] for case in cases),
              "worst_scaled_error":worst,
              "largest_even_scaled_error":max((case for case in cases if case["q"]%2 == 0),
                                             key=lambda case:case["q_times_max_error"]),
              "margin_below_proposed_constant_3":3-worst["q_times_max_error"],
              "selected_cases":[case for case in cases if case["q"] in
                                (1,2,3,4,5,10,20,50,100,200,500,1000)]}
    if args.full:
        report["all_cases"] = cases
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()

