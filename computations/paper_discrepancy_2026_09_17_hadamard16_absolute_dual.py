"""Finite complete-eigencode absolute-column dual for the alternating Walsh core.

This is n=16 evidence, not an asymptotic assertion. The final rational
weights are checked on every physical sign column using integer overlaps.
"""

from __future__ import annotations

import argparse
import json

import numpy as np
from scipy.optimize import linprog


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full",action="store_true")
    args = parser.parse_args()
    p = 4
    n = p*p
    labels = [(u,v) for u in range(p) for v in range(p)]
    fourier = np.asarray([[(-1)**(bin(u&b).count("1")+bin(v&a).count("1"))
                           for a,b in labels] for u,v in labels],dtype=np.int64)
    assert np.array_equal(fourier@fourier,n*np.eye(n,dtype=np.int64))
    # Fix the first spin +1; absolute response identifies antipodal words.
    indices = 2*np.arange(1 << (n-1),dtype=np.int64)+1
    words = 2*((indices[:,None] >> np.arange(n))&1)-1
    images = words@fourier
    positive = np.all(images == p*words,axis=1)
    negative = np.all(images == -p*words,axis=1)
    code = words[positive | negative]
    sectors = np.where(positive[positive | negative],1,-1)
    overlaps = np.abs(words@code.T)
    count = len(code)
    result = linprog(np.r_[np.zeros(count),-1.0],
                     A_ub=np.column_stack((-overlaps,np.ones(len(words)))),
                     b_ub=np.zeros(len(words)),
                     A_eq=np.r_[np.ones(count),0.0][None,:],b_eq=[1.0],
                     bounds=[(0,None)]*count+[(None,None)],method="highs")
    assert result.success
    denominator = 165
    weights = np.rint(np.maximum(result.x[:count],0)*denominator).astype(np.int64)
    assert int(weights.sum()) == denominator
    numerators = overlaps@weights
    numerator = int(numerators.min())
    worst = int(np.argmin(numerators))
    report = {"status":"PASS every physical sign query",
                      "order":n,"old_absolute_cap":n*(p+1)//2,
                      "included_deficit_window":n,
                      "projective_positive_eigenwords":int(positive.sum()),
                      "projective_negative_eigenwords":int(negative.sum()),
                      "uniform_complete_code_lower":float(np.min(overlaps.mean(axis=1))),
                      "numerical_lp_lower":float(result.x[-1]),
                      "certified_lower_numerator":numerator,
                      "certified_lower_denominator":denominator,
                      "certified_lower_over_sqrt_n":numerator/(denominator*p),
                      "certified_positive_sector_mass":int(weights[sectors > 0].sum())/denominator,
                      "worst_query":words[worst].tolist()}
    if args.full:
        report["positive_weight_atoms"] = [{"word":word.tolist(),"sector":int(sector),
                                               "numerator":int(weight)}
                                              for word,sector,weight in zip(code,sectors,weights)
                                              if weight]
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()
