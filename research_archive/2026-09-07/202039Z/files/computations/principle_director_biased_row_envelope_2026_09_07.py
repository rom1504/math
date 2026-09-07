"""Floating falsifier/diagnostic for the marked centered-row certificate.

Reuses the existing rate-distortion dual and reproduction-grid error.
Bounds use floating arithmetic, NOT interval certificates. The analytic
formula being tested is stated explicitly in every output record.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from continued_convergence_latent_supersolution_test_2026_09_06 import SourceEnvelope, ent


class RawEnvelope(SourceEnvelope):
    def __init__(self, p, a, t, grid):
        self.x = np.array([0., (1-a)/np.sqrt(p), (-1-a)/np.sqrt(p)])
        self.p = np.array([1-p, p*(1+a)/2, p*(1-a)/2])
        self.t = t
        self.y = np.linspace(self.x.min(), self.x.max(), grid)
        self.h = self.y[1]-self.y[0]
        self.cache = {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--p', type=float, default=.96)
    ap.add_argument('--biases', default='0,.1,.25,.5,.75,.9,.97,.99')
    ap.add_argument('--temperatures', default='2,4.85,8,16,32')
    ap.add_argument('--grid', type=int, default=151)
    ap.add_argument('--tolerance', type=float, default=.002)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    records = []
    for a in map(float, args.biases.split(',')):
        for t in map(float, args.temperatures.split(',')):
            env = RawEnvelope(args.p, a, t, args.grid).upper(args.tolerance, 160)
            variance = 1-a*a
            slice_entropy = args.p*ent([(1+a)/2, (1-a)/2])
            denom = 2*t*np.sqrt(args.p)*variance
            ratios = [(t*variance+slice_entropy+env[k])/denom for k in ['lower', 'upper']]
            record = dict(p=args.p, bias=a, t=t, variance=variance,
                          slice_entropy=slice_entropy, envelope=env,
                          variance_ratio_lower=ratios[0], variance_ratio_upper=ratios[1])
            records.append(record)
            print(json.dumps(record), flush=True)
            args.output.write_text(json.dumps(dict(
                status='FLOATING DIAGNOSTIC, not an interval certificate or actual cap lower bound',
                formula='(t*(1-a^2)+p*h((1+a)/2)+E_t(nu_p,a))/(2*t*sqrt(p)*(1-a^2))',
                records=records), indent=2)+'\n')


if __name__ == '__main__':
    main()
