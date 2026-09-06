"""Floating discriminator only: paired-query Haar-involution state evolution.

This is not an interval certificate and does not give a universal lower bound
for arbitrary hollow sign matrices. It tests fixed inertial local updates on
the flat-involution/Haar comparator class.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc


def run(power, steps, inertia, rule, seed, schedule=None):
    samples = 2**power
    draws = qmc.Sobol(d=steps+2, scramble=True, seed=seed).random_base2(power)
    spin = np.where(draws[:, 0] < .5, -1., 1.)
    normals = ndtri(draws[:, 1:])
    q, p = [], []
    u = spin.copy()
    trace = []
    for step in range(steps+1):
        aq = np.array([np.mean(x*u) for x in q])
        ap = np.array([np.mean(x*u) for x in p])
        residual = u.copy()
        field = np.zeros(samples)
        for x, y, a, b in zip(q, p, aq, ap):
            residual -= a*x+b*y
            field += a*y+b*x
        variance = np.mean(residual**2)
        sigma = np.sqrt(max(0., variance))
        z = normals[:, step].copy()
        # Finite-quadrature orthogonalization only; disappears in the
        # population recursion. It reduces spurious accumulated overlaps.
        for x in q+p:
            z -= np.mean(z*x)*x
        if sigma > 1e-9:
            newq = residual/sigma
            z -= np.mean(z*newq)*newq
            z /= np.sqrt(np.mean(z*z))
            field += sigma*z
            q.append(newq)
            p.append(z)
        energy = float(np.dot(aq, ap))
        stability_gap = float(np.mean(np.abs(field)-u*field))
        trace.append(dict(step=step, energy=energy,
                          residual_variance=float(variance),
                          empirical_energy=float(np.mean(u*field)/2),
                          field_norm=float(np.mean(field*field)),
                          stability_gap=stability_gap,
                          norm=float(np.mean(u*u))))
        step_inertia = inertia if schedule is None else schedule[min(step, len(schedule)-1)]
        if rule == 'sign':
            u = np.sign(field+step_inertia*u)
        elif rule == 'clip':
            u = np.clip(u+step_inertia*field, -1., 1.)
        else:
            raise ValueError(rule)
    return dict(status='floating only, not a certificate', power=power,
                steps=steps, inertia=inertia, schedule=schedule, rule=rule, seed=seed,
                trace=trace)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--power', type=int, default=17)
    parser.add_argument('--steps', type=int, default=8)
    parser.add_argument('--inertia', type=float, default=.5)
    parser.add_argument('--rule', choices=['sign', 'clip'], default='sign')
    parser.add_argument('--seed', type=int, default=617)
    parser.add_argument('--schedule', help='comma-separated inertias; last is repeated')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    schedule = None if args.schedule is None else [float(x) for x in args.schedule.split(',')]
    result = run(args.power, args.steps, args.inertia, args.rule, args.seed, schedule)
    encoded = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(encoded)
    else:
        print(encoded, end='')


if __name__ == '__main__':
    main()
