#!/usr/bin/env python3
"""Numerical variational probes, not a partition-sum upper certificate."""

import argparse
import json

import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import minimize_scalar


def profile_pressure(values, weights, t):
    kernel = 0.5 * (np.exp(-t * (values[:, None] - values[None, :])**2)
                    + np.exp(-t * (values[:, None] + values[None, :])**2))
    logscale = np.zeros(len(values))
    for iteration in range(20000):
        update = -np.log(kernel @ (weights * np.exp(logscale)))
        error = np.max(np.abs(update - logscale))
        logscale = 0.5 * (logscale + update)
        if error < 1e-11:
            break
    return float(-2 * weights @ logscale), float(error), iteration


def mixture(p, block_size, order):
    nodes, weights = hermgauss(order)
    positive = nodes > 0
    nodes = np.sqrt(2) * nodes[positive]
    weights = 2 * weights[positive] / np.sqrt(np.pi)
    if block_size == 1:
        return nodes, weights
    values = np.concatenate((np.sqrt(1-p) * nodes,
                             np.sqrt(1-p+p*block_size) * nodes))
    probability = np.concatenate(((1-1/block_size)*weights, weights/block_size))
    return values, probability


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--p', type=float, default=15/16)
    parser.add_argument('--order', type=int, default=160)
    parser.add_argument('--blocks', type=int, nargs='+', default=[1,2,4,8,16,32,64])
    args = parser.parse_args()
    t = np.sqrt(args.p)/(2*(1-args.p))
    for block_size in args.blocks:
        values, weights = mixture(args.p, block_size, args.order)
        entropy = (1-(1-args.p)**block_size)/block_size*np.log(2)
        pressure, error, iterations = profile_pressure(values, weights, t)
        print(json.dumps({'status':'quadrature probe, not certified bound',
                          'p':args.p, 'block_size':block_size, 'quadrature_order':args.order,
                          't':t, 'row_entropy':entropy, 'permanent_pressure':pressure,
                          'tilted_exponent':entropy+pressure/2+t*(1-np.sqrt(args.p)),
                          'scaling_error':error, 'iterations':iterations}), flush=True)


if __name__ == '__main__':
    main()
