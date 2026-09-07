"""Diagnostics of the non-signing pressure countermodel; proof is analytic."""
import math
import numpy as np


def amplitude(n):
    return 0.45+0.01*math.sin(math.log(math.log(n+math.e)))


def profile(n, beta):
    a = amplitude(n)
    return (1-1/n)*a*beta*beta/(np.sqrt(4*a*a+beta*beta)+2*a)


def main():
    temperatures = np.logspace(-4, 3, 81)
    kg = math.pi/(2*math.asinh(1))
    count = 0
    for total in range(2, 129):
        if total % 2 == 0:
            assert np.all(profile(total, temperatures) >=
                          .5*profile(total//2, math.sqrt(2)*temperatures)-1e-12)
        for left in range(1, total):
            right = total-left
            r = left/total
            s = right/total
            beta = temperatures
            parent = profile(total, beta)
            child = r*profile(left, beta*math.sqrt(r))+s*profile(right, beta*math.sqrt(s))
            gap = parent-child
            elementary = .18*beta*beta*r*s/(beta+2)
            cap = np.minimum(2*math.log(2)+beta*beta/2,
                             2*math.log(2)+2*beta*math.sqrt(math.log(2)))
            budget = np.minimum(4*cap, 2*cap+1)+4*kg*cap+np.sqrt(8*kg*cap/math.pi)
            credit = .5*beta*beta*r*s*np.exp(-2*budget/min(r, s))
            assert np.all(gap >= elementary-1e-12*beta*beta)
            assert np.all(gap >= credit-1e-12*beta*beta)
            assert np.all(parent <= (1-1/total)*beta*beta/4+1e-12)
            assert np.all(parent <= (1-1/total)*amplitude(total)*beta+1e-12)
            count += len(beta)
    print('PASS: %d split/temperature cases; countermodel is not an actual signing model' % count)


if __name__ == '__main__':
    main()
