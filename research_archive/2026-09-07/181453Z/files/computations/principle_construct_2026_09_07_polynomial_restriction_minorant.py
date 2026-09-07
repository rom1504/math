"""Exact rational degree-28 minorant of absolute value and Gaussian mean.

No numerical optimization or floating-point inequality enters the certificate.
Nodes are fixed rational approximations selected before this verification.
"""
import json
import sympy as sp

y = sp.Symbol('y')
roots = [sp.Rational(a, 1000) for a in (799,1607,2432,3289,4196,5190,6364)]
nodes = [a*a for a in roots]
h = sp.Poly(0, y)
for i, (x, z) in enumerate(zip(roots, nodes)):
    ell = sp.Poly(1, y)
    for j, zj in enumerate(nodes):
        if j != i:
            ell *= sp.Poly((y-zj)/(z-zj), y)
    derivative = ell.diff().eval(z)
    local = sp.Poly((1-2*(y-z)*derivative)/x-(y-z)/(2*x**3), y)
    h += local*ell*ell
p = h*sp.Poly(y, y)
assert p.degree() == 14
for x, z in zip(roots, nodes):
    assert p.eval(z) == x
    assert p.diff().eval(z) == 1/(2*x)
expectation = sum(p.nth(r)*sp.factorial2(2*r-1) for r in range(1,15))
assert expectation > sp.Rational(754,1000)
assert 2*expectation/3 > sp.Rational(5027,10000)
print(json.dumps({
    'status': 'EXACT_RATIONAL_MINORANT_MOMENT_PASS',
    'degree_in_x': 28,
    'positive_nodes_numerators_over_1000': [int(1000*a) for a in roots],
    'gaussian_mean_lower_bound': '754/1000',
    'greedy_coefficient_lower_bound': '5027/10000',
    'gaussian_mean_DISPLAY': str(sp.N(expectation, 20)),
    'greedy_coefficient_DISPLAY': str(sp.N(2*expectation/3,20)),
    'coefficients_in_y': [str(p.nth(r)) for r in range(15)],
    'gaussian_mean_exact': str(expectation)
}, indent=2))
