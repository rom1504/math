"""Finite algebra and pointwise replay for external-field regularization."""

from fractions import Fraction as F
import itertools
import json
import math

import numpy as np
from scipy.special import ndtr, ndtri
from scipy.stats import binom


SEED = 2026091711
rng = np.random.default_rng(SEED)
kappa = math.sqrt(2/math.pi)


def entropy(p):
    return 0.0 if p in (0,1) else -p*math.log(p)-(1-p)*math.log1p(-p)


def envelope(alpha):
    return 0.0 if alpha == 0 else .5*math.log1p(alpha**(2/3))+entropy(alpha**(2/3)/2)


def sharp_envelope(alpha):
    if alpha == 0:
        return 0.0
    if alpha > math.exp(-1):
        return envelope(alpha)
    ell = math.log(1/alpha)
    sigma = math.sqrt(2/ell)
    return alpha/sigma+entropy(sigma*alpha/2)


# Families with exactly known Gaussian widths and entropies.
entropy_checks = 0
for n in (2,4,8,16,64,256):
    for d in range(n+1):
        # One affine Boolean subcube with d unconstrained coordinates.
        alpha = kappa*d/n
        assert d*math.log(2)/n <= envelope(alpha)+1e-12
        assert d*math.log(2)/n <= sharp_envelope(alpha)+1e-12
        entropy_checks += 1
    alpha = kappa/math.sqrt(n)  # An antipodal pair.
    assert math.log(2)/n <= envelope(alpha)+1e-12
    assert math.log(2)/n <= sharp_envelope(alpha)+1e-12
    entropy_checks += 1

pointwise_checks = 0
antipodal_checks = 0
pattern_checks = 0
for n in (4,6,8):
    words = np.array(list(itertools.product((-1.,1.),repeat=n)))
    for _ in range(20):
        a = rng.choice((-1.,1.),size=(n,n))
        a = np.triu(a,1)
        a += a.T
        energy = np.einsum('bi,ij,bj->b',words,a,words)/2/math.sqrt(n)
        offset = np.abs(energy)
        t,s,eta = .3,.07,.2
        g = rng.normal(size=n)
        field = math.sqrt(t)*(words@g)
        base = offset+field
        m = base.max()
        code = base >= m-eta*n
        z = rng.normal(size=(80,n))
        perturb = math.sqrt(s)*(words@z.T)
        increments = np.max(base[:,None]+perturb,axis=0)-m
        near_max = np.max(perturb[code],axis=0)
        assert np.all(near_max <= eta*n+increments+1e-10)
        pointwise_checks += len(z)
        abs_code = np.abs(energy+field) >= np.max(np.abs(energy+field))-eta*n
        assert abs(np.max(np.abs(energy+field))-m) < 1e-10
        # Lexicographic sign words have antipodal index equal to reversed order.
        assert np.all(~abs_code | code | code[::-1])
        antipodal_checks += 1
        k = 3
        labels = rng.choice((-1.,1.),size=(n,k))
        patterns = np.array(list(itertools.product((-1.,1.),repeat=k)))
        values = offset[:,None]+.2/math.sqrt(k)*(words@labels@patterns.T)
        joint = values >= values.max()-eta*n
        individual = values >= values.max(axis=0)[None,:]-eta*n
        assert np.all(~joint | individual)
        pattern_checks += 1

# Exact rational exponents, separate from floating diagnostics.
K,r = F(2,3),F(1,24)
q = K+r
delta = r+K/2-F(1,2)
eta = -F(1,4)
s = eta+delta
comparison = delta+1-K/4
increment = 1+s-delta
width = increment-s/2
failure = 2*eta+1-2*delta
new_child = 3*q/2
physical_window = eta+F(3,2)
code_entropy = 1+(width-1)*F(2,3)
assert (q,delta,s,comparison,increment,width,failure,new_child,physical_window,code_entropy) == (
    F(17,24),-F(1,8),-F(3,8),F(17,24),F(3,4),F(15,16),F(3,4),F(17,16),F(5,4),F(23,24))
assert comparison < increment and K < failure and new_child < physical_window

# One-dimensional W1, computed from exact binomial masses and Gaussian
# CDF antiderivatives. Floating diagnostics, not a replacement for Stein.
wasserstein = []
for q_test in (1,2,3,4,8,16,32,64,128,512):
    sd = math.sqrt(q_test)
    density = lambda z: math.exp(-z*z/2)/math.sqrt(2*math.pi)
    primitive = lambda x: x*ndtr(x/sd)+sd*density(x/sd)
    distance = 2*(sd*density(q_test/sd)-q_test*ndtr(-q_test/sd))
    for j in range(q_test):
        left,right = 2*j-q_test,2*j-q_test+2
        cdf = binom.cdf(j,q_test,.5)
        crossing = sd*ndtri(cdf)
        mid = min(right,max(left,crossing))
        distance += cdf*(mid-left)-(primitive(mid)-primitive(left))
        distance += primitive(right)-primitive(mid)-cdf*(right-mid)
    assert 0 <= distance <= 1+1e-9
    wasserstein.append(dict(q=q_test, W1=float(distance), proven_bound=1))

print(json.dumps(dict(status='PASS',seed=SEED,
    exact_width_entropy_family_checks=entropy_checks,
    deterministic_secondary_field_inequalities=pointwise_checks,
    absolute_polarity_reductions=antipodal_checks,
    joint_nearcode_subset_checks=pattern_checks,
    binomial_gaussian_W1_diagnostics=wasserstein,
    exact_exponents={name:str(value) for name,value in
                    dict(q=q,delta=delta,s=s,comparison=comparison,
                         increment=increment,width=width,failure=failure,
                         child_cap=new_child,physical_window=physical_window,
                         entropy=code_entropy).items()},
    scope='No empirical concentration claim; proof is analytic and deterministic identities are replayed'),indent=2))
