# Exact half-normal profile pressure and the fixed-tilt weave obstruction

Date: 2026-09-06. Independent reconstruction of the convergence agent's
Gaussian-profile lower obstruction. The variational optimization in
Section 1 is exact. Its application to the actual Boolean partition sum
is only a LOWER bound, never a Gaussian replacement for that partition sum.

## 0. Director's complementary finite-type upper tool

The complete file
`continued_director_permanent_profile_variation_2026_09_06.md` was read
and independently reconstructed. A contingency table N_ab with margins
n_a is realized by exactly product_a(n_a!)^2/product_ab N_ab!
permutations. Dividing by d! proves its finite formula. Applying
n!<=[e(n+1)](n/e)^n only to the numerator factorials and
n!>=(n/e)^n to the denominator gives the prefactor
e^(2l)(d+1)^(2l); the at most (d+1)^(l^2) tables give precisely the
displayed quantitative upper bound.

For fixed l and positive K, round down an approximately optimizing
coupling with the CURRENT exact margins. The remaining integral row
and column deficits total O(l^2) and admit a nonnegative bipartite
transport, so Stirling also gives the fixed-type matching lower limit.
Positivity of K matters for this last filling argument; zeros cause
no problem for the upper bound. Binwise supremum domination of a
nonnegative permanent is valid without requiring the bin matrix to
remain positive semidefinite. The maximum over the removed diagonal
coordinate must still be taken, as the director states.

The final profile-count bound is conditional on a genuine uniform
count of actual Boolean rows. It does not infer such a count from a
typical Gaussian spectral law. This upper tool and the lower obstruction
below therefore have complementary, correctly separated quantifiers.

## 1. Exact optimizer of the limiting profile problem

Let nu be the law of |N(0,1)| and define

```math
K_t(a,b)=\tfrac12\{e^{-t(a-b)^2}+e^{-t(a+b)^2}\}
=e^{-t(a^2+b^2)}\cosh(2tab),\qquad a,b\ge0.
```

For t>0 choose the unique rho in (0,1) satisfying
2t=rho/(1-rho^2). Let pi_rho be the law of (|G_1|,|G_2|) for a
standard bivariate normal pair with correlation rho. Its density ratio
relative to nu tensor nu is

```math
r_\rho(a,b)=\frac1{\sqrt{1-\rho^2}}
 \exp\left\{-\frac{\rho^2(a^2+b^2)}{2(1-\rho^2)}\right\}
 \cosh\left(\frac{\rho ab}{1-\rho^2}\right).
```

Consequently log K_t-log r_rho is a sum of marginal potentials:

```math
\tfrac12\log(1-\rho^2)
 +\left(-t+\frac{\rho^2}{2(1-\rho^2)}\right)(a^2+b^2).
```

Every coupling pi with marginals nu has E(a^2+b^2)=2. Therefore

```math
\int\log K_t\,d\pi-D(\pi\Vert\nu\otimes\nu)
=-F(t)-D(\pi\Vert\pi_\rho),
```

```math
F(t)=2t(1-\rho)-\tfrac12\log(1-\rho^2).                 (1)
```

This proves that the permanent profile variational value is EXACTLY
-F(t), uniquely attained at pi_rho. The equality is not merely a trial
coupling bound. Infinite entropy couplings cause no exception: log r_rho
has at most quadratic growth and the marginal second moments are fixed.
At t=0 the statement extends with rho=0 and F(0)=0.

## 2. Typical actual partial-Hadamard profiles

Let H be ANY real Hadamard matrix of order m. Independently choose a
uniform k-subset T of its rows and uniform Boolean spins x on T, with
k/m tending to p in (0,1). Define

```math
a_j=|(H[T,:]^Tx)_j|/\sqrt k.
```

Every single coordinate is the absolute value of a normalized sum of
k independent signs. For any distinct column pair, exactly m/2 original
rows have equal column signs. The number M of selected agreements is
hypergeometric and M/k tends in probability to 1/2. Conditional on M,
the signed pair has the law

```math
(U+V,U-V)/\sqrt k,
```

where U and V are independent sign sums of lengths M and k-M. Its
characteristic function converges to that of independent standard
normals. This argument is identical for every distinct column pair,
regardless of the Hadamard family.

For each fixed bounded continuous test function, the empirical profile
average therefore has convergent mean and variance tending to zero.
The empirical law converges in probability to nu. In addition,

```math
\mathbb E\frac1m\sum_j a_j^4=3-2/k,\qquad
\mathbb E\frac1m\sum_j(a_j-R)_+^2\le3/R^2.             (2)
```

These are actual Boolean-seed statements, not a statement that the m
spectral coordinates are independent Gaussians at exponential scale.

## 3. Tail control and finite types give the needed lower bound

For every epsilon>0, writing a_R=min(a,R), direct scalar squaring gives

```math
K_t(a,b)\ge
\exp\{-2t(1+1/\epsilon)[(a-R)_+^2+(b-R)_+^2]\}
K_{t(1+\epsilon)}(a_R,b_R).                              (3)
```

Thus the square root of the normalized permanent loses at most
2t(1+1/epsilon) times the SUM of the squared tails in its logarithm.
Both the row and column occurrence of each coordinate are counted.

The clamped Gaussian coupling gives a useful uniform lower comparison.
Clamping decreases |a-b| and a+b, so it increases K_s(a,b). Pushing
pi_rho and nu through the clamp cannot increase relative entropy.
Consequently the clamped profile variational value at K_s is at least
-F(s) for every fixed R.

Partition [0,R] into finitely many bins and use the INFIMUM of K_s on
each bin pair. The exact contingency-table permanent formula, with
Stirling bounds from below and integer tables of the exact empirical
margins, gives the profile variational lower bound. The kernel is
strictly positive and continuous on this compact square, so binning
error tends uniformly to zero with the mesh. Bin frequencies converge
in probability by Section 2. The clamped endpoint atom can be placed
inside the last bin, with all other boundaries continuity points.

For any fixed target tail tolerance choose R large enough that (2)
and Markov's inequality leave a positive-probability good-tail event.
Intersect it with the high-probability empirical-bin event. This event
has probability bounded away from zero, which is enough for an
exponential LOWER bound on the expected permanent; no exponentially
strong concentration claim is needed. Fix all truncation/bin parameters
before taking m to infinity, remove the desired error afterward, and
finally send epsilon to zero in (3).

Deleting one spectral coordinate changes the bounded empirical-bin
frequencies by O(1/m) and cannot increase the unnormalized tail sum.
Thus this argument applies to any one specified deletion, and hence to
the maximum-over-deletions square-root permanent L_t in the actual soft
weave theorem. Since E_T Z_T=2^k E_(T,x)L_t, it proves

```math
\liminf_m\frac1m\log\mathbb E_T Z_T(t)
\ge p\log2-\frac12F(t).                                 (4)
```

This bound leaves room for rare non-Gaussian spectral profiles to make
the actual partition sum exponentially LARGER.

## 4. Consequence for the fixed-tilt sufficient certificate

The soft weave certificate at restriction proportion p adds the tilt
cost t(1-sqrt(p)) before any desired positive cap gap. By (4), its
limiting exponent is bounded below by

```math
g_p(t)=p\log2-\tfrac12F(t)+t(1-\sqrt p).
```

Implicit differentiation of (1) gives F'(t)=2(1-rho), so
g_p'(t)=rho-sqrt(p). The global minimum is attained at

```math
\rho=\sqrt p,\qquad t=\frac{\sqrt p}{2(1-p)},\qquad
\inf_{t>0}g_p(t)=p\log2+\tfrac14\log(1-p).              (5)
```

This is nonnegative for 0<p<=p_*, where

```math
p_*=0.9225232669048273\ldots,\qquad 1-p_*=2^{-4p_*}.
```

Therefore no FIXED t can satisfy the current strict negative-exponent
sufficient certificate in that proportion range. At p=7/8 the lower
obstruction is (log2)/8>0. At p=15/16 it is -(log2)/16, so this particular
obstruction does not rule the certificate out; it does not prove that
the certificate succeeds there.

The result concerns the specified uniformly random selectors and the
fixed-tilt soft first-moment route. It is not a lower bound on every
optimized restriction, an upper cap certificate, or an exclusion of all
possible n-dependent tilts without further uniform estimates.

## 5. Matching input signs obstruct Jensen across the square root

The complete Section 7 of
`continued_convergence_restricted_weave_2026_09_06.md` has been read
and independently verified. For any Hadamard basis F, random input
signs D, and uniform selector T, let P be the maximum deleted normalized
permanent. The true one-row quantity is 2^k E sqrt(P). Replacing it by
the larger Jensen quantity 2^k sqrt(E P) loses the needed exponent.

Require the effective input signs D x on T to match one fixed row of F.
This costs exactly 2^-k, independently of both F and T. The matching
output spike is sqrt(k) and can be removed. Every remaining output
coordinate is a balanced finite-population sum, and distinct pairs
have four equal sign-class populations by Hadamard orthogonality.
Their typical magnitude profile is half-normal of variance 1-p, with
retained squared norm exactly m-k. The previously audited clamp/type
argument gives

```math
\liminf\frac1m\log[2^k\sqrt{\mathbb EP}]
\ge\frac p2\log2-\frac12F((1-p)t).
```

After the half-cap tilt cost, the optimized lower exponent is

```math
\frac p2\log2+
\frac14\log\left(1-\frac p{(1+\sqrt p)^2}\right)>0.
```

The optimum has rho=sqrt(p)/(1+sqrt(p)). Positivity follows directly
from log(1-u)>=-u/(1-u) and
u/(1-u)=p/(1+2sqrt(p)): the displayed value exceeds
p(log2/2-1/4)>0. Thus this annealed Jensen simplification fails at
every fixed 0<p<1, even if the underlying random Hadamard basis is
otherwise ideal. The rare matching event is harmless at this entropy
scale before Jensen; it is the misplaced square root that pays only
half its probability cost.

This obstruction also survives the HIGH-PROBABILITY selector filters
currently proposed. The matching sign event remains independent of
(F,T), and its residual profile is typical on a set of selectors of
probability tending to one. Conditioning on another event of probability
tending to one preserves that conclusion. This is distinct from the
rare-selector obstruction, which can disappear under such a filter.
No claim is made for arbitrary conditioning on a low-probability set.
