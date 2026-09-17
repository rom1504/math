# Kernel coloring, an isotropy obstruction, and an exact non-Gaussian slice repair

2026-09-17. This answers a bounded mechanism question from the
Gram--Schmidt/partial-coloring track. Bare kernel extreme-point rounding
cannot be exactly isotropic in a diffuse projection. A structured repair
can nevertheless give an exactly isotropic, uniformly subGaussian
physical law with response below the paired-Gaussian spectral floor.
No covering theorem for actual minimizing nearcodes is asserted.

## 1. What kernel extreme-point rounding really guarantees

Let P be an orthogonal projector of rank r in R^n and let
gamma=max_i P_ii. Every extreme point z of

```
[-1,1]^n intersect ker(P)
```

has at most r strictly fractional coordinates. Otherwise those columns
of P have a nonzero linear dependency supported on the fractional
coordinates, which gives a feasible segment through z.

Round the fractional coordinates independently to signs with means z_i,
leaving already integral coordinates unchanged. The resulting physical
sign vector h satisfies

```
E h=z,
E ||Ph||_2^2=sum_i (1-z_i^2)P_ii <= r gamma.          (1)
```

Arbitrary mixtures of these rounded vertices, including a global fair
sign to center them, retain (1). But exact isotropy would require
E||Ph||^2=tr(P)=r. Thus if r>0 and gamma<1, NO such mixture alone is
exactly isotropic. This is a second-moment obstruction, not merely a
missing choice of mixture weights.

More quantitatively, put Lambda(P)=max_h ||Ph||^2. If a mixture of
weight 1-alpha on laws satisfying (1), and weight alpha on an arbitrary
repair law, is exactly isotropic, then

```
alpha >= r(1-gamma)/(Lambda(P)-r gamma),             (2)
```

whenever the denominator is positive. A sharper known projection-energy
bound for the first component can be substituted into the same identity.
This shows exactly which covariance resource must be restored.

For example, when n is even and P=11^T/n, the kernel's vertices are
balanced Boolean words. Uniform balanced signs have covariance
n(I-P)/(n-1). Mixing weight 1/n on the two all-aligned words restores
exact isotropy, and makes E|h.1|=1. However the normalized projection
then has atoms of size sqrt(n) with probability 1/(2n), forcing any
subGaussian proxy K I to have

```
K >= n/[2 log(2n)].                                 (3)
```

Indeed its MGF is at least exp(t sqrt(n))/(2n); optimize the resulting
necessary quadratic-exponential inequality at t=2log(2n)/sqrt(n).
Exact isotropy by itself does not make this rare-event repair acceptable
for a dimension-free discrepancy theorem.

## 2. Every fixed slice has a dimension-free centered MGF

Let h be uniform on any nonempty slice sum_i h_i=m. If sum_i t_i=0,
then F=t.h has mean zero and

```
E exp(F)<=exp(||t||_2^2).                            (4)
```

Here is a self-contained exchangeable-pair proof, including constants.
Choose uniformly an ordered pair of distinct coordinates and exchange
their signs, obtaining h' and F'. With lambda=2/(n-1),

```
E[F-F' | h]=lambda F,
E[(F-F')^2 | h]
 <=4/[n(n-1)] sum_(i!=j)(t_i-t_j)^2
 =8||t||^2/(n-1).                                  (5)
```

For theta>=0, exchangeability and the exponential trapezoid inequality
give

```
E F exp(theta F)
 =E[(F-F')(exp(theta F)-exp(theta F'))]/(2lambda)
 <=theta E[(F-F')^2 exp(theta F)]/(2lambda)
 <=2theta ||t||^2 E exp(theta F).
```

Integrate from zero to one to obtain (4); replacing t by -t covers
the other sign. No independence of coordinates within the slice is
assumed. The bound is uniform in the slice magnetization m, including
the degenerate endpoint slices.

## 3. A one-dimensional repair with proxy 2I-P

Take n=12s^2, s a positive integer, and let M=6s=sqrt(3n). Sample h
uniformly from the slice of magnetization zero with probability 2/3,
and from each of the slices of magnetization +M or -M with probability
1/6. The resulting law nu is globally sign symmetric. Exchangeability
and E(sum h_i)^2=n imply exactly

```
E h=0,             E hh^T=I.                         (6)
```

In the flat direction it has

```
E |h.1|/sqrt(n)=1/sqrt(3).                           (7)
```

This restores variance at a CONSTANT magnetization scale sqrt(n),
rather than through the rare size-n atoms in (3).

Let J=(sum h_i)/sqrt(n), so J has probabilities 2/3 at zero and
1/6 at each of +-sqrt(3). It is exactly 1-subGaussian:

```
E exp(uJ)=(2+cosh(sqrt(3)u))/3 <=exp(u^2/2).         (8)
```

One proof compares even Taylor coefficients: for k>=1,
E J^(2k)=3^(k-1)<=(2k-1)!!, the corresponding standard Gaussian
moment. Decompose t=t_perp+bar(t)1 and apply (4) conditionally on J,
then (8). With P=11^T/n,

```
E exp(t.h)
 <=exp(||t_perp||^2+n bar(t)^2/2)
 =exp(t^T(2I-P)t/2).                                (9)
```

Thus the law is exactly isotropic and has the explicit anisotropic
subGaussian proxy 2I-P. In particular its projection to the cheap
flat direction has the SHARP variance proxy one, despite the absolute
response 1/sqrt(3)<sqrt(2/pi). The non-Gaussian distribution is essential:
it is not obtained by assuming a Gaussian absolute-moment formula from
its covariance.

To obtain full physical support, mix weight 1/8 of independent fair
signs into nu. Exact isotropy and (9) survive, since I<=2I-P. The
response at the flat word is at most

```
7/(8sqrt(3))+1/8 < .65.                              (10)
```

The same construction can be switched coordinatewise to any fixed
Boolean center u. For any query x at Hamming distance at most rho*n
from u or -u, exact isotropy and Cauchy--Schwarz give the uniform bound

```
E |h.x|/sqrt(n) <=7/(8sqrt(3))+1/8+2sqrt(rho).        (11)
```

This pays every query in the specified tube; it does not assert that
an actual minimizing nearcode lies in that tube.

## 4. A whole block-constant code with proxy 2I+P

Partition N coordinates into blocks of sizes n_b=12s_b^2. Let P project
onto the subspace of blockwise constant vectors. Use ONE common active
indicator A with probability 1/3. If A=0, sample independently a uniform
balanced slice in every block. If A=1, choose independent fair poles
epsilon_b and independently sample the block slice with magnetization
epsilon_b sqrt(3n_b).

Within each block, the unconditional magnetization second moment is
n_b; between different blocks it is zero because their active poles
are independent. Hence this complete physical sign law is centered
and exactly isotropic on all N coordinates.

For t_perp=(I-P)t, conditional application of (4) costs
exp(||t_perp||^2). Put a_b=sqrt(n_b)bar(t_b). The remaining MGF is

```
2/3+(1/3) product_b cosh(sqrt(3)a_b)
 <=exp((3/2)sum_b a_b^2).
```

Consequently

```
E exp(t.h)<=exp(t^T[2(I-P)+3P]t/2)
           =exp(t^T(2I+P)t/2).                      (12)
```

Every block-constant Boolean query x has, EXACTLY,

```
E |h.x|=(1/sqrt(3)) E |sum_b sqrt(n_b)epsilon_b|
          <=sqrt(N/3).                              (13)
```

Thus one centered isotropic 3-subGaussian law protects the ENTIRE
block-constant code. The rank of P may grow; no Gaussian limit or
fixed-rank condition was used. For r equal-sized blocks, the normalized
response is E|S_r|/sqrt(3r), tending to sqrt(2/(3pi)). Adding 1/8
independent-cube mass gives full physical support and the same upper
bound (10), simultaneously for every code word. The Hamming-tube
extension (11) applies to the whole code as well.

## 5. What this does and does not replace

The block law is an exact non-Gaussian physical replacement with a
dimension-free MGF and a genuinely small mean response, not a weighted
Gaussian rounding surrogate. It shows that the sharp Gaussian spectral-
band floor is not a barrier for all isotropic subGaussian sign laws.

On the other hand, Section 1 prevents simply declaring a kernel
extreme-point mixture isotropic. The successful repair above exploits
exact slice symmetry and explicit magnetization distributions. It does
not prove such a repair for an arbitrary high-covariance or signed-energy
spectral projector, and the projector need not be a block projector
in the actual optimization problem. There is also no claim that a
whole minimizing nearcode has the required block-mode or Hamming-tube
structure. Those are substantive missing hypotheses, not consequences
of the covariance identity.

The construction and swap-MGF proof were independently reconstructed
by the Bernoulli researcher. The director's separate diffuse finite-rank
tilting theorem generalizes the non-Gaussian variance-allocation idea
using smooth invariance and a small exact-isotropy repair; it has
different fixed-rank/diffuseness hypotheses and does not subsume this
growing block-rank finite statement.

## 6. Exact replay

```
.venv/bin/python computations/paper_discrepancy_2026_09_17_kernel_slice_repair.py
```

The replay enumerates all 4096 words at n=12 and verifies exact
centering, covariance I, the slice response 2, and full support after
the 1/8 repair. The latter has exact unnormalized center response
4277/2048. It also checks 1920 swap increments and one hundred exact
three-point moment inequalities. The block response formulas are
evaluated exactly before separately labeled decimal normalization.
These checks accompany the proofs above; no finite grid of MGF tests
is substituted for the exchangeable-pair argument.
