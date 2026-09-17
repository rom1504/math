# Cheap physical response or concentration of signed energy in a small covariance subspace

2026-09-17. A structural alternative combining the audited paired
Gaussian-sign response mechanism with a Frobenius cross-term estimate.
It assumes NO operator-norm bound on A and NO covariance bound on the
query law. It does not prove that the exceptional alternative is absent
for actual minimizers. External novelty is unestablished.

**Subsequent closure in a low-cap regime:** the director's
[half-range argument](paper_director_low_cap_uniform_response_2026_09_17.md)
uses the cap budget to exclude this exception when
C<c_half+c0/(2K_G). The independent
[explicit-constant audit](paper_discrepancy_low_cap_response_audit_2026_09_17.md)
proves a uniform isotropic physical response gap for Q(A)<=n^(3/2)/2
and the entire code |H_A|>=(3/10)n^(3/2). The general alternative below
remains valid outside that regime. Its final residual-obligation wording
records the scope before this additional cap-budget argument.

## 1. Finite asymptotic-response theorem

Let A be any hollow symmetric full signing of order n. Let mu be any
probability law supported on |H_A(x)|>=c n^(3/2), where c>0 is fixed.
Write

```
kappa=sqrt(2/pi),
sigma(x)=sign H_A(x),
Sigma=E_mu xx^T,   K=E_mu sigma(x)xx^T,
m=E_mu sigma(x),  K0=K-mI,
Q=1_(Sigma>2),    P=I-Q,
M=tr(Sigma Q),    J=E_mu sigma(x)(Qx)^T A(Qx).
```

Choose any fixed 0<a<=min{1/16,c^2/200}. There is a positive explicit
delta(c,a), independent of n,A,mu, with the following alternative:

1. There is a centered, EXACTLY isotropic physical sign law nu, with
   subGaussian covariance proxy (3/2)I, such that
   E_mu E_nu |h.x|/sqrt(n) <= kappa-delta(c,a)+epsilon_n; or
2. M<an and J>c n^(3/2).

Here epsilon_n=O(n^(-1/6)sqrt(log(en))) is uniform in A and mu. The law
in the first branch is an equal mixture of signs of two centered
Gaussians with unit diagonal and spectra contained in [1/2,3/2].
Thus the theorem does not replace physical signs by a Gaussian query.

One valid explicit choice is

```
t=min{1/12,c/sqrt(1728)},
d_H=kappa^2 (1/2-pi/12)a,
d_L=kappa^2 t c^2/8,
delta(c,a)=(kappa/8) min{d_H^2,d_L^2}.              (1)
```

The exceptional quantity J is a SIGNED projected quadratic energy.
The vectors Qx need not be Boolean, and Q need not be a coordinate
projection. Small M means small average squared Euclidean norm in this
subspace, not that its contribution to the quadratic energy is small.

## 2. Imported response lemma and exact isotropy

We use the bounded-spectrum Gaussian-sign scalar comparison reconstructed
in [the quenched universality proof](flatify_independent_2026_09_07_gaussian_sign_quenched_universality.md),
and the paired-law calculation in
[the signed-covariance response theorem](paper_bernoulli_signed_covariance_response_2026_09_17.md).
Its hypotheses are uniform unit diagonal and spectrum in [1/2,3/2];
it imposes no spectral assumption on the original signing A.

Specifically, if H is symmetric and hollow and ||tH||_op<=1/2, set
R_+=I+tH and R_-=I-tH. For the equal mixture of sign N(0,R_+) and
sign N(0,R_-), the Gaussian arcsine identity gives exactly

```
E_nu hh^T=I,
v_x=(kappa^2/n) x^T asin(tH) x,
E_nu |h.x|/sqrt(n)=kappa f(v_x)+epsilon_(n,x),
f(v)=(sqrt(1+v)+sqrt(1-v))/2,
sup_x |epsilon_(n,x)|<=epsilon_n.                  (2)
```

Arcsine is entrywise, and the diagonal of asin(tH) is zero. The
paired covariance matrices are I plus or minus kappa^2 asin(tH), so
|v_x|<=1. Since f(v)<=1-v^2/8, ANY lower bound
|E_mu s(x)v_x|>=d with |s(x)|<=1 yields

```
E_mu E_nu |h.x|/sqrt(n)<=kappa-(kappa/8)d^2+epsilon_n. (3)
```

The subGaussian proxy (3/2)I follows from the Gaussian product Holder
inequality at covariance R_+ or R_-<= (3/2)I, applied to exponential
functions of the coordinate signs. The mixture retains the same proxy.
These source mechanisms were independently audited in the campaign;
neither the scalar comparison nor Gaussian Holder is newly proved here.

## 3. Large covariance mass gives the first branch directly

Suppose M>=an. Take H=Q-diag(Q), so ||H||_op<=1, and t=1/2.
All coefficients of the entrywise power series asin(z)-z are
nonnegative. Schur products therefore give the matrix inequality
asin(Q/2)>=Q/2. Also asin(q/2)<=q*pi/6 for 0<=q<=1.
Using diag(Sigma)=1 and rank(Q)<=M/2,

```
E_mu v_x
 = (kappa^2/n)[tr(Sigma asin(Q/2))-sum_i asin(Q_ii/2)]
 >= (kappa^2/n)[M/2-(pi/6)rank(Q)]
 >= kappa^2(1/2-pi/12)a=d_H.                       (4)
```

Apply (3) with s(x)=1. Notice that the original signed energy is not
used in this branch; appreciable covariance mass above level two is
already enough to produce a response discount.

## 4. The Frobenius estimate isolates the only remaining energy term

Suppose M<an. Since -Sigma<=K<=Sigma and P Sigma P<=2P, the matrix
B=PK0P has operator norm at most three. Because A is a full signing,

```
||A||_F^2=n(n-1).
```

The cross term is bounded without ||A||_op:

```
|E_mu sigma(x)(Px)^T A(Qx)|
 <= [E_mu ||A Px||_2^2 * E_mu ||Qx||_2^2]^(1/2)
 <= [2 n(n-1) M]^(1/2)
 <= sqrt(2a)n^(3/2).                               (5)
```

The centering correction is also small. Since tr(A)=0, |m|<=1 and
rank(Q)<=an/2,

```
|m<A,Q>|<=||A||_F ||Q||_F<=sqrt(a/2)n^(3/2).        (6)
```

Expand K in its P/Q blocks, retaining the actual sign sigma inside
every expectation. Then

```
<A,B>
 = E_mu sigma x^T A x
   -2E_mu sigma(Px)^T A(Qx)-J+m<A,Q>
 >= [2c-J/n^(3/2)-(5/sqrt(2))sqrt(a)] n^(3/2).       (7)
```

If J<=c n^(3/2), our choice a<=c^2/200 gives
<A,B>>=(3c/4)n^(3/2), in particular at least (c/2)n^(3/2).
Frobenius Cauchy--Schwarz then gives ||B||_F^2>=c^2 n/4.

Now put H=B-diag(B), so ||H||_op<=6. Crucially K0 has zero diagonal
and B=PK0P, whence

```
<K,H>=<K0,H>=<K0,B>=||B||_F^2.                     (8)
```

The Schur tensor-compression inequality gives
||H^(entrywise j)||_op<=||H||_op^j for positive integers j. For
t<=1/12, the arcsine remainder consequently obeys

```
||asin(tH)-tH||_op<=216t^3.
```

Since ||K||_*<=E_mu ||sigma xx^T||_*=n, equations (8) and (1) imply

```
E_mu sigma(x)v_x
 >= kappa^2[t ||B||_F^2/n-216t^3]
 >= kappa^2 t c^2/8=d_L.                            (9)
```

Apply (3). This proves the stated alternative. There is no invocation
of an upper bound on ||A||_op in (5)--(9).

## 5. An asymptotically sharp formulation of the remaining obstruction

Let A_n be any sequence of full signings, and let mu_n be supported
on |H_(A_n)(x)|>=c n^(3/2). Suppose every paired physical law of the
form (2), with spectra in [1/2,3/2], has mu_n-average normalized
response at least kappa-o(1). Then

```
tr(Sigma_n Q_n)=o(n),
rank(Q_n)=o(n),
||P_n K0_n P_n||_F^2=o(n),
E_mu sigma(Q_n x)^T A_n(Q_n x)
       =2 E_mu |H_(A_n)(x)|+o(n^(3/2)).              (10)
```

Indeed (4) forces M/n->0. To obtain the third assertion, set
b_n=||B||_F^2/n<=9. If b_n stays above a positive constant along a
subsequence, choose a sufficiently small FIXED t in (8)--(9); (3)
then gives a fixed response discount, a contradiction. Finally (5)
and (6) are o(n^(3/2)), while
|<A,B>|<=||A||_F ||B||_F=o(n^(3/2)); the exact expansion (7) proves
the final identity in (10).

This conclusion applies in particular to any putative almost-maximal
dual query laws for the unrestricted full-column response game. Uniform
independent columns always have normalized response kappa+o(1), so
almost-maximal here means approaching this independent-column value.
It also applies to the corresponding minimax game restricted to paired
isotropic column laws.

Conversely, if for a complete high-energy code every query law avoids
the exceptional alternative in Section 1, finite-dimensional minimax
over the compact convex hull of the paired laws yields ONE centered
isotropic (3/2)-subGaussian column law with the same fixed discount
at EVERY code word. The exception must be excluded for every query
law; checking one convenient ground-state average is insufficient.

## 6. Exact scope, and why this is not the earlier bounded-spectrum claim

The earlier [adaptive bounded-spectrum proof](paper_director_adaptive_energy_response_2026_09_17.md)
was valid but redundant: for ||A||_op=O(sqrt(n)), an already archived
explicit radial Gaussian pair gives a stronger response discount.
The present statement has no such assumption. Its new content is
the exact alternative (7), and the localization conclusion (10), not
another bounded-spectrum corollary.

A thin coherent subspace CAN carry substantial signed energy, so (10)
cannot be discarded just because its covariance mass tends to zero.
The [near-minimizer covariance-spike examples](paper_discrepancy_nearminimizer_covariance_spike_2026_09_17.md)
also show that exact support on a small nearlevel does not itself force
bounded covariance. Those examples do not prove that the obstruction
in (10) occurs for dual-optimal laws of exact minimizers.

In particular, this theorem does not establish a universal cheap law
for actual minimizing nearcodes, a favorable parent value, or an
asymptotic recurrence. It identifies a precise residual obligation:
rule out, exploit, or cheaply remove leading MARKED energy carried by
a subspace of vanishing covariance mass, while preserving the entire
absolute nearcode and the original full-sign problem.
