# Bounded-response nuclear comparison by a sparse/dense dual split

2026-09-06. **Proved finite-polynomial-frame implication under the explicit
uniform cut and alias hypotheses below.** This is a new way to pass from
local Boolean Stein comparisons to the required nuclear comparison without
polynomial approximation of bounded functions of the literal coherent law.
It does not, by itself, verify these hypotheses for every richer history.

The intended application is the seed agent's fixed rich-frame probe. Its
remaining returned-query problem BC is separate from this mixed covariance
theorem. In particular the conclusion below alone is not a new cap gain.

## 1. Exact finite-frame hypotheses

All polynomial degrees and vector dimensions below are fixed before n.
The underlying seeds are independent unbiased signs. Constants in the
following GLOBAL bounds are uniform in n, not merely polylogarithmic.

Let X_i be a centered homogeneous Walsh polynomial of degree P. Its exact
diagonal-free kernels have all proper fixed-root cuts at most epsilon_n,
where epsilon_n tends to zero, and every global cut, including the root-map
cut, at most C. Write

```
K=E X X^T, q_i=K_ii;  ||K||op<=C_1, max_i q_i<=C_1.
```

Let W_i be a fixed finite list of literal Boolean polynomials of degrees
at most M<P, with uniformly bounded row second moments. Let R_i be another
literal polynomial of fixed degree D, with the same row-moment bound,
whose homogeneous degree-P component is EXACTLY zero. For every q>P,
assume the degree-q component R_q has all proper fixed-root cuts at most
epsilon_n. Its lower-degree components need no small-cut hypothesis.
Uniform fixed moments follow from fixed-degree hypercontractivity.

Fix k>=2 and put p=kP. Let U_i be the exact degree-p squarefree main of

```
q_i^(k/2) h_k(X_i/sqrt(q_i)),
```

with the variance-Hermite convention at q_i=0. Here h_k is orthonormal.
Assume the standard polynomial noise-product surgery and covariance
conclusion, which concern polynomial kernels only:

```
max_i E|U_i-q_i^(k/2)h_k(X_i/sqrt(q_i))|^2 -> 0,
||Cov(U)-K^{circ k}||_*/n -> 0.                         (A)
```

The main U_i is a product of k equal X branches with all original labels
distinct. The X global-cut bounds therefore imply the following weighted
two-factor-transport bound, also for these exact distinct projections:

```
max proper cut of sum_j m_j U_j <= C_2 max_j |m_j|.      (B)
```

This follows from the banked two-factor flat-transport factorization;
its weights need not be rows of a signing. Constants depend on P,k,C.

The only additional alias condition is a ROW-EUCLIDEAN scalar covariance
condition at degree p:

```
A_ij=E[R_{p,i} U_j],       max_i sum_j |A_ij|^2 -> 0.  (C)
```

It is not the desired response covariance in disguise. For the intended
same-old-source application, P>M and k>=3 give (C) by an elementary Hall
bound described in Section 6. No condition is needed if R_p=0.

Let H be fixed bounded continuous on the finite-dimensional W space,
and let psi be fixed bounded smooth, with bounded derivatives (a fixed
sufficient finite differentiability order is enough). Set

```
C_i=H(W_i) psi(R_i+X_i),
d_i=E[H(W_i) psi^(k)(R_i+sqrt(q_i)N)]/sqrt(k!),
```

where N is independent of the ACTUAL joint law of (W_i,R_i). Then

```
||Cov(C,U)-D_d K^{circ k}||_*/n -> 0.                   (1)
```

No exponential moment or polynomial density of the coherent law is used.
The coefficients d_i are uniformly bounded, including at small q_i.
Bounded masks with already justified ordered L2 approximations can be
handled after this smooth/continuous theorem.

## 2. The exact uniform two-root local theorem needed

For ANY roots i,j, apply the finite-cube stable-noise characteristic-
function argument to the two-dimensional noise (X_i,X_j), retaining
the ACTUAL vector (W_i,R_i). The means and moments of that retained
vector are uniformly controlled; it need not be moment-determinate.

All mixed derivative contractions are small. Against W the noise degree
P is larger. Against R_q with q<P they are proper for X; with q>P,
the full P-contraction is proper for R_q. There is no q=P term. The
proper-cut bounds also control the noise-noise derivative fluctuations.
The proof works for a singular covariance K restricted to {i,j}; it
uses characteristic functions, not an inverse covariance matrix.

Thus (X_i,X_j) can be replaced jointly by the centered Gaussian pair
with its ACTUAL covariance, independent of the literal (W_i,R_i),
uniformly in i,j in bounded continuous tests. Fixed-degree moment bounds
allow the degree-k polynomial test in X_j by truncation. The main error
in (A) is uniform as well. Conditional Gaussian Hermite integration by
parts consequently proves

```
eta_n := max_(i,j) |E[C_i U_j]-d_i K_ij^k| -> 0.       (2)
```

This is the only sparse-part local theorem used below. It is a direct
Boolean Stein consequence of the listed proper cuts. It does NOT infer
joint independence from separate marginal limits or use moment density
to approximate H(W).

## 3. Split an arbitrary nuclear dual test

Let M be ANY real n-by-n matrix with ||M||op<=1. Fix tau>0 independently
of n, and write

```
M^large_ij=M_ij 1{|M_ij|>tau},
M^small=M-M^large.
```

Every row and column of M has squared Euclidean norm at most one.
Therefore

```
sum_(i,j)|M^large_ij| <= n/tau,
||M^large||op <= 1/tau,
||M^small||op <= 1+1/tau,
||M^small||F^2 <= n,
max_(i,j)|M^small_ij| <= tau.                           (3)
```

The operator bound follows from row AND column l1 bounds <=1/tau.
Entrywise thresholding is NOT asserted to preserve operator norm.
Equation (2) now controls the large-entry error, after division by n,
by eta_n/tau, which tends to zero at this fixed tau.

For the target matrix, k>=2 and bounded K give

```
sum_(i,j)|K_ij|^k
 <= (max |K_ij|)^(k-2) tr(K^2) <= C_3 n.
```

Consequently the small-entry target contribution is at most
`tau C_3 ||d||infinity n`. There is no dimension loss from this term.

## 4. The dense dual probe is new local Gaussian noise

Put Y=M^small U. It is homogeneous of degree p. Its proper fixed-root
cuts are at most C_2*tau by (B), even though M^small need not have small
operator norm. The latter issue is handled in the CORRECT order:
at each fixed tau, (A) and (3) give

```
average_i E Y_i^2
 = tr(M^small K^{circ k}(M^small)^T)/n+o_tau(1)
 <= C_4+o_tau(1),                                      (4)
```

where C_4 does not depend on tau. The error is bounded by
`(1+1/tau)^2 ||Cov(U)-K^{circ k}||_*/n`, so n must precede tau.

Retain roots with E Y_i^2<=R_0, with R_0 fixed. At those roots use the
same exact Boolean stable-noise argument for Y_i against the literal
list (W_i,R_i,X_i). Lower query degrees give proper Y cuts; higher
degrees give proper right cuts. At the single equal degree p, (C) and
the row l2 bound on M^small imply

```
|E[Y_i R_{p,i}]| <= (sum_j |A_ij|^2)^(1/2) -> 0.       (5)
```

All the resulting mixed derivative bounds are O_R0(tau)+o_n(1).
The characteristic-function proof makes Y_i approximately a centered
Gaussian, independent of the WHOLE actual retained list. Its variance
is its actual E Y_i^2; zero variance is harmless.

More precisely, the bounded row moments and fixed-degree
hypercontractivity imply tightness and uniform integrability. Uniform
characteristic-function bounds, followed by a compactness argument,
give a modulus omega_R0(tau) tending to zero such that

```
limsup_n sup_M |average_(retained i) E[C_i Y_i]|
 <= omega_R0(tau).                                    (6)
```

One can use quantitative smooth-test Stein bounds instead, but no
unproved numerical rate for this modulus is required. The Gaussian
comparison has zero cross, since C_i is a bounded function of the
retained actual list and the comparison Gaussian is independent and
centered. The entire joint law, not just its separate marginals, is
retained.

The discarded roots cost at most

```
||C||infinity sqrt(fraction_discarded * average E Y_i^2)
 <= C_5/sqrt(R_0)+o_tau(1).                            (7)
```

Combine (2)--(7), then take n to infinity at fixed tau and R_0, next
tau down to zero at fixed R_0, and finally R_0 to infinity. Uniformity
in the original operator-norm dual test M proves (1).

## 5. Why this bypasses the difficult Boolean high-order chain rule

Only FIRST-order finite-cube Stein is used. The nonlinear right source
U is handled in two different ways: individual entries use the joint
Gaussian pair (X_i,X_j), whereas a small-entry root combination of U
is itself a new approximately Gaussian homogeneous field by (B).
No derivative of order kP of the bounded response is expanded across
label-dependent flipped configurations. No bounded random coefficient
is pulled out of a closed expected diagram.

The two parts need both local input types. Entrywise convergence alone
would not prove a nuclear conclusion; conversely, the dense transport
cannot be treated as an operator contraction after thresholding. The
fixed-tau order and (4) are what make the argument valid.

## 6. Checking the alias condition in the intended same-source model

Suppose R_p is a flat bounded-op transport of an exact old-source forest
whose primitive degrees are at most M<P. The right U_j consists of k
equal degree-P noise branches. Pull R_p through its single flat B row.
At the source, any collection of k' right noise branches has P*k'
distinct labels but each old primitive has at most M slots. Hall
matching assigns k distinct old primitives to the k right noises.
Each selected merge is proper on the degree-P noise, or retains a
label and uses its small influence. Assuming the established fixed-root
noise cuts are O(n^-1/2), the source covariance is O(n^-k/2).
The one pulled-back B row costs sqrt(n), so

```
|A_ij| = O(n^((1-k)/2)),
sum_j |A_ij|^2 = O(n^(2-k)) -> 0      for k>=3.
```

Internal old coherent collisions are retained by the same finite
equality-diagram argument. This is a polynomial scalar alias bound,
not a bounded-response nuclear comparison. It does not require moments
of the coherent law to determine that law.

This source hypothesis is essential. Arbitrary collections of noise
channels with small proper cuts need not satisfy it: a later channel
formed by transporting U itself has a leading covariance with U.
Hence this theorem cannot be generalized by omitting the same-old-
source condition or (C).

## 7. Exact scope of the rich-frame use

The suggested application takes X as the ENTIRE exact original-P block
of BF and R=BF-X, so the absence of degree P in R is exact. For fixed
polynomial source stages the theorem above is literal. Extending it to
bounded sources requires their separately justified ordered source-L2
approximations, with source structural cuts at each fixed stage. It
does not require polynomial approximation of psi(R) in the law of R.

The uniform-global-cut premise must actually be verified for the
selected old-tree diagrams. The generic all-global-cut theorem giving
polylogarithmic bounds alone does not meet the fixed-tau hypothesis.
A shrinking tau_n would then need quantitative rates in (A), (2), and
the alias estimates; that variant is not silently asserted here.

Even after this mixed covariance is available, a full feedback gain
requires the additional literal returned-query comparison involving
BC-beta E. That remains a separate task in the rich-frame extension.
