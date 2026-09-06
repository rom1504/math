# Local full-noise contractions against bounded coherent functions

2026-09-06. Finite-degree lemma, independently reconstructed in messages
by the director and adversarial audit agent. Section 4's weighted-tail
and multiplication bounds were also independently audited. It concerns LOCAL
contractions. It does not promote one small factor to a two-root open
Frobenius bound after summing over two output indices.

## 1. Bounded chain rule controls the local Stein product

Let X be a centered homogeneous Boolean polynomial of fixed degree P,
with bounded variance and all proper cuts at most epsilon_n, where
epsilon_n times any fixed power of log(n+1) tends to zero. Let W be a
fixed finite vector of polynomials of degrees at most M<P, with bounded
second moments. Let A:R^d->R be bounded C^2 with bounded first and second
derivatives. Define the EXACT Boolean product

```
Gamma(X,A(W))=(1/P) sum_a Delta_a X Delta_a[A(W)].
```

Then

```
||Gamma(X,A(W))||_2 <= C epsilon_n log(n+1)^C -> 0.     (1)
```

Indeed exact first-order Taylor expansion of the finite difference gives

```
Delta_a[A(W)] = sum_l partial_l A(W) Delta_a W_l + R_a,
|R_a| <= C ||Delta_a W||^2.
```

For each l, Gamma(X,W_l) has L2 norm O(epsilon_n), by the standard
proper-contraction argument: W_l has degree less than P. It is itself
a fixed-degree polynomial, so its fixed higher moments obey the same
bound by hypercontractivity. The bounded first derivatives control the
main term. For the remainder,

```
|sum_a Delta_a X R_a|
 <= C max_a |Delta_a X| sum_a ||Delta_a W||^2.
```

Proper X cuts give max_a Inf_a(X)<=C epsilon_n^2. Fixed-degree
hypercontractivity and a union bound give
`||max_a |Delta_a X|||_4<=C epsilon_n log(n+1)^C`.
Also, by Minkowski and hypercontractivity of each Delta_a W_l,

```
||sum_a ||Delta_a W||^2||_4
 <= C_M sum_(a,l) E(Delta_a W_l)^2
 <= C_M sum_l Var(W_l).
```

Hölder proves (1). High own-spin influences in W are allowed. No
Gaussianization or moment determinacy of W appears in this argument.

The same proof permits W to have any fixed degree, not necessarily below
P, provided each primitive mixed product Gamma(X,W_l) tends to zero in
L2 uniformly. Its fixed-degree hypercontractivity then gives the higher
moments used above. The strict inequality M<P is one sufficient way to
verify that primitive mixed hypothesis, not an essential part of the
bounded chain rule. This extension does NOT assert the primitive mixed
hypothesis for a new rich frame without checking it.

There is also a RATE-FREE averaged version, useful when collision
deletion supplies only averaged small influences. Put

```
delta_i=max_a Inf_a(X_i),
gamma_i^2=sum_l ||Gamma(X_i,W_li)||_2^2.
```

Assume the degrees, row variances and hence total influences are
uniformly bounded. Without a union bound or a logarithmic rate,

```
||max_a |Delta_a X_i|||_4^4
 <= sum_a E|Delta_a X_i|^4
 <= C sum_a Inf_a(X_i)^2 <= C delta_i.
```

The same chain argument therefore gives

```
||Gamma(X_i,A(W_i))||_2 <= C(gamma_i+delta_i^(1/4)),
average_i ||Gamma(X_i,A(W_i))||_2^2
 <= C(average_i gamma_i^2 + sqrt(average_i delta_i)).   (1a)
```

In particular averaged primitive mixed-product smallness and averaged
maximum-influence smallness imply averaged bounded-response mixed-product
smallness. No quantitative logarithmic influence rate is needed for
this weaker but often sufficient local conclusion. All constants here
are for a fixed bounded smooth A and a fixed finite polynomial frame.

## 2. Isolating each full contraction by ORIGINAL output degree

Let A_q be the exact original-degree-q Walsh projection of A(W). For
every fixed q>P,

```
||K_X star_P K_(A_q)||F -> 0.                          (2)
```

Project Gamma(X,A(W)) onto original output degree d=q-P. In the exact
Boolean product formula, a contraction sharing r labels has output
degree `P+e-2r` when the right homogeneous degree is e. Therefore only
the FINITE list `e=d-P+2r`, 1<=r<=P, can contribute. The single full-X
term r=P comes from e=d+P=q. Every r<P term is a proper X contraction,
whose norm is at most C epsilon_n ||A_e||_2. These right norms are at
most ||A(W)||_2. Equation (1) and orthogonal degree projection isolate
the full term and prove (2), up to fixed nonzero factorial conventions.

There is no symmetrization gap in this extraction. In the full term
ALL remaining q-P marks belong to the symmetric right tensor A_q;
the full contraction is already symmetric in those remaining marks.
The other terms may be symmetrized/projection-compressed, which only
reduces their Hilbert norm. Infinite higher degrees of A cannot feed
the fixed output d because contracting the degree-P X removes at most
P original labels.

At q=P the same reasoning isolates the scalar covariance. For local
product surgery at fixed degrees, (2) supplies precisely the full
noise-to-coherent contractions that degree comparison alone did not
control. Proper overlaps remain controlled by X's own cuts.

## 3. Norm and approximation limitations

This is a fixed-root estimate, uniform when its stated row bounds are
uniform. An n-by-n matrix of such errors may have entries O(n^-1/2)
and Frobenius norm O(sqrt(n)); that alone is NOT an o(n) nuclear
comparison or an o(sqrt(n)) open two-root tensor estimate.

The result also does not silently assert that ordinary finite-degree
Walsh truncations approximate products X A(W) in L2. Such a passage
requires a weighted degree-tail or stronger approximation estimate.
The next section supplies one for fixed trigonometric A; (1)--(2)
themselves do not require it.

## 4. Trigonometric coefficients have enough ORIGINAL Sobolev tails

For a fixed trigonometric polynomial A(W), every fixed s satisfies

```
sum_q binom(q,s) ||A_q||_2^2
 = E sum_(|U|=s) |Delta_U A(W)|^2 <= C_s,              (3)
```

uniformly under the stated fixed-degree, bounded-row-variance hypotheses
on W. To prove it, use the exact finite cover formula for exp(it dot W).
Each cover factor satisfies
`|exp(it(-2)^|V| chi_V Delta_V W)-1|<=2^|V||t||Delta_V W|`.
For each relative cover, the full U column embeds injectively into the
Cartesian product of its subtuple columns. Its root-row Hilbert norm
is therefore at most the product of the corresponding derivative row
Hilbert norms. For every fixed p,

```
|| (sum_(|V|=r) |Delta_V W|^2)^(1/2) ||_p
 <= C_(M,p) (sum_(|V|=r) E|Delta_V W|^2)^(1/2)
 <= C_(M,p,r) sqrt(Var(W)).
```

This is scalar hypercontractivity followed by Minkowski; no dimension
factor is introduced. Hölder over the finitely many cover factors proves
(3). Constants grow with s and the fixed frequencies, but not with n.
In particular for any fixed r and s>r,

```
sum_(q>Q) (q+1)^r ||A_q||_2^2 <= C_(s,r) Q^(r-s).     (4)
```

There is also a dimension-free multiplication bound for any fixed
degree-P homogeneous Boolean polynomial X and any polynomial Y:

```
||X Y||_2^2 <= C_P ||X||_2^2
                     sum_q (q+1)^P ||Y_q||_2^2.       (5)
```

One deliberately generous explicit choice for P>=1 is
`C_P=(P+1)^2 P! P^P`.

For completeness, decompose a product X Y_q according to the r shared
seed labels, 0<=r<=P. The exact symmetric-tensor contraction has Hilbert
norm at most the product of the two input Hilbert norms, and its
normalizing binomial factors are bounded by C_P(q+1)^(P/2). Distinct
remaining labels are a Hilbert projection. Each resulting output degree
receives at most P+1 different q contributions. Cauchy--Schwarz over
this fixed number proves (5). This is the ordinary fixed-degree product
bound; it does not use a hypercontractive factor exponential in q.
More explicitly, the squared coefficient for r shared labels is

```
P! q! (P+q-2r)! / [r!^2 (P-r)!^2 (q-r)!^2]
 <= P! P^P (q+1)^P.
```

There are at most P+1 overlaps for each input degree and at most P+1
input degrees for each output degree, giving the stated constant.

Apply (5) to the finite tails of A, then use (4) with s>P. Those products
are Cauchy in L2 and converge in probability to X times the L2 limit A;
hence the L2 limit is indeed X A. Thus

```
sup_n ||X(A-A_(<=Q))||_2 -> 0                         (6)
```

when ||X||_2 is uniformly bounded. The same argument applies to a fixed
finite product of polynomial noise factors, using its fixed total degree
and finitely many homogeneous components.

Together, (2) and (6) give an ordered local noise/coherent product surgery
for fixed trigonometric coefficients: truncate original degree first,
use the fixed-degree contraction bounds as n grows, then remove the
weighted tail. This is still a LOCAL statement, not a two-root open
contraction estimate.

## 5. Finite algebra verification

`computations/transfer_seed_local_sobolev_product_checks_2026_09_06.py`
passes 90 direct finite-difference versus Walsh Sobolev identities,
90 finite-cube weighted multiplication tests, and 8726 checks of the
displayed factorial coefficient bound (P<=12, q<=100). The script reads
no input files and writes no generated artifacts. These checks support
the finite formulas, not an asymptotic two-root theorem. Its initial
run used Python's unavailable `int.bit_count`; the checked version uses
an elementary population-count helper compatible with this workspace.
