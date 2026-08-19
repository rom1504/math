# Inhomogeneous bridge discrepancy: exact boundary and method no-go

**Status.**  The exact inhomogeneous convex-body reduction and its
power-scale implication are proved below for actual exact child minimizers.
The direct Lovett--Meka/entropy-threshold implementation is proved to fail
its own hypothesis exponentially on balanced target-scale instances.  A
universal near-maximum shell theorem and a frozen order-`8+3` computation
show why unlabeled slack histograms do not contain the missing geometry.

There is **no unconditional exponent improvement** here.  In particular,
this note does not prove

```math
b_{m+n}\le b_m+b_n+o(m+n).
```

The no-go is only for slack/histogram-only coordinate partial coloring and
its standard triangle-additive recursion.  It is not a no-go for all
discrepancy arguments, correlated walks using the rank-one row geometry, or
selection of a favorable actual minimizer class.

## 1. Setup and the exact inhomogeneous body

For a hollow sign matrix `A`, put

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
Q(A)=\max_x|H_A(x)|,
\qquad
M_k=\min_AQ(A),
\qquad
b_k=M_k^{2/3}.
```

Fix `N=m+n` and actual exact minimizers

```math
A\in\operatorname*{argmin}_{A'}Q(A'),
\qquad
D\in\operatorname*{argmin}_{D'}Q(D').
```

Write `p=M_m`, `q=M_n`, and, for `\epsilon\in\{+1,-1\}`,

```math
U_\epsilon(x,y)=H_A(x)+\epsilon H_D(y).
```

For a sign bridge `B\in\{\pm1\}^{m\times n}`, let

```math
P_{\epsilon,B}=
\begin{pmatrix}A&B\\B^{\mathsf T}&\epsilon D\end{pmatrix}.
```

For fixed `x,y`, replacing `y` by `-y` leaves `U_\epsilon` unchanged and
reverses `x^{\mathsf T}By`.  Since

```math
\max\{|u+v|,|u-v|\}=|u|+|v|,
```

one obtains the exact identity

```math
\boxed{
Q(P_{\epsilon,B})
=\max_{x,y}
 \left\{|U_\epsilon(x,y)|+|x^{\mathsf T}By|\right\}.}       \tag{1}
```

Define the ideal power target

```math
T_0=(b_m+b_n)^{3/2}
   =(p^{2/3}+q^{2/3})^{3/2}.                         \tag{2}
```

Since `(u+v)^{3/2}\ge u^{3/2}+v^{3/2}` for `u,v\ge0`, (2) gives
`T_0\ge p+q`, so all widths below are nonnegative.  For `\Delta\ge0`, define the
centrally symmetric convex body in `R^(mn)`

```math
\boxed{
\mathcal K_{\epsilon,\Delta}(A,D)=
\left\{Z\in\mathbb R^{m\times n}:
 |x^{\mathsf T}Zy|
 \le T_0+\Delta-|U_\epsilon(x,y)|
 \quad\hbox{for every }x,y\right\}.}                 \tag{3}
```

Thus the offsets do not shift the bridge intervals: the exact `y -> -y`
symmetry converts them into nonuniform symmetric widths.

### Property `P_gs(Delta;m,n)`

There exist actual exact child minimizers `A,D`, an orientation `epsilon`,
and

```math
B\in\{\pm1\}^{m\times n}
  \cap\mathcal K_{\epsilon,\Delta}(A,D).              \tag{4}
```

This property is existential in the child classes because the recurrence
only needs one favorable pair of actual minimizers.  Requiring (4) for every
pair of actual minimizers is a strictly stronger statement.

By (1)--(4), property `P_gs` implies

```math
M_N\le T_0+\Delta.                                    \tag{5}
```

No surrogate child, conference restriction, or near-minimizer enters this
implication.

## 2. Exact power-scale arrow

Applying the concavity of `s -> s^(2/3)` at `T_0` to (5) gives

```math
\begin{aligned}
b_N-b_m-b_n
&\le (T_0+\Delta)^{2/3}-T_0^{2/3}\\
&\le {2\Delta\over3T_0^{1/3}}
 ={2\Delta\over3\sqrt{b_m+b_n}}.
\end{aligned}                                         \tag{6}
```

Consequently, for any split family satisfying

```math
b_m+b_n\ge\kappa N                                    \tag{7}
```

with a fixed `kappa>0`, the quantitative property

```math
P_{\rm gs}(C N^{3/2-\delta};m,n)
```

implies

```math
\boxed{
b_N-b_m-b_n
\le {2C\over3\sqrt\kappa}\,N^{1-\delta}.}            \tag{8}
```

In particular, `Delta=o(N^(3/2))` implies an `o(N)` power defect whenever
(7) holds.  On comparable large splits, the current lower frontier
`c_*=0.336493364431...` gives

```math
\kappa=c_*^{2/3}-o(1),
\qquad
{2\over3\sqrt\kappa}
={2\over3c_*^{1/3}}+o(1)
=0.9584804067\ldots+o(1).                             \tag{9}
```

This is the exact `P => E_b` conversion.  All missing work is therefore in
producing the sign point (4) with `Delta=o(N^(3/2))`.

### A precise convex-balancing property that would suffice

Let `\gamma_L` be standard Gaussian measure on `\mathbb R^L`, `L=mn`, and
write `(1/5)K=\{z:5z\in K\}`.  The property

```math
\exists A,D,\epsilon:\qquad
\gamma_L\!\left({1\over5}
 \mathcal K_{\epsilon,\Delta}(A,D)\right)\ge{1\over2} \tag{10}
```

is sufficient for (4).  Indeed, apply the standard `5K` form of
Banaszczyk's vector-balancing theorem to the `L` coordinate vectors.  It
produces a sign vector in

```math
5\left({1\over5}\mathcal K_{\epsilon,\Delta}\right)
=\mathcal K_{\epsilon,\Delta}.
```

No known result supplies the optimizer-specific Gaussian estimate (10).
The apparently weaker condition `gamma_L(K)>=1/2` only produces a sign point
in `5K`; that constant dilation costs order `N^(3/2)` here and does not
improve the power-defect exponent.

## 3. Direct Lovett--Meka threshold failure

Specialize to the balanced split `m=n=r`, so `N=2r` and the number of bridge
variables is

```math
L=r^2.
```

Index the absolute rank-one query rows by

```math
q_{x,y}=xy^{\mathsf T}\quad\hbox{modulo }q\sim-q.
```

Independent global flips of `x` and `y` preserve the width and give exactly

```math
K=2^{2r-2}                                             \tag{11}
```

distinct absolute rows.  Every row has Euclidean norm

```math
\|q_{x,y}\|_2=r.                                      \tag{12}
```

Let `A,D` be arbitrary actual order-`r` ground-state minimizers and suppose
`Delta_r=o(r^(3/2))`.  A direct application from zero of the main
Lovett--Meka partial-coloring lemma would use the target widths

```math
d_{x,y}=T_0+\Delta_r-|U_\epsilon(x,y)|,
\qquad
c_{x,y}={d_{x,y}\over r},                             \tag{13}
```

and would require

```math
\sum_{x,y}e^{-c_{x,y}^2/16}\le {L\over16}.            \tag{14}
```

The established asymptotic construction gives

```math
M_r\le(1/2+o(1))r^{3/2}.
```

Consequently

```math
T_0=2\sqrt2\,M_r
\le(\sqrt2+o(1))r^{3/2},                              \tag{15}
```

and every width in (13) is at most `T_0+Delta_r`.  Hence

```math
\begin{aligned}
\sum_{x,y}e^{-c_{x,y}^2/16}
&\ge
2^{2r-2}
\exp\left\{-{(T_0+\Delta_r)^2\over16r^2}\right\}\\
&=\exp\left\{
 \left(2\log2-{1\over8}+o(1)\right)r\right\}.
\end{aligned}                                         \tag{16}
```

The exponential rate is

```math
2\log2-{1\over8}=1.2612943611\ldots>0,                \tag{17}
```

whereas the right side of (14) is only `r^2/16`.  Thus (14) fails
exponentially, uniformly over all actual order-`r` child minimizers and both
orientations.

Giving the first phase only part of the final slack makes (14) still worse.
Therefore a recursive Lovett--Meka proof which certifies the final error by
the triangle sum of nonnegative phase errors cannot begin at the desired
scale.  The same calculation for the classical Beck entropy criterion uses
its `exp(-c^2/9)` branch and has rate

```math
2\log2-{2\over9}=1.1640721389\ldots>0.                \tag{18}
```

so that hypothesis also fails exponentially.

For the Beck criterion, rows with `c\le0.1` lie in its logarithmic branch,
which is still larger than a fixed multiple of the lower bound used in
(18); thus no omitted small-threshold case affects the conclusion.

This proves failure only of the published rowwise threshold certificates.
Large theorem-allowed thresholds do not force the algorithm's realized
discrepancies to be large, and (16) does not exclude a new output-specific
analysis or a correlated walk exploiting relations among the rank-one rows.

## 4. Universal near-maximum shell theorem

The low-slack rows cannot be declared subexponential at any fixed fractional
distance from the child extrema.

### Theorem

Let `A` be any hollow sign matrix of order `k\ge2`, put `p=Q(A)>0`, and choose
`\sigma\in\{\pm1\}` and `x_*` with

```math
\sigma H_A(x_*)=p.
```

For `0<=ell<=k`, define the Hamming sphere

```math
S_\ell(x_*)=\{x:d_H(x,x_*)=\ell\}
```

and

```math
\rho_{k,\ell}
=1-{4\ell(k-\ell)\over k(k-1)}.                      \tag{19}
```

Then for every `\theta<\rho_{k,\ell}`,

```math
\boxed{
\#\{x\in S_\ell(x_*):\sigma H_A(x)\ge\theta p\}
\ge {\rho_{k,\ell}-\theta\over1-\theta}
       {k\choose\ell}.}                              \tag{20}
```

### Proof

Switch by `x_*`, so the maximizing spin is the all-one vector.  For a
uniform `ell`-subset of flipped coordinates and every pair `i<j`,

```math
\mathbb E[x_ix_j]
=1-{4\ell(k-\ell)\over k(k-1)}
=\rho_{k,\ell}.
```

Linearity and `sigma H_A(x_*)=p` give

```math
\mathbb E[\sigma H_A(X)]=\rho_{k,\ell}p.             \tag{21}
```

Let `alpha` be the fraction of the sphere satisfying the event in (20).
Since `sigma H_A(x)<=p` everywhere and the complement has value at most
`theta p`,

```math
\rho_{k,\ell}p
\le\alpha p+(1-\alpha)\theta p.
```

Rearranging proves (20).  `square`

If `ell/k -> delta in (0,1/2)`, then

```math
\rho_{k,\ell}\longrightarrow(1-2\delta)^2.
```

Thus, for every fixed `theta in (0,1)` and

```math
0<\delta<{1-\sqrt\theta\over2},                      \tag{22}
```

(20) and Stirling's formula give at least

```math
\exp\{(h(\delta)+o(1))k\},
\qquad
h(u)=-u\log u-(1-u)\log(1-u),                        \tag{23}
```

spins with energy at least `theta p` in the maximizing orientation.

Apply (20) independently to two actual minimizing children, after the
allowed global matrix orientations have aligned their extrema.  There are
at least

```math
\exp\{m h(\delta_A)+n h(\delta_D)-o(N)\}              \tag{24}
```

paired rows with

```math
U(x,y)\ge\theta_A M_m+\theta_D M_n.                  \tag{25}
```

Therefore any property asserting that every fixed-fraction near-extreme
slack shell of every actual minimizer is subexponential is false.  The
theorem does not settle an `o(N^(3/2))`-thin boundary window: there
`theta -> 1`, the admissible `delta -> 0`, and (23) may itself be only
`exp(o(N))`.

## 5. Same slack histogram, different optimized continuation

The two exhaustive order-eight exact-minimizer classes have `M_8=10` and
the identical projective signed-energy histogram

```text
-10:4, -8:10, -6:12, -4:16, -2:16, 0:12,
  2:16,  4:16,  6:12,  8:10, 10:4.
```

Let `C` be either order-three triangle signing.  It is an actual exact
minimizer with `M_3=3`.  The frozen complete bridge optimization gives

```math
F_C(A_0)=17,
\qquad
F_C(A_1)=19,                                          \tag{26}
```

for both triangle orientations, where

```math
F_C(A)=\min_{B\in\{\pm1\}^{8\times3}}
Q\!\begin{pmatrix}A&B\\B^{\mathsf T}&C\end{pmatrix}.
```

On the other hand the ideal power target is

```math
T_0=(10^{2/3}+3^{2/3})^{3/2}
=17.4267354345015\ldots .                             \tag{27}
```

It follows that

```math
\mathcal K_{\epsilon,0}(A_0,C)
 \cap\{\pm1\}^{8\times3}\ne\varnothing,
\qquad
\mathcal K_{\epsilon,0}(A_1,C)
 \cap\{\pm1\}^{8\times3}=\varnothing.              \tag{28}
```

The common signed-energy histogram makes the complete multiset

```math
\{T_0-|H_{A_i}(x)+\epsilon H_C(y)|:x,y\}              \tag{29}
```

identical for `i=0,1`, for each orientation.  Hence the unlabeled slack
histogram does not determine sign feasibility.  More precisely, any
histogram-only sufficient criterion sharp enough to accept `A_0` at this
target would be unsound if it also accepted every actual minimizer having
that histogram.  A conservative criterion may reject both, and the
cross-order recurrence may select the favorable class `A_0`; neither
possibility is contradicted by (28).

The order-three optimization is a solver-certified exact finite
computation: CP-SAT and an independent HiGHS model both certify `17` versus
`19`, and all feasible witnesses are independently evaluated over every
parent spin.  Neither solver emitted a standalone infeasibility proof
object, so (29) retains that explicit computational qualification.  It is a
finite falsifier, not an asymptotic gap.

## 6. Reproducibility

The frozen files and current SHA-256 hashes are:

```text
620eb3941aad86c70079a52fc77d3af3bf037aabd4a05cb04c593199219cddc1
  extremal_information/experiments/exact_minimizer_optimized_bridge_response.py

24047fb3563b51c983eea3e66907393319dae463fdbaaa63f570452a607e84a5
  extremal_information/experiments/results/exact_minimizer_optimized_bridge_response.json

996d8b34cefbb2ea80e32d47bb4842d1c5f78286dc709b37be1af74069009375
  computations/enumerate_minimizer_orbits.py

119633851ebe024924331b9015354fcc2989eaedde415f49e466376dd39ddbf4
  computations/results/m8_minimizer_orbits.json

10618f45db6a146ebe3949ca0902ac814a2ad8254d28208a3cc160fdc706268b
  computations/results/m3_minimizer_orbits.json
```

The two order-eight canonical orbit hashes are

```text
class 0  a5bbc9a3785f85e929367e670d5e7e0bf6bc46cec302dc054c01b5eee2d07fc9
class 1  75397bf3565083fefc1f5be5b402c0f0bb92c607871ebad3d2a01d53f39c7c5d
```

The reproduction command is

```bash
.venv/bin/python -u \
  extremal_information/experiments/exact_minimizer_optimized_bridge_response.py
```

The protocol and independent audit are

```text
extremal_information/experiments/exact_minimizer_optimized_bridge_response_protocol.md
extremal_information/drafts/exact_minimizer_optimized_bridge_response_independent_audit.md
```

## 7. Primary discrepancy references and final boundary

The exact threshold condition used in (14) is Theorem 4 of
Shachar Lovett and Raghu Meka,
[*Constructive Discrepancy Minimization by Walking on the Edges*](https://arxiv.org/abs/1203.5747).
The convex-body implication (10) uses the standard `5K` normalization of
Wojciech Banaszczyk,
[*Balancing vectors and Gaussian measures of n-dimensional convex bodies*](https://doi.org/10.1002/(SICI)1098-2418(199807)12:4%3C351::AID-RSA3%3E3.0.CO;2-S).

The proved boundary is therefore:

1. a sign point in the exact body (3) with `Delta=o(N^(3/2))` would give the
   desired `o(N)` power defect by (6);
2. standard nonuniform coordinate partial-coloring hypotheses fail
   exponentially before producing such a point;
3. fixed-width near-extreme shells are universally exponential; and
4. even the complete unlabeled slack histogram forgets feasibility-relevant
   labeled geometry on actual minimizers.

A surviving discrepancy route must prove a new optimizer-specific labeled
statement, such as a sufficiently strong Gaussian-measure/section theorem
for (3), or analyze a correlated rounding trajectory beyond its rowwise
threshold guarantees.  No such theorem, and hence no exponent improvement,
is established here.
