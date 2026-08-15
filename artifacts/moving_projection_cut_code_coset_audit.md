# Moving projections for augmented cut-code cosets

Date: 2026-08-15.

Status: **proved coset inequality, verified source mapping, finite numerical
falsification of the canonical rank-one family, and a scalable no-go for a
rank-one transversal escape**.  No bound on `M_n` is improved.

The source mechanism is Chapter 2 of OpenAI's
[*Ten Advances in Mathematics and Theoretical Computer Science*](https://cdn.openai.com/pdf/ten-proofs-oai.pdf),
with its construction history in the
[*Mathematical Discovery Notes*](https://cdn.openai.com/pdf/reasoning-walkthroughs.pdf).
The derivations below use the paper's *complete* mixed-term Gram remainder;
they do not declare its individual harmonic channels positive or pay for
those channels separately.

## 1. Exact cut-code mapping and the packing mismatch

Put

```math
E=\binom n2,
\qquad
\tau(z)={1\over E}\sum_{e\in E(K_n)}z_e,
\qquad z\in\{\pm1\}^E.
\tag{1.1}
```

The augmented cut code is the subgroup

```math
\mathcal C_n^+
=\{(\sigma x_ix_j)_{i<j}:\sigma,x_i\in\{\pm1\}\}.
\tag{1.2}
```

It has `2^n` words.  For an edge signing `a`, define

```math
\mu(a)=\max_{c\in\mathcal C_n^+}\tau(ac).
\tag{1.3}
```

The global sign in (1.2) gives the exact identities

```math
\mu(a)={Q(a)\over E},
\qquad
d(a,\mathcal C_n^+)={E-Q(a)\over2},
\qquad
M_n=E-2\rho(\mathcal C_n^+).
\tag{1.4}
```

For `n>=5`, the minimum distance of `C_n^+` is `n-1`.  A nonconstant cut
has weight `r(n-r)>=n-1`, its complement has weight at least `n-1`, and a
singleton cut attains equality.  Therefore the paper's packing theorem sees
only

```math
s=1-{2(n-1)\over E}=1-{4\over n}.
\tag{1.5}
```

Every translate `a C_n^+` has the same internal distances.  The packing
theorem is consequently coset-blind when applied directly; it cannot locate
the translate relative to the origin and does not itself bound (1.3).

This distinction agrees with the older exact coding formulation of signed
graphs by Solé and Zaslavsky,
[*A Coding Approach to Signed Graphs*](https://doi.org/10.1137/S0895480189174374):
cosets of the cocycle code are switching classes, while the covering radius
is the maximum frustration/imbalance.  That paper supplies the mapping, not
the moving-representation inequality below.

## 2. The unsplit moving Gram remainder

For any of the paper's Boolean moving projections, let `P_z` be the
equal-rank projection attached to `z`.  Translation and permutation
equivariance make

```math
K(z,w)=\operatorname{tr}(P_zP_w)=\kappa(zw)\ge0.
\tag{2.1}
```

Equations (25)--(27) of the paper construct maps `Theta_z` satisfying

```math
\langle\Theta_z,\Theta_w\rangle
=\bigl(\tau(zw)-\lambda\bigr)\kappa(zw).
\tag{2.2}
```

Thus both

```math
\kappa(z)
\quad\hbox{and}\quad
F(z)=\bigl(\tau(z)-\lambda\bigr)\kappa(z)
\tag{2.3}
```

are positive-definite functions on the edge cube, while `kappa` is also
pointwise nonnegative.  The positivity of `F` is precisely where the two
mixed terms in the paper cancel.  It is not a collection of separately
bounded scalar responses.

## 3. A joint-Gram covering inequality

Define

```math
T_a=\sum_{c\in\mathcal C_n^+}\kappa(ac),
\qquad
J_n=\sum_{c\in\mathcal C_n^+}F(c).
\tag{3.1}
```

> **Theorem 3.1 (moving-Gram coset inequality).**  For every edge signing
> `a`,
>
> ```math
> \boxed{(\lambda-\mu(a))_+T_a\le J_n.}
> \tag{3.2}
> ```
>
> If `T_a>0`, this gives
>
> ```math
> \boxed{\mu(a)\ge\lambda-{J_n\over T_a}.}
> \tag{3.3}
> ```

**Proof.**  Use normalized Fourier coefficients on the edge cube.  Since
`F` is positive definite, `F_hat(S)>=0` for every edge subset `S`.  If
`D=(C_n^+)^perp`, subgroup orthogonality gives

```math
\begin{aligned}
S_a:=\sum_{c\in\mathcal C_n^+}F(ac)
 &=|\mathcal C_n^+|\sum_{S\in D}\widehat F(S)\chi_S(a),\\
J_n
 &=|\mathcal C_n^+|\sum_{S\in D}\widehat F(S).
\end{aligned}
\tag{3.4}
```

Hence `-J_n<=S_a<=J_n` and `J_n>=0`.  On the other hand, (2.1) and the
definition of `mu(a)` imply

```math
S_a
=\sum_c(\tau(ac)-\lambda)\kappa(ac)
\le(\mu(a)-\lambda)T_a.
\tag{3.5}
```

If `mu(a)<lambda`, combine the lower bound in (3.4) with (3.5).  If
`mu(a)>=lambda`, (3.2) is automatic.  This proves the theorem. `square`

The inequality is rooted: the internal double sum

```math
\sum_{c,d\in a\mathcal C_n^+}F(cd)
=|\mathcal C_n^+|J_n
\tag{3.6}
```

is translation-blind, but the cross term `S_a` in (3.4) retains the root
`a`.  Fourier positivity is applied only after the paper's complete Gram
remainder (2.2) has been formed.

## 4. Explicit numerator and exact remaining obligation

For a vertex sign vector with `r` negative entries, the normalized edge
correlation of its cut is

```math
q_r={(n-2r)^2-n\over n(n-1)}.
\tag{4.1}
```

Therefore the numerator in (3.3) is completely explicit:

```math
J_n={1\over2}\sum_{r=0}^n\binom nr
\left[(q_r-\lambda)\kappa(q_r)
+(-q_r-\lambda)\kappa(-q_r)\right].
\tag{4.2}
```

The denominator is the one joint weighted coset enumerator

```math
T_a={1\over2}\sum_{x\in\{\pm1\}^n}
\left[
\kappa\!\left({H_a(x)\over E}\right)
+\kappa\!\left(-{H_a(x)\over E}\right)
\right].
\tag{4.3}
```

It is not the maximum in (1.3).  However, a uniform lower bound for (4.3)
is not presently known to be simpler than controlling the original coset
tail.  Theorem 3.1 is therefore a rigorous new interface, not yet primary
progress on convergence.

The dual code is

```math
D=\{S\subseteq E(K_n):
\deg_S(v)\equiv0\pmod2\ \forall v,
\ |S|\equiv0\pmod2\}.
\tag{4.4}
```

Consequently,

```math
{T_a\over|\mathcal C_n^+|}
=\sum_{S\in D}\widehat\kappa(S)\chi_S(a).
\tag{4.5}
```

For the whole-cube construction using levels `k,...,L`, this includes signed
even-Eulerian coefficients through degree `2L`.  The growing hierarchy has
been compressed to one scalar, but it has not disappeared.

## 5. Correct scale and the canonical finite test

The paper's corrected transition weights satisfy

```math
c_i^{(k)}
={(i-k+1)(E-i-k)\over E\sqrt{(i+1)(E-i)}}
\le\sqrt{L+1\over E}.
\tag{5.1}
```

Thus its top eigenvalue obeys

```math
\lambda\le2\sqrt{L+1\over E}.
\tag{5.2}
```

For the rank-one fiber `k=0`, the last `2 by 2` principal block also gives
`lambda=Theta(sqrt(L/E))` when `L=o(E)`.  Hence the desired
`lambda=Theta(n^(-1/2))` is possible only at, and is attained spectrally at,

```math
L=\Theta(n).
\tag{5.3}
```

This is the correct leading scale because

```math
E\,{1\over\sqrt n}
=\left({1\over2}+o(1)\right)n^{3/2}.
\tag{5.4}
```

The accompanying computation evaluated every canonical rank-one level
`1<=L<=min(E-1,3n)` on saved exact/good signing histograms and the order-18
conference histogram.

| vertex order `n` | observed cap | best certified expression `E(lambda-J/T)` | best level |
|---:|---:|---:|---:|
| 6 | 5 | 3.872983 | 1 |
| 8 | 10 | 5.291503 | 1 |
| 10 | 13 | 6.708204 | 1 |
| 12 | 18 | 8.124038 | 1 |
| 14 | 21 | 9.539392 | 1 |
| 18 | 33 | 12.369317 | 1 |

The best value is always `sqrt(E)`, exactly the elementary RMS bound.  At
level one the dual distance four makes `J_n=0`.  At the first level whose
`lambda` exceeds the observed cap ratio, the denominator has already
collapsed; for example,

```math
\begin{array}{c|c|c|c|c|c}
n&L&\lambda&J_n/|\mathcal C|&T_a/|\mathcal C|&\lambda-J_n/T_a\\ \hline
14&3&0.242814&6.3352\cdot10^{-5}&1.0440\cdot10^{-5}&-5.82536\\
18&4&0.229312&5.9675\cdot10^{-6}&3.9254\cdot10^{-8}&-151.793
\end{array}
\tag{5.5}
```

This is floating-point finite evidence, not an asymptotic no-go.  Reproduce
it with

```text
.venv/bin/python computations/audit_moving_projection_cut_code.py
```

## 6. A tempting scalar transversal escape is too small

There is a natural attempt to make (4.5) constant rather than estimate an
unknown coset.  It produces an exact theorem, but a hypercube spectral bound
shows that it cannot reach (5.4).

Write the edge cube additively as `G=F_2^E`, let `D=(C_n^+)^perp`, and let
`A_E` be normalized adjacency on `Q_E`.  Suppose `v>=0` is supported on a
partial `D`-transversal and

```math
A_Ev\ge\lambda v.
\tag{6.1}
```

Define

```math
g(z)=\sum_Rv_R\chi_R(z),
\qquad
\kappa(z)=g(z)^2.
\tag{6.2}
```

If `q=v*v` is the Fourier coefficient vector of `kappa`, then

```math
q\ge0,
\qquad
A_Eq-\lambda q=(A_Ev-\lambda v)*v\ge0.
\tag{6.3}
```

Thus `kappa` and `(tau-lambda)kappa` are positive definite, while `kappa` is
pointwise nonnegative.  Transversality gives

```math
T_a=|\mathcal C_n^+|\,\|v\|_2^2
\quad\hbox{for every }a.
\tag{6.4}
```

Let `A_bar` be normalized adjacency of the quotient Cayley graph `G/D`, and
let `v_bar` be the push-forward of `v`.  Direct fiber summation gives

```math
{J_n\over|\mathcal C_n^+|\|v\|_2^2}
=\rho_{\rm quot}(v)-\lambda,
\qquad
\rho_{\rm quot}(v)
={\langle v_{\rm bar},A_{\rm bar}v_{\rm bar}\rangle
  \over\|v\|_2^2}.
\tag{6.5}
```

Theorem 3.1 would therefore give the optimizer-free bound

```math
\mu(a)\ge2\lambda-\rho_{\rm quot}(v).
\tag{6.6}
```

The quotient leakage in (6.5) measures edges that exist between selected
quotient states but are not aligned by the chosen section.  This is a real
compressed state; it never queries a parent signing.

It nevertheless has a scalable obstruction.  A partial transversal has at
most `2^n` support points.  Theorem 4 of Bollobás--Lee--Letzter,
[*Eigenvalues of subgraphs of the cube*](https://arxiv.org/abs/1605.06360),
states that Hamming balls of radius `i=o(E)` asymptotically maximize the
spectral radius among `O(|B_E(i)|)`-vertex cube subgraphs.  Taking
`i=ceil(2n/log n)` makes `|B_E(i)|>=2^n` for large `n`, and their Theorem 3
gives

```math
\rho(Q_E[\operatorname{supp}v])
=O\!\left(\sqrt{iE}\right)
=O\!\left({n^{3/2}\over\sqrt{\log n}}\right).
\tag{6.7}
```

Since (6.1) implies
`lambda<=rho(Q_E[supp v])/E`,

```math
\boxed{
\lambda=O\!\left({1\over\sqrt{n\log n}}\right)
=o(n^{-1/2}).}
\tag{6.8}
```

Thus no rank-one partial-transversal autocorrelation can reach the leading
scale, even with zero quotient leakage.  Inverting the same estimate shows
that a correct-scale cube construction needs hidden support
`exp(Omega(n log n))`.  This does **not** rule out the paper's high-rank
moving projections or a nonabelian `S_n`-adapted construction; it explains
quantitatively why bounded or `exp(O(n))` scalar support is insufficient.

## 7. Exact sufficient theorem and falsification criterion

A genuinely matrix-valued moving-representation route would prove
convergence if it supplied, independently of `M_n`, pointwise-nonnegative
kernels `kappa_n` and complete Gram remainders

```math
F_n=(\tau-\lambda_n)\kappa_n\succeq_{\rm PD}0
\tag{7.1}
```

such that

```math
\lambda_n={1-o(1)\over\sqrt n},
\qquad
\sup_a{J_n\over T_a}=o(n^{-1/2}).
\tag{7.2}
```

Then (3.3)--(5.4) give uniformly

```math
Q(a)\ge\left({1\over2}-o(1)\right)n^{3/2},
\tag{7.3}
```

and the known conference upper bound forces
`M_n/n^(3/2)->1/2`.

The remaining lemma must be more than a repackaged coset histogram.  A valid
candidate must provide an `S_n`-equivariant algebraic lower bound for `T_a`,
or an annular pointwise lower bound for `kappa_n`, using the complete
matrix-valued remainder.  Reject a family if a scalable signing makes
`J_n/T_a` exceed `lambda_n` at fixed `n^(-1/2)` scale, if `T_a` collapses,
or if its compressed state is invertibly equivalent to the signed Eulerian
histogram.

The published whole-cube rank-one family and the scalar-transversal family
are now inactive.  A cut-specific, high-rank moving representation remains
open.
