# Banach-space / harmonic-analysis audit for the quadratic signing problem

## 1. Exact translation

Let
\[
\Lambda_{2,n}=\{\chi_{\{i,j\}}(x)=x_ix_j:1\le i<j\le n\},
\qquad L_n=|\Lambda_{2,n}|=\binom n2 .
\]
Define the equal-modulus, or unimodular, Sidon constant
\[
U_n:=
\max_{\varepsilon_{ij}\in\{\pm1\}}
\frac{L_n}
{\left\|\sum_{i<j}\varepsilon_{ij}\chi_{\{i,j\}}\right\|_\infty}.
\]
Then exactly
\[
\boxed{U_n=\frac{\binom n2}{F_n}},\qquad
\frac{F_n}{n^{3/2}}
=\frac{1-1/n}{2}\frac{\sqrt n}{U_n}.
\tag{1}
\]
Thus the original limit exists if and only if \(U_n/\sqrt n\) converges.

The ordinary Sidon constant is
\[
\operatorname{Sid}(\Lambda_{2,n})
=\sup_{a\ne0}
\frac{\sum_{i<j}|a_{ij}|}
{\left\|\sum_{i<j}a_{ij}\chi_{\{i,j\}}\right\|_\infty}.
\tag{2}
\]
Consequently
\[
U_n\le \operatorname{Sid}(\Lambda_{2,n}),
\tag{3}
\]
but (2) allows arbitrary coefficient magnitudes whereas (1) forces every
magnitude to be one.  This distinction is precisely where the available
Banach-space theorems stop short of the problem.

## 2. Closest exact theorem: projection constants

Defant, Galicer, Mansilla, Mastyło and Muro,
[“Asymptotic insights for projection, Gordon–Lewis and Sidon constants in
Boolean cube function spaces”](https://arxiv.org/abs/2302.00233)
(IMRN 2024), prove
\[
\lambda(\mathcal B^n_{=2})
=2^{-n}\sum_{x\in\{\pm1\}^n}
\left|\sum_{i<j}x_ix_j\right|
\]
and
\[
\boxed{\lim_{n\to\infty}\frac{\lambda(\mathcal B^n_{=2})}{n}
=\sqrt{\frac{2}{\pi e}}
=0.4839414490\ldots .}
\tag{4}
\]
Indeed,
\[
\sum_{i<j}x_ix_j
=\frac12\left[\left(\sum_i x_i\right)^2-n\right],
\]
so (4) is a central-limit calculation.

The numerical proximity of (4) to the ROM candidate
\(\sqrt{15}/8=0.4841229183\ldots\) is striking but does not give an identity:
their difference is \(1.81469\cdot10^{-4}\), and the invariants have different
definitions and scalings.  Formula (4) is the average absolute value of the
**all-positive** Walsh kernel; \(F_n\) is a minimum over coefficient signings
of a supremum.

The same paper proves only dimension-free comparisons for Sidon constants:
for fixed \(d\),
\[
\operatorname{Sid}(\mathcal B^n_{=d})
\asymp_{C^d} \binom nd^{(d-1)/(2d)}
\]
(equivalently, for \(d=2\),
\(\operatorname{Sid}(\Lambda_{2,n})\asymp\sqrt n\)).
Theorem 5.1 bounds the Sidon constant by a degree-dependent multiple of the
projection constant one level lower.  The multiplier is fixed away from one;
neither the theorem nor its proof yields an asymptotic constant for
\(\operatorname{Sid}(\Lambda_{2,n})/\sqrt n\), much less for its equal-modulus
subproblem \(U_n/\sqrt n\).

A July 2026 extension by Defant et al.,
[“Support-Sensitive Bohnenblust–Hille Inequalities and Local Invariants on
Hamming Schemes”](https://arxiv.org/abs/2607.05594), again obtains
\[
\operatorname{Sid},\ \chi,\ \mathrm{gl}
\asymp_{c(q)^d}(n/d)^{(d-1)/2}
\]
for spherical/homogeneous/tetrahedral spaces.  The constants are
dimension-free, not \(1+o(1)\), so this does not supply stabilization.

## 3. Strongest theorem directly about the minimum over signs

Astashkin and Lykov,
[“Random unconditional convergence of Rademacher chaos in \(L_\infty\) and
sharp estimates for discrepancy of weighted graphs and
hypergraphs”](https://arxiv.org/abs/2412.20107), Theorem 3, prove that for
arbitrary edge weights \(w_{ij}\), with universal comparison constants,
\[
\begin{aligned}
&\min_{\varepsilon_{ij}=\pm1}
\left\|\sum_{i<j}\varepsilon_{ij}w_{ij}r_ir_j\right\|_\infty\\
&\quad\asymp
\mathbb E_\varepsilon
\left\|\sum_{i<j}\varepsilon_{ij}w_{ij}r_ir_j\right\|_\infty\\
&\quad\asymp
\max\left\{
\sum_{i=1}^{n-1}\left(\sum_{j=i+1}^nw_{ij}^2\right)^{1/2},
\sum_{j=2}^{n}\left(\sum_{i=1}^{j-1}w_{ij}^2\right)^{1/2}
\right\}.
\tag{5}
\end{aligned}
\]
For \(w_{ij}=1\), the last quantity is exactly
\[
\sum_{k=1}^{n-1}\sqrt{k}
=\frac23n^{3/2}+O(\sqrt n).
\tag{6}
\]
The minimum on the first line of (5) is exactly \(F_n\).  Hence (5) is the
closest published functional-analytic theorem to the problem.

It still proves only
\[
c\sum_{k<n}\sqrt k\le F_n\le C\sum_{k<n}\sqrt k
\]
with fixed universal \(c,C\).  In the displayed proof, an explicit lower
constant is already degraded by decoupling/polarization to
\(1/(16\sqrt2)\) relative to the row functional.  There is no assertion that
either comparison ratio converges or tends to one.  Even though the random
expectation in (5) has a spin-glass limit, a fixed-factor RUC comparison cannot
transfer that limit to the minimum.

## 4. Bohnenblust–Hille and Kahane–Salem–Zygmund

Defant, Mastyło and Pérez,
[“On the Fourier spectrum of functions on Boolean
cubes”](https://doi.org/10.1007/s00208-018-1756-y), prove the dimension-free
Boolean Bohnenblust–Hille inequality
\[
\left(\sum_{|S|=d}|\widehat f(S)|^{2d/(d+1)}\right)^{(d+1)/(2d)}
\le C_d\|f\|_\infty .
\tag{7}
\]
For \(d=2\) and \(|a_{ij}|=1\), (7) gives only a constant-times
\(n^{3/2}\) lower bound.  Its constant is independent of \(n\), but it is not
sharp in this equal-modulus quadratic problem and contains no cross-dimension
relation.

The Boolean Kahane–Salem–Zygmund construction gives a signing with an
\(O(n^{3/2})\) norm.  Again it is an existence theorem at each dimension, not
a compatible family or a submultiplicative inequality.

For genuinely multilinear forms with independent input vectors,
Albuquerque and Rezende,
[“Asymptotic estimates for unimodular multilinear forms on sequence spaces
with small norms”](https://arxiv.org/abs/1710.09711), prove the sharp order
\[
\inf\|A\|\asymp
\left(\sqrt{n_1}+\cdots+\sqrt{n_d}\right)
\prod_jn_j^{1/2-1/p_j}.
\tag{8}
\]
Pellegrino and Raposo,
[“Constants of the Kahane–Salem–Zygmund inequality asymptotically bounded by
\(1\)”](https://arxiv.org/abs/2006.12892), sharpen the upper KSZ constant in
the large rectangular multilinear regime.

These do not solve the present problem.  In degree two their norm is
\[
\max_{x,y\in\{\pm1\}^n}|x^\top Ay|,
\]
whereas ours is
\[
\max_{x\in\{\pm1\}^n}|x^\top Ax|/2
\]
with \(A\) symmetric and zero diagonal.  Decoupling or polarization compares
the two only up to a fixed constant.  That fixed loss is at the leading
\(n^{3/2}\) scale, so a multilinear limit or a KSZ constant tending to one
does not imply the desired diagonal limit.

## 5. Why Gordon–Lewis / summing factorization does not close the gap

The Gordon–Lewis and \(1\)-summing machinery is useful because those constants
have an ideal property under factorizations.  The 2024 and 2026 papers exploit
it to compare:

- the ordinary Sidon constant;
- the unconditional basis constant over **all** coefficient vectors; and
- the Gordon–Lewis constant of the whole polynomial space.

The equal-modulus set
\[
\{(a_{ij}):|a_{ij}|=1\}
\]
is nonlinear and is not preserved by the averaging, projection, or
factorization maps in those proofs.  The output coefficient vector generally
has unequal magnitudes.  Therefore the ideal property does not produce a
map between \(U_n\) and \(U_m\), nor a Fekete-type inequality.

Permutation transitivity of the weight-two Walsh set also does not justify
flattening: symmetrizing a coefficient vector averages its **sign pattern**
as well as its magnitudes and can collapse the polynomial.  No cited theorem
shows
\[
\operatorname{Sid}(\Lambda_{2,n})\le(1+o(1))U_n.
\tag{9}
\]

## 6. Verdict and precise missing theorem

No theorem found in this literature implies convergence of
\(F_n/n^{3/2}\).  The strongest exact reformulation is (1), and the strongest
direct general theorem is the RUC equivalence (5).

The Banach-space route would become decisive if one could prove both:

1. **asymptotic flattening**
   \[
   \operatorname{Sid}(\Lambda_{2,n})=(1+o(1))U_n;
   \]
2. **stabilization of the ordinary invariant**
   \[
   \operatorname{Sid}(\Lambda_{2,n})/\sqrt n
   \quad\text{has a limit}.
   \]

Neither statement is present in the primary sources checked.  Existing
projection, Gordon–Lewis, \(p\)-summing, BH, KSZ, RUC, tensor, and decoupling
results all lose a fixed multiplicative constant or leave the equal-modulus
slice.  Those losses are exactly of leading order here.

