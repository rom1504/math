# Independent audit of the conference quartic basin reduction

**Verdict: PASS, with two nonfatal scope qualifications.**  This audit is
against the following frozen task-local inputs:

```text
draft   03158cff923d7630b80e9bdfffa316aa7b37640e2730e8d29f92c53e3e4f63f3
script  116f4d6b175160ecfe8dffdf5180c1bbc5353e7d6c79148f4706b5f4ad2f8b23
result  6f5cac6d45b9ccac76026da66f707be1a97773177eca2429e15d70f7bf359723
```

The theorem-level claims QR.1--QR.5 are correct at those hashes.  In
particular, the pressure-to-quartic implication does not silently assert an
unconditioned bridge theorem: it is conditional on the explicitly stated
uniform operator and fixed-power delocalization hypotheses.

Two wording qualifications should be retained if the result is promoted.

1. QR.2 proves a speed-`r^2` **upper** bound for both orientations.  The
   matching elementary lower bound obtained by assigning mass `2^{-r^2}` to
   `B=A+I` is supplied only for `epsilon=-1`.  Thus the introductory phrase
   “has probability `exp(-Theta(r^2))`” is unconditional for that orientation,
   but the draft has not proved nonemptiness, and hence a matching lower
   bound, for every `epsilon=+1` instance.
2. The KL corollary requested after the draft was frozen is valid with an
   additive constant, and therefore as `c r^2` for all sufficiently large
   `r`.  Without adjusting constants separately, it should not be stated as a
   strictly positive `c r^2` lower bound at every small order.

## 1. Algebra and normalization audit

Put `D=BB^T` and `K_epsilon=AB+epsilon BA`.  Direct block multiplication
gives

```math
S^2=
\begin{pmatrix}
(r-1)I+D&K_\epsilon\\
K_\epsilon^T&(r-1)I+B^TB
\end{pmatrix}.
```

Using `Tr D=r^2` and
`||BB^T||_F=||B^TB||_F` gives exactly

```math
\operatorname{Tr}S^4
=2r(r-1)^2+4r^2(r-1)+2\|D\|_F^2+2\|K_\epsilon\|_F^2
=6r^3-8r^2+2r+2J_\epsilon(B).
```

For a complete hollow signing of order `N`, the repeated-edge/four-cycle
expansion independently gives

```math
\kappa_4\!\left(\sum_{i<j}s_{ij}z_iz_j\right)
=3\operatorname{Tr}S^4-2N(N-1)(3N-4).
```

Substitution of `N=2r` produces

```math
\kappa_4(H)=6J_\epsilon(B)-30r^3+32r^2-10r,
```

and multiplying by
`t^4/24=beta^4/(96r^2)` verifies every coefficient in QR.10.
The parent-minus-two-children second-order term is `beta^2/4`.  The child
conference cumulant is

```math
-3r^3+8r^2-5r,
```

so subtraction verifies QR.10a, including the constants `-1/3` and
`5/(16r)`.

The exact expectation checks in QR.17b are also correct:

```math
\mathbb E\|BB^T\|_F^2=2r^3-r^2,
\qquad
\mathbb E\|AB+\epsilon BA\|_F^2=2r^2(r-1).
```

The cross term in the second identity has expectation zero because `A` is
hollow.

## 2. Gauge invariance and the equality case

For independently switched children,

```math
(A,C,B)\mapsto(DAD,ECE,DBE)
```

sends

```math
BB^T\mapsto D(BB^T)D,
\qquad
AB+\epsilon BC\mapsto D(AB+\epsilon BC)E.
```

Both Frobenius norms are therefore invariant.  This verifies that `J` is
information on the switching quotient rather than a choice of child gauge.

Since `BB^T` is positive semidefinite with trace `r^2`,

```math
\|BB^T\|_F^2\ge r^{-1}(\operatorname{Tr}BB^T)^2=r^3.
```

Equality in `J>=r^3` forces both `BB^T=rI` and
`AB+epsilon BA=0`.  Because `B` is square, the first identity also gives
`B^TB=rI`; substituting in the displayed block square proves equivalence with
`S^2=(2r-1)I`.  Conversely that conference identity gives both conditions.
There is no missing equality branch.

For `epsilon=-1`, `B=A+I` has zero intertwiner and

```math
\|(A+I)(A+I)^T\|_F^2=r^3+4r(r-1),
```

as stated.

## 3. Projection and small-ball audit

A symmetric conference matrix has eigenvalues
`+sqrt(r-1)` and `-sqrt(r-1)`.  Its zero trace makes both multiplicities
`r/2`.  In the tensor eigenbasis, the coefficient of the `(i,j)` entry under
`B -> AB+epsilon BA` is `lambda_i+epsilon lambda_j`.  Exactly `r^2/2`
coordinates are active, and every active singular value is
`2sqrt(r-1)`.  Hence QR.20--QR.21 hold with an orthogonal projection of rank
`d=r^2/2`.

If `J<=(1+delta)r^3`, the Gram lower bound leaves

```math
b^TP_\epsilon b\le {\delta r^3\over4(r-1)}.
```

For `delta<=1/2` and every admissible `r>=2`, the right side is at most
`r^2/4=d/2`.  Rademacher `b` has identity covariance, so
`E b^TPb=Tr P=d`, while `||P||_F^2=d` and `||P||_op=1`.  Theorem 1.1 of
Rudelson--Vershynin therefore yields

```math
\Pr\{b^TPb\le d/2\}
\le2\exp\{-c\min(d/4,d/2)\}
\le2e^{-c'r^2}.
```

This exactly verifies QR.2; no independence in the spectral basis is being
assumed.

## 4. Uniform free-cumulant remainder

Moment/free-cumulant inversion gives

```math
r_j(\nu)=\sum_{\pi\in NC(j)}\mu(\pi,1_j)m_\pi(\nu).
```

The bounds `|NC(j)|<=4^j`,
`|mu(pi,1_j)|<=4^j`, and `|m_pi|<=K^j` justify the deliberately crude
`|r_j(nu)|<=(16K)^j`.  Thus for `16K beta<=1/2`, termwise integration is
uniformly valid and

```math
\mathfrak f(\beta\nu)
=\sum_{j\ge1}{\beta^jr_j(\nu)\over2j}.
```

Averaging the two signs retains the `j=2,4` terms with coefficients `1/4`
and `1/8`.  The remaining even tail is bounded by

```math
\sum_{\substack{j\ge6\\j\ {\mathrm{even}}}}
{(16K\beta)^j\over2j}
\le{(16K\beta)^6\over6}.
```

Hence both the order and the explicit constant in QR.25 are valid.  Strictly
speaking, `mathfrak f(nu)` for an arbitrary unscaled compact law need not be
represented by an `R`-series all the way to `u=1`; the lemma only uses the
small scaled laws `+-beta nu`, for which its proof establishes the required
domain.  This does not affect QR.3 or QR.4.

## 5. Primary-source audit of the FMW transfer

The primary source was downloaded independently from
[Fan--Misiakiewicz--Wang--Wen, arXiv:2607.10102](https://arxiv.org/abs/2607.10102)
(PDF SHA-256
`9cf559df9fd416c53996a5aaaa5845fe083b555bd78789699d5c282ba729b452`).
Its Assumption 2.9 requires exactly:

- a deterministic symmetric matrix with a uniform operator bound;
- weak convergence of its empirical spectral law to a centered,
  nondegenerate compact law; and
- for each fixed power and every positive exponent slack, the maximum
  diagonal deviation and maximum off-diagonal entry to be
  `o` in the stated `n^{-1/2+eta}` sense.

Corollary 2.10(c) imports Corollary 2.8 at zero field when the uniform
operator bound is strictly below `1/2`.  Applying it to `X=+-beta Y_r` is
legitimate because `beta K<1/2`; scaling preserves QR.29.  At zero external
field global spin-flip symmetry gives every one-spin Gibbs mean zero at every
finite order, so the two-replica parameter in Corollary 2.8 is exactly
`q_*=0`.  With the paper's normalization this leaves

```math
{1\over n}\log\sum_x e^{x^TXx/2}
\longrightarrow
\log2+{1\over2}\int_0^1R_\nu(u)\,du.
```

After normalizing by `2^{-n}`, this is precisely the draft's
`mathfrak f(nu)`.  The identity

```math
\overline Z_{\rm cosh}(X)
={Z(X)+Z(-X)\over2^{n+1}}
```

and arithmetic--geometric mean verify QR.34.  No equality between the two
orientation limits is assumed for the parent.

For the child, `A_r/sqrt r` has limiting Bernoulli law on `{-1,+1}`, and
the same theorem gives `T_r/(2r)=psi(beta)+o(1)`.  QR.33 follows directly by
dividing QR.7 by `8r^3`.  Its limiting fourth free cumulant is `(j-5)/4`,
which verifies the coefficient `(j-5)/32` in QR.35 and all constants in
QR.36--QR.37.

## 6. Compactness and the finite regular class

The compactness step is sound.  The bound `||Y_r||_op<=K` puts every
empirical law on `[-K,K]`; QR.33 then also bounds `J/r^3`.  A subsequence can
therefore simultaneously realize the limsup of `J/r^3` and a weak spectral
limit.  Compact support upgrades weak convergence to convergence of the
second and fourth moments, so the selected value `j` is retained in QR.35.

For `D_r(K)`, any selected sequence is FMW-regular: a fixed `p` is eventually
at most `q_r`, and for each fixed `eta>0`, eventually `eta_r<eta`.  If QR.38
failed infinitely often, selecting those failures would satisfy QR.31 and
QR.29 while keeping `J/r^3>3/2`, contradicting QR.4's bound at most `5/4`.
Thus the sequential argument really does yield an eventual uniform
containment; it is not merely a pointwise subsequence statement.

## 7. Reproducibility

The frozen script was rerun with the documented seed and sample count, with
output written under `/home/math/quadra/tmp`.  The rerun was byte-identical to
the checked-in result and had the same SHA-256
`6f5cac6d45b9ccac76026da66f707be1a97773177eca2429e15d70f7bf359723`.
The script exhausts all order-two bridges and checks 64 seeded order-six
bridges in both orientations.  It verifies QR.7, QR.9, and QR.21; as the draft
states, this is regression evidence, not support for the asymptotic step.

## 8. Entropic cost corollary

The proposed corollary is valid in the following precise form.  Fix an
orientation, let `U` be uniform on all sign bridges, and let `q` be any bridge
law.  If `0<delta<=1/4` and

```math
\mathbb E_qJ_\epsilon(B)\le(1+\delta)r^3,
```

then, writing

```math
E=\{J_\epsilon(B)\le(1+2\delta)r^3\},
```

the deterministic lower bound `J>=r^3` and Markov's inequality give
`q(E)>=1/2`.  QR.2, now with parameter `2delta<=1/2`, gives
`U(E)<=2e^{-c_0r^2}`.  Binary data processing yields

```math
\begin{aligned}
D(q\|U)
&\ge d_{\rm bin}(q(E)\|U(E))\\
&\ge {1\over2}\log{1\over U(E)}-\log2\\
&\ge {c_0\over2}r^2-{3\over2}\log2.
\end{aligned}
```

Consequently `D(q||U)>=c_1r^2` for all sufficiently large `r`, with a
universal `c_1>0` (in fact independent of `delta` throughout the stated
range).  The same proof works for a joint law on orientation and bridge if
the reference orientation is uniform, because QR.2 holds in each slice.
This corollary is a genuine entropy consequence of QR.2 and does not require
the pressure regularity theorem.
