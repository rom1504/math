# Hypercube flatness of exposed spherical optimizers: a focused literature scout

**Status.** Source-grounded negative assessment, one precise imported theorem,
one derived random-sign consequence, and one exact deterministic falsifier.  This
note does not propose a new framework.

## 1. Decision

No theorem found in the audited Kashin--Dvoretzky, discrepancy/vector
balancing, inverse Littlewood--Offord, eigenvector-universality, or sign-rounding
literature forces

```math
\phi(u):=1-\frac{\|u\|_1}{n}=o(1),
\qquad \|u\|_2^2=n,                                  \tag{HF.1}
```

for an exposed spherical optimizer of an **adversarial** symmetric hollow sign
matrix.  In fact:

1. a completely generic spherical vector has

   ```math
   \phi(u)\longrightarrow
   1-\sqrt{\frac2\pi}=0.202115\ldots;                \tag{HF.2}
   ```

2. eigenvector universality transfers precisely this Gaussian constant to an
   iid symmetric sign-Wigner model; and
3. an explicit deterministic hollow sign matrix has a simple exposed top
   eigenvector with

   ```math
   \phi(u)\longrightarrow1-\cos(\pi/8)
   =0.076120\ldots.                                  \tag{HF.3}
   ```

Thus delocalization, isotropy, and almost-Euclidean `l_1` sections point in the
wrong quantitative direction: they give a constant distance from the cube, not
`o(1)` distance.  A positive theorem would have to use special algebraic or
near-minimizer structure that synchronizes one exposed optimizer with the cube.

The exact geometric identity behind this distinction is

```math
\min_{x\in\{+-1\}^n}\frac{\|u-x\|_2^2}{n}
=2\left(1-\frac{\|u\|_1}{n}\right)=2\phi(u).         \tag{HF.4}
```

No global or dependent choice of signs can improve (HF.4) in Euclidean
distance; coordinatewise `sgn(u)` is already optimal.

## 2. One precise imported theorem

The cleanest applicable import is Bourgade and Yau's eigenvector-universality
theorem, [*The Eigenvector Moment Flow and local Quantum Unique
Ergodicity*](https://arxiv.org/abs/1312.1301), Communications in Mathematical
Physics 350 (2017), 231--278,
[DOI](https://doi.org/10.1007/s00220-016-2627-6).

### Imported theorem HF.1 (Bourgade--Yau, specialized to the real case)

Let `H_N` be a real symmetric generalized Wigner matrix.  Its upper-triangular
entries are independent and centered; writing
`sigma_ij^2=E h_ij^2`, assume

```math
\sum_i\sigma_{ij}^2=1,
\qquad C^{-1}N^{-1}\le\sigma_{ij}^2\le CN^{-1},      \tag{HF.5}
```

and, for every fixed `p`,

```math
\sup_{i,j,N}\mathbb E|\sqrt N h_{ij}|^p<\infty.     \tag{HF.6}
```

Order the eigenvalues increasingly and choose an orthonormal eigenbasis
`(v_k)`.  In the cited version there is a fixed `delta>0`, depending only on
the ensemble assumptions, and the allowed index set is exactly

```math
\mathbb T_N=
[1,N^{1/4}]
\mathbin\cup[N^{1-\delta},N-N^{1-\delta}]
\mathbin\cup[N-N^{1/4},N],                          \tag{HF.7a}
```

with integer endpoints understood.  In particular it includes both spectral
edges and the stated central bulk interval, but this note makes no claim for
the two omitted mesoscopic index ranges.  Bourgade and Yau prove, uniformly
for `k in T_N` and for every fixed coordinate set `J`, that

```math
(\sqrt N v_k(j))_{j\in J}
\Longrightarrow (g_j)_{j\in J}                     \tag{HF.7}
```

in moments, modulo the irrelevant global eigenvector sign, where the `g_j` are
independent standard real Gaussians.  Their theorem is stronger: it treats
projections onto arbitrary deterministic unit vectors and proves local quantum
unique ergodicity as well.

The variance lower bound in (HF.5), including on the diagonal, is part of this
precise cited version.  The paper resolves the arbitrary global eigenvector
sign by multiplying it by an independent uniform sign; absolute-coordinate
observables below are invariant under that convention.  Therefore the
statement below is made only for the exact ensemble covered by the paper; it
is not silently advertised as a theorem for a hollow ensemble.

### Corollary HF.2 (Gaussian flatness deficit for sign-Wigner eigenvectors)

Let `S_N` be symmetric with independent Rademacher upper-triangular entries,
including the diagonal, and put `H_N=S_N/sqrt(N)`.  Fix an allowed eigenvector
index `k=k(N)`, let `v_k` have Euclidean norm one, and set

```math
u_k=\sqrt N v_k,
\qquad L_N=\frac{\|u_k\|_1}{N}.                     \tag{HF.8}
```

Then

```math
L_N\longrightarrow\sqrt{\frac2\pi}
\quad\hbox{in }L^2,
\qquad
\phi(u_k)\longrightarrow1-\sqrt{\frac2\pi}.       \tag{HF.9}
```

#### Proof

The one- and two-coordinate cases of the **uniform** convergence (HF.7),
together with moment convergence, give uniformly in distinct `i,j`

```math
\mathbb E|u_k(1)|\to\mu,
\qquad
\mathbb E|u_k(1)u_k(2)|\to\mu^2,
\qquad \mu=\mathbb E|g|=\sqrt{2/\pi}.               \tag{HF.10}
```

Moment convergence supplies uniform integrability for the absolute-value
observables.  Averaging the uniform one-coordinate limit, and using the
deterministic identity `sum_i u_k(i)^2=N` for the diagonal part of the second
moment, gives

```math
\mathbb E L_N
=\frac1N\sum_i\mathbb E|u_k(i)|\to\mu,
\qquad
\mathbb E L_N^2
=\frac1{N^2}\sum_i\mathbb E u_k(i)^2
 +\frac1{N^2}\sum_{i\ne j}\mathbb E|u_k(i)u_k(j)|
\to\mu^2.                                           \tag{HF.11}
```

This proves `L^2` convergence and (HF.9). `square`

**Source audit: PASS after two scope repairs.**  First, HF.2 is asserted only
for iid Rademacher entries **including the diagonal**, exactly satisfying
(HF.5)--(HF.6); no hollow-matrix extension is inferred.  Second, `k(N)` is
restricted to the explicit set (HF.7a).  The derivation uses the paper's
uniform one- and two-coordinate convergence rather than assuming an
exchangeable eigenbasis on the exceptional event of a repeated eigenvalue.
No part of this import concerns adversarial near-minimizers.

**Model where the import helps.**  The iid symmetric sign-Wigner ensemble is
the closest standard random-matrix benchmark for the signing matrices.  HF.2
settles the genericity question: even asymptotically Gaussian, isotropic,
microscopically delocalized eigenvectors are not Boolean-flat.  The theorem
does **not** address adversarially selected near-minimizers; that extra
selection is exactly where any positive result would have to enter.

For comparison, the same constant for a Haar spherical vector is elementary.
If `g_i` are iid standard Gaussians and
`u=sqrt(N)g/||g||_2`, then the strong law gives

```math
\frac{\|u\|_1}{N}
=\frac{N^{-1}\sum_i|g_i|}
       {(N^{-1}\sum_i g_i^2)^{1/2}}
\longrightarrow\sqrt{\frac2\pi}\quad\hbox{a.s.}   \tag{HF.12}
```

Thus (HF.2) is the baseline for generic directions, not a defect of a weak
delocalization estimate.

## 3. Exact hollow-sign falsifier

The following elementary family rules out an all-adversarial theorem without
any random-matrix hypothesis.

Let `n=2m`.  Split the coordinates into two sets of size `m` and define the
symmetric hollow signing

```math
A_m=
\begin{pmatrix}
J_m-I_m & J_m\\
J_m & -(J_m-I_m)
\end{pmatrix}.                                      \tag{HF.13}
```

On the block-constant subspace, in the orthonormal basis of normalized block
indicators, `A_m` is

```math
K_m=
\begin{pmatrix}
m-1&m\\ m&-(m-1)
\end{pmatrix},
\qquad
\operatorname{spec}(K_m)=\{+-r_m\},
\quad r_m=\sqrt{(m-1)^2+m^2}.                       \tag{HF.14}
```

The remaining eigenvalues are `-1` and `+1`, each with multiplicity `m-1`.
Thus the top eigenvalue `r_m` is simple.  Its norm-`sqrt(n)` eigenvector is
constant with values `alpha_m,beta_m>0` on the two blocks, where

```math
t_m:=\frac{\beta_m}{\alpha_m}
=\frac{r_m-(m-1)}m,
\qquad
\alpha_m=\sqrt{\frac2{1+t_m^2}},
\qquad \beta_m=t_m\alpha_m.                         \tag{HF.15}
```

It is therefore an exposed optimizer of the positive spherical quadratic
channel, but

```math
\frac{\|u_m\|_1}{2m}
=\frac{1+t_m}{\sqrt{2(1+t_m^2)}}
\longrightarrow \cos(\pi/8)<1.                     \tag{HF.16}
```

This proves (HF.3).  The bottom channel has the same limiting absolute
coordinate profile.

The finite spectrum, eigenvector residual, formula (HF.16), its limit, and the
Gaussian baseline (HF.12) are checked reproducibly by
[`verify_hypercube_flatness_scout.py`](../experiments/verify_hypercube_flatness_scout.py).

**Scope of the falsifier.**  This family has a quadratic extremum of order
`n^2`, so it is not a near-minimizer for the original `n^(3/2)` signing
problem.  It kills only the unrestricted claim

```math
\text{``sign matrix + exposed eigenvector''}
\Longrightarrow \phi=o(1).                         \tag{HF.17}
```

It does not kill a theorem whose hypotheses use a verified rigidity property
of near-minimizers.

## 4. What the other audited literatures actually provide

| Source/mechanism | Precise output relevant here | Why it does not yield `phi=o(1)` |
|---|---|---|
| Milman's concentration proof of Dvoretzky ([1971 paper](https://doi.org/10.1007/BF01086740)); Kashin/Garnaev--Gluskin almost-Euclidean sections ([primary 1984 source](https://www.mathnet.ru/eng/dan9506)) | On suitable subspaces, `||.||_1` is uniformly comparable, or nearly proportional after normalization, to `||.||_2`.  For the `l_1^N` norm, the spherical concentration center is `E||theta||_1 ~ sqrt(2/pi)sqrt(N)`. | Near-Euclideanity controls variation about the Gaussian center.  At norm `sqrt(N)`, that center is `||u||_1/N~sqrt(2/pi)`, a fixed distance below one.  Kashin's constant-equivalence version is weaker still. |
| Banaszczyk vector balancing ([1998 primary paper](https://doi.org/10.1002/%28SICI%291098-2418%28199807%2912%3A4%3C351::AID-RSA3%3E3.0.CO;2-S)); constructive Gram--Schmidt walk ([Bansal--Dadush--Garg--Lovett](https://theoryofcomputing.org/articles/v015a021/)) | Given short vectors and a convex body of Gaussian measure at least `1/2`, choose signs whose **linear image** lies in a constant dilation of the body. | This balances `sum epsilon_i v_i`; it does not approximate a prescribed spherical vector by a coordinate sign vector.  For that latter task (HF.4) is already the exact optimum. |
| Tao--Vu inverse Littlewood--Offord ([Annals 2009](https://annals.math.princeton.edu/2009/169-2/p06), [sharp version](https://arxiv.org/abs/0902.2357)) | Abnormally high concentration of `sum epsilon_i a_i` forces most coefficients `a_i` into a low-rank generalized arithmetic progression. | Exposure supplies no requisite small-ball concentration.  Additive structure of the coordinate multiset does not force `|a_i|=1+o(1)`; even the inverse conclusion is about progression structure, not cube distance. |
| Random hyperplane/sign embeddings, e.g. Plan--Vershynin ([primary preprint](https://arxiv.org/abs/1111.4452)) | Encode directions by signs of random measurements while approximately preserving angular/Hamming distance on a prescribed set. | The generated signs live in measurement space.  This is an embedding theorem, not an assertion that the original coordinate vector lies near `{+-1}^N`. |
| Eigenvector delocalization/universality | Individual coordinates are small and, under the hypotheses of HF.1, jointly Gaussian after `sqrt(N)` scaling. | Gaussian coordinates have mean absolute value `sqrt(2/pi)`, giving the fixed deficit (HF.9).  Bounds such as `||v||_infty=O(sqrt(log N/N))` are much weaker and cannot imply equal coordinate magnitudes. |

There is therefore no hidden upgrade in these sources from “mass is spread over
many coordinates” to “almost every coordinate has magnitude one.”  The latter
is a rigidity statement at the equality case of Cauchy--Schwarz.

## 5. The exact missing statement and stopping judgment

A theorem strong enough for the exposed-rounding carrier would need a
near-minimizer-specific hypothesis, for example the following schematic form:

> For every sequence of signing instances in the required compositional class
> with normalized Boolean cap within `o(1)` of the optimum, and every bounded
> admissible port mass, at least one `o(1)`-near-optimal spherical response
> channel has an exposed optimizer `u_N` satisfying `phi(u_N)=o(1)`.

This statement is genuinely stronger than delocalization and genuinely weaker
than storing the full Boolean response table: it asks for one checkable exposed
witness.  None of the imported theorems supplies its near-minimizer selection
step.  The deterministic family (HF.13) shows that some such structural
hypothesis is logically necessary.

**Recommendation.**  Do not pursue generic Kashin/Dvoretzky, discrepancy,
inverse Littlewood--Offord, or Wigner-delocalization arguments as a way to prove
exposed Boolean flatness.  Keep the carrier only for algebraically synchronized
or rigorously characterized near-minimizer classes.  A future positive route
must prove an equality-case rigidity/synchronization theorem, not another form
of ordinary delocalization.
