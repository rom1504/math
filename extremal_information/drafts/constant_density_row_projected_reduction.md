# Constant-density row fibres reduce to one sharp projected-edge lemma

**Status.** Task-local theorem report.  This continues the projected-coupling
campaign for an arbitrary centrally symmetric row event.  It proves the
low-rank covariance peel, its nuclear cost, and the required Frobenius
transport.  It does **not** prove the remaining sharp operator edge.

## 1. Setup

Let `U_r` be uniform on `{+-1}^r`, let `E_r=-E_r` satisfy

```math
p_r=U_r(E_r)\ge p_0>0,
\tag{CR.1}
```

and let `mu_r=U_r(.|E_r)`.  Write `R~mu_r`, `W~U_r`, and

```math
\Sigma_r=\mathbb E RR^T.
\tag{CR.2}
```

Central symmetry gives `E R=0`.  No permutation symmetry, magnitude
description, or bounded-degree definition of `E_r` is assumed.

## 2. Fourier information forces a low-rank covariance peel

### Lemma CR.1 (constant density has only bounded total covariance defect)

Uniformly in `r`,

```math
\boxed{
\|\Sigma_r-I\|_F^2\le2(p_0^{-1}-1).}
\tag{CR.3}
```

**Proof.**  Let `g=1_E/p_r` be the density of `mu_r` relative to `U_r`.
For every nonempty coordinate set `S`, write

```math
m_S=\mathbb E_{\mu_r}\prod_{j\in S}R_j
=\langle g,\chi_S\rangle_{U_r}.
```

Parseval gives

```math
\sum_{S\ne\varnothing}m_S^2
=\mathbb E_Ug^2-1=p_r^{-1}-1.
\tag{CR.4}
```

The diagonal of `Sigma_r` is one, while its `(i,j)` entry for `i!=j` is
`m_{\{i,j\}}`.  Summing the ordered off-diagonal entries and discarding all
other Fourier levels proves (CR.3).  `square`

For any sequence `delta_r downarrow0`, let `P_r` be the spectral projection
of `Sigma_r` onto eigenvalues outside `[1-delta_r,1+delta_r]`.  Lemma CR.1
gives

```math
\boxed{
\operatorname{rank}P_r
\le {2(p_0^{-1}-1)\over\delta_r^2},
\qquad
\|(I-P_r)\Sigma_r(I-P_r)\|_{op}\le1+\delta_r.}
\tag{CR.5}
```

Taking, for example, `delta_r=r^(-1/4)` makes the rank `O_(p_0)(sqrt r)=o(r)`.

## 3. The removed response subspace has negligible nuclear price

Let `B_r` have `r` independent rows with law `mu_r`.  Conditioning and
`p_r>=p_0` imply for every deterministic projection `P` of rank `k` that

```math
\mathbb E\|PR\|_2^2
\le p_0^{-1}\mathbb E\|PW\|_2^2={k\over p_0}.
\tag{CR.6}
```

Consequently

```math
\begin{aligned}
\mathbb E\|B_rP\|_*
&\le\sqrt{k}\,\mathbb E\|B_rP\|_F\\
&\le\sqrt{k}\left(r\mathbb E\|PR\|_2^2\right)^{1/2}
\le k\sqrt{r/p_0}.
\end{aligned}
\tag{CR.7}
```

For the projection in (CR.5), this is `O_(p_0)(r)=o(r^(3/2))`.
The iid bridge obeys the same estimate with `p_0=1`.  Thus the nuclear-cost
hypothesis of the projected-coupling theorem is automatic.

## 4. Entropy transport gives the Frobenius coupling

### Lemma CR.2 (a constant-density row event is `O(sqrt r)` edits away)

There is a coupling `(R,W)` with

```math
\boxed{\mathbb E d_H(R,W)\le C_{p_0}\sqrt r.}
\tag{CR.8}
```

**Proof.**  Every real function `phi` which is one-Lipschitz for Hamming
distance satisfies the bounded-difference MGF inequality

```math
\log\mathbb E_Ue^{s(\phi-\mathbb E_U\phi)}\le {s^2r\over8}.
\tag{CR.9}
```

Entropy duality and optimization in `s` therefore give

```math
\mathbb E_{\mu_r}\phi-\mathbb E_{U_r}\phi
\le\sqrt{{r\over2}D(\mu_r\|U_r)}
\le\sqrt{{r\over2}\log(1/p_0)}.
\tag{CR.10}
```

Apply the same estimate to `-phi` and use Kantorovich--Rubinstein duality on
the finite Hamming cube.  The resulting `W_1` bound is exactly the infimum
of `E d_H(R,W)` over couplings, proving (CR.8).  `square`

Couple the rows independently and call the resulting matrices `B_r,W_r`.
Every changed sign contributes four to squared Frobenius distance, so

```math
\mathbb E\|B_r-W_r\|_F^2
=4r\mathbb E d_H(R,W)=O_{p_0}(r^{3/2}).
\tag{CR.11}
```

Projection is a Frobenius contraction, and hence

```math
\boxed{
\mathbb E\|(B_r-W_r)(I-P_r)\|_F=O_{p_0}(r^{3/4})=o(r).}
\tag{CR.12}
```

Thus the Frobenius-transport hypothesis of the projected-coupling theorem is
also automatic.

## 5. Exact remaining lemma

Combining (CR.5), (CR.7), and (CR.12) with the audited projected-coupling
criterion leaves only the following statement.

### `L_edge` (sharp edge after the Fourier covariance peel)

For every fixed `p_0>0`, every centrally symmetric `E_r` satisfying
(CR.1), and the projection `P_r` from (CR.5) with, say,
`delta_r=r^(-1/4)`, prove that

```math
\boxed{
\|B_r(I-P_r)\|_{op}\le(2+o_\Pr(1))\sqrt r.}
\tag{CR.13}
```

For the pressure theorem it is enough, at each fixed
`beta<sqrt(2)/6`, to replace `o(1)` by any constant strictly below
`1/(sqrt(2)beta)-3`.

If `L_edge` holds, every constant-density centrally symmetric row-product
law has

```math
\mathbb E[(h_\beta-f(B_r)/r)_+]\longrightarrow0
\tag{CR.14}
```

throughout the strict conference interval.  This implication uses no
pressure information in the definition of `P_r`.

The covariance statement alone does not prove (CR.13).  The rows have
bounded likelihood ratio relative to Rademacher measure and are uniformly
subgaussian in every linear direction, but their coordinates can have
arbitrary high-order dependence.  Standard nonasymptotic independent-row
theorems give a constant multiple of `sqrt r`, not the sharp coefficient
two.  Therefore importing a Bai--Yin edge without verifying its dependence
hypotheses would be invalid.

## 6. Frontier movement

The previous class-level question asked simultaneously for a low-rank peel,
a low-cost coupling, and a sharp regular bulk.  Lemmas CR.1--CR.2 prove the
first two uniformly over the entire constant-density row class.  The SML is
now the standalone dependent-row random-matrix statement (CR.13).

A counterexample to `L_edge` would be valuable only if its excess top
singular direction remains after the covariance peel.  A covariance spike,
a row magnitude, or finitely many exceptional directions is already removed
by construction and cannot falsify the lemma.
