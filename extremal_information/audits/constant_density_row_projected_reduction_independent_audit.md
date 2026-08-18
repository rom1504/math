# Independent audit: constant-density row projected reduction

**Frozen source:**
`extremal_information/drafts/constant_density_row_projected_reduction.md`

**SHA-256:**
`748bf8f42ad65eda569c573c6a177353e15f854cc38de062f2ec523a0f30ef67`

**Verdict:** **PASS.**  The Fourier normalization, the factor two in the
covariance defect, the spectral-rank estimate, both nuclear estimates, the
Hamming transport constants, the matrix Frobenius scaling, and the reduction
to `L_edge` are correct.  No repair to the frozen source is required.

This audit does **not** certify `L_edge`; the source states it as the sole
remaining lemma.  High-order constant-density events with no covariance
defect show that this is a genuinely stronger random-matrix assertion, not a
consequence hidden in CR.1--CR.2.

## 1. Fourier normalization and the ordered-pair factor

With normalized uniform measure on the cube and

```math
g={1_E\over p},
```

one has

```math
\widehat g(\varnothing)=1,
\qquad
\mathbb E_U g^2={p\over p^2}={1\over p}.
```

Consequently normalized Parseval gives exactly

```math
\sum_{S\ne\varnothing}\widehat g(S)^2=p^{-1}-1.
```

For `i != j`,

```math
(\Sigma-I)_{ij}=\mathbb E_\mu R_iR_j
=\widehat g(\{i,j\}),
```

whereas the diagonal vanishes.  The Frobenius sum counts `(i,j)` and `(j,i)`
separately, while the Fourier sum contains the set `{i,j}` once.  Hence

```math
\|\Sigma-I\|_F^2
=2\sum_{i<j}\widehat g(\{i,j\})^2
\le 2(p^{-1}-1)
\le 2(p_0^{-1}-1).
```

Thus the factor two in CR.3 is necessary and correct.  Central symmetry gives
`E_mu R=0`, so `Sigma` is both the second-moment matrix used in the proof and
the centered covariance matrix advertised by the text.

## 2. Spectral peel

Let `lambda_i` be the eigenvalues of `Sigma`.  Since

```math
\sum_i(\lambda_i-1)^2=\|\Sigma-I\|_F^2,
```

the number outside `[1-delta_r,1+delta_r]` is at most

```math
{2(p_0^{-1}-1)\over\delta_r^2}.
```

The spectral projection `P_r` commutes with `Sigma`; on its orthogonal
complement all eigenvalues lie in the displayed interval.  Therefore

```math
\|(I-P_r)\Sigma(I-P_r)\|_{op}\le1+\delta_r.
```

For `delta_r=r^(-1/4)`, the exceptional rank is `O_(p_0)(sqrt r)=o(r)`.
There is no missing dimensional factor here.

## 3. Conditioned nuclear cost

For every deterministic rank-`k` projection `P`, positivity of
`||P R||_2^2` and the density bound `d mu/dU <= p_0^(-1)` give

```math
\mathbb E_\mu\|PR\|_2^2
\le p_0^{-1}\mathbb E_U\|PW\|_2^2
={k\over p_0}.
```

Since `rank(BP) <= k`, Cauchy--Schwarz for singular values and then Jensen
give

```math
\begin{aligned}
\mathbb E\|BP\|_*
&\le \sqrt{k}\,\mathbb E\|BP\|_F\\
&\le \sqrt{k}\,(\mathbb E\|BP\|_F^2)^{1/2}\\
&\le \sqrt{k}\,(rk/p_0)^{1/2}
=k\sqrt{r/p_0}.
\end{aligned}
```

The same calculation for the iid bridge has `p_0=1`.  With
`k=O_(p_0)(sqrt r)`, each expectation is `O_(p_0)(r)`, which is indeed
`o(r^(3/2))`.  The projection is deterministic because it is defined from
the deterministic row law, as required by the projected-coupling theorem.

## 4. Hamming transport and Frobenius scaling

For the unnormalized Hamming metric, a one-Lipschitz function has coordinate
oscillation at most one.  Hoeffding's bounded-difference MGF therefore has
the exact variance proxy used in CR.9:

```math
\log\mathbb E_U e^{s(\phi-\mathbb E_U\phi)}
\le {s^2r\over8}.
```

The entropy variational inequality gives, for `D=D(mu||U)=log(1/p)`,

```math
\mathbb E_\mu\phi-\mathbb E_U\phi
\le {D\over s}+{sr\over8}.
```

Optimizing at `s=sqrt(8D/r)` yields `sqrt(rD/2)`.  Applying the same bound
to `-phi` and using finite-space Kantorovich--Rubinstein duality proves

```math
W_1(\mu,U)\le
\sqrt{{r\over2}\log(1/p_0)}.
```

Thus CR.8 has a valid constant depending only on `p_0`.

Taking independent copies of an optimal row coupling preserves both row
product marginals.  A sign mismatch has numerical difference two and hence
contributes four to squared Frobenius distance.  Therefore

```math
\mathbb E\|B-W\|_F^2
=4r\,\mathbb E d_H(R,W)
=O_(p_0)(r^{3/2}).
```

Right multiplication by an orthogonal projection is a Frobenius contraction,
and Jensen then gives

```math
\mathbb E\|(B-W)(I-P)\|_F
\le (\mathbb E\|B-W\|_F^2)^{1/2}
=O_(p_0)(r^{3/4})=o(r).
```

The factor four, the number `r` of rows, and the final exponent `3/4` are all
correct.

## 5. Exact implication from `L_edge`

Fix `beta<sqrt(2)/6` and put

```math
c_\beta={1\over\sqrt2\,\beta}-3>0.
```

`L_edge`, together with the standard iid Bernoulli edge, lets one choose a
fixed `delta in (0,c_beta)` and events of probability tending to one on
which

```math
\max\{\|W_r\|_{op},\|B_r(I-P_r)\|_{op}\}
\le(2+\delta)\sqrt r.
```

This is PC.5, and `beta(3+delta)/sqrt(2)<1/2` is exactly PC.4.  Sections 3
and 4 above give PC.6--PC.7.  The audited projected-coupling theorem then
implies CR.14.  Intersecting the two operator-regular events loses only
`o(1)` probability.  On their complement, nonnegativity of the cosh pressure
bounds the normalized positive shortfall, so no additional endpoint
integrability is being assumed.

More generally, if the right side of CR.13 is
`(2+c+o_Pr(1))sqrt(r)`, the same argument works whenever
`c<c_beta`: choose the fixed PC margin strictly between `c` and `c_beta`.
This verifies the source's relaxed fixed-temperature formulation.

Therefore CR.3, CR.7, CR.12, and `L_edge` imply the pressure conclusion for
the entire declared class.  No definition in the reduction uses pressure or
the child optimum.

## 6. Attempted high-order falsifiers

Two stress tests make clear both the strength and the exact limitation of
the reduction.

### Full parity

For even `r>2`, let

```math
E_r=\left\{x:\prod_{j=1}^r x_j=1\right\}.
```

This event is centrally symmetric and has density `1/2`.  Its density is
`g=1+chi_[r]`; hence every degree-one and degree-two Fourier coefficient
vanishes.  Thus `Sigma=I` and `P_r=0`, even though the row coordinates have
an exact order-`r` constraint.  All proved steps above remain valid, while
`L_edge` becomes the genuinely nontrivial assertion that a square matrix
with independent parity-conditioned rows has sharp edge two.

### Four-block majority parity

Partition the coordinates into four odd blocks and condition on the product
of the four block-majority signs being `+1`.  The event is centrally
symmetric and has density `1/2`.  Any collection of at most three whole
blocks retains its uniform product marginal under this conditioning.  In
particular every coordinate pair is uniform, so again `Sigma=I` and `P_r=0`,
despite a nonlinear macroscopic dependency involving all four blocks.

These examples do not falsify the reduction: its Fourier, transport, and
nuclear arguments deliberately make no assumption on higher Fourier levels.
They do rule out interpreting CR.3 as if it already controlled the operator
edge.  A counterexample among them, or among more complicated
constant-density fibres, would falsify `L_edge` and hence stop this route,
but it would not expose a gap in the claimed reduction to `L_edge`.

## 7. Scope of the verdict

The result is a strict and correctly stated reduction for **row-product**
laws obtained by conditioning each uniform row on the same centrally
symmetric constant-density event.  It does not cover dependencies between
different rows.  It also does not import a sharp Bai--Yin theorem under
insufficient hypotheses: the source explicitly leaves that sharp edge as
the SML.  Subject to those stated boundaries, the whole-class reduction is
rigorous.
