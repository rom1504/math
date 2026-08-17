# Independent audit: PSD compatibility fibres

**Audited files.** `psd_compatibility_fibre.md` and
`verify_psd_compatibility_fibre.py`.

## Verdict

**PASS WITH REQUIRED REPAIRS.**  Theorems PF.1--PF.3, the multi-piece PSD
correlation carrier, truncation constants, net cardinality, and response
scaling are correct.  One displayed sharpness example is false at the
allowed endpoint `tau=1`: PF.21--PF.22 require `0<tau<1`.  Two scope
clarifications are also required:

1. “operator ball modulo gauges” means a presented associated fibre; one
   must not quotient `W` alone down to its singular values, which would lose
   real cross-block information;
2. the response consequence is compressive for growing `p` only at
   `mp/r=O(1)`, i.e. `m=O(r/p)`.  At the anti-pin choice `m=r`, its metric
   accuracy and net complexity become strongly `p`-dependent.

An independent verifier adds arbitrary multi-piece recovery, a zero-rank
Douglas case, the `tau=1` endpoint diagnostic, and the net/response
arithmetic.

## 1. PF.1, including singular supports

For a PSD block matrix, positivity of the quadratic polynomial in a scalar
parameter gives

```math
|x^TCy|^2\le (x^TK_1x)(y^TK_2y).                     \tag{APF.1}
```

If `Y_1^Tx=0`, then `x^TK_1x=0`, so APF.1 makes
`x^TCy=0` for every `y`; the same holds on the second kernel.  Hence

```math
(Y_1^Tx,Y_2^Ty)\longmapsto x^TCy                    \tag{APF.2}
```

is well-defined even when either marginal matrix is singular.  Because
`Y_i` has full column rank, `Y_i^T` maps onto `R^(r_i)`.  APF.1 says that
APF.2 is a bilinear contraction on those Euclidean support spaces.  Riesz
representation therefore gives `W` with `||W||op<=1` and
`C=Y_1WY_2^T`.

Uniqueness follows by applying left inverses of both skinny factors.  The
rank-zero case is included: if `K_1=0`, APF.1 forces `C=0`, and the unique
coordinate contraction has shape `0 by r_2`.  The independent verifier
checks this edge case and confirms that any nonzero cross block makes the
joined matrix indefinite.

Conversely,

```math
\begin{pmatrix}I&W\\W^T&I\end{pmatrix}\succeq0
\quad\Longleftrightarrow\quad \|W\|_{op}\le1,        \tag{APF.3}
```

by a Schur complement or singular-value decomposition.  Congruence by the
block factor matrix proves the reverse direction.  PF.1 is therefore valid
without invertibility of `K_1` or `K_2`.

The gauge formula is also correct.  The precise intrinsic object is a
contraction between the two marginal support Hilbert spaces.  Once bases
are chosen it has `r_1r_2` coordinates, and simultaneous changes of the
factor bases act as `W -> O_1^T W O_2`.  Taking the double quotient of the
matrix ball while forgetting the factors would retain essentially only
singular values and would **not** parameterize the cross blocks.  Thus the
word “modulo” in the draft must be read as gauge equivalence of the whole
presentation `(Y_1,Y_2,W)`, not as a smaller standalone state for `W`.

## 2. Multi-piece correlation carrier

Let `D=diag(Y_1,...,Y_s)`.  Given a joined PSD matrix `K`, every vector
supported in piece `i` and lying in `ker(Y_i^T)=ker(K_i)` has zero quadratic
form, hence lies in `ker K`.  Therefore

```math
\operatorname{range}K\subseteq\operatorname{range}D. \tag{APF.4}
```

Set

```math
\Omega=D^\dagger K(D^\dagger)^T.                    \tag{APF.5}
```

It is PSD by congruence, and APF.4 gives `K=D Omega D^T`.  Its diagonal
blocks are

```math
\Omega_{ii}=Y_i^\dagger(Y_iY_i^T)(Y_i^\dagger)^T=I_{r_i}. \tag{APF.6}
```

Conversely any PSD `Omega` with these identity diagonal blocks produces a
PSD joined matrix with exactly the prescribed marginals.  This proves the
claimed equivalence, including singular marginals.  Principal restriction
and PSD completion are associative at the level stated in Section 2.

Pairwise contractions are necessary but not sufficient for `s>=3`: the
full block matrix `Omega` must be PSD.  The canonical verifier builds one
forward three-piece example.  The independent verifier additionally starts
from two arbitrary singular global Gram presentations, recovers `Omega` by
APF.5, checks every identity block, and reconstructs the joined matrix.

## 3. PF.2 truncation constants

The trace hypothesis implies that the number of eigenvalues exceeding
`tau p_i` is strictly less than `1/tau`, hence PF.13.  The discarded
marginal block has operator norm at most `tau p_i`, so for Boolean `x`

```math
x^T(K_i-K_i^h)x\le \tau p_i^2.                       \tag{APF.7}
```

Use a Douglas contraction in square-root coordinates.  Then

```math
C-P_1CP_2=Q_1C+P_1CQ_2.                              \tag{APF.8}
```

For the first term,

```math
|x^TQ_1Cy|
\le\|K_1^{1/2}Q_1x\|\,\|K_2^{1/2}y\|
\le \sqrt\tau,p_1p_2,                               \tag{APF.9}
```

and the symmetric estimate holds for the second.  The joined quadratic
counts the cross block twice, giving

```math
{\tau(p_1^2+p_2^2)+4\sqrt\tau p_1p_2\over(p_1+p_2)^2}. \tag{APF.10}
```

Writing `a=p_1/(p_1+p_2)` reduces this to

```math
\tau+(4\sqrt\tau-2\tau)a(1-a),                      \tag{APF.11}
```

whose maximum for `0<tau<=1` is
`sqrt(tau)+tau/2`.  The factor two in PF.16 then follows from
`d_q=2 max(q_p^+,q_p^-)`.  All constants are correct.

### Required endpoint repair in PF.21--PF.22

The sharpness construction works for `0<tau<1`.  Its first eigenvalue is
exactly `tau p_1` and is discarded by the strict cutoff, while the second
eigenvalue `p_2` is retained because `p_2>tau p_2`.  The lost all-positive
query is then exactly

```math
{\tau+2\sqrt\tau\over4}.                             \tag{APF.12}
```

At `tau=1`, however, the second eigenvalue is also exactly at the strict
threshold and is discarded.  Then `K^h=0` and the normalized all-positive
loss is `1`, not `3/4`.  The theorem PF.2 remains true (`1<=3/2`); only the
claimed exact example needs the explicit restriction `0<tau<1`.  The
independent verifier records this discrepancy exactly.

For `tau` tending to zero, APF.12 is `Theta(sqrt(tau))`, so the conclusion
that error-`eta` reusable truncation may require cutoff `Theta(eta^2 p_i)`
and rank `O(eta^(-2))` is sound.

## 4. PF.3 stability

The factor hypotheses give

```math
\|Y_iY_i^T-Z_iZ_i^T\|_{op}\le2\delta p_i.           \tag{APF.13}
```

Thus the two diagonal Boolean contributions total at most
`2delta(p_1^2+p_2^2)`.  The three-term expansion of the cross difference is
algebraically exact and has bilinear response at most
`(2delta+zeta)p_1p_2`; it is counted twice.  Dividing by
`p^2=(p_1+p_2)^2` gives

```math
2\delta+2\zeta{p_1p_2\over p^2}
\le2\delta+{\zeta\over2}.                            \tag{APF.14}
```

PF.25 and, after the sector factor two, PF.26 are correct.  The canonical
verifier tests 200 random perturbations satisfying the inward factor and
contraction bounds.

## 5. Net cardinality

With `tau=(eta/4)^2` and `delta=zeta=eta/16`,

```math
(2\sqrt\tau+\tau)+(4\delta+\zeta)
={13\eta+\eta^2\over16}
\le {7\eta\over8}.                                   \tag{APF.15}
```

The retained marginal ranks are at most `16/eta^2`.  An `r_1 by r_2`
operator contraction has Frobenius norm at most `sqrt(min(r_1,r_2))`.  A
maximal internal Frobenius `zeta`-net therefore has at most

```math
\left(1+{2\sqrt{\min(r_1,r_2)}\over\zeta}\right)^{r_1r_2} \tag{APF.16}
```

points.  Internal centres remain contractions, and Frobenius accuracy
implies the required operator accuracy.  Substitution gives base
`1+128/eta^2`, exponent at most `256/eta^4` per sector, and exactly the
two-sector logarithmic bound in PF.32.

This is a two-piece compatibility bound.  For many pieces the globally PSD
correlation carrier may require a growing joint net; independently netting
all pairwise contractions would ignore PF.9 and is not justified.  The draft
does not make that extension.

## 6. Response scaling and critical scope

PF.29 is a `d_q` error bound.  GE.20 converts it to normalized spherical
response error

```math
c\sqrt{\eta/2}+{c^2\eta\over8},
\qquad c={mp\over r}.                                 \tag{APF.17}
```

This is correct.  For `c=O(1)`, response accuracy `epsilon` requires
`eta=Theta(epsilon^2)` near the hard edge.  Consequently the reusable
marginal rank from this architecture is `O(epsilon^(-4))`, and the
compatibility-coordinate upper bound is `O(epsilon^(-8))`; both remain
independent of the total number of ports.

The compressible scaling is

```math
mp=O(r),\qquad m=O(r/p).                              \tag{APF.18}
```

At the original anti-pin choice `m=r`, one has `c=p`.  Constant spherical
response error then needs roughly `eta=O(p^(-2))`.  Feeding that accuracy
into PF.30--PF.32 produces rank bounds of order `p^4` (capped only by the
actual marginal sizes) and compatibility exponents of order `p^8`, before
logarithmic factors.  Thus PF.32 is not a growing-port compression theorem
for the anti-pin regime.  Section 6 states the positive `c=O(1)` assumption,
but this contrasting limitation should be made as explicit as it is in the
metric-entropy draft.

Under a uniform trust margin, GE.24 makes the response dependence linear in
metric error, but it still contains the `c^2` amplification.  The same
scaling distinction remains.

Nothing here controls the exact Boolean old-spin response or proves that
pairwise contraction states can be iterated without the global carrier.
The draft states both limitations honestly.

## 7. Verifier assessment

The canonical verifier passes 450 checks.  It verifies forward Douglas
constructions with rank-deficient marginals, random truncation, the
`0<tau<1` sharp examples used numerically, factor stability, and one forward
multi-piece correlation construction.  It does not test rank-zero
marginals, multi-piece converse recovery, PF.31--PF.34 arithmetic, or the
literal `tau=1` endpoint.

The independent verifier fills those gaps and reports:

```text
PSD compatibility-fibre independent audit: PASS
multi_recovery=2 zero_rank=2 tau_endpoint_diagnostic=1 arithmetic=3
```

Run both:

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_psd_compatibility_fibre.py

./.venv/bin/python \
  extremal_information/experiments/verify_psd_compatibility_fibre_independent_audit.py
```
