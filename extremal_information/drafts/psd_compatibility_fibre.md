# PSD compatibility fibres for collective Gram gluing

**Status.** Rigorous task-local exact gluing theorem, approximate truncation
law, and finite verifier.  This concerns the positive Gram sectors
`K^+-=(G+-R)/2`.  It is an exact statement about the spherical carrier, not
about the old-spin Boolean integrality gap.

Marginal PSD states do not determine their joint state.  The missing object
is a contraction between their support spaces.  This is the continuous
Gram analogue of a relative orientation fibre.  Unlike a scalar orientation
bit, it has `r_1 r_2` real coordinates at marginal ranks `r_1,r_2`.

## 1. Exact two-piece compatibility

Let `K_1\succeq0` and `K_2\succeq0` have sizes `p_1,p_2`, and consider

```math
K=\begin{pmatrix}K_1&C\\C^T&K_2\end{pmatrix}.          \tag{PF.1}
```

Choose skinny full-column-rank factors

```math
K_i=Y_iY_i^T,
\qquad Y_i\in\mathbb R^{p_i\times r_i},
\qquad r_i=\operatorname{rank}K_i.                   \tag{PF.2}
```

### Theorem PF.1 (Douglas compatibility fibre)

The block matrix (PF.1) is positive semidefinite if and only if

```math
\boxed{C=Y_1WY_2^T\quad\text{for a contraction}
       \quad W\in\mathbb R^{r_1\times r_2},\ \|W\|_{op}\le1.} \tag{PF.3}
```

For the fixed full-column factors in (PF.2), `W` is unique.  Under factor
gauge changes

```math
Y_i\longmapsto Y_iO_i,
\qquad O_i\in O(r_i),                                \tag{PF.4}
```

the same cross block is represented by

```math
W\longmapsto O_1^TWO_2.                              \tag{PF.5}
```

For fixed marginal factor gauges, the exact compatibility fibre is the
operator-norm unit ball.  Intrinsically, it is the space of presented
triples `(Y_1,Y_2,W)` modulo the simultaneous action (PF.4)--(PF.5).  One
must not double-quotient `W` alone: changing its singular vectors while
holding the labelled marginal factors fixed generally changes `C`.

#### Proof

Positivity of (PF.1), applied to `(x,ty)` and optimized over `t`, gives

```math
|x^TCy|^2\le(x^TK_1x)(y^TK_2y).                      \tag{PF.6}
```

Consequently the bilinear form

```math
(Y_1^Tx,Y_2^Ty)\longmapsto x^TCy
```

is well-defined and contractive on the two support spaces.  Riesz
representation gives a contraction `W` satisfying (PF.3).  Full column rank
makes it unique.

Conversely,

```math
K=
\begin{pmatrix}Y_1&0\\0&Y_2\end{pmatrix}
\begin{pmatrix}I&W\\W^T&I\end{pmatrix}
\begin{pmatrix}Y_1^T&0\\0&Y_2^T\end{pmatrix},        \tag{PF.7}
```

and the middle matrix is PSD exactly when `||W||op<=1`.  Equations
(PF.4)--(PF.5) are immediate. `square`

For the Gram--Rayleigh state, apply PF.1 independently to `K^+` and `K^-`.
Thus exact two-piece gluing requires contractions `(W^+,W^-)`.  At fixed
marginal ranks the compatibility dimension is at most

```math
r_{1,+}r_{2,+}+r_{1,-}r_{2,-}.                      \tag{PF.8}
```

It is created by composition: neither marginal PSD matrix contains it.

## 2. Multi-piece form and associativity warning

For pieces `i=1,...,s`, exact joint PSD compatibility is equivalently a
block correlation operator

```math
\Omega=(\Omega_{ij})_{i,j\le s}\succeq0,
\qquad \Omega_{ii}=I_{r_i},                          \tag{PF.9}
```

with

```math
K_{ij}=Y_i\Omega_{ij}Y_j^T.                         \tag{PF.10}
```

The gauge action is block orthogonal conjugation.  This gives an exact
associative presented carrier: restriction takes principal block matrices,
and gluing is PSD completion of the common presentation.

Pairwise contractions cannot be chosen independently once `s>=3`.
Condition (PF.9), not merely `||Omega_ij||<=1` for every pair, is the joint
compatibility law.  Storing independent pairwise `W_ij` without the global
PSD constraint would therefore be a false gluing theorem.

## 3. Truncating marginal spectra

Let `p=p_1+p_2`.  Assume the natural port normalization

```math
\operatorname{tr}K_i\le p_i.                       \tag{PF.11}
```

Fix `0<tau<=1`.  Let `P_i` be the spectral projection of `K_i` onto
eigenvalues strictly greater than `tau p_i`, and put `Q_i=I-P_i`.  Define

```math
K_i^h=P_iK_iP_i,
\qquad C^h=P_1CP_2,
\qquad
K^h=\begin{pmatrix}K_1^h&C^h\\(C^h)^T&K_2^h\end{pmatrix}. \tag{PF.12}
```

Then `K^h\succeq0` and

```math
\operatorname{rank}K_i^h\le {1\over\tau}.           \tag{PF.13}
```

The high cross block has the same form (PF.3), with `W` compressed to the
two retained support spaces.

For a symmetric matrix on `p` port coordinates, write

```math
q_p(M,M')={1\over p^2}
 \max_{\epsilon\in\{+-1\}^p}
 |\epsilon^T(M-M')\epsilon|.                         \tag{PF.14}
```

### Theorem PF.2 (cross-safe marginal truncation)

Every PSD joined state obeys

```math
\boxed{q_p(K,K^h)\le\sqrt\tau+{\tau\over2}.}         \tag{PF.15}
```

For a pair of Gram sectors, the query metric of GE.3 therefore satisfies

```math
d_q((G,R),(G^h,R^h))\le2\sqrt\tau+\tau.              \tag{PF.16}
```

#### Proof

Use PF.1 in positive-square-root coordinates,
`C=K_1^(1/2)WK_2^(1/2)`.  For Boolean vectors
`x in {+-1}^{p_1}`, `y in {+-1}^{p_2}`,

```math
x^T(K_i-K_i^h)x\le\tau p_i^2.                       \tag{PF.17}
```

Moreover

```math
C-P_1CP_2=Q_1C+P_1CQ_2.                             \tag{PF.18}
```

Contractivity of `W` and the spectral cutoff give, for either term,

```math
|x^TQ_1Cy|\le\sqrt\tau\,p_1p_2,
\qquad
|x^TP_1CQ_2y|\le\sqrt\tau\,p_1p_2.                 \tag{PF.19}
```

Remembering that the cross block occurs twice in the joined quadratic form,

```math
q_p(K,K^h)
\le {\tau(p_1^2+p_2^2)+4\sqrt\tau\,p_1p_2
       \over(p_1+p_2)^2}.                            \tag{PF.20}
```

Writing `a=p_1/p`, `1-a=p_2/p`, the numerator ratio is

```math
\tau+(4\sqrt\tau-2\tau)a(1-a),
```

whose maximum for `tau<=1` is `sqrt(tau)+tau/2` at `a=1/2`.  This proves
(PF.15).  Apply it separately to the two sectors and use
`d_q=2max(q_p^+,q_p^-)` to obtain (PF.16). `square`

The square root is unavoidable for marginal truncation.  For `0<tau<1`, take
`p_1=p_2=s`, unit all-ones directions `u,v`, and

```math
K_1=\tau s\,uu^T,
\qquad K_2=s\,vv^T,
\qquad C=\sqrt\tau\,s\,uv^T.                        \tag{PF.21}
```

This is PSD with `W=1`.  If the first eigenvalue is discarded at the strict
threshold `tau p_1`, the lost all-positive quadratic response is

```math
{\tau+2\sqrt\tau\over4}=\Theta(\sqrt\tau).           \tag{PF.22}
```

Thus retaining only eigenvalues above `eta p_i` can leave fixed
`Theta(sqrt(eta))` joined error.  To guarantee `O(eta)` after arbitrary PSD
gluing, the marginal cutoff must be of order `eta^2 p_i`, and the retained
rank is in general `O(eta^(-2))`, not `O(eta^(-1))`.

## 4. Netted factors and contractions

Write the retained blocks as

```math
K_i^h=Y_iY_i^T,
\qquad C^h=Y_1WY_2^T,
\qquad \|Y_i\|_F\le\sqrt{p_i},
\qquad \|W\|_{op}\le1.                              \tag{PF.23}
```

Suppose presented approximants satisfy

```math
\|Y_i-Z_i\|_F\le\delta\sqrt{p_i},
\quad \|Z_i\|_F\le\sqrt{p_i},
\quad \|W-V\|_{op}\le\zeta,
\quad \|V\|_{op}\le1.                              \tag{PF.24}
```

Let `\widetilde K^h` be reconstructed from `(Z_1,Z_2,V)`.

### Proposition PF.3 (presented-fibre stability)

```math
q_p(K^h,\widetilde K^h)\le2\delta+{\zeta\over2}.     \tag{PF.25}
```

Consequently, for two sectors,

```math
d_q((G^h,R^h),(\widetilde G^h,\widetilde R^h))
\le4\delta+\zeta.                                    \tag{PF.26}
```

#### Proof

The diagonal factor errors have operator norm at most `2delta p_i`, so
their total Boolean quadratic contribution is at most
`2delta(p_1^2+p_2^2)`.  Expanding the cross difference into

```math
(Y_1-Z_1)WY_2^T+Z_1(W-V)Y_2^T+Z_1V(Y_2-Z_2)^T
```

bounds its bilinear response by `(2delta+zeta)p_1p_2`.  The joined
quadratic counts this twice.  Division by `p^2` gives

```math
2\delta{p_1^2+p_2^2\over p^2}
+2(2\delta+\zeta){p_1p_2\over p^2}
=2\delta+2\zeta{p_1p_2\over p^2}
\le2\delta+\zeta/2.                                 \tag{PF.27}
```

This proves (PF.25), then (PF.26). `square`

## 5. Fixed-accuracy compatibility complexity

Choose, for `0<eta<=1`,

```math
\tau=(\eta/4)^2,
\qquad \delta=\zeta=\eta/16.                        \tag{PF.28}
```

Then PF.2--PF.3 give total joined error

```math
d_q\le(2\sqrt\tau+\tau)+(4\delta+\zeta)
\le {7\eta\over8}<\eta.                             \tag{PF.29}
```

Each marginal factor rank is at most

```math
r_i\le {16\over\eta^2}.                             \tag{PF.30}
```

Once the marginal factor states are given, each sector's compatibility
contraction therefore has at most `256/eta^4` real coordinates.  A maximal
internal Frobenius `zeta`-net of the contraction ball has cardinality at
most

```math
\left(1+{2\sqrt{\min(r_1,r_2)}\over\zeta}\right)^{r_1r_2}. \tag{PF.31}
```

For the two sectors, (PF.28)--(PF.31) give the explicit compatibility bound

```math
\log N_{\rm compat}
\le {512\over\eta^4}
 \log\left(1+{128\over\eta^2}\right).               \tag{PF.32}
```

This is `O_eta(1)`, independent of `p_1,p_2`.  The marginal factor nets
themselves still cost `O_eta(p_1+p_2)` logarithmic states, as required by the
linear Gram metric-entropy theorem.  Compatibility adds only a fixed number
of accuracy-dependent parameters.

The `eta^(-4)` displayed dependence is a robust upper bound, not an
optimality claim.  PF.21 proves only that the `eta^2` marginal spectral
threshold is unavoidable for this truncation architecture.

## 6. Consequence for collective spherical response

Combine PF.29 with the hard-edge response modulus GE.20.  If the total port
mass obeys

```math
c={mp\over r}=O(1),                                  \tag{PF.33}
```

then an `eta`-accurate joined compatibility carrier changes the normalized
SA.3 spherical response by at most

```math
c\sqrt{\eta/2}+{c^2\eta\over8}.                      \tag{PF.34}
```

Under a trust margin the loss is linear in `eta`.  This transfers the exact
PSD gluing law to a genuine approximate collective response law.

It still does not control the exact Boolean old-spin response.  Nor does it
say that pairwise compatibility matrices can be reused indefinitely without
retaining the joint PSD carrier (PF.9).

## 7. Theory interpretation

1. **Marginal state:** low-rank factors `Y_i^+-,` defined modulo orthogonal
   gauge.
2. **Composition-created state:** contractions `W^+-` in the chosen
   marginal frames; changing a frame acts simultaneously on its factor and
   every incident contraction.
3. **Exact gluing:** Douglas factorization PF.1, or the block correlation
   carrier PF.9 for several pieces.
4. **Approximate gluing:** marginal rank `O_eta(1)` plus `O_eta(1)`
   compatibility parameters, with explicit response error.
5. **New obstruction:** marginal spectral truncation has a square-root
   cross penalty.  A rank `O(eta^(-1))` marginal carrier, although sufficient
   for one closed query metric, is not uniformly reusable under arbitrary
   PSD composition at error `eta`.

This is a strict, quantitative instance of the general law

```math
\text{reusable state}
=\text{marginal response image}
+\text{compatibility fibre}.                         \tag{PF.35}
```

The response conclusion is compact only when `c=mp/r=O(1)`, equivalently
`m=O(r/p)` for growing `p`.  At the one-port anti-pin scaling `m=r`, one has
`c=p`; constant response accuracy then forces `d_q` error of order
`p^(-2)`, so neither the fixed-rank truncation nor its net remains
dimension-free.

## 8. Verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_psd_compatibility_fibre.py
```

The verifier checks Douglas reconstruction including singular marginals,
multi-piece block correlations, the truncation inequality on random and
extremal instances, the square-root example, and the factor/contraction
stability constant.
