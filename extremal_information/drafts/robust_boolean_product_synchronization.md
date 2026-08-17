# Robust Boolean product synchronization

**Status.** Rigorous theorem.  Approximate spectral closure of the Fourier
products of one majority selector gives a uniform Boolean trust-response
bound.  The defect is stable under corresponding tensor products with an
additive, dimension-free law.  Combining this with the projective-row
response coreset gives a strict approximate carrier at fixed total port mass.

The result is deliberately narrower than arbitrary Boolean bridge
optimization.  It assumes one public selector and controls the spectral
defect of all products in that selector's Fourier support.  When the port
arity grows, certifying the resulting operator-norm condition may itself be
expensive unless the active products have a generated algebraic
presentation.

## 1. Normalization and the selector defect

Let `H` be a real symmetric `n times n` matrix and fix `r>0` such that

```math
\|H\|_{op}\le r.                                  \tag{RS.1}
```

Let `w_1,...,w_p in {+-1}^n` be Boolean ports.  Fix an antipodally odd
tie-broken majority selector

```math
\tau:\{+-1\}^p\longrightarrow\{+-1\},
\qquad \tau(-a)=-\tau(a),                         \tag{RS.2}
```

which equals `sgn(sum_i a_i)` off the zero-sum layer.  On that layer any
antipodally consistent tie rule is allowed.  Use the normalized Boolean
Fourier convention

```math
\widehat\tau(S)=2^{-p}\sum_{a\in\{+-1\}^p}
 \tau(a)\prod_{i\in S}a_i.                       \tag{RS.3}
```

Let `A_tau={S:hat tau(S)!=0}`.  Antipodal oddness makes every active set
odd.  For `S in A_tau`, put

```math
z_S=\bigodot_{i\in S}w_i,
\qquad Z=(z_S)_{S\in\mathcal A_\tau}.             \tag{RS.4}
```

Thus `Z` is `n times q`, where `q=|A_tau|`; repeated product columns, if
there are any, remain separately labelled.  Define

```math
G={Z^TZ\over n},\qquad
R={Z^THZ\over rn},\qquad
D=G-R.                                            \tag{RS.5}
```

The matrix `D` is automatically positive semidefinite:

```math
D={1\over n}Z^T(I-H/r)Z\succeq0.                 \tag{RS.6}
```

For an endpoint label `epsilon`, define `a^epsilon` by (RS.10) below and put

```math
\Delta_\tau(H,W)
=\max_{\epsilon\in\{+-1\}^p}(a^\epsilon)^TD a^\epsilon,
\qquad \delta=\|D\|_{op}.                         \tag{RS.6a}
```

Parseval will give `||a^epsilon||_2=1`, so `0<=Delta_tau<=delta`.
`Delta_tau` is the intrinsic selector defect; `delta` is a convenient
checkable sufficient certificate.  Both are **joint** defects.  They are not
the sum of the Rayleigh deficits of the active products and allow cross terms
to cancel before the quadratic form or norm is taken.

For `epsilon in {+-1}^p` and `m>=0`, set

```math
u_\epsilon=\sum_{i=1}^p\epsilon_iw_i             \tag{RS.7}
```

and define the labelled Boolean trust response

```math
\mathcal B_\epsilon(H;W,m)
=\max_{x\in\{+-1\}^n,\ \sigma\in\{+-1\}}
 \left\{{\sigma\over2}x^THx+m u_\epsilon^Tx\right\}.
                                                               \tag{RS.8}
```

## 2. Uniform robust selector theorem

### Theorem RS.1 (joint product-synchronization bound)

Under (RS.1)--(RS.8), every labelled channel obeys

```math
0\le {rn\over2}+m\|u_\epsilon\|_1
       -\mathcal B_\epsilon(H;W,m)
\le {rn\over2}(a^\epsilon)^TD a^\epsilon
\le {\Delta_\tau rn\over2}
\le {\delta rn\over2}.                            \tag{RS.9}
```

In particular, `Delta_tau=o(1)` gives `o(rn)` error simultaneously in all
`2^p` endpoint channels; `delta=o(1)` is sufficient.  Exact product closure
is the case `Delta_tau=0` (and in particular follows from `D=0`).

#### Proof

For each endpoint label define the Fourier coefficient vector

```math
a^\epsilon_S=\widehat\tau(S)\prod_{i\in S}\epsilon_i,
\qquad S\in\mathcal A_\tau.                      \tag{RS.10}
```

Parseval and `tau^2=1` give the dimension-free identity

```math
\|a^\epsilon\|_2^2
=\sum_S\widehat\tau(S)^2=1.                      \tag{RS.11}
```

Use the **same selector witness in all channels**:

```math
x_\epsilon(j)
=\tau(\epsilon_1w_1(j),\ldots,\epsilon_pw_p(j))
=Za^\epsilon.                                     \tag{RS.12}
```

It is Boolean.  Where `u_epsilon(j)` is nonzero, the majority rule makes
`x_epsilon(j)=sgn(u_epsilon(j))`; at a tie the field is zero.  Hence

```math
u_\epsilon^Tx_\epsilon=\|u_\epsilon\|_1.         \tag{RS.13}
```

Moreover, using `||x_epsilon||_2^2=n`, (RS.5), and (RS.11),

```math
\begin{aligned}
{rn\over2}-{1\over2}x_\epsilon^THx_\epsilon
 &= {rn\over2}(a^\epsilon)^T(G-R)a^\epsilon\\
 &\le {\Delta_\tau rn\over2}
 \le {\delta rn\over2}.                          \tag{RS.14}
\end{aligned}
```

Evaluating (RS.8) at `(x_epsilon,+1)` proves its lower bound.  Conversely,
for every Boolean `x` and either `sigma`,

```math
{\sigma\over2}x^THx\le {rn\over2},
\qquad
u_\epsilon^Tx\le\|u_\epsilon\|_1,              \tag{RS.15}
```

which proves the upper bound and (RS.9). `square`

The use of Parseval is essential: the estimate is independent of the number
`q` of active products.  Bounding the product channels separately would
replace (RS.11) by an `l_1` charge and can lose an exponential factor in
`p`.

There is also a diagonal-only sufficient bound which makes that possible
loss explicit.  For each active product let

```math
d_S=1-{z_S^THz_S\over rn}=D_{SS}\ge0.             \tag{RS.15a}
```

Since `D` is positive semidefinite,
`|D_{ST}|<=sqrt(d_Sd_T)`.  Therefore

```math
\boxed{
(a^\epsilon)^TD a^\epsilon
\le\left(\sum_{S\in\mathcal A_\tau}
 |\widehat\tau(S)|\sqrt{d_S}\right)^2.}          \tag{RS.15b}
```

Combining (RS.14) with (RS.15b) gives the same response theorem with
`delta` replaced by the right-hand side of (RS.15b).  In particular, if
every active product has Rayleigh deficit at most `d`, the loss parameter
is at most

```math
d\left(\sum_S|\widehat\tau(S)|\right)^2\le dq.    \tag{RS.15c}
```

This corollary is convenient when individual product relations are known,
but it exposes the growing-arity danger sharply: Fourier Parseval controls
the joint operator defect dimension-freely, whereas diagonal information
alone pays the Fourier `l_1` norm.  We make no unproved claim here about a
uniformly bounded majority Fourier `l_1` norm.

## 3. Tensor stability

Take two systems `(H_t,W_t,n_t,r_t)`, `t=1,2`, with the same port labels and
selector `tau`.  Form

```math
H_{12}=H_1\otimes H_2,
\qquad w_i^{12}=w_i^1\otimes w_i^2.               \tag{RS.16}
```

Then `||H_12||_op<=r_1r_2`, and every active product column is

```math
z_S^{12}=z_S^1\otimes z_S^2.                      \tag{RS.17}
```

Consequently, with `circ` denoting the entrywise Schur product,

```math
G_{12}=G_1\circ G_2,
\qquad R_{12}=R_1\circ R_2,                       \tag{RS.18}
```

and therefore

```math
\boxed{
D_{12}=D_1\circ G_2+R_1\circ D_2.}               \tag{RS.19}
```

There is a sharper semantic law before taking operator norms.

### Theorem RS.2 (intrinsic selector defects are tensor-subadditive)

For corresponding tensor products,

```math
\boxed{
\Delta_\tau(H_1\otimes H_2,W_1\boxtimes W_2)
\le\Delta_\tau(H_1,W_1)+\Delta_\tau(H_2,W_2).}   \tag{RS.19a}
```

Hence the `L`-fold intrinsic defect is at most `sum_t Delta_tau,t`.

#### Proof

For a row pattern `b in {+-1}^p`, let its active character vector be

```math
\chi(b)_S=\prod_{i\in S}b_i.                     \tag{RS.19b}
```

The Gram matrix is the row average

```math
G={1\over n}\sum_b N_W(b)\,\chi(b)\chi(b)^T.     \tag{RS.19c}
```

Also

```math
a^\epsilon\circ\chi(b)=a^{\epsilon\circ b}.     \tag{RS.19d}
```

Using `R_t=G_t-D_t`, expand the tensor defect symmetrically:

```math
D_{12}=D_1\circ G_2+G_1\circ D_2-D_1\circ D_2.  \tag{RS.19e}
```

By (RS.19c)--(RS.19d), for every `epsilon`,

```math
(a^\epsilon)^T(D_1\circ G_2)a^\epsilon
=\mathbb E_{b\sim\mu_{W_2}}
 (a^{\epsilon\circ b})^TD_1a^{\epsilon\circ b}
\le\Delta_{\tau,1},                              \tag{RS.19f}
```

and the second positive term is at most `Delta_(tau,2)`.  Finally
`D_1 circ D_2` is positive semidefinite by the Schur product theorem, so its
subtracted quadratic form is nonnegative.  This proves (RS.19a). `square`

The required norm estimate is not a positivity shortcut: `R_1` can be
indefinite.  We isolate the exact Schur-multiplier fact.

### Lemma RS.3 (Gram and contraction kernels are Schur contractions)

Suppose a matrix `K=(K_ij)` has a Hilbert-space factorization

```math
K_{ij}=\langle a_i,b_j\rangle,
\qquad \sup_i\|a_i\|\le1,
\qquad \sup_j\|b_j\|\le1.                       \tag{RS.20}
```

Then for every square matrix `X`,

```math
\|K\circ X\|_{op}\le\|X\|_{op}.                 \tag{RS.21}
```

#### Proof

Define diagonal embeddings on the standard basis by

```math
V_a e_i=e_i\otimes a_i,
\qquad V_b e_j=e_j\otimes b_j.                    \tag{RS.22}
```

Their norms are at most one, and direct inspection of entries gives

```math
K\circ X=V_a^*(X\otimes I)V_b.                   \tag{RS.23}
```

Taking operator norms proves (RS.21). `square`

For `G_t`, take `a_S=b_S=z_S^t/sqrt(n_t)`.  For `R_t`, put
`T_t=H_t/r_t`, take `a_S=z_S^t/sqrt(n_t)` and
`b_S=T_tz_S^t/sqrt(n_t)`.  Define

```math
\kappa_t=\max_{S\in\mathcal A_\tau}
 {\|H_tz_S^t\|_2\over r_t\sqrt{n_t}}\le1.         \tag{RS.23a}
```

The first factor family has norm one and the second has norm at most
`kappa_t`.  Thus the Gram Schur multiplier is contractive and the Rayleigh
Schur multiplier has norm at most `kappa_t`, even though only `G_t` is
necessarily positive semidefinite.

### Theorem RS.4 (dimension-free additive tensor defect)

For the corresponding tensor product (RS.16),

```math
\boxed{
\delta_{12}=\|D_{12}\|_{op}
\le\min\{\delta_1+\kappa_1\delta_2,
          \kappa_2\delta_1+\delta_2\}
\le\delta_1+\delta_2.}                            \tag{RS.24}
```

More generally, for `L` corresponding factors,

```math
\delta_{1\cdots L}\le\sum_{t=1}^L\delta_t.       \tag{RS.25}
```

#### Proof

Apply Lemma RS.3 to the two terms in (RS.19).  The `G_2` multiplier has norm
at most one and the `R_1` multiplier has norm at most `kappa_1`, giving the
first refined bound.  The symmetric identity

```math
D_{12}=D_1\circ R_2+G_1\circ D_2
```

gives the second.  Dropping the `kappa` factors and iterating proves
(RS.25). `square`

This is dimension-free but not automatically summable.  Tensoring one
factor of fixed positive defect `L` times gives only `L delta`; a long-depth
theorem needs exact factors, summable local defects, or an additional
forgetting mechanism.  Nor does (RS.24) make the active product list cheap:
for a generic growing-arity majority selector it can contain exponentially
many columns.  The theorem becomes a strict state reduction only when those
products or their defect admit a generated sub-landscape presentation.

## 4. Approximate projective-histogram carrier

Let `mu_W` be the probability histogram of the rows of `W`, modulo global
sign, on

```math
G_p=\{+-1\}^p/\{s\sim-s\}.                       \tag{RS.26}
```

For probability measures on `G_p`, recall

```math
d_p(\mu,\nu)=\max_{\epsilon\in\{+-1\}^p}
\left|\mathbb E_\mu{|s\cdot\epsilon|\over p}
     -\mathbb E_\nu{|s\cdot\epsilon|\over p}\right|.
                                                               \tag{RS.27}
```

Put `c=mp/r`, the total port-mass ratio.

### Corollary RS.5 (robust fixed-scale sub-landscape decoder)

If the intrinsic selector defect is at most `Delta_tau` and
`d_p(mu_W,nu)<=eta`, then the decoder

```math
\widehat{\mathcal B}_\epsilon(\nu)
={rn\over2}+mn\,\mathbb E_\nu|s\cdot\epsilon|     \tag{RS.28}
```

satisfies, simultaneously for all endpoint labels,

```math
\boxed{
\left|\mathcal B_\epsilon-
       \widehat{\mathcal B}_\epsilon(\nu)\right|
\le rn\left({\Delta_\tau\over2}+c\eta\right)
\le rn\left({\delta\over2}+c\eta\right).}        \tag{RS.29}
```

The dimension-free empirical response theorem supplies such a `nu` on at
most

```math
k=\left\lceil{16\over\eta^2}\right\rceil         \tag{RS.30}
```

projective row types.  An external fixed-scale cover therefore indexes the
decoder with at most

```math
(p-1)\left\lceil{16\over\eta^2}\right\rceil
                                                               \tag{RS.31}
```

bits, apart from the public scalars.  For bounded `c`, `Delta_tau=o(1)`, and
`eta=o(1)`, this gives `o(rn)` error with a state strictly smaller than the
full `2^n` Boolean landscape whenever `p/eta^2=o(n)`.  The scalar
`Delta_tau` is a certified class parameter, not part of the response table;
if it is not otherwise available, `||D||_op` supplies the public bound.

#### Proof

Because a row of `W` has projective type `s`,

```math
\|u_\epsilon\|_1
=n\mathbb E_{\mu_W}|s\cdot\epsilon|.              \tag{RS.32}
```

Equation (RS.27) bounds the change in this support by `np eta`.  Multiplying
by `m` gives `rn c eta`; combine with Theorem RS.1.  The coreset and bit
count are RC.1. `square`

Under corresponding tensor products the exact histograms convolve on
`G_p`, while Theorem RS.2 controls the intrinsic selector defect.  If factor histograms
`mu_t` have approximants `nu_t` with errors `eta_t`, then convolution
nonexpansiveness gives

```math
d_p(\mu_1*\cdots*\mu_L,\nu_1*\cdots*\nu_L)
\le\sum_{t=1}^L\eta_t.                            \tag{RS.33}
```

For the tensor product, whose final mass ratio is `c`, the convolved decoder
therefore has error at most

```math
rn\left\{ {1\over2}\sum_t\delta_t
                  +c\sum_t\eta_t\right}.         \tag{RS.34}
```

The sharper `sum_t Delta_(tau,t)` can replace `sum_t delta_t` in (RS.34).

As with the spectral defect, independent fixed errors accumulate.  Equation
(RS.34) is an honest reusable law, not a claim that a fixed coreset can be
silently resparsified at every depth.

## 5. Exact-sign completion consequence

Assume now that deleting the diagonal of `H` gives a hollow exact signing
and that `tr(H)=0`; this includes the regular symmetric Hadamard setting.
Append `p` auxiliary shores, each of width `m`, and give every vertex in
shore `i` the old--new column `w_i`.  Before filling auxiliary--auxiliary
edges, its cap is exactly

```math
\max_\epsilon\mathcal B_\epsilon(H;W,m).          \tag{RS.35}
```

Indeed, after fixing the old spin, every repeated shore sign can be chosen
to align its field with either sign of the old quadratic term.

Fill the `pm` auxiliary vertices by any public hollow sign matrix `C`.  If
`Q(C)` is its Boolean quadratic cap, pointwise triangle inequalities give

```math
|Q(\operatorname {Parent}_C)-
  \max_\epsilon\mathcal B_\epsilon|
\le Q(C)\le {pm\choose2}.                         \tag{RS.36}
```

Thus the histogram decoder for the completed exact signing has uniform
error at most

```math
rn\left({\delta\over2}+c\eta\right)+Q(C).         \tag{RS.37}
```

The sharper intrinsic `Delta_tau` may again replace `delta` in (RS.37).

If `pm=O(r)` and `r=o(n)`, then `Q(C)=O(r^2)=o(rn)`.  In particular, at the
Hadamard scale `r=sqrt(n)`, arbitrary public completion costs only `O(n)`
against the leading `n^(3/2)` response.  The total parent order is
`n+O(r)=n+o(n)`, so the leading normalization is unchanged.

## 6. Sharp scope and falsifiers

1. The equal-`(G,R)` four-port collision does not contradict RS.1: at least
   one active majority product has a fixed spectral defect, so `delta` does
   not vanish.
2. The theorem controls one selector's Fourier span jointly.  It does not
   assert that all Boolean optimizers lie in that span or reconstruct the
   parent maximum.
3. For `p=Theta(log n)`, full odd-product closure can have `Theta(sqrt n)`
   generated products and is a strict carrier.  For generic `p=Theta(n)`,
   an explicit active matrix can be exponentially wide; RS.1 alone then
   supplies accuracy but not information compression.
4. Fixed positive local defects are not washed out by tensor product.
   Summability in (RS.25) is a real hypothesis, not a consequence of
   dimension-free Schur contraction.

## 7. Verification

The companion script checks exact antipodal-selector examples, saturation of
the `delta rn/2` bound, tensor identities, random Schur contractions, the
histogram normalization, and the completion Lipschitz estimate:

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_robust_boolean_product_synchronization.py
```
