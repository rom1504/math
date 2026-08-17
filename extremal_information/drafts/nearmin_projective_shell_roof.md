# Projective shell covers give local-field response roofs

Date: 2026-08-17.

Status: **proved draft; finite verifier supplied separately**.

This note proves the compression half of the exact-projective alternative.
A small projective cover of the positive shell determines the whole absolute
trust response on the corresponding `l_1` query ball, up to a quadratic
function of the cover radius.  Each chart stores only a signed centre, its
baseline, and its `n` oriented local fields.  It does **not** store the
internal edge matrix.

The result is static and one-block.  Constructing the shell cover may itself
require hard optimization, and no update or cross-order congruence is proved.

## 1. Signed cuts and their projective metric

Let `n>=3` and

```math
E={n\choose2},\qquad
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
Q=Q(A)=\max_x|H_A(x)|.
```

Write `c(x)=(x_ix_j)_(i<j)` and let

```math
\mathcal Z_n=\{\sigma c(x):\sigma\in\{\mathord\pm1\},
                         x\in\{\mathord\pm1\}^n\}
```

be the augmented-cut family.  For `z,z' in mathcal Z_n`, define

```math
d_{\rm P}(z,z')
=\min\{d_E(z,z'),d_E(z,-z')\}
={E-|\langle z,z'\rangle|\over2}.                 \tag{PR.1}
```

For `G>=0`, the positive signed-cut shell is

```math
\mathcal S_G^+(A)
=\{z\in\mathcal Z_n:\langle a,z\rangle\ge Q-G\}. \tag{PR.2}
```

Fix an integer radius `0<=R<E-floor(n^2/4)` and put

```math
k_R=\max\{0\le d\le\lfloor n/2\rfloor:d(n-d)\le R\}. \tag{PR.3}
```

The strict upper bound on `R` excludes the nearly-antipodal branch of the
signed-edge metric.  It is harmless in the asymptotic regime used below,
where `R/E` tends to zero.

### Lemma PR.1 (a projective signed ball is two local spin charts)

Let `z=\sigma c(x)` and `z_0=\sigma_0c(u)`.  If

```math
d_{\rm P}(z,z_0)\le R<E-\lfloor n^2/4\rfloor,       \tag{PR.4}
```

then there are `eta,tau in {+-1}` and a set `S subseteq[n]` with
`|S|<=k_R` such that

```math
\sigma=\eta\sigma_0,
\qquad x=\tau u^S.                                  \tag{PR.5}
```

Here `u^S` is obtained from `u` by flipping the coordinates in `S`.
Thus one projective ball splits into the two signed charts
`eta=+1` and `eta=-1`; `tau` is only the irrelevant global spin gauge.

#### Proof

Choose `eta` so that

```math
d_E(z,\eta z_0)=d_{\rm P}(z,z_0).
```

Let `d=min(d_H(x,u),n-d_H(x,u))` and globally negate `x` if needed, so
`d<=floor(n/2)`.  The unsigned cut distance is exactly

```math
d_E(c(x),c(u))=d(n-d).                              \tag{PR.6}
```

If `sigma ne eta sigma_0`, the distance in (PR.4) would be
`E-d(n-d)>=E-floor(n^2/4)>R`, a contradiction.  Hence the signed
orientations agree, (PR.6) is at most `R`, and `d<=k_R`.  This is (PR.5).
`square`

## 2. The local-field atlas

Let

```math
\mathcal C=\{z^r=\sigma_r c(u^r):1\le r\le L\}
\subseteq\mathcal S_G^+(A)                         \tag{PR.7}
```

be any shell-centred `R`-cover: every atom of `S_G^+(A)` has projective
distance at most `R` from some member of `C`.  One important way to obtain
such a cover is to take an inclusion-maximal `R`-separated family, where
separated means `d_P(z^r,z^s)>R`; maximality then gives the cover property.
The theorem below needs only the cover property, not separation.

For each centre retain only

```math
h_r=\sigma_rH_A(u^r),
\qquad
\ell_{r,i}=\sigma_r u_i^r(Au^r)_i\quad(1\le i\le n), \tag{PR.8}
```

together with `(sigma_r,u^r)`.  The values `ell_(r,i)` are the oriented
local fields at the centre.  No coefficient `a_ij` is retained.

For a real query `g in R^n`, define the atlas roof

```math
\widehat{\mathcal B}_{\mathcal C,R}(g)
=\max_{\substack{1\le r\le L,\ \eta,\tau\in\{\mathord\pm1\}\\
                  S\subseteq[n],\ |S|\le k_R}}
\left\{
 \eta h_r-2\eta\sum_{i\in S}\ell_{r,i}
 +\tau g\mathbin\cdot u^r-2\tau\sum_{i\in S}g_i u_i^r
\right\}.                                          \tag{PR.9}
```

The target response is

```math
\mathcal B_A(g)
=\max_{\sigma,x}\{\sigma H_A(x)+g\mathbin\cdot x\}. \tag{PR.10}
```

### Theorem PR.2 (uniform projective-shell response roof)

For every `g` with `||g||_1<=G`, one has the exact sandwich

```math
\boxed{
\left|\widehat{\mathcal B}_{\mathcal C,R}(g)
             -\mathcal B_A(g)\right|
\le 2k_R(k_R-1).}                                  \tag{PR.11}
```

The shell width needed in (PR.7) is exactly `G`, not `2G`.

#### Proof

First, choose a signed ground state `(sigma_*,x_*)`.  Since globally
negating `x_*` leaves its quadratic energy unchanged,

```math
\mathcal B_A(g)\ge Q+|g\mathbin\cdot x_*|\ge Q.    \tag{PR.12}
```

If `(sigma,x)` is an optimizer in (PR.10), then

```math
\sigma H_A(x)
=\mathcal B_A(g)-g\mathbin\cdot x
\ge Q-\|g\|_1\ge Q-G.                              \tag{PR.13}
```

Thus its augmented cut lies in `S_G^+(A)`.  This proves the claimed shell
width and is the only place the query-radius hypothesis is used.

For a centre `r`, orient and switch its edge signs by putting

```math
d^r_{ij}=\sigma_ra_{ij}u_i^ru_j^r.
```

The exact flip identity is

```math
\sigma_rH_A((u^r)^S)
=h_r-2\sum_{i\in S}\ell_{r,i}
   +4\sum_{\{i,j\}\subseteq S}d^r_{ij}.            \tag{PR.14}
```

After multiplying the quadratic orientation by `eta`, the part omitted by
(PR.9) has absolute value at most

```math
4{ |S|\choose2}\le2k_R(k_R-1).                     \tag{PR.15}
```

Every tuple `(r,eta,tau,S)` in (PR.9) corresponds to the genuine competitor

```math
(\eta\sigma_r,\ \tau(u^r)^S)
```

in (PR.10).  Equations (PR.14)--(PR.15) therefore show that every atlas
value is at most `B_A(g)+2k_R(k_R-1)`.

Conversely, cover the optimizer from (PR.13) by a centre `r`.  Lemma PR.1
gives precisely an `eta,tau,S` with `|S|<=k_R` representing it.  Its atlas
value is at least the optimizer's exact value minus (PR.15).  This proves
the other half of (PR.11). `square`

### Corollary PR.3 (sorting evaluator and information cost)

For fixed `(r,eta,tau)`, put

```math
w_i=-2\eta\ell_{r,i}-2\tau g_i u_i^r.              \tag{PR.16}
```

The optimum over `S` in (PR.9) is obtained by adding the largest at most
`k_R` positive values among the `w_i`.  Hence the atlas roof is evaluated
by four sorts per centre, in `O(Ln log n)` arithmetic comparisons (or by
linear-time order selection if desired).

The centre costs `n+1` bits, `|h_r|<=E`, and
`|ell_(r,i)|<=n-1`.  Thus the complete exact atlas has description length

```math
O(Ln\log n)\quad\hbox{bits}.                       \tag{PR.17}
```

This excludes the query `g`, which is input rather than retained state.
The baseline is actually redundant because `2h_r=sum_i ell_(r,i)`, but
keeping it makes the response formula transparent.

## 3. Normalized radius and the strict compression regime

Set `R=floor(gamma E)` and suppose

```math
0\le\gamma<{E-\lfloor n^2/4\rfloor\over E}.        \tag{PR.18}
```

Since `k_R(n-k_R)<=gamma E` and `n-k_R>=n/2`,

```math
k_R\le {2\gamma E\over n}=\gamma(n-1).             \tag{PR.19}
```

Consequently

```math
2k_R(k_R-1)\le2\gamma^2(n-1)^2.                    \tag{PR.20}
```

For a sequence `gamma_n=o(n^(-1/4))`, the uniform response error in
(PR.11) is `o(n^(3/2))`.  If at the same time

```math
L_n=o(n/\log n),                                   \tag{PR.21}
```

then (PR.17) is `o(n^2)` bits.  This is a strict static response
compression of the declared `l_1` query ball as a data representation, but
its **target-scale response content requires a qualification**.  Uniformly
on that query ball one always has the zero-state estimate

```math
Q\le\mathcal B_A(g)\le Q+G.                        \tag{PR.21a}
```

Hence if `G=o(n^(3/2))`, target-scale approximation is already trivial
without an atlas.  The atlas is genuinely finer than this scalar roof when
`k_R^2=o(G)`, and it is target-scale nontrivial when `G` itself is
macroscopic while (PR.20) is subleading.  It remains an exact quantitative
statement outside those regimes, but should not be advertised as new
`n^(3/2)`-scale compression there.  Nor is it an efficient encoder: finding
the required shell cover is not claimed to be easier than evaluating the
shell.

## 4. Packing or compression, and the critical `n^(-1/4)` scale

For the packing-versus-covering alternative, now choose `C` specifically to
be an inclusion-maximal `R`-separated family.  The same family gives an exact
finite alternative.

* Its balls cover the whole shell, so Theorem PR.2 gives an
  `O(Ln log n)`-bit response atlas with error at most `2k_R(k_R-1)`.
* Its centres satisfy

  ```math
  { |\langle z^r,z^s\rangle|\over E}<1-2\gamma
  \qquad(r\ne s),                                  \tag{PR.22}
  ```

  (The strict inequality is exact despite the floor, because the metric is
  integer valued.)  Thus Theorem AO.2 applies to the `G`-near-top centres
  with projective-overlap gap `Gamma=2gamma`.
  It turns `L` centres into `L` low-cap physical contextual states whenever
  its quantitative gap remains positive.

Here `L` is the size of the chosen maximal net, not the minimum number of
ordinary radius-`R` balls needed to cover the shell.  These quantities must
not be conflated.  If `P_R` denotes the maximum cardinality of an
`R`-separated subset, then every inclusion-maximal net has `L<=P_R`.  Also,
an ordinary radius-`R/2` cover of size `N` implies `L<=N`, since two
`R`-separated points cannot lie in the same such ball.  A bound only on the
minimum radius-`R` cover does not bound `L`.  Conversely, any explicitly
given small radius-`R` cover may be used directly in Theorem PR.2, but it
does not itself provide the separated AO.2 branch.

This is a genuine packing-versus-covering theorem, but it does **not** close
the near-minimizer frontier at one asymptotic scale.  With fixed sparse-flip
intensity `alpha` and shore ratio `lambda`, AO.2 has leading separation

```math
\asymp\gamma\min\{\alpha,\lambda\}n^{3/2},          \tag{PR.23}
```

while its simultaneous concentration error is

```math
O(\sqrt\alpha\,n^{5/4}+n).                         \tag{PR.24}
```

Even after optimizing `0<alpha<=1`, that compiler resolves the packing
uniformly only above the scale

```math
\gamma\gg n^{-1/4}
```

(and also requires `G=o(gamma n^(3/2))`).  In contrast, the local-field
atlas has subleading error precisely in the strict regime

```math
\gamma=o(n^{-1/4}).                                 \tag{PR.25}
```

At `gamma asymp n^(-1/4)`, both the omitted internal-edge term and the AO.2
sampling fluctuation are of leading order.  Thus AO.2 supplies the exact
opposite **finite** branch and a physical branch above the critical scale,
but not a target-scale physical conclusion for a large packing below it.
Removing this critical gap would require either a sharper joint physical
compiler or more information than the centre local fields.

## 5. Scope and frontier judgment

1. PR.2 is universal; near-minimality enters only through an explicit small
   projective cover (for compression), or through a bound on the appropriate
   packing/maximal-net number (for the same-family dichotomy).  An ordinary
   radius-`R` covering number and an `R`-packing number are not interchangeable.
2. The atlas answers one common continuum of trust queries uniformly.  It
   has no transition law, no reusable physical bridge, and no all-order
   realization theorem.
3. Small atlas size is an existential structural statement.  The retained
   state is strictly smaller than the full edge matrix only under (PR.21),
   but the theorem does not prove that exact minimizers satisfy it.  When
   `G=o(n^(3/2))`, (PR.21a) already gives target-scale scalar compression;
   the atlas adds information only at the finer scale described above.
4. Large fixed-scale packing is already handled physically by AO.2.  Large
   packing only at radii below `n^(-1/4)` is the unresolved middle case.
5. The theorem explains why the exponent `1/4` is intrinsic to the present
   pair of tools: local quadratic Taylor error is `Theta((gamma n)^2)`,
   while exact-sign sparse-flip concentration is `Theta(n^(5/4))`.

The result therefore sharpens `L_projective` into a multiscale question:
determine whether exact-minimizer shells have a small projective cover below
`n^(-1/4)`, a resolvable packing above `n^(-1/4)`, or genuinely critical
entropy concentrated at that scale.
