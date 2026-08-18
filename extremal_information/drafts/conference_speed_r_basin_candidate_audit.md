# Speed-`r` conference basin audit: affine fibres and weak template biases do not work

**Status.**  Task-local theorem report.  This note searches for an explicit
conference bridge family of probability `exp(-O(r))` whose pressure has a
fixed linear improvement over the uniform-bridge value.  It does not edit
canonical sources.

No favorable basin was found.  Two broad mechanisms can be ruled out
rigorously:

1. every affine parity fibre of codimension `O(r)`--including prescribed
   row/column products and every sparse planted-coordinate family of the
   right entropy--has the same leading pressure as a uniform bridge;
2. every finite-block weak template bias with total squared bias `O(r)`--
   including rank-one and block biases whose exact type shells have
   probability `exp(-O(r))`--fails to lower the leading pressure.

The first obstruction is deterministic.  The second uses a moment-matched
Lindeberg comparison and the convex-even bridge symmetry.  Together they
show that a speed-`r` favorable basin, if it exists, must use a nonlinear
correlation not reducible to sparse repair, affine parity information, or
weak independent/template bias.  They do not rule out extracting an
additional exponentially small pressure-correlated subset from a type
shell; doing that is the original lower-deviation problem again.

## 1. Setup and the exact edit modulus

Fix a symmetric conference signing `A_r`, an orientation
`epsilon in {+-1}`, and

```math
S_{\epsilon,B}
=\begin{pmatrix}A_r&B\\B^T&\epsilon A_r\end{pmatrix},
\qquad
t={\beta\over\sqrt{2r}},
\tag{CB.1}
```

```math
f_{\epsilon,r}(B)
=\log\left[2^{-2r}\sum_{x,y}
 \cosh\{t(H_{A_r}(x)+\epsilon H_{A_r}(y)+x^TBy)\}\right].
\tag{CB.2}
```

Throughout

```math
0<\beta<{\sqrt2\over6}.
\tag{CB.3}
```

The audited conference theorem gives, separately for the two orientations
and a uniform sign bridge `U_r`,

```math
{f_{\epsilon,r}(U_r)\over r}\longrightarrow h_\beta
\quad\hbox{in probability},
\qquad
{\mathbb Ef_{\epsilon,r}(U_r)\over r}\longrightarrow h_\beta,
\tag{CB.4}
```

where

```math
h_\beta=2\psi(\beta/\sqrt2)+{\beta^2\over4}.
\tag{CB.5}
```

The same-temperature child target is

```math
\tau_\beta=2\psi(\beta)=h_\beta-\gamma(\beta),
\qquad \gamma(\beta)>0.
\tag{CB.6}
```

For arbitrary real bridges `B,B'`, the pointwise Boltzmann comparison gives

```math
\boxed{
|f_{\epsilon,r}(B)-f_{\epsilon,r}(B')|
\le t\|B-B'\|_1.}
\tag{CB.7}
```

In particular, for sign bridges,

```math
\boxed{
|f_{\epsilon,r}(B)-f_{\epsilon,r}(B')|
\le2t\,d_H(B,B').}
\tag{CB.8}
```

This identifies the deterministic edit threshold.  A linear pressure
change cannot be forced by `o(r^(3/2))` sign edits.  Fixing only `O(r)`
coordinates costs the desired `O(r)` entropy, but can change pressure by
only `O(sqrt(r))`.  Conversely, the planted construction in
`conference_sublevel_gradient_audit.md` operates at
`Theta(r^(3/2))` edits and pays that superlinear entropy scale.

## 2. A general Hamming-retraction barrier

The right formulation is not restricted to coordinate subcubes.

### Theorem CB.1 (small-repair images retain the conference pressure)

Let `F_r` be a family of sign bridges.  Suppose there is a map

```math
\pi_r:\{+-1\}^{r\times r}\longrightarrow F_r
\tag{CB.9}
```

such that

1. if `U_r` is uniform on the full bridge cube, then `pi_r(U_r)` is uniform
   on `F_r`;
2. for every bridge `B`,

   ```math
   d_H(B,\pi_r(B))\le s_r=o(r^{3/2}).
   \tag{CB.10}
   ```

Then, separately for both orientations,

```math
\boxed{
{f_{\epsilon,r}(B_r)\over r}\longrightarrow h_\beta
\quad\hbox{in probability for uniform }B_r\in F_r.}
\tag{CB.11}
```

In particular, for every fixed `eta>0`,

```math
{1\over|F_r|}
\#\{B\in F_r:f_{\epsilon,r}(B)\le(h_\beta-\eta)r\}
\longrightarrow0.
\tag{CB.12}
```

Hence `F_r` itself cannot be a fixed-lower-deviation basin and cannot reach
the target `tau_beta r+o(r)`.

**Proof.**  Couple a uniform full bridge `U_r` to `pi_r(U_r)`.  By (CB.8),

```math
{1\over r}|f(\pi_r(U_r))-f(U_r)|
\le{2ts_r\over r}=o(1).
\tag{CB.13}
```

Now use (CB.4) and the pushforward assumption. `square`

This theorem is an exact response statement, not an entropy heuristic.

### Corollary CB.2 (all codimension-`O(r)` affine parity fibres fail)

Encode bridge signs as bits.  Let `F_r` be any nonempty affine solution set
of `q_r` independent binary parity equations.  Then

```math
|F_r|=2^{r^2-q_r}.
\tag{CB.14}
```

Choose `q_r` pivot coordinates for the parity-check matrix.  Given an
arbitrary bridge, retain its `r^2-q_r` free coordinates and overwrite the
pivots with the unique values satisfying the syndrome.  Every output has
exactly `2^{q_r}` preimages, so this map sends the uniform cube to the
uniform fibre and changes at most `q_r` entries.

Consequently, if `q_r=o(r^(3/2))`, Theorem CB.1 applies.  In particular,
every affine fibre with `q_r=O(r)` has the exact desired cardinality

```math
|F_r|=2^{r^2}\exp{-O(r)\},
\tag{CB.15}
```

but retains pressure `h_beta r+o(r)` rather than the smaller target.

This includes:

1. fixing `O(r)` arbitrary bridge coordinates;
2. planting a block containing `O(r)` prescribed entries;
3. adding any `O(r)` parity checks, regardless of their support sizes;
4. prescribing all row and column products.  The last system has rank
   `2r-1` when its syndromes satisfy the one consistency relation, hence

   ```math
   |F_r|=2^{(r-1)^2}.
   \tag{CB.16}
   ```

Thus neither sparse planting, row/column parity, nor an affine algebraic
fibre can be the required speed-`r` conference basin.

## 3. Weak template biases: a leading-order no-gain theorem

Affine fibres do not cover a fluctuation-scale rank-one bias.  Such a bias
can have relative entropy `Theta(r)` while its conditional mean has a
leading `sqrt(r)` singular value.  Nevertheless it cannot lower the
conference pressure.

The key fact is the exact symmetry

```math
\boxed{f_{\epsilon,r}(-B)=f_{\epsilon,r}(B).}
\tag{CB.17}
```

Indeed, substitute `y -> -y`; its internal quadratic energy is unchanged
and the bridge term changes sign.  The function is also convex in every
real bridge entry.

### Theorem CB.3 (independent weak biases cannot lower pressure)

Let `q_r` be a product law on sign bridges with deterministic coordinate
means

```math
m_e=\mathbb E_{q_r}B_e,
\qquad
\sum_{e=1}^{r^2}m_e^2\le Cr,
\tag{CB.18}
```

where `C` is fixed.  Then, separately for both orientations,

```math
\boxed{
\mathbb E_{q_r}f_{\epsilon,r}(B)
\ge h_\beta r-o(r).}
\tag{CB.19}
```

Moreover, for every fixed `eta>0`, all large `r` satisfy

```math
\boxed{
q_r\{f_{\epsilon,r}(B)\le(h_\beta-\eta)r\}
\le\exp\left\{-{\eta^2\over4\beta^2}r\right\}.}
\tag{CB.20}
```

**Proof.**  Let `V_e` be independent unbiased signs, put

```math
a_e=\sqrt{1-m_e^2},
\qquad
Y_e=m_e+a_eV_e,
\tag{CB.21}
```

and write `M=(m_e)` and `W=(a_eV_e)`.  Coordinatewise, `B_e` and `Y_e`
have the same first two raw moments.  For fixed values of the other bridge
coordinates,

```math
{\partial^3 f\over\partial B_e^3}
=t^3\kappa_3(\sigma x_i y_j),
\tag{CB.22}
```

where the observable is sign-valued under the current Gibbs law.  Hence
the third derivative has absolute value at most `2t^3`.  The ordinary
one-coordinate Lindeberg replacement, summed over the `r^2` entries, gives

```math
|\mathbb E_{q_r}f(B)-\mathbb E f(M+W)|
\le C_1r^2t^3=O_\beta(\sqrt r).
\tag{CB.23}
```

For completeness, the variables in a replacement are bounded by `2`, so
the Taylor remainders are uniformly summable; no high-temperature or
operator-norm hypothesis enters this step.

With the already chosen mean profile `M` fixed, also freeze
`a_e=sqrt(1-m_e^2)` and write `W_M=(a_eV_e)`.  The translated convolution

```math
G_M(Z)=\mathbb E_V f(Z+W_M)
\tag{CB.24}
```

is convex as a function of `Z`.  Since `W_M` and `-W_M` have the same law
and `f` is even, `G_M` is
even.  Therefore

```math
G_M(M)\ge G_M(0)=\mathbb E_Vf(W_M).
\tag{CB.25}
```

Finally, (CB.7) and `1-sqrt(1-u^2)<=u^2` give the pointwise comparison

```math
|f(W)-f(V)|
\le t\sum_e(1-a_e)
\le t\sum_em_e^2
=O_{\beta,C}(\sqrt r).
\tag{CB.26}
```

Combining (CB.4) and (CB.23)--(CB.26) proves (CB.19).

Changing one sign coordinate changes `f` by at most `2t`.  McDiarmid's
inequality under the product law gives

```math
q_r\{f\le\mathbb E_{q_r}f-u\}
\le\exp\left\{-{u^2\over\beta^2r}\right\}.
\tag{CB.27}
```

Take `u=eta r/2` and use (CB.19). `square`

This is stronger than the generic entropy-transport estimate for this
class.  A product tilt with `D(q_r||U_r)=Theta(r)` is allowed by transport
to gain `Theta(r)` pressure; CB.3 shows that its actual leading gain is
nonpositive, for every bias geometry satisfying (CB.18).

## 4. Exact rank-one and block type shells

The preceding product theorem has an exact microcanonical consequence.
Fix a set `E_0` of

```math
N_r=\theta r^2+O(r),
\qquad 0<\theta\le1,
\tag{CB.28}
```

and any sign template `R_e` on it.  Fix `c>0`, let

```math
m_r={c\over\sqrt r},
\qquad
k_r=\left\lfloor{1+m_r\over2}N_r\right\rfloor,
\tag{CB.29}
```

and define the exact type shell

```math
\mathcal T_r(R,E_0,c)
=\left\{B:\#\{e\in E_0:B_e=R_e\}=k_r\right\}.
\tag{CB.30}
```

Outside `E_0` the entries are free.  Its cardinality is exactly

```math
|\mathcal T_r|=2^{r^2-N_r}{N_r\choose k_r},
\tag{CB.31}
```

and Stirling's formula gives

```math
\boxed{
\log|\mathcal T_r|
=r^2\log2-{\theta c^2\over2}r+O(\log r).}
\tag{CB.32}
```

Thus this is precisely a speed-`r` family.  It covers a full rank-one
correlation shell by taking `R=uv^T`, an arbitrary algebraic template by
taking general `R`, and a block bias by restricting `E_0`.

### Corollary CB.4 (exact weak-template shells are typically not favorable)

For every fixed `eta>0`,

```math
\boxed{
{1\over|\mathcal T_r|}
\#\{B\in\mathcal T_r:
 f_{\epsilon,r}(B)\le(h_\beta-\eta)r\}
\le e^{-c_{\beta,\eta}r}}
\tag{CB.33}
```

for all large `r`, after harmlessly changing the positive constant.

**Proof.**  Apply CB.3 to the product law with mean `m_rR_e` on `E_0`
and zero outside.  Its squared-mean sum is
`N_rm_r^2=theta c^2r+O(1)`.  Conditional on having exactly `k_r`
template agreements, this law is uniform on `T_r`.  The conditioning event
is a central binomial point event of probability `Theta(1/r)`.  Divide
(CB.20) by that polynomial probability and absorb the factor into the
exponent. `square`

The result is deliberately stronger than saying that one shell member is
bad: an overwhelming fraction remains on or above the uniform pressure
rate.  It also permits finitely many disjoint blocks and prescribed type in
each block.  The exact count then has leading entropy loss

```math
{r\over2}\sum_{\ell}\theta_\ell c_\ell^2+O(\log r),
\tag{CB.34}
```

and conditioning costs only a fixed power of `r`.

## 5. One-sided row/block biases also contain manifestly high outputs

A still coarser candidate is a halfspace such as

```math
u^TBv\ge c r^{3/2},
\tag{CB.35}
```

or the `r` row-majority constraints
`u_i sum_j B_ijv_j>=0`.  Such families have speed at most `r`, but they do
not force low pressure.  They contain the rank-one bridge `B=uv^T`, whose
cross energy at `(x,y)=(u,v)` is `r^2`.

More generally there is an exact factorization, obtained by pairing `y`
with `-y`,

```math
\overline Z_{2r}(S_{\epsilon,B},t)
=\mathbb E_{x,y}
 \left[
 \cosh\{t(H_A(x)+\epsilon H_A(y))\}
 \cosh(t x^TBy)
 \right].
\tag{CB.36}
```

Thus any family containing a bridge with a large prescribed cross witness
also contains a high-pressure member; signed cancellation with the child
energies cannot rescue that member after the global-spin pairing.  For the
rank-one bridge, for example,

```math
f_{\epsilon,r}(uv^T)
\ge {\beta\over\sqrt2}r^{3/2}-2r\log2-O(1).
\tag{CB.37}
```

The point is not that every halfspace member is high, but that one-sided
row or rank-one constraints by themselves do not define a favorable
basin.  Selecting only their low-pressure members adds exactly the
unresolved pressure-dependent criterion.

## 6. Adaptive row/column gauge selectors also fail on the regular sector

The row-majority family has more structure than the preceding observation
uses.  It is a cross-section for row switching.  This permits a much
stronger uniform argument.

Let

```math
\mathcal G_r
=\big(\{+-1\}^r\times\{+-1\}^r\big)/\{(1,1),(-1,-1)\}
\tag{CB.38}
```

act on bridges by

```math
(s,u)\mathbin\cdot B=D_sBD_u.
\tag{CB.39}
```

It has `2^(2r-1)` elements and acts freely.  A cross-section therefore has
exactly `2^(r^2-(2r-1))` elements.  Row-only switching similarly has
`2^r` elements and cross-sections of size `2^(r^2-r)`.  Choosing the sign
of every row by its majority against a fixed vector, with a deterministic
tie rule that chooses one of each antipodal pair, is one such row
cross-section.

### Theorem CB.5 (adaptive gauge selection cannot create a regular basin)

Let `g_r(B)` be an arbitrary bridge-dependent row/column switching and put
`\widehat B=g_r(B)\mathbin\cdot B`.  Then, for every fixed `eta>0`, there
are positive constants `c_0,c_1` such that

```math
\boxed{
\Pr\{f_{\epsilon,r}(\widehat B)\le(h_\beta-\eta)r\}
\le e^{-c_0r}+e^{-c_1r^2}}
\tag{CB.40}
```

for all large `r`, when `B` is uniform.  Consequently the uniform law on
any row, column, or row--column switching cross-section has vanishing
lower-deviation fraction.  Its exact speed-`r` constraint cannot by itself
reach the child target.

**Proof.**  Choose `delta>0` and `kappa<1/2` as in the regular-sector
conference theorem.  The rectangular norm tail gives

```math
\Pr\{\|B\|_{op}>(2+\delta)\sqrt r\}\le e^{-c_0r}.
\tag{CB.41}
```

Every row/column switching preserves `||B||_op`.  On the complementary
event, the triangle estimate RC.5 puts the parent formed with every
possible switched bridge in the same operator-regular class.

For fixed `(s,u)`, conjugating the parent by `diag(D_s,D_u)` gives

```math
f_{A_r,\epsilon A_r}(D_sBD_u)
=f_{D_sA_rD_s,\epsilon D_uA_rD_u}(B).
\tag{CB.42}
```

Equivalently, change variables to `D_sBD_u`, which is again a uniform
bridge.  The regular-sector `exp(-c_eta r^2)` lower-tail theorem therefore
holds with the same constants for every fixed gauge.  Union over at most
`2^(2r-1)` gauges costs only `exp(O(r))`, leaving `exp(-c_1r^2)`.  Add the
operator-irregular probability (CB.41). `square`

The theorem permits the gauge to depend arbitrarily on the entire bridge;
it is not restricted to a local or computable selector.  Thus row-majority,
column-majority, simultaneous row/column canonicalization, and every other
switching-orbit transversal are closed as regular-sector basin mechanisms.
As before, the speed-`r` remainder is exactly the operator-irregular sector,
not a proved favorable family.

## 7. What the no-go theorems do and do not establish

The tested speed-`r` mechanisms now separate as follows.

| proposed basin | exact entropy | rigorous pressure verdict |
|---|---:|---|
| fix `O(r)` entries / sparse planted block | `r^2 log2-O(r)` | deterministic `O(sqrt r)` edit; rate remains `h_beta` |
| row/column products | `(r-1)^2 log2` | affine retraction changing at most `2r-1` bits; rate remains `h_beta` |
| any `O(r)` affine parity/algebraic constraints | `(r^2-O(r))log2` | affine retraction; rate remains `h_beta` |
| rank-one weak type shell | `r^2log2-c^2r/2+O(log r)` | an `1-e^{-Omega(r)}` fraction is not a fixed lower deviation |
| finitely many block/template weak type shells | `r^2log2-Theta(r)` | same Lindeberg/convex-even obstruction |
| rank-one or row-majority halfspace | at least speed `r` | contains manifestly very high bridges; constraint alone is not favorable |
| any adaptive row/column switching cross-section | `r^2log2-(2r-1)log2` or larger | regular lower deviations retain speed `r^2`; only operator-irregular outputs remain |
| planted `r^(3/4)`-by-`r^(3/4)` block | entropy loss `Theta(r^(3/2))` | can change pressure by `Theta(r)` but is too rare |

The new positive information is the structural boundary:

```text
speed-r linear constraints
    + a small Hamming repair, or
speed-r independent weak bias
    + convex-even bridge pressure

cannot create a leading favorable response.
```

The theorems do **not** prove that the full favorable set has probability
`e^{-omega(r)}`.  In CB.33, an `e^{-Theta(r)}` exceptional subset could
still have cardinality `exp(r^2 log2-O(r))`; defining that subset by the
pressure itself would merely rename the desired basin.  Nor do the results
cover nonlinear rowwise magnitude constraints not reducible to a gauge and
whose natural repair changes `Theta(r^(3/2))` entries, a joint spin--bridge
condition, or a genuinely
operator-irregular algebraic construction with an independent pressure
upper certificate.

## 8. Research judgment

No explicit speed-`r` target-reaching conference basin emerges from row or
column products, sparse planted blocks, rank-one/block fluctuation bias, or
ordinary affine algebraic fibres.  These are not inconclusive numerical
failures: CB.1--CB.4 give scalable analytic obstructions with exact entropy
counts.

Together with the regular-sector and Hamming-collar theorems, the only live
possibility is now narrower:

```math
\boxed{
\begin{array}{c}
\text{a nonlinear, deeply operator-irregular family of size}
\quad2^{r^2}e^{-O(r)},\\[2mm]
\text{with a pressure upper certificate not obtained by selecting on }f
\text{ itself.}
\end{array}}
\tag{CB.43}
```

The most discriminating next candidate is a rowwise nonlinear **magnitude**
fibre, rather than a sign/majority selector: one constant-probability
constraint per row on `|sum_jB_ijv_j|` has the correct entropy, is
switching-invariant, and its nearest repair can lie at the critical
`r^(3/2)`-edit scale.  To count as progress it must yield either a direct
upper bound `f<=(h_beta-eta)r` for all its members, or a theorem showing
that its conditioned pressure still has rate `h_beta`; another mean or
operator-norm calculation alone is insufficient.

## 9. Archive comparison

1. The global conference transport inequality allows a product law with
   `Theta(r)` relative entropy to gain `Theta(r)` pressure.  CB.3 is
   strictly sharper on the independent-bias class: convex evenness removes
   the entire possible leading gain.
2. The biased selector in `bounded_cap_boundary_roof_selector.md` uses a
   weak rank-one mean to *increase* a query-owned boundary response and
   localize an optimizer.  It does not study the conference bridge pressure
   and does not imply CB.3's no-lowering theorem.
3. The operator-regular theorem gives a speed-`r^2` lower tail after an
   analytic norm conditioning.  CB.1 instead handles arbitrary affine
   parity fibres without any norm hypothesis, by an exact uniform
   retraction.  CB.3 likewise does not assume operator regularity.
4. The planted mesoscopic increment in
   `conference_sublevel_gradient_audit.md` shows that
   `Theta(r^(3/2))` edits can change pressure by `Theta(r)`.  CB.1 proves the
   complementary deterministic statement that fewer-order edits cannot.
5. None of the archived conference notes found in the repository states
   the general affine-fibre retraction theorem, the weak-template
   Lindeberg/convex-even barrier, their exact type-shell consequence, or the
   adaptive-gauge regular-sector union theorem.
