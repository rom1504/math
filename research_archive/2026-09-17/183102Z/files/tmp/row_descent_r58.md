# Wave 58: relative bare-fibre transfer and the failure of migration-only row descent

## Status and relation to the ledger

The orientation-unified energy and row formulas below are **Verified**, but
their row part is not a new mechanism: in the same-orientation cut sector it
is exactly the shore identity (10.914).  The averaged orbit consequence
(10.1158)--(10.1159) is also already in the ledger.  The new useful object is
the exact transfer of the **bare favourable fibre** between two arbitrary
oriented cuts, expressed through a relative cut/uncut masked row.  Its exact
fixed-slice variance was checked by enumeration, and its exponential form is
the Bernoulli-conditioning/Hanson--Wright argument already verified in
(10.1124), now applied to the relative mask.

Current localized migration does **not** force strict descent of
`h+lambda R`, even after randomized or dyadic selection of the positive
field block.  A stored exact order-eight minimizer gives a finite bare-fibre
wall in the correct energy sector and even supplies exact-ground migrated
responses.  This is a finite universal-implication wall, not an asymptotic
counterexample to a project-scale row theorem.  It is closely related to,
but does not merely restate, the empty-neighbour walls immediately after
(10.1159): those concern hard complement columns and one-flips, whereas the
example here concerns the actual bare incidence and a nonlocal exact-ground
response.

Reproducible audit: `tmp/row_descent_r58_check.py`.

## 1. Exact relative cut/uncut calculus

Let

```math
d=(\tau,x),\qquad \omega=(\upsilon,z),\qquad
\kappa=\tau\upsilon,\qquad y_i=x_i z_i.
```

Gauge at `d` by

```math
s_{ij}=\tau a_{ij}x_ix_j,qquad
r_i=\sum_{j\ne i}s_{ij},qquad
R(d)=\sum_i r_i^2.
```

Define the relative disagreement set and its masked degrees by

```math
\mathcal D=\{ij:\kappa y_i y_j=-1\},\qquad
w_{ij}=s_{ij}{\bf1}_{\{ij\in\mathcal D\}},
\qquad b_i=\sum_{j\ne i}w_{ij},
\qquad W=\sum_{ij\in\mathcal D}s_{ij}.
\tag{R58.1}
```

If `kappa=+1`, `mathcal D` is the cut between the two level sets of `y`.
If `kappa=-1`, it is the complementary **uncut**, consisting of the edges
within the two level sets.  In both sectors

```math
s_{ij}(\omega)=s_{ij}(d)(1-2{\bf1}_{\{ij\in\mathcal D\}}).
```

Consequently the exact global identities are

```math
\boxed{
\begin{aligned}
\Delta(\omega)-\Delta(d)&=4W=2\sum_i b_i,\\
r_i(\omega)&=r_i(d)-2b_i,\\
R(\omega)-R(d)
&=4\lVert b\rVert_2^2-4\langle r(d),b\rangle\\
&=-4\sum_i b_i\{r_i(d)-b_i\}.
\end{aligned}}
\tag{R58.2}
```

The last line is (10.914) when the relative sector is a physical shore; the
only extension here is the uniform treatment of the opposite-orientation
uncut sector.  It exposes why counting reversed positive edges is
insufficient: row change depends on the complete signed degree vector `b`,
including every unselected positive edge and every negative edge of
`mathcal D`.

## 2. The exact bare-fibre transfer

For a selector `S`, put

```math
W_S=\sum_{ij\in\mathcal D,\ i,j\in S}s_{ij},
\qquad p_j=\frac{(m)_j}{(n)_j}.
```

The principal norm cancels when the two local deficits are subtracted, while
the parent-deficit term uses the first line of (R58.2).  Thus

```math
\boxed{
\widehat\ell(S,\omega)-\widehat\ell(S,d)
=4\{W_S-p_2W\}.}
\tag{R58.3}
```

This is the exact fibre-transfer variable missing from a migration argument.
It is centered under `S sim U_m`.  Equal, adjacent, and disjoint edge-pair
counting gives

```math
\boxed{
\begin{aligned}
\operatorname{Var}(W_S-p_2W)
={}&(p_2-2p_3+p_4)|\mathcal D|\\
&+(p_3-p_4)\lVert b\rVert_2^2
 +(p_4-p_2^2)W^2.
\end{aligned}}
\tag{R58.4}
```

Since `p_4-p_2^2<=0`, this yields

```math
\operatorname{Var}\{\widehat\ell(S,\omega)-\widehat\ell(S,d)\}
\le C_p\{n^2+\lVert b\rVert_2^2\},
\qquad
\lVert b\rVert_2
\le\frac{\sqrt{R(d)}+\sqrt{R(\omega)}}2.
\tag{R58.5}
```

The sharper statement uses the masked symmetric matrix `G=(w_ij)`.  With
`Y=diag(y)`,

```math
G=\frac12\{S_d-\kappa YS_dY\},
\qquad \lVert G\rVert_{\rm op}\le\lVert A\rVert_{\rm op},
\qquad G\mathbf1=b.
\tag{R58.6}
```

Applying the already verified Bernoulli-conditioning and Hanson--Wright proof
of (10.1124) to `G` gives, uniformly on a compact density window,

```math
\boxed{
\Pr_{U_m}\{|\widehat\ell(S,\omega)-\widehat\ell(S,d)|\ge a\}
\le C\sqrt n\exp\left[-c\min\left\{
\frac{(a-Cq_n/n)_+^2}{n^2+\lVert b\rVert_2^2},
\frac{(a-Cq_n/n)_+}{\lVert A\rVert_{\rm op}}
\right\}\right].}
\tag{R58.7}
```

Here the independent-to-slice mean correction is only `O(q_n/n)`: indeed
`2W` is the total quadratic energy of `G`, and
`|2W|=|E(d)-E(omega)|/2<=q_n`.  Also exact minimality gives the stored
spectral estimate

```math
\lVert A\rVert_{\rm op}\le\sqrt{2q_n}=O(n^{3/4}).
```

For fibres `F_t^d={S:widehat ell(S,d)<=t}`, (R58.7) implies the concrete
one-sided transfer

```math
\boxed{
u_{t+a}(\omega)\ge u_t(d)-\varepsilon_a(d,\omega),}
\tag{R58.8}
```

where `epsilon_a` is the right-hand side of (R58.7).  Notice the widened
threshold.  At a fixed threshold, a fibre may sit entirely at its boundary;
neither centering nor two-sided concentration gives the sign needed to keep
it.

## 3. Exact project-scale cost: fibre stability needs the relative row

Put

```math
H=n^{3/4-c},\qquad T_n=n^{3/2-c},\qquad
R_*=n^{9/4-c}.
```

If

```math
\lVert b\rVert_2^2\le C_bR_*,
```

then (R58.7), with `a=K T_n`, becomes

```math
\boxed{
\varepsilon_{KT_n}(d,\omega)
\le C\sqrt n\exp\{-c\min(K^2/C_b,K)H\}.}
\tag{R58.9}
```

The two exponents are exact at power level:

```math
\frac{T_n^2}{R_*}=n^{3/4-c}=H,
\qquad
\frac{T_n}{n^{3/4}}=H.
```

Thus if `u_t(d)>=e^{-C_0H}`, choosing a sufficiently large **fixed** `K`
in terms of `C_0,C_b` gives

```math
u_{t+KT_n}(\omega)\ge\frac12u_t(d),
\qquad h_\omega(t+KT_n)\le h_d(t)+\log2.
\tag{R58.10}
```

Both endpoint row bounds `R(d),R(omega)<=C R_*` are a simple sufficient
condition for the relative-row premise, but the exact requirement is only
the masked relative row `||b||^2`.  This is a real conditional bridge: once
the response and input have controlled relative row, favourable mass is
stable at the affordable entropy speed after spending another fixed multiple
of `T_n`.

It is also exactly circular for the proposed high-row migration.  If only
the generic response bound `R(omega)=O(n^{5/2})` is known, (R58.5) permits
`||b||^2=O(n^{5/2})`, and the Gaussian exponent in (R58.7) drops to

```math
\frac{T_n^2}{n^{5/2}}=n^{1/2-2c}=o(H).
\tag{R58.11}
```

More generally, if `R(d)` is a factor `L` above `R_*`, the available first
exponent is only `O(H/L)` unless cancellation makes `b` much smaller.  A
dyadic decomposition names this loss but does not remove it.  The elementary
variance bound is even less useful for a rare fibre: at target rows it loses
`O(H^{-1})` **additive** mass, much larger than `e^{-C_0H}`.

## 4. Why randomized positive blocks do not force scalar descent

For a correct-sector column `d`, (10.1371) extracts a block `E` of
`d`-positive edges from positive signed fields and localized migration
produces some near-ground `omega` with a constant fraction of `E` in
`mathcal D`.  This says only

```math
|E\cap\mathcal D|\ge c|E|.
```

Neither a random choice within a field dyadic level nor randomization over
levels controls:

1. the negative `d`-edges of `mathcal D`, which enter every coordinate of
   `b` with the opposite sign;
2. the unselected positive edges of `mathcal D`;
3. `||b||^2`, hence the sign of the row formula (R58.2); or
4. the correlation of the old favourable fibre with the signed transfer
   variable in (R58.3).

The response is existential and can depend adversarially on the sampled
block.  Averaging the blocks does not average one fixed response unless a new
minimax theorem controls the joint pair `(D,b)`.  The common active-face law
(10.1340) still controls only block marginals and has the same omission.

There is also a direct finite wall at a scalar optimum.  Use the stored exact
minimizer

```math
A_8=\begin{pmatrix}
0&1&1&1&1&1&1&1\\
1&0&1&-1&1&1&-1&-1\\
1&1&0&1&-1&1&-1&-1\\
1&-1&1&0&-1&-1&-1&1\\
1&1&-1&-1&0&-1&1&-1\\
1&1&1&-1&-1&0&1&1\\
1&-1&-1&-1&1&1&0&1\\
1&-1&-1&1&-1&1&1&0
\end{pmatrix},\qquad Q(A_8)=20.
```

Take `m=6,t=0`, and

```math
d=(-1;(1,-1,-1,-1,1,1,-1,-1)).
```

Exact enumeration gives

```math
E(d)=16,\quad \Delta(d)=4,\quad
r(d)=(3,5,1,-1,-1,3,5,1),\quad R(d)=72,
\quad u_0(d)=18/28.
\tag{R58.12}
```

Thus `d` is in the correct sector.  For `lambda=1/1000`, it is a global
minimizer of `h_d+lambda R(d)` over all 256 oriented projective states; all
eight optimizers have `(u,R,Delta)=(18/28,72,4)`, and the next objective gap
is `log(9/8)=0.117783...`.

Choose the field-positive block

```math
E=\{12,13\}.
```

It has `L=|E|=2` and satisfies the finite high-row inequality
`R(d)>=2(n-1)L`.  The exact parent ground

```math
\omega=(-1;(1,-1,1,-1,1,-1,-1,1))
```

reverses one of the two edges, stronger than the required one-quarter
escape, and has

```math
(u_0(\omega),R(\omega),\Delta(\omega))=(11/28,64,0),
\qquad |F_0^d\cap F_0^\omega|=8.
```

It really lowers row, but its scalar change is

```math
\boxed{
[h_\omega+\lambda R(\omega)]-[h_d+\lambda R(d)]
=\log(18/11)-8/1000
=0.484476...>0.}
\tag{R58.13}
```

A second response

```math
\omega'=(+1;(1,-1,-1,-1,1,1,1,1))
```

reverses both block edges, remains on the same deficit-four face, and has
`(u,R)=(16/28,72)`: even its row does not decrease.  Since `d` is already a
global scalar minimizer, **no** choice of response, randomized block, or
dyadic mixture can yield strict scalar descent at this price.

The small block does not satisfy the asymptotic numerical rounding condition
`2 eta_E/lambda_0<t<=4|E|` used to prove (10.1323).  The example instead
grants the desired qualitative conclusion directly, with an exact-ground
response.  It therefore falsifies an algebraic implication from the current
migration conclusion to scalar descent; it does not falsify an asymptotic
theorem that adds a new project-scale aggregate or stability hypothesis.

## 5. Resulting judgment

The migration-only row-descent proposal is closed in its present form.
Equations (R58.2)--(R58.3) isolate two independent missing controls:

```math
\langle r(d),b\rangle-\lVert b\rVert^2>0
```

at enough scale to lower row, and favourable-fibre mass on the correct side
of the centered transfer.  Randomized selection of `E` controls neither.

The surviving positive result is (R58.7)--(R58.10).  A future theorem that
constructs a response with `R(omega)=O(R_*)` and either `R(d)=O(R_*)` or the
sharper `||b||^2=O(R_*)` automatically obtains entropy-speed fibre stability
after a fixed tolerance enlargement.  It still does not give descent of the
same-threshold scalar objective without a margin statement.  For genuinely
high rows `R(d)/R_* -> infinity`, a new relative-row susceptibility or a
multi-stage theorem is indispensable; a dyadic relabeling alone loses the
factor `R(d)/R_*` in the exponent.

This also clarifies the status of the older orbit route.  Formula (10.1158)
can lower row on average, and (10.1159) would finish if average log-fibre
inflation were controlled.  The exact relative transfer above identifies
that missing inflation, while the A8/A9 empty-neighbour examples and the new
bare A8 response show why it is not a consequence of row decrease itself.
