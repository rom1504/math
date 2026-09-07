# Wave 16: the common-mosaic response dual and its localization wall

## Status and conclusion

- **Verified:** the finite pure response problem is a minimum of support
  functions.  Its exact lower certificate is tuple-dependent.  Convexifying
  the local minimizers gives a genuine LP whose exact dual is a probability
  distribution on signed global states plus one lower-support variable per
  block.
- **Verified:** the dual has exact block-local *signed* deficit marginals, but
  neither their sign nor a temporal/endpoint location is supplied by duality.
  Nonnegative block totals can always be manufactured from a positive scalar
  dual value, but this is only the root/full-set inequality.
- **Verified finite counterexample:** on the order-nine `6+1+1+1` capture wall,
  the pure response gap is `12`, whereas the convexified/mixed LP gap is only
  `8`.  Its common dual distribution has one block marginal equal to the
  already-defined excess `X=12` and average global slack `4`.  Thus the local
  marginal literally restates `X` and does not even certify the pure gap.
- **Verified finite counterexample:** the endpoint-neutral `A_8` admits a
  `4+4` partition with pure response gap `8`, mixed gap `4`, and nonnegative
  dual block marginals `(4,4)`, while every **root** tolled residual is zero.
  Conservation leaves zero descendant allocation, so no universal map from
  response marginals to the existing allocation-weighted endpoint resource
  can exist, even with an arbitrary finite dilation.
- **Falsified:** an unconditional, zero-error inference

  ```math
  \text{forced common-mosaic response}
  \Longrightarrow
  \text{actual endpoint service or all causal Hall cuts}.
  ```

  The LP proves only a scalar separation.  A compatible nonnegative flow is
  an additional Hall theorem, and A5/A7/A8 violate respectively the relevant
  prefix, suffix, and endpoint-capacity implications.
- **Open:** a bridge with an independently proved temporal cocycle, localized
  terminal/negative credits, structural exclusion of neutral mosaics, window
  contraction, or a tail-summable additive error.  The finite examples do not
  rule out such an asymptotic theorem; they specify constraints it must meet.

All new finite assertions below are checked with exact integer enumeration by
`tmp/check_response_dual_r16.py`.  No claim relies on an LP solver tolerance.

## 1. The finite response profile

Let an order-`n` global minimizer be partitioned into nonempty blocks `V_i`:

```math
A=C+\bigoplus_{i=1}^kD_i,
\qquad Q(A)=q_n.
```

Here `C` is the fixed cross-only mosaic.  Write

```math
s_i=Q(D_i),\qquad b_i=q_{|V_i|},\qquad
x_i=s_i-b_i\ge0,
```

and `S=\sum_i s_i`, `B=\sum_i b_i`, `X=S-B=\sum_i x_i`.
Let `\mathcal M_i` be the finite set of all labelled order-`|V_i|`
minimizers and let

```math
\Gamma=\prod_i\mathcal M_i.
```

For `\mathbf G=(G_i)\in\Gamma`, put

```math
H_{\mathbf G}=C+\bigoplus_iG_i.
```

Use the finite signed projective state space

```math
\Omega=\{(\sigma,z):\sigma\in\{\pm1\},\ z\in\{\pm1\}^n,
\ z_1=1\}.
```

For `\omega=(\sigma,z)` define the global ground slack

```math
a_\omega=q_n-\sigma z^{\mathsf T}Az\ge0
```

and the block deficit difference

```math
d_{i,G}(\omega)
=\delta_{D_i}(\sigma,z_i)-\delta_G(\sigma,z_i),
\qquad
\delta_M(\sigma,y)=Q(M)-\sigma y^{\mathsf T}My.
```

The profile coordinate relevant to response is

```math
\begin{aligned}
f_{\mathbf G}(\omega)
&=\sum_i d_{i,G_i}(\omega)-a_\omega\\
&=X+\sigma z^{\mathsf T}H_{\mathbf G}z-q_n.
\end{aligned}
\tag{R16.1}
```

The second equality is exact: the `D_i` energies cancel against the `A`
part of `a_\omega`.  Since

```math
\Phi_{\mathbf D}=q_n-S,
\qquad
\Phi_{\mathbf G}=Q(H_{\mathbf G})-B,
```

(R16.1) gives

```math
\max_{\omega\in\Omega}f_{\mathbf G}(\omega)
=\Phi_{\mathbf G}-\Phi_{\mathbf D}
=X+Q(H_{\mathbf G})-q_n.
\tag{R16.2}
```

Thus the fixed-partition pure value is

```math
\boxed{
\rho_{\rm pure}
=\min_{\mathbf G\in\Gamma}\max_{\omega\in\Omega}
f_{\mathbf G}(\omega)
=X+\min_{\mathbf G\in\Gamma}[Q(H_{\mathbf G})-q_n].
}
\tag{R16.3}
```

This is exactly (10.575)--(10.576) for one partition.  It is a minimum of
finite support functions, not a linear program in the local blocks.

### 1.1 Exact pure separation certificate

Let `v_{\mathbf G}\in\mathbb R^\Omega` have coordinates
`f_{\mathbf G}(\omega)`.  For a proposed lower bound `r`,

```math
\rho_{\rm pure}\ge r
\quad\Longleftrightarrow\quad
\forall\mathbf G\in\Gamma\ \exists\omega\in\Omega:
f_{\mathbf G}(\omega)\ge r.
\tag{R16.4}
```

Equivalently, the finitely many tuples are covered by the state halfspaces

```math
\Gamma=\bigcup_{\omega\in\Omega}
\{\mathbf G:f_{\mathbf G}(\omega)\ge r\}.
```

An exact probability form is

```math
\boxed{
\rho_{\rm pure}
=\max_{(p^{\mathbf G})_{\mathbf G\in\Gamma}}
\min_{\mathbf G\in\Gamma}
\sum_\omega p^{\mathbf G}_\omega f_{\mathbf G}(\omega),
\qquad p^{\mathbf G}\in\Delta(\Omega).
}
\tag{R16.5}
```

For each tuple, take `p^{\mathbf G}` concentrated on one maximizing state.
Conversely, every expectation is bounded by the corresponding maximum.
Formula (R16.5) is an LP if one introduces a common lower variable and a
separate simplex `p^{\mathbf G}` for every tuple.  Its decisive feature is
that the separating state distribution may depend on the *entire* tuple.
There is no common state measure and no block factorization.

## 2. Convexified/mixed minimizers and the exact LP dual

Now let the minimizer choose a probability distribution on local minimizers
and let the state player respond only to the *expected* profile.  Correlation
between blocks is immaterial: only the block marginals enter the expectation,
and a product distribution realizes any list of marginals.  Write

```math
\lambda_{i,G}\ge0,
\qquad \sum_{G\in\mathcal M_i}\lambda_{i,G}=1,
\qquad K_i=\sum_G\lambda_{i,G}G.
```

The mixed value is

```math
\boxed{
\rho_{\rm mix}
=X+\min_{K_i\in\operatorname{conv}\mathcal M_i}
\left[Q\!\left(C+\bigoplus_iK_i\right)-q_n\right].
}
\tag{R16.6}
```

It has the primal LP

```math
\begin{array}{ll}
\text{minimize}&t\\[2mm]
\text{subject to}&
t\ge X-q_n+\sigma z^{\mathsf T}Cz
+\displaystyle\sum_i\sum_{G\in\mathcal M_i}
\lambda_{i,G}\,\sigma z_i^{\mathsf T}Gz_i
\quad(\omega=(\sigma,z)\in\Omega),\\[2mm]
&\sum_G\lambda_{i,G}=1,\quad \lambda_{i,G}\ge0.
\end{array}
\tag{R16.7}
```

Its exact dual is

```math
\boxed{
\begin{array}{ll}
\text{maximize}&
X-q_n+\displaystyle\sum_\omega p_\omega\sigma z^{\mathsf T}Cz
+\sum_i\eta_i\\[2mm]
\text{subject to}&p_\omega\ge0,\quad\sum_\omega p_\omega=1,\\[1mm]
&\displaystyle
\eta_i\le\sum_\omega p_\omega\sigma z_i^{\mathsf T}Gz_i
\quad(G\in\mathcal M_i).
\end{array}
}
\tag{R16.8}
```

This also follows directly from finite minimax:

```math
\rho_{\rm mix}
=\max_{p\in\Delta(\Omega)}
\left\{
X-q_n+\mathbb E_p[\sigma z^{\mathsf T}Cz]
+\sum_i\min_{G\in\mathcal M_i}
\mathbb E_p[\sigma z_i^{\mathsf T}Gz_i]
\right\}.
\tag{R16.9}
```

### 2.1 Separation/KKT certificate

Let `K_i^*` and `p^*` be primal/dual feasible.  They certify equality in
(R16.6) exactly when

1. `p^*` is supported on signed states active for
   `Q(C+\bigoplus_iK_i^*)`;
2. every `G` used with positive `\lambda_{i,G}` minimizes
   `\mathbb E_{p^*}[\sigma z_i^{\mathsf T}Gz_i]` over `\mathcal M_i`;
3. the two displayed objective values agree.

Geometrically, `p^*` is a probability subgradient of the coordinatewise
maximum at the optimal convexified profile and

```math
\langle p^*,v-v^*\rangle\ge0
```

for every convexified profile `v`.  This is the exact common separating
hyperplane sought by the support-function formulation.

### 2.2 The exact block deficit marginals

Put

```math
\boxed{
r_i(p)
=x_i+min_{G\in\mathcal M_i}
\mathbb E_p\!\left[
\sigma z_i^{\mathsf T}(G-D_i)z_i
\right].
}
\tag{R16.10}
```

Equivalently,

```math
r_i(p)
=\mathbb E_p\delta_{D_i}
-\max_{G\in\mathcal M_i}\mathbb E_p\delta_G.
```

Then (R16.9) becomes the exact deficit form

```math
\boxed{
\rho_{\rm mix}
=\max_{p\in\Delta(\Omega)}
\left[\sum_i r_i(p)-\mathbb E_p a_\omega\right].
}
\tag{R16.11}
```

An entirely linear version introduces unrestricted variables `r_i` and uses

```math
r_i\le x_i+\sum_\omega p_\omega
\sigma z_i^{\mathsf T}(G-D_i)z_i
\quad(G\in\mathcal M_i),
\tag{R16.12}
```

while maximizing `\sum_i r_i-\sum_\omega p_\omega a_\omega`.

These `r_i` are the only block-local marginals furnished by the exact dual.
They have three limitations.

1. **They are signed.**  For example, if `D_i` is itself a minimizer with
   `q_{|V_i|}>0` and `p` is concentrated on a positive ground of `D_i`, then
   `-D_i\in\mathcal M_i` and (R16.10) gives `r_i\le-2q_{|V_i|}`.  No sign is
   built into LP duality.
2. **Positive totals can be made nonnegative only by bookkeeping.**  If a
   dual certificate has value `\rho\ge0`, then

   ```math
   \sum_i r_i=\rho+\mathbb E_pa\ge\rho.
   ```

   Hence one may shave the positive parts to numbers
   `0\le s_i\le[r_i]_+` with `\sum_i s_i=\rho`.  This construction uses only
   the scalar total; it has no canonical state support or endpoint location.
3. **Decomposing `X` is already trivial.**  The numbers `x_i\ge0` sum to `X`
   by definition.  A dual for which `r_i=x_i` has not extracted new matrix
   service unless it also proves a compatibility/capacity statement.  The
   A9 and A8 certificates below have exactly this tautological form.

### 2.3 Why this mixed game is not a deterministic replacement theorem

If the adversary observes the realized tuple, randomization gives

```math
\mathbb E_{\mathbf G}\max_\omega f_{\mathbf G}(\omega),
```

whose minimum over distributions is attained by a pure tuple and equals
`\rho_{\rm pure}`.  The LP instead uses

```math
\max_\omega\mathbb E_{\mathbf G}f_{\mathbf G}(\omega)
```

and may be strictly smaller.  Therefore replacing the tuple-dependent
certificates in (R16.5) by one common `p` silently changes the pure problem
to its convexified relaxation.

## 3. Exact order-nine test: pure `12`, mixed `8`

Use `A_9` from (10.550) and the partition from (10.579):

```math
V=(0,1,2,4,5,6),\qquad\{3\},\{7\},\{8\}.
```

As in the ledger,

```math
Q(A_9)=24,quad Q(C)=22,quad Q(D_V)=22,quad q_6=10,quad X=12.
\tag{R16.13}
```

Complete enumeration of all 384 labelled order-six minimizers gives 40
hybrids of norm `24` and no smaller hybrid.  Therefore

```math
\boxed{\rho_{\rm pure}=12.}
\tag{R16.14}
```

### 3.1 Exact mixed primal certificate

The two order-six minimizers

```math
G^{(a)}=
\begin{pmatrix}
0&-1&1&1&-1&-1\\
-1&0&-1&1&-1&1\\
1&-1&0&1&1&1\\
1&1&1&0&-1&1\\
-1&-1&1&-1&0&1\\
-1&1&1&1&1&0
\end{pmatrix},
```

```math
G^{(b)}=
\begin{pmatrix}
0&-1&1&-1&1&-1\\
-1&0&-1&1&1&-1\\
1&-1&0&1&-1&-1\\
-1&1&1&0&-1&-1\\
1&1&-1&-1&0&-1\\
-1&-1&-1&-1&-1&0
\end{pmatrix}
```

satisfy

```math
Q\!\left(C\oplus_{\rm diag}\frac{G^{(a)}+G^{(b)}}2\right)=20.
\tag{R16.15}
```

The checker verifies (R16.15) by doubling the matrix and finding exact norm
`40`.

### 3.2 Exact mixed dual certificate

Give mass `1/2` to each of

```math
\begin{aligned}
\omega_+&=\left(+1,(1,1,-1,-1,-1,1,-1,-1,-1)\right),\\
\omega_-&=\left(-1,(1,1,-1,1,-1,1,-1,1,1)\right).
\end{aligned}
\tag{R16.16}
```

Their restrictions to `V`, in the displayed block order, are identical:

```math
(1,1,-1,-1,1,-1).
```

The orientations are opposite.  Consequently, for **every** symmetric local
matrix `G`, not only for a minimizer,

```math
\mathbb E_p[\sigma z_V^{\mathsf T}Gz_V]=0.
\tag{R16.17}
```

Their oriented cross energies are respectively `18` and `22`, so every
fractional hybrid has expected oriented energy `20`.  Hence its norm is at
least `20`, proving equality in (R16.15).  Therefore

```math
\boxed{\rho_{\rm mix}=X+20-q_9=8<12=\rho_{\rm pure}.}
\tag{R16.18}
```

The two oriented `A_9` energies are `16` and `24`.  Thus

```math
\mathbb E_pa=\frac{(24-16)+(24-24)}2=4.
```

By (R16.17), the unique nontrivial block marginal is

```math
r_V=x_V+0-0=12=X.
```

The dual identity is exactly

```math
\rho_{\rm mix}=12-4=8.
\tag{R16.19}
```

This is the requested stress test.  A common separating distribution does
produce a nonnegative block number, but that number is precisely the input
excess `X`; sign cancellation makes it independent of every replacement.
It is not a new endpoint or temporal resource, and the global slack prevents
it from certifying the pure response.

### 3.3 Persistent grounds explain the pure gap

For each of the 40 norm-24 pure completions, exact enumeration finds a signed
projective state that is a ground of both `A_9` and the hybrid.  At every such
persistent ground the signed profile is

```math
(C,D,G)=(14,10,10)\quad\text{or}\quad(18,6,6).
\tag{R16.20}
```

Thus the internal state energy does not change.  The old/new deficit pairs
are `(12,0)` or `(16,4)`, so `d_V=12` solely because the local baseline falls
from `Q(D)=22` to `q_6=10`.

More explicitly, for

```math
A_t=C+(1-t)D+tG
```

convexity gives `Q(A_t)\le24`, while the persistent ground gives the reverse
inequality.  Hence `Q(A_t)=24` for every `t\in[0,1]`, but the linearly based
response is

```math
\Phi_t=24-(22-12t)=2+12t.
\tag{R16.21}
```

This also explains (10.581): the response shift is a falling deficit
baseline, not work done by the global norm.  Subtracting it as scalar credit
would falsely give `24\le22+10-12=20`.

## 4. Exact A8 no-go for endpoint-capacity localization

Use the endpoint-neutral matrix `A_8` from (10.445) and partition it into

```math
V_1=(0,3,4,5),\qquad V_2=(1,2,6,7).
\tag{R16.22}
```

Exact enumeration gives

```math
Q(A_8)=20,quad Q(D_1)=Q(D_2)=12,quad
q_4=8,quad X=(4+4)=8,quad Q(C)=16.
\tag{R16.23}
```

Among all `48^2` labelled minimizing replacements, the least hybrid norm is
`20`, attained by 16 pairs.  Hence

```math
\boxed{\rho_{\rm pure}=8.}
\tag{R16.24}
```

For the mixed problem, `0\in\operatorname{conv}\mathcal M_4` because `G`
and `-G` are both minimizers, so the cross matrix itself is feasible and has
norm `16`.  The matching dual puts mass `1/2` on

```math
\omega_+=(+1,(1,1,1,1,1,1,1,1))
```

and

```math
\omega_-=(-1,(1,-1,-1,1,1,1,-1,-1)).
```

On each block the two spin restrictions agree up to an overall sign, while
the orientations are opposite.  Every local quadratic energy therefore
cancels.  Both oriented cross energies equal `16`, proving the mixed lower
bound.  The oriented `A_8` energies are `12` and `20`.  Consequently

```math
\boxed{
\rho_{\rm mix}=8+(16-20)=4,qquad
(r_1,r_2)=(4,4),qquad \mathbb E_pa=4.
}
\tag{R16.25}
```

Again, the nonnegative marginals are exactly the already-known local
excesses.

Now enumerate the four projective positive and four projective negative
grounds of `A_8`.  Every one of their 16 pairs has two four-vertex shores.
On every directed shore,

```math
Q(X)=8,qquad h_X=0,qquad L_X=10,qquad
\partial=Q(A_8)-Q(X)=12.
```

For either orientation the endpoint capacity and tolled residual are

```math
c=[2L_X-(Q(X)-\sigma h_X)]_+=12,
\qquad [c-\partial]_+=0.
\tag{R16.26}
```

Thus the root tolled residual is zero for every endpoint pair.  Conservation
leaves zero descendant allocation, so the allocation-weighted tolled resource
is zero on every endpoint tree, exactly as in (10.538).  Some unused
descendant buckets do have positive raw residual; no claim is made that all
raw descendant residuals vanish.

### Endpoint no-go theorem

**Verified.**  There is no universal rule which, for every common mosaic,

1. extracts nonnegative service of total at least `\rho_{\rm mix}` (and hence
   none of total at least `\rho_{\rm pure}` or `X`), and
2. assigns it either to root endpoint buckets or through the existing
   conserved endpoint allocation, while keeping bucket load bounded by any
   finite multiple of the corresponding tolled resource.

Indeed, (R16.25) demands positive total `4`, while (R16.26) makes every root
bucket zero and conservation makes every allocation-weighted descendant
bucket zero.  This counterexample is stronger than a sign failure: the dual
marginals are already nonnegative.

The theorem concerns the actual residual resource entering (10.543).  It
does **not** rule out inventing a new temporal response resource through an
additional, separately proved conservation identity.  The response LP by
itself supplies no such identity.

## 5. The missing localization is exactly another Hall problem

Suppose nonnegative block totals `s_i` have somehow been selected and a
matrix-derived relation says which endpoint buckets `b` can receive each
block's mass.  Let the actual bucket capacities be `R_b`.  A capacity-respecting
assignment `y_{ib}` must satisfy

```math
\sum_b y_{ib}=s_i,qquad
\sum_i y_{ib}\le R_b,qquad
y_{ib}=0\ \text{off the compatibility relation}.
\tag{R16.27}
```

Max-flow/min-cut gives the exact separation condition

```math
\boxed{
\sum_{i\in I}s_i\le R(N(I))
\quad\text{for every }I\subseteq[k].
}
\tag{R16.28}
```

The uncovered mass is

```math
\max_{I\subseteq[k]}[s(I)-R(N(I))]_+.
```

The response dual provides neither the compatibility relation nor any of the
proper-subset inequalities.  Its scalar objective can at most imply a
full-set/root total.  A8 violates even that full-set endpoint condition.

After supplies are placed at temporal nodes, honest descendant service is a
second flow problem.  Its exact cuts are the subtree/antichain cuts
(10.540)--(10.543).  Therefore saying “place the response marginal where it
is needed” without proving (R16.28) and every temporal subtree inequality is
precisely the desired response-to-temporal theorem restated in flow language.

## 6. A5/A7/A8 timing and neutrality audit

### 6.1 A5: endpoint absence and the opposite timing wall

For the order-five minimizer in (10.518), the endpoint pair with positive
ground `\mathbf1` and negative ground `(1,1,1,1,-1)` has zero tolled residual
on both the `4`-vertex and singleton shores.  Nevertheless the deterministic
singleton deletion has

```math
8-\frac{64\sqrt5}{25}>0
```

centered demand and zero terminal excess.  Thus an endpoint-capacity map has
no resource at that node.

The separate deletion in (10.563)--(10.564) has positive residual `r=2` and
a zero-deficit negative-ground completion, but its temporal demand is still

```math
6-\frac{24\sqrt{15}}{25}>0.
```

So a favorable signed completion, like a favorable static deficit marginal,
cannot be reclassified as negative temporal credit by sign alone.

Under ancestor-resource timing, the first positive prefix from (10.518)
cannot use a later credit even though the complete normalized gap is
negative.  This is the prefix obstruction (10.545).  A static response dual
has no timestamp and therefore cannot resolve it.

### 6.2 A7: a proper descendant cut fails although the root succeeds

For the deterministic `7\to6\to5` chain in (10.544),

```math
x_1=10-\frac{108\sqrt{42}}{49}<0,
```

```math
x_2=\frac{108\sqrt{42}-90\sqrt{35}}{49}-2>0,
```

but `x_1+x_2<0`.  The order-six suffix has no descendant negative credit,
zero terminal excess, and zero matching endpoint residual.  Hence its proper
descendant Hall cut has positive deficiency `x_2`, although the root/global
sum is harmless.  Any block decomposition that proves only a total scalar
inequality is therefore insufficient for (10.543).

The A5 and A7 chains have opposite timing: A5 needs a later resource to pay an
earlier demand, while A7 would need an earlier resource to pay a later demand.
Choosing one causal direction does not remove both walls.  Contracting both
events into one window removes the proper cut, but then one has returned to a
global sink on that window and must prove tail-summable boundary errors.

### 6.3 A8: neutrality survives nonnegative response marginals

The original A8 wall already has a positive parent-averaged singleton demand,
exact order-seven children (zero terminal excess), and zero root endpoint
residual for every endpoint choice.  The new partition calculation
(R16.22)--(R16.26)
adds that a positive common-mosaic response and nonnegative dual block
marginals coexist with this complete endpoint neutrality.  Therefore the
failure is not caused by choosing a dual with the wrong sign.

The full-shore event in (10.539) alone does not rule out cancellation from
other partial-shore outcomes, as the ledger notes.  The singleton-eligible
parent average removes that particular escape.  An asymptotic theorem may
still pay a finite A8 loss as additive error or exclude this configuration
structurally; no globally minimizing amplification is known.

## 7. Exact scope of the no-go and surviving target

The support-function dual gives a useful diagnostic, but its mathematical
content stops at separation in signed-state profile space:

```math
\text{common state distribution}
\Longrightarrow
\text{signed block totals minus global slack}.
```

It does not give

```math
\text{endpoint pair}
+\text{bucket capacity}
+\text{temporal location}
+\text{proper Hall cuts}.
```

The order-nine wall shows that a common distribution loses part of the pure
response and can return `r_i=x_i` purely by sign cancellation.  Persistent
grounds show that even the exact pure response may be entirely a falling
local baseline while `Q` is constant.  The order-eight wall then proves that
nonnegative margins of this kind need not dominate any actual endpoint
resource.

A viable positive theorem must therefore add at least one genuinely new
ingredient:

1. an explicit nonnegative response measure with a proved matrix-derived
   compatibility kernel to endpoint/temporal nodes;
2. capacity domination and all subset/subtree Hall inequalities, not only
   the root total;
3. a temporal cocycle explaining how response mass enters with the correct
   sign, possibly together with localized terminal excess and negative
   credits; and
4. treatment of A5 prefixes, A7 suffixes, and A8 zero-capacity parents, or an
   explicit tail-summable error/structural hypothesis excluding them.

Whether such an enriched, windowed asymptotic bridge exists is **Open**.  The
finite LP dual alone cannot supply it.
