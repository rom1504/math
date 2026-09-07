# Wave 55: exact-minimizer localized curvature becomes a near-ground migration theorem

## Outcome

The hidden-curvature configuration is not pointwise excluded.  Instead, exact
discrete minimality forces a sharp and genuinely subproject-scale migration of
the parent witness.

Let `H` be any edge set, with `e=|H|` and `v` incident vertices, and put

```math
\eta_H=\sqrt{8ev\log 2}.
```

If all parent states in a deficit-`t` layer have positive signed correlation
at least `lambda e` with the current signs on `H`, attenuating precisely those
edges and then sign-rounding gives a concrete integral replacement whose cap
is smaller than the exact minimum.  Consequently an exact minimizer must have
a deficit-`t` state whose `H`-correlation is smaller than `lambda e` whenever

```math
0<\lambda\le1,qquad \frac{2\eta_H}{\lambda}<t\le4e.
```

For a child ground, one can choose `H` to contain only child-positive edges
and to have size comparable to the sum of any selected internal fields.  A
Lorentz-prefix argument shows that `I >> n^(9/4-c)` produces such an `H` with
`eta_H=o(e)` and `eta_H=O(n^(5/4))=o(T_n)`.  Thus a hidden spike in an exact
minimizer necessarily produces an `o(T_n)`-deficit parent state which
disagrees with the child on a fixed positive fraction of the selected excess
edges.

This is the exact obstruction as well as the exact positive theorem.  The
small quantities `X`, `d`, and `H_T` constrain the single child word `y`; they
do not constrain the migrated parent state.  A box proof now needs a
common-witness/fibre-compatibility theorem, or an aggregate way to charge the
migrated states.  Stored A9/A10 hidden spikes already take the migration branch
at parent deficit zero, so universal pointwise compatibility is false.

The derivations and finite checks are reproduced by
`tmp/localized_curvature_r55_check.py`; saved output is in
`tmp/localized_curvature_r55_check.out`.

## 1. Exact attenuated-block replacement theorem

Let `A` be an exact order-`n` minimizer and `q=Q(A)=q_n`.  On the lifted
oriented state space write

```math
E_A(\omega)=\tau x^{\mathsf T}Ax,qquad
\Delta_A(\omega)=q-E_A(\omega),qquad
s_e(\omega)=\tau a_ex_ix_j.
```

For an edge set `H`, define

```math
S_H(\omega)=\sum_{e\in H}s_e(\omega),qquad
e=|H|,qquad v=|V(H)|,qquad
\eta_H=\sqrt{8ev\log2}.
```

The following statement is exact and uses global discrete minimality.

> **Attenuated-block replacement theorem.**  Fix `0<lambda<=1` and
> `0<t<=4e`.  If
>
> ```math
> \Delta_A(\omega)\le t\quad\Longrightarrow\quad
> S_H(\omega)\ge\lambda e
> \tag{W55.1}
> ```
>
> for every oriented parent state, then there is an integral replacement of
> the signs on `H`, leaving every other edge fixed, whose cap obeys
>
> ```math
> \boxed{
> Q(A')\le q-\frac{\lambda t}{2}+\eta_H.}
> \tag{W55.2}
> ```

In particular, if `lambda t/2>eta_H`, (W55.2) contradicts exact
minimality.  Equivalently, every exact minimizer satisfies the **Verified
localized migration theorem**

```math
\boxed{
\frac{2\eta_H}{\lambda}<t\le4e
\quad\Longrightarrow\quad
\exists\omega:\
\Delta_A(\omega)\le t,\quad S_H(\omega)<\lambda e.}
\tag{W55.3}
```

### Proof and constant audit

Replace the signs on `H` provisionally by the fractional coefficients

```math
z_e=\left(1-\frac{t}{4e}\right)a_e.
\tag{W55.4}
```

The restriction `t<=4e` keeps (W55.4) on the line segment from `a_H` to
zero.  The oriented energy becomes

```math
E_z(\omega)
=q-\Delta_A(\omega)-\frac{t}{2e}S_H(\omega).
\tag{W55.5}
```

If `Delta<=t`, (W55.1) and (W55.5) give
`E_z<=q-lambda t/2`.  If `Delta>t`, the universal bound `S_H>=-e` gives

```math
E_z(\omega)<q-t+\frac t2=q-\frac t2
\le q-\frac{\lambda t}{2}.
```

Hence the fractional cap is at most `q-lambda t/2`.

Independently round each `z_e` to a sign `zeta_e` with mean `z_e`.  For a
fixed local oriented pattern, the rounding error is

```math
2\sum_{e\in H}(\zeta_e-z_e)\tau x_ix_j.
```

Each summand has range length four.  There are at most `2^v` different
oriented edge patterns on the `v` incident vertices.  Hoeffding and a union
bound therefore give, for every `epsilon>0`, a rounding whose cap is at most
the fractional cap plus `eta_H+epsilon`.  Choose `epsilon` smaller than
`lambda t/2-eta_H`; this proves (W55.2).  The rounding can also be found by
finite enumeration or conditional expectation, so the replacement is an
actual complete signing, not a convex relaxation.

For example, if `eta_H/e<2/3`, take `lambda=1/2` and `t=6eta_H`.  Then an
exact minimizer has a state satisfying

```math
\Delta_A(\omega)\le6\eta_H,qquad S_H(\omega)<\frac e2.
\tag{W55.6}
```

If all edges of `H` are positive in a reference child orientation and
`N_H(omega)=|{e in H:s_e(omega)=-1}|`, then
`S_H=e-2N_H`; (W55.6) says `N_H>e/4`.

## 2. A field spike supplies a localized excess-edge block

Fix `S`, `|S|=m`, and an oriented child ground `(sigma,y)`.  Write

```math
c_{ij}=\sigma a_{ij}y_iy_j,qquad
r_i=\sum_{j\in S\setminus\{i\}}c_{ij}\ge0,qquad
I=\sum_{i\in S}r_i^2.
```

For any vertex set `U subseteq S`, put `R_U=sum_(i in U) r_i`.  Vertex `i`
has exactly `(m-1+r_i)/2` child-positive incident edges.  Since
`0<=r_i<=m-1`, choose any `r_i` of them and call this set `H_i`.  With

```math
H=\bigcup_{i\in U}H_i
```

every edge of `H` is child-positive, and each edge is selected at most
twice.  Therefore the **Verified excess-edge extraction** is

```math
\boxed{
\frac{R_U}{2}\le |H|\le R_U,qquad |V(H)|\le m.}
\tag{W55.7}
```

This removes the positive-star baseline: only `r_i`, the positive-minus-
negative excess, determines how many edges are retained.

There is always a concentrated choice of `U`.  Sort
`r_(1)>=...>=r_(m)` and write `R_h=sum_(j<=h)r_(j)`.  If

```math
K=\max_{1\le h\le m}\frac{R_h}{\sqrt h},
```

then `r_(h)<=R_h/h<=K/sqrt(h)`.  Hence, with
`H_m=sum_(h<=m)1/h`,

```math
I\le K^2H_m.
```

Choosing a maximizing prefix gives the **Verified Lorentz localization**

```math
\boxed{
R_U\ge\sqrt{\frac{Ih}{H_m}}
\ge\sqrt{\frac I{H_m}}.}
\tag{W55.8}
```

Combine (W55.7)--(W55.8).  Along a fixed-density sequence with
`0<c<1/4` and

```math
I\gg n^{9/4-c},
```

the extracted set obeys

```math
e\ge\frac12\sqrt{I/H_m}\gg n,qquad
\frac{\eta_H}{e}le
\sqrt{\frac{8m\log2}{e}}=o(1).
\tag{W55.9}
```

Also `sum_i r_i=Q(A[S])<=q_n=O(n^(3/2))`: the principal norm is at most
the full norm by uniform completion.  Thus `e<=q_n`, and

```math
\boxed{
\eta_H\le\sqrt{8nq_n\log2}=O(n^{5/4})
=o(n^{3/2-c})=o(T_n).}
\tag{W55.10}
```

For every fixed `lambda>0`, (W55.9) eventually permits
`t=3eta_H/lambda` in (W55.3).  Therefore an exact minimizer with the
postulated spike has a parent state with

```math
\boxed{
\Delta_A=o(T_n),qquad
N_H>\frac{1-\lambda}{2}e.}
\tag{W55.11}
```

In words: at project-negligible parent deficit, the response must reverse a
constant fraction of a child-positive edge set whose cardinality is forced
by the localized curvature.  This uses exact minimality in an essential way
and is much sharper in deficit than applying one simultaneous flip to all of
`H`, which would give only `Delta=O(e)`.

## 3. What (W55.11) proves, and the exact remaining compatibility gap

The theorem does **not** prove that the box value is small.  It proves the
following precise conditional exclusion.

> If the low-cross hidden-spike incidence additionally had the property that
> every parent state of deficit at most `t` retained signed correlation at
> least `lambda e` on its extracted excess-edge set, then (W55.4), followed
> by sign rounding, would lower the parent cap by at least
> `lambda t/2-eta_H`.  Such an incidence cannot occur in an exact minimizer
> when that quantity is positive.

The missing hypothesis is genuinely about a **common parent witness**.  The
conditions

```math
X=\lVert By\rVert_2^2\ \hbox{small},\qquad
d=BPy+CBy\ \hbox{small},\qquad
\mathsf H_T\ \hbox{small}
```

refer only to the fixed word `y`.  They say nothing about `s_e(omega)` for a
different parent near-ground state.  Thus neither (10.1298)--(10.1303) nor
the new replacement theorem bridges the child fibre to the full active
parent face.

A viable continuation must do at least one of the following:

1. select a low-cross child incidence whose excess-edge pattern is stable on
   the relevant parent near-ground face;
2. charge every migrated state in (W55.11) to another favorable child fibre
   and prove that repeated migration terminates at a box-canceling incidence;
3. aggregate many choices of the excess sets `H_i` so that witness migration
   forces saved near-ground mass rather than merely one migrating state.

A pointwise claim that every low-cross hidden spike has the first property is
already false finitely.

## 4. Exact finite diagnostics and the scalable obstruction

### Stored exact A9 and A10 hidden spikes

For the stored `A9,m=8` and `A10,m=8` incidences with

```math
I=120,qquad X=0,qquad d=0,qquad\mathsf H_T=0,
```

the deterministic top-prefix/excess-edge construction gives respectively

| instance | `R_U` | `e` | `eta_H` | minimum escape deficit | ground-state negative-edge range |
|---|---:|---:|---:|---:|---:|
| A9 | 20 | 15 | 25.795761150576 | 0 | 0--7 |
| A10 | 20 | 17 | 27.461684807186 | 0 | 0--9 |

These small orders do not satisfy `eta_H=o(e)`, so they are not tests of the
asymptotic constant in (W55.3).  They are nevertheless **Numerical exact
finite diagnostics** of the missing compatibility: the parent ground face
contains both a state with zero disagreements and states with more than one
quarter of `H` reversed.  Witness migration is present even with
`X=d=H_T=0`.

### Wave 54 scalable hidden-spike template

The scalable pair-duplicated hub construction remains a valid obstruction to
all generic block algebra, but it is not an exact minimizer and hence does not
contradict (W55.3).  The finite order-22 realization makes this distinction
concrete.  Exact enumeration gives

```math
I=1840,qquad X=d=\mathsf H_T=0,qquad Q(A_{\rm hub})=170.
```

Keeping its cross block and conference exterior fixed, replace the internal
pair-duplicated hub block by the pair duplication of the stored exact `A8`.
This is a concrete internal signing replacement, and exhaustive projective
enumeration gives

```math
\boxed{Q(A_{\rm replacement})=118<170.}
\tag{W55.12}
```

Thus the finite template fails exact minimality spectacularly.  Equation
(W55.12) is not claimed as a uniform replacement theorem for the full
scalable family; that stronger assertion is unnecessary and unproved.

## 5. Research judgment

The requested generic implication

```math
I\gg n^{9/4-c},\quad X,d,\mathsf H_T\ \hbox{small}
\quad\Longrightarrow\quad\hbox{box cancellation}
```

is not established.  What is now exact is its minimizer-specific replacement
content: a spike at that scale forces subproject-deficit migration on a
field-excess edge block, and it can persist only because the active parent
witness leaves the child fibre.  This identifies the missing theorem without
reusing trace or operator arguments.

The most promising next form is an **aggregate compatibility theorem**, not
a pointwise one: use the many admissible choices of `H_i` in (W55.7), count
how many can share the same escaping near-ground parent state, and transfer
the resulting mass back to low-cross child incidences.  The stored finite
examples require such a theorem to tolerate migration inside the parent
ground face.
