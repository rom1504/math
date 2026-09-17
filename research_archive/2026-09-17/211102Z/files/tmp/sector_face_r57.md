# Wave 57 sector-face audit: exact disagreement identity and a pointwise wall

## Scope and conclusion

I attacked the gap between the common active-face law (10.1340)--(10.1343)
and bare favorability (10.1352), keeping the two orientations separate.  There
is a useful exact identity: an escaping parent differs from the chosen child
ground on at least one quarter of the selected child-positive edges, in either
orientation sector, and its complete local deficit is exactly four times the
**signed** child weight of the whole sector-disagreement set.  Thus the missing
property is cancellation by child-negative disagreement edges, not merely a
larger positive-edge reversal count.

This identity also exposes a no-go.  Flipping the complete disagreement set is
the canonical state-adapted second perturbation, but it raises that parent's
energy by its local deficit.  It therefore supports, rather than contradicts,
exact signing minimality.  A stored exact order-nine minimizer gives a stronger
finite wall: one parent ground escapes more than one third of all selector
blocks while none of its canonical escaped incidences is bare-favorable at
zero slack.  It satisfies the full-edge energy identity, has a rank-one PSD
orientation matrix, uses exact child grounds, and has no parent-sector
cancellation.

This is an **exact finite obstruction to pointwise and identity-only
implications**.  It is **not** an asymptotic counterexample at the target
`T_n` scale.  A surviving theorem must be aggregate and scale-sensitive.

## 1. Unified identity in the two relative orientation sectors

Let `A` be a signing, `q=Q(A)`, and let the oriented parent state be

```math
\omega=(\tau,x),\qquad
E_A(\omega)=\tau x^{\mathsf T}Ax=q-\Delta.
```

For a selector `S`, choose an exact oriented child ground `(sigma,y)`, so

```math
q_S=Q(A[S])=\sigma y^{\mathsf T}A[S]y.
```

Gauge the internal edges and the parent restriction by

```math
c_{ij}=\sigma a_{ij}y_iy_j,\qquad
z_i=x_iy_i,\qquad
\kappa=\tau\sigma.
```

Thus `2 sum_(i<j in S)c_ij=q_S`.  The parent edge sign is

```math
s_{ij}=\tau a_{ij}x_ix_j=c_{ij}\kappa z_iz_j.
```

Define the **sector-disagreement set**

```math
D_S(\omega,y)=\{ij\subset S:\ \kappa z_iz_j=-1\}.
```

When `kappa=+1`, this is the cut of `z`.  When `kappa=-1`, it is the
uncut of `z`.  On `D`, the parent sign is `s_e=-c_e`; off `D`, it is
`s_e=c_e`.  Because quadratic forms count every undirected edge twice,
the complete local deficit is exactly

```math
\boxed{
\delta_S(\omega)
=q_S-\tau x_S^{\mathsf T}A[S]x_S
=4\sum_{e\in D_S(\omega,y)}c_e.
}
\tag{SF.1}
```

This audits the orientation and the factor four simultaneously.  Exact child
optimality gives the nonnegativity of the final signed sum.

Let `E_S` be any legal Wave 56 excess block: every edge in it is
child-positive (`c_e=+1`), and write `e_S=|E_S|`.  Put
`a_S=|E_S cap D_S|`.  Then

```math
\boxed{
S_{E_S}(\omega)=\sum_{e\in E_S}s_e=e_S-2a_S.
}
\tag{SF.2}
```

Consequently the common-law escape threshold has the same interpretation in
both sectors:

```math
S_{E_S}(\omega)\le e_S/2
\quad\Longleftrightarrow\quad
a_S\ge e_S/4.
\tag{SF.3}
```

For `kappa=+1`, at least a quarter of the selected positive edges cross the
relative cut.  For `kappa=-1`, at least a quarter lie in its relative uncut.
The earlier apparently asymmetric description disappears after the correct
orientation gauge.

If `P_D` and `N_D` denote the numbers of child-positive and child-negative
edges in `D`, respectively, then (SF.1) is

```math
\delta_S=4(P_D-N_D),\qquad P_D\ge a_S\ge e_S/4
\quad\hbox{on escape}.
\tag{SF.4}
```

Bare favorability at threshold `t` is exactly

```math
4(P_D-N_D)
\le B_{n,m}+p_2\Delta+t.
\tag{SF.5}
```

Thus escape becomes favorable only when sufficiently many child-negative
edges join the sector-disagreement set.  Counting positive reversals alone
cannot establish (SF.5).

## 2. The canonical second flip has the wrong sign

Let `A^D` be obtained by flipping precisely the internal edges in `D`.  On
those edges `s_e=-c_e`, and (SF.1) gives

```math
\sum_{e\in D}s_e=-\frac{\delta_S}{4}.
```

Flipping one undirected signing edge changes the oriented quadratic energy by
`-4s_e`.  Therefore the same parent state has

```math
\boxed{
E_{A^D}(\omega)
=q-\Delta-4\sum_{e\in D}s_e
=q-\Delta+\delta_S.
}
\tag{SF.6}
```

In particular, an unfavorable large `delta_S` makes `A^D` visibly *less*
competitive.  Exact minimality only says `Q(A^D)>=q`, so (SF.6) satisfies it
with positive slack.  The state-dependent perturbation naturally suggested
by the exposed child face cannot provide the proposed contradiction.  A
useful second perturbation would have to be common to many incidences or have
an independent cap upper bound; neither follows from (10.1340).

The orientation-conditioned PSD matrices do not repair this pointwise gap.
For a law `mu`, they are

```math
M_\tau=\mathbb E[xx^{\mathsf T}\mid\tau],\qquad M_\tau\succeq0,
\qquad \operatorname{diag}M_\tau=\mathbf1.
```

They linearly encode mixtures of cut matrices and the global exposed parent
face.  The finite wall below already has one orientation only and
`M_-=xx^T`, a rank-one correlation matrix.  Hence PSD, the full-energy
identity, and exact exposed child grounds cannot by themselves imply a
pointwise favorable alternative.

## 3. Exact `A_9,m=7` obstruction

Use the stored exact order-nine minimizer `A_9`, whose enumerated cap is
`q=24`.  Take the parent ground

```text
tau = -1,        x = +-++-----
```

(coordinates are indexed `0,...,8`).  Exact enumeration gives

```math
\tau x^{\mathsf T}A_9x=24,\qquad
\Delta=0,\qquad
x^{\mathsf T}A_9^2x=96,
\qquad \sum_es_e=12=q/2.
\tag{SF.7}
```

For each of the `binom(9,7)=36` selectors, choose the lexicographically first
projective exact child ground and, at each vertex, the lexicographically first
`r_i` child-positive incident edges, then take their union.  Every choice is a
legal block from (10.1324).  The parent ground has

```text
canonical escaped selectors        28 / 36
bare-favorable selectors at t=0     5 / 36
escaped and bare-favorable          0 / 36
```

The escaped incidences split as follows:

| relative sector `kappa` | local deficit `delta_S` | count | `hat ell=delta_S-B` |
|---:|---:|---:|---:|
| `-1` | `4`  | 8  | `18-56 sqrt(7)/9` |
| `-1` | `12` | 3  | `26-56 sqrt(7)/9` |
| `+1` | `8`  | 14 | `22-56 sqrt(7)/9` |
| `+1` | `16` | 3  | `30-56 sqrt(7)/9` |

Indeed, for `p=7/9` and `p_2=7/12`,

```math
B_{9,7}
=24\left\{(7/9)^{3/2}-7/12\right\}
=\frac{56\sqrt7}{9}-14
=2.462452602\ldots.
\tag{SF.8}
```

Here `Delta=0`, so (10.1352) is exactly
`hat ell=delta_S-B_(9,7)`.  Since `0<B_(9,7)<4`, every escaped value in the
table is strictly positive.  This explicitly checks the full effective loss,
not just `delta_S`, in both relative orientation sectors.

The example is not fragile with respect to choosing the child witness or the
`r_i` positive edges.  Exhausting **every** exact child ground and every legal
per-vertex edge selection shows that 15 of the 36 selectors escape for this
same parent under every possible choice; all 15 have `delta_S>=4`.  Thus a
Fubini-size (`>1/3`) subfamily of bad escapes is choice-robust.  The canonical
rule sharpens this to 28 bad escapes and separates all five favorable
selectors from them.

This example is stronger than the old formal all-negative response:

- it is an actual parent ground of an exact finite minimizer;
- it obeys `sum_e s_e=(q-Delta)/2` exactly;
- its orientation matrix is PSD rank one with diagonal one;
- every child witness is an exact exposed ground; and
- it uses only one parent orientation, so cross-orientation cancellation is
  not the explanation.

It nevertheless does **not** supply a fixed-density asymptotic sequence, and
at order nine a target-scale allowance `t=C T_n` could exceed the first
positive deficit step.  It therefore falsifies pointwise/identity-only
bridges, not an aggregate asymptotic theorem at scale `T_n`.

## 4. Precise surviving aggregate lemma

Let `mu` be the particular common active-face law after the near-ground
conditioning in (10.1341), and sample `(S,omega)` from `U_m tensor mu`.
For each selector use an exact child ground and its legal excess block, and
form `D_S`, `P_D`, and `N_D` above.  The least additional assertion that
actually uses the common escape mechanism is the following entropy-scale
**sector-disagreement cancellation lemma**: for some fixed constants `C,C'`,

```math
\boxed{
(U_m\otimes\mu)\left\{
S_{E_S}\le e_S/2,
\quad
4(P_D-N_D)\le B_{n,m}+p_2\Delta+C T_n
\right\}
\ge e^{-C'H}.
}
\tag{SF.9}
```

By (SF.1), the second event is exactly bare favorability.  No constant
conditional fraction is needed; an `e^{-O(H)}` fraction of the already
constant-size escape incidence suffices.  Taking `nu=mu` then gives the mass
part of the selected-prior lemma.  The remaining independent row requirement
can be written in its weakest captured form as

```math
\mathbb E\left[R(\omega)\mathbf1_{\mathrm{event\ in\ (SF.9)}}\right]
\le C n^{9/4-c}\,
(U_m\otimes\mu)(\mathrm{event\ in\ (SF.9)}).
\tag{SF.10}
```

Equations (SF.9)--(SF.10), together with the already conditioned
`Delta=o(T_n)`, feed directly into (10.1353).  The finite obstruction does not
falsify them because they are existential in the full minimax law, aggregate
over selectors and states, and allow the asymptotic `T_n` window.

What it does prove is that a derivation of (SF.9) must add genuinely joint
information: for example, an entropy bound on the family of disagreement
sets, a selector-aggregate theorem forcing child-negative cancellation, or a
common perturbation that covers many bad `D_S` at once.  Sector splitting,
PSD positivity, the total parent energy, and the child exposed face are all
already present in the counterexample and are insufficient on their own.

## Reproduction

Run:

```bash
.venv/bin/python tmp/sector_face_r57_check.py
```

The checker enumerates all projective parent and child cuts used above,
verifies (SF.1), (SF.2), and (SF.6) edge by edge, checks every factor of two,
audits the exact `hat ell` threshold, verifies the rank-one PSD and full-energy
identities, and exhausts all legal excess-block choices for the robust count.
