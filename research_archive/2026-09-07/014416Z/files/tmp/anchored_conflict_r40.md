# Wave 40: anchored-conflict four-corner audit

## Verdict

Exact signing minimality gives useful signed four-corner inequalities, but none
contains a positive multiple of anchored Hamming conflict.  Even exact
complement maximizers can have nonzero common-overlap disagreement and zero
cyclic/normal-fan gap.  Thus no universal positive conflict modulus follows
from the displayed parent deficits, complement cap slacks, row costs, or
ordinary cyclic monotonicity alone.

This is a statement about arbitrary selected witnesses.  The finite examples
do **not** rule out choosing different witnesses from the same fibers so as to
obtain a globally coherent low-conflict kernel.

## Exact parent four-corner identity

Let `A` be an exact minimizer with `Q(A)=q`.  Fix a common response orientation
`sigma`, selectors `S,T` containing an anchor `v`, and gauge
`x_v=y_v=1`.  Put

```math
E_x=\sigma x^TAx=q-\Delta_x,
\qquad E_y=\sigma y^TAy=q-\Delta_y.
```

Let

```math
K=\{i:x_i\ne y_i\},\qquad
H=(K\cap S\cap T)\setminus\{v\},\qquad
R=K\setminus H,
```

and write `w_{ij}=sigma a_{ij}x_ix_j`.  Define

```math
J=\sum_{i\in H,\,j\in R}w_{ij},
\qquad
P=\sum_{i\in H,\,j\notin K}w_{ij}.
```

For the signed parent cuts `C_x(H),C_y(H)`, exactly

```math
C_x(H)=J+P,\qquad C_y(H)=J-P.
```

If `x^H,y^H` denote crossover words obtained by flipping `H`, then

```math
E_A(x^H)+E_A(y^H)=E_x+E_y-8J.
```

The four individual parent-cap inequalities are equivalently

```math
-\frac{\Delta_x+\Delta_y}{8}
\le J\le
\frac q2-\frac{\Delta_x+\Delta_y}{8},
```

and

```math
\left|P-\frac{\Delta_y-\Delta_x}{8}\right|
\le J+\frac{\Delta_x+\Delta_y}{8},
```

```math
\left|P-\frac{\Delta_y-\Delta_x}{8}\right|
\le \frac q2-J-\frac{\Delta_x+\Delta_y}{8}.
```

These are sharp signed cones.  They have no term depending on `|H|` or on
the corrected conflict `sum_(i in H) 1/p_i`.  Singleton flips similarly give

```math
\frac{C_x(\{i\})+C_y(\{i\})}{2}
\ge-\frac{\Delta_x+\Delta_y}{8}\qquad(i\in H),
```

so weighting by `1/p_i` only places conflict beside the deficits in a lower
bound; it does not give positive Hamming coercivity.

## Exact complement four-corner identity

Use the complement signing

```math
B_S=-A+2P_SAP_S,
```

and set

```math
L_x=\langle B_S,d_x\rangle,\quad
L_y=\langle B_T,d_y\rangle,\quad
M_S=Q(B_S),\quad M_T=Q(B_T),
```

```math
\Gamma_x=M_S-L_x,\qquad \Gamma_y=M_T-L_y.
```

Let

```math
D_x=C_{B_S,x}(H),\qquad D_y=C_{B_T,y}(H),
\qquad
\Xi=\frac{D_x+D_y}{2},\quad
\Omega=\frac{D_x-D_y}{2}.
```

The four complement-cap inequalities are exactly

```math
-\frac{\Gamma_x}{4}\le D_x\le
\frac{M_S}{2}-\frac{\Gamma_x}{4},
\qquad
-\frac{\Gamma_y}{4}\le D_y\le
\frac{M_T}{2}-\frac{\Gamma_y}{4}.
```

Equivalently, their two centered forms are

```math
\left|\Omega-\frac{\Gamma_y-\Gamma_x}{8}\right|
\le \Xi+\frac{\Gamma_x+\Gamma_y}{8},
```

```math
\left|\Omega-
\left(\frac{M_S-M_T}{4}+\frac{\Gamma_y-\Gamma_x}{8}\right)\right|
\le
\frac{M_S+M_T}{4}-\Xi-\frac{\Gamma_x+\Gamma_y}{8}.
```

The cross-overlap quantity in the first cone has the exact expansion

```math
\Xi=
\sum_{i\in H,\,j\in R}w_{ij}
  \bigl(\mathbf1_{j\in S}+\mathbf1_{j\in T}-1\bigr)
+\sum_{i\in H,\,j\notin K}w_{ij}
  \bigl(\mathbf1_{j\in S}-\mathbf1_{j\in T}\bigr).
```

Again the inequalities are signed and can be saturated at `Xi=0` with
nonzero `H`, even when `Gamma_x=Gamma_y=0`.

## Parent-row crossover

Write `u=x 1_H`, `r=x 1_R`, and `c=x 1_(K^c)`.  Thus
`x=c+u+r`, `y=c-u-r`, `x^H=c-u+r`, and `y^H=c+u-r`.  Exactly,

```math
R_2(x^H)+R_2(y^H)
=R_2(x)+R_2(y)-8\langle Au,Ar\rangle,
```

and

```math
\|A(u+r)\|_2^2\le\frac{R_2(x)+R_2(y)}2.
```

Although

```math
\|Au\|_2^2=(n-1)|H|
+2\sum_{i<j\in H}x_ix_j(A^2)_{ij},
```

the Gram correction and `\langle Au,Ar\rangle` have no favorable sign.
Consequently the row identities also provide no Hamming modulus.

## Cyclic monotonicity and its exact failure mode

Suppose a canonical rule chooses `d_S` by maximizing
`<B_S,d>-Phi(d)` with one selector-independent penalty `Phi`.  Adding the
two optimality inequalities gives only

```math
\langle B_S-B_T,d_S-d_T\rangle\ge0.
```

For the same orientation and `K={i:x_i ne y_i}`, the left side is exactly

```math
8\sum_{\substack{i<j\\|\{i,j\}\cap K|=1}}
\sigma a_{ij}x_ix_j
\left(\mathbf1_{\{i,j\}\subset S}
-\mathbf1_{\{i,j\}\subset T}\right).
```

It sees only edges that both cross the global disagreement cut and lie in
the selector-support difference; it does not see `|H|` directly.

More generally define the two-way cross regret

```math
G_{S,T}=
\bigl[M_S-\langle B_S,d_T\rangle\bigr]
+\bigl[M_T-\langle B_T,d_S\rangle\bigr].
```

Then exactly

```math
G_{S,T}=\Gamma_S(d_S)+\Gamma_T(d_T)
+\langle B_S-B_T,d_S-d_T\rangle.
```

For exact maximizers the cap slacks vanish, so cross regret is precisely the
cyclic gap.  The order-nine example below has `G_(S,T)=0` and nonzero
overlap conflict.

## Exact finite witnesses

At order five, use

```text
A5 =
 0 -1 -1 -1 -1
-1  0 -1 -1  1
-1 -1  0  1 -1
-1 -1  1  0  1
-1  1 -1  1  0
```

Exhaustive enumeration gives `q_5=8`, so this is an exact minimizer.  Take
`S=(0,1,2,3)`, `T=(0,1,2,4)`, anchor `0`, and `sigma=+1`.

1. With
   `x=(1,1,-1,-1,-1)` and `y=(1,-1,1,-1,-1)`, one has
   `Delta=(0,0)`, `Gamma=(0,0)`, complement caps and energies `(8,8)`,
   rows `(24,24)`, `H=K={1,2}`, and `J=0`.  Thus the parent lower cone is
   flat despite maximal disagreement on the nonanchor overlap.  Here
   `Xi=4` and the cyclic gap is `32`.

2. With
   `x=(1,1,-1,-1,1)` and `y=(1,-1,1,1,-1)`, the same deficits, slacks,
   caps, energies, rows, and `H={1,2}` occur, while
   `K={1,2,3,4}`, `Xi=0`, `J=4`, and the cyclic gap is `16`.
   Thus the complement lower cone is flat despite maximal overlap
   disagreement.

The checker exhausts all exact minimizers at `n=3,4` on the proper
high-ratio slice `m=n-1` and finds no pair of distinct selectors with common
anchor and orientation, simultaneous parent exact grounds and complement
exact maximizers, and nonanchor-overlap disagreement.  Order two is
structurally vacuous.  Therefore `n=5` is minimal for this **narrowly scoped
zero-slack/exact-ground pattern**, not for every conceivable conflict
counterexample.

The cyclic obstruction needs order nine in the present audit.  On the known
exact minimizer `A9`, take

```text
S=(0,1,2,3,4,5,6),  T=(0,1,2,3,4,5,7),
x=(1,-1,1,1,1,1,-1,1,1),
y=(1, 1,1,1,1,1,-1,1,1),  sigma=+1, anchor=0.
```

The exact data are

```text
q=24; complement caps=(28,28); complement energies=(28,28);
Gamma=(0,0); parent energies=(8,-8); Delta=(16,32);
rows=(112,112); H=K={1}; J=Xi=row-cross=0; cyclic gap=0.
```

Hence no inequality with a universal `kappa>0` of the form

```math
\kappa\,c_{S,T}
\le \langle B_S-B_T,d_S-d_T\rangle
+\text{(parent/complement scalar slacks)}
```

can hold for arbitrary exact selected witnesses when the displayed slacks
vanish.  This pair does not show that all possible choices from the two
maximizer fibers conflict.

## Entropy comparison and the surviving added hypothesis

For a deterministic anchored selection, put
`p_i=Pr(i in S)` and
`m_i=E[y_i^S | i in S]`.  Its expected corrected conflict is exactly

```math
\mathcal C_v
=\frac12\sum_{i\ne v}p_i(1-m_i^2).
```

For two independent binary responses with means `m_S,m_T`, exactly

```math
\Pr(Y_S\ne Y_T)
=\frac14\bigl[(1-m_S^2)+(1-m_T^2)+(m_S-m_T)^2\bigr].
```

Thus an entropy/KL control of marginal drift alone misses within-selector
depolarization; conversely, polarization alone misses drift.  Both pieces
must be priced for randomized kernels.

The weakest honest quantitative addition exposed by this audit is a
normal-fan/cross-overlap modulus **for the selected response**, together with
an explicit flat-face price:

```math
\kappa_n c_{S,T}\le G_{S,T}+F_{S,T},
\qquad F_{S,T}\ge0,
```

followed by an affordable expectation bound on `G+F`.  The `A9` witness
shows that either `F` must charge disagreement inside a common exposed face,
or the selection rule must forbid such a pair.  A merely
selector-independent penalty proves nonnegative cyclic monotonicity and is
insufficient.  A strongly convex response on the convexified cut polytope
could provide a modulus, but one would still have to preserve complement
thresholds, row cost, and control integral rounding.  Alternatively, one may
couple the exact-maximizer fibers and lexicographically minimize global
anchored conflict; proving that this coupled optimum is affordable remains
the open coherent-selection problem relevant to (10.1089).

## Verification

`tmp/anchored_conflict_r40_check.py` verifies all displayed finite energies,
slacks, cuts, crossover identities, row identities, and cyclic formulas with
integer arithmetic.  It also performs the stated exhaustive `n=3,4`
minimality check.  Current output ends with:

```text
PASS anchored_conflict_r40_check
```
