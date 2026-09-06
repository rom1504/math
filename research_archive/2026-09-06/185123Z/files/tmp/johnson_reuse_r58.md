# Wave 58: common escape does not force Johnson reuse from local sector data

## Conclusion and scope

The Fubini step does not by itself repair the collision gap.  Even after one
center has been fixed and made to escape on every selector, an
`e^{-O(H)}`-sized signed-cancellation subfamily can be Johnson-random.  The
complete weighted histogram, not only ordinary degrees, then stays below its
degree-two baseline.

There is new signing content behind this obstruction.  A single order-six
child signing admits two actual relative cuts which have local deficits zero
and four, while the same legal child-positive excess block escapes in both
cases.  Thus the obstruction preserves exact child optimality, the two sector
identities, legality of the Wave-56 block, and a constant child-cap profile.
It is stronger than the Wave-57 balanced-partition wall, which only assigned
different center labels to selectors.

The construction is nevertheless **local/abstract**.  The selector-wise
switched edge patterns are not restrictions of one global signing.  It does
not preserve the parent row, a global parent deficit, or exact-minimizer
status.  It therefore falsifies deductions from common escape plus the listed
local child/sector data; it does not falsify a genuinely global
signing-specific codegree theorem.

## 1. One common center can have a random cancellation family

Let `Omega=binom([n],m)`, `N=|Omega|`, and fix one center `z`.  Suppose it
escapes on all selectors, but is bare-favorable exactly on a family
`F subset Omega` of size `r`.  For the exact normalizations of (10.1362),
write

```math
P_\ell(F)=\frac1r\sum_{S,T\in F}K_\ell(S,T),
\qquad
\Pi_s(F)=\frac{E_s(F)}{rD_s}.
```

If `F` is a uniformly random `r`-subset of `Omega`, pair counting gives

```math
\boxed{
\mathbb E P_\ell
=h_\ell+(1-h_\ell)\frac{r-1}{N-1},
\qquad
\mathbb E\Pi_s=\frac{r-1}{N-1}.}
\tag{R58J.1}
```

This differs from the Wave-57 partition wall in the important quantifier:
there is only one center, and it has already survived the common-state Fubini
step.  The random label is now whether that same state receives the signed
negative-edge cancellation in (10.1359).

At the mesoscopic scale, take `r/N=e^{-Theta(H)}`.  Since
`h_ell=e^{-Theta(n)}` and `lambda_2(ell)` is polynomial, (R58J.1) implies the
existence of a deterministic `F` with

```math
P_\ell(F)=e^{-\Theta(H)}+e^{-\Theta(n)}\ll\lambda_2(\ell).
\tag{R58J.2}
```

The exact threshold inequality (10.1363) then also yields

```math
\Pi_{s_*}(F)
\le \frac{P_\ell(F)-h_\ell}{R_{\ell,s_*}}
=e^{-\Theta(H)+O(\log n)},
\tag{R58J.3}
```

exponentially below the `n^{-1/16}` term in (10.1364).  The full histogram is
already evaluated in (R58J.2); combining shells cannot repair this example.

## 2. Exact child-ground-compatible sector templates

The following matrix is used for both local records:

```math
C=\begin{pmatrix}
0&1&1&1&1&1\\
1&0&-1&-1&1&1\\
1&-1&0&1&-1&1\\
1&-1&1&0&1&-1\\
1&1&-1&1&0&-1\\
1&1&1&-1&-1&0
\end{pmatrix}.
```

Exact enumeration gives `Q(C)=10`.  The positive orientation of
`y=mathbf 1` is a child ground, with nonnegative fields

```math
C\mathbf1=(5,1,1,1,1,1).
```

The Wave-56 per-vertex rule may therefore choose the five-edge star

```math
E=\{01,02,03,04,05\}.
```

Vertex zero chooses all five positive edges and every other vertex chooses
its positive edge to zero, so this is a legal union and
`Q(C)/2=|E|=5`.

Use the positive relative sector `kappa=+1` and the two actual cut vectors

```math
z_0=(1,-1,1,1,-1,1),
\qquad
z_1=(1,-1,-1,1,1,-1).
```

Writing `D_a={ij:z_{a,i}z_{a,j}=-1}`, direct counting gives

```math
\boxed{
\sum_{e\in D_0}c_e=0,
\quad \sum_{e\in E}c_ez_{0,i}z_{0,j}=1;
\qquad
\sum_{e\in D_1}c_e=1,
\quad \sum_{e\in E}c_ez_{1,i}z_{1,j}=-1.}
\tag{R58J.4}
```

Consequently the exact sector identities give

```math
\delta_0=0,
\qquad \delta_1=4,
\qquad
2S_E(z_0)=2\le5,
\qquad
2S_E(z_1)=-2\le5.
\tag{R58J.5}
```

Both records escape, while any bare threshold `b in [0,4)` accepts exactly
the first.  Assign record zero to `S in F` and record one otherwise.  This
produces one common center escaping everywhere whose cancellation family is
the arbitrary `F` in §1.  All child caps are the same (`10`), every child
ground and disagreement set is genuine, and every block is the same legal
block.  In particular the construction does not exploit cap variation or an
illegal choice of positive edges.

For a concrete full-histogram audit take `(n,m,ell)=(11,6,3)` and the eight
selectors recorded in the checker.  Then

```math
N=462,quad d_3=1120,quad h_3=\frac1{56},
\quad\lambda_2(3)=\frac1{14},
```

and exact enumeration gives

```math
\boxed{P_3(F)=\frac3{112}<\frac1{14}.}
\tag{R58J.6}
```

The best threshold is `s_*=4`, for which

```math
D_4=180,quad R_{3,4}=\frac9{14},quad
\Pi_4(F)=\frac1{120}<\frac1{12}
=\frac{\lambda_2-h_3}{R_{3,4}}.
\tag{R58J.7}
```

Thus the exact `Pi_s` diagnostic and the complete histogram fail, despite
common escape and locally exact child-sector structure.

## 3. What child-ground exchange can and cannot add

For an actual global signing and adjacent selectors
`S=U cup {i}`, `T=U cup {j}`, extending an oriented child ground across the
new vertex with the better sign proves

```math
|q_S-q_T|\le2(m-1).
```

A fixed parent's oriented induced energy changes by at most `4(m-1)`, hence

```math
\boxed{|\delta_S(d)-\delta_T(d)|\le6(m-1).}
\tag{R58J.8}
```

This only transports a threshold through
`u=O(T_n/n)=O(n^{1/2-c})` exchanges before using the whole target slack.
For `ell=o(n)`, a `K_ell` partner and the optimizing threshold shell have
intersection `p^2n+o(n)`, hence Johnson distance
`p(1-p)n+o(n)=Theta(n)`.  The exchange bound therefore cannot reach the
shell used by (10.1364).  The local templates even keep `q_S` constant and
change `delta` by only four, so they obey these scalar exchange constraints.

What remains genuinely open is exactly the global glue omitted above: prove
that restrictions of one exact minimizing signing, one parent cut, and one
coherent collection of child grounds cannot place the cancellation labels
Johnson-randomly.  Such a theorem must use overlapping-edge consistency,
parent row/global energy, or exact minimality across many selectors.  Neither
child optimality one selector at a time nor the full histogram as a mere
accounting identity supplies it.

## Reproduction

Run:

```bash
.venv/bin/python tmp/johnson_reuse_r58_check.py
```

