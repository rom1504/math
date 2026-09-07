# Compatible replacement: exact metric wall and boundary-aware remainder

Scratch memo for Wave 12.  All finite claims are checked by
`tmp/verify_compatible_replacement_r12.py`; no tracked file is edited.

## 1. Uniform deficit stability is exactly a directed `Q`-metric

For same-order symmetric zero-diagonal matrices `X,Y`, put

```math
\delta_X(\sigma,z)=Q(X)-\sigma z^{\mathsf T}Xz,
\qquad \sigma\in\{\pm1\},\quad z\in\{\pm1\}^m.
```

Then the error used in (10.506) has the exact closed form

```math
\boxed{
\varepsilon(X,Y)
:=\sup_{\sigma,z}[\delta_X(\sigma,z)-\delta_Y(\sigma,z)]_+
=Q(X)-Q(Y)+Q(X-Y).
}
\tag{R12.1}
```

Indeed, before the positive part the difference is
`Q(X)-Q(Y)-\sigma z^T(X-Y)z`; maximizing over `(sigma,z)` adds
`Q(X-Y)`.  Triangle inequality gives
`Q(Y)\le Q(X)+Q(X-Y)`, so the resulting number is already nonnegative.
Equivalently,

```math
\boxed{Q(Y)+\varepsilon(X,Y)=Q(X)+Q(X-Y).}
\tag{R12.2}
```

Thus paying the uniform error erases the entire decrease from `Q(X)` to
`Q(Y)` and adds the `Q`-distance.  For example, if `Y` is an order-`m`
minimizer in an orbit `\mathcal O(G)`, then

```math
\boxed{
\min_{Y\in\mathcal O(G)}\varepsilon(X,Y)
=Q(X)-q_m+\operatorname{dist}_Q(X,\mathcal O(G)).
}
\tag{R12.3}
```

In particular the terminal error already contains the whole terminal
optimality excess.  It vanishes only when `X` itself is the chosen
minimizer.  If two signings differ on `k` undirected edges, the Rademacher
second-moment identity gives

```math
Q(X-Y)^2\ge
\mathbb E_z(z^{\mathsf T}(X-Y)z)^2=16k.
\tag{R12.4}
```

Consequently orbit proximity here is a genuine strong requirement, not a
formal choice of gauge.

Let `R` be a signed permutation matrix.  Exactly,

```math
\delta_{RGR^{\mathsf T}}(\sigma,z)
=\delta_G(\sigma,R^{\mathsf T}z),
\qquad
\varepsilon(X,RGR^{\mathsf T})
=Q(X)-Q(G)+Q(X-RGR^{\mathsf T}).
\tag{R12.5}
```

Uniform switching average does not repair this.  Uniform vertex switching
has `\mathbb E_R RGR^T=0`, so convexity of `Q` yields

```math
\boxed{
\mathbb E_R\varepsilon(X,RGR^{\mathsf T})
\ge Q(X)-Q(G)+Q(X)=2Q(X)-Q(G).
}
\tag{R12.6}
```

Hence a random-orbit argument cannot make the uniform errors `o(m^{3/2})`
on a macroscopic scale without additional orbit-closeness information.

For one replacement level, inserting (R12.1) into the transfer from a
hybrid successor `Y` to the original successor `X` gives

```math
Q(Y)+\mathcal L_C(Y)
\le Q(X)+\mathcal L_C(X)+Q(X-Y).
\tag{R12.7}
```

Thus (10.506) cannot exploit the lower norm of the replacement.  Along a
tower it charges

```math
\sum_t\varepsilon_t
=\sum_t\{Q(A_{t+1})-Q(\widetilde A_{t+1})
          +Q(A_{t+1}-\widetilde A_{t+1})\},
\tag{R12.8}
```

so geometric block sizes alone provide no cancellation.  The proposed
uniform-deficit route is therefore closed unless one first proves precisely
the strong suffix `Q`-proximity that the replacement was meant to create.

## 2. The exact boundary-aware quantity

For a fixed cross mosaic `C`, write

```math
F_C(z)=2\|Cz\|_1,
\qquad
\mathcal L_C(X)=\max_{\sigma,z}[F_C(z)-\delta_X(\sigma,z)]_+.
```

Since a `Q(X)`-ground makes the expression nonnegative,

```math
\boxed{
J_C(X):=Q(X)+\mathcal L_C(X)
=\max_z\{2\|Cz\|_1+|z^{\mathsf T}Xz|\}.
}
\tag{R12.9}
```

Therefore the layer difference has the exact centered decomposition

```math
\boxed{
\mathcal L_C(Y)-\mathcal L_C(X)
=Q(X)-Q(Y)+J_C(Y)-J_C(X).
}
\tag{R12.10}
```

The weakest local scalar error left after removing the compulsory norm
shift is

```math
\kappa_C(X,Y)=[J_C(Y)-J_C(X)]_+.
\tag{R12.11}
```

It is genuinely weaker than (10.506).  Orbit choice is exactly boundary
column choice:

```math
\boxed{J_C(RGR^{\mathsf T})=J_{CR}(G).}
\tag{R12.12}
```

Thus a compatible orbit representative should minimize the extension
functional `J`, not the full deficit distance.

There is also a universal singleton wall.  If an order-`m+1` global
minimizer is the singleton extension of `X` by row `c`, then
`q_{m+1}=J_c(X)`.  Replacing `X` by any order-`m` minimizer `G` gives another
signing, whence `J_c(G)\ge q_{m+1}`.  Therefore

```math
\boxed{
\mathcal L_c(G)-\mathcal L_c(X)
\ge Q(X)-q_m.
}
\tag{R12.13}
```

So even perfect boundary optimization must pay the entire child excess at a
singleton split.  This is an exact theorem, not a finite-order accident.

## 3. Exact orbit certificate: cross compatibility can evade `Q`-distance

Use the order-six/order-five matrices `A_6,D,D'` from (10.505), with
`C=(1,1,1,-1,-1)`.  The signed-permutation orbit of `D'` has 192 distinct
members.  Their exact `(Q(D-G),mathcal L_C(G))` counts are

| pair | count |
|---|---:|
| `(0,2)` | 1 |
| `(8,6)` | 5 |
| `(8,10)` | 10 |
| `(16,2)` | 11 |
| `(16,6)` | 55 |
| `(16,10)` | 110 |

Thus eleven orbit elements have `Q(D-G)=16` but
`mathcal L_C(G)=mathcal L_C(D)=2`.  For them the uniform error is 16 while
the centered boundary error is zero.  Conversely the uniform orbit average
has

```math
\mathbb E\mathcal L_C(G)=\frac{33}{4}>2,
\qquad
\mathbb E Q(D-G)=\frac{367}{24}.
\tag{R12.14}
```

This simultaneously proves that (i) uniform `Q`-distance is needlessly
strong for a fixed boundary, and (ii) averaging the orbit is not a selector;
one must minimize a boundary-dependent functional.

## 4. Exact order-seven no-go for the best boundary selector

Exhaustive enumeration of all `2^15` gauge-fixed order-seven signings gives
`q_7=18` and 3,240 gauge minimizers.  One is

```math
A_7=\begin{pmatrix}
0&1&1&1&1&1&1\\
1&0&1&1&-1&-1&1\\
1&1&0&1&-1&1&-1\\
1&1&1&0&1&-1&-1\\
1&-1&-1&1&0&-1&-1\\
1&-1&1&-1&-1&0&-1\\
1&1&-1&-1&-1&-1&0
\end{pmatrix}.
\tag{R12.15}
```

Split off vertex 4.  In tail order `(1,2,3,5,6,7)`, the boundary row is
`C=(1,1,1,1,-1,-1)`.  Its original child satisfies

```math
Q(D_6)=18,
\qquad \mathcal L_C(D_6)=0.
```

Enumeration of all 384 labeled order-six minimizers (twelve gauge classes
and all switchings) gives

```math
\boxed{
\min_{Q(G)=q_6=10}\mathcal L_C(G)=8.
}
\tag{R12.16}
```

This attains equality in (R12.13): the best lexicographic/orbit-compatible
choice still shifts the complete child excess eight into the boundary layer.

More sharply, take tail `(1,2,3,4,5)` and head `(6,7)`.  Then

```math
C=\begin{pmatrix}
1&-1&1&-1&-1\\
1&1&-1&-1&-1
\end{pmatrix},
\qquad
D_5=\begin{pmatrix}
0&1&1&1&1\\
1&0&1&1&-1\\
1&1&0&1&-1\\
1&1&1&0&1\\
1&-1&-1&1&0
\end{pmatrix}.
```

Exactly,

```math
Q(D_5)=12,
\quad \mathcal L_C(D_5)=4,
\quad J_C(D_5)=16,
```

whereas enumeration of all 192 labeled order-five minimizers gives

```math
\boxed{
\min_{Q(G)=q_5=8}\mathcal L_C(G)=12,
\qquad
\min_{Q(G)=8}J_C(G)=20>16.
}
\tag{R12.17}
```

Hence even the centered boundary functional (R12.11) cannot always be made
nonincreasing: its best possible error here is four.  This falsifies a
universal greedy or lexicographic `J_C`-monotonicity lemma, even when the
parent is globally minimizing and all minimizer switching classes are
available.

The verifier fixes the first row by switching, enumerates exactly
`2^{binom(n-1,2)}` gauges, then expands every lower-order minimizing gauge
class through all vertex switchings.  It also checks all displayed norms,
layers, orbit counts, and minima.

## 5. What survives

The exact common-mosaic identity explains the obstruction.  Partition the
original minimizer into diagonal blocks `D_j`, replace them by minimizers
`G_j`, and keep all cross blocks.  If `Phi_D` and `Phi_G` denote the two
single-witness cross responses after subtracting their respective diagonal
norms, then

```math
\boxed{
\Phi_G-\Phi_D
=\sum_j\bigl(Q(D_j)-q_{|D_j|}\bigr)
 +\bigl(Q(\widetilde A)-q_n\bigr).
}
\tag{R12.18}
```

The first term is the full internal optimality excess and the second is a
nonnegative replacement excess.  Thus no compatible choice can make the raw
hybrid response close to the original response when the blocks carry
substantial excess: moving that excess into the hybrid cross response is
exactly what global minimality forces.

The weakest surviving replacement target is therefore global rather than
layerwise: choose all block minimizers jointly so that
`Q(tilde A)-q_n=o(n^{3/2})`, then use the resulting hybrid itself as a
near-minimizer with improved diagonal structure.  Alternatively compare
only the signed mismatch on the one common maximizing mosaic state and
retain/subtract the compulsory total block excess in (R12.18).  Neither goal
follows from orbit averaging or local lexicographic choice, and neither has
been proved here.  Any summable geometric-scale statement must be phrased
for this centered global replacement excess (or coupled directly to the
same-window endpoint allocation), not for the uniform errors in (10.506).

