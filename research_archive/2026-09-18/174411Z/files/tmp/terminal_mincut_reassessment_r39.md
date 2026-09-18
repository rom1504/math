# Wave 39: terminal-min-cut reassessment

## Verdict

The Johnson-entry wall does not apply to the exact coordinate min-cut.  There
is an exact, entry-free subtree operation: on an active tree edge, flip every
disagreeing coordinate throughout the descendant subtree.  It deletes that
entire coordinate separator and changes child energy only through a laminar
sum of signed shores.  This gives a rigorous Euler inequality for any
flip-closed penalized selection problem.

The inequality does **not** yet prove the required likely cut compression.
The favorable fibers are not closed under the subtree flips, and the laminar
signed-shore terms neither telescope nor have a known project-scale upper
budget.  A complete-signing construction below shows that exact child energy,
zero shore, project-scale `Q`, and low row do not control conflicts between
arbitrarily prescribed grounds.  It is not an actual-minimizer counterexample
and therefore leaves an existential minimizer-specific selection theorem
open.

## 1. Exact descendant flip

Fix a rooted tree `T`, coherent oriented representatives `y^v` on `S_v`, and
optimal full Hamming completions `x^v` in (10.922).  For an edge
`e=(p,c)`, let `H_e` be the component below `c` and put

```math
U_e=\{i:x_i^p\ne x_i^c\}.
\tag{R39.1}
```

For every `v` in `H_e`, flip the coordinates in `U_e` in `x^v`, and flip
`U_e\cap S_v` in its partial label.  Every edge internal to `H_e` has both
ends flipped, every edge outside has neither end flipped, and `e` loses all
of its disagreements.  Therefore, writing `y^{(e)}` for the new partial-label
configuration,

```math
D_T(y^{(e)})\le D_T(y)-|U_e|.
\tag{R39.2}
```

This operation pays no `|S_v\setminus S_u|` term.  It acts on coordinates
already carried by the optimal completion and hence directly exploits the
free-coordinate screening absent from (10.953).

There is also an exact sparse-code description.  Given the root `z`, the
sets `U_e` reconstruct every completion by

```math
x_i^v=z_i(-1)^{\#\{e\in[0,v]_T:i\in U_e\}},
\qquad
\sum_e|U_e|=D_T.
\tag{R39.3}
```

Thus a forest of cost at most `k_0` is precisely a root plus at most `k_0`
coordinate--edge toggles, followed by restriction to the selectors.  For a
fixed tree with `s` non-root vertices the number of toggle sets of size at
most `k_0` is at most

```math
\sum_{j\le k_0}\binom{ns}{j}
\le \exp\!\left\{O\!\left(k_0\log\frac{ens}{k_0}\right)\right\}
=\exp\{O(L_0)\}
\tag{R39.4}
```

at the project scales.  This is exponentially cheaper than forcing the
selectors themselves into a Johnson tree.

## 2. Exact laminar Euler inequality

Attach to node `v` an energy sector `epsilon_v` in `{+1,-1}` and define

```math
\delta_v(y)=Q(A[S_v])-\epsilon_v y^T A[S_v]y.
\tag{R39.5}
```

For `F\subseteq S_v`, use the signed shore

```math
C_v(F;y)=\epsilon_v
\sum_{a\in F,\ b\in S_v\setminus F}A_{ab}y_ay_b.
\tag{R39.6}
```

Direct expansion gives

```math
\delta_v(y^F)-\delta_v(y)=4C_v(F;y).
\tag{R39.7}
```

Let `gamma>0`, and suppose `y` minimizes

```math
\Psi_\gamma(y)=D_T(y)+\gamma\sum_{v\ne0}\delta_v(y)
\tag{R39.8}
```

over a domain closed under every simultaneous descendant flip above.  Apply
(R39.2) and (R39.7) to one edge.  Minimality gives the exact Euler bound

```math
\boxed{
|U_e|\le4\gamma\sum_{v\in H_e}
C_v(U_e\cap S_v;y).
}
\tag{R39.9}
```

In particular the aggregate descendant shore sum for each edge is positive;
individual terminal shores may have either sign.  Summing the tree edges
yields

```math
\boxed{
D_T(y)\le4\gamma\,
\mathcal L_T(y,x),\qquad
\mathcal L_T=
\sum_e\sum_{v\in H_e}C_v(U_e\cap S_v;y).
}
\tag{R39.10}
```

The same statement holds for arbitrary additive soft potentials `pi_v`:
if `y` minimizes `D_T+sum_v pi_v(y^v)` over all labels, then

```math
|U_e|\le
\sum_{v\in H_e}
[\pi_v((y^v)^{U_e\cap S_v})-\pi_v(y^v)].
\tag{R39.11}
```

Equations (R39.9)--(R39.11) are the clean minimizer-specific interface.  They
use coherent orientations and the exact terminal separator, not a Johnson
surrogate.

## 3. Exact reduced lemma at the required exponent

Put `L_0=n^(3/4-c_0)`, `k_0=Theta(L_0/log n)`,
`r=Theta(n^(1/4+c_0)log n)`, and allow at most `T<=n^eta` groups, with
`eta<c_0`.  It is sufficient to prove the following against every selector
law, on an event of probability at least
`exp{-O(rL_0)}=exp{-O(n log n)}`:

For each group, there are a low-row root, a tree, additive soft potentials,
and a minimizer of their penalized completion problem such that

1. every selected partial label lies in its actual favorable fiber; and
2. the total positive laminar potential increment in (R39.11) is `O(k_0)`.

Then summing (R39.11) gives `D_T=O(k_0)` in every group, which is exactly the
forest input (10.923)--(10.924).  For the energy penalty (R39.8), item 2 is
the concrete estimate

```math
\mathcal L_T=O(k_0/\gamma).
\tag{R39.12}
```

This reduced lemma is genuinely weaker than Johnson concentration: (R39.4)
allows arbitrary selectors.  The unresolved content is simultaneous
**favorability plus laminar-curvature control**, not coordinate entry.

## 4. Why the present exchange identities do not close it

There are two exact gaps.

First, a favorable fiber is a hard, selector-dependent subset of projective
words.  The flip `U_e\cap S_v` can leave it.  Minimizing only over favorable
labels therefore does not permit the comparison in (R39.9).  Replacing the
hard constraint by a soft potential is useful only after proving that a
global penalized minimizer is still favorable.

Second, a terminal `v` occurs in every ancestor shore.  The sets `U_e` vary
with `e`, and all shores are evaluated at the same terminal word rather than
at the successive intermediate flips which would produce a telescoping
energy difference.  Existing exchange identities control signed first-order
totals; they do not upper-bound the unweighted, repeatedly reused laminar sum
in (R39.10).  Taking absolute values loses the needed scale.

## 5. Complete-signing scoped obstruction

The older wall used zero/real weights.  It can be upgraded to a complete
signing with the correct `n^(3/2)` scale, though not to a proved global
minimizer.

Partition `n=gb` vertices into `g` blocks of size `b`.  Take a symmetric
zero-diagonal signing `K` with `||K||=O(sqrt(g))`.  Let `b=2h`, choose a
symmetric sign matrix `R` with `||R||=O(sqrt(h))`, and put

```math
H=\begin{pmatrix}R&-R\\-R&R\end{pmatrix},
\qquad H\mathbf1=0,
\qquad \|H\|=O(\sqrt b).
\tag{R39.13}
```

Give each diagonal block of `A` the signing `J_b-I_b` and each off-diagonal
block `(a,c)` the signing `K_(ac)H`.  Then `A` is a symmetric complete
signing.  Decompose a spin in block `a` as `x_a=mu_a 1+v_a`, with
`v_a` perpendicular to `1`, and put `V=sum_a||v_a||^2`.  If

```math
kappa=\frac{\|K\|\|H\|}{b}\le1-\frac2b,
\tag{R39.14}
```

then, with `E_0=n(b-1)`,

```math
E_0-b(1+\kappa)V\le x^TAx
\le E_0-b(1-\kappa)V,
\qquad 0\le V\le n.
\tag{R39.15}
```

Consequently `Q(A)=E_0`.  Every block-uniform spin is an exact positive
ground and satisfies

```math
Ax=(b-1)x,qquad R_2(x)=n(b-1)^2=O(n^2).
\tag{R39.16}
```

Choose `b/g` a sufficiently large fixed constant.  Then (R39.14) holds and
`b=Theta(sqrt(n))`, so `Q(A)=Theta(n^(3/2))`.  A union of blocks can be
flipped from one block-uniform ground to another at projective distance
`Theta(n)`, while its signed shore is exactly zero because both endpoint
energies equal `E_0`.  Fixed-density selectors which are unions of blocks
inherit the same statement.  Two prescribed exact child grounds can
therefore force `D_T=Theta(n)>>k_0` while all relevant ground-to-ground shore
and deficit changes vanish.

This falsifies any deterministic claim that arbitrary exact grounds are cut
compressed merely from their deficits, shores, parent `Q`, and row scale.  It
does **not** falsify the desired existential theorem: the same fibers also
contain compatible common block grounds, and the construction is not known
to attain the global minimum `q_n`.  A successful theorem must use global
joint label selection and some additional consequence of exact global
minimality.

## 6. Verification

`tmp/terminal_mincut_reassessment_r39_check.py` verifies exactly:

- deletion of precisely one coordinate-edge disagreement set by every
  descendant flip on a random small tree;
- the deficit/shore factor `4` in (R39.7);
- reconstruction from the sparse cut code (R39.3); and
- an explicit order-16 instance of the complete-signing wall, with
  `Q(A)=112`, ground row-square `784`, zero inter-block shore, and two grounds
  at projective distance `8`.

It passes under the repository virtual environment.

## 7. Frontier recommendation

Retain terminal min-cut as an independent fallback, now in the sharper form
(R39.11)--(R39.12).  Do not spend another wave on Johnson transport or on a
bare signed-shore summation.  Revisit it only with one of two new inputs:

1. a soft favorable potential whose global minimizers remain favorable and
   whose laminar increments have an `O(k_0)` budget; or
2. an actual-minimizer theorem giving a canonical jointly compatible choice
   from the favorable fibers.

Without one of these, the exact subtree calculus stops at signed cancellation
and does not prove the required likely event.
