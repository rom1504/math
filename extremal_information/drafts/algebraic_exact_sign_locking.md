# Algebraic exact-sign locking: a metric compiler and two sharp lock ceilings

**Status.**  Rigorous task-local draft.  The exact finite checks are in
[`../experiments/verify_algebraic_exact_sign_locking.py`](../experiments/verify_algebraic_exact_sign_locking.py).

This note freezes the weakest useful exact-sign compiler before testing
particular lock architectures.  That distinction matters: a query-dependent
coordinate pin already transfers the complete response metric with only a
linear number of vertices, whereas a single query-independent Hamiltonian
cannot robustly identify a duplicate of every Boolean state.  The former is a
positive contextual theorem.  The latter obstruction explains why several
more ambitious equality/Hadamard-lock constructions fail.

There are three different targets.

1. A **metric compiler** only asks that every pairwise response distance be
   exposed by some predeclared disjoint exact-sign continuation.
2. A **pointwise overlay compiler** asks for one continuation indexed by `T`
   whose response approximates `max_x|H_B(x)-H_T(x)|` for every `B`.
3. A **low-cap compiler** additionally requires the complete parent itself to
   have Boolean cap `O(N^(3/2))`, without subtracting a quadratic baseline.

The first, weakest target is solved exactly below.  The second and third are
not.  In particular, subtracting a common calibration is legitimate for a
contextual response metric but not for the original signing minimization.

## 1. The weakest compiler and an exact coordinate pin

For a hollow symmetric signing `A` of order `k`, write

```math
H_A(x)={1\over2}x^TAx.
```

Fix a query `u in {+-1}^k`.  Append `k` new spins `y`, use the complete
rank-one sign bridge

```math
R_u=u\mathbf1^T,                                      \tag{EL.1}
```

and put a positive clique on the new block,

```math
C=J-I.
```

The resulting order-`2k` matrix is hollow and has a sign on every
off-diagonal entry.  Its one-sided parent response is

```math
F_u(A)=\max_{x,y\in\{\pm1\}^k}
 \{H_A(x)+x^TR_uy+H_C(y)\}.                           \tag{EL.2}
```

### Theorem EL.1 (exact disjoint coordinate compiler)

For every hollow symmetric signing `A` and every Boolean `u`,

```math
\boxed{
F_u(A)={3k^2-k\over2}+H_A(u).}                        \tag{EL.3}
```

Consequently, for any two children `A,A'`,

```math
\sup_u|F_u(A)-F_u(A')|
=\max_u|H_A(u)-H_{A'}(u)|
=Q(A-A').                                             \tag{EL.4}
```

Thus coordinate pins are an isometric embedding of the same-support
quadratic response metric into disjoint complete exact-sign contexts on twice
as many vertices.

#### Proof

Put

```math
a=x\cdot u,
\qquad b=y\cdot\mathbf1.
```

Then the locking part of (EL.2) is

```math
ab+{b^2-k\over2}.                                     \tag{EL.5}
```

Because `H_A(-x)=H_A(x)`, change the global sign of `x` and then of `y` so
that `a,b>=0`.  For fixed `a`, (EL.5) is a convex function of `b` and its
larger endpoint value is at `b=k`.  Hence an optimizer may be taken with
`y=mathbf1` and with `x` at Hamming distance `d<=k/2` from `u`, so
`a=k-2d`.  Relative to `(u,mathbf1)`, the lock loses exactly `2kd`.

Switching `d` coordinates of a quadratic signing changes only the
`d(k-d)` crossing edges.  Therefore

```math
H_A(x)-H_A(u)\le2d(k-d).                              \tag{EL.6}
```

The total change from the planted configuration is at most

```math
-2kd+2d(k-d)=-2d^2\le0.                              \tag{EL.7}
```

At `d=0`, (EL.5) equals
`k^2+(k^2-k)/2`, proving (EL.3).  Subtracting two copies of (EL.3) and taking
the supremum over `u` proves (EL.4). `square`

### Corollary EL.1a (exact-sign transfer of the short-seed Gram packing)

Apply EL.1 to the `2^h` children in Theorem SG.4/21.26.  At child order

```math
k=256r^2,
\qquad h={r(r-1)\over2}\ge{k\over1024},
```

every distinct pair has pin-response distance at least

```math
{\sqrt2\over16}k^{3/2}={1\over32}(2k)^{3/2}.          \tag{EL.8}
```

Hence at exact parent order `N=2k`, response error below
`N^(3/2)/64` requires at least

```math
h\ge {N\over2048}                                    \tag{EL.9}
```

bits.  The child block, bridge, and auxiliary block are all exact signs and
the future vertices are disjoint from the child.

This closes the **metric** exact-sign/disjoint scope cost of SG.4.  It does
not compile one negative-clone query `-H_T` into one context.  More
importantly, every parent in EL.1 has the common locking baseline

```math
{3k^2-k\over2}=\Theta(N^2).                           \tag{EL.10}
```

The theorem is therefore a response-information result, not a low-cap
construction near the original signing problem.

## 2. Classification of one-layer duplicate locks

The exact pin above is query dependent.  A tempting stronger construction
is one fixed complete sign Hamiltonian on paired blocks

```math
X=(x_1,\ldots,x_k),\qquad Y=(y_1,\ldots,y_k)
```

whose ground states include every duplicate `(u,u)`.  Exact signs force this
architecture to have a large hidden tie set.

### Theorem EL.2 (all duplicate ground states force one balance constraint)

Let `L` be a hollow symmetric signing on `X union Y`.  Every `(u,u)`,
`u in {+-1}^k`, is a global maximizer of `H_L` if and only if there is a sign
vector `s in {+-1}^k` such that

```math
L_(x_i,y_i)=1                                         \tag{EL.11}
```

and, for every `i<j`, the block between the two coordinate pairs is

```math
\begin{pmatrix}
L_(x_i,x_j)&L_(x_i,y_j)\\
L_(y_i,x_j)&L_(y_i,y_j)
\end{pmatrix}
=-s_is_j
\begin{pmatrix}1&-1\\-1&1\end{pmatrix}.              \tag{EL.12}
```

Equivalently, with

```math
d_i={x_i-y_i\over2}\in\{-1,0,1\},
```

the complete energy is

```math
\boxed{H_L(x,y)=k-2\left(\sum_i s_id_i\right)^2.}     \tag{EL.13}
```

In particular, for `k>=2` the duplicate relation is never isolated.  Every
duplicate ground state has a nonduplicate ground state obtained by
mismatching two coordinate pairs.

#### Proof

At a global maximizer, flipping one spin cannot increase the energy.  Flip
`y_i` at `(u,u)`.  Multiplying its local field by `u_i` gives

```math
L_(x_i,y_i)+
\sum_{j\ne i}
 \bigl(L_(y_i,x_j)+L_(y_i,y_j)\bigr)u_iu_j\ge0       \tag{EL.14}
```

for every `u`.  The products `u_i u_j` may be chosen independently.  The
first coefficient is `+-1`, while every parenthesized coefficient belongs to
`{-2,0,2}`.  Thus (EL.14) for all signs forces

```math
L_(x_i,y_i)=1,
\qquad L_(y_i,x_j)=-L_(y_i,y_j).                     \tag{EL.15}
```

Flipping `x_i`, and applying the same argument with `i,j` interchanged,
forces both row and column sums of every pair block to vanish.  Hence the
block has the form

```math
a_(ij)\begin{pmatrix}1&-1\\-1&1\end{pmatrix},
\qquad a_(ij)\in\{\pm1\}.                            \tag{EL.16}
```

Now mismatch a set `S` of coordinate pairs, leaving all other pairs
duplicated.  Relative to a duplicate, the energy change is

```math
-2|S|+4\sum_{i<j\in S}a_(ij)u_iu_j.                 \tag{EL.17}
```

For `|S|=3`, a triangle with
`a_(ij)a_(i\ell)a_(j\ell)=+1` admits signs making all three summands
positive.  Equation (EL.17) would then equal `6>0`, contradicting global
optimality.  Thus every triangle product is `-1`.  A complete signed graph
with this property is antibalanced: fix vertex `1`, take `s_1=-1` and
`s_i=a_(1i)`; then

```math
a_(ij)=-s_is_j.                                      \tag{EL.18}
```

(For `k=2`, choose either factorization directly.)  Substituting (EL.18) and
`x_i y_i=1-2d_i^2` into (EL.11)--(EL.12) gives (EL.13).

Conversely, (EL.13) is at most `k` and equals `k` on every duplicate, so the
displayed form is sufficient.  Finally choose distinct `i,j` with
`s_i d_i=-s_jd_j=1` and all other `d` zero.  This is a nonduplicate state
with the same energy `k`. `square`

### Consequence

Ordinary polarization by a universal complete sign equality lock cannot
implement the same-support overlay: the exact-sign constraint collapses the
putative `k` coordinate equalities to one signed magnetization equation.
This is a theorem about a **single fixed pair lock**.  It does not contradict
EL.1, whose rank-one bridge changes with the exposed query and pins one
chosen state rather than all `2^k` states simultaneously.

## 3. Why parallel balance locks need growing algebraic rank

One can try to repair EL.2 by superposing balance penalties with sign rows
`s^(1),...,s^(q)`:

```math
-\sum_{a=1}^q\left(s^{(a)}\cdot d\right)^2.           \tag{EL.19}
```

This is an idealized weighted architecture; it is more permissive than one
exact sign layer.

### Proposition EL.3 (ternary-kernel counting ceiling)

If (EL.19) vanishes only at `d=0` among
`d in {-1,0,1}^k`, then

```math
q\ge {k\log2\over\log(2k+1)}.                        \tag{EL.20}
```

#### Proof

Let `c_i in {+-1}^q` be the `i`th column of the row matrix.  If two distinct
subsets of columns had the same sum, their difference would give a nonzero
ternary vector in the kernel.  Thus all `2^k` subset sums are distinct.  Each
of their `q` coordinates is an integer in `[-k,k]`, so there are at most
`(2k+1)^q` possible sums.  Therefore `2^k<=(2k+1)^q`. `square`

A literal parallel implementation that spends one length-`k` duplicate
layer per balance square consequently uses
`Omega(k^2/log k)` spin occurrences.  EL.3 is not a lower bound for every
possible two-layer Ising gadget: a nonlinear shared auxiliary layer might
encode several constraints at once.  It does rule out the straightforward
Hadamard-code repair in which finitely many independent aggregate penalties
are expected to isolate all coordinate mismatches.

## 4. High-rate Boolean eigenspaces have only square-root lock margin

A different proposal uses a symmetric full sign matrix `W` with exponentially
many Boolean vectors in one extremal eigenspace (self-dual bent vectors are
the canonical example).  This gives many exact relations, but Frobenius mass
limits their robustness.

### Theorem EL.4 (universal repeated-eigenspace lock ceiling)

Let `W in {+-1}^{k times k}` be symmetric.  Suppose its `lambda`-eigenspace,
`lambda>0`, contains a Boolean code `mathcal C` with

```math
|\mathcal C|\ge2^{\alpha k}.                          \tag{EL.21}
```

Then

```math
\lambda\le\sqrt{k/\alpha}.                           \tag{EL.22}
```

Consider the following universal repeated-lock architecture with `s` blocks
of `k` spins.  Every cross block is `W`; all unexposed internal blocks use
the hollow part of `W`; one exposed block is replaced by an arbitrary hollow
sign child `A`.  If the architecture is required to keep every repeated
codeword `(u,...,u)`, `u in mathcal C`, even as a local maximum for every
exact hollow sign child with operator norm at most `9sqrt(k)`, then
necessarily

```math
s\ge1+{k-1\over\lambda}
\ge1+(k-1)\sqrt{\alpha/k}.                            \tag{EL.23}
```

In particular it needs `Omega(k^(3/2))` total vertices at positive code rate,
not `O(k)`.

#### Proof

If a real subspace has dimension `d`, some projection onto `d` coordinates
is injective on it.  Its intersection with the Boolean cube therefore has at
most `2^d` points.  Equation (EL.21) implies that the `lambda`-eigenspace has
dimension at least `alpha k`.  Frobenius accounting gives

```math
\alpha k\lambda^2
\le\|W\|_F^2=k^2,
```

which proves (EL.22).

Fix `u in mathcal C` and a coordinate `i`.  There exists a hollow symmetric
signing `A` such that

```math
u_i(Au)_i=-(k-1),
\qquad\|A\|_(2\to2)\le9\sqrt k.                      \tag{EL.24}
```

Indeed, prescribe the `i`th star by
`A_(ij)=-u_i u_j`.  Its operator norm is `sqrt(k-1)`.  Fill the remaining
principal block by a symmetric random signing.  A `1/4`-net and the standard
Rademacher quadratic-form tail show that some completion has norm at most
`8sqrt(k)`; the triangle inequality gives (EL.24).

At the repeated codeword, the local signed field of coordinate `i` in the
exposed block is

```math
(s-1)\lambda+u_i(Au)_i=(s-1)\lambda-(k-1).           \tag{EL.25}
```

It must be nonnegative at a local maximum.  This proves (EL.23). `square`

The theorem is deliberately scoped.  It kills a **universal** compiler that
relies only on an exponential Boolean eigenspace and a constant number of
repeated Hadamard locks.  It does not show that the particular alternating-
form children contain the adversarial star completion (EL.24).  To use a
Boolean eigenspace for that narrower family one must prove additional
state-dependent local-margin or restricted-witness structure.  The operator
bound `||A_B||=O(sqrt(k))` alone is insufficient, because one row field may
still have order `k`.

## 5. The exact Hadamard pullback obligation

Let `W` be a symmetric Hadamard matrix with `W^2=kI`.  On Boolean pairs
satisfying

```math
y={1\over\sqrt k}Wx\in\{\pm1\}^k,                   \tag{EL.26}
```

an exact sign query `C` on the `y` block pulls back to

```math
H_C(y)={1\over2k}x^TWCWx.                            \tag{EL.27}
```

Thus a Hadamard lock can reproduce a target signing `A_T` on its exact
Boolean relation only if the conjugate

```math
{1\over k}WCW                                         \tag{EL.28}
```

has the required off-diagonal signs (diagonal values are only a Boolean
calibration).  This is a rigid bi-unimodular closure condition, not a
consequence of spectral flatness.  The short-seed construction proves
simultaneous bounds for `A_u odot chi_B`; it proves neither (EL.28) nor that
the contrast witnesses lie in one Boolean eigenspace.

Accordingly, the minimal positive lemma for this algebraic route is:

> Find a predeclared Boolean regularizer code `U` and exact sign pullbacks
> `C_T` such that every nonzero alternating-form contrast has a
> `Theta(k^(3/2))` witness in `U`, and prove stability of that witness set
> under the complete lock.

This is strictly less than full parent maximization, but EL.4 shows that it
cannot be replaced by the bare assertion that `U` is exponentially large.

## 6. `l_1` elimination and completion accounting

An isolated auxiliary spin gives the exact identity

```math
\max_{z=\pm1}z(x_i-T_(ij)x_j)
=|x_i-T_(ij)x_j|
=1-T_(ij)x_ix_j.                                     \tag{EL.29}
```

One auxiliary per edge therefore compiles `-H_T` exactly, but uses
`Theta(k^2)` vertices before completing all missing edges.  With `m=O(k)`
dense auxiliaries and no internal interaction, elimination instead gives

```math
\max_y x^TRy=\sum_{a=1}^m|(R^Tx)_a|.                \tag{EL.30}
```

No representation of all SG.4 contrasts by such a positive `l_1` sum is
proved here.

It would be incorrect to dismiss exact completion merely by paying the
absolute cap of the auxiliary block separately.  A negative clique has

```math
H_{-(J-I)}(y)={m-(\sum_a y_a)^2\over2}\le {m\over2}, \tag{EL.31}
```

so its **one-sided** contribution is only linear.  Combined with a complete
sign bridge it yields the legitimate joint balanced-`l_1` response

```math
\max_y\left\{x^TRy+{m-(\sum_a y_a)^2\over2}\right\}. \tag{EL.32}
```

Whether a query-dependent dense bridge can make (EL.32) reproduce the
alternating-form pointwise overlay remains open.  It must be analyzed
jointly; independently bounding the bridge and completion would lose the
same leading cancellation the program is seeking.

## 7. Research judgment

The exact-sign question splits cleanly.

* **Solved:** contextual response distance and its linear information rate
  survive disjoint exact-sign composition, by EL.1 and Corollary EL.1a.
* **Falsified:** one universal two-block equality lock; every such lock is
  only a signed-magnetization constraint and has exact two-defect ties.
* **Falsified under a precise architecture:** a constant-replica extremal-
  eigenspace lock robust to all bounded-operator exact sign children.
* **Open:** a pointwise overlay compiler, a restricted alternating-family
  Hadamard pullback, or any exact-sign compiler whose parent cap remains
  `O(N^(3/2))` rather than hiding the response under a `Theta(N^2)` pinning
  baseline.

Thus algebraic locking does produce a real theorem, but it does not yet
provide a near-original composition law.  The next discriminating target is
not another generic Hadamard variant.  It is either a bounded-cap coordinate
exposure theorem with `o(N^(3/2))` calibration cost, or a proof that every
exact sign exposure of an exponential coordinate family necessarily carries
a quadratic common baseline.
