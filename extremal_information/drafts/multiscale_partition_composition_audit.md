# What the multiscale partition shell does and does not compose

Date: 2026-08-17.

Status: proof draft.  This note assumes MP.1--MP.7 from
`multiscale_partition_affine_shell.md`.  It isolates the strongest direct
cross-child consequence of the partition budget, gives a genuinely
multi-anchor exact-sign compiler, and proves a target-scale lower bound for
the whole linear-anchor/Lipschitz use of that compiler.  It does **not**
claim an impossibility for a theorem in which the omitted bridge channel is
cancelled jointly against the child and shore energies.

## 1. A common favorable support exists for a small family

Let `A^(1),...,A^(t)` be hollow signings on the same labelled vertex set
`[n]`, and put `Q_r=Q(A^(r))`.  Fix one **equitable** partition

```math
[n]=J_1\sqcup\cdots\sqcup J_q,                  \tag{MC.1}
```

whose cell sizes differ by at most one.

For each `r`, independently orient and switch an absolute ground state as
in MP.3.  Denote the resulting matrix by `D^(r)`, its nonnegative ground
row sums by `ell_i^(r)`, and set

```math
F_(r,a)=2\sum_(i\in J_a)\ell_i^(r)
        +4Q_-(D^(r)[J_a]).                       \tag{MC.2}
```

### Theorem MC.1 (simultaneous partition-shell support)

There is one block `I=J_a`, common to all `t` children and of size at least
`floor(n/q)`, such that for every `r` and every `S subseteq I`,

```math
\boxed{
 \rho_rH_(A^(r))((x^(r))^S)
 \ge Q_r-{8\sum_(u=1)^tQ_u\over q}.}            \tag{MC.3}
```

In particular, if `Q_r<=C n^(3/2)` uniformly, then a common coordinate
support has vanishing normalized shell width whenever `t=o(q)`.

#### Proof

The proof of MP.3 gives, separately for every child,

```math
\sum_(a=1)^qF_(r,a)\le8Q_r.                     \tag{MC.4}
```

Sum (MC.4) over `r` and choose a block for which the summed cost is at most
its average.  Since every summand is nonnegative,

```math
F_(r,a)\le\sum_uF_(u,a)
          \le {8\sum_uQ_u\over q}               \tag{MC.5}
```

for every `r`.  The exact flip calculation in MP.3 now gives (MC.3).
`square`

This is a useful correction to the naive claim that the favorable block is
necessarily child-dependent.  **Coordinate support** can be synchronized
across any fixed, or more generally `o(q)`-sized, family.  What MC.1 does
not synchronize is the ground-state gauge `x^(r)`, the energy orientation,
or the response of child `r` to mixtures of fields designed for the other
children.

The `t=o(q)` range is a real boundary of the budget argument, rather than
just an artifact of choosing the minimum separately.

### Proposition MC.2 (a target-scale many-child selection obstruction)

There are hollow real matrices of cap `Theta(n^(3/2))` for which every
child separately has an exact favorable partition block, while no common
block gives a vanishing-width affine shell for all children.

More explicitly, let `n=r^4`, `q=r`, `m=r^3`, and partition `[n]` into
`q` blocks of size `m`.  For `a in [q]`, let `D^(a)` be the matrix which is
the all-positive hollow clique on `J_a` and zero elsewhere.  Then

```math
Q(D^(a))={m(m-1)\over2}=Theta(n^(3/2)).          \tag{MC.6}
```

Every block other than `J_a` gives a zero-defect affine cube for child `a`.
But for every proposed common block `J_b`, child `b` has a mask
`S subseteq J_b`, `|S|=floor(m/2)`, for which

```math
Q(D^(b))-|H_(D^(b))(1^S)|
   ={m(m-1)\over2}-O(m)=Theta(n^(3/2)).          \tag{MC.7}
```

#### Proof

Only the clique block contributes.  Its energy at a spin of magnetization
`z` is `(z^2-m)/2`; this proves both (MC.6) and (MC.7). `square`

MC.2 is a sharp falsifier for deriving an arbitrarily deep common selector
from the partition budget alone: here `t=q` and the guarantee in MC.1 is no
longer small.  The example is weighted (`0/1`), not an exact signing and not
a near-minimizer.  It therefore does not disprove a signing-specific
synchronization theorem.  It says precisely that such a theorem would need
new input beyond MP.1--MP.7.

## 2. More than one scalar anchor can be compiled jointly

The scalar microcanonical compiler is not intrinsically one-dimensional.
The clean extension is obtained by conditioning cell sums.

Partition the `s` shore coordinates into nonempty cells

```math
[s]=C_1\sqcup\cdots\sqcup C_d,
\qquad h_c=|C_c|,                                \tag{MC.8}
```

and let `P` be the orthogonal block-averaging projection: `P eta` is
constant on each `C_c`, with value

```math
\bar\eta_c={1\over h_c}\sum_(j\in C_c)\eta_j.   \tag{MC.9}
```

For every old row `i` and cell `c`, prescribe a feasible integer

```math
u_(i,c)\in\{-h_c,-h_c+2,\ldots,h_c\}.           \tag{MC.10}
```

### Theorem MC.3 (block-microcanonical joint compiler)

There is `B in {+-1}^{n times s}` with

```math
\sum_(j\in C_c)b_(i,j)=u_(i,c)                  \tag{MC.11}
```

for every `i,c`, and an absolute constant `C` such that

```math
\boxed{
 \max_(eta\in\{+-1\}^s)
 \|B eta-BP eta\|_1
 \le C\sqrt{n(s-d)(n+s)}.}                     \tag{MC.12}
```

Here

```math
(BP eta)_i=\sum_(c=1)^d u_(i,c)\bar\eta_c.      \tag{MC.13}
```

Thus all `2^d` block-constant Boolean endpoints are preserved exactly, and
every endpoint is reduced, up to (MC.12), to one point in an at-most
`d`-dimensional field channel.  This is a joint compiler, not `d`
separately paid scalar compilers.

#### Proof

Independently for every row and every cell, sample the row signs uniformly
from the Boolean slice of sum `u_(i,c)`.  The cell sums, hence `BP eta`, are
deterministic.  For fixed `eta`, put

```math
X_i=b_i^T(I-P)eta.                               \tag{MC.14}
```

Each cell contribution is a centred sample-without-replacement sum of the
numbers `eta_j-bar eta_c`.  Hoeffding comparison and Hoeffding's lemma give
a subgaussian variance proxy bounded by

```math
C_0\sum_(c:h_c>1)h_c
 \le2C_0\sum_c(h_c-1)=2C_0(s-d).                \tag{MC.15}
```

For every fixed `(z,eta) in {+-1}^n times {+-1}^s`, independence of the
rows makes `sum_i z_iX_i` subgaussian with proxy `C_1n(s-d)`.  A union bound
over the `2^(n+s)` pairs gives, with positive probability,

```math
\max_(z,eta)|\sum_i z_iX_i|
 \le C\sqrt{n(s-d)(n+s)}.                       \tag{MC.16}
```

Maximizing over `z` for fixed `eta` is exactly `sum_i|X_i|`, proving
(MC.12). `square`

For fixed `d`, (MC.12) has the same order as the one-anchor compiler.  To
make the residual smaller at fixed shore ratio, the retained dimension has
to occupy a nonvanishing fraction of the shore.  The following theorem
makes that necessity quantitative and does not depend on this random
construction.

## 3. A projection no-go at fixed ratio

### Theorem MC.4 (multi-anchor residual lower bound)

Let `B in {+-1}^{n times s}` and let `P` be any rank-`d` orthogonal
projection on `R^s`.  Then

```math
\boxed{
 \max_(eta\in\{+-1\}^s)\|B(I-P)eta\|_1
 \ge {ns-\|BP\|_F^2\over\sqrt{2s}}.}           \tag{MC.17}
```

Suppose in addition that the range of `P` has an orthonormal basis
`u_a=eta^(a)/sqrt s`, where the `eta^(a)` are Boolean, pairwise orthogonal
anchor endpoints, and put `g^(a)=B eta^(a)`.  If

```math
\|g^(a)\|_1\le L n
\qquad(1\le a\le d),                             \tag{MC.18}
```

then

```math
\boxed{
 \max_eta\|B(I-P)eta\|_1
 \ge n\sqrt{s/2}
       \left(1-{dL^2n\over s^2}\right).}        \tag{MC.19}
```

Consequently, if `s/n -> theta in (0,infinity)`, `L=O(1)`, and `d=o(n)`,
the omitted channel has size

```math
(1-o(1))\sqrt{theta/2}\,n^(3/2).                \tag{MC.20}
```

Any scheme that wants the right side of (MC.19) to be `o(n^(3/2))` must
retain

```math
d\ge(1-o(1)){s^2\over L^2n}=Omega(n).           \tag{MC.21}
```

Linear rank is not by itself an information lower bound.  For example, a
Hadamard bridge can have a full orthogonal anchor basis with extremely
simple algebraic presentation.  MC.21 says that a fixed-ratio escape cannot
be a *few-channel* extension of the scalar compiler; it leaves open a
linearly growing but algebraically closed channel system.

#### Proof

Write `r_i=b_i(I-P)`.  The sharp elementary `p=1` Khintchine lower bound
gives

```math
E_eta|r_i^Teta|\ge{1\over\sqrt2}\|r_i\|_2.     \tag{MC.22}
```

Since `||r_i||_2<=sqrt s`,

```math
\sum_i\|r_i\|_2
 \ge {1\over\sqrt s}\sum_i\|r_i\|_2^2
 ={\|B(I-P)\|_F^2\over\sqrt s}.                \tag{MC.23}
```

Average over `eta`, then maximize, and use the orthogonal Frobenius
decomposition

```math
\|B(I-P)\|_F^2=ns-\|BP\|_F^2                  \tag{MC.24}
```

to get (MC.17).  Under the additional hypothesis,

```math
\|BP\|_F^2
 =\sum_(a=1)^d\|Bu_a\|_2^2
 \le {1\over s}\sum_(a=1)^d\|g^(a)\|_1^2
 \le {dL^2n^2\over s}.                         \tag{MC.25}
```

Substitution proves (MC.19)--(MC.21). `square`

MC.4 is an actual scalable no-go, but its scope matters.  It rules out the
following certificate architecture:

1. retain a sublinear-dimensional linear field channel;
2. require its orthogonal Boolean anchor fields to be balanced in `l_1`;
3. approximate every free shore endpoint by projection to that channel;
4. pay the omitted field by child trust-response `l_1` Lipschitzness.

At fixed ratio, this architecture incurs a leading `Theta(n^(3/2))`
certificate error.  MC.4 does **not** say that the actual parent cap differs
by this amount.  A theorem evaluating the retained and omitted channels
jointly against the child quadratic energy, before absolute values, could
escape it.  Such a theorem is exactly the joint-cancellation ingredient
that MP.1--MP.7 do not provide.

## 4. Why the positive pieces do not make a fixed-ratio recurrence

Combining MC.1 and MC.3 synchronizes more than was previously explicit:

* a common coordinate block can serve `t=o(q)` children;
* a finite-dimensional family of exact anchor fields can be compiled into
  one balanced exact-sign cross block;
* all block-constant shore words, not merely one selected word, are then
  represented jointly.

There are nevertheless two distinct unresolved arrows.

First, MC.1 gives child-dependent gauges `x^(r)`.  Feeding their anchor
fields into MC.3 makes a generic endpoint produce a *mixed* field

```math
\sum_(c=1)^d\bar\eta_c u_c.                     \tag{MC.26}
```

The MP response theorem for child `r` controls its own star-selector
language; it does not control (MC.26) when the other `u_c` came from other
children.  Bounding those channels separately throws away precisely the
joint cancellation under investigation.

Second, a fixed ratio `s=Theta(n)` makes the MC.3 residual
`Theta(n^(3/2))` unless `s-d=o(n)` at a quantitative rate.  MC.4 shows that
this is not a weakness of the probabilistic proof when only `d=o(n)`
balanced anchors are retained.  Near-order insertion (`s=o(n)`) makes the
same residual subleading, but that is the already-known scale and does not
accumulate into a fixed-ratio recurrence.

There is also no cross-**order** selector in MC.1.  Its conclusion concerns
children on one labelled vertex set.  For matrices of different orders,
one first has to specify an injection, restriction, or blow-up identifying
their coordinates.  Producing such a map with cap control is already the
missing transfer problem; the partition budget contains no rule for it.
MC.2 shows, at the level of the budget algebra, how a favorable coordinate
can cycle through a growing family even though every individual member has
many perfect choices.

The exact implication graph is therefore

```text
MP.1--MP.7
   |
   +--> common coordinate support for t=o(q) children       [MC.1]
   |       |
   |       X--> common gauge / mixed-field response
   |
   +--> growing one-shot affine selector language
           |
           +--> d-dimensional joint exact compilation       [MC.3]
                   |
                   X--> o(n^(3/2)) fixed-ratio error
                        for d=o(n), separately paid residual [MC.4]

Missing escape arrow:
joint child + bridge + shore cancellation before absolute values,
or a nonlinear/cross-level congruence not represented by a fixed projection.
```

Thus no fixed-ratio composition law with `o(n^(3/2))` loss follows from the
partition shell plus balanced linear compilation.  The strongest honest
positive statement is MC.1: finite-family **support** compatibility is
available.  The strongest honest negative statement is MC.4: turning that
support into a uniformly reusable fixed-ratio physical interface requires
linear anchor rank, unless one proves a genuinely joint cancellation
theorem.  Neither conclusion uses or assumes the full endpoint language.
