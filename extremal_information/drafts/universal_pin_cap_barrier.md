# Universal coordinate pins have quadratic cap

**Status.** Rigorous task-local draft.  This theorem closes the low-cap
version of the most direct exact metric compiler.  It applies to arbitrary
joint interaction among the auxiliary spins; its essential hypothesis is
that one fixed future makes one prescribed old configuration optimal for
**every** complete sign child.

It does not rule out a response-metric compiler whose optimizing old state
depends on the child while the resulting response differences nevertheless
encode the desired coordinate.

## Setup

For a hollow complete signing `A` on `k` old spins, write

```math
H_A(x)=\sum_(i<j)A_(ij)x_ix_j.
```

Fix an old--new interaction matrix `B in R^(k times m)` and any hollow
symmetric auxiliary matrix `C`.  The appended quadratic future induces the
even effective landscape

```math
g(x)=\max_(y in {+-1}^m)
 \{x^TBy+H_C(y)\}.                                  \tag{UP.1}
```

The completed parent for child `A` is

```math
P_A(x,y)=H_A(x)+x^TBy+H_C(y).                       \tag{UP.2}
```

We use the Boolean operator convention

```math
||B||_(infinity->1)=
 \max_(x\in{+-1}^k,y\in{+-1}^m)|x^TBy|.
```

No sign or norm assumption is needed for the first margin conclusion.  The
cap conclusion is most relevant when `B,C` are exact signs and (UP.2) is a
complete signing.

## Theorem UP.1 (universal pinning forces quadratic cap)

Fix `u in {+-1}^k` and `eta>=0`.  Suppose the same future `(B,C)` makes
`u` projectively `eta`-optimal for every hollow complete sign child:

```math
H_A(x)+g(x)
 <=H_A(u)+g(u)+eta                                  \tag{UP.3}
```

for every complete signing `A` and every `x in {+-1}^k`.  If `d` is the
projective Hamming distance from `x` to `u`, so `0<=d<=floor(k/2)`, then

```math
g(u)-g(x)>=2d(k-d)-eta.                             \tag{UP.4}
```

Consequently, with `d_*=floor(k/2)`,

```math
osc(g)>=2d_*(k-d_*)-eta,                            \tag{UP.5}
```

and

```math
||B||_(infinity->1)
 >=d_*(k-d_*)-eta/2.                               \tag{UP.6}
```

Every parent (UP.2), for every `A`, therefore has Boolean cap

```math
Q(P_A)>=d_*(k-d_*)-eta/2.                           \tag{UP.7}
```

In particular, an exact universal coordinate pin (`eta=0`) has

```math
Q(P_A)>= {k^2over4}-O(1).                          \tag{UP.8}
```

If its total order is `N=k+m=O(k)`, it cannot belong to a bounded
`Q(P_A)=O(N^(3/2))` family.

### Proof

Because both terms in (UP.1) are unchanged under `(x,y)->(-x,-y)`, `g` is
even.  Replace `x` by `-x` if necessary, and let `S` be the set of `d`
coordinates on which `x` and `u` differ.

Choose a child signing `A=A_(x,u)` by putting

```math
A_(ij)=-u_i u_j
```

on every edge crossing `(S,S^c)`; choose all remaining signs arbitrarily.
Only crossing edges change their Boolean monomial between `u` and `x`, and
each chosen edge gains two.  Hence

```math
H_A(x)-H_A(u)=2d(k-d).                              \tag{UP.9}
```

Substitute this `A` and `x` in (UP.3) to obtain (UP.4).  Choosing
`d=d_*` proves (UP.5).

For any two old states,

```math
g(x)-g(x')
 <=\max_y (x-x')^TBy
 <=2||B||_(infinity->1),                            \tag{UP.10}
```

which proves (UP.6).  Finally fix arbitrary `x,y`.  The two parent energies
at `(x,y)` and `(x,-y)` are

```math
H_A(x)+H_C(y)+x^TBy,
\qquad
H_A(x)+H_C(y)-x^TBy.                                \tag{UP.11}
```

At least one has absolute value at least `|x^TBy|`.  Maximizing gives

```math
Q(P_A)>=||B||_(infinity->1),                        \tag{UP.12}
```

and (UP.7)--(UP.8) follow. `square`

## Scope and consequence

The rank-one compiler in `algebraic_exact_sign_locking.md` satisfies the
hypothesis with `eta=0`, and its actual common baseline is `Theta(k^2)`.
UP.1 proves that this is not an accident of its ferromagnetic clique: every
fixed quadratic future that robustly pins one Boolean coordinate against
the entire complete-sign child class has quadratic cap, even with arbitrary
auxiliary--auxiliary interaction.

The hypothesis is deliberately stronger than exact contextual metric
embedding.  A max of many child-dependent affine pieces can agree with a
linear response on all coefficient-cube vertices without one common old
witness being active.  Excluding that possibility would require a new
exposed-face theorem; it is not silently assumed here.
