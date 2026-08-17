# Hadamard bridge synchronization: an exact deficit and a cap-scale defeat

**Status.**  Rigorous task-local draft.  The exact checks are in
`experiments/verify_hadamard_bridge_synchronization.py`.

An orthogonal sign bridge has the correct `k^(3/2)` cap scale and an exact
Boolean relation on bent inputs.  The question is whether that fact alone
supplies a bounded-cap exact compiler.  It does not.  The cross deficit is a
squared Euclidean distance to the normalized Hadamard image, and a
`Theta(sqrt(k))`-coordinate departure costs only `Theta(k)`.  An exact child
whose own Boolean cap is `O(k^(3/2))` can gain `Theta(k^(3/2))` along the
same departure.  This destroys the planted relation by a fixed leading
amount already at one layer.

The obstruction is deliberately sharper than saying that arbitrary internal
terms might be large: all three blocks are exact sign matrices, the complete
parent has `O(N^(3/2))` cap, and the loss is computed exactly.

## 1. Exact cross-deficit identities

Let `W in {+-1}^{k times k}` satisfy

```math
W^TW=kI.                                               \tag{HS.1}
```

For a Boolean old spin `x`, put

```math
v_x={1\over\sqrt k}W^Tx,
\qquad ||v_x||_2=\sqrt k.                             \tag{HS.2}
```

### Lemma HS.1 (pair and eliminated deficits)

For every Boolean pair `(x,y)`,

```math
k^{3/2}-x^TWy
={\sqrt k\over2}||y-v_x||_2^2.                       \tag{HS.3}
```

Consequently

```math
\begin{aligned}
\max_y x^TWy&=||W^Tx||_1,\\
delta_W(x)
&:=k^{3/2}-||W^Tx||_1\\
&={1\over2\sqrt k}\sum_{a=1}^k
       \left(|(W^Tx)_a|-\sqrt k\right)^2\\
&={\sqrt k\over2}\operatorname{dist}
       (v_x,{+-1}^k)^2.                              \tag{HS.4}
\end{aligned}
```

Thus `delta_W(x)=0` exactly when the normalized Hadamard transform is
Boolean.  On that regularizer code, the unique optimizing signs are
`y=v_x` unless a transform coordinate vanishes.

#### Proof

Both `y` and `v_x` have squared norm `k`, while
`y dot v_x=x^TWy/sqrt(k)`.  Expanding their squared distance proves
(HS.3).  Minimizing coordinatewise over Boolean `y` gives the first and last
lines of (HS.4).  Finally, `||W^Tx||_2^2=k^2`, so

```math
\sum_a(|(W^Tx)_a|-\sqrt k)^2
=2k^2-2\sqrt k||W^Tx||_1,
```

which is the middle identity. `square`

This is an exact Moreau-type formula, not a spectral estimate.  For any
future energy `K` on the new shore,

```math
\max_y\{x^TWy+K(y)\}
=k^{3/2}+\max_y\left\{
 K(y)-{\sqrt k\over2}||y-v_x||_2^2\right\}.         \tag{HS.5}
```

It makes the required regularity transparent.

### Proposition HS.2 (a quantitative sufficient condition)

Suppose `K` has an extension to the radius-`sqrt(k)` sphere which is
`L`-Lipschitz in Euclidean distance.  Then

```math
\max_y\{x^TWy+K(y)\}
\le k^{3/2}+K(v_x)+{L^2\over2\sqrt k}.               \tag{HS.6}
```

If `v_x` is Boolean, evaluating at `y=v_x` gives the matching zero-error
lower bound

```math
\max_y\{x^TWy+K(y)\}\ge k^{3/2}+K(v_x).
```
For `K(y)=H_C(y)` and a real symmetric `C`, one may take

```math
L=\sqrt k||C||_(2->2),
```

so the error in (HS.6) is at most

```math
{\sqrt k\over2}||C||_(2->2)^2.                      \tag{HS.7}
```

#### Proof

Put `r=||y-v_x||_2`.  Lipschitz continuity bounds the braces in (HS.5) by
`K(v_x)+Lr-sqrt(k)r^2/2`; its real maximum is
`K(v_x)+L^2/(2sqrt(k))`.  If `v_x` is Boolean, choose it.  For the quadratic,

```math
|H_C(y)-H_C(v)|
\le {1\over2}||C||\,||y-v||\,||y+v||
\le\sqrt k||C||\,||y-v||.
```

`square`

The proposition gives a genuine positive compiler criterion: a Boolean
pullback witness plus `||C||=o(sqrt(k))` has `o(k^(3/2))` synchronization
error.  It does not solve exact-sign closure, because every dense sign block
has Frobenius mass forcing `||C||>=sqrt(k-1)`.  A restricted exact-sign
theorem must exploit alignment beyond the generic Lipschitz estimate.

## 2. Exact departure from a Walsh regularizer

Take

```math
k=2^{2m},\qquad q=\sqrt k=2^m,                       \tag{HS.8}
```

and take the symmetric Sylvester--Walsh matrix, indexed by pairs in
`F_2^m times F_2^m`.  Let

```math
x_0(u,v)=(-1)^{u\cdot v}.                            \tag{HS.9}
```

The elementary Walsh summation gives

```math
Wx_0=q x_0.                                          \tag{HS.10}
```

Thus `(x_0,x_0)` is an exact Boolean Hadamard relation.  Let `S` be any set
of `d<q/2` old coordinates and let `x_S` flip `x_0` on `S`.

### Lemma HS.3 (small Hamming departures are exactly cheap)

The bridge optimizer for `x_S` is still `y=x_0`, and

```math
\begin{aligned}
x_S^TWx_0&=k^{3/2}-2d\sqrt k,\\
delta_W(x_S)&=2d\sqrt k.                             \tag{HS.11}
\end{aligned}
```

In particular `W x_S/sqrt(k)` is not Boolean when `d>0`.

#### Proof

For each transform coordinate `a`, put

```math
c_a=x_0(a)\sum_{i\in S}x_0(i)W_(ai).
```

Then

```math
{(Wx_S)_a\over\sqrt k}
=x_0(a)\left(1-{2c_a\over\sqrt k}\right).           \tag{HS.12}
```

Since `|c_a|<=d<sqrt(k)/2`, every factor in parentheses is positive, so the
coordinatewise optimizing sign remains `x_0(a)`.  Moreover

```math
\sum_a c_a
=\sum_{i\in S}x_0(i)(Wx_0)_i=d\sqrt k.             \tag{HS.13}
```

Summing the absolute transform coordinates proves (HS.11).  Its positive
deficit rules out Booleanity. `square`

The bridge charges only `2sqrt(k)` per flipped old coordinate throughout
this window.  A cap bound does not impose a comparable local charge on an
internal landscape.

## 3. A bounded-cap exact parent defeats the relation

We use one elementary probabilistic fact.  For every `n`, there is a hollow
complete signing `G_n` with

```math
Q(G_n)\le2n^{3/2}.                                   \tag{HS.14}
```

Indeed, for a fixed Boolean spin a random signing is a Rademacher sum over
`binom(n,2)` edges.  Hoeffding at `2n^(3/2)` and a union bound over `2^n`
spins have total probability below one (small orders may be checked or
absorbed into the constant).

### Theorem HS.4 (cap-scale Hadamard synchronization obstruction)

For all sufficiently large `k=2^{2m}`, there is a hollow complete signing
`A` on the old shore and a hollow complete signing `C` on the new shore such
that

```math
Q(A)\le3k^{3/2},
\qquad Q(C)={1\over2}k^{3/2}.                        \tag{HS.15}
```

The exact complete order-`N=2k` parent

```math
P=\begin{pmatrix}A&W\\W&C\end{pmatrix}              \tag{HS.16}
```

has

```math
Q(P)\le {9\over2}k^{3/2}=O(N^{3/2}).                \tag{HS.17}
```

Nevertheless, for `d=floor(sqrt(k)/4)` and some `S` of size `d`, the
off-relation configuration `(x_S,x_0)` beats the exact relation
configuration `(x_0,x_0)` by

```math
H_P(x_S,x_0)-H_P(x_0,x_0)
=2d(k-d-\sqrt k)
\ge {1\over8}k^{3/2}.                               \tag{HS.18}
```

Thus bounded cap, exact signs, and an exact Boolean Hadamard witness do not
give even one-layer robust synchronization.  Any depth-independent theorem
requires an additional alignment or pullback hypothesis.

#### Proof

On every edge crossing `(S,S^c)`, prescribe

```math
A_(ij)=-x_0(i)x_0(j).                                \tag{HS.19}
```

Fill the two principal shores independently by signings satisfying
(HS.14).  Flipping all spins in `S` preserves both principal energies and
changes every prescribed crossing contribution from `-1` to `+1`.  Hence

```math
H_A(x_S)-H_A(x_0)=2d(k-d).                           \tag{HS.20}
```

For every Boolean `z`, the absolute crossing contribution is at most
`d(k-d)`.  The two fillers therefore give

```math
Q(A)\le d(k-d)+2d^{3/2}+2(k-d)^{3/2}
\le3k^{3/2}                                          \tag{HS.21}
```

for all sufficiently large `k` and this choice of `d`.

Take `C` to be the hollow part of `W`.  The Walsh trace is zero at the even
tensor powers in (HS.8), so

```math
H_C(y)={1\over2}y^TWy.
```

The spectral upper bound and (HS.10) give the second equality in (HS.15).
The cross block obeys

```math
\max_{x,y}|x^TWy|
\le||W||\,||x||_2||y||_2=k^{3/2}.                   \tag{HS.22}
```

Triangle inequality proves (HS.17).

Compare the two displayed parent configurations.  They use the same new
spin `x_0`, so the `C` energy cancels.  Equations (HS.11) and (HS.20) give
the equality in (HS.18).  If `q=sqrt(k)` is large, then
`d>=q/8` and `k-d-q>=k/2`, proving its final bound. `square`

The conclusion is a stability obstruction, not the assertion that
`(x_S,x_0)` is the global parent maximizer.  That stronger assertion would
require controlling all exact relation states.  What is ruled out is the
universal compiler step “the Hadamard bridge keeps a planted Boolean
pullback witness valid under arbitrary `O(k^(3/2))` internal energies.”  It
fails by a fixed leading amount at depth one.

## 4. Research consequence

Hadamard synchronization has a clean exact state:

```math
\text{state}(x)=W^Tx/\sqrt k,
```

and a clean quadratic penalty to the Boolean shore.  But this state is not
closed under bounded-cap exact-sign interaction.  The failure is local in
Hamming geometry: bridge stiffness is `Theta(sqrt(k))` per coordinate,
whereas a cap-scale signing may coherently expose `Theta(k)` per coordinate
over a `Theta(sqrt(k))` cut.

Therefore a positive restricted compiler must prove at least one of:

1. a pullback identity making child/query gains cancel before the deficit is
   paid;
2. an operator/local-field bound, not merely a Boolean cap bound;
3. witness rigidity excluding the coherent cuts used in HS.4.

The generic Lipschitz theorem HS.2 quantifies the second option.  The exact
bi-unimodular conjugacy obligation in EL.27--EL.28 quantifies the first.  A
bare regularizer code, regardless of its size, supplies neither.
