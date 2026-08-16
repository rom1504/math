# Bounded-composition response compression and error accumulation

**Status.** Proved director draft.  The estimate is an elementary
quantization/error-propagation theorem.  Its role is to connect the one-shot
response net and exact-congruence regimes without pretending that a rounded
net is an associative feature algebra.

## 1. A rounded response program

Let `(M,star,1)` be a commutative monoid, let `F:M->R`, and put

```math
d_F(x,y)=\sup_{c\in M}|F(x\star c)-F(y\star c)|.  \tag{BDR.1}
```

This is a translation-contractive pseudometric.  Let `C subseteq M` be an
actual `delta`-net and fix a projection

```math
\pi:M\longrightarrow C,
\qquad d_F(x,\pi(x))\le\delta,                    \tag{BDR.2}
```

with `pi(c)=c` on `C`.  Define the rounded binary operation

```math
u\widehat\star v=\pi(u\star v),
\qquad u,v\in C.                                 \tag{BDR.3}
```

The operation in (BDR.3) need not be associative.  A binary evaluation tree
is therefore part of the program.

### Theorem BDR.1 (finite-composition law)

Label the `ell` leaves of any full binary tree by exact objects
`x_1,...,x_ell`.  Encode a leaf by `pi(x_i)` and evaluate every internal node
with (BDR.3).  If `s_T` is the root state, then

```math
d_F\left(s_T,x_1\star\cdots\star x_\ell\right)
\le(2\ell-1)\delta.                              \tag{BDR.4}
```

The bound is independent of the shape of the tree.  If the leaves already
belong to `C` and are not re-encoded, it improves to

```math
d_F\left(s_T,x_1\star\cdots\star x_\ell\right)
\le(\ell-1)\delta.                               \tag{BDR.5}
```

Consequently two bracketings of the same `ell` centers give states at
response distance at most `2(ell-1)delta`.

#### Proof

For a subtree `T`, let `e_T` be the response distance between its computed
state and the exact product of its leaves.  At an encoded leaf, `e_T<=delta`.
If `T` has children `T_1,T_2`, translation contraction and (BDR.2) give

```math
e_T\le e_{T_1}+e_{T_2}+\delta.                   \tag{BDR.6}
```

A full tree has `ell` leaves and `ell-1` internal nodes, proving (BDR.4).
Starting with zero leaf error proves (BDR.5).  Apply (BDR.5) to both
bracketings and use the triangle inequality for the last assertion.
`square`

### Corollary BDR.2 (depth-indexed state bound)

Let `N_F(delta)` be the least size of an actual `delta`-net of `(M,d_F)`.
Every product with at most `ell` leaves can be evaluated, for a declared
binary tree, to terminal response error `epsilon` with at most

```math
N_F\left({\epsilon\over2\ell-1}\right)           \tag{BDR.7}
```

states.

This is a bounded-composition statement, not an exact closed summary.  The
computed root need not equal the designated encoding of the exact product,
and two bracketings need not give the same state.  Equation (BDR.4) is exactly
the price of allowing that distinction.

## 2. Prime-cycle interpolation

For

```math
M=\mathbb Z/p\mathbb Z,
\qquad F_p(x)=\cos(2\pi x/p),                    \tag{BDR.8}
```

the phase-mesh proof of Theorem CSC.2 gives

```math
N_{F_p}(\delta)\le
\min\left\{p,\left\lceil{2\pi\over\delta}\right\rceil+1\right\}.
                                                            \tag{BDR.9}
```

Thus products of at most `ell` phases have a rounded response program with

```math
O\left(\min\{p,\ell/\epsilon\}\right)            \tag{BDR.10}
```

states at error `epsilon`.  In contrast, Theorem CSC.2 proves that for every
fixed `epsilon<1` and all sufficiently large prime `p`, an exact associative
homomorphic summary valid for arbitrary composition depth requires all `p`
states.

The same smooth response orbit therefore exhibits three regimes:

1. one raw future query: `O(1/epsilon)` metric states;
2. at most `ell` summarized factors: `O(ell/epsilon)` rounded states; and
3. arbitrary-depth exact closure: `p` congruence states.

This does not prove that (BDR.10) is the optimal bounded-depth rate.  It does
prove that low response-metric entropy alone controls only a declared finite
composition budget.

## 3. Consequence for syndrome landmark states

Corollary SL.2 gives actual syndrome response nets by choosing one support
from every nonempty landmark-summary cell.  Let `ell` be fixed and seek final
radius error `epsilon*w`, where `0<epsilon<1`.  Choose a landmark radius

```math
r=\left\lfloor{\epsilon w\over2(2\ell-1)}\right\rfloor.          \tag{BDR.11}
```

The resulting actual response net has radius at most `2r`, so Theorem BDR.1
turns it into a declared-tree rounded union program with final error at most
`epsilon*w`.  For fixed `epsilon,ell` its description length is bounded by

```math
2^{\left(1-h_2(\epsilon/(2(2\ell-1)))+o(1)\right)w}.             \tag{BDR.12}
```

The product is allowed unbounded computation: union the two representative
supports and project the result back to a selected representative.  It is
not the canonical landmark summary of the exact union and is generally
nonassociative.  The hard-core/Rees quotient supplies the complementary
arbitrary-depth result: a larger-error region is collapsed by a true
congruence, so no re-encoding error accumulates.

## 4. What this adds

Metric nets, error propagation, and congruences are classical ingredients.
The combined response law separates three information prices that the phrase
“composable summary” otherwise conflates:

- **contextual metric entropy** for one stored object against a raw future;
- **bounded-depth rounded entropy** at the finer scale forced by (BDR.4);
- **congruence entropy** for exact summary-only composition at arbitrary
  depth.

The prime cycle proves that the first and third can have unboundedly different
state counts even when the response orbit is smooth and one-dimensional.
Absorbing response ideals are a mechanism that defeats accumulation: their
Rees quotient is both metric-safe and an exact congruence.  Selector cubes
and Grassmann carriers are the opposite mechanism: legal contexts turn many
latent directions into a macroscopic response packing before the decoder
pays its error.

This is a general law about feature-algebra growth, but not a universal
compression theorem.  The arbitrary-table union-semilattice construction
in Theorem CRL.2 shows that the metric entropy in (BDR.7) can itself be as
large as an unconstrained response landscape.
