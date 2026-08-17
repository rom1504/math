# Independent audit: algebraic exact-sign locking

**Verdict:** PASS, with one wording repair recommended for EL.4 and important
response-scope qualifications.  I checked EL.1--EL.4 and the Hadamard
pullback calculation independently.  I also ran
`experiments/verify_algebraic_exact_sign_locking.py`; all 120,225 finite
checks pass.

## 1. EL.1 is an exact complete-sign disjoint compiler

With `R_u=u1^T`, the cross energy is

```math
x^TR_uy=(x\cdot u)(y\cdot1)=ab,
```

and the positive clique contributes `(b^2-k)/2`.  Both blocks may be
globally flipped without changing their internal quadratic energies, so an
optimizer may have `a,b>=0`.  For fixed `a`, the lock is convex in the
allowed scalar `b`; its larger endpoint is `b=k`.  Writing
`a=k-2d`, the lock loses `2kd` from its value at `(u,1)`.

Changing `d` old spins affects exactly `d(k-d)` old edges and changes each
edge contribution by at most two.  Hence the child can regain at most
`2d(k-d)`, leaving the exact deficit `2d^2`.  This proves

```math
F_u(A)=(3k^2-k)/2+H_A(u).
```

There is no factor-of-two error under `H_A=x^TAx/2`.  The order-`2k` parent
is hollow and complete: `A` and `J-I` supply exact signs within the two
blocks and every entry of `u1^T` is an exact sign across them.  Subtraction
and maximization over `u` give EL.4 exactly.

The scale conversion in EL.1a is also correct:

```math
{(\sqrt2/16)k^{3/2}\over(2k)^{3/2}}={1\over32},
```

so half-gap error is `N^(3/2)/64`, and
`h>=k/1024=N/2048`.

**Scope.**  `F_u` is the one-sided maximum of the complete parent, not its
absolute Boolean cap.  The bridge signs depend on the declared query `u`.
Thus EL.1 gives an exact contextual metric embedding with disjoint appended
vertices; it does not give one fixed `T`-indexed negative-clone overlay, and
it does not imply an identity for `Q` of the parent.  The common
`Theta(N^2)` positive calibration is harmless for response differences but
is not subtractable in the original minimization problem.  The draft states
these distinctions correctly.

## 2. EL.2 classification is iff

At a duplicate `(u,u)`, the one-spin local maximum condition for `y_i` is

```math
L_(x_i,y_i)+\sum_{j\ne i}
 (L_(y_i,x_j)+L_(y_i,y_j))u_iu_j\ge0
```

for every independent choice of the products `u_i u_j`.  The constant is a
sign and every other coefficient is in `{-2,0,2}`.  Therefore the constant
must be `+1` and every other coefficient zero.  Repeating this for `x_i`
and with `i,j` interchanged forces every `2 by 2` coordinate-pair block to
be

```math
a_(ij)[[1,-1],[-1,1]].
```

For a ternary mismatch vector `d`, its exact change from duplicate energy is

```math
-2|supp d|+4\sum_{i<j}a_(ij)d_id_j.
```

If a triangle product were positive, signs of three nonzero `d_i` could
make all three pair terms positive, increasing the energy by six.  Hence all
triangle products are negative.  This is equivalent to
`a_ij=-s_i s_j` (with the direct factorization used at `k=2`).  Substitution
gives exactly

```math
H_L=k-2(s\cdot d)^2.
```

This proves both necessity and sufficiency.  Choosing two nonzero mismatch
coordinates with opposite signed contributions gives the asserted
nonduplicate tie.  There is no hidden assumption that the duplicate maxima
are strict; indeed the theorem proves that they cannot be.

## 3. EL.3 counting is sharp at its stated level

If two different column subsets have the same sum, their signed incidence
difference is a nonzero vector in `{-1,0,1}^k` in the row kernel.  Thus the
`2^k` subset sums must be distinct.  Each of the `q` integer coordinates has
at most `2k+1` possible values, giving

```math
2^k<=(2k+1)^q
```

and EL.20.  The resulting `Omega(k^2/log k)` spin-occurrence cost applies to
a literal stack of `q` length-`k` balance layers.  It is not a lower bound
for a nonlinear shared auxiliary gadget, as the draft explicitly notes.

## 4. EL.4 eigenspace and star-row bounds

If a real eigenspace has dimension `d`, some projection onto `d` coordinate
positions is injective on it.  Its Boolean intersection therefore has at
most `2^d` points.  Positive-rate Boolean code size forces
`d>=alpha k`, and Frobenius accounting gives

```math
alpha k lambda^2<=||W||_F^2=k^2,
```

hence EL.22.

For a fixed codeword and row, prescribing
`A_ij=-u_i u_j` on the star gives signed local field `-(k-1)` and a star
matrix of norm `sqrt(k-1)`.  A symmetric random completion of the remaining
principal block has norm at most `8sqrt k` for some realization (a
`1/4`-net and the Rademacher quadratic-form tail suffice), so the full exact
hollow signing has norm below `9sqrt k`.  At a repeated codeword the exposed
coordinate's signed field is exactly

```math
(s-1)lambda-(k-1).
```

Nonnegativity is necessary even for a local maximum, yielding EL.23 and
`Omega(k^(3/2))` total vertices.

**Recommended wording repair.**  Replace “for every child with operator
norm `O(sqrt k)`” by the quantified statement “for every exact hollow sign
child with operator norm at most `9sqrt k`” (or with any fixed bound at least
that constant).  Big-O inside a universal hypothesis leaves the class
formally ambiguous.  The proof supplies and only needs the explicit
`9sqrt k` class.

The theorem assumes a universal repeated-block architecture and is not a
claim about the narrower alternating-form family.  It also uses a full
signed symmetric `W` (diagonal included) for cross blocks; internal copies
use its hollow part.  These choices still make the assembled parent hollow
and complete.  Only the exposed block's local field is used, so the missing
diagonal of internal copies causes no normalization error.

## 5. Hadamard pullback

For symmetric `W` with `W^2=kI` and a Boolean pair satisfying
`y=W x/sqrt k`, direct substitution gives

```math
H_C(y)={1\over2k}x^TWCWx.
```

Equality with an old quadratic up to Boolean calibration therefore requires
the off-diagonal part of `WCW/k` to equal the target signing; diagonal terms
are constant on the cube.  Spectral flatness alone supplies neither this
closure nor enough exact Boolean witnesses.  The stated regularizer/pullback
lemma is consequently a genuine remaining obligation rather than an
implicit consequence of the earlier Gram construction.

Subject to the one EL.4 quantifier clarification, the draft is rigorous and
its hierarchy is accurate: exact contextual metric compilation is achieved,
whereas pointwise negative-clone compilation and bounded-cap closure remain
open.
