# Exact closed summaries versus one-shot future-response nets

**Status.** Adversarial theorem draft.  The quotient/congruence equivalence is
classical universal algebra.  The project-level consequence is the exact
distortion characterization below and the separation showing that a small
all-context response net need not be a small reusable feature algebra.

## 1. Decoder model

Let `(M,star,1)` be a commutative monoid and `F:M->R`.  An **exact closed
`epsilon`-summary** consists of

1. a finite monoid `(S,odot)`;
2. a homomorphism `sigma:M->S`; and
3. a terminal decoder `fhat:S->R`

such that

```math
|F(x)-\widehat f(\sigma(x))|\le\epsilon
\qquad(x\in M).                                   \tag{CSC.1}
```

Then summaries compose without returning to the objects:

```math
\sigma(x\star y)=\sigma(x)\odot\sigma(y),          \tag{CSC.2}
```

and every future answer is decoded from
`sigma(x) odot sigma(c)`.  Both associativity and the summary of a product
are exact; only the scalar answer is approximate.

This is stronger than a one-shot sketch from which a decoder, given the raw
context `c`, approximates `F(x star c)`.  It is also stronger than an
approximately multiplicative update whose error is reintroduced at each
composition.

## 2. Exact characterization

For a monoid congruence `theta` on `M`, let

```math
\operatorname{osc}_F(C)
=\sup_{x\in C}F(x)-\inf_{x\in C}F(x)              \tag{CSC.3}
```

for each congruence class `C`.

### Theorem CSC.1 (minimum closed-summary index)

The smallest number of states in an exact closed `epsilon`-summary is

```math
\boxed{
\min_\theta
\left\{|M/\theta|:
 \operatorname{osc}_F(C)\le2\epsilon
 \text{ for every }C\in M/\theta
\right\}.}                                       \tag{CSC.4}
```

The same statement holds with cardinal infima when no finite quotient
exists.

#### Proof

Given a summary, restrict `S` to the image of `sigma`.  Its kernel is a
monoid congruence.  If `x,y` lie in one kernel class, (CSC.1) and the triangle
inequality give `|F(x)-F(y)|<=2epsilon`, proving necessity.

Conversely, quotient by any congruence in (CSC.4).  On a class `C`, decode by
the midpoint of `inf_C F` and `sup_C F`; every value in the class is within
`epsilon`.  The quotient operation is well defined because `theta` is a
congruence. `square`

For a noncommutative monoid, the same proof uses a two-sided congruence.  If
the decoder is allowed the raw future context but the state has an exact
transition under every context, the kernel need only be a right congruence;
commutativity makes this distinction disappear.

### What the theorem does and does not say

- It is an exact law for **homomorphic feature algebras**.
- It does not claim that finding the minimizing congruence is easy.
- It does not characterize one-shot response metric entropy.
- It does not cover a product rule that re-rounds to a net and accumulates
  error, nor a nonassociative approximate product.
- For vector-valued outputs, diameter at most `2epsilon` need not imply a
  radius-`epsilon` center in the same normed space.  The scalar real midpoint
  is essential to the converse as stated.

Theorem HRC.1 gives one family of congruences in (CSC.4): Rees congruences
which collapse a response ideal.  Theorem CSC.1 shows that Rees collapse is
not automatically the globally smallest quotient; other congruence classes
may also have small `F`-oscillation.

## 3. Prime-cycle separation

Let

```math
M=\mathbb Z/p\mathbb Z,
\qquad
F_p(x)=\cos(2\pi x/p),                            \tag{CSC.5}
```

where `p` is an odd prime and composition is addition.

### Theorem CSC.2 (small response net, rigid closed algebra)

For every `0<epsilon<1`:

1. the complete future-response family

   ```math
   R_x(c)=F_p(x+c)
   ```

   has an `epsilon`-net in uniform norm of size at most
   `ceil(2pi/epsilon)+1`, independently of `p`;
2. for all sufficiently large primes `p`, the minimum exact closed
   `epsilon`-summary has `p` reachable states (equivalently, every reduced
   such summary has `p` states).

#### Proof

Put `k=ceil(2pi/epsilon)`.  If `p<=k`, choose every residue.  Otherwise take
the residue nearest each point of the regular `k`-mesh on the circle.
Every phase is within `pi/k` of a mesh point and that point is within `pi/p`
of its selected residue.  Because `p>k>=2pi/epsilon`, their sum is at most
`epsilon`.  Since cosine is one-Lipschitz in its argument,

```math
\sup_c|F_p(x+c)-F_p(y+c)|\le\epsilon,             \tag{CSC.6}
```

giving the claimed one-shot response net.

Congruences on a group are coset congruences of normal subgroups.  The prime
cyclic group has only the equality and universal congruences.  Equality has
`p` classes.  The oscillation on the universal class is

```math
1+\cos(\pi/p),                                    \tag{CSC.7}
```

because the closest residues to phase `pi` give the minimum
`-cos(pi/p)`.  For every fixed `epsilon<1`, (CSC.7) exceeds `2epsilon` once
`p` is large enough.  The universal congruence is then forbidden by
Theorem CSC.1, leaving only equality. `square`

The response maps form a smooth one-dimensional phase circle, so metric
covering sees only `O(1/epsilon)` states.  Exact algebraic closure sees a
simple group and permits no intermediate quotient at all.  This is a
scalable separation, not a finite residue artifact.

## 4. Consequence for the syndrome landmark theorem

The landmark summary in `phase3_syndrome_landmark_quotient.md` is a valid
small **one-shot all-context sketch**:

- its decoded profile is uniformly close to the true word profile;
- min-plus convolution by the raw future profile is nonexpansive; and
- one sketch answers every adversarial appended-fragment radius query.

What has not been proved is a binary operation on landmark summaries which
is associative and returns the designated summary of the exact support
union.  The state-dependent basis charts make this a substantive missing
step.  Calling the construction an “approximate congruence” is safe only in
the metric sense that response-close states remain response-close after a
common continuation.  It is not an exact closed feature algebra in the sense
of (CSC.1)--(CSC.2).

This does not weaken its positive answer to the stated one-shot
all-future-response net problem.  It does change the feature-algebra lesson:

```math
\text{low interface metric entropy + nonexpansive raw-context action}
```

is sufficient for one-shot future-query compression, while

```math
\text{a low-index small-oscillation congruence}
```

is necessary and sufficient for exact reusable algebraic compression.  The
prime-cycle example proves that neither condition implies the other with
comparable complexity.

## 5. Relation to the hard-core/Rees theorem

For a nonnegative antitone deficit with an absorbing terminal state,
Theorem HRC.2 collapses the maximal terminal interval `{F<=2epsilon}`.
That set is an ideal, so its Rees relation is one explicit admissible
congruence in (CSC.4).  The cyclic example has no order, terminal ideal, or
proper congruence; its smooth response geometry is invisible to exact
quotient algebra.

The combined law is therefore a genuine dichotomy:

1. **metric compression** controls a system used once against a raw future;
2. **congruence compression** controls a state repeatedly composed using
   summaries alone.

Any claim about feature-algebra growth must state which resource is meant.
Without this distinction, a small response net can be mistaken for a closed
algebra, or a congruence lower bound can be misreported as a one-shot
information lower bound.

## 6. Director judgment

Theorem CSC.1 is classical quotient logic, not a new semigroup theorem.  Its
extremal-information value is diagnostic: it gives the exact obstruction to
turning a response sketch into a reusable algebra.  Theorem CSC.2 shows the
obstruction can be maximal even for a one-dimensional smooth response
family.

The strongest next theorem is now more precise than “find a smaller net”:

> For syndrome-support union at distortion `epsilon*w`, either construct a
> low-index congruence with classwise radius oscillation at most
> `2epsilon*w`, or prove that every such congruence has much larger index
> than the landmark response net.

That question measures feature-algebra growth rather than one-shot source
coding and cannot be answered by the landmark restriction data alone.
