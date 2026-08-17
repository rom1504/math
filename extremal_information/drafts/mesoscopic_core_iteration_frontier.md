# Iterating mesoscopic agreement cores: a conditional theorem and the exact stall

Date: 2026-08-17.

Status: **proved conditional theorem and decisive scope memo; no canonical edit**.

MB.2 is frozen.  This note asks whether its third witness can be iterated to
a growing packing at scale `cM_n`.  There is a clean multi-anchor iteration
theorem, with no accumulation of shell width.  Its exact state variable is
the intersection of the anchors' common-correct coordinates.  PP.4 supplies
that variable for the first two anchors but gives no lower bound after the
third.  A finite union of `o(E)` projective balls does not repair this gap,
because `M_n=o(E)`.

Thus repeated pairwise bipartite cores do not presently imply a growing
packing.  The next missing statement is a genuinely higher-order
common-correct-mass or non-recycling theorem, not another shell-cardinality
estimate.

## 1. A multi-anchor localized-flip theorem

Let `A` be an exact order-`n` minimizer, put

```math
E={n\choose2},\qquad M=Q(A)=M_n,
```

and let `z_1,...,z_k` be positive augmented cuts with

```math
\langle a,z_i\rangle\ge M-2s.                     \tag{MI.1}
```

Their common-correct reservoir is

```math
R(z_1,...,z_k)
=\{e:a_e(z_i)_e=+1\text{ for every }i\},
\qquad p_k=|R(z_1,...,z_k)|.                       \tag{MI.2}
```

### Theorem MI.1 (one reservoir adds one jointly separated state)

Fix `0<theta<1/2`.  If `1<=r<=p_k`, `2r<M`, `s+r<M`, and

```math
2^n\exp\{-2(1/2-\theta)^2r\}<1,                   \tag{MI.3}
```

then there is a positive augmented cut `z_(k+1)` of deficit `d<=2r` such
that, simultaneously for every `i<=k`,

```math
\boxed{
d_{\rm P}(z_{k+1},z_i)
\ge\min\{\theta p_k,M-s-r\}.}                     \tag{MI.4}
```

All `k+1` words lie in the common deficit-`2max{s,r}` shell.  The new word
is obtained by flipping actual edges of `A`; no auxiliary bridge is used.

#### Proof

Choose `F` uniformly from the `r`-subsets of the common reservoir.  For a
fixed augmented cut `w` with

```math
|\{e\in R:a_ew_e=-1\}|\le\theta p_k,
```

the hypergeometric Hoeffding bound makes the probability that at least
`r/2` sampled edges have `a_ew_e=-1` at most the right side of (MI.3)
without its factor `2^n`.  Union over all at most `2^n` augmented cuts.
There is therefore an `F` for which every cut having at least `r/2`
negative correlations on `F` has more than `theta p_k` negative
correlations on the whole reservoir.

Flip `F`.  Exact order-`n` minimality supplies an augmented cut `z` with

```math
M\le\langle a^F,z\rangle
=M-d-2\sum_{e\in F}a_ez_e.                        \tag{MI.5}
```

As in MB.1, at least `r/2+d/4` of the sampled correlations are negative,
and `d<=2r`.  Hence more than `theta p_k` reservoir coordinates have
`a_ez_e=-1`.  Every old word has `a_e(z_i)_e=+1` there, so its actual signed
distance from `z` exceeds `theta p_k`.

If that actual distance is `h_i`, positivity gives the complementary
bound

```math
E-h_i\ge {\langle a,z\rangle+\langle a,z_i\rangle\over2}
\ge M-s-d/2\ge M-s-r.                              \tag{MI.6}
```

Taking `min{h_i,E-h_i}` proves (MI.4). `square`

### Corollary MI.2 (packing-or-higher-order-collapse iteration)

Assume `M>=c_0n^(3/2)` and `s=o(M)`.  Fix `0<beta<=1` and an integer target
`K`.  Take `theta=1/4`, `r=ceil(9n log2)`, and
`s_*=max{s,r}`.  Starting from a positive deficit-`2s_*` list whose mutual
projective distances are at least `(beta/4-o(1))M`, repeat Theorem MI.1
(with `s_*` in place of `s`) whenever its common-correct reservoir has size
at least `beta M`.

Then, without any cumulative shell loss, one of the following happens:

1. the deficit-`2max{s,r}` shell contains `K` words at mutual projective
   distance at least

   ```math
   (\beta/4-o(1))M;                                \tag{MI.7}
   ```

2. at some `k<K`, the already constructed `k`-point packing has the
   explicit higher-order collapse certificate

   ```math
   |R(z_1,...,z_k)|<\beta M.                       \tag{MI.8}
   ```

Every new witness is obtained by a fresh perturbation of the same base
signing `A`, so the `O(n)` deficit in (MI.7) is paid once, not once per
iteration.

For the initial PP.4 pair, its common-correct set `Z` obeys
`|Z|>=M-2s`; hence, for every fixed `beta<1`, MI.2 recovers an
`Omega(M)`-separated third state.  (MB.2's sharper `(1/2-o(1))M` constant
uses a varying `theta` and the larger `O(n log^2n)` sample.)  Nothing in
PP.4 or MB.2 prevents that third state from making the triple intersection
in (MI.8) small.  The endpoint `beta=1` in the abstract iteration requires
its stronger reservoir premise and is not supplied by PP.4.

## 2. Why pairwise cores do not imply reservoir persistence

The missing implication is not a set-theoretic formality.  Even exact
augmented-cut words can have large pairwise common-correct sets and a tiny
joint one.

Take disjoint vertex sets `S,T` of size `k=Theta(sqrt n)` and the three
actual augmented cuts

```math
z_0=\mathbf1,qquad z_1=-c(v_S),qquad z_2=-c(v_T). \tag{MI.9}
```

Partition the edges into the four sign-pattern cells

```math
A_0=\delta(S)\cap\delta(T),
B_0=\delta(S)\setminus\delta(T),
C_0=\delta(T)\setminus\delta(S),
D_0=E\setminus(A_0\cup B_0\cup C_0).              \tag{MI.10}
```

Their `(z_0,z_1,z_2)` patterns are respectively

```text
(+,+,+), (+,+,-), (+,-,+), (+,-,-).
```

Choose an exact `+-1` signing whose cell sums, up to parity-one rounding, are

```math
\sum_{A_0}a=0,qquad
\sum_{B_0}a=\sum_{C_0}a=m,qquad
\sum_{D_0}a=-m,                                   \tag{MI.11}
```

where `m=Theta(n^(3/2))` is at most the two middle cell sizes.  Then

```math
\langle a,z_0\rangle
=\langle a,z_1\rangle
=\langle a,z_2\rangle=m+O(1).                    \tag{MI.12}
```

The pairs `(z_0,z_1)` and `(z_0,z_2)` each have
`Omega(n^(3/2))` common-correct coordinates in `B_0` and `C_0`, while all
three are simultaneously correct only on the `a=+1` part of `A_0`, of
size `O(k^2)=O(n)`.  Thus pairwise mesoscopic reservoirs do not satisfy a
Helly principle at the required scale, even for genuine cut words and an
exact sign coefficient vector.

This is a **cut-realizable pattern obstruction**, not an exact-minimizer
counterexample.  By choosing the signs randomly subject to the four cell
sums, a union bound over all `2^n` cuts gives `Q(A)=O(n^(3/2))`: every cut's
mean is at most the sum of the absolute cell sums, `O(n^(3/2))`, and its
centred fluctuation is `O(n^(3/2))` uniformly.  It does not prove that the
three displayed words are in an `o(n^(3/2))` shell of a true minimizer.
Its role is narrower and exact: cut algebra plus pairwise core mass alone
cannot prove persistence of (MI.2).

## 3. Why finitely many `o(E)` balls do not force growth at scale `M`

The scales in the proposed hypothesis are mismatched.  Since

```math
M=Theta(n^(3/2))=o(E),                             \tag{MI.13}
```

one projective cut ball of radius `o(E)` can contain a very large
`cM`-packing.  For example, reserve `L=Theta(sqrt n)` vertices and take a
constant-rate binary code of length `L` and fixed relative distance.  Every
mask lies in one projective cut ball of radius

```math
L(n-L)=Theta(n^{3/2})=o(E),
```

every pair has projective cut distance `Theta(n^(3/2))`, and the ball
contains `exp(Theta(sqrt n))` such points.  Therefore a finite `o(E)`-ball
cover neither bounds nor forces response entropy at the `M` scale.

If the assumed balls had radius `o(M)` instead, every newly produced word
at distance `cM` from all previous representatives would occupy a new ball.
But proving “from all previous representatives” is precisely MI.1's
common-reservoir hypothesis; pairwise PP.4 gives it only at the first step.
The ball cover does not create the missing labelled intersection.

## 4. Relation to FB.3 and AO.20

MI.1 is genuinely geometric but uses no new optimality engine: it is the
multi-anchor localized form of the exact finite-flip certificate underlying
FB.3.  MI.6 is exactly AO.20.  What is new is the higher-order state
`R(z_1,...,z_k)` and the no-accumulation iteration alternative MI.2.

Trying to bypass (MI.8) with shell cardinality returns to the known ceiling.
At projective edge radius `Theta(M)`, a cut ball has

```math
\exp(Theta(\sqrt n\log n))
```

words.  FB.3 at any `o(n^(3/2))` shell width supplies too little logarithmic
cardinality to beat that volume uniformly.  Thus this route should stop if
its next input is merely “the shell is large.”

## 5. Frontier verdict

There is a rigorous iteration theorem, but it is conditional on a genuinely
higher-order invariant.  The assumptions “finite `o(E)` ball cover” and
“every encountered near-antipodal pair has a PP.4 core” do not imply that
invariant and hence do not yet imply a growing `cM` packing.

The smallest non-cardinality next lemma is one of:

```math
|R(z_1,...,z_k)|\ge\beta M
\quad\text{for every fixed }k,                    \tag{MI.14}
```

after a suitable choice of the witnesses, or a replacement theorem proving
that failure of (MI.14) itself exposes a new projective state.  Either would
be genuinely beyond FB.3/AO.20.  Without such a higher-order linkage, MB.2
can recycle among finitely many `Theta(M)`-separated caps, and no further
iteration is justified.
