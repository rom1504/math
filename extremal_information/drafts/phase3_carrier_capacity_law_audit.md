# Adversarial audit: carrier capacity law

**Verdict.**  CC.1--CC.6 are mathematically correct, including the constant
`p`, all packing/covering directions, the rate--distortion thresholds, and
the rank-metric constants.  The verifier passes independently.  One scope
qualification is needed: CC.4 is a bona fide rank-metric Cayley realization,
but its packing proof uses only an equilateral host inside that metric.  It
does not yet exhibit a capacity phenomenon intrinsically special to rank
geometry, so the phrase "not a disguised Hamming example" is stronger than
the proved distinction.

## 1. Constant and metric-entropy transfers

From `0<=pi_theta<=p`, write

```math
F_\theta=d_{C_\theta}+e_\theta,
\qquad 0\le e_\theta\le p.
```

The important point is that

```math
\|e_\theta-e_{\theta'}\|_\infty\le p,
```

not `2p`.  Combining this with the exact identity

```math
\|d_C-d_{C'}\|_\infty=d_H(C,C')
```

proves CC.1.  The constant is sharp: the same singleton carrier with
presentation costs zero and `p` has Hausdorff distance zero and response
distance `p`.  Alternatively, CC.6 realizes a carrier separation of size
`p=diam(X)` with identical responses.

Every inequality in CC.2 has the correct direction.  Explicitly,

```math
d_H>s+p\Longrightarrow d_{resp}>s,
```

```math
d_{resp}>s\Longrightarrow d_H>s-p,
```

which gives the two packing comparisons.  A Hausdorff `(s-p)`-net is a
response `s`-net, while a response `(s-p)`-net is a Hausdorff `s`-net,
which gives the two covering comparisons.  Duplicate indexed carriers cause
no failure because their response distance is at most `p`; the upper packing
comparison is used only for `s>p`.

For deterministic reconstruction to uniform error `epsilon`, profiles in
one summary cell are at distance at most `2epsilon`.  Hence the lower
threshold is `2epsilon+p` in carrier Hausdorff distance.  A carrier
`(epsilon-p)`-net supplies the upper code.  The displayed sandwich assumes
`epsilon>p`, which is sufficient for both sides; the lower bound alone would
only require `2epsilon>p`.  The direct `Delta`-packing threshold
`epsilon<(Delta-p)/2` and the Fano bound are correct.

Thus CC.1 is sharp as a uniform theorem under only the hypotheses
`0<=pi<=p`.  At scales much larger than `p`, no stronger universal packing
comparison can be inferred without additional structure on the costs or the
carrier family.  This supports the claim that it is the strongest law
currently justified **inside the bounded-presentation carrier model**.  It
should not be read as a universal characterization of responses that have no
such carrier representation.

### Query-mass extension CC.1a

The added `L^s(mu)` statement is correct.  The pointwise perturbation

```math
(F_\theta-F_{\theta'})-(d_C-d_{C'})
```

has absolute value at most `p`, so its `L^s(mu)` norm is at most `p` for a
probability law.  The reverse triangle inequality gives CC.6b with the same
constant.

At a Hausdorff witness `x_0`, put `g=d_C-d_{C'}`.  Both terms are
one-Lipschitz, so `g` is two-Lipschitz and

```math
|g(x)|\ge\Delta-2t\qquad(x\in B(x_0,t)).
```

Subtracting the presentation perturbation and integrating over this ball
gives exactly

```math
(\Delta-2t-p)_+\,\mu(B(x_0,t))^{1/s}.
```

The constant two cannot be improved in a general metric space.  On a line,
take singleton carriers at `0` and `Delta`; moving distance `t` from the
first toward the second changes the difference of distance functions from
`Delta` to `Delta-2t`.

The scope distinction is important and correctly stated.  Uniform response
separation can be witnessed by one endpoint, whereas fixed-`s` average-query
separation also needs a neighborhood with nonnegligible `mu`-mass.  Under a
uniform law on a Hamming or rank-metric ambient space, the relevant ball may
have exponentially small mass, making CC.6c exponentially weaker.  The ball
factor is a rigorous local exposure lower bound, not necessarily the exact
global `L^s` distance: other query regions can contribute additional mass.
For a restricted sup-query set, the `p` comparison remains valid relative to
the restricted distance-transform norm, but Hausdorff entropy transfers only
when the restricted set still contains suitable witnesses.

## 2. Rank-metric application

For `E=F_{q^D}`, every nonzero multiplication map `M_a` is invertible and
therefore has rank `D`.  If `U ne U'` and `a in U\U'`, then for every
`a' in U'`,

```math
d_{rk}(M_a,M_{a'})=rank(M_{a-a'})=D.
```

The carrier Hausdorff distance is consequently the full diameter `D`.
The coefficient-support presentation cost is at most `k`, so CC.1 gives
response separation at least `D-k`.

The Cayley interpretation is exact.  A rank-`r` residual is a sum of `r`
rank-one matrices and no fewer, while the chosen multiplication maps form a
basis of `C_U`, so an element with coefficient support `s` needs exactly `s`
scalar-closed shortcut letters.  Minimizing over the shortcut contribution
is precisely CC.15.

The Gaussian product satisfies

```math
{D\brack k}_q
=\prod_{i=0}^{k-1}{q^D-q^i\over q^k-q^i}
\ge q^{k(D-k)}.
```

For `k<=D/4`, the response gap is at least `3D/4`; two reconstructions with
error `epsilon D` cannot coincide when `epsilon<3/8`.  Thus CC.18 has the
correct coefficient and strict threshold.

The supplied verifier independently exhausts all 29,161 pairs of weighted
carriers on the five-cycle, including sharp slack one for `p=1`.  In the
rank example it checks all 35 two-subspaces of `F_16`, their 595 pairs, and
all 65,536 matrix queries for one representative carrier.  It passes.

### Qualification on "genuinely non-Hamming"

The ambient landscape is genuinely rank-metric: its primitive moves are
rank-one updates and the query space is the full endomorphism algebra.  It is
therefore a valid non-Hamming *realization* of CC.1.  However, the lower-bound
proof sees the multiplication host `C_0` only through the fact that all of
its distinct points are at common distance `D`.  That restricted host is an
equilateral metric, and the same Grassmannian subset argument could be placed
inside a Hamming model with an equidistant host.  No rank-ball volume,
rank-intersection number, or other specifically rank-metric geometry enters
the capacity estimate.

Accordingly, CC.4 establishes portability beyond a Hamming ambient space,
but not yet an irreducibly rank-metric information law.  A stronger third
model would use carriers whose Hausdorff entropy depends on nontrivial rank
geometry rather than only an equilateral subspace.

## 3. Collapse examples

CC.5 is correct: all surjections have the same image carrier `W`, so with
zero presentation cost their distance response is identically zero despite
the large parametrization space.  This directly refutes raw algebraic
dimension as a capacity proxy.

CC.6 is also correct.  With `A=diam(X)`, every costly point satisfies

```math
d(x,c)+A\ge A\ge d(x,c_0),
```

while `c_0` attains `d(x,c_0)`.  Arbitrarily rich carriers containing
`c_0` therefore have the same response.  This both validates the need for a
submacroscopic presentation radius and realizes the full `p` loss in CC.1.
