# Audit of the proposed atom-net composition theorem

## Verdict

The proposed theorem is correct **provided one common net and one fixed
atom-to-centre quantizer are used throughout each target merge tree**. It is
a genuine extension of the *upper-bound/approximation* part of Theorem 16.15
from a finite atom alphabet to any atom-response image with finite metric
entropy. It does not subsume Theorem 16.15's arithmetic-rank exponent,
coarsest-state result, or robust lower bounds.

External centres and projective norms cause no problem when formulated in
the ambient normed vector space (or intrinsically in its quotient). The only
serious trap is allowing the net to vary locally with the current block
mass: then histograms need not be compatible and the claimed exact update is
false.

## Precise valid statement

Let `X` be a normed real vector space, let `Z` be an atom set, and let
`phi:Z->X`. Suppose `N_eta={v_1,...,v_D} subset X` is an **external**
`eta`-net of `phi(Z)`. Choose once and for all a deterministic map

```math
Q_eta:Z->[D],\qquad ||phi_z-v_(Q_eta(z))||_X<=eta.
```

For a word/multiset `A=(z_1,...,z_n)`, define

```math
F_A=sum_i phi_(z_i),
qquad c_j(A)=#{i:Q_eta(z_i)=j},
qquad Ftilde_A=sum_j c_j(A)v_j.
```

Then:

1. `c(A concatenated C)=c(A)+c(C)` exactly.
2. `||F_A-Ftilde_A||_X<=n eta`.
3. On any merge tree of total leaf mass `N`, if the same `Q_eta` is used at
   every leaf, the root error is at most `N eta`, independently of depth and
   bracketing.
4. The number of attainable mass-`n` summaries is at most

   ```math
   {n+D-1 choose D-1}.
   ```

   Equality is not asserted: some net centres can have empty preimage or
   different histograms can decode to the same response.
5. Hence

   ```math
   log_2 #states
   <=log_2 {n+D-1 choose D-1}
   <=D log_2(n+1).
   ```

   A useful symmetric alternative is

   ```math
   log_2 #states<=n log_2(D+1),
   ```

   so the minimum of the two bounds may be used.

The proof is just histogram addition, the triangle inequality, and stars
and bars. Importantly, the triangle inequality is paid once per microscopic
atom, not once per internal merge.

## Varying `D` and the root-scale-net requirement

For a sequence of target masses, one may choose a different pair
`(eta_N,N_(eta_N))` for each **whole target problem**. To compose a target
of total mass `N`, every leaf and intermediate summary must use that same
root-scale codebook. A child summary previously encoded using its own
mass-dependent net generally cannot be added to a sibling summary encoded
under another net.

A two-atom counterexample makes this unavoidable. Let `X=R` and
`phi(Z)={0,1}`. A coarse child codebook of radius one may map both atoms to
the centre zero, so its histogram loses the number of `1` atoms. A later
parent codebook of radius below `1/2` must separate `0` and `1`, but no update
from the coarse child histogram can recover the missing count. Thus the
family `eta_n` is not by itself a compatible multiscale state system.

Nested nets plus explicit refinement data could provide a different
theorem, but one would need either enough child information to refine or a
separate bound for every coarsening/requantization. Without that structure,
repeated requantization can accumulate error with depth.

## Asymptotic corollary

Let `D(eta)` denote an external covering number of the fixed atom-response
image, and make the root-scale convention above. If

```math
eta_n -> 0,
qquad D(eta_n)log(n+1)=o(n),
```

then the response error `n eta_n=o(n)` and the per-instance summary length
is `o(n)` bits.

If, for all sufficiently small `eta`,

```math
D(eta)<=C eta^(-p)
```

with constants `C<infinity` and `p>0`, then choosing
`eta_n=n^(-a)` for any fixed `0<a<1/p` gives

```math
n eta_n=n^(1-a)=o(n),
```

and

```math
D(eta_n)log(n+1)
<=C n^(ap)log(n+1)=o(n).
```

This calculation is correct. The codebook is treated as public side
information for the mass-`n` scheme; if the description cost or
computability of a nonconstructive net is part of the resource model, it
must be added separately.

## External-centre issue

External centres are legitimate here. The summary decoder returns an
element of `X`, not necessarily the response of a realizable atom or
landscape. Addition of centres is defined in `X`, and the error proof needs
only

```math
||sum_i(phi_(z_i)-v_(Q(z_i)))||
<=sum_i||phi_(z_i)-v_(Q(z_i))||.
```

This would fail only under an extra requirement that every intermediate
decoded state itself be a realizable microscopic object, or if future
continuations acted nonlinearly on atom identities rather than by addition
in `X`. Neither is part of the proposed theorem. Its scope should therefore
say explicitly that concatenation is represented by vector addition of
responses and that no mass-dependent normalization or cross-atom
interaction is present.

## Projective-response issue

The theorem works without changed constants in a genuine quotient normed
space such as

```math
X=l_infinity(Q)/R 1,
qquad ||[f]||=osc(f)/2.
```

Equivalence classes add linearly, an external net consists of classes, and
the same triangle inequality proves `n eta`. One must not independently
choose arbitrary representatives for atoms and centres and then interpret
their sum as a literal response: inconsistent additive baselines can grow
with `n`. Either work intrinsically in the quotient or choose one common
linear anchoring convention. A fixed-query anchor is generally only
bi-Lipschitz up to a factor two relative to half-oscillation, so literal and
projective radii should not be interchanged silently.

## Relation to Theorem 16.15

This proposal is strictly broader in one direction and strictly weaker in
another.

* **Broader atom class.** Theorem 16.15 assumes finitely many exact atoms.
  The net theorem permits infinite or continuous `Z`, provided
  `phi(Z)` is totally bounded at the selected scales. For instance a
  Lipschitz one-parameter curve of response atoms with
  `D(eta)=O(eta^(-1))` immediately has an `o(n)`-bit, `o(n)`-error additive
  carrier. This includes continuous heterogeneous-field quantization as a
  special case.
* **Weaker structural conclusion.** For finite atoms, the net histogram can
  be far from coarsest. Theorem 16.15 quotients exact additive relations and
  proves `Theta(n^(r_Z))` states, whereas the net theorem merely gives the
  crude `binom(n+D-1,D-1)` upper bound. Irrational precision, lattice margin,
  exposed conditioning, and lower bounds are absent from the new theorem.

At `eta=0` for a finite response image, the proposed construction recovers
an exact histogram upper bound, but not the arithmetic-rank theorem. It is
best presented as a **metric-entropy extension of the additive upper
carrier**, not as a replacement or strict strengthening of Theorem 16.15.

## Scope counterexamples to mention

1. **Changing normalization.** If a site's effective response depends on
   final mass, e.g. through a coefficient `1/n`, one fixed atom map `phi`
   does not describe all merge levels. The theorem does not apply without a
   common final-mass encoding.
2. **Interacting atoms.** A cross term depending on pairs of atom types is
   not determined by the summed atom response unless that interaction has a
   separately proved quotient closure.
3. **Nonlinear continuation.** An external centre is safe for additive
   readout, but a nonlinear future map may amplify its error or require
   information not present in the response vector.
4. **Mass-dependent local nets.** As above, exact histogram addition fails
   unless codebooks are shared or connected by certified compatible maps.

## Bottom line

Accept the theorem with the root-scale common-net hypothesis and intrinsic
projective formulation. The asymptotic choice `eta_n=n^(-a)`,
`0<a<1/p`, is correct. The result is a useful and genuinely more general
upper theorem for additive continuous atom families, but it should be
described as complementary to Theorem 16.15 rather than stronger than it.
