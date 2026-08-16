# Independent verification of the phase-3 algebra claims

**Status.** Independent proof audit.  I checked
`phase3_hard_core_rees_compression.md`,
`phase3_closed_summary_congruence.md`, and Propositions SN.3--SN.5 of
`phase3_selective_neutralization.md`.  I also exhaustively tested the matroid
response formula over all binary flats through width four, and the
multi-carrier pointwise inequality over all binary instances with
`w<=3`, `m<=3`, and every future dictionary.

## Verdicts

### 1. Maximal midpoint terminal cell: verified

For a nonnegative antitone deficit with absorbing terminal state,

```math
d_F(q,0)=F(q).
```

The identity context attains `F(q)`, and every future value lies in
`[0,F(q)]`.  Thus all states with `F(q)<=2eta` can share the constant decoder
`eta`, uniformly after every continuation.  Conversely, if `q` and `0`
share one deterministic uniformly `eta`-accurate message, their response
vectors are at distance at most `2eta`, hence `F(q)<=2eta`.  Antitonicity
makes this entire sublevel set an ideal, so its Rees collapse is closed under
arbitrarily many compositions.  The maximality claim and factor two are
correct.

This verdict assumes the summary model stated in the drafts: scalar real
output, deterministic uniform error, and either an exact homomorphic update
or the same raw context supplied to the decoder.  It would not automatically
extend to vector outputs in a norm without a radius/diameter hypothesis.

### 2. Matroid quotient, response metric, and rate: verified

For flats `X,Y`, submodularity gives

```math
|r(X\cup T)-r(Y\cup T)|
<=max{r(X\vee Y)-r(X),r(X\vee Y)-r(Y)}.
```

The contexts `Y` and `X` attain the two directed increments, respectively.
Consequently HRC.10 is exact, and equality of all future responses is
equivalent to equality of closure.  Exhaustive binary tests through `w=4`
found exact agreement for every pair of flats and all tested contexts,
including the two attaining witnesses.

The strengthened HRC.5 with fixed `epsilon<1/4` is also correct.  For
`w=2d`, a Grassmann packing at injection distance
`s=ceil(gamma*w)`, where `2epsilon<gamma<1/2`, has logarithmic size

```math
((1/2-gamma)^2-o(1))w^2.
```

Its response distance is at least `s>2epsilon*w`, so common messages are
impossible.  A basis description gives the matching `O(w^2)` upper bound.
For complete formal precision, the odd-width sentence can specify the
embedding `X=U direct-sum W` with a fixed line `U` and `W` ranging over a
middle-dimensional Grassmann packing in a fixed even-dimensional
hyperplane; the injection distances are unchanged.

### 3. Prime-cycle metric/algebra separation: verified

The regular phase mesh gives a one-shot response net with at most
`ceil(2pi/epsilon)` representatives (so the stated `+1` bound is safe): the
angular error is at most `pi/k+pi/p<epsilon`.  Exact closed summaries are
different.  The image of a group under a monoid homomorphism is a group, so
the kernel is a group congruence.  A prime cyclic group has only equality
and the universal congruence, while the universal class has oscillation

```math
1+cos(pi/p)>2epsilon
```

for every fixed `epsilon<1` and all sufficiently large `p`.  Hence exactly
`p` reachable states are necessary.  This genuinely separates a raw-context
response net from an exact reusable feature algebra; it does not claim a
lower bound for approximate/nonassociative updates.

### 4. Multi-carrier collapse under arbitrary futures: verified, with
minor domain corrections

The proof of SN.17 is valid uniformly for arbitrary background `B`, future
dictionary `E`, and target `x`.  All atoms in `P(L)` used by a shortest
representation can be combined to one vector `z in L`; decomposing
`z=z_1+...+z_m`, `z_i in W_i`, replaces that atom by at most `m` carrier
atoms.  Nothing about `B` or `E` is separately bounded, so cancellation with
future atoms does not invalidate the argument.  Maximization over `x`
preserves the same `m-1` bound.  Exhaustive tests verified the pointwise
claim in every binary case with `w<=3`, `m<=3`, and every `E`.

Three hypotheses should be made explicit before promotion:

1. SN.1/SN.7 requires `1<=d<=w`; for `d=0`, the displayed radius
   `w-d+1` is false.
2. SN.4 requires `m>=1` and should assume `D` (equivalently `bar D`) spans
   `G`, unless synthesis lengths are explicitly extended by infinity.
3. In SN.5, “a bounded number of finite labels” must mean labels from a
   fixed or otherwise `q^{O_m(w^2)}`-sized alphabet.  Unrestricted finite
   labels can encode arbitrarily much information and defeat the count.

With those corrections, SN.3--SN.5 are sound.  The statement for `m=o(w)`
correctly concerns the **dense-carrier part**: the common background `B` is
retained on both sides and may itself contain information.

## Model-distinction judgment

The coding and sparse-synthesis readings of SN.1 are operationally useful,
but they are not mathematically independent validations: both are exactly
the same finite-field projective-dictionary response algebra (and the binary
code case is literally the same word-length function).  They should be
described as two interpretations, not two distinct model classes.  The
matroid residual-rank theorem is the genuinely different third model.

## Minor terminology correction

After CSC.2, the landmark construction should consistently be called a
one-shot all-context **sketch** rather than a quotient or closed feature
algebra.  Its nonexpansive raw-context action is proved; an associative
homomorphic product on landmark summaries is not.

## Overall judgment

No requested central claim was falsified.  The strongest conclusions are:

1. response-net complexity and exact closed-algebra complexity can differ
   maximally even for a smooth one-dimensional response orbit;
2. the Rees midpoint cell is exactly maximal for scalar antitone deficits;
3. projective-matroid future-response complexity is genuinely
   `Theta(w^2)` for every fixed `epsilon<1/4`; and
4. bounded or `o(w)` many linear carriers collapse uniformly to their span
   before the extremum is taken, even under arbitrary adversarial futures.

The three explicit domain corrections above are required before theorem-file
promotion; none changes the asymptotic or conceptual conclusions.

## Addendum: bounded-composition rounded nets

I subsequently audited `phase3_bounded_depth_response.md`.

### 5. BDR.1 and the syndrome factor two: verified

For two child subtrees, commutativity and translation contraction give

```math
d_F(s_1\star s_2,x_1\star x_2)
<=d_F(s_1,x_1)+d_F(s_2,x_2).
```

Projection back to the net costs one further `delta`.  Counting one error at
each encoded leaf and each internal node gives `(2ell-1)delta`; if the leaves
are already centers, only the `ell-1` internal-node errors remain.  No tree
shape is being hidden in the estimate.

For syndrome landmarks, a summary cell has response diameter at most `2r`,
so selecting one actual support from each nonempty cell gives an actual net
of radius `delta=2r`, not `r`.  Therefore

```math
r=floor(epsilon*w/(2(2ell-1)))
```

does give final error at most `epsilon*w`.  The factor two in BDR.11 and the
entropy argument in BDR.12 are correct.  The projection at an internal node
may use unbounded computation on the union of its two representative
supports, exactly as stated.

The draft correctly refuses to call this an exact closed summary: the root
state need not be the designated encoding of the exact union, the rounded
operation need not be associative, and a binary evaluation tree is part of
the program.  CSC.2 therefore does not contradict the prime-cycle bounded-
composition upper bound.

Two scope clarifications are advisable:

1. The quantitative parameter is the number `ell` of leaves (equivalently,
   the total number of rounded compositions), not tree depth alone.  A tree
   of depth `D` can have `2^D` leaves, and the theorem gives no `O(D*delta)`
   guarantee.  “Bounded-composition” or “bounded-size” is more precise than
   “bounded-depth” unless this convention is stated.
2. BDR.12 should state `0<epsilon<1` (or at least
   `epsilon/(2(2ell-1))<1/2`) when invoking the binary-entropy covering
   asymptotic.  Larger relative errors are trivial or require a separate
   endpoint sentence.

No correction to a theorem constant is required.

## Addendum: tropical defect saturation

I also audited `phase3_tropical_defect_saturation.md` after it appeared.
The central identity TDS.8 is exact: all coarse factors combine to
`(star_i f_i) star b^(star m)`, nonexpansiveness gives the upper bound, and
`f_i=delta_0` attains it.  Since the kernel powers decrease to their
shortest-path closure, the arbitrary-depth defect saturates at
`||b-b_*||_infinity`.  Reapplying the blur at internal nodes only changes
the total kernel power, so it introduces no uncounted tree dependence.

The fixed-chart syndrome quotient TDS.4 is also correct.  Subgroup min-filtering
commutes exactly with min-plus convolution on `G/H`, and the common contained
basis bounds each coset's word-length oscillation by `r`.  Thus the true
radius lies in an interval of length `r`, giving midpoint error `r/2`; the
state count `(w+1)^(2^(w-r))` follows directly.  This is an exact algebra of
coarse profiles with an approximate scalar decoder, not an exact encoding
of the full word profile.

Three scope edits would make the statement fully precise:

1. TDS.2 should quantify over proper extended profiles (at least one finite
   value), or simply finite-valued profiles.  Convolution with a finite
   kernel does not make the identically-infinite profile finite.
2. The indicator-kernel corollary extends beyond the finite-kernel hypotheses
   of Section 1.  Say explicitly that the same algebraic proof is now used
   in the extended reals; the resulting defect may be infinite.
3. The asymptotic after TDS.19 should state `0<epsilon<1/2`, as required for
   `r=floor(2epsilon*w)<w` and the displayed nontrivial exponent.

The Hamming-ball sentence correctly measures `f-P_(mB_r)f`; it should not be
read as the exact one-blur-versus-`m`-blur defect, which is a different
quantity (and is infinite on the complete extended profile class).
