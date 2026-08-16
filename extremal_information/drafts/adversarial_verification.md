# Adversarial verification of the extremal-information drafts

Date: 2026-08-16.

## Scope and verdict convention

This audit used only:

- `drafts/query_response_body.md`;
- `drafts/overlap_entropy_report.md`;
- `drafts/rate_distortion_report.md`; and
- the scripts and generated results in `experiments/`.

The long ledger and main state files were not read or modified.  `ACCEPT`
means that the stated claim and proof survive the audit.  `CORRECT` means the
core result survives but the indicated wording, convention, or boundary must
be repaired.  `REJECT` means the statement is false under its written
quantifiers; a replacement statement is supplied.

The main result of the audit is:

1. QR-A, QR-B, QR-C, QR-D, and QR-E are mathematically sound, subject to the
   convention and scope corrections below.
2. The zero-entropy support theorem, quadratic noise-cloud theorem, global
   response support identity, Curie--Weiss block pair, homometric code pair,
   and pinned-query information theorem survive exact reconstruction.
3. Proposition 4.4 of `overlap_entropy_report.md` is false as written for an
   arbitrary fixed positive bin width.  It becomes correct when the
   first-energy resolution is required to separate the two limiting edges.
4. The Ising response family is `2a`-separated, not a `2a`-packing under the
   drafts' stated strict-packing convention.  It is an `r`-packing for every
   `r<2a`.  The lossless lower bound for error strictly below `a` is still
   valid and sharp.

## 1. Query-response body

### QR-A: response duality — ACCEPT

For a finite landscape, extend the roof by `-infinity` off `P_phi`.  It is a
proper closed concave polyhedral function.  Its concave conjugate is exactly

\[
V_H(\theta)=\max_u\{\widehat H_\phi(u)+\langle\theta,u\rangle\},
\]

and concave biconjugacy gives QR.5, including at relative-boundary points of
`P_phi`.  Equivalently, finite-dimensional linear-programming duality gives
the same formula.  No missing relative-interior hypothesis is needed.

Required precision patch: when comparing landscapes with different feature
polytopes, say that their roofs are compared as extended functions on the
common ambient `R^d`, with value `-infinity` off their respective domains.
Then equality of all responses also recovers the domain.

The phrase "minimal quotient" is valid only in the equivalence-relation
sense: the roof is a complete invariant for equality of all declared linear
query answers.  It is not a minimal-bit representation or a compression
theorem.  This distinction should be stated beside that phrase.

### QR-B: sup-convolution — ACCEPT

The generator set of the composed lifted polytope is the set sum of the two
generator sets, and

\[
\operatorname{conv}(A+B)=\operatorname{conv}(A)+\operatorname{conv}(B).
\]

Taking the upper fiber at `u` yields QR.6, while support functions yield
QR.7.  The maximum in QR.6 exists because the feasible decomposition set is
compact.

Editorial patch only: write `K_{phi_oplus}(H_oplus)` rather than the ambiguous
`K_phi(H_oplus)` in the displayed equivalence.

### QR-C: one-step bilinear maximum — ACCEPT

There is a shorter adversarial proof that makes the independence assumption
transparent.  If distributions `lambda` and `mu` attain the two roofs at
`u,v`, then the right-hand objective in QR.9 equals

\[
\mathbb E_{x\sim\lambda,y\sim\mu}
 [H_1(x)+H_2(y)+\phi_1(x)^TB\phi_2(y)],
\]

which is at most the maximum over pure pairs.  Point masses give the reverse
inequality.  Thus no hidden joint-correlation assumption is being made for
the scalar one-step maximum.

The limitation stated after the theorem is real: QR-C does not produce the
parent roof for a new feature, and therefore is not an iterative closure
theorem.

### QR-D: packing/covering sandwich — ACCEPT, with two corrections to its use

The sandwich is correct for deterministic worst-case summaries if:

- `Pack(r)` uses pairwise distance strictly greater than `r`; and
- `Cov(epsilon)` is an internal cover by landscapes in `H` (an external cover
  in response-function space would also work if its decoder centers are
  explicitly allowed).

Two response functions decoded from the same state are at distance at most
`2 epsilon`, and an internal `epsilon`-net gives the upper bound.

The subsequent Ising sentence must be corrected.  The proof in
`rate_distortion_report.md` establishes

\[
d_\Theta(H_A,H_B)\ge 2a\quad(A\ne B),
\]

not strict inequality.  Under the draft's explicit convention this is not a
`2a`-packing.  Already at `n=2`, the two landscapes have query distance
exactly `2a`.  Replace "a `2a`-packing" by either:

- "a `2a`-separated family"; or
- "an `r`-packing of size `2^N` for every `r<2a`."

The lower bound for uniform error strictly below `a` is unaffected, either by
the direct Walsh decoder or by choosing an `r` strictly between twice the
actual error and `2a`.

The Hausdorff paragraph also needs scope precision.  QR.11 samples support
functions only in directions `(theta,1)`.  Its supremum is exactly a
restricted-support-function pseudometric, not the ordinary Hausdorff metric
unless the query experiment is enlarged to all normalized support
directions.  Positive last-coordinate directions recover the upper roof;
negative directions would additionally inspect the imposed lower body.

### QR-E: fixed-ambient compactness and finite realization — ACCEPT

The compactness proof survives the two likely failure points.

1. Downwardness is closed under Hausdorff limits.  If `(u_n,t_n)` tends to
   `(u,t)` and `s<=t`, use `(u_n,min(s,t_n))`; these points belong to the
   approximating downward sets and tend to `(u,s)`.
2. If `S_delta` is a finite `delta`-net in `U`, then
   `conv(S_delta) subset U` and is still `delta`-dense in `U`.  Its downward
   extension is also contained in `U`, contains `conv(S_delta)`, and hence
   remains `delta`-dense.  Treating the net points `(u,t)` as finite states
   makes this extension exactly the truncated response body of that finite
   landscape.  Arbitrary points below the target roof cause no problem: the
   upper hull and downward extension remain inside `U`.

QR.14 is the standard support-function bound, but the statement should name
the norm used to define `d_H`; `||.||_*` is then its dual norm.

Terminology patch: the theorem fixes the ambient feature dimension and region
`C`, not a fixed feature alphabet, state set, or fixed projection polytope.
"Fixed-ambient-interface compactness" would prevent a stronger
model-specific interpretation.  The draft correctly excludes constrained
Ising/sign-matrix realization.

### Compression and tautology audit

- Proposition 5.1 of the overlap report and QR-A's sufficiency direction are
  exact changes of coordinates: an objective that factors through retained
  data can be optimized from that data.  They are correct, but by themselves
  are not compression results.
- QR-D is an operational metric-entropy reduction.  It does not estimate the
  covering numbers; its structural content begins only when a class-specific
  packing or cover is proved.
- With `phi(x)=x` on the Boolean cube, QR.10 recovers `H(x)` at every vertex.
  The roof therefore retains the full labelled landscape.  The draft
  explicitly acknowledges this, so there is no concealed compression claim.

## 2. Entropy and overlap report

### Zero-entropy definitions — CORRECT

The separation between support-sensitive `log 0=-infinity` data and clipped
`log(1+N)` data is necessary and correctly described.  Theorem 4.1 uses the
former; Theorem 4.2 is what permits edge recovery from the strictly positive
part for homogeneous quadratics.

One definition should be aligned with the later theorem.  Section 3.2 defines
positive-tail identity using "limiting rates at continuity points", while
Corollary 4.3 defines the operative function with a `liminf`.  Choose one of
these formulations explicitly.  The strongest clean version for the stated
corollary is to compare the loci

\[
\{t:\liminf_n n^{-1}\log N^\uparrow_{H,n}(t)>0\}.
\]

### Theorem 4.1, support remembers the maximum — ACCEPT

If `m_H>m_G`, a rational interval `(r,s)` with
`m_G<r<m_H<s` is eventually nonempty for the first sequence and empty for the
second.  Its two capacities are respectively at least zero and `-infinity`.
The use of `limsup` causes no gap because convergence of the maxima makes
these eventual statements, not merely subsequential ones.

### Theorem 4.2, universal quadratic high-energy cloud — ACCEPT

For a uniformly chosen Hamming sphere of radius `r`, the exact multiplier is

\[
\lambda_{n,r}=\frac{(n-2r)^2-n}{n(n-1)}.
\]

The expectation identity follows term by term for every homogeneous
quadratic.  Since every sampled energy is at most `M_n`, the displayed lower
bound on the high-energy fraction follows whenever
`lambda_{n,r}>theta`.  The binomial-sphere exponent and limiting constant in
4.1 are correct.

### Corollary 4.3, positive tail determines the edge — ACCEPT

For every strict sub-edge threshold, Theorem 4.2 supplies an exponential
cloud; every super-edge threshold is eventually empty.  This proves the
endpoint identity.  Add one boundary sentence: if `max H_n=0` for some `n`,
then a homogeneous quadratic with maximum zero is identically zero (its cube
average is zero), so the required tail count is trivially `2^n`.  This closes
the only case not literally covered by Theorem 4.2's `M_n>0` hypothesis.

### Proposition 4.4, pair profile detects a rare high state — REJECT AS WRITTEN

The literal quantifier "fix positive energy and overlap bin widths" is too
strong.  A fixed energy grid can be coarser than the gap between the two
limiting maxima.  In the extreme, one energy bin contains the entire common
interval `[-C,C]`; the remaining overlap counts are then landscape-independent,
so different maxima need not be detected.

The proof silently chooses a first-energy bin lying wholly above the smaller
edge and containing the larger edge.  Such a bin need not exist in an
arbitrary pre-fixed grid.

Replacement proposition: if `m_H>m_G` and the pair profile is retained for
all rational open boxes, choose a rational interval `I` with

\[
m_G<\inf I<m_H<\sup I.
\]

Alternatively, require a grid/algebra containing such an `I` (a sufficiently
fine, suitably aligned grid or an allowed union of cells).  Pair an
`H`-maximizer with all `2^n` second replicas and pigeonhole over the finite
second-energy/overlap partition.  One fixed cell recurs on a subsequence with
exponent at least `log 2`, while the corresponding `G` cell is eventually
empty.  With this resolution hypothesis, the proof is correct.

### Proposition 5.1, global two-replica response — ACCEPT (TAUTOLOGICAL)

The response is the maximum of a function over the image support
`T_{H,n}`.  Exact support therefore determines it.  This is a useful boundary
on possible counterexamples, but it is deliberately a factorization
tautology, not an entropy-compression theorem.  The asymptotic Hausdorff
version is valid for continuous objectives on a common compact domain.

### Curie--Weiss left/right block pair — ACCEPT

All constants and feasibility conditions check exactly for even `m`:

\[
\max Q_L=m-1,\qquad \min Q_L=-1,
\]

\[
\Gamma^L_{H^L,n}(u)=\frac{mu^2-1}{m-1},\qquad
\Gamma^L_{H^R,n}(u)=1,
\]

and, under `R_L(x,y)=0`, independent sign flips allow `u,v>=0`, while the
joint-type constraint gives `u+v<=1`.  Hence the exact maximum of
`Q_L(x)+Q_L(y)` is `m-2`, attained by a uniform and a balanced block.  Thus

\[
\Theta^L_{H^L,n}(0)=\frac{m-2}{2(m-1)},\qquad
\Theta^L_{H^R,n}(0)=1.
\]

The soft perturbation constants also check:

\[
\mathcal B^L_{H^L,n}=\frac1{m-1},\qquad
\mathcal B^L_{H^R,n}=\frac{m+1}{m-1}.
\]

The coordinate swap proves exact equality of every global
energy--energy--overlap fiber.  The construction is valid only for a labelled,
anchored block experiment; modulo simultaneous relabelling of the landscape
and apparatus it is the same object.  The draft states this limitation
correctly.  It also correctly does not claim a dense hollow `+/-1` pair.

The direct-sum response-gap calculation is correct:

\[
\frac{a_n}{2(a_n+b_p)}\frac{m}{m-1}
=\frac{a_n}{2(a_n+b_p)}+o(1).
\]

Assume explicitly that `a_n+b_p>0` when introducing that normalization.

### Homometric code tensor pair — ACCEPT

For

\[
C=\{0000,0011,0101,0110\},\qquad
D=\{0000,0011,0101,1001\},
\]

both ordered pair-distance enumerators are `{0:4,2:12}`.  Their covering
radii are respectively `2` and `3`; `1110` is distance `3` from every word of
`D`.  For Cartesian concatenation, pair-distance enumerators convolve and

\[
r(C^{\times k})=k r(C),\qquad r(D^{\times k})=k r(D),
\]

because distance-to-code separates blockwise and maxima of the resulting sum
separate.  The normalized covering-radius gap is therefore exactly `1/4` at
every tensor power.

The experiment script computes the convolved enumerator but records the
tensor radii by this product theorem rather than exhaustively recomputing
them.  That is mathematically valid, but a one-line assertion or proof comment
in the script would make the verification boundary explicit.

## 3. Rate--distortion report

### Deterministic packing converse and Bayesian/deficiency bounds — ACCEPT

The deterministic collision argument, mutual-information reduction, and
Fano bound have the correct directions and quantifiers.  Apply the same
strict-packing terminology correction described under QR-D.

### Theorem 4.1, finite Ising pinned-query rate--distortion — ACCEPT

The flip-set bound is exact.  Flipping `k` coordinates loses `2Mk` in field
value, while at most `k(n-k)` quadratic edges change, for interaction gain at
most `2ak(n-k)`.  Thus `M>a(n-1)` makes the queried `u` the unique maximizer
and gives 4.6.

The constant shift `c_A` has Walsh degree zero, so every degree-two
coefficient is exactly `a A_ij`.  Bessel's inequality and sign thresholding
give, pointwise,

\[
\frac{d_{Ham}(A,\widehat A(Z))}{N}\le d_Q(A,Z).
\]

Entropy subadditivity plus binary Fano and concavity then give
`I(A;Z)>=N[1-h_2(D)]`.  No independence of posterior edge errors is assumed.

The query quantifier is essential and is correctly stated: one common
transcript must decode the entire function `u -> V(h^u)`, with distortion
averaged over a fresh uniform `U`.  The theorem does not apply to a separate
sketch for each preselected query.

For uniform error `<a`, each recovered Walsh coefficient has error `<a`, so
every sign is exact.  Storing `A` gives the matching `N`-bit upper bound.
Replace only the phrase "`2a`-packing" by "`2a`-separated" (or `r`-packing
for every `r<2a`).  At error exactly `a`, the `n=2` pair admits a midpoint
decoder, so the strict threshold is substantive.

### Bounded rank-one coupling variant — ACCEPT

For a nonconstant `z_i=u_i x_i`, the rank-one term loses
`2Lk(n-k)` and the unknown quadratic gains at most `2ak(n-k)`.
Thus `L>a` pins the pair `x=+/-u` and 4.13 follows.  The fact that only
`2^(n-1)` distinct rank-one queries remain does not weaken the Walsh decoder,
because `q_A(u)=q_A(-u)`.

### State-information comparison — ACCEPT

The coordinate queries return exactly `a s_i`, yielding the stated Hamming
reduction and `n[1-h_2(D)]` lower bound.  The common-optimizer fiber has size
at least `2^(N-n)` by pigeonhole (informative once `N>=n`), and the pinned
response remains injective on that fiber.

### Rare-success lemma — ACCEPT

Data processing to the success indicator gives a binary divergence.  With
joint success at least `1-delta` and product success at most `p`, dropping the
nonnegative failure term gives exactly

\[
I(\Theta;Z)\ge(1-\delta)\log_2(1/p)-h_2(\delta).
\]

### Posterior-polarization converse — ACCEPT

For each edge,

\[
\frac{1-|w_e|}{2}\le\frac{1-w_e^2}{2}\le\frac12.
\]

Monotonicity and concavity of binary entropy give 6.3; averaging and applying
concavity once more gives 6.4.  Conditioning the uniform source on an event
then gives 6.5 exactly.  The result requires small expected posterior
variance `E V_Z`, as stated.

### Vertex-prize Max-Cut corollary — ACCEPT

A flip set changes at most `k(n-k)` cut edges, while the field loses `2Mk`,
so `M>(n-1)/2` pins `u`.  The Walsh coefficient is exactly `-B_ij/2`.
Uniform value error `<1/4` therefore recovers every edge by thresholding, and
the worst-case lossless rate is exactly `N` bits.

## 4. Exact computation reruns

All generated files were written to `/tmp`; no repository result was
overwritten.  The repository virtual environment was used.

1. `entropy_overlap_lab.py` reran successfully.  Its output was byte-for-byte
   identical to `entropy_overlap_results.json`, with SHA-256
   `b7efd0139af06d054f4d1b66a389e0509d33e477d4686abda8e0cdacdaabcf80`.
   It reproduces the code radii `2,3`, the order-eight equal-energy-histogram
   caps `16,20`, and the census through order eight with no pair-signature
   collision having different one-vertex-response multisets.
2. `pinned_query_rate_verify.py` reran successfully.  Its output was
   byte-for-byte identical to `pinned_query_rate_results.json`, with SHA-256
   `82afc47acf84db5bd2a6bb3c2c7cca9f35a76a3a44bba891e38a14b560b7e6ce`.
   It checks `2,8,64,1024` landscapes at orders `2,3,4,5`, respectively, and
   recovers every edge sign from every complete response vector.
3. `build_quadratic_landscape_dataset.py` reran successfully.  Its output was
   byte-for-byte identical to `quadratic_landscape_order8.json`, with
   SHA-256
   `949987affbe23f341c3c8abc47c01ce8aa233eeb5b6a034d1a8a41dba87cbd92`.
   The `1044` rooted-gauge graph representatives partition into `243` exact
   pair-signature classes; bucket sizes sum to `1044`, all `243` signature
   hashes are distinct, the minimum absolute cap is `10`, and exactly two
   signature classes attain it.

Independent brute-force checks, separate from the supplied scripts, gave:

```text
codes {0: 4, 2: 12} {0: 4, 2: 12} radii 2 3
tensor powers 1,2,3: equal enumerators; radii (2,3), (4,6), (6,9)
Curie--Weiss m=2,4,6,8: constrained sums m-2; soft maxima m+1
n=2 pinned-query distance between the two landscapes: exactly 2a
```

The experiment census is finite evidence only.  It proves its exhaustive
claims through order eight after switching the first row positive and
quotienting the remaining graph by isomorphism; it does not establish an
all-order completeness theorem for pair signatures.

## 5. Required patch list

1. Replace Proposition 4.4 by the resolution-aware version above; do not
   claim detection for an arbitrary coarse fixed grid.
2. Everywhere replace strict "`2a`-packing" by "`2a`-separated", or by
   "`r`-packing for every `r<2a`."
3. Compare QR-A roofs as extended functions when feature polytopes differ,
   and qualify "minimal quotient" as equivalence-minimal, not bit-minimal.
4. Define whether the QR-D cover is internal and describe QR.11 as a
   restricted support-function pseudometric unless all support directions are
   queried.
5. In QR-E name the norm defining Hausdorff distance and call the interface
   fixed in ambient dimension/region, not in feature alphabet.
6. Align the positive-tail identity definition with the `liminf` locus used
   by Corollary 4.3, and mention the trivial zero quadratic boundary case.

After these corrections, no remaining theorem in the audited drafts is
falsified by the boundary cases or exact computations above.
