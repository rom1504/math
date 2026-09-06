# Independent audit of `meanfield_response_quantization.md`

**Scope.** I checked the draft against Corollary 4.4 and independently
rederived the binary top-\(k\), response-duality, histogram, quadratic,
and quotient claims. I did not use other repository material.

## Bottom-line verdict

The mathematical core is sound. In particular, the following nontrivial
claims check out:

1. the uniform-field response is injective on the binary local top-\(k\)
   profiles;
2. fixed-grid histograms compose exactly and have the stated stars-and-bars
   count;
3. the \(n\eta/2\) value error is independent of merge depth when rounding
   happens once per microscopic site;
4. for fixed pair coefficient \(J\), equality of linear responses is a
   congruence for the bilinear merge, so the concave-envelope quotient really
   is associative;
5. the two-site strict-quotient example is correct (when \(B>0\)); and
6. \(J\ge 2B\) is a valid size-uniform sufficient condition for endpoint
   collapse in the draft's maximization convention.

I found no fatal counterexample. I would classify the draft as **correct
subject to revisions**, mainly to distinguish several notions of state
count, repair one asymptotic bound outside its implicit parameter regime,
make the maximization and normalization conventions unmistakable, and avoid
claiming sharpness or novelty beyond what is proved.

## Revisions that should be made before treating it as a theorem statement

### 1. Repair the \(M\)-dependence in (1.6)

With

\[
M=1+\left\lceil\frac{2B}{\eta}\right\rceil,
\]

the unconditional bound is

\[
\log_2|\mathcal S_{n,\eta}|
 \le (M-1)\log_2(n+1)
 =O\!\left((1+B/\eta)\log(n+1)\right).
\]

The draft's \(O((B/\eta)\log(n+1))\) is valid only after assuming, for
example, \(0<\eta\le 2B\) (and \(B>0\)). If \(\eta\gg B\), \(M\) remains at
least two under the chosen endpoint-grid convention while \(B/\eta\) tends
to zero.

The degenerate case \(B=0\) also needs a separate sentence: take the
singleton grid \(\{0\}\), one histogram at each fixed mass, and zero error.
As written, (4.1) becomes \(0/0\). The construction also assumes
\(\eta>0\); \(\eta=0\) does not yield a finite grid for arbitrary real
fields.

### 2. Say that the full chemical-potential range is essential, not merely convenient

For known mass and no quadratic term, \(\lambda\in[-B,B]\) is indeed enough.
For unknown mass it is not: blocks consisting entirely of fields \(-B\) all
have zero response on that closed interval, regardless of mass. With the
quadratic term, supporting fields can have magnitude of order \(|J|n\); for
the endpoint-collapsed profile the transition occurs at

\[
\lambda_*=-\frac1n\left(\sum_i h_i+J\binom n2\right).
\]

Thus the sentence in Section 2 calling unbounded \(\lambda\) "convenient"
should say it is required unless mass is separately stored and a sufficiently
large model-dependent interval is declared. This also matters to the
coarsest-quotient claim: restricting the terminal coefficient range generally
produces a coarser observational equivalence.

### 3. Make the max-energy convention explicit

Every proof in the draft **maximizes** \(H\). The signs in Sections 6.2--6.3
are correct for a reward/score, or for the negative of a conventionally
minimized Hamiltonian. In standard physics language a ground-state energy is
usually minimized. Under negation, the roles of positive and negative
curvature reverse. The document should either call \(H\) a score, say that
the ground state of \(-H\) is being computed, or give the min-plus translation.

This is especially important before calling a sign of \(J\) ferromagnetic,
attractive, or repulsive; none of those physical labels is normalization- and
sign-free.

### 4. Qualify the \(J\ge2B\) threshold and equality case

For the precise model

\[
q(k)=p(k)+J\binom{k}{2},
\]

the calculation in (6.18) is correct: \(J\ge2B\) makes every raw sequence
discrete convex, hence its least concave majorant is its endpoint chord.
This is a size-uniform sufficient condition, sharp if blocks of size two must
be covered.

It is not the sharp threshold for endpoint response collapse at a fixed mass
\(n\ge2\). Uniformly over all fields in \([-B,B]\), the exact fixed-\(n\)
threshold is

\[
J\ge \frac{4B}{n}.
\]

Indeed, writing \(P_k=\sum_{i\le k}a_i\) and \(S=\sum_i a_i\),

\[
P_k-\frac{k}{n}S\le \frac{2B\,k(n-k)}{n},
\]

with equality for \(k\) fields equal to \(B\) and the remainder equal to
\(-B\). The quadratic endpoint chord exceeds \(q(k)\) exactly when
\((J/2)k(n-k)\) exceeds that heterogeneity term.

At \(J=2B\), intermediate occupancies can tie (certainly for \(n=2\), and
also in degenerate cases), so the safe statement is that the **response is
determined by the endpoints**, not that only endpoint configurations can
ever maximize.

### 5. Do not transport the number \(2B\) across quadratic normalizations

The draft's \(J\) in (6.5) is the coefficient of each unordered occupied
pair. For

\[
\Phi(k)=\alpha k^2+\beta k,
\]

the effective cross/pair coefficient is \(2\alpha\), as (6.20) correctly
states. Consequently:

* for \(Jk^2/N\), the effective pair coefficient is \(2J/N\) (and
  \(k^2=2\binom{k}{2}+k\)); at final mass \(n=N\), the sharp endpoint-roof
  condition is \(J\ge2B\), whereas making every partial block of every size
  collapse under that same coefficient requires the much stronger
  \(J\ge BN\);
* for \((J/N)\binom{k}{2}\), the corresponding final-mass threshold is
  \(J\ge4B\), while size-uniform collapse of all partial blocks requires
  \(J\ge2BN\).

These statements use the draft's maximization convention. They illustrate
why the unnormalized condition \(J\ge2B\) must always be quoted together
with (6.5), not as a generic Curie--Weiss threshold.

For Ising variables, the field conversion in Section 2 is correct, but pair
constants should also be mentioned if thresholds are compared. Specifically,

\[
\sum_{i<j}\sigma_i\sigma_j
=4\binom K2-2(n-1)K+\text{constant}.
\]

Thus an Ising pair coefficient acquires a factor four in the occupancy pair
coefficient, plus a known mass-dependent linear term and a constant.

### 6. Separate finite table width, representation length, and number of possible states

Corollary 4.4's "polynomial-state" assertion counts attainable feature
fibres in a **single** block table. In this binary specialization that table
has \(n+1\) entries. With arbitrary real fields there are uncountably many
possible exact response summaries; no finite exact-state cardinality has
been proved.

For an \(M\)-level grid, the draft correctly counts response-equivalence
classes at fixed mass:

\[
\#\mathcal S_n=\binom{n+M-1}{M-1}.
\]

If all masses \(0\le n\le N\) are admitted, the exact count is

\[
\sum_{n=0}^N\binom{n+M-1}{M-1}=\binom{N+M}{M}.
\]

In the strong-coupling endpoint-collapse regime, an arithmetic \(M\)-point
grid has only

\[
n(M-1)+1
\]

distinct fixed-mass quotient states, because only the total field sum
remains. With real fields the "single number" in Section 6.3 still ranges
over a continuum. These distinctions prevent "number of table coordinates,"
"bits used to encode one table," and "number of equivalence classes" from
being conflated.

The sentence saying concavity "reduces" an \((n+1)\)-entry profile to an
empirical multiset is structural rather than a size reduction: the multiset
still contains \(n\) real slopes.

### 7. Retained additive baselines affect finiteness and error accounting

Section 2 correctly notes the baseline
\(c_A=\sum_i e_i(0)\). If it is retained and arbitrary real-valued, then the
quantized slope histogram alone is not a finite state: \(c_A\) must be stored
exactly, quantized with its own error budget, or explicitly quotiented out.
The histogram state counts in Sections 4--5 apply to the anchored
\(H(0)=0\) convention.

### 8. Narrow the broad aggregate-term wording

The configurationwise estimate (4.7) remains true in the presence of any
identically evaluated interaction that does not itself change when fields
are rounded. That is an error-stability statement. It does not by itself
prove that one anonymous root histogram contains enough information to
evaluate an arbitrary term depending separately on occupancies of old named
blocks. The fixed total-occupancy quadratic and scalar terminal field are
safe; more general labelled or hierarchical terms require the additional
state needed to evaluate them. Section 2's label caveat and Section 7.2
already point in the right direction, but Theorem 4.1 should use the narrower
wording.

### 9. Separate affinity is sufficient, not logically necessary

Separate affinity of the cross term is the clean general condition used by
the mixture/conjugacy proof. It is not a theorem that roof closure occurs
*only if* the kernel is separately affine; special non-biaffine kernels or
restricted realizable families can close for other reasons. Section 6.4
should say "is guaranteed by this argument when" rather than a literal
"only when." Section 7.2's invitation to prove a model-specific closure
identity is the correct qualification.

### 10. Small presentation defects

There are several vertical-tab/control-character corruptions before the word
`arepsilon` in Sections 5 and the opening summary. They should be replaced
by proper \(\varepsilon\). The packing should also define \(K_\delta\)
formally. Its constants themselves check out, including the need for strict
error below \(\Delta/2\) in the microscopic lower bound.

## Claim-by-claim mathematical audit

### Coarsest context state: verified with the declared scope

For additive union,

\[
\mathcal R_A(C,\lambda)=R_A(\lambda)+R_C(\lambda).
\]

The empty future makes equality of \(R_A\) necessary, and equality of
\(R_A\) plainly suffices for every declared context. Since \(p_A\) is closed
discrete concave, conjugacy recovers it and its slopes. Thus equality of
responses, profiles, and field multisets is the same contextual equivalence.
Calling any injective encoding of that equivalence class "coarsest" is
correct in the usual deterministic/Myhill--Nerode sense.

It is not coarsest for futures with separately addressable old blocks,
site-labelled fields, nonlinear terminal lookup/pinning terms, restricted
\(\lambda\) ranges, or a changing interaction coefficient. The draft mostly
states these scope limits correctly.

Also note that Corollary 4.4 uses "roof over a feature" for the fibrewise
conditional maximum. In Section 6 of the draft, "roof" means the least
concave majorant of that raw conditional table. Those are different objects
when \(J>0\). The draft defines both, but consistently saying "raw fibre
roof" versus "concave envelope" would avoid attributing the strict quotient
directly to Corollary 4.4.

### Biconjugacy and exact range: verified

For any finite sequence \(q(0),\ldots,q(n)\),

\[
L_q(\lambda)=\max_k(q(k)+\lambda k),\qquad
q^{\mathrm{cav}}(u)=\inf_\lambda(L_q(\lambda)-\lambda u)
\]

give the least upper-semicontinuous concave majorant on \([0,n]\).
Raw coordinates are recovered exactly iff the sequence is discrete
concave. This validates (1.3) and (6.3)--(6.4), including tied slopes.

For clarity, the exact realizable range in the local binary model is:

* \(p(0)=0\), and \(p(k)-p(k-1)\) is a nonincreasing sequence in
  \([-B,B]\); every such sequence is realized by its slopes;
* equivalently,
  \(R(\lambda)=\sum_i(h_i+\lambda)_+\), so \(R\) is convex piecewise linear,
  equals zero for \(\lambda\le-B\), equals
  \(n\lambda+\sum_i h_i\) for \(\lambda\ge B\), and its distributional
  second derivative is an integer atomic measure of mass \(n\), supported
  on \([-B,B]\), with atoms at \(-h_i\).

For fixed quadratic coefficient \(J\), a raw profile is realizable exactly
when

\[
b_k:=q(k)-q(k-1)-J(k-1)
\]

is nonincreasing and lies in \([-B,B]\). Not every abstract concave roof is
necessarily the roof of a bounded-field microscopic block, although the
\(\star_J\) operation can consistently be defined on a larger roof class.

### Exact and approximate counts: verified subject to the distinctions above

The fixed-mass grid count (1.5) is exact: every weak composition is realized
by taking the corresponding numbers of fields at the grid values, and
different histograms have different sorted slopes and responses. The
microscopic lower bound is also exact: response distances between distinct
arithmetic-grid histograms are at least \(\Delta\), so error strictly below
\(\Delta/2\) cannot identify two of them.

The macroscopic packing in Section 5.2 is algebraically valid. Under
\(q\le n/8\) and \(q^2\le B/(64\varepsilon)\),

\[
s\ge\frac{n}{8q},\qquad
sw\ge\frac{Bn}{16q^2}\ge4\varepsilon n,
\]

and the disjoint triangular response perturbations give \(2^q\) pairwise
separated systems. The stated
\(\Omega(\min\{n,\sqrt{B/\varepsilon}\})\) conclusion follows when a
positive integer \(q\) of that order exists. It does **not** match the
histogram entropy at macroscopic error, exactly as the draft says.

For asymptotic claims such as (4.8)--(4.9), state explicitly that \(B>0\) is
fixed (or give the corresponding assumptions if \(B=B_N\)).

### Nonaccumulating quantization error: verified

If \(|h_i-Q(h_i)|\le\eta/2\), then for every configuration with occupancy
set \(S\),

\[
|H(x)-\widetilde H(x)|\le |S|\eta/2\le N\eta/2.
\]

Maximizing preserves the same bound even when a fixed interaction is
present in both systems. Exact addition of histograms rounds no site a
second time. Therefore the absolute error grows with the number of leaves,
not with merge depth or bracketing; the normalized error is at most
\(\eta/2\). "Nonaccumulating" should continue to be qualified this way,
because the absolute error does accumulate linearly in total mass, and
re-quantizing after merges would be a different algorithm.

The constant \(1/2\) is tied to nearest rounding on an endpoint grid of
spacing at most \(\eta\). Floor/ceiling rounding gives an \(n\eta\) bound
instead. For a target **absolute** root error \(\epsilon\), one needs
\(\eta\le2\epsilon/N\); a fixed grid only gives a fixed per-site error.

### Quadratic signs: verified for maximization

With decreasing fields \(a_k\),

\[
\Delta q(k)=a_k+J(k-1),\qquad
\Delta q(k+1)-\Delta q(k)=a_{k+1}-a_k+J.
\]

Hence discrete concavity is exactly
\(a_k-a_{k+1}\ge J\). Every \(J\le0\) preserves concavity, while no
\(J>0\) preserves it uniformly over all field multisets because ties already
violate the condition. These signs are all correct under maximization.

### Roof congruence and associativity: verified

A particularly direct independent proof of the key congruence is

\[
L_{A\sqcup C}(\lambda)
=\max_\ell\left\{q_C(\ell)+\lambda\ell
                  +L_A(\lambda+J\ell)\right\}.
\]

Thus replacing \(q_A\) by any other profile with the same linear response
cannot change the parent's response; the symmetric formula proves the same
for \(C\). Equality of response functions is therefore a two-sided
congruence. The raw cocycle

\[
Jk\ell+J(k+\ell)r=J\ell r+Jk(\ell+r)
\]

makes raw merging associative, and quotienting by a congruence makes the
induced roof law associative. This validates Theorem 6.1 for every sign of
fixed \(J\). The independent-mixture proof in the draft is also sound:
the outer concavification is essential because a fixed mean total may use
correlated pure occupancy pairs even though product mixtures suffice to
identify every linear support value.

### Strict quotient: verified

For \(B>0\), \(J>0\), and
\(0<a<\min(B,J/2)\), the profiles

\[
(0,0,J),\qquad (0,a,J)
\]

are distinct while their middle coordinates lie strictly below the common
endpoint chord. Both responses are

\[
\max\{0,J+2\lambda\}.
\]

The example is therefore a genuine strict contextual quotient, not just a
non-unique optimizer example. Add the hypothesis \(B>0\); no such choice of
\(a\) exists when \(B=0\).

## Classical content and possible overclaim

Most ingredients are classical:

* top-\(k\) order statistics and the identity
  \(R(\lambda)=\sum_i(h_i+\lambda)_+\) (a stop-loss/hinge transform);
* Legendre--Fenchel biconjugacy and upper concave hulls;
* max-plus convolution of cardinality profiles and sorted-list union;
* histogram monoids, stars-and-bars counting, and deterministic scalar
  quantization error;
* the all-to-all cardinality reduction and the bilinear pair-count cocycle;
* contextual equivalence/coarsest deterministic quotients, which are the
  standard observational or Myhill--Nerode construction once the context
  class is fixed.

The useful contribution of the draft is the clean synthesis: it declares a
nonlookup context family, proves that the linear response quotient is a
congruence for fixed bilinear mean-field merging, and keeps exact-grid and
macroscopic-error lower bounds separate. Corollary 4.4 already supplies the
exact polynomial-width fibre DP; Theorem 6.1 is an additional minimality
result caused by narrowing future observations to linear terminal fields.

Without a literature review, phrases such as "missing rigorous benchmark,"
"stronger," or "optimality remains open" should be scoped to this report.
In particular, the report proves that its packing does not settle the
macroscopic metric entropy; it does not establish that the question is open
in the literature. A safe formulation is "not resolved by the present
argument." Likewise, no asymptotic information-theoretic novelty should be
claimed from (1.5)--(1.7): those are standard histogram and rounding bounds,
albeit correctly deployed here.

## Recommended final disposition

Keep the central theorems and proofs. Before publication or use as a formal
benchmark, make the parameter-domain, maximization, baseline, and coupling
normalizations explicit; change (1.6) to include the additive \(1\); call
\(J\ge2B\) a size-uniform sufficient threshold and optionally state the
sharp fixed-mass \(4B/n\) threshold; distinguish Corollary 4.4's fibrewise
roof from the concave-envelope quotient; and scope all novelty/open-problem
language. With those changes, the proposed benchmark theorem is rigorous.
