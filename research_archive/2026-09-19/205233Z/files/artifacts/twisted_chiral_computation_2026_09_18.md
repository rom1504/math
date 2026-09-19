# Twisted chiral computational track, 2026-09-18

Original campaign window: 17:09:11--20:09:11 UTC, September 18. A usage-limit
interruption suspended work after about 50 active minutes; authorized work
resumed September 19 at 18:58 UTC with about 2h10 active time remaining.
This track distinguishes exact
fixed-matrix caps, exhaustive fixed-child family optima, and separately stored
global minima. The convention is Q(A)=max_x |sum_{i<j} A_ij x_i x_j|.

## Inputs and algorithm

Stored minimizer class counts for orders 3 through 8 are 1,1,1,1,3,2.
The initial orders 9 and 10 seeds were the concrete matrices in
`exact_m9.json` and `exact_m10.json`. Complete independent censuses now
establish all nine minimizing classes at order9 and both classes at order10,
with globally optimal caps12 and13. The selected optimal-child class
representatives are therefore complete through order10.

We use B_ij=s_i s_j A_(p_i,p_j), with p mapping destination to source.
The parent is D=[[A,B+diag(d)],[B+diag(d),-A]]. Its energy is

    E(x,y)=q_A(x)-q_A(y)+x^T B y+d·(x*y).

The chiral map (x,y)->(y,-x) negates E. Together with global spin negation,
this allows absolute-value enumeration with x_0=y_0=+1. For each B, the
search computes the exact minimum/maximum base energy at every t=x*y
with t_0=+1, then exhausts all 2^r matching vectors d. Search fitness uses
the integer cap and a heuristic tie-break from profile extremizer counts;
the tie-break is not a mathematical bound. Exhaustive mode additionally
enumerates all r! 2^(r-1) signed permutations, deduplicating equal B matrices.

Every saved final witness is reconstructed and independently checked by
the pre-existing `exact_fixed_signing_gray.cpp`, which visits all 2^(2r-1)
projective parent configurations using integer arithmetic and no chiral
reduction. Sources and random seeds are retained with the generated JSON.

## Results

| Child r | Stored optimal-child class | Q(A) | Exact fixed-child family minimum | Distinct B tested |
|---|---|---|---|---|
| 3 | 0 | 3 | 5 | 4 |
| 4 | 0 | 4 | 10 | 48 |
| 5 | 0 | 4 | 13 | 192 |
| 6 | 0 | 5 | 18 | 384 |
| 7 | 0 | 9 | 21 | 53,760 |
| 7 | 1 | 9 | 21 | 23,040 |
| 7 | 2 | 9 | 25 | 26,880 |
| 8 | 0 | 10 | 30 | 215,040 |
| 8 | 1 | 10 | 30 | 322,560 |

Each B is tested against every matching d. These exhaustive fixed-child
optima do not claim global optimality beyond the separately stored global
values through order 14. In particular the global M16 is not established
by this experiment (the archive has the bracket 28 <= M16 <= 30).

The class-2 optimal child at r=7 cannot reach global parent optimum M14=21;
classes 0 and 1 can. Thus the quantifier "some optimal child" is strictly
weaker than "every optimal child" even in this small regime.

For every optimal child at r=5,6,8, the best parent ratio
Q(D)/(2^(3/2)Q(A)) is respectively 1.149048519, 1.272792206, 1.060660172.
This rules out exact normalized nonincrease without an additive error; it
does not rule out an asymptotic subleading-error theorem.

An independent NumPy checker (`twisted_chiral_independent_profile_2026_09_18.py`)
has independently exhausted every B and d through r=6, using direct
integer matrix products over all x_0=+1 and unrestricted y. It reproduces
the same family minima and distinct B counts. Family-cap histograms are
retained in `twisted_chiral_2026_09_18_independent.json`.

### Reproduction

The driver automatically compiles missing/stale scratch executables from
the tracked C++ sources with `g++ -O3 -march=native -std=c++17`. For example:

    .venv/bin/python computations/twisted_chiral_driver_2026_09_18.py --orders 3 4 5 6 7 8 --mode exhaustive --tag exhaustive_histograms
    .venv/bin/python computations/twisted_chiral_independent_profile_2026_09_18.py --orders 3 4 5 6

The first run takes about a minute on the campaign machine. The exact
family minimization compares caps lexicographically before heuristic
tie-breaks; no floating-point score decides the mathematical cap.

## Order-nine census and larger parents

The independent census `twisted_chiral_classify9_2026_09_18.cpp` exhausts
all 268,435,456 root-gauged order-nine signings and proves M9=12 with
1,607,760 minimizing root-gauged signings in exactly nine equivalence
classes. Root-gauged permutation/global-negation orbit sizes are

    362880, 362880, 362880, 120960, 181440, 45360, 120960, 45360, 5040.

The exhaustive counts agree exactly with independently generated orbit
sizes. Every surviving signing is checked on all 256 projective spins.
For a root-signing bit mask b and spin z, the direct popcount energy is

    q_b(z)=q_all_plus(z)-2 popcount(b)+4 popcount(b & negative_product_mask(z)).

Only a witnessed violation may stop the spin loop early. Global negation
followed by root regauging complements every free edge bit, so both the
permuted root code and its 28-bit complement enter each equivalence orbit.
The source and full census are retained in `twisted_chiral_2026_09_18_classify9.json`.

For archived `exact_m9.json`, the entire twisted family has minimum 33:
46,448,640 distinct B matrices, reduced to 23,226,624 classes under its
order-two signed automorphism group. The exact bound-mode run completed
in 432.641 seconds. Its witness is

    p=(0,1,2,3,5,4,6,7,8),
    s=d=(1,1,1,1,1,-1,-1,-1,-1).

The independent Gray verifier gives parent energies min=-33, max=33,
with 186 maximizers and 186 minimizers among 131,072 projective states.
The root agent's independent full Walsh transform also verified all
262,144 states. This is a fixed-child family optimum, not global M18.

All nine order-nine family traversals are complete: census classes
0,2,5,7,8 have family minimum33; classes1,3,4,6 have family minimum37.
The order-ten census and its arithmetically independent replay both
completed, exhausting all 2^36 root-gauged signings. Both identify exactly
two minimizing classes, canonical root codes225334962 and1294879660,
each with362880 root-gauged signings. Their exact twisted-family minima
are both44. Each class has185794560 distinct B matrices and9292808
automorphism-quotient representatives. Original-order and accelerated
state-order family runs independently agree at44 for the second class.

The order-nine independent census uses signed half-edge sums rather than
the popcount identity above and reproduces all class counts. The analogous
order-ten independent replay completed in two disjoint 2^35-sized shards;
each contains181440 members of each minimizing class. Results are
`twisted_chiral_2026_09_18_classify10_independent.json` and the associated
shard logs. No census was rerun after the September19 resumption.

The standardized complete finite table and machine-readable evidence index
are `artifacts/twisted_chiral_computation_table_2026_09_18.md` and
`computations/results/twisted_chiral_2026_09_18_summary.json`.

An independent direct-spin untwisted audit, with x_0=+1 and unrestricted y
and all matching d, separates the effects of twisting and child selection.
For r7 the untwisted minima are25,25,25 versus full-family21,21,25.
For r8 both untwisted minima are32 versus full-family30. Thus twists
genuinely help at parent14 and16. At r9 the conference-deletion class8
already reaches33 without twisting, whereas the other classes give39 or41
untwisted. Both optimal r10 classes already reach44 untwisted; twisting
does not improve either one. Source and results are
`twisted_chiral_untwisted_audit_2026_09_19.py` and
`twisted_chiral_untwisted_audit_2026_09_19.json`.

The exact switching-only subclass (p=id, all2^(r-1) switches and all2^r
matchings) is also complete. Permutations are indispensable for fixed
optimal-child classes r7 IDs0,1 (25 versus21), r8 both (32 versus30),
and r9 IDs0,2 (37 versus33), ID5 (39 versus33). Other stored classes
match their full-family minimum by switching alone. The table is
`artifacts/twisted_chiral_switching_table_2026_09_19.md`.

For J={i:s_i=-1}, rotating the coordinate pair at each i in J by
(x_i,y_i)=(yprime_i,-xprime_i) transforms the switching-only lift into
an untwisted lift of A with ONLY its internal-J edges reversed, and
matching d*s. This matrix identity passed120 random exact tests, and
every saved best and transformed untwisted witness passed the independent
Gray verifier. All1024 partial reversals of the conference10 child have
cap histogram {15:82,17:762,19:180}, so none is a globally optimal Q13
child. Results: `twisted_chiral_switching_audit_2026_09_19.json`.

### Selectable optimal children fail at order20

The conference matrix in `conference_order10_gf9.json` is a nonoptimal
child with Q(A)=15. Taking p=id, s=d=(1,...,1) yields a parent with exact
cap40, independently checked over all524288 projective states by the
Gray verifier and all1048576 states by the root agent's full Walsh code.
Both extreme energies have96 projective occurrences. Moreover the full
fixed-conference-child twisted family has exact minimum40, from all
2580480 distinct B matrices,1926 quotient representatives.

Consequently the strict finite separation is

    min_{Q(A)=M10=13, p,s,d} Q(D) =44 >40 >= M20.

This disproves the assertion that some optimal child always suffices for
the globally best parent. It does not identify M20. The conference-child
witness and full-family result are retained in
`twisted_chiral_2026_09_18_alternative_conference.json`.

The upper cap40 now also has an algebraic certificate independent of
enumeration: H=D+diag(-I10,+I10) is symmetric Hadamard, H^2=20I, and
the symmetric-Hadamard20 quadratic bound is |z^T H z|<=80. The saved
spin attains80. A standard-library-only verifier with the complete child
and spin is `twisted_chiral_hadamard20_certificate_2026_09_19.py`.
The alternate sparse-vector proof is in
`artifacts/twisted_chiral_hadamard20_quadratic_certificate_2026_09_19.md`;
it references the root's canonical Hadamard20 theorem.

One bounded joint search, started September19 after the completed census
audit, permits Q(A)<=17 and simultaneously varies A,p,s while exactly
optimizing all matching d. Its target is parent cap38, seeded at the
conference cap40 witness; the budget is3600 seconds, seed20260919.
Source `twisted_chiral_joint_search_2026_09_19.cpp` includes the frozen
integer profile evaluator. The companion driver rebuilds its executable,
records both source hashes, logs all best matrices, and independently
checks each best with the full projective Gray verifier. Failure of this
bounded heuristic is not an exclusion certificate.

At elapsed2024.52 seconds, this search found a second cap40 parent with
exactly the original conference child of cap15, but a different twist.
Its two projective extreme multiplicities
are60 each rather than the original96, proving these two parents are
inequivalent under signed permutations/global negation. The new parent
has no sign diagonal making it Hadamard: choosing h_0=+/-1 uniquely
forces all h_j from row0 of D^2+(h_i+h_j)D_ij=0, and both candidates
fail. The exact exclusion algebra and multiplicities are retained by
`twisted_chiral_joint_witness_structure_2026_09_19.py` and its JSON result.

The director subsequently authorized exactly one additional1200-second
continuation, AFTER the original3600-second run, expanding the child cap
limit to23 and starting at the best saved child/twist. The necessary bound
Q(parent)>=2Q(child)-10 means a target38 parent can still have odd child
cap19,21 or23; limiting the search to17 omits those possibilities.
The extended source/driver have separate filenames; the original source
is unchanged. A zero-second input/reconstruction smoke test passed, and
the actual continuation began September19 at about20:07UTC, with no
authorization for any further extension.

The original3600-second run completed with6387870 proposals,4139221 exact
evaluations,446149 accepted moves and1064 restarts. Its best cap remained40;
no cap38 was found. Candidate child-cap counts were
{13:1780,15:775487,17:3361953}. The parent-cap histogram and child-cap
histogram each total4139220, exactly one less than the evaluation count
because the initial seed evaluation is not in those histograms. The
proposal count minus2248650 child-bound rejections also equals4139220.
Repeated visits are included: none of these counts is a distinct-state
coverage certificate. Complete results and negative-search statistics are
`twisted_chiral_2026_09_19_joint_target38.json` and its JSONL log.

## Matching-independent finite obstructions

For the third order-seven minimizing class, all 26,880 bridge matrices B
have maximum conditional half-width exactly 24. The matching contribution
is constant within each t=x*y sector, so even arbitrary real sector
recentering cannot lower the cap below 24. Integral parent energies are
odd, proving the family lower bound 25 before optimizing matching signs.

For the optimal order-six child, independent direct enumeration gives

    core-cap histogram Q(D0): {18:104,22:280},
    conditional-half-width histogram: {16:104,18:280}.

Thus min_B Q(D0)=18 and min_B max_t half-width=16. The full signed
matching family minimum is also18. Results are in
`twisted_chiral_2026_09_18_independent_r6_core.json`.

At order9, the exact minimum conditional half-width is32 for both bad
classes3 and6, despite their full matching family minima37. Each width
traversal covers15482880 bridge matrices via the signed symmetry quotient;
the identity bridge already attains32. Therefore a width-only obstruction
of36 is false here: restrictions on simultaneous matching-realizable
sector centers are material. Result:
`twisted_chiral_2026_09_18_width_bad_census9.json`.

## Preserved implementation checkpoint

The search implementation is frozen as
`computations/twisted_chiral_search_frozen_2026_09_18.cpp`, SHA256
`e79f22aa4374851fcf8093462d8e1e46645b72b362608f7023214837f84e1813`.
Drivers automatically build that frozen source and record its checksum.
The initial bound-mode binary source, preceding extra width diagnostics,
is separately retained as `twisted_chiral_search_bound_initial_frozen_2026_09_18.cpp`.
Generated-ELF provenance is in
`twisted_chiral_2026_09_18_build_products_candidate.json`; this is only a
candidate for the root-owned archive manifest. The initial checkpoint
did not modify that shared manifest; the later explicitly authorized merge
is recorded below.

On September19, `twisted_chiral_build_replay_2026_09_19.py` rebuilt all
12 distinct retained C++ source files. The executable `.text` and `.rodata`
sections agree for all14 original scratch ELF products, including the
captured historical initial-bound executable. The candidate manifest covers
all26 scratch ELFs, including the12 independent replay builds; there are
no uncovered executables. The included frozen evaluator dependency of the
joint search is recorded explicitly. Replay evidence is
`twisted_chiral_build_replay_2026_09_19.json`.
At the director's explicit request, these26 reviewed entries were merged
into canonical `research_archive/reviewed_build_products.json`, preserving
unrelated entries. No archive snapshot files were edited. The subsequent
switching-only evaluator and prepared extended-search evaluator, together
with their independent replay builds, bring the current reviewed count
to30 ELF products from14 distinct sources. All code/data comparisons pass,
and these four additional entries were also merged under the director's
provenance authorization.

The slower redundant full-histogram runs for r9/r10 were intentionally
interrupted after preserving their verified best witnesses in
`twisted_chiral_stopped_histogram_runs_2026_09_18.json`. Those incomplete
runs carry no family-optimality claim. The complete bound-mode tests
check all d simultaneously with integer bitsets, stopping only when every
d has an explicit violating spin pair. They agree with full-profile
exhaustion for every minimizing child through order eight.

No result here proves a uniform construction estimate or convergence.

## Closing bounded search and stability diagnostic, September19

The sole authorized wider-child continuation completed exactly1200 seconds
at about20:27 UTC. It made1423529 proposals,1415920 exact evaluations,
7610 child-cap rejections and237 restarts, retaining cap40. Evaluated child
cap counts were {13:210,15:157707,17:485727,19:497767,21:222216,23:52292}.
They sum to calls minus the initial evaluation and include repeated visits.
No cap38 witness was found; this is not an exclusion certificate. There
were no further search extensions. Complete matrices and failure statistics
are preserved in `twisted_chiral_2026_09_19_joint_target38_child23.json`.

After this run stopped, a separately authorized3.48-second wind-tunnel
evaluated151 parent matrices:24 seeded uniform twists for each of six
child groups, plus seven known minima. Exact Gray histograms count positive
energies and positive strict one-spin local maxima, checking every parent
coordinate. Six independent direct NumPy histogram replays agree. At the
2sqrt2 Q(child) threshold the stability filter reduces average counts by
roughly6.4 to37.2 times, but all144 random twists violate the threshold
and every sampled stable-count mean remains above1. Sample averages are
heuristic diagnostics, not expectation certificates. Full data and the
other two thresholds are described in
`artifacts/twisted_chiral_stability_wind_tunnel_2026_09_19.md`.

The stability executable and its independent replay bring the canonical
reviewed build-product count to32 from15 source files. All executable
code/read-only-data comparisons pass. The original two cap40 parent
witnesses remain the only deduplicated best matrices; the wider run did
not produce a third one.

## Independent near-order proof audit, September19

`twisted_chiral_near_order_audit_2026_09_19.py` independently checks the
director's near-order cycle-repair theorem on105 actual full signings,
15 random cases at each n=4,...,10, with signed permutations and admissible
deleted sets. Every case passes the changed-column and both image-set
containments, exact row-norm identities, beta(B0-B1)<=6 L_union, incident
matrix norm bounds, and principal Q/beta restriction inequalities.
For45 cases n<=6, direct enumeration of all2^(2n) parent spins checks both
repaired deletion and completion from an independently chosen child twist.
These are pointwise checks, not assumed optimizers or a new search for F.

Every case also exhausts all Bernoulli reservoir subsets. The expected
row norm and event probabilities are rational numbers; the square-root
expectation bound is checked exactly by squaring its nonnegative excess.
The paired expectation, probability of size at least r, and existence of
a simultaneously good reservoir also pass. All input matrices, twists,
deleted sets, norms and rational expectations are retained in
`twisted_chiral_near_order_audit_2026_09_19.json`. This finite independent
audit supports the separate proof but does not replace that proof.

## Exactness of the accelerated family traversal

The quotient engine first visits all r! permutations of A and root-gauges
each resulting B. It deduplicates the root-gauged matrix codes, then
attaches every one of the 2^(r-1) switching vectors with s_0=1. These are
all distinct full B matrices: their root rows determine the switching
vector, and regauging determines the root-gauged B. No B is omitted.

Let H be the signed permutations U, modulo ±I, with U^T A U=±A.
For the plus sign, simultaneous conjugation of both parent halves acts on
B and permutes d, preserving Q. For the minus sign, additionally exchange
the two parent halves. Thus H acts by cap-preserving transformations of
the B family. The engine keeps the least edge-bit code in each orbit.
It explicitly counts its stabilizer and assigns weight |H|/|Stab(B)|;
the final weight sum must equal the full distinct-B count. Complete
weighted family-cap histograms agree with unquotiented enumeration for
every optimal-child class through order8.

For each surviving B, bound mode tests whether any matching d can have
cap at most T=best_cap-2. T has the correct energy parity. Write the
matching sign mask as b, and t's mask (including its fixed positive
coordinate0) as a. An encountered base energy e imposes exactly

    (r+e-T)/2 <= HammingDistance(a,b) <= (r+e+T)/2.

All endpoints are integers in these runs. Precomputed Hamming-distance
bitsets remove exactly the matching vectors violating that inequality.
The run abandons a B only when every matching has a witnessed violating
spin pair. A surviving matching causes full exact profile evaluation.
If the incumbent later improves, all earlier exclusions remain valid
because they excluded an even larger threshold. Consequently a completed
weighted traversal proves the recorded fixed-child family minimum,
independently of the heuristic fitness used for tie-breaking or search.
