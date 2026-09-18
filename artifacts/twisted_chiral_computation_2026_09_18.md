# Twisted chiral computational track, 2026-09-18

Campaign window: 17:09:11--20:09:11 UTC. This track distinguishes exact
fixed-matrix caps, exhaustive fixed-child family optima, and separately stored
global minima. The convention is Q(A)=max_x |sum_{i<j} A_ij x_i x_j|.

## Inputs and algorithm

Stored minimizer class counts for orders 3 through 8 are 1,1,1,1,3,2.
Orders 9 and 10 use the concrete matrices in `exact_m9.json` and
`exact_m10.json`; those files record globally optimal caps 12 and 13.
The selected child representatives are only complete through order 8.

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

The stored order-ten optimal child has a verified parent20 witness cap44.
Its exact family traversal remains ongoing. A separate CP-SAT sampler
found two inequivalent optimal order-ten children; a bounded two-shard
exhaustive 2^36-signing census is running with a 45-minute budget.
All nine order-nine child classes are being covered by exact family runs.

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

## Preserved implementation checkpoint

The search implementation is frozen as
`computations/twisted_chiral_search_frozen_2026_09_18.cpp`, SHA256
`e79f22aa4374851fcf8093462d8e1e46645b72b362608f7023214837f84e1813`.
Drivers automatically build that frozen source and record its checksum.
The initial bound-mode binary source, preceding extra width diagnostics,
is separately retained as `twisted_chiral_search_bound_initial_frozen_2026_09_18.cpp`.
Generated-ELF provenance is in
`twisted_chiral_2026_09_18_build_products_candidate.json`; this is only a
candidate for the root-owned archive manifest, which this agent does not edit.

The slower redundant full-histogram runs for r9/r10 were intentionally
interrupted after preserving their verified best witnesses in
`twisted_chiral_stopped_histogram_runs_2026_09_18.json`. Those incomplete
runs carry no family-optimality claim. The complete bound-mode tests
check all d simultaneously with integer bitsets, stopping only when every
d has an explicit violating spin pair. They agree with full-profile
exhaustion for every minimizing child through order eight.

No result here proves a uniform construction estimate or convergence.

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
