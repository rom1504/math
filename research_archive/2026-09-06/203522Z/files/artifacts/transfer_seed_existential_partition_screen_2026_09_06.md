# Existential partitions of exact minimizers: a frozen finite screen

Date: 2026-09-06. Finite evidence and a falsifiable sufficient criterion,
not a partition theorem or a convergence claim.

## 1. What survives the averaged restriction obstructions

An average or high-probability lower bound on random child caps does
not exclude ONE favorable partition. The quantifier that a convergence
argument could actually use is weaker still: for each prescribed pair
of comparable orders m,n, at least ONE exact minimizing parent of
order `N=m+n` admits the required partition. The parent may depend on
the prescribed pair; no uniform claim over all near-minimizers is needed.

A concrete sufficient target is: there are absolute `C,N0` such that,
for every `m,n>=N0`, `1/2<=m/n<=2`, an exact minimizing A of order N
has an m-subset T satisfying

```math
Q(A_T)^{2/3}+Q(A_{T^c})^{2/3}
\le M_N^{2/3}+C\sqrt N.                                  (P)
```

Since the two child caps dominate `M_m,M_n`, (P) immediately gives the
actual-minimum reverse-Fekete inequality. The balanced-tree argument
in `transfer_fresh_reverse_fekete_2026_09_06.md` then proves convergence.
This is an implication, not a proof of (P). Equal-half partitions alone
do not supply all prescribed comparable sizes and are not being
advertised as an all-order convergence criterion.

For a finite falsification screen, freeze `C=1` and require the bound
at every available prescribed comparable split. A parent/split for
which EVERY partition violates it falsifies the ALL-exact-minimizer
version, but not the weaker existential-parent version unless all
exact minimizing classes at that order have also been covered.
An asymptotic `O(sqrt N)` assertion with unspecified constant C would
require an unbounded sequence of excess-to-`sqrt N` ratios to refute;
one small counterexample cannot do that.

Two simpler finite screens were also frozen before looking at results:

- BP-half: at `N>=10`, some floor/ceil-half partition has BOTH child
  normalized caps strictly below one half.
- BP-exact: some floor/ceil-half partition has BOTH children exactly
  minimizing at their respective orders.

Neither simpler statement alone compares the parent cap to the sum of
powered child caps, and neither alone proves reverse Fekete.

## 2. Exact finite results

`computations/transfer_seed_balanced_minimizer_screen_2026_09_06.py`
replays all 15 cases of the archived isotropy screen: the stored
minimizing orbit representatives for orders three through eight and
one stored representative at each order nine through fourteen. It
checks 13,139 prescribed-partition instances. Boolean energies are
computed as integers. Powered-cap and square-root comparisons in the
`C=1` screen use exact rational enclosing intervals, not floating-point
inequality decisions.

Every stored case passes BP-exact and the `C=1` screen at every
prescribed comparable split. Every stored case of order ten through
fourteen passes BP-half. The balanced counts at these larger orders are:

| Order | Partitions | Both exact children | Both strictly subhalf |
| --- | ---: | ---: | ---: |
| 10 | 126 | 20 | 20 |
| 11 | 462 | 7 | 42 |
| 12 | 462 | 7 | 112 |
| 13 | 1,716 | 52 | 364 |
| 14 | 1,716 | 624 | 624 |

Global minimality is imported from each case's recorded provenance,
including the solver-dependent lower certificates at the larger
orders; those lower certificates were NOT rerun by this screen.
The new computation independently rechecks caps and partitions only.
Above order eight this is a sample of representatives, not a complete
classification of minimizing matrices. No theorem at growing order
is inferred from the observed counts.

The separate original order-14 signing in
`heuristic_m14_from_conference.json` gives the same balanced histogram:
624 pairs of caps `(9,9)` and 1,092 pairs `(11,11)`.
It already falsifies the ZERO-error version of (P) at that order,
because `21^2<8*9^2`. It does not falsify (P) with `C=1`.

One distinction appears at prescribed non-half sizes: the stored
order-12 exact minimizer has ZERO partitions of sizes four and eight
with both children exactly minimizing. Thus the all-parent BP-exact
statement cannot simply be extended to every prescribed comparable
split. Its four/eight powered `C=1` test nevertheless passes. The
global minimum status of the parent remains the imported certificate
status already specified above.

### 2.1 Supplementary order-15 and order-16 witness screen

The optional `--external-witnesses` run additionally reads the independently
replayed integer matrices in
`computations/results/transfer_adversary_external_order15_16_witness_verify_2026_09_06.json`.
These have caps 27 and 30 at orders fifteen and sixteen, respectively.
Their caps are rechecked again here. The order-15 global lower certificate
has NOT been replayed in this campaign; the order-16 witness is NOT
claimed exactly minimizing. They are therefore supplementary actual-sign
data, not newly certified exact-minimizer cases.

The order-15 witness has 104 of its 6,435 balanced partitions with child
caps `(9,10)`, exactly the known two child minima and strictly subhalf.
The order-16 witness has only THREE of its 6,435 balanced partitions
with child caps `(10,10)`; these are precisely its strict-subhalf
balanced pairs. Both witnesses pass the `C=1` inequality for every
prescribed comparable split. The complete extended run performs
53,465 exact partition evaluations. The low count at order sixteen
illustrates why existential and average assertions are materially
different, but supplies no growing-order rarity theorem.

## 3. Why the planted clique does not already decide this test

The independently audited near-minimizer construction changes at most
`binom(r,2)` signs, with `r=N^(2/3)`. For every prescribed T,

```math
|Q(A'_T)-Q(A_T)|\le r(r-1)=o(N^{3/2}).                   (1)
```

Thus, if a parent already has a balanced partition whose two child
normalized caps are bounded below one half by a FIXED margin, the
same partition survives the clique change for all large N. The clique
may be distributed across both shores; no argument forces an
existentially chosen partition to put the whole clique on one side.
Consequently the clique theorem alone is not a counterexample to
BP-half or to the existential exact-parent assertion (P).

There is also no license to extend (P) from exact to arbitrary
normalized near-minimizers. At balanced sizes, the universal positive
lower cap scale and the derivative of `u^(2/3)` bound the effect of
the replacement by `O(r^2/sqrt N)=O(N^(5/6))` on the powered scale.
This is much larger than the desired `O(sqrt N)` defect. Normalized
cap closeness `o(N^(3/2))` therefore has insufficient quantitative
precision to transfer (P). This sensitivity calculation is not itself
a proof that every favorable partition is destroyed.

The recorded scope remains: averaged failures do not decide an
existential selector; all-near-minimizer failures do not decide selected
exact minimizers; finite passing samples do not establish a recurrence.

## 4. Bounded cap-preserving neighbor and orientation screen

At the director's request, the companion probe also tests the immediately
falsifiable exact-minimizer assertion
`|max H_A-max(-H_A)|<=2`. A gap at least four at an exact minimizing
cap would be decisive. No such gap is assumed impossible in advance.

`computations/transfer_seed_balanced_neighbor_probe_2026_09_06.py`
first replays every stored minimizing orbit class at orders three
through eight. Their oriented gaps are respectively `2,0,0,0`,
then two for each of the three order-seven classes and zero for both
order-eight classes.

The probe then exhausts every one-, two-, and three-edge modification
of the stored minimizing witnesses at orders nine through twelve:

| Order | Changed signings tested | Same-cap neighbors | Largest gap among parent and retained neighbors |
| --- | ---: | ---: | ---: |
| 9 | 7,806 | 68 | 0 |
| 10 | 15,225 | 0 | 0 |
| 11 | 27,775 | 0 | 2 |
| 12 | 47,971 | 0 | 0 |

All 68 retained order-nine neighbors have a balanced pair of exact
minimizing children, with at least 21 such partitions in each case.
No orientation-gap or BP-exact counterexample was found. This is a
complete local Hamming-radius-three test around these selected
witnesses, not a classification of all exact minimizers at those
orders. Minimality remains imported from the recorded global lower
certificates. The source preserves the preliminary radius-two test,
which likewise retained no neighbors at orders eleven and twelve;
the final radius-three run strictly includes those candidates.
