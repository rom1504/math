# Blind structural audit of quadratic near-minimizers

Status: **FINITE EXPERIMENTAL AUDIT / HYPOTHESIS GENERATION ONLY / NOT A THEOREM**

Date: 2026-08-17 UTC

## Bottom line

Rigid conference flatness, exact cap-active isotropy, and hereditary minimization all fail on genuine exact or one-parity-step near-minimizers in the available orders. The two quantities that best survived deliberate falsification were instead:

1. a **soft near-cap shell** rather than the exact active set; and
2. a **coarse sixth-moment tail bound** rather than conference-like spectral flatness.

Both are only candidate implications. Orders at most 14, strong construction bias at orders 12--14, and explicit low-order counterexamples make an asymptotic claim unjustified.

The strongest negative witnesses are explicit and machine-verifiable:

- an order-9 signing with `Q=14=M_9+2`, conference defect `D4=0.95833`, normalized sixth moment `R6=12.94444`, and no improving edge flip; and
- an order-11 signing with `Q=19=M_11+2`, exactly one cap-active projective spin, maximal active-frame defect `10`, and no improving edge flip.

Thus near-minimality does not force a rigid algebraic matrix or a rich *exact* optimizer code, even at the next parity-compatible cap.

## Blind protocol and reproducibility

Before freezing the observables I used only the naked objective, the supplied values

`M_3,...,M_14=(3,4,4,5,9,10,12,13,17,18,20,21)`,

and matrix-bearing files and computation scripts under `computations/`. I did not inspect the near-minimizer prompt, favored synchronization hypothesis, theory/portfolio drafts, or PC.3 drafts. The freeze is recorded in

`extremal_information/experiments/nearmin_blind_observable_freeze.md`.

I read the campaign prompt only after that file was written. I did not use the theory/portfolio or PC.3 drafts in this audit.

Reproduce the run with

```bash
python3 extremal_information/experiments/nearmin_blind_structural_audit.py
```

The deterministic default run (seed `20260817`) takes about 37 seconds in the current environment and writes

`extremal_information/experiments/nearmin_blind_structural_results.json`.

Every reported cap is recomputed over all `2^(n-1)` projective spins. The run independently reproduced the archived root-gauged minimizer counts at orders 3--7: `2, 6, 12, 12, 3240`. Twelve random switching/permutation/global-sign transformations preserved every tested scalar observable to at most `3.2e-15`, and preserved the exact deletion and edge-response multisets.

Data generated or ingested in the default run:

| stratum | saved matrices / records | status |
|---|---:|---|
| repository matrices with cap `M_n` | 61 byte-distinct | exact cap recomputation; many are symmetry/construction duplicates |
| repository matrices with cap `M_n+2` | 9 | exact cap, heuristic or certified discovery provenance retained |
| all single-edge neighbors of one exact seed, `n=7,...,14` | 420 | exact cap; all satisfy `Q<=M_n+2` |
| cap-constrained adversarial-walk samples | 116 | heuristic discovery, exact cap |
| independently seeded greedy low-cap samples | 28 | heuristic discovery, exact cap |
| uniform controls | 48 per order | exact cap |
| uniform draws that happened to have cap at most `M_n+2` | 182 | mostly orders 3--9; one at order 11 |
| all cyclic-distance controls | `2^floor(n/2)` per order | exact cap |
| exhaustive root-gauged populations | all signings through order 7 | exact enumeration |

“Exact minimizer” below means the recomputed cap equals the supplied exact `M_n`. This audit does not independently reproduce every lower-bound certificate establishing those supplied values.

## Frozen observables

For `e_A(x)=sum_{i<j}a_ij x_i x_j`, `C=Q(A)`, and `N=binom(n,2)`, the main quantities were:

- conference defect
  `D4=||A^2-(n-1)I||_F^2/[n(n-1)^2]`;
- absolute two-walk correlations `|(A^2)_ij|/(n-2)` and spectral effective rank;
- absolute triangle bias;
- normalized landscape moments, especially
  `R6=E_x[e_A(x)^6]/N^3`;
- exact and soft boundary masses
  `b_j=Pr(|e_A(x)|>=C-j)` for `j=0,2,4`;
- exact active-code overlap/frame geometry;
- extremal local fields
  `z_i=sign(e_A(x))x_i(Ax)_i/C`;
- principal-deletion caps; and
- all single-edge-flip cap responses.

These are invariant under diagonal switching, vertex permutation, and global matrix sign, either directly or as unlabeled/absolute multisets.

## Exact minimizer inventory

The repository contains an authoritative exhaustive orbit classification only through order 8.

| `n` | `M_n` | minimizing root-gauged labeled signings | signed-permutation/global-sign classes | conference splits cover all classes? |
|---:|---:|---:|---:|---|
| 3 | 3 | 2 | 1 | yes |
| 4 | 4 | 6 | 1 | yes |
| 5 | 4 | 12 | 1 | yes |
| 6 | 5 | 12 | 1 | yes |
| 7 | 9 | 3240 | 3 | yes |
| 8 | 10 | 4200 | 2 | no |

The failure of conference coverage already at order 8 is an early warning against equating minimization with a single algebraic family.

At orders 9--14 there is no exhaustive orbit claim in the inspected results. I therefore inventory only available representatives. “Signature classes” below are separated by absolute spectrum, energy histogram up to reflection, deletion-cap multiset, `D4`, and triangle bias. Different signatures prove distinction; a signature collision is **not** asserted to prove switching equivalence.

| `n` | `M_n` | byte-distinct available exact representatives | distinguishable frozen-observable signatures |
|---:|---:|---:|---:|
| 9 | 12 | 6 | 5 |
| 10 | 13 | 9 | 2 |
| 11 | 17 | 9 | 2 |
| 12 | 18 | 8 | 1 |
| 13 | 20 | 2 | 1 |
| 14 | 21 | 3 | 1 |

The sources include exact MILP witnesses, certified bridge constructions, exhaustive extensions/restrictions, conference completions, and heuristic discoveries with exhaustively recomputed Boolean profiles. Source file, nested key, and original classification are retained beside every matrix in the JSON result.

The byte-distinct repository one-step-near inventory is:

| `n` | cap | representatives | discovery provenance |
|---:|---:|---:|---|
| 8 | 12 | 1 | exact finite arithmetic/cospectral-collision audit |
| 10 | 15 | 2 | proved GF(9) construction; exhaustive extension of an exact child |
| 11 | 19 | 1 | exhaustive extension of an exact child |
| 12 | 20 | 4 | exhaustive conference restriction; exact algebraic generator; heuristic search; exhaustive extension of a heuristic child |
| 13 | 22 | 1 | heuristic search |
| 14 | 23 | 0 | none in the pre-existing result files |

“Heuristic” here describes how the matrix was found. Its saved cap is still exhaustively recomputed in this audit. The independently generated greedy and cap-constrained-walk matrices are labelled the same way: heuristic discovery, exact finite cap.

## Exact cap stratification at order 7

The complete order-7 population gives the cleanest search-independent signal. Medians by cap are:

| cap | population | median `D4` | median `R6` | median `b_2` | median active-frame defect |
|---:|---:|---:|---:|---:|---:|
| 9 (`M_7`) | 3240 | 0.357 | 5.794 | 0.2031 | 2.571 |
| 11 | 16884 | 0.611 | 9.630 | 0.0781 | 4.286 |
| 13 | 9660 | 0.865 | 14.709 | 0.0312 | 6.000 |
| 15 | 2520 | 1.373 | 26.423 | 0.0156 | 6.000 |
| 17 | 420 | 2.008 | 48.917 | 0.0156 | 6.000 |
| 19 | 42 | 2.897 | 86.753 | 0.0156 | 6.000 |
| 21 | 2 | 4.167 | 151.023 | 0.0156 | 6.000 |

So cap, high moments, and exact-active degeneracy are strongly associated in the complete finite population. The order-7 one-step stratum is nevertheless already broad enough to reach the random median in several observables; association is not rigidity.

## Exact, one-step-near, and random comparison

The next table uses byte-distinct exact/near matrices from all repository and generated low-cap strata. `E` is cap `M_n`; `N` is cap `M_n+2`; `R` is the 48-sample uniform random control. Maxima are shown for defects/moments and minima for soft-shell mass. The near samples at orders 12--14 are heavily influenced by one-edge perturbations of structured exact seeds.

| `n` | `E max D4` | `N max D4` | `R median D4` | `E max R6` | `N max R6` | `R median R6` | `E min b_2` | `N min b_2` | `R q90 b_2` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 7 | .611 | .865 | .865 | 8.08 | 12.84 | 12.84 | .2031 | .0313 | .1094 |
| 8 | .245 | .735 | .816 | 5.16 | 10.41 | 12.02 | .2188 | .0625 | .0938 |
| 9 | .569 | .958 | .847 | 8.30 | 12.94 | 13.24 | .0781 | .0195 | .0363 |
| 10 | .395 | .553 | .948 | 7.24 | 8.91 | 14.51 | .2344 | .0234 | .0137 |
| 11 | .478 | .711 | .769 | 8.54 | 10.85 | 12.67 | .0400 | .0137 | .0103 |
| 12 | .165 | .386 | .871 | 6.45 | 8.22 | 13.59 | .0449 | .0146 | .0057 |
| 13 | .083 | .442 | .887 | 6.69 | 9.20 | 13.80 | .0190 | .0071 | .0013 |
| 14 | 0 | .081 | .899 | 6.63 | 7.18 | 14.22 | .0635 | .0249 | .0017 |

Two robust descriptive facts emerge:

- Exact minimizers have substantially suppressed sixth moments in every observed order. One-step-near samples generally interpolate toward random, and at order 9 almost reach the random median.
- Exact boundary multiplicity is fragile, but the one-lattice-step shell `b_2` remains enhanced over random controls at orders 10--14, including after adversarial searches aimed at reducing it.

These are rankings for further work, not asymptotic evidence. The large apparent separation at orders 12--14 may mainly be inherited from conference/double constructions.

## Counterexample-first verdict on the preregistered hypotheses

### H1: conference spectral/two-walk flatness — rigid version falsified

Exact minimizers themselves range from `D4=0` to `0.6111` in the authoritative small-order orbits. More decisively, the independently generated order-9 matrix with internal hash

`4f33c84ef851b44d2c80f33debad02e988f8694b94e4bb37ecdb4944c0696986`

has

`Q=14=M_9+2`, `D4=0.958333`, maximum singular deviation `0.961850`, effective-rank fraction `0.585460`, and `R6=12.944444`.

The uniform order-9 controls have median `D4=0.847222`, q10 effective-rank fraction `0.57631`, and median `R6=13.2428`. Thus this genuine one-step near-minimizer is random-like or worse on the strongest rigid spectral summaries. It also has no improving single-edge flip, so adding edgewise local optimality does not repair H1.

Its upper-triangle row encoding (`+` means `+1`, rows separated by `/`) is

```text
+-+--+--/+-+--++/++-+--/-++--/-+-+/-+-/+-/+
```

The full matrix is stored beside the hash in the JSON result.

### H2: exact-boundary anti-spikiness/isotropy — falsified; soft shell survives

The cap-constrained adversarial search found an order-11 matrix with hash

`1dc0388681c4720b9fd6f723ef4814d8e4653882897e3ecb669e9e07819e2ae7`

and

`Q=19=M_11+2`, one cap-active projective spin out of 1024, active-frame defect `10=n-1` (the maximum singleton value), and no improving edge flip.

Its upper-triangle encoding is

```text
++-+--++++/--++++-++/-++-+-+-/---+-++/-+--+-/-+++-/----/---/+-/+
```

This kills any implication based on multiplicity, overlap dispersion, or isotropy of the *exact* active code. But 25 projective spins lie within two energy units of the cap and 54 lie within four. That is why the surviving candidate uses a soft shell rather than exact maximizers.

### H3: hereditary robustness — falsified even for an exact minimizer

The certified bridge witness from `computations/results/bridge_5_6_sign1.json`, internal hash

`bb6584568fb39a4a870bbe8544d764e288e838d55cb5f85580cc7b565a6f76dc`,

has `Q=17=M_11`, while **all eleven** principal deletions have cap 17. Since `M_10=13`, none is exact or one-step-near at order 10. Therefore “most deletions of an exact minimizer are minimizers/one-step-near-minimizers” is false.

The weaker normalized observation `Q(A_{-v})-M_{n-1}=o(n^(3/2))` is not tested by these orders and is too close to the original cross-order problem to count as a discovered invariant here.

### H4: edgewise stationarity — true for exact minimizers but not informative

No exact minimizer can have an improving edge flip, simply because its cap is globally minimal. At cap `M_n+2`, an improving flip must land at `M_n`, so the improving fraction mostly measures adjacency to the exact minimizer set. An order-7 one-edge witness has 10 improving flips out of 21, while most higher-order sampled one-step witnesses have zero or one. This is neither uniform enough nor independent enough of `Q` to promote.

## The two strongest surviving candidate implications

These statements are deliberately schematic and are not claimed.

### Candidate 1: soft-shell response entropy

Let

`S_2(A)={x projective: |e_A(x)|>=Q(A)-2}`

and define its deficit rate

`P_shell(A)=-(1/n) log_2(|S_2(A)|/2^(n-1))`.

A serious target would be a bound of the form

`Q(A)<=M_n+epsilon n^(3/2)  =>  P_shell(A)<=1-c+g(epsilon)+o_n(1)`

for some `c>0` and `g(epsilon)->0`, or a weaker quantitatively useful polynomial lower bound on `|S_2(A)|`.

Evidence: at orders 10--14, the smallest observed one-step-near `b_2` exceeds the random-control q90 by factors about `1.7, 1.3, 2.6, 5.5, 14.3`, respectively. The exact-active singleton falsifier still has `|S_2|=25`.

Strongest limitation: the smallest observed one-step-near shell counts at orders 7--14 are `2,8,5,12,14,30,29,204`. These tiny orders cannot distinguish exponential, polynomial, or construction-specific growth. The order-9 count 5 is an especially strong warning. The next useful experiment is therefore a scale study, not a proof attempt based on the present constants.

Why it remains interesting: it is collective, switching-invariant, stable under the exact-active collapse above, and directly phrased as near-top response entropy rather than as another calculation of `Q`.

### Candidate 2: strict sixth-moment tail suppression

With `N=binom(n,2)`, let

`R6(A)=E_x[e_A(x)^6]/N^3`.

A coarse target is

`Q(A)<=M_n+epsilon n^(3/2)  =>  R6(A)<=15-delta+g(epsilon)+o_n(1)`

for some universal `delta>0`, rather than `A^2` being close to `(n-1)I`.

The reference value 15 is the Gaussian sixth moment and is close to the observed random medians at the larger orders. Exact minimizer maxima are at most `8.55` in the inspected orders; one-step-near maxima at orders 7--14 are at most `12.95`. The complete order-7 enumeration shows a clean monotone cap stratification of `R6`.

Strongest limitation: the order-9 witness above has `R6=12.9444`, essentially the random median, so any strict gap is at most modest. Order 4 even has a one-step-near example with `R6=27.11`, ruling out an unqualified finite-order statement. Moreover, a sixth-moment bound is useful only if a downstream response/composition argument genuinely consumes it; this audit does not establish that arrow.

Why it ranks above `D4`: the strong conference-defect claim is directly falsified by the order-9 witness, while the coarser high-tail moment envelope survived all order-7 exhaustive and order-8--14 adversarial samples.

## Other observables that did not earn promotion

- Triangle bias varies across exact minimizer orbits and adds no clean separation.
- Spectral effective rank overlaps the random lower tail at order 9.
- Exact active-frame isotropy is maximally false for the order-11 singleton witness.
- Extremal local-field effective-support fractions overlap random controls and partly move in the opposite direction from naive “balance”. The universal checks `z_i in [0,1]` and `sum_i z_i=2` held, but they are local-optimality identities/inequalities rather than discovered near-minimizer structure.
- Deletion profiles are highly non-hereditary in absolute finite gaps.
- Edge response is mostly a restatement of adjacency to the minimizer set at one parity step.

## Interpretation and limitations

The data favor a softened collective landscape statement, not a rigid matrix classification. A low cap can coexist with random-like `A^2` and a singleton exact optimizer. Any promising observable should average over a near-top band and/or use a higher landscape moment, so that one or a few exact-sign edits cannot erase it.

This conclusion has strict limits:

- only orders 3--14 were used;
- orbit classification is exhaustive only through order 8;
- random controls number only 48 per order;
- orders 12--14 are dominated by related conference/double witnesses and their one-edge neighborhoods;
- the constrained walks explore connected local components, not the entire low-cap set;
- byte-distinct matrices overcount switching/permutation duplicates;
- frozen-observable signature counts undercount or equal the true number of available orbits, and are not orbit certificates; and
- no finite trend is promoted to an asymptotic theorem.

The correct frontier label is therefore **NO THEOREM; TWO NARROWED EXPERIMENTAL TARGETS; THREE RIGID ROUTES FALSIFIED**.
