# Near-minimizer blind observable freeze

Status: **FROZEN BEFORE READING THE NEAR-MINIMIZER PROMPT OR THEORY DRAFTS**

Freeze date: 2026-08-17 UTC

The only problem-specific inputs used before this freeze were

- the naked objective `Q(A)=max_x |sum_{i<j} a_ij x_i x_j|`;
- the supplied exact values `M_3,...,M_14=(3,4,4,5,9,10,12,13,17,18,20,21)`; and
- machine-readable matrices, profiles, orbit enumerations, and certification labels under `computations/results/`, plus computation scripts needed to interpret those files.

The audit treats diagonal switching, vertex permutation, and global matrix sign as nuisance symmetries. Every frozen observable is invariant under them, either exactly or after taking an absolute value / an unlabeled multiset.

## Frozen primary observables

For a symmetric zero-diagonal signing `A` of order `n`, put `e_A(x)=sum_{i<j}a_ij x_i x_j`, `C=Q(A)`, and `N=binom(n,2)`.

1. **Conference spectral defect**
   `D4 = ||A^2-(n-1)I||_F^2 / (n(n-1)^2)`.
   Also retain the normalized maximum singular deviation
   `Dop=max_i ||lambda_i(A)|/sqrt(n-1)-1|`.
2. **Pair-correlation / two-walk defect**
   The unlabeled multiset `| (A^2)_ij |/(n-2)`, `i<j`; summarize by mean, RMS, maximum, and zero fraction. This is the local decomposition of `D4`.
3. **Triangle bias**
   `T3=|tr(A^3)|/[n(n-1)(n-2)]`, the absolute mean signed triangle product.
4. **Boolean landscape moments and entropy**
   Over projective spins (`x_1=1`), retain the exact energy histogram, `E[e^4]/N^2`, `E[e^6]/N^3`, Shannon entropy of `e/C`, and the cap-normalized absolute-energy quantiles.
5. **Boundary mass**
   Fractions of projective spins with `|e|=C`, `|e|>=C-2`, and `|e|>=C-4` (using the actual parity lattice rather than treating these bands as continuous).
6. **Extremal-code geometry**
   For `X_*={x:x_1=1, |e_A(x)|=C}`, retain `|X_*|/2^(n-1)` and the pairwise multiset `|x.y|/n`; summarize its mean, maximum below the diagonal, and collision histogram.
7. **Extremal local-field profile**
   For every cap-active `x`, form `z_i=sign(e_A(x)) x_i(Ax)_i/C`. Retain the pooled unlabeled multiset and per-state min/max/variance. At a genuine absolute extremum these lie in `[0,1]` and sum to `2`; deviations therefore expose implementation errors.
8. **Principal-deletion robustness**
   The unlabeled multiset `{Q(A_{-v})/C:v in [n]}` and, where supplied exact `M_(n-1)` is relevant, `{Q(A_{-v})-M_(n-1)}`.
9. **Single-edge perturbation response**
   The unlabeled multiset `{Q(A^(ij))-C:i<j}`, where one edge sign is flipped; summarize improving / neutral / worsening fractions and extrema.

## Frozen secondary summaries

- normalized cap `C/sqrt(N)` and asymptotic-scale cap `C/n^(3/2)`;
- absolute determinant on the natural scale `(n-1)^(n/2)` and spectral effective rank;
- number and balance of positive versus negative cap states (reported only up to swapping);
- orbit/source/certification metadata, never silently pooled across exhaustive, certified-witness, heuristic, or sampled populations.

## Pre-registered candidate implications to try to falsify

These are hypotheses for finite-data testing, not claims.

- **H1 (spectral/two-walk flatness):** a cap within one parity step of `M_n` forces small `D4`, equivalently collectively small off-diagonal entries of `A^2`.
- **H2 (landscape anti-spikiness):** a cap within one parity step of `M_n` forces non-negligible near-boundary mass and a dispersed cap-active code, rather than one isolated spike pair.
- **H3 (hereditary robustness):** most principal deletions of a near-minimizer remain near-minimal at order `n-1`.
- **H4 (edgewise stationarity):** exact minimizers have no improving edge flip and near-minimizers have few improving edge flips, with a characteristic neutral/worsening response profile.

Controls and falsification strata were also frozen: exhaustive root-gauged signings where feasible, cap-stratified random signings, uniformly random controls, cyclic/Toeplitz controls, the all-one switching class, all available exact minimizer orbit representatives, and independently generated local-search low-cap samples. Finite patterns will not be promoted to theorems.
