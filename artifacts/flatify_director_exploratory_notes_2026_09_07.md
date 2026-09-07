# Director exploratory notes — preservation, not theorem claims

Status: unfinished/rejected routes, 2026-09-07 constructive campaign.
Keep these notes separate from the verified upper theorem and classification.

## Rediscoveries checked against the archive

- Closing audit: deriving a small-coordinate deletion with bounded operator
  norm from SDP/Grothendieck weights only rediscovers the existing
  Grothendieck--Pietsch regular-core theorem. At vanishing deleted fraction
  its operator constant diverges; it does not repair the opposite-spectral
  normalization loss. No new progress is claimed from this rediscovery.
- The full quenched comparison suggests possible fourth-order refinements
  using global sign symmetry. No fourth-cumulant theorem was proved, and
  improving the transfer exponent would not itself pay the Gaussian parent
  optimum. This was not promoted over the actual convergence obligation.

- H2 conjugation exposes all clique-flipped child objectives; already proved
  in `decisive_independent_h2_exact_profiles_and_algebra_2026_09_07.md`.
  The new pair-rotation note classifies the proposed additional angle freedom.
- Full rank-one weave column signs erase arbitrary off-diagonal seed signs;
  matching fibre modes force sqrt(p)/2 after arbitrary balanced restriction.
  Both are already in `transfer_adversary_seed_loss_2026_09_06.md`.
  Trying to balance diagonal fibre signs does not remove matching witnesses.
- With unequal fibre sizes k_i, the same two-polarity matching argument gives
  Q >= max_matching sum_(i,j matched) k_i k_j. Sorting and pairing adjacent
  sizes makes this at least .5 sum_i k_i^2-O(m^2) when 0<=k_i<=m.
  This elementary extension is not a solution mechanism.
- A common-amplitude source is required for unequal-retention row estimates.
  E_t is convex in its source (rate-distortion J is an infimum of linear
  source functionals), so heterogeneous retention cannot improve that scalar
  certificate. The independent agent derived and wrote the precise version.
- Grothendieck/Pietsch deletion gives bounded operator norm after deleting a
  fixed small vertex fraction; random refill pays a small leading constant.
  This was already proved repeatedly in the archive, including
  `grothendieck_pietsch_spectral_regularization.md`. It does not yield the
  missing all-order action realization or exact flatification.

## Concrete but unresolved operations

1. One-vertex completion minimizes
   max_x (|H_A(x)|+|v dot x|) over full sign rows v. Choosing v from old rows
   relates the cost to local fields. Their average at one extremizer is
   O(sqrt(n)) for low-cap children, but exchanging the extremizer maximum
   with row selection is unjustified. No O(sqrt(n)) insertion theorem proved.
2. MUB conjugation preserves a seed's spectral law while giving full sign
   cross blocks. Spherical relaxation pays an SDP rather than Q(A); it cannot
   be used as the desired cap bound. The constructive and adversarial tracks
   are testing the actual rotated cube, not only this relaxation.
3. A hidden-partition variance martingale could start at the uniform profile
   and end at a true balanced block profile, avoiding illegal free reselection
   of a partition in the middle of a deterministic interpolation. What is
   missing is the sign/size of the optimized-pressure generator. The existing
   fixed-path PSD-plus-flip-slack formula does not supply that sign; width
   versus absolute cap is an additional gap. No new state is being promoted.
4. Majority/nonlinear combinations of several independent weave channels
   might remove individual matching witnesses, but their first-degree signal
   is diluted and higher-degree noise may restore a large cap. No useful
   inequality was derived. Do not treat this as a viable construction yet.
5. Orthogonally invariant spectral ensembles have a natural variational
   language, but spectrum alone does not preserve Boolean maxima, and known
   high-temperature TAP results do not give all-temperature adversarial sign
   realization. No imported theorem closes the target.

## Literature reconnaissance, no theorem imported from these searches

- Volberg, *An estimate of Sidon constant for complex polynomials with
  unimodular coefficients*, https://arxiv.org/abs/2205.04936 . Complex-domain
  and coefficient-norm results must not be identified with this real absolute
  quadratic optimum. Full mapping not established.
- Fan–Li–Sen, *TAP equations for orthogonally invariant spin glasses at high
  temperature*, https://arxiv.org/abs/2202.09325 . The stated regime is high
  temperature, with random orthogonal disorder; it is not a theorem about
  globally optimized full signs at zero temperature.
- *The Gale–Berlekamp game for complex Hadamard matrices*,
  https://arxiv.org/abs/1310.1810 . Complex bilinear glow is not the real
  same-spin cap. Existing project complex-phase obstructions should be read
  before treating this as a new connection.

These are retained to prevent repeated optimistic reinterpretation. They
are not an exhaustive literature audit and not evidence of impossibility.

## Second-block independent attempts and their actual gaps

- Small Frobenius distance to the sign cube does not directly give favorable
  rounding when coefficients exceed one. Clipping the exceptional large
  entries is unpaid: even a sparse set of edge deletions can increase cap at
  leading scale. Low-rank exceptional structure might help, but coherent
  block-diagonal penalties can perform real midpoint centering, so they
  cannot simply be declared irrelevant. No near-sign stability theorem was
  proved from this idea.
- A scalar/nuclear stationary law is not enough to characterize optimality.
  Bent-convolution Hadamards have a complete Boolean eigenbasis, giving a
  signed covariance of operator norm exactly 1/n despite their suboptimal
  asymptotic cap 1/2. The independent researcher wrote the exact counterexample.
- The root noticed a fixed-nuclear-budget rate improvement via interpolation:
  complex 1-to-infinity norm is one and complex infinity-to-one norm is at
  most four times the real bilinear norm, hence ||A||op<=4sqrt(Q(A)). The
  rank-one Schur envelope then has bias <=2sqrt(Q(A)) Tr|G|. This improves
  the old mask proof for fixed budgets, but the subsequently discovered
  target-contraction operation supersedes the entire bias calculation.
- A self-weave using the seed itself as a row basis looked genuinely dense.
  Root's exact matching witness kills it: pair opposite diagonal signs,
  choose partner columns as fibre spins, and average component reversals.
  Unmatched majority-diagonal fibres use their own columns. The result is a
  half-floor without any orthogonality assumption. The constructive agent
  independently checked this and is preserving the full formula.
- Random orthogonal spin-glass results found in the primary literature
  cover a different coupling law or strict high temperature. They do not
  justify rounding Haar coefficients to full signs without leading error.
  The 2024 Boolean-cube projection/Sidon paper gives order estimates and
  projection limits, not the required complete-unimodular optimum limit:
  https://doi.org/10.1093/imrn/rnae083 . No new implication was imported.
- Root independently reconstructed the adversary's stopped-Brownian Parisi
  lower test from Auffinger--Chen arXiv:1606.05335, Theorem 1 and Corollary 2.
  This DOES have an exact application to iid bridges; the full proof and
  its Bernoulli comparison are being saved separately. It is not a lower
  bound for adversarial full signings in the original optimization.
