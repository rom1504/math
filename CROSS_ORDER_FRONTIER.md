# Cross-order frontier

## Permanent campaign target

**Top-level SML (fixed):** Prove a genuinely sublinear composition defect
for actual optimizing children.

The primary quantity is the minimized scaled-temperature pressure

```math
P_n(\beta)=\min_A\log\!\left(2^{-n}\sum_x
\cosh\!\left({\beta H_A(x)\over\sqrt n}\right)\right).
```

## Checkpoint 0 — inherited baseline

- **Best proved bound:** for `N=m+n`, contracted-temperature minimizing
  children and every split,

  ```math
  P_N(\beta)\le P_m(\beta)+P_n(\beta)
  +mn\log\cosh(\beta/\sqrt N)
  \le P_m(\beta)+P_n(\beta)+{\beta^2mn\over2N}.
  ```

- **Previous/current exponent:** `1 / 1`; the comparable-split defect is
  `Theta_beta(N)`.
- **Assumptions:** fixed `beta>0`; exact contracted-temperature pressure
  minimizers are used in the annealed bridge identity.
- **Actual optimizing children:** yes.
- **Order coverage:** all splits and all orders, but the estimate is not
  sublinear on comparable splits.
- **Top-level SML:** **UNCHANGED**.

Every subsequent auxiliary statement is counted only together with an
explicit implication `P => E_N <= [quantitative bound]` (or the analogous
bound for another quantity already shown to imply convergence).

## Checkpoint 1 — fractional-cardinality coefficient improvement

- **Best proved bound:** for every `0<q<=1`, all orders, and exact
  own-scale pressure-minimizing children,

  ```math
  E_{m,n}(\beta)
  \le {mn\log\cosh(q\beta/\sqrt N)\over q}
   +{(N-1)(1-q)\log2\over q}-\Delta_A-\Delta_D.
  ```

  At an equal split and `beta>=sqrt(8log2)`, the optimized positive
  coefficient is `beta sqrt(log2/2)-log2`; it is strictly below the
  annealed `beta^2/8` when the inequality on `beta` is strict.
- **Previous/current exponent:** `1 / 1`; the coefficient improves, but the
  comparable-split defect is still `Theta_beta(N)`.
- **Assumptions:** fixed `beta>0`; exact child pressure minimizers;
  rank-one cut-word support and the exact radial payments
  `Delta_A,Delta_D`.
- **Actual optimizing children:** yes.
- **Order coverage:** every split and every order; the coefficient statement
  is asymptotic at comparable equal splits.
- **Top-level SML:** **UNCHANGED**.

## Checkpoint 2 — switching-orbit reduction and scalar floors

- **Best proved bound:** in addition to the fractional-cardinality bound,
  every pair of exact own-scale minimizing children satisfies

  ```math
  E_{m,n}(\beta)
  \le \Psi_{m,n}(\beta/\sqrt N)-\Delta_A-\Delta_D
      +\log(1-|u_Au_D|),
  ```

  where `Psi` is the minimum pure bipartite pressure.  The cruder
  orientation-averaged version drops the logarithm.  At equal splits that
  orientation-averaged certificate has a positive linear floor for every
  fixed `beta>4.0515964866...`; the fractional-cardinality certificate,
  even with its exact orientation term, also has a positive linear floor
  for `beta>7.814070149...`.
- **Previous/current exponent:** `1 / 1`; taking
  `Psi<=mn log cosh(beta/sqrt N)` recovers an `O_beta(N)` bound, and no
  proved estimate makes either exact right side sublinear.
- **Assumptions:** fixed `beta>0`; exact own-scale child pressure minimizers;
  the floor statements use only the rigorous all-order cap upper frontier
  and are lower bounds on certificate values, not on the true defect.
- **Actual optimizing children:** yes.
- **Order coverage:** the switching inequality holds for every split and
  every order; its positive floor is proved for all sufficiently large
  balanced orders in the stated temperature range.
- **Top-level SML:** **UNCHANGED**.

## Checkpoint 3 — joint quotients and separable-method barriers

- **Best proved bound:** the actual-child recurrence is still
  `E_(m,n)(beta)=O_beta(N)`.  Two scalable barriers now delimit possible
  improvements: every deterministic positive Hölder allocation that pays
  the two children and bridge separately has a positive linear certificate
  on balanced splits for `beta>2.554944594406...`; and an exact labelled
  switch-convolution quotient of dimension `ell` can gain at most

  ```math
  2\beta {\max(m,n)\over\sqrt N}\,\ell.
  ```

  Hence cancelling a residual `cN+o(N)` term by this genuinely joint
  quotient requires `ell=Omega(sqrt(N))`.
- **Previous/current exponent:** `1 / 1`; no actual-child exponent
  improvement has been proved.
- **Assumptions:** fixed `beta>0`; exact own-scale pressure minimizers.  The
  Hölder floor concerns deterministic three-factor positive separation.
  The Fourier lower bound is conditional only on a positive linear
  pre-cancellation term and leaves `sqrt(N)<=ell<N-1` open.
- **Actual optimizing children:** yes; both method barriers are uniform in
  the selected child minimizers.  The separate marginal-entropy collision
  is an abstract switching-group no-go, not an actual-child example.
- **Order coverage:** all orders for the exact inequalities; comparable
  balanced large orders for the stated linear floors.
- **Top-level SML:** **UNCHANGED**.

## Checkpoint 4 — Hadamard replicas isolate the core

- **Best proved bound:** the unconditional actual-child bound remains
  `E_(m,n)(beta)=O_beta(N)`.  For every admissible symmetric Hadamard order
  `k` and every `r`, however, annealed completion gives the exact
  multiplicative-order inequality

  ```math
  P_{kr}(\beta)-kP_r(\beta)
  \le \mathcal D^\star_{k,r}(\beta)+{\beta^2(k-1)\over4},
  ```

  where `D*` minimizes the exact rotated-cube core defect over the actual
  order-`r` minimizing fibre and over admissible outer matrices.  Thus the
  entire sign-completion payment is `O_beta(k)=o(kr)` when `r->infinity`;
  only the core can retain a linear defect.  At zero temperature the proved
  arrow is

  ```math
  b_{kr}-kb_r
  \le {2\over3}(k^{3/2}M_r)^{-1/3}(R_cap+C)_{+},
  ```

  with `C=sqrt(k(k-1)r(kr+2)log2)`.  Hence
  `R_cap=O((kr)^(3/2-delta))` implies
  `E=O((kr)^(1-delta)+kr/sqrt(r))`.
- **Previous/current exponent:** `1 / 1`; the completion component improved
  from `1/2` to `0` for fixed `k`, but no bound on the actual minimizing
  core improves the total exponent.
- **Assumptions:** fixed `beta`; actual pressure or cap minimizer at order
  `r`; a symmetric Hadamard outer matrix.  Exact data prove
  `P_16(beta)<=4P_4(beta)+3beta^2/4`, but growing replica order is
  nonmonotone and can have a linear core floor for a fixed actual child.
- **Actual optimizing children:** yes; the child is selected from the exact
  minimizing fibre.  The growing-`k` floor uses an actual all-temperature
  order-four minimizer but keeps its child order fixed.
- **Order coverage:** every child order `r` and every available symmetric
  Hadamard multiplier `k`; this is not an arbitrary comparable split, and
  the core estimate required for sublinearity is open.
- **Top-level SML:** **UNCHANGED**.

## Checkpoint 5 — full marginal accounting and bounded replica selection fail

- **Best proved bound:** the best unconditional comparable-order inequality
  remains `E_(m,n)(beta)=O_beta(N)`, with exponent one.  Two stronger
  certificates are now quantitatively excluded.  An unbiased iid bridge,
  even with the complete child energy tables and exact one-constraint
  binomial probabilities, cannot certify a parent energy radius below

  ```math
  \sqrt{2mn(m+n)\log2}-o(N^{3/2}).
  ```

  At an equal split the `1/2` cap upper frontier forces this method to pay
  at least `0.0724620948...N` in the `b=M^(2/3)` defect.  Separately, the
  unique exact order-five cap minimizer satisfies
  `M(H_2 tensor A)=12>8sqrt(2)` and
  `M(H_4 tensor A)=40>32`, so selecting between the two bounded canonical
  replica orders does not make the core favorable.
- **Previous/current exponent:** `1 / 1`; neither result improves the
  unconditional exponent, while both prove a fixed linear floor for their
  respective certificate classes.
- **Assumptions:** iid-bridge floor: independent unbiased bridge and
  marginal union/first-moment accounting, with arbitrary exact
  energy-dependent thresholds.  Replica falsifier: exact cap minimizers,
  Sylvester outer orders two and four; its dyadic amplification produces
  one bad symmetric outer at every dyadic multiplier.
- **Actual optimizing children:** yes.  The iid result is uniform over all
  children and hence over optimizers; the replica result uses the exhaustive
  unique order-five minimizing class.  Neither is a lower bound on the true
  optimized parent defect.
- **Order coverage:** the iid obstruction covers every comparable split and
  is asymptotically sharp for its method class.  The replica result is a
  finite child-order falsifier with a fixed-seed amplification, not an
  all-large-order obstruction or a recurrence.
- **Top-level SML:** **UNCHANGED**.

## Checkpoint 6 — Gaussian soft minima do not yet calibrate the endpoint

- **Best proved bound:** the unconditional actual-child inequality remains
  `E_(m,n)(beta)=O_beta(N)`.  For every fixed disorder exponent `s`, the
  physical Rademacher bridge moment has an all-split Gaussian replacement
  error `O_(beta,s)(sqrt(N))`, and hence

  ```math
  P_N(\beta)\le \Phi_s^{\rm Gau}+O_{\beta,s}(\sqrt N).
  ```

  This does not improve the defect because no proved upper calibration puts
  `Phi_s^Gau` within `o(N)` of `P_m+P_n`.  After resolving one ordinary
  exponential sector, the natural child-variance-to-cross-variance
  interpolation moves in the wrong direction for every `s<1`; retaining
  the physical `cosh` sector makes its signed square indefinite and leaves
  that route open.
- **Previous/current exponent:** `1 / 1`; the `1/2` universality exponent
  belongs only to the replacement error, not to the total cross-order
  defect.
- **Assumptions:** fixed `beta` and fixed `s`; actual own-scale minimizing
  children for the bridge reduction.  The sector-resolved sign theorem is
  uniform in the deterministic base measure, while its exact order-two
  optimizer witness is finite and does not provide an asymptotic floor.
- **Actual optimizing children:** yes for the replacement and direct
  `P_N` arrow; no actual-child theorem controls the remaining Gaussian
  endpoint.  A one-spin joint-sector example proves that sector resolution
  cannot be restored for free.
- **Order coverage:** the replacement covers every order and split.  It
  yields no sublinear inequality on comparable large orders.
- **Top-level SML:** **UNCHANGED**.  This checkpoint is a **STRIKE**, not a
  criterion-four reset, because the same covariance architecture with the
  actual joint-sector sign remains live.

## Checkpoint 7 — campaign verdict after outward review

- **Best proved bound:** the actual-child inequality remains
  `E_(m,n)(beta)=O_beta(N)`.  No `o(N)` comparable-order defect was proved.
  Separately, a log-periodic hereditary countermodel has
  `beta D_(N,floor(N/2))-S_(N,floor(N/2))
  =(0.02 beta+o(1))N^2`; this persists under every
  `O(n^gamma)`, `gamma<3/2`, scalar perturbation, including subleading
  genuine cut-cap data and an affine `2^n`-cut-query representation.
- **Previous/current exponent:** `1 / 1` for the actual defect.
- **Assumptions:** fixed `beta`; exact own-scale pressure minimizers for the
  actual bound.  The no-go concerns the explicitly listed scalar
  heredity/extension/entropy/canonical resources, not the unit-leading
  joint incidence of the actual cut characters.
- **Actual optimizing children:** yes for the unchanged `O_beta(N)` bound;
  no for the proof-class countermodel.
- **Order coverage:** every order and split for the actual linear bound;
  an infinite sequence of equal-ratio pairs for the scalable no-go.
- **Top-level SML:** **UNCHANGED**.  The struck actual-child implementation
  remains frozen.  Campaign accounting is **RESET under criterion 4 only**:
  the scalar canonical-pressure proof class is quantitatively insufficient
  and any continuation must use leading joint cut incidence.
