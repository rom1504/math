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
