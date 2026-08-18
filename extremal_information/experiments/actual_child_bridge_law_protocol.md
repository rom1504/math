# Preregistered finite audit: actual child-induced bridge laws

Date frozen: 2026-08-18, before running the audit.

## Question

For the finite-temperature composition identity, do not substitute a
conference or Paley child.  At total order `N=m+n`, target inverse-temperature
parameter `beta`, and raw temperature

```math
t={\beta\over\sqrt N},
```

enumerate the signings which actually minimize the normalized child pressure

```math
p_t(A)=\log\left(2^{-|A|}\sum_x\cosh(tH_A(x))\right).
```

For each minimizing child pair `(A,C)` and each relative child orientation
`epsilon in {+-1}`, define

```math
f_{A,C,\epsilon}(B)=
\log\left(2^{-N}\sum_{x,y}
 \cosh(t[H_A(x)+\epsilon H_C(y)+x^TBy])\right)
```

on the complete bridge cube and audit the canonical negative-disorder law

```math
q_\lambda(B)=
{e^{-\lambda f(B)}\over\mathbb E_Ue^{-\lambda f}}U(B).
```

This is an audit of the actual contracted-temperature child law.  Ground-state
minimizers, conference matrices, Paley matrices, and heuristic children are
not admissible substitutes.

## Frozen finite scope

- total orders: every `N=4,...,9`;
- splits: every `2<=m<=n` with `m+n=N`;
- `beta` grid: `0.25, 0.5, 1, 2, 4`;
- `lambda` grid: `0.25, 0.5, 1, 2, 4`;
- every switching-orbit representative which minimizes each child pressure;
- every resulting child equivalence class/pair and both relative
  orientations.

Order nine is the preregistered computational ceiling because its balanced
bridge has `20` signs and hence `2^20` bridge points.  The next balanced order
has `25` signs and `2^25` points.  The implementation may terminate an order
only for a recorded resource failure; it may not replace it with sampling.

The child signing enumeration and every bridge enumeration are exact finite
enumerations.  Energies and histograms are integer-exact.  Transcendental
pressure comparisons, Gibbs weights, entropies, and divergences are evaluated
numerically; they are not formal certificates.  Child winners and their
nearest distinct pressure competitor will be re-evaluated at high precision,
and the numerical separation will be reported.

## Frozen measurements

For every law, record:

1. `D(q_lambda || U)` and its densities per parent vertex and bridge sign;
2. both the mean pressure gain `E_U f-E_q f` and the negative-moment gain
   `E_U f+lambda^{-1}log E_U exp(-lambda f)`, normalized by `N`;
3. the exact chain-rule effective support `s_*` when feasible, otherwise an
   explicitly labelled upper proxy obtained from row-major, column-major,
   reverse, snake, and marginal-information orderings;
4. row total correlation, individual row marginal `D_2` values, and the
   corresponding column quantities as a transpose check;
5. for every row permutation, the `q`-weighted distribution of
   `D_2(q(R_j | R_{<j}) || U_row)` along the row filtration, including
   weighted means, quantiles, maxima, and threshold masses;
6. for every latent subset of zero or one complete rows, the residual total
   correlation of the other rows conditioned on that latent state and the
   conditional marginal-row `D_2` tail.  This is the preregistered natural
   latent/product decomposition; conditioning on all but one row is excluded
   because it trivially reconstructs the bridge.

The audit will search for, but will not infer asymptotics from, a finite
falsifier to exact row-product structure, bounded-size coordinate support, or
vanishing conditional dependence.  Any claimed falsifier must be a complete
bridge enumeration and must report a numerical robustness margin.

## Companion row-product shadow audit

Added and frozen before running this companion calculation.  For every law
at `N<=8`, and for the balanced `4+5` split at `N=9`, numerically minimize the
exact row-product variational objective

```math
E_{p_1\otimes\cdots\otimes p_m}f
+\lambda^{-1}\sum_iD(p_i\Vert U_n).
```

Use cyclic exact coordinate Gibbs updates, starting from the uniform product,
the product of the exact escort row marginals, and `16` seeded random softmax
starts (`seed=20260818`).  Stop only after a complete sweep changes the
objective by at most `10^-11`, or after `100` sweeps.  Preserve the complete
monotone objective trace and best-response residual.  This is a nonconvex
calculation: the best evaluated value is a heuristic upper bound on
`V_row`, while its improvement over the uniform product is a rigorous
numerical lower bound (up to reported floating error) on the optimal
row-product gain.  It is not a global-optimality certificate.

After the exact target-threshold calculation was inspected, one explicitly
post hoc discriminating run repeated the shadow calculation at the balanced
order-eight `beta=4` target threshold.  It uses the same starts and update
rule but is not part of the preregistered grid; its purpose is to compare
forward total correlation with the directed product-shadow value at a law
that actually reaches the finite target.

## Interpretation discipline

- Variation across actual minimizing child classes is part of the result and
  may not be optimized away after inspection.
- A finite increase of conditional `D_2` is evidence, not a proof of escaping
  asymptotic mass.
- A positive total correlation is a finite falsifier to exact row-product
  structure, not a theorem of extensive dependence.
- No statistic may be defined by selecting a parent optimum or a pressure
  sublevel after seeing the bridge.
