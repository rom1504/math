# Protocol: raw negative-tilt overlap of actual minimizing children

Status: **finite falsification audit protocol**.  This experiment
evaluates the exact observable in `L_raw-negative-overlap`; it does not use a
conference/Paley child, a ground-state surrogate, or a proxy row statistic.

## Question

For actual contracted-temperature pressure-minimizing children at
`t=beta/sqrt(N)`, compute

```math
\widehat\rho_N^-(\lambda)
={1\over\lambda mn}\int_{-\lambda}^0
 E_{\widehat\Pi_s}\sum_e r_e(B_{-e})^2\,ds,
\qquad d\widehat\Pi_s\propto e^{sL(B)}dU(B).       \tag{P.1}
```

The complete bridge-cube audit covers balanced splits for `N=4,...,9`,
`beta in {1,2,4}`, both relative child orientations, and

```text
lambda in {0.5, 1, 2, 4, 5.382104195764755}.
```

The last value is the recorded target-reaching inverse-disorder value for the
`N=8,beta=4` case.  An additional complete `3+7` audit at `N=10` is run from
the same program and is clearly labelled as a nonbalanced but comparable
split.

## Child selection

For each `(child order,beta,N)`, enumerate every root-gauged signing, compute
its exact integer absolute-energy histogram, and compare the finitely many
histogram pressures with `mpmath` at 80 decimal digits.  Retain every
signed-permutation/global-sign class in the minimizing histogram set.  The
output records:

- the complete signing and histogram counts;
- the minimizing class counts and representative hashes;
- the high-precision optimum and gap to the next histogram type.

This corrects the loose-tolerance classification in the archived early
tilted-influence experiment.

## Bridge and cavity evaluation

For each selected child-class pair and each orientation, enumerate all
`2^(mn)` bridges.  The Walsh-convolution implementation computes `L(B)` and
checks three masks against direct finite spin sums.

For edge `e`, the exact half-flip cavity identity is

```math
r_e(B_{-e})^2
={\tanh^2((L(B)-L(B^e))/2)\over\tanh^2 t}.         \tag{P.2}
```

The program checks that the right side is invariant within every
`{B,B^e}` pair.  On three deterministic bridge masks it independently
computes the full Gibbs response `m_e` and uses

```math
r_e={m_e-B_e\tanh t\over1-B_e\tanh(t)m_e}          \tag{P.3}
```

to verify the normalized sum of cavity squares against (P.2).

## Tilt integration and recorded quantities

For every declared `lambda`, 32-node Gauss--Legendre quadrature evaluates
(P.1).  The output keeps separate:

- the fair-bridge value at `s=0`;
- the path average (P.1);
- the midpoint at `s=-lambda/2`;
- the endpoint at `s=-lambda`;
- the minimum effective fraction of the exact finite weighted sum;
- the pointwise minimum and maximum over the complete bridge cube.

The bridge sum is exact; pressures, Gibbs responses, and one-dimensional
quadrature are floating evaluations.  No interval/asymptotic certificate is
claimed.

## Reproduction

From the repository root:

```bash
.venv/bin/python \
  extremal_information/experiments/actual_child_negative_overlap_exact.py

.venv/bin/python \
  extremal_information/experiments/actual_child_negative_overlap_exact.py \
  --min-total-order 10 --max-total-order 10 --left-order 3 \
  --output computations/results/actual_child_negative_overlap_exact_n10_3x7.json
```

The primary output is
[`../../computations/results/actual_child_negative_overlap_exact.json`](../../computations/results/actual_child_negative_overlap_exact.json).

The held-out extension through equal splits at `N=14` uses exhaustive child
selection but sampled bridges.  It is separately classified and reproduced
by
[`actual_child_negative_overlap_sample.py`](actual_child_negative_overlap_sample.py).

## Decision rule

- Values bounded away from zero across these orders are only a **finite
  falsifier** to an easy decay proof.  They do not establish an actual
  minimizing sequence with positive limiting overlap.
- Visible decay is only motivation for a theorem and is not fitted to an
  exponent from this short range.
- A RESET requires a uniform optimizer theorem or a scalable certified
  sequence; neither can be inferred from this protocol alone.
