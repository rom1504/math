# Finite audit of affine closure in near-top shells

Date: 2026-08-17.

Status: reproducible finite evidence only.  The computation is
[`experiments/nearmin_shell_parallelogram_audit.py`](../experiments/nearmin_shell_parallelogram_audit.py)
with frozen output
[`experiments/nearmin_shell_parallelogram_results.json`](../experiments/nearmin_shell_parallelogram_results.json).

## Frozen observable

For every available exact minimizer, one-edge near-minimizer, and sampled
uniform signing through order 14, enumerate projective spins and their
absolute energy deficits.  Record the smallest shell width containing four
distinct projective spins whose coordinatewise product is the all-positive
spin.  This is the first nontrivial affine parallelogram in the projective
spin cube.  The statistic was chosen before inspecting its values.

## Result

At orders 9--14, almost every stored exact minimizer already has such a
parallelogram in its exact active set; the only exception among these
records has threshold 2.  Uniform random controls usually first acquire one
at deficit 4, with observed thresholds between 0 and 10 (one order-8 sample
has 14).  One-step near-minimizers are mixed.  Thus exact activity shows a
finite separation, but the normalized thresholds of both target and random
classes are already small and decrease with order.

The experiment supplied a useful conjecture--near-top affine structure--but
not a near-minimizer-specific asymptotic invariant.  The subsequent
low-local-field theorem proves a much stronger statement for **every**
bounded-cap signing: an affine cube of dimension `Theta(sqrt n)` lies within
`O(n)` of the cap.  Therefore the existence of one `o(n^(3/2))`-shell
parallelogram cannot distinguish near-minimizers from generic bounded-cap
controls.

## Verdict

* Exact-shell parallelograms are a real finite feature of the stored
  minimizers.
* Mere affine closure in a vanishing normalized shell is universally forced
  by a ground state's low local fields, so it is **falsified as a selective
  near-minimizer law**.
* Any surviving candidate must concern how an affine shell responds to
  declared contexts or composes, not whether it exists.

No asymptotic inference is made from the finite table.
