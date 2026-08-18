# Finite actual-child audit of the signed sector-mode product trial

Classification: **80-digit optimizing-child selection, complete finite
bridge enumeration, and floating eigensystem/transcendental evaluation**.
This is reproducible numerical evidence, not an interval certificate and not
an asymptotic theorem.

The experiment
[`actual_child_signed_sector_mode_falsifier.py`](actual_child_signed_sector_mode_falsifier.py)
tests the one-parameter product path in Theorem SM.2 of
[`../drafts/actual_child_signed_sector_mode_product_certificate.md`](../drafts/actual_child_signed_sector_mode_product_certificate.md).
It covers every balanced order `N=4,...,9`, `beta=1,2,4`, and `lambda=1`.

## Protocol

For each order and temperature, the driver first sets `mp.mp.dps=80` and
selects the genuine contracted-temperature minimizing class.  It then:

1. chooses the bias-canceling row direction and orientation of Theorem
   37.32;
2. forms the signed cross-row quadratic matrix `M` from the exact zero-bridge
   child law;
3. removes the arbitrary basis choice in a multiple bottom eigenspace by
   projecting the first available coordinate vector into that eigenspace;
4. rounds each nonzero row block by the odd hyperplane bit
   `phi_i(b)=sgn_*(<v_i,b>)`;
5. computes `w_i=E_r[B phi_i(B)]` under the exact canonical row law; and
6. whenever

   ```math
   R_{\rm quad}:={-t^2w^{\mathsf T}Mw\over k}>1,
   ```

   evaluates both the exact quadratic product gain and the exact
   physical-interaction gain at the prescribed SM.24 amplitude.

Thus no product coordinate descent, external-field grid, or bridge surrogate
enters the calculation.  Every physical bridge in the declared finite cube
is enumerated.

## Result: quadratic instability is canceled at physical amplitude

Four genuine minimizing-child cases cross the rounded quadratic threshold:

| `N` | `beta` | active rows `k` | `R_quad` | SM.24 `a` | quadratic gain | exact actual-`h` gain | remainder contribution |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | 4 | 2 | `3.717397` | `.612372` | `.792592` | `-.071170` | `-.863762` |
| 5 | 4 | 2 | `2.977547` | `.612372` | `.572196` | `-.192136` | `-.764332` |
| 6 | 4 | 3 | `1.105939` | `.268035` | `.009683` | `-.064272` | `-.073955` |
| 7 | 4 | 2 | `2.274495` | `.612372` | `.362763` | `-.278556` | `-.641319` |

In every case the signed quadratic mode gives a strictly positive coherent
product gain, including the explicit positive lower bound in SM.25, while
the same declared product law has **negative** gain against the full actual
bridge interaction.  The physical order-at-least-four remainder reverses
the sign.

This is a sharp finite actual-law falsifier of the following tempting
promotion rule:

```text
an unstable sector quadratic mode by itself implies coherent physical
product retuning along its rounded principal direction.
```

It does not falsify Theorem SM.2, whose exact statement subtracts the
physical remainder, and it does not rule out a different coherent direction
or a full optimal-product gain.  It also does not prove that the cancellation
persists at fixed `beta` as `N` tends to infinity.  Its exact lesson is that
the remainder premise in SM.22--SM.25 cannot be deleted using only the
present finite actual-minimizer evidence.

The other fourteen cases have `R_quad<=1` for the canonical eigenspace
projection and therefore do not trigger this trial.  In particular, a large
Frobenius sector--Gram mass is not enough: the certificate depends on a
signed negative mode, its row rounding, and the physical remainder.

## Internal checks and scope

- The corrected selector finds exactly one minimizing class at every child
  order used here.  All nontrivial gaps to the next histogram are stored in
  the output.
- Both canonical orientations are not enumerated: the experiment deliberately
  uses the bias-canceling orientation requested by the frontier.
- Canonical row log-likelihoods agree across rows to floating error before
  the iid product is formed.
- Every odd rounded feature has canonical mean zero to floating error.
- Every declared product trial normalizes to error below `2e-11`.
- Eigenvalue multiplicities, eigengaps, feature tie counts, canonical `J`,
  quadratic gains, actual gains, and remainder contributions are all stored
  in the machine-readable output
  [`../../computations/results/actual_child_signed_sector_mode_falsifier.json`](../../computations/results/actual_child_signed_sector_mode_falsifier.json).

Because eigenvectors and transcendental weights are evaluated in floating
arithmetic, the displayed strict inequalities are numerical.  Their margins
are large except for the positive quadratic gain at `N=6`; an interval
certificate would be required before promoting any decimal here to a finite
theorem.
