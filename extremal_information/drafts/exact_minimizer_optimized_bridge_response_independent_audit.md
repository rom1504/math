# Independent audit: exact-minimizer optimized bridge responses

**Verdict: PASS, with the advertised finite-computation qualification.**  The
complete table is reproducible.  An independent HiGHS mixed-integer model,
which uses no first-bridge-row symmetry reduction, certifies every displayed
optimum through order eight.  In particular, it independently reproduces the
order-eight separation `17` versus `19`.  Neither the CP-SAT computation nor
the independent HiGHS computation supplies a standalone proof object, so the
order-three results remain solver-certified exact finite computations rather
than paper proofs.

The audited frozen files have SHA-256 hashes

```text
program 620eb3941aad86c70079a52fc77d3af3bf037aabd4a05cb04c593199219cddc1
result  24047fb3563b51c983eea3e66907393319dae463fdbaaa63f570452a607e84a5
```

## 1. Independent input and reproduction checks

I reran `computations/enumerate_minimizer_orbits.py` from the certified
representatives for every order `3 <= n <= 8`.  The run independently
reproduced the following exhaustive root-gauged counts and
signed-permutation/global-sign class counts:

| `n` | minimizers | classes | class sizes |
|---:|---:|---:|---:|
| 3 | 2 | 1 | `2` |
| 4 | 6 | 1 | `6` |
| 5 | 12 | 1 | `12` |
| 6 | 12 | 1 | `12` |
| 7 | 3240 | 3 | `1680,720,840` |
| 8 | 4200 | 2 | `1680,2520` |

All class rows, canonical orbit hashes, representative matrices, exact
profiles, self-complementarity flags, and class sizes agree with the six
authoritative input files.  The only fields intentionally absent from this
rerun were the optional conference-partition annotations.  I also reran
`computations/certify_m11_m12.py`; its `M_11=17`, `M_12=18` certificate JSON
was byte-identical to the committed file.  This preserves the existing
qualification that the `M_11` lower bound is solver-certified without a
standalone proof object.

I then reran the audited optimized-response program in full.  After ignoring
wall-clock fields, the following all agree exactly with the frozen result:

- input hashes and declared exact values;
- query-orbit counts;
- every response and canonical response partition;
- every exhaustive bridge count and optimal-bridge count for `k <= 2`;
- every saved witness;
- every CP-SAT feasibility status, branch count, and conflict count; and
- every switching, permutation, and global-sign transport count.

## 2. Independent optimization model

As a solver-independent formulation check, I rebuilt every instance directly
for SciPy 1.13.1/HiGHS, without importing the experiment and without its
first-row orbit restriction.  Write each bridge sign as `b_ij=1-2z_ij` with
`z_ij` binary.  For every `x` with `x_0=1` and every unrestricted query spin
`y`, I imposed

```math
-T\le
H_A(x)+H_C(y)
+\sum_{i,j}x_i y_j
-2\sum_{i,j}x_i y_j z_{ij}
\le T.                                                     \tag{IA.1}
```

Thus this model contains every projective parent spin and all `2^(nk)`
bridges, but none of the source model's bridge symmetry reduction.  For all
nine exact-minimizer classes and all four displayed queries, HiGHS found the
reported target `T` feasible and certified the parity-predecessor `T-2`
infeasible.  This independently verifies every entry in the response table,
not only its two nonconstant partitions.

For the decisive order-eight positive- and negative-triangle instances the
independent results are:

| class | triangle product | target 17 | target 19 |
|---:|:---:|:---:|:---:|
| 0 | `+1` | feasible | not needed |
| 0 | `-1` | feasible | not needed |
| 1 | `+1` | infeasible | feasible |
| 1 | `-1` | infeasible | feasible |

Independent exhaustive spin evaluation of the four HiGHS witnesses gives
caps `17,17,19,19`, respectively.  The class-1 target-17 infeasibility runs
took about 23 seconds each with no symmetry restriction.  HiGHS, like CP-SAT,
reported a final infeasibility status but emitted no separately checkable
proof certificate.

## 3. Formulation and quotient audit

The source formulation is exact.

1. Fixing the first old-child spin to `+1` removes precisely the global spin
   reversal.  Keeping every query spin unrestricted therefore enumerates all
   projective parent contexts exactly once.
2. Each parent energy is an integer sum of
   `binom(n+k,2)` signs, so all possible caps have that parity.  Starting at
   the certified global lower bound and advancing by two cannot skip an
   attainable objective.
3. If `U` is a signed permutation with `U^T C U=C`, then `B -> BU` is a
   bijection of the bridge fibre preserving the parent cap.  Also `B -> -B`
   is conjugation by switching every old-child vertex.  Hence restricting
   the first bridge row to one representative of these genuine orbits is
   valid.  The independent model above confirms that it removes no optimum.
4. Hollow sign queries have one switching/permutation class for `k=1,2` and
   two for `k=3`, distinguished by the triangle product.  Direct exhaustive
   orbit enumeration reproduces `1,1,2`.
5. The authoritative minimizer files identify `A` with `-A`, while

   ```math
   F_C(-A)=F_{-C}(A).
   ```

   Negating a three-vertex query reverses its triangle product.  Sorting the
   two order-three coordinates is therefore exactly the residual global-sign
   quotient; comparing the oriented pair would not be class-invariant.

The frozen signatures and partitions are consequently correct.  At order
seven, class 1 is separated from classes 0 and 2 already by `k=1`.  At order
eight, both classes have responses `(12,15)` for `k=1,2`, but their
order-three pairs are `(17,17)` and `(19,19)`.

## 4. Archive comparison and scope

The order-seven separation recovers the one-row insertion obstruction in
`artifacts/scale_transfer_profile_no_go.md`: two minimizers with the same
radial signed-extrema profile have optimized extension caps 12 and 10.  The
archive also distinguishes the two order-eight classes under the *fixed*
universal-double construction, with caps 40 and 32.  That fixed composition
does not imply the present result, because here every bridge is optimized
separately and child switching and global-sign gauge are removed.  I found no
earlier archived table that tests these optimized query fibres and locates
the first separation at query order three.

The result has a deliberately narrow scope:

- it is a fixed-order separation by two, not a macroscopic gap;
- two classes do not yield a packing-rate lower bound;
- there is no evidence that query order three remains discriminating in a
  scalable family;
- evaluating `F_C(A)` still carries out full bridge optimization; and
- it provides neither a compressed composition law nor an asymptotic
  convergence statement.

What is rigorously supported is the finite falsifier: a gauge-invariant state
that identifies the two order-eight minimizer classes cannot predict every
optimized three-vertex continuation with additive error strictly below one.
Optimizing the bridge erases pure switching labels, but does not erase all
information distinguishing inequivalent exact minimizers.
