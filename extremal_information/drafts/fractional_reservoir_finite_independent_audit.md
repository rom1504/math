# Independent audit: finite shell-wide fractional reservoirs

**Verdict:** **PASS for the stated finite numerical and exact-certificate
claims.**  The experiment does not establish an asymptotic reservoir theorem,
and its corpus is not an ensemble sample.

## Frozen objects checked

| object | SHA-256 |
|---|---|
| input corpus | `2c086cf7523ead804942948e800c6231eac33d954e049b5aa113c9fb0cca47a5` |
| protocol | `0210db40945b61fd4caf50d602578b94c5c0e9c83a7f9590cd0e69bd9cff41e6` |
| solver script | `a0b3db41b8169bbc249dcfe5bf48b82f15002b398f2fca457b229dfbcf75149d` |
| stored result | `de1d74c8716fc99413343c3db69882b06f149dba7e8d31bc4ddf8115fa6cdcc8` |

I reran the frozen script with its output redirected under
`/home/math/quadra/tmp/fractional_audit_verify/`.  The recomputed JSON was
byte-for-byte identical to the stored result.

## Independent reconstruction

I separately reconstructed the corpus from the input JSON rather than trusting
the stored inventory.  The stated priority, SHA ordering, and global
deduplication produce 66 candidates, two duplicates, and 64 distinct matrices:

| stratum | count |
|---|---:|
| exact | 21 |
| one-step-near | 9 |
| independently generated heuristic low-cap | 6 |
| uniform-random low-cap | 12 |
| structured control | 16 |

The two removed duplicates are a size-13 `control_extremes` matrix already
selected from the cyclic controls and a size-3 control already selected as an
exact representative.  All 64 selected arrays are symmetric, hollow, have
off-diagonal entries in `{-1,1}`, and have distinct serialized-matrix hashes.
Independent enumeration of all `2^(n-1)` projective spin configurations
reproduced every recorded cap and cap delta; all 21 exact-labelled matrices
have cap delta zero.

I then independently rebuilt each oriented active and deficit-2 response
matrix.  All 128 shell sizes, response extrema, common-correct counts, column
pattern counts, and exact uniform-certificate integer residuals agree with the
stored records.

## LP and dual checks

For response matrix `R`, threshold `m`, and `E` edges, the solved primal is

```math
\min \mathbf 1^T w
\quad\text{such that}\quad
-Rw\le -m\mathbf 1,
\qquad 0\le w\le\mathbf 1.
```

Under SciPy/HiGHS's marginal convention, the dual objective is exactly

```math
(-m\mathbf 1)^T y+\mathbf 1^T u,
```

where `y` is `ineqlin.marginals` and `u` is `upper.marginals`; the lower-bound
term vanishes because the lower bound is zero.  Thus the reconstruction used
by the script has the correct signs and includes every nonzero bound term.

I re-solved all 128 dual-simplex problems and checked stationarity and
complementary slackness, which the stored script did not explicitly record.
The worst residuals were:

| check | maximum residual |
|---|---:|
| stationarity | `3.00e-13` |
| inequality complementary slackness | `4.50e-12` |
| lower-bound complementary slackness | `0` |
| upper-bound complementary slackness | `0` |
| primal-dual objective gap | `4.55e-13` |
| recomputed versus stored objective | `2.84e-14` |

The stored cross-method and primal checks also reproduce: maximum
dual-simplex/interior-point discrepancy `4.15e-12`, constraint violation
`3.07e-11` (the independently read value is `3.0603e-11` before rounding),
box violation `5.77e-12`, and duality gap `4.55e-13`.

## Exact full-mass certificate

Let `K` be the number of shell rows and let
`c_e=sum_z R_{z,e}`.  Every feasible `w` obeys

```math
m
\le {1\over K}\sum_z\sum_e R_{z,e}w_e
=\sum_e {c_e\over K}w_e.
```

If `E c_e <= K m` for every edge, nonnegativity of `w` gives

```math
m\le {m\over E}\sum_e w_e,
```

so `sum_e w_e >= E`.  The box constraints give the reverse inequality, and
full edge weight is feasible because every oriented shell row has total
response at least `m`.  Therefore the only feasible weight vector is
`w=1`, not merely one numerical optimum.  This verifies both the logic and
the integer implementation of (FA.1).

Independent exact integer evaluation finds 18 certificates, all on active
shells.  Seven are in the exact stratum: orders 3, 5, 6, both selected
order-13 representatives, and both selected order-14 representatives.  There
are no deficit-2 uniform-average certificates.  Every stored certificate flag
agrees with the independent calculation.

## Recomputed headline summaries

| shell | `C_inst` range | median | mean | numerical full mass | integral reported optimum |
|---|---:|---:|---:|---:|---:|
| active | `1`--`4.333333` | `2.333333` | `2.323287` | 36/64 | 62/64 |
| deficit 2 | `1`--`5` | `3` | `2.942564` | 31/64 | 47/64 |

For the 21 exact records, the corresponding means are `2.814544` and
`3.568767`, with full-mass counts 13/21 and 15/21 and no-literal-reservoir
counts 17/21 and 21/21.  Spearman correlations of shell size with `C_inst`
are `0.934920` and `0.858160`.  The within-order cap-delta comparison also
reproduces exactly: active means favor delta zero in 9/10 orders with one tie
and mean difference `0.719353`; deficit-2 means do so in 10/10 orders with
mean difference `1.124631`.

The large-order statements check exactly against the frozen records.  Both
selected order-14 exact matrices have active `(K,m,W)=(156,21,91)` and
deficit-2 `(520,19,91)`.  Both order-13 exact matrices have active
`(78,20,78)` and deficit-2 `(78,18,70.2)`; in the latter shell all included
responses actually have total 20, so the uniform weight `0.9` is feasible.

## Scope cautions (not failures)

1. The embedded hashes certify which files were used, but the present snapshot
   cannot independently prove the temporal assertion that the protocol was
   written before any LP was inspected.
2. The corpus is deliberately low-cap and highly selected.  SHA-first choice
   is deterministic, not an ensemble-sampling theorem; exact orders 9--14 use
   at most two repository representatives and are not claimed to be distinct
   switching-permutation orbits.  The reported correlations and paired
   contrasts therefore remain descriptive and confounded by order, shell size,
   and source composition.
3. Dual simplex and interior point are two HiGHS algorithms sharing the same
   implementation and presolve stack.  Their agreement is strong numerical
   validation, not an exact certificate.  Only the 18 integer
   uniform-average cases have the exact full-mass proof recorded here.
4. The script currently selects repository exact representatives for every
   input order `n>=9`, whereas the protocol says orders 9--14.  The frozen
   input has maximum order 14, so this causes no discrepancy in this run, but
   a future reuse should impose `9 <= n <= 14` or revise the protocol.
5. Full-shell LPs have anchor counts as large as 520.  They do not contradict
   a fixed-anchor theorem whose constant may depend on the anchor count, and
   no finite trend here supplies a scalable asymptotic family.

Subject to these scope limits, the draft's central research judgment is
supported: whole-shell fractionalization frequently retains the full physical
edge interface, while sequential finite-anchor fractional reservoirs are not
falsified by this audit.
