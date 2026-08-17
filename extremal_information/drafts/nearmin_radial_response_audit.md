# Exact finite audit of the near-minimizer radial response roof

Date: 2026-08-17.

Status: **exact finite computation, except for four explicitly bounded
set-cover entries**.  This is a solution-hidden hypothesis audit, not an
asymptotic theorem.  The observables and comparison rule were frozen in
`experiments/nearmin_radial_response_observable_freeze.md` before any response
data were computed.

## 1. Question and population

For a signing `a`, an edge edit `F`, and an augmented cut witness `z`, the
audited roof is

```math
g_z(F)=\langle a^F,z\rangle-Q(a),
\qquad
R_a(F)=\max_z g_z(F)=Q(a^F)-Q(a).
```

Every context in `B_r={F:|F|<=r}` was enumerated for `r=1,2,3`.  At each
order `7<=n<=11`, the population contains three byte-distinct repository
exact minimizers, three cap-`M_n+2` one-edge perturbations, and three fresh
unconditioned random signings.  Thus there are 45 matrices.  The three saved
exact representatives need not be different switching/permutation orbits;
beyond order eight the sample is not an orbit census.

For each matrix and radius the computation records:

* the cardinality `S` of the guaranteed RS shell
  `S_a(Q(a)-M_n+2r)`;
* the number `E` of augmented cut witnesses exposed by some context;
* the minimum number `C_delta` of exposed witnesses whose upper envelope is
  within `delta in {0,2,4}` on every context;
* the mutual information `I(F;Z)` when `F` is uniform and an optimizer `Z`
  is selected uniformly among all ties;
* scalar response entropy and the contextual affine-function metric.

The finite set-cover programs proved 401 of 405 optima.  The four unresolved
entries are only the tolerance-four covers of the three selected order-ten
exact minimizers, each certified in `[10,12]`, and one order-ten/radius-two
tolerance-two cover certified in `[19,22]`.  Every `C_0` value below is exact.

## 2. Radius-three data

The table gives medians across the three matrices in each cell.  The tuple is

```text
S / E / C_0 / C_2 / C_4 ; I(F;Z)/n.
```

For the order-ten exact class, `C_4=10--12` denotes the certified interval.

| `n` | exact minimizer | one-step near-minimizer | uniform random |
|---:|:---|:---|:---|
| 7 | `49 / 49 / 29 / 17 / 10 ; .585` | `43 / 43 / 25 / 11 / 7 ; .522` | `39 / 39 / 24 / 13 / 5 ; .491` |
| 8 | `84 / 84 / 42 / 20 / 12 ; .555` | `80 / 74 / 39 / 22 / 12 ; .532` | `73 / 54 / 28 / 18 / 7 ; .401` |
| 9 | `85 / 85 / 25 / 25 / 11 ; .410` | `132 / 129 / 42 / 42 / 9 ; .383` | `82 / 82 / 24 / 24 / 6 ; .357` |
| 10 | `120 / 120 / 120 / 40 / 10--12 ; .511` | `200 / 130 / 60 / 33 / 21 ; .472` | `186 / 75 / 36 / 16 / 10 ; .319` |
| 11 | `212 / 117 / 68 / 31 / 20 ; .430` | `198 / 119 / 71 / 42 / 17 ; .385` | `175 / 85 / 37 / 20 / 8 ; .238` |

The corresponding exact-cover description rates `log_2(C_0)/n` are:

| `n` | exact | one-step near | random |
|---:|---:|---:|---:|
| 7 | .694 | .663 | .655 |
| 8 | .674 | .661 | .601 |
| 9 | .516 | .599 | .509 |
| 10 | .691 | .591 | .517 |
| 11 | .553 | .559 | .474 |

At radii one and two the same qualitative ordering is noisy at the smallest
orders but stabilizes in the radius-three sample: near-optimal signings are
not visibly simpler than random controls.  At radius three, both the exact
cover rate and optimizer mutual information of exact and one-step-near
matrices are at least the random median at every audited order.  In the
order-ten exact sample every one of the 120 exposed witnesses is required by
the exact roof.

This comparison is especially adverse to the pre-registered compression
hypothesis because the random controls have larger cap excess and hence a
weakly broader *allowed* RS shell.  Nevertheless they generally expose and
need fewer witnesses.

## 3. Scalar simplicity does not imply roof simplicity

The scalar response distribution tells the opposite story.  At radius three,
the median response entropy divided by `log_2|B_3|` is

| `n` | exact | one-step near | random |
|---:|---:|---:|---:|
| 7 | .144 | .163 | .165 |
| 8 | .131 | .140 | .155 |
| 9 | .105 | .100 | .118 |
| 10 | .062 | .116 | .142 |
| 11 | .100 | .112 | .128 |

Thus the exact minimizers have a comparatively simple histogram of scalar cap
changes while using a comparatively rich family of optimizer identities to
realize it.  The sharpest instance is order ten: scalar response entropy is
less than half the random median, but the exact cover needs 120 witnesses
versus 36 for the random median.  An unrooted response histogram therefore
misses precisely the local roof information used by edit composition.

There is a second distinction.  For every audited matrix at radii two and
three, all exposed affine functions form a pairwise `d_r>4` packing.  Hence
their exact metric covering number at tolerance four is `E`.  Yet their
*upper envelope* often has a tolerance-four cover of only 5--21 witnesses.
Approximate extremal sufficiency here comes from global domination in the
upper envelope, not from clustering the affine response functions.  This is a
finite example where response-image metric entropy and approximate
query-sufficient envelope complexity differ sharply.

## 4. Director judgment

### What is falsified at these orders

The data reject the finite hypothesis that near-minimality itself forces a
small exact local edit-response roof.  They also reject using scalar response
entropy as a proxy for rooted optimizer complexity.  Exact and one-step-near
matrices are usually *more* contextually complex than unconditioned random
controls after the frozen normalizations.

### What remains open

Fixed radius has only polynomially many contexts:
`|B_r|=n^{O(r)}`.  Therefore every fixed-`r` exact response cover carries only
`O(r log n)` bits regardless of the signing, and the displayed rates must
eventually vanish if `r` is fixed.  This experiment cannot decide the regime
where `r=r(n)` grows enough to support a physical `n^(3/2)`-scale compiler.
Nor does a small tolerance-four envelope cover establish a reusable
composition congruence.

The surviving precise question is consequently narrower than “is the shell
small?”:

> When the edit radius grows with order, can the approximate envelope cover
> remain sub-landscape in size **and** be refreshed under composition, even
> though its exposed affine functions are not metrically clustered?

This finite audit supplies no positive evidence for that statement.  It does
show that any proof should target envelope domination or a dynamic
congruence, rather than shell cardinality, scalar cap-change histograms, or a
static metric net of exposed functions.

## 5. Reproduction and verification

```bash
.venv/bin/python \
  extremal_information/experiments/nearmin_radial_response_audit.py \
  --solver-seconds 30

.venv/bin/python \
  extremal_information/experiments/verify_nearmin_radial_response_audit.py
```

The second command independently reconstructs every cap, response histogram,
RS shell and exposed set; checks all 401 certified covers; and compares the
closed contextual-metric formula with direct context enumeration on 319
deterministically sampled witness pairs.  Its recorded terminal result is:

```text
PASS: 45 records; 401 certified covers; 319 metric pairs
```

Primary data:
`experiments/nearmin_radial_response_results.json`.
