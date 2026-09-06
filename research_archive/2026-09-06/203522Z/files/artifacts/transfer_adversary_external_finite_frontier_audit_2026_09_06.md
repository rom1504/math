# Attributed external finite frontier: independently verified witnesses and certificate boundaries

Date: 2026-09-06. External author: Rob Sneiderman, GitHub user Robby955.
Repository: [mo-413935-research](https://github.com/Robby955/mo-413935-research),
pinned commit `61ac2268b6f9234c3ea268e7a7072ac72141f36b`.
This note does not import that repository's asymptotic bounds or
composition claims into the present project's frontier.

## Independently verified mathematical data

Our checker
`computations/transfer_adversary_external_order15_16_witness_verify_2026_09_06.py`
initially extracted only literal witness data from upstream syntax; it
did not import or run the upstream energy evaluators. A separately written
integer Gray-cube scan checks every projective spin, exact first and
second moments, and periodic direct-energy/field reconstructions.
The matrices and complete projective energy histograms are preserved in
`computations/results/transfer_adversary_external_order15_16_witness_verify_2026_09_06.json`.

The DEFAULT checker mode, also selected by `--replay-preserved`, now
re-enumerates those preserved matrices, histograms and deletion caps
without reading the excluded external cache. Only explicit
`--extract-upstream` regenerates the data from upstream syntax, after
checking its source hashes against the durable dependency manifest.
The preserved-data replay therefore remains standalone after snapshots
exclude the downloaded sources. The separately named final-level
assessment below intentionally retains a documented download dependency.

| Witness | Independently verified Q | Positive / negative projective grounds | Single-deletion caps |
|---|---:|---:|---|
| `/cases/0/matrix`, order 15 | 27 | 66 / 0 | 21 twice, 27 thirteen times |
| `/cases/1/matrix`, order 16 | 30 | 19 / 19 | 29 at all sixteen deletions |

The order-16 matrix digest and its entire full-cube histogram also match
the upstream pins. These are unconditional UPPER witnesses. This replay
does not independently establish `M_15=27` or `M_16>=28`.

The seed agent separately checked the matrices' balanced partitions:
104 of 6435 order-15 splits have child caps 9 and 10, and 3 of 6435
order-16 splits have child caps 10 and 10. These match the known smaller
minimum values, but do not establish optimality of either parent.
See `transfer_seed_existential_partition_screen_2026_09_06.md` for the
independent partition source and precise scope.

## What the order-15 lower certificate would require

The [pinned upstream certifier](https://github.com/Robby955/mo-413935-research/blob/61ac2268b6f9234c3ea268e7a7072ac72141f36b/verification/research_order15_certify.py)
separates an explicit cap-27 witness from its lower certificate. The
lower chain enumerates root-normalized residual graphs at order twelve,
then extends twice and canonically deduplicates, producing 1,313,164
order-fourteen records. The final task rules out every cap-at-most-25
one-row extension of those records. Odd order-fifteen energy parity
then yields 27, CONDITIONAL on completeness of that catalogue.

The extension identity is exact:
`Q(extension by v)=max_x (|H_A(x)|+|v dot x|)`.
Both x and v can be represented modulo global sign. Rejecting v against
any subset of x is a safe necessary-condition filter; survivors require
a full-cube test. The final checker implements this correctly. Its
binary32 products contain only small exact integers; the relevant
products and partial sums stay far below `2^24`, so this use does not
introduce a floating-point threshold ambiguity on the stated CPU path.

We independently verified the downloaded catalogue's 19,697,460-byte
decompressed payload, record count, and SHA-256:

```text
a0f0c0d2e21f382bb73923f014b3936cb60f59ac0a7b454324fdcb5c142d4ecc
```

A digest establishes the identity of the data, not its completeness.
The author explicitly retains that boundary. Full reconstruction of the
lower levels is documented as approximately one day and 15 GB peak disk
in the [pinned sizing note](https://github.com/Robby955/mo-413935-research/blob/61ac2268b6f9234c3ea268e7a7072ac72141f36b/verification/order15_tower/F15_chain.txt).
It was not launched here. The tower driver also contains ordinary shell
pipelines with `set -e` but no `pipefail`; a documentary receipt alone
is not an independent certificate that every producer completed. This
is an audit boundary, not evidence that the recorded mathematical
result is false.

## Bounded replay assessment, not a full lower replay

Our separate assessment driver is
`computations/transfer_adversary_external_order15_replay_assessment_2026_09_06.py`.
It verifies the pinned source hash before importing the already-inspected
external final-level functions. It does not call the upstream main or
run any tower shell script. Its result JSON has the same basename under
`computations/results/`.

Exactly 256 evenly spaced catalogue records were tested: every one had
cap 25, and no cap-at-most-25 extension survived. Sixteen of their caps
were independently rechecked by our integer Gray-cube implementation.
With standard BLAS thread controls fixed to one, the timed sample took
1.7675 seconds. Linear extrapolation is about 151 minutes for one worker
or an IDEAL 37.8 minutes for four. These are estimates, not completion
guarantees; sample variation and shared-machine contention matter.
No full 1,313,164-record replay was launched during this bounded task.
Even completing that final replay would still leave the lower-level
catalogue-completeness dependency.

One preliminary benchmark attempt stopped immediately because
`threadpoolctl` was unavailable. No package was installed. The corrected
driver sets the standard BLAS thread environment controls before NumPy
is imported; its complete bounded rerun passed. No failed attempt wrote
a result file.

## Order sixteen and the restricted order-fifteen census

The [pinned reduction program](https://github.com/Robby955/mo-413935-research/blob/61ac2268b6f9234c3ea268e7a7072ac72141f36b/verification/research_order16_reduction.py)
does NOT claim `M_16=30`. Its asserted reduction depends on the complete
order-fourteen catalogue and `M_15=27`. Within that framework a
hypothetical cap-28 order-sixteen signing must have every order-fourteen
pair-deletion at cap 27. The committed 30 rooted order-fifteen records
cover only minimizers having SOME order-fourteen deletion of cap at most
25. They are not a complete list of all order-fifteen minimizers.
Neither that census's completeness nor the reduction's full scan was
replayed here. The pinned-frame lower bound 30 is a restricted-family
statement and must not be treated as a global lower bound at order 16.

## Provenance, license, and preservation

The complete upstream tree contains no LICENSE file, and the downloaded
GitHub repository metadata reports `license: null`. No reuse or
redistribution grant has been inferred. The 25 downloaded upstream
files remain local under
`external_research/2026-09-06/robby955_mo_413935_61ac2268b6f9/`.
They are reproducible external dependencies and are NOT to be copied
into project snapshots or commits merely because they were inspected.
Nothing was deleted and no upstream environment or build product was
downloaded or created.

`research_archive/reviewed_external_dependencies.json` records each
individual local path, byte SHA-256, exact-commit URL, source commit,
reason, and UTC. Every raw file was additionally checked against its
Git blob SHA in the pinned API tree. Exclusion applies only to matching
reviewed hashes, not arbitrary later edits in that directory. The
independently authored code, this note, timing results, and mathematical
witness data remain ordinary preservable project research.
