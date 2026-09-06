# Excluded build-product preservation audit

Audited 2026-09-06, approximately 19:07 UTC, against
`research_archive/2026-09-06/185123Z/manifest.json` and its recorded Git
head `96a5f1ebcf0c400eb093c62c3a2827c01464a121`.

**Result: 77 entries = 27 ELF executables + 50 Python bytecode files.**
Of the executables, 26 have an identified preserved source; one remains
unresolved and should be preserved as a binary pending clarification.
All bytecode entries have identified source. No excluded file was
deleted, copied, executed, rebuilt, or reclassified during this audit.
Only this report was created.

## Evidence and status definitions

`T` means the source is tracked at the manifest's recorded Git head;
this was checked with `git ls-tree`, not inferred from current existence.
`A` means the source is in the manifest as research and its archived
snapshot was checked against both its manifest SHA-256 and live source.
`P` means installed package source exists and is independently verified.

For ordinary ELF rows, the source basename also appears in the binary's
ELF `FILE` symbol. The renamed coset binaries have stronger section-level
matching evidence described below. These findings establish source
preservation, **not bit-for-bit reproducibility of every historical
compiler invocation or source revision**. Builds were not run.

Build shorthand in the table, from `/home/math/quadra`:

- `C17`: `g++ -O3 -std=c++17 SOURCE -o BINARY`.
- `C20`: `g++ -O3 -std=c++20 SOURCE -o BINARY`.
- `C20a`: `g++ -O3 -std=c++20 -include array SOURCE -o BINARY`.
  The preserved `moment_min.cpp` uses `std::array` but omits its include;
  this explicit forced include is a proposed reproducible remedy, not
  a recovered historical flag or a build verified in this audit.

`SOURCE` and `BINARY` mean the exact corresponding paths in that row.
All identified C++ sources use standard-library headers, with no
unpreserved project-local include dependency found. Where no historical
build command was located, the command above is a source-inspection
recipe, not a claim that that exact command produced the old executable.

## ELF executable table: all 27 entries

| Excluded path | Preserved source / build recipe | Status |
|---|---|---|
| `tmp/audit_growing_hadamard_core` | `computations/audit_growing_hadamard_core.cpp`; C17 | T; historical command in `artifacts/cross_order_growing_hadamard_replica.md:473` |
| `tmp/audit_growing_hadamard_core_check` | `computations/audit_growing_hadamard_core.cpp`; C17 | T; matching ELF source name; check-build flags not recovered |
| `tmp/audit_k4_strong_cyclic_nonexistence` | `computations/audit_k4_strong_cyclic_nonexistence.cpp`; C17 | T; build command in source header |
| `tmp/audit_k4_strong_cyclic_nonexistence_warnings` | `computations/audit_k4_strong_cyclic_nonexistence.cpp`; C17, optionally warning flags | T; matching ELF source name; exact warning flags not recovered |
| `tmp/blank_coding/cospectral` | `tmp/blank_coding/cospectral.cpp`; C20 | A; `std::popcount` requires C++20 |
| `tmp/blank_coding/enumerate` | `tmp/blank_coding/enumerate.cpp`; C20 | A; `std::popcount` requires C++20 |
| `tmp/blank_coding/moment_min` | `tmp/blank_coding/moment_min.cpp`; C20a | A; source preserved; missing `<array>` caveat above |
| `tmp/continued_convergence_weave_ascent_2026_09_06` | `computations/continued_convergence_weave_ascent_2026_09_06.cpp`; C17 | T; exact build recipe in `computations/continued_convergence_restricted_weave_2026_09_06.py:44` |
| `tmp/exact_fixed_signing_gray` | `computations/exact_fixed_signing_gray.cpp`; C17 | T |
| `tmp/fresh_limit_algebra_anneal` | `computations/fresh_limit_algebra_anneal.cpp`; C17 | T |
| `tmp/fresh_limit_algebra_core` | `computations/fresh_limit_algebra_core.cpp`; C17 | T |
| `tmp/fresh_limit_tree_pair_partition_audit` | `computations/fresh_limit_tree_pair_partition_audit.cpp`; C17 | T; historical command in `artifacts/fresh_limit_leading_diagram_and_unmarked_audit_2026_09_05.md:62` |
| `tmp/fresh_limit_tree_pair_partition_root_audit` | `computations/fresh_limit_tree_pair_partition_audit.cpp`; C17 | T; same source name and BuildID as preceding binary |
| `tmp/panel_verify/coset_exact` | `computations/coset_terminal_drift_exact.cpp`; C17 | T |
| `tmp/panel_verify/coset_sample` | `computations/coset_terminal_drift_sample.cpp`; C17 | T |
| `tmp/phase2_profile_collision_n8` | `computations/phase2_profile_collision_n8.cpp`; C17 with `-Wall -Wextra -pedantic` | T; command in `computations/results/phase2_profile_collision_n8.json` |
| `tmp/phase2_profile_collision_n8_verify` | `computations/phase2_profile_collision_n8.cpp`; C17 | T; matching ELF source name; exact verify-build flags not recovered |
| `tmp/phase2_profile_collision_search` | `computations/phase2_profile_collision_n8.cpp`; C17 with `-Wall -Wextra -pedantic` | T; command in `artifacts/phase2b_phi6_collision_report.md:43` and result JSON |
| `tmp/phase2_subset_caps_gray` | `computations/phase2_subset_caps_gray.cpp`; C17 | T; historical command in `artifacts/phase2_falsification_report.md:14` |
| `tmp/phase2_subset_caps_gray_verify` | `computations/phase2_subset_caps_gray.cpp`; C17 | T; same source name and BuildID as preceding binary |
| `tmp/phase2b_exact_fixed` | `computations/exact_fixed_signing_gray.cpp`; C17 | T; same source name and BuildID as `tmp/exact_fixed_signing_gray` |
| `tmp/resumed_power16` | `computations/resumed_director_power16_exhaustive.cpp`; C17 | T; source header records C17 build under another output name |
| `tmp/retrieval_panel_2026_08/verifications/coset_drift_enum` | `computations/coset_terminal_drift_exact.cpp`; C17 | T; old filename resolved by exact code/data/relocation section comparison |
| `tmp/retrieval_panel_2026_08/verifications/coset_drift_sample` | `computations/coset_terminal_drift_sample.cpp`; C17 | T; old filename resolved by exact code/data/relocation section comparison |
| `tmp/retrieval_panel_2026_08/verifications/coset_drift_sample_new` | `computations/coset_terminal_drift_sample.cpp`; C17 | T; byte-identical to preceding binary |
| `tmp/search_r33` | No source or documented build command found | **UNRESOLVED: preserve binary pending clarification** |
| `tmp/transfer_seed_restricted_tensor_probe_2026_09_06` | `computations/transfer_seed_restricted_tensor_probe_2026_09_06.cpp`; C17 | A at this snapshot; build command in `artifacts/transfer_seed_restricted_tensor_2026_09_06.md:164` |

### Archived-source identities

The following sources are preserved under the snapshot's `files/`
directory with their original relative paths. Live and archived hashes
both matched the manifest:

| Source | SHA-256 |
|---|---|
| `tmp/blank_coding/cospectral.cpp` | `a0185631f765211b7993237dc4fd5c578655e43f6ab855f0a099c5104941daa9` |
| `tmp/blank_coding/enumerate.cpp` | `b81531707019fd0afec1e3d897ceefc23218a9d56633ee9e30213d17115a40b5` |
| `tmp/blank_coding/moment_min.cpp` | `f7a7b32753b94c5e3ad3469dcfa8dc0fe1caf477fb034ecf5b431994451725a3` |
| `computations/transfer_seed_restricted_tensor_probe_2026_09_06.cpp` | `493a6fe65c21b04964330d8e5bd96275b42163849511957439b3575dd6317378` |

### Renamed coset-source verification

The old source paths named by the three retrieval-panel binaries no
longer exist. Merely finding similar names would not establish provenance.
Read-only `readelf -x` comparisons instead showed identical `.text`,
`.rodata`, `.data`, `.rela.dyn`, and `.rela.plt` sections between:

- `coset_drift_enum` and `tmp/panel_verify/coset_exact`;
- `coset_drift_sample` and `tmp/panel_verify/coset_sample`.

The latter panel binaries identify the tracked renamed sources in their
ELF `FILE` symbols. Those sources were introduced in tracked commit
`3d668e8` and are present at the snapshot head. The sample and sample-new
retrieval binaries are also byte-identical, SHA-256
`4df696d2291902a9331355af2a3fe12c91eb3e625b806d01a9b6030d3ba8a715`.
This resolves their source mapping without rebuilding or running them.

### Unresolved binary: explicit preservation candidate

`tmp/search_r33` is 21,808 bytes, SHA-256
`b02af09feaac26707ff606feaec819d71434d119fb8ee3af6ced4bf884915743`,
ELF BuildID `c80c8ec125e2bf3947aef597f3653851f114aa0f`.
Its non-runtime ELF `FILE` entry is blank. Symbols include `dfs(int,int)`,
`distc(int,int)`, `tar`, `sel`, `cnt`, `rc`, and `nodes`; they do not identify
a preserved source. Searches of tracked/ignored research, archived text,
source names, these distinctive symbols, and Git path history found no
source or complete build command. Similar `r33` mathematical notes do
not establish that they contain this program.

The candidate and hash were reported to the root before this report.
**ELF format is not sufficient grounds to exclude this binary from
preservation.** Preserve it pending clarification or source recovery.
No classification change or binary copy was made by this audit.

## Python bytecode table: all 50 entries

Four files cache tracked research scripts. In the first table, the
excluded-path prefix `E/` expands exactly to
`tmp/eit_pycache/home/math/quadra/extremal_information/experiments/`.

| Excluded path | Source | Status |
|---|---|---|
| `E/build_quadratic_landscape_dataset.cpython-39.pyc` | `extremal_information/experiments/build_quadratic_landscape_dataset.py` | T |
| `E/entropy_overlap_lab.cpython-39.pyc` | `extremal_information/experiments/entropy_overlap_lab.py` | T |
| `E/pinned_query_rate_verify.cpython-39.pyc` | `extremal_information/experiments/pinned_query_rate_verify.py` | T |
| `E/verify_code_replica_hierarchy.cpython-39.pyc` | `extremal_information/experiments/verify_code_replica_hierarchy.py` | T |

The other 46 are standard-library/package caches. In the table below,
`P/` expands exactly to `tmp/eit_pycache/usr/lib/python3.9/`, and `L/`
expands exactly to `/usr/lib/python3.9/`. Every listed source exists and
is owned by Debian package `libpython3.9-minimal:amd64`, installed version
`3.9.2-1+deb11u7`. A read-only `dpkg -V libpython3.9-minimal` returned no
differences, including for `sitecustomize.py`. They are identified vendor
source, not missing research code; they are not claimed to be separately
archived in this repository.

| Excluded path | Source | Status |
|---|---|---|
| `P/_bootlocale.cpython-39.pyc` | `L/_bootlocale.py` | P |
| `P/_collections_abc.cpython-39.pyc` | `L/_collections_abc.py` | P |
| `P/_sitebuiltins.cpython-39.pyc` | `L/_sitebuiltins.py` | P |
| `P/_weakrefset.cpython-39.pyc` | `L/_weakrefset.py` | P |
| `P/abc.cpython-39.pyc` | `L/abc.py` | P |
| `P/codecs.cpython-39.pyc` | `L/codecs.py` | P |
| `P/collections/__init__.cpython-39.pyc` | `L/collections/__init__.py` | P |
| `P/collections/abc.cpython-39.pyc` | `L/collections/abc.py` | P |
| `P/contextlib.cpython-39.pyc` | `L/contextlib.py` | P |
| `P/copyreg.cpython-39.pyc` | `L/copyreg.py` | P |
| `P/encodings/__init__.cpython-39.pyc` | `L/encodings/__init__.py` | P |
| `P/encodings/aliases.cpython-39.pyc` | `L/encodings/aliases.py` | P |
| `P/encodings/latin_1.cpython-39.pyc` | `L/encodings/latin_1.py` | P |
| `P/encodings/utf_8.cpython-39.pyc` | `L/encodings/utf_8.py` | P |
| `P/enum.cpython-39.pyc` | `L/enum.py` | P |
| `P/functools.cpython-39.pyc` | `L/functools.py` | P |
| `P/genericpath.cpython-39.pyc` | `L/genericpath.py` | P |
| `P/heapq.cpython-39.pyc` | `L/heapq.py` | P |
| `P/importlib/__init__.cpython-39.pyc` | `L/importlib/__init__.py` | P |
| `P/importlib/abc.cpython-39.pyc` | `L/importlib/abc.py` | P |
| `P/importlib/machinery.cpython-39.pyc` | `L/importlib/machinery.py` | P |
| `P/importlib/util.cpython-39.pyc` | `L/importlib/util.py` | P |
| `P/io.cpython-39.pyc` | `L/io.py` | P |
| `P/keyword.cpython-39.pyc` | `L/keyword.py` | P |
| `P/linecache.cpython-39.pyc` | `L/linecache.py` | P |
| `P/operator.cpython-39.pyc` | `L/operator.py` | P |
| `P/os.cpython-39.pyc` | `L/os.py` | P |
| `P/pkgutil.cpython-39.pyc` | `L/pkgutil.py` | P |
| `P/posixpath.cpython-39.pyc` | `L/posixpath.py` | P |
| `P/py_compile.cpython-39.pyc` | `L/py_compile.py` | P |
| `P/re.cpython-39.pyc` | `L/re.py` | P |
| `P/reprlib.cpython-39.pyc` | `L/reprlib.py` | P |
| `P/runpy.cpython-39.pyc` | `L/runpy.py` | P |
| `P/site.cpython-39.pyc` | `L/site.py` | P |
| `P/sitecustomize.cpython-39.pyc` | `L/sitecustomize.py` | P |
| `P/sre_compile.cpython-39.pyc` | `L/sre_compile.py` | P |
| `P/sre_constants.cpython-39.pyc` | `L/sre_constants.py` | P |
| `P/sre_parse.cpython-39.pyc` | `L/sre_parse.py` | P |
| `P/stat.cpython-39.pyc` | `L/stat.py` | P |
| `P/token.cpython-39.pyc` | `L/token.py` | P |
| `P/tokenize.cpython-39.pyc` | `L/tokenize.py` | P |
| `P/traceback.cpython-39.pyc` | `L/traceback.py` | P |
| `P/types.cpython-39.pyc` | `L/types.py` | P |
| `P/typing.cpython-39.pyc` | `L/typing.py` | P |
| `P/warnings.cpython-39.pyc` | `L/warnings.py` | P |
| `P/weakref.cpython-39.pyc` | `L/weakref.py` | P |

For any bytecode row, a source-level recreation recipe is
`PYTHONPYCACHEPREFIX=/home/math/quadra/tmp/eit_pycache .venv/bin/python -m py_compile SOURCE`,
using the expanded absolute source path and Python 3.9. The configured
interpreter reports Python 3.9.2. This command is documented, **not run**.
Cache header timestamps and historical source revisions can differ;
the claim is that the source is accounted for, not that old `.pyc` bytes
would be reproduced identically.

## Requested action

Preserve `tmp/search_r33` as unresolved pending clarification. The
remaining exclusions have the source provenance recorded above, subject
to the explicit historical-build and `moment_min.cpp` caveats. This audit
does not authorize deleting any live artifact or changing the archiver's
classification policy.
