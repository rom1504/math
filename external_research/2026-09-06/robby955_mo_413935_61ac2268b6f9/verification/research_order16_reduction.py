#!/usr/bin/env python3
"""Computer-assisted reduction for F(16), and a minimizer census at order 15.

Banked facts used: F(15) = 27 (this repository's order-15 certificate);
F(16) <= 30 (framed-Hadamard construction); even order-16 parity and
monotonicity, so F(16) is 28 or 30.

Any order-16 signing W with M(W) <= 28 has order-15 deletions of odd
maximum at most 28, hence exactly 27.  If additionally some pair-deletion
of W has maximum at most 25, that order-14 signing lies, up to switching
and relabeling, in the committed catalogue of all 1,313,164 order-14
signings with maximum at most 25, and both incident rows of W must extend
the base to order-15 maximum exactly 27.  This script checks, exactly:

* the Paley special case: the order-14 conference matrix has exactly 56
  extension rows of order-15 maximum 27, and all 56 x 56 x 2 completions
  to order 16 have maximum at least 32;
* the general case: over the whole catalogue, exactly 163 bases admit
  such rows, 316 rows in all, and every completion (u, v, mutual sign)
  has maximum at least 30, with 30 attained;
* the census: the 316 order-15 minimizers collapse under nauty labelg to
  exactly 30 rooted canonical classes, each re-verified at maximum 27;
  the committed list ``order15_minimizers_rooted.g6`` is their file.

Consequences, modulo the catalogue completeness chain pinned by the
order-15 certificate: every order-16 signing with maximum at most 28 has
all 120 order-14 pair-deletions at maximum exactly 27, so its order-15
deletions are minimizers all of whose deletions have maximum 27; no such
"locally-27" minimizer is currently known, and the committed census shows
none among the minimizers with a sub-25 deletion.  Deciding F(16) is
thereby reduced to the locally-27 case.

Float32 matrix products are used only as exact containers for small
integers (all values bounded far below 2^24); every threshold decision is
made on exactly recovered integers.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

DETERMINISTIC_SEED = 413935
RESIDUAL_ORDER = 13
CATALOGUE_NAME = "order15_level14_catalogue.g6.gz"
CATALOGUE_SHA256 = (
    "a0f0c0d2e21f382bb73923f014b3936cb60f59ac0a7b454324fdcb5c142d4ecc"
)
CATALOGUE_RECORDS = 1_313_164
CENSUS_NAME = "order15_minimizers_rooted.g6"
CENSUS_SHA256 = (
    "dab8a70fb5840f53d5ef636ef93fd33cebff125c3eb59f78126c0d689fd98b5b"
)

ORDER15_MAX = 27
TARGET16 = 28
PRUNE_K = 160

PALEY_Q = 13
PALEY_GOOD_ROWS = 56
PALEY_EXTENSION_DIST = {27: 56, 29: 3640, 31: 9464, 33: 2912, 35: 312}
PALEY_FAMILY_MIN = 32

BASES_WITH_ROWS = 163
MINIMIZER_ROWS = 316
CATALOGUE_FAMILY_MIN = 30
CENSUS_CLASSES = 30


def spin_block(order: int, projective: bool) -> np.ndarray:
    free = order - 1 if projective else order
    masks = np.arange(1 << free, dtype=np.int64)
    X = np.ones((1 << free, order), dtype=np.float32)
    for bit in range(free):
        X[:, bit] = np.where((masks >> bit) & 1, -1.0, 1.0)
    return X


def exact(values: np.ndarray, bound: int) -> np.ndarray:
    integers = np.rint(values).astype(np.int64)
    if np.abs(integers).max(initial=0) > bound:
        raise AssertionError(("integer bound", int(np.abs(integers).max())))
    return integers


def energies(matrix: np.ndarray, spins: np.ndarray) -> np.ndarray:
    doubled = exact((spins @ matrix.astype(np.float32) * spins).sum(axis=1),
                    4 * matrix.shape[0] ** 2)
    if np.any(doubled % 2):
        raise AssertionError("odd doubled energy")
    return doubled // 2


def graph6_signing(record: str) -> np.ndarray:
    order = ord(record[0]) - 63
    if order != RESIDUAL_ORDER:
        raise ValueError(("graph6 order", order))
    edges = order * (order - 1) // 2
    if len(record) != 1 + (edges + 5) // 6:
        raise ValueError("graph6 length")
    bits: list[int] = []
    for ch in record[1:]:
        value = ord(ch) - 63
        if not 0 <= value < 64:
            raise ValueError("graph6 byte")
        bits.extend((value >> s) & 1 for s in range(5, -1, -1))
    if any(bits[edges:]):
        raise ValueError("nonzero graph6 padding")
    n = order + 1
    matrix = np.ones((n, n), dtype=np.int64)
    k = 0
    for j in range(1, order):
        for i in range(j):
            matrix[i, j] = matrix[j, i] = -1 if bits[k] else 1
            k += 1
    np.fill_diagonal(matrix, 0)
    return matrix


def graph6_encode(bits: list[int], order: int) -> str:
    edges = order * (order - 1) // 2
    padded = bits + [0] * ((-len(bits)) % 6)
    chars = [chr(order + 63)]
    for k in range(0, len(padded), 6):
        value = 0
        for b in padded[k : k + 6]:
            value = (value << 1) | b
        chars.append(chr(value + 63))
    if len(bits) != edges:
        raise AssertionError("edge count")
    return "".join(chars)


def verify_paley() -> None:
    residues = {(a * a) % PALEY_Q for a in range(1, PALEY_Q)}
    n = PALEY_Q + 1
    C = np.ones((n, n), dtype=np.int64)
    for i in range(PALEY_Q):
        for j in range(PALEY_Q):
            if i != j:
                C[i + 1, j + 1] = 1 if (i - j) % PALEY_Q in residues else -1
    np.fill_diagonal(C, 0)
    if not ((C @ C.T) == PALEY_Q * np.eye(n, dtype=np.int64)).all():
        raise AssertionError("conference identity")
    X = spin_block(n, projective=False)
    E = energies(C, X)
    if int(np.abs(E).max()) != 21:
        raise AssertionError("M(C14)")
    A = exact(X @ X.T, n)
    M15 = (np.abs(E)[None, :] + np.abs(A)).max(axis=1)
    values, counts = np.unique(M15, return_counts=True)
    if dict(zip(values.tolist(), counts.tolist())) != PALEY_EXTENSION_DIST:
        raise AssertionError("extension distribution")
    good = np.flatnonzero(M15 == ORDER15_MAX)
    best = 10 ** 9
    for i in good:
        for j in good:
            for e in (1, -1):
                best = min(
                    best,
                    max(
                        int(np.abs(E + A[i] + A[j] + e).max()),
                        int(np.abs(E + A[i] - A[j] - e).max()),
                    ),
                )
    if len(good) != PALEY_GOOD_ROWS or best != PALEY_FAMILY_MIN:
        raise AssertionError(("paley family", len(good), best))
    print(f"paley_stage=M21 good_rows={len(good)} family_min_M16={best}")


def load_catalogue(directory: Path) -> list[str]:
    digest = hashlib.sha256()
    chunks: list[bytes] = []
    with gzip.open(directory / CATALOGUE_NAME, "rb") as stream:
        while True:
            chunk = stream.read(1 << 22)
            if not chunk:
                break
            digest.update(chunk)
            chunks.append(chunk)
    if digest.hexdigest() != CATALOGUE_SHA256:
        raise AssertionError("catalogue digest")
    lines = b"".join(chunks).decode("ascii").split("\n")
    if lines and not lines[-1]:
        lines.pop()
    if len(lines) != CATALOGUE_RECORDS:
        raise AssertionError("catalogue records")
    return lines


_RECORDS: list[str] | None = None
_SPINS: np.ndarray | None = None


def _init(records: list[str]) -> None:
    global _RECORDS, _SPINS
    _RECORDS = records
    _SPINS = spin_block(14, projective=True)


def scan_chunk(bounds: tuple[int, int]) -> tuple[int, list, int]:
    start, stop = bounds
    X = _SPINS
    hits = []
    best = 10 ** 9
    for idx in range(start, stop):
        B = graph6_signing(_RECORDS[idx])
        E = energies(B, X)
        absE = np.abs(E)
        if int(absE.max()) > 25:
            raise AssertionError("catalogue member above 25")
        top = np.argpartition(-absE, PRUNE_K)[:PRUNE_K]
        inner = exact(X @ np.ascontiguousarray(X[top].T), 14)
        cand = np.flatnonzero(
            (np.abs(inner) + absE[top][None, :]).max(axis=1) <= ORDER15_MAX
        )
        rows = []
        for c in cand:
            column = exact(X @ X[int(c)], 14)
            if int((np.abs(column) + absE).max()) <= ORDER15_MAX:
                rows.append(int(c))
        if not rows:
            continue
        hits.append((idx, rows))
        A = np.stack([exact(X @ X[r], 14) for r in rows])
        for i in range(len(rows)):
            for j in range(len(rows)):
                for e in (1, -1):
                    best = min(
                        best,
                        max(
                            int(np.abs(E + A[i] + A[j] + e).max()),
                            int(np.abs(E + A[i] - A[j] - e).max()),
                            int(np.abs(E - A[i] + A[j] - e).max()),
                            int(np.abs(E - A[i] - A[j] + e).max()),
                        ),
                    )
    return stop - start, hits, best


def verify_catalogue(directory: Path, jobs: int) -> list[tuple[int, list]]:
    records = load_catalogue(directory)
    step = 8192
    chunks = [
        (s, min(s + step, len(records))) for s in range(0, len(records), step)
    ]
    processed = 0
    hits: list[tuple[int, list]] = []
    best = 10 ** 9
    with ProcessPoolExecutor(
        max_workers=jobs, initializer=_init, initargs=(records,)
    ) as pool:
        for count, chunk_hits, chunk_best in pool.map(scan_chunk, chunks):
            processed += count
            hits.extend(chunk_hits)
            best = min(best, chunk_best)
    row_total = sum(len(rows) for _, rows in hits)
    if (
        processed != CATALOGUE_RECORDS
        or len(hits) != BASES_WITH_ROWS
        or row_total != MINIMIZER_ROWS
        or best != CATALOGUE_FAMILY_MIN
    ):
        raise AssertionError(
            ("catalogue stage", processed, len(hits), row_total, best)
        )
    print(
        f"catalogue_stage=records{processed} bases_with_rows={len(hits)} "
        f"minimizer_rows={row_total} family_min_M16={best} witnesses_at_28=0"
    )
    return [(idx, rows, records[idx]) for idx, rows in hits]


def verify_census(directory: Path, hits) -> None:
    X14 = spin_block(14, projective=True)
    X15 = spin_block(15, projective=True)
    rooted = []
    for idx, rows, record in hits:
        B = graph6_signing(record)
        for r in rows:
            u = exact(X14[r], 1)
            A15 = np.zeros((15, 15), dtype=np.int64)
            A15[:14, :14] = B
            A15[:14, 14] = u
            A15[14, :14] = u
            if int(np.abs(energies(A15, X15)).max()) != ORDER15_MAX:
                raise AssertionError(("census maximum", idx, r))
            switched = (u[:, None] * u[None, :]) * B
            bits = [
                1 if switched[i, j] < 0 else 0
                for j in range(1, 14)
                for i in range(j)
            ]
            rooted.append(graph6_encode(bits, 14))
    labelg = shutil.which("labelg") or shutil.which("nauty-labelg")
    if labelg is None:
        raise FileNotFoundError("labelg (nauty)")
    canon = subprocess.run(
        [labelg, "-q"],
        input=("\n".join(rooted) + "\n").encode("ascii"),
        capture_output=True,
        check=True,
    ).stdout.decode("ascii").split()
    unique = sorted(set(canon))
    committed = (directory / CENSUS_NAME).read_text().split()
    digest = hashlib.sha256(
        (directory / CENSUS_NAME).read_bytes()
    ).hexdigest()
    if (
        len(rooted) != MINIMIZER_ROWS
        or len(unique) != CENSUS_CLASSES
        or unique != sorted(committed)
        or digest != CENSUS_SHA256
    ):
        raise AssertionError(
            ("census stage", len(rooted), len(unique), digest)
        )
    print(
        f"census_stage=rows{len(rooted)} rooted_classes={len(unique)} "
        "committed_list_matches=yes"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=4)
    arguments = parser.parse_args()
    if arguments.jobs < 1 or arguments.jobs > 16:
        parser.error("--jobs must be between 1 and 16")
    directory = Path(__file__).resolve().parent
    verify_paley()
    hits = verify_catalogue(directory, arguments.jobs)
    verify_census(directory, hits)
    print("reduction=every order-16 signing with M<=28 has all pair-deletions"
          " at exactly 27")
    print("certificate_boundary=catalogue_completeness (order-15 certificate"
          " chain); census counts rooted classes")
    print(f"deterministic_seed={DETERMINISTIC_SEED}")
    print("order16_reduction_verification=PASSED")


if __name__ == "__main__":
    main()
