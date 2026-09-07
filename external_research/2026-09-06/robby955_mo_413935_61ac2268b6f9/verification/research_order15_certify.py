#!/usr/bin/env python3
"""Computer-assisted certificate for F(15)=27.

The value is decided by a threshold tower.  A complete nauty ``geng`` pass
over all 1,018,997,864 graphs on 11 vertices produced the 82,502,142
root-normalized order-12 signings with maximum at most 24; two one-vertex
expansion levels with ``labelg`` canonical deduplication produced the
282,202,131 order-13 signings with maximum at most 24 and then the 1,313,164
order-14 signings with maximum at most 25.  The committed catalogue
``order15_level14_catalogue.g6.gz`` is that final level.  The tower driver,
its run logs, and the sizing notes are committed under ``order15_tower/``.

This script certifies the top of the chain exactly:

* the order-15 witness has maximum exactly 27 (pure-integer spin cube), so
  F(15) <= 27;
* order-15 energies are odd, so F(15) is odd;
* the four order-14 signings with maximum at most 23 have extension minima
  27, 29, 29, 29, reproving that the threshold-23 tower level at order 15
  is empty;
* the committed catalogue has the pinned SHA-256 digest and record count,
  contains those four signings, and (in ``--full`` mode) every one of its
  1,313,164 records has maximum at most 25 and no one-vertex extension of
  any record has maximum at most 25.

Given completeness of the catalogue, the empty extension level means no
order-15 signing has maximum at most 25, so F(15) >= 27 and hence
F(15) = 27.  Completeness of the catalogue is the trust boundary: it rests
on the recorded tower run (complete geng pass at order 12, two expansion
levels, labelg canonical dedup), whose logs are digest-pinned here and
whose driver is committed for full re-derivation.  The lower tower levels
are not replayed by this script.

Bulk energy evaluations use float32 matrix products as an exact container
for small integers: every intermediate value is an integer of magnitude at
most 196, far below the 2^24 exactness threshold of binary32, and every
threshold decision is made on integers recovered exactly from those
products.  A deterministic subsample is re-checked in pure Python integer
arithmetic.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np


DETERMINISTIC_SEED = 413935
ORDER = 15
SIGNING_ORDER = 14
RESIDUAL_ORDER = 13

CATALOGUE_NAME = "order15_level14_catalogue.g6.gz"
CATALOGUE_RECORDS = 1_313_164
CATALOGUE_SHA256 = (
    "a0f0c0d2e21f382bb73923f014b3936cb60f59ac0a7b454324fdcb5c142d4ecc"
)
TOWER25_LOG_SHA256 = (
    "d2fb2f52e78bc7e97e7e33950b0b6c9cc5de5d9cdaae1b596b24bcfe5293f3b7"
)
TOWER23_LOG_SHA256 = (
    "166fe2f6f86e7083cf4947658ca4d41df48a914156a0520a217d239e2ea239bf"
)
LEVEL13_RECORDS = 282_202_131
LEVEL12_RECORDS = 82_502_142
GENG_TOTAL_RECORDS = 1_018_997_864

SIGNING_THRESHOLD = 25
CERTIFIED_VALUE = 27
PRUNE_K = 128

WITNESS_ROWS = (
    "0++++++++++++++",
    "+0+-++----++-++",
    "++0+-++----++-+",
    "+-+0+-++----+++",
    "++-+0+-++----++",
    "+++-+0+-++----+",
    "+-++-+0+-++---+",
    "+--++-+0+-++--+",
    "+---++-+0+-++-+",
    "+----++-+0+-+++",
    "++----++-+0+-++",
    "+++----++-+0+-+",
    "+-++----++-+0++",
    "++-++----++-+0+",
    "++++++++++++++0",
)
WITNESS_MAXIMUM = 27
WITNESS_EDGE_FLIP_MAXIMUM = 29

THRESHOLD23_RECORDS = (
    ("L}akqXXWomdULJ", 21, 27),
    ("L?MRL\\]lBczHyK", 23, 29),
    ("LoSoCurpp]iuYj", 23, 29),
    ("LrluUGNgQcbNHn", 23, 29),
)


def witness_matrix() -> list[list[int]]:
    matrix = []
    for row_index, row in enumerate(WITNESS_ROWS):
        if len(row) != ORDER:
            raise ValueError(("witness row length", row_index))
        entries = []
        for column_index, character in enumerate(row):
            if character == "0":
                if column_index != row_index:
                    raise ValueError(("witness off-diagonal zero", row_index))
                entries.append(0)
            elif character == "+":
                entries.append(1)
            elif character == "-":
                entries.append(-1)
            else:
                raise ValueError(("witness character", character))
        matrix.append(entries)
    for row in range(ORDER):
        if matrix[row][row] != 0:
            raise ValueError(("witness diagonal", row))
        for column in range(ORDER):
            if matrix[row][column] != matrix[column][row]:
                raise ValueError(("witness symmetry", row, column))
    return matrix


def graph6_signing(record: str) -> list[list[int]]:
    """Decode a graph6 record on 13 vertices into an order-14 signing.

    Edge bit one means coefficient -1; the appended root row is all +1.
    """
    if not record:
        raise ValueError("empty graph6 record")
    residual_order = ord(record[0]) - 63
    if residual_order != RESIDUAL_ORDER:
        raise ValueError(("graph6 order", residual_order, RESIDUAL_ORDER))
    edge_count = residual_order * (residual_order - 1) // 2
    data_characters = (edge_count + 5) // 6
    if len(record) != 1 + data_characters:
        raise ValueError(("graph6 length", len(record), 1 + data_characters))

    bits: list[int] = []
    for character in record[1:]:
        value = ord(character) - 63
        if not 0 <= value < 64:
            raise ValueError(("graph6 byte", character))
        bits.extend((value >> shift) & 1 for shift in range(5, -1, -1))
    if any(bits[edge_count:]):
        raise ValueError("nonzero graph6 padding")

    matrix = [[0] * (residual_order + 1) for _ in range(residual_order + 1)]
    bit_index = 0
    for column in range(1, residual_order):
        for row in range(column):
            coefficient = -1 if bits[bit_index] else 1
            matrix[row][column] = matrix[column][row] = coefficient
            bit_index += 1
    for vertex in range(residual_order):
        matrix[vertex][residual_order] = matrix[residual_order][vertex] = 1
    return matrix


def projective_spins(order: int) -> list[tuple[int, ...]]:
    return [
        tuple(
            -1 if mask >> vertex & 1 else 1 for vertex in range(order - 1)
        )
        + (1,)
        for mask in range(1 << (order - 1))
    ]


def energy(matrix: list[list[int]], spins: tuple[int, ...]) -> int:
    return sum(
        matrix[row][column] * spins[row] * spins[column]
        for row in range(len(matrix))
        for column in range(row + 1, len(matrix))
    )


def maximum_absolute_energy(matrix: list[list[int]]) -> int:
    return max(
        abs(energy(matrix, spins)) for spins in projective_spins(len(matrix))
    )


def spin_block(order: int) -> np.ndarray:
    """All sign vectors of the given order with last entry +1, as float32."""
    count = 1 << (order - 1)
    masks = np.arange(count, dtype=np.uint32)
    block = np.ones((count, order), dtype=np.float32)
    for bit in range(order - 1):
        block[:, bit] = np.where((masks >> bit) & 1, -1.0, 1.0)
    return block


def exact_int(array: np.ndarray, bound: int) -> np.ndarray:
    """Recover exact integers from a float32 product and assert the bound."""
    integers = np.rint(array).astype(np.int32)
    if np.abs(integers).max(initial=0) > bound:
        raise AssertionError(("integer bound", int(np.abs(integers).max()), bound))
    return integers


def signing_energies(matrix: list[list[int]], spins: np.ndarray) -> np.ndarray:
    coefficients = np.asarray(matrix, dtype=np.float32)
    doubled = (spins @ coefficients * spins).sum(axis=1)
    doubled_int = exact_int(doubled, 4 * len(matrix) * len(matrix))
    if np.any(doubled_int % 2):
        raise AssertionError("odd doubled energy")
    return doubled_int // 2


def extension_profile(
    matrix: list[list[int]], spins: np.ndarray, block: int = 1024
) -> np.ndarray:
    """Exact M of every one-vertex extension of an order-14 signing."""
    energies = signing_energies(matrix, spins)
    absolute = np.abs(energies)
    profile = np.empty(len(spins), dtype=np.int32)
    for start in range(0, len(spins), block):
        rows = spins[start : start + block]
        inner = exact_int(rows @ spins.T, SIGNING_ORDER)
        profile[start : start + block] = (
            np.abs(inner) + absolute[None, :]
        ).max(axis=1)
    return profile


def verify_witness() -> None:
    matrix = witness_matrix()
    energies = [
        energy(matrix, spins) for spins in projective_spins(ORDER)
    ]
    maximum = max(abs(value) for value in energies)
    if maximum != WITNESS_MAXIMUM:
        raise AssertionError(("witness maximum", maximum))
    if any(value % 2 == 0 for value in energies):
        raise AssertionError("even order-15 energy")
    maximizers = sum(abs(value) == maximum for value in energies)

    corrupted = [row[:] for row in matrix]
    corrupted[0][1] *= -1
    corrupted[1][0] *= -1
    flipped = maximum_absolute_energy(corrupted)
    if flipped == WITNESS_MAXIMUM or (
        WITNESS_EDGE_FLIP_MAXIMUM is not None
        and flipped != WITNESS_EDGE_FLIP_MAXIMUM
    ):
        raise AssertionError(("edge-flip corruption control", flipped))
    print(
        f"order15_witness_M={maximum} maximizers={maximizers} "
        f"edge_flip_control_M={flipped}"
    )
    print("order15_energy_parity=odd")


def verify_threshold23_records(spins: np.ndarray) -> None:
    for record, expected_maximum, expected_extension in THRESHOLD23_RECORDS:
        matrix = graph6_signing(record)
        pure_maximum = maximum_absolute_energy(matrix)
        energies = signing_energies(matrix, spins)
        maximum = int(np.abs(energies).max())
        if maximum != pure_maximum or maximum != expected_maximum:
            raise AssertionError((record, maximum, pure_maximum))
        extension_minimum = int(extension_profile(matrix, spins).min())
        if extension_minimum != expected_extension:
            raise AssertionError((record, extension_minimum))
        print(
            f"order14_record={record} M={maximum} "
            f"extension_minimum={extension_minimum}"
        )
    print("threshold23_extension_floor=27 => level15_threshold23_empty")


def load_catalogue(directory: Path) -> list[bytes]:
    path = directory / CATALOGUE_NAME
    digest = hashlib.sha256()
    records: list[bytes] = []
    with gzip.open(path, "rb") as stream:
        while True:
            chunk = stream.read(1 << 22)
            if not chunk:
                break
            digest.update(chunk)
            records.append(chunk)
    payload = b"".join(records)
    if digest.hexdigest() != CATALOGUE_SHA256:
        raise AssertionError("catalogue digest mismatch")
    lines = payload.decode("ascii").split("\n")
    if lines and not lines[-1]:
        lines.pop()
    if len(lines) != CATALOGUE_RECORDS:
        raise AssertionError(("catalogue records", len(lines)))
    membership = set(lines)
    for record, _, _ in THRESHOLD23_RECORDS:
        if record not in membership:
            raise AssertionError(("catalogue membership", record))
    print(
        f"catalogue_records={len(lines)} sha256={CATALOGUE_SHA256[:16]}... "
        "threshold23_records_present=4"
    )
    return [line.encode("ascii") for line in lines]


def verify_tower_receipts(directory: Path) -> None:
    for name, expected in (
        ("order15_tower/tower25.log", TOWER25_LOG_SHA256),
        ("order15_tower/tower23.log", TOWER23_LOG_SHA256),
    ):
        digest = hashlib.sha256((directory / name).read_bytes()).hexdigest()
        if digest != expected:
            raise AssertionError(("tower receipt digest", name))
    log = (directory / "order15_tower/tower25.log").read_text()
    for line in (
        f"|T12(M<=24)| =  {LEVEL12_RECORDS}",
        f"|T13(M<=24)| =  {LEVEL13_RECORDS}",
        f"|T14(M<=25)| =  {CATALOGUE_RECORDS}",
        "T15(M<=25) EMPTY  =>  F(15) = 27 (witness already verified)",
        "TOWER25_DONE",
    ):
        if line not in log:
            raise AssertionError(("tower receipt line", line))
    print(
        f"tower_receipts=pinned levels {GENG_TOTAL_RECORDS}->"
        f"{LEVEL12_RECORDS}->{LEVEL13_RECORDS}->{CATALOGUE_RECORDS}->0"
    )


_WORKER_SPINS: np.ndarray | None = None


def _worker_init() -> None:
    global _WORKER_SPINS
    _WORKER_SPINS = spin_block(SIGNING_ORDER)


def replay_chunk(chunk: list[bytes]) -> tuple[int, int, int, int]:
    spins = _WORKER_SPINS if _WORKER_SPINS is not None else spin_block(
        SIGNING_ORDER
    )
    transpose = np.ascontiguousarray(spins.T)
    processed = 0
    survivors = 0
    minimum_seen = 10**9
    maximum_seen = 0
    for raw in chunk:
        matrix = graph6_signing(raw.decode("ascii"))
        energies = signing_energies(matrix, spins)
        absolute = np.abs(energies)
        maximum = int(absolute.max())
        minimum_seen = min(minimum_seen, maximum)
        maximum_seen = max(maximum_seen, maximum)
        if maximum > SIGNING_THRESHOLD:
            raise AssertionError(("catalogue record above threshold", raw))
        top = np.argpartition(-absolute, PRUNE_K)[:PRUNE_K]
        inner = exact_int(spins @ spins[top].T, SIGNING_ORDER)
        bounds = (np.abs(inner) + absolute[top][None, :]).max(axis=1)
        candidates = np.flatnonzero(bounds <= SIGNING_THRESHOLD)
        for candidate in candidates:
            row = spins[candidate]
            inner_full = exact_int(spins @ row, SIGNING_ORDER)
            exact_maximum = int((np.abs(inner_full) + absolute).max())
            if exact_maximum <= SIGNING_THRESHOLD:
                survivors += 1
        processed += 1
    return processed, survivors, minimum_seen, maximum_seen


def verify_full(records: list[bytes], jobs: int, sample_stride: int) -> None:
    chunk_size = 4096
    chunks = [
        records[start : start + chunk_size]
        for start in range(0, len(records), chunk_size)
    ]
    processed = 0
    survivors = 0
    minimum_seen = 10**9
    maximum_seen = 0
    with ProcessPoolExecutor(
        max_workers=jobs, initializer=_worker_init
    ) as executor:
        for count, found, minimum, maximum in executor.map(
            replay_chunk, chunks, chunksize=1
        ):
            processed += count
            survivors += found
            minimum_seen = min(minimum_seen, minimum)
            maximum_seen = max(maximum_seen, maximum)
            if processed % (chunk_size * 64) < chunk_size:
                print(
                    f"  replay {processed}/{len(records)}",
                    file=sys.stderr,
                    flush=True,
                )
    if processed != CATALOGUE_RECORDS or survivors != 0:
        raise AssertionError(("full replay", processed, survivors))

    spins = spin_block(SIGNING_ORDER)
    for index in range(0, len(records), sample_stride):
        matrix = graph6_signing(records[index].decode("ascii"))
        pure_maximum = maximum_absolute_energy(matrix)
        numpy_maximum = int(np.abs(signing_energies(matrix, spins)).max())
        if pure_maximum != numpy_maximum or pure_maximum > SIGNING_THRESHOLD:
            raise AssertionError(("pure python cross-check", index))
    print(
        f"catalogue_signing_maxima=[{minimum_seen},{maximum_seen}] "
        f"records={processed} extension_survivors={survivors}"
    )
    print(
        "pure_python_cross_check="
        f"{len(range(0, len(records), sample_stride))}_records"
    )
    print("level15_threshold25=EMPTY")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--sample-stride", type=int, default=65536)
    arguments = parser.parse_args()
    if arguments.jobs < 1 or arguments.jobs > 16:
        parser.error("--jobs must be between 1 and 16")

    directory = Path(__file__).resolve().parent
    verify_witness()
    spins = spin_block(SIGNING_ORDER)
    verify_threshold23_records(spins)
    verify_tower_receipts(directory)
    records = load_catalogue(directory)
    if arguments.full:
        verify_full(records, arguments.jobs, arguments.sample_stride)
        print("certified_value=F(15)=27")
    else:
        print("certified_value=F(15)<=27 and F(15) odd; run --full for the bound F(15)>=27")
    print(
        "certificate_boundary=tower_catalogue_completeness"
        " (geng pass + two labelg-deduplicated expansion levels, receipts pinned)"
    )
    print(f"deterministic_seed={DETERMINISTIC_SEED}")
    print("order15_verification=PASSED")


if __name__ == "__main__":
    main()
