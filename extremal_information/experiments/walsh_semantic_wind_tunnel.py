#!/usr/bin/env python3
"""Exact small-order semantic wind tunnel for linear-label Walsh tuples.

The rooted relation form classifies tuples up to Walsh-coordinate relabeling.
This experiment asks a different question: how many of those algebraic orbit
states are distinguished by several *scalar extremal* query families?

All reported maxima are obtained by exhaustive Boolean enumeration.  At
``m=1`` we enumerate complete multi-block landscapes.  At ``m=2`` we use a
two-block meet-in-the-middle enumeration: one spin in each block is fixed,
and the two missing orientations are optimized analytically.  Floating point
is used only as an exact BLAS carrier for integers of magnitude below 2**24;
every block product is checked to be integral.

The output is evidence about small orders, not an asymptotic theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import time
from collections import defaultdict
from itertools import combinations, product
from pathlib import Path
from typing import Iterable

import numpy as np


def parity_dot(a: int, b: int) -> int:
    return bin(a & b).count("1") & 1


def sylvester(q: int) -> np.ndarray:
    return np.asarray(
        [[(-1) ** parity_dot(a, u) for u in range(q)] for a in range(q)],
        dtype=np.int16,
    )


def modulation(m: int, label: int) -> np.ndarray:
    q = 1 << m
    return np.asarray(
        [(-1) ** parity_dot(label, v) for u, v in product(range(q), repeat=2)],
        dtype=np.int16,
    )


def pole(m: int, label: int) -> np.ndarray:
    """The Boolean Walsh pole s_label(u,v)=(-1)^(u.v+label.v)."""

    q = 1 << m
    return np.asarray(
        [
            (-1) ** (parity_dot(u, v) ^ parity_dot(label, v))
            for u, v in product(range(q), repeat=2)
        ],
        dtype=np.int16,
    )


def walsh_data(m: int) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
    q = 1 << m
    w = np.kron(sylvester(q), sylvester(q)).astype(np.int16)
    children = []
    for label in range(q):
        d = modulation(m, label)
        children.append(((d[:, None] * w) * d[None, :]).astype(np.int16))
    return w, tuple(children)


def relation_state(labels: tuple[int, ...], m: int) -> tuple[object, ...]:
    k = len(labels)
    gram = tuple(parity_dot(labels[i], labels[j]) for i in range(k) for j in range(k))
    relations: list[int] = []
    rooted: list[int] = []
    omega = (1 << m) - 1
    for mask in range(1 << k):
        value = 0
        for i, label in enumerate(labels):
            if (mask >> i) & 1:
                value ^= label
        if value == 0:
            relations.append(mask)
        if value == omega:
            rooted.append(mask)
    return gram, tuple(relations), tuple(rooted)


def state_key(state: tuple[object, ...]) -> str:
    return json.dumps(state, separators=(",", ":"))


def projective_spins(n: int) -> np.ndarray:
    """All Boolean spins with first coordinate +1."""

    count = 1 << (n - 1)
    values = np.arange(count, dtype=np.uint64)[:, None]
    bits = ((values >> np.arange(n - 1, dtype=np.uint64)) & 1).astype(np.int8)
    return np.concatenate([np.ones((count, 1), dtype=np.int8), 1 - 2 * bits], axis=1)


def graph_edges(k: int) -> tuple[tuple[int, int], ...]:
    return tuple(combinations(range(k), 2))


def exhaustive_components(m: int, labels: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return child, edge, and omega-field component arrays on projective spins."""

    w, children = walsh_data(m)
    n = len(w)
    k = len(labels)
    spins = projective_spins(k * n).astype(np.int16)
    blocks = spins.reshape(len(spins), k, n)
    child = np.zeros(len(spins), dtype=np.int64)
    for i, label in enumerate(labels):
        # Each child's diagonal sums to zero, so the half quadratic is integral.
        child += np.einsum("bi,ij,bj->b", blocks[:, i], children[label], blocks[:, i], optimize=True) // 2
    edge_values = []
    for i, j in graph_edges(k):
        edge_values.append(np.einsum("bi,ij,bj->b", blocks[:, i], w, blocks[:, j], optimize=True))
    edge = np.stack(edge_values, axis=1) if edge_values else np.zeros((len(spins), 0), dtype=np.int64)
    omega = (1 << m) - 1
    root = pole(m, omega)
    # q*s_omega is exactly the field W*y_omega from the rooted Walsh query.
    fields = (1 << m) * np.einsum("bkn,n->bk", blocks, root, optimize=True)
    return child, edge, fields.astype(np.int64)


def max_rows(values: np.ndarray) -> tuple[int, ...]:
    return tuple(int(x) for x in np.max(values, axis=0))


def m1_signatures(k: int) -> dict[tuple[int, ...], dict[str, tuple[int, ...]]]:
    """Complete graph/weight/root signatures for every m=1 label tuple."""

    m = 1
    edges = graph_edges(k)
    unweighted = np.asarray(list(product((0, 1), repeat=len(edges))), dtype=np.int16)
    signed = np.asarray(list(product((-1, 0, 1), repeat=len(edges))), dtype=np.int16)
    root_masks = np.asarray(list(product((0, 1), repeat=k)), dtype=np.int16)
    answer: dict[tuple[int, ...], dict[str, tuple[int, ...]]] = {}
    for labels in product(range(1 << m), repeat=k):
        child, edge, fields = exhaustive_components(m, labels)
        graph_values = child[:, None] + edge @ unweighted.T
        weighted_values = child[:, None] + edge @ signed.T
        # The root family uses every unweighted graph and every subset of
        # canonical omega fields.  Absolute value accounts for the omitted
        # global Boolean orientation after projectivization.
        root_values = child[:, None, None]
        root_values = root_values + (edge @ unweighted.T)[:, :, None]
        root_values = root_values + np.abs((fields @ root_masks.T)[:, None, :])
        answer[tuple(labels)] = {
            "graphs": max_rows(graph_values),
            "signed_weights": max_rows(weighted_values),
            "omega_roots": tuple(int(x) for x in np.max(root_values, axis=0).ravel()),
        }
    return answer


def extension_signatures_m1(
    base_k: int,
    full: dict[tuple[int, ...], dict[str, tuple[int, ...]]],
) -> dict[tuple[int, ...], tuple[int, ...]]:
    """One-block futures: append 0, omega, or any repeated base label."""

    answer: dict[tuple[int, ...], tuple[int, ...]] = {}
    for labels in product((0, 1), repeat=base_k):
        extension_labels = sorted({0, 1, *labels})
        signature: list[int] = []
        for new_label in extension_labels:
            signature.extend(full[tuple(labels) + (new_label,)]["graphs"])
        answer[tuple(labels)] = tuple(signature)
    return answer


def pair_representatives_m2() -> list[tuple[int, int]]:
    states: dict[str, tuple[int, int]] = {}
    for labels in product(range(4), repeat=2):
        key = state_key(relation_state(tuple(labels), 2))
        states[key] = min(states.get(key, tuple(labels)), tuple(labels))
    return sorted(states.values())


def m2_pair_signatures(
    chunk: int = 128,
    evaluate_roots: bool = True,
) -> dict[tuple[int, int], dict[str, tuple[int, ...]]]:
    """Exact two-block signatures by meet-in-the-middle enumeration."""

    m = 2
    q = 1 << m
    w, children = walsh_data(m)
    n = len(w)
    spins = projective_spins(n)
    spins64 = spins.astype(np.int64)
    child_energy = np.stack(
        [np.einsum("bi,ij,bj->b", spins64, c.astype(np.int64), spins64, optimize=True) // 2 for c in children]
    ).astype(np.int16)
    root = pole(m, (1 << m) - 1).astype(np.int64)
    root_field = (q * (spins64 @ root)).astype(np.int16)
    left_walsh = spins.astype(np.float32) @ w.astype(np.float32)

    reps = pair_representatives_m2()
    graph_maxima = np.full((len(reps), 3), -32000, dtype=np.int16)

    chunks = list(range(0, len(spins), chunk))
    for ci, start in enumerate(chunks):
        stop = min(start + chunk, len(spins))
        bridge_float = left_walsh[start:stop] @ spins.T.astype(np.float32)
        bridge = bridge_float.astype(np.int16)
        if not np.array_equal(bridge_float, bridge.astype(np.float32)):
            raise AssertionError("nonintegral BLAS bridge product")
        base = np.stack(
            [child_energy[a, start:stop, None] + child_energy[b, None, :] for a, b in reps]
        ).astype(np.int16)
        absolute_bridge = np.abs(bridge)[None, :, :]
        for weight in (0, 1, 2):
            values = base + weight * absolute_bridge
            graph_maxima[:, weight] = np.maximum(
                graph_maxima[:, weight], np.max(values, axis=(1, 2)).astype(np.int16)
            )
        if ci % max(1, len(chunks) // 8) == 0:
            print(f"m=2 graph sweep {ci + 1}/{len(chunks)}", flush=True)

    # Only states still colliding after weighted scalar graph queries need
    # rooted queries to decide the nested semantic quotient.  This adaptive
    # subfamily is frozen by the collision partition, not by observed rooted
    # values, and avoids evaluating irrelevant 2**30 pair tables.
    by_weighted: defaultdict[tuple[int, ...], list[int]] = defaultdict(list)
    for ri in range(len(reps)):
        by_weighted[tuple(int(x) for x in graph_maxima[ri])].append(ri)
    unresolved = sorted({ri for values in by_weighted.values() if len(values) > 1 for ri in values})

    root_queries = tuple((weight, r1, r2) for weight in (0, 1, 2) for r1, r2 in ((1, 0), (0, 1), (1, 1)))
    root_maxima = np.full((len(reps), len(root_queries)), -32000, dtype=np.int16)
    if unresolved and evaluate_roots:
        for ci, start in enumerate(chunks):
            stop = min(start + chunk, len(spins))
            bridge_float = left_walsh[start:stop] @ spins.T.astype(np.float32)
            bridge = bridge_float.astype(np.int16)
            if not np.array_equal(bridge_float, bridge.astype(np.float32)):
                raise AssertionError("nonintegral BLAS bridge product")
            hx = root_field[start:stop, None]
            hy = root_field[None, :]
            for ri in unresolved:
                a, b = reps[ri]
                base = child_energy[a, start:stop, None] + child_energy[b, None, :]
                for qi, (weight, r1, r2) in enumerate(root_queries):
                    # Representatives have x_0=y_0=+1. sigma restores the
                    # relative orientation; tau restores common orientation.
                    plus = base + weight * bridge + np.abs(r1 * hx + r2 * hy)
                    minus = base - weight * bridge + np.abs(r1 * hx - r2 * hy)
                    value = int(max(np.max(plus), np.max(minus)))
                    if value > root_maxima[ri, qi]:
                        root_maxima[ri, qi] = value
            if ci % max(1, len(chunks) // 8) == 0:
                print(f"m=2 root sweep {ci + 1}/{len(chunks)} unresolved={len(unresolved)}", flush=True)

    answer: dict[tuple[int, int], dict[str, tuple[int, ...]]] = {}
    for labels in reps:
        ri = reps.index(labels)
        none = tuple(int(v) for v in graph_maxima[ri])
        roots = (
            tuple(int(v) for v in root_maxima[ri])
            if evaluate_roots and ri in unresolved
            else ()
        )
        answer[labels] = {
            "graphs": none[:2],  # weights 0 and 1
            "weighted_bridges": none,
            "adaptive_omega_roots": roots,
        }
    return answer


def extension_signatures_m2(
    pairs: dict[tuple[int, int], dict[str, tuple[int, ...]]]
) -> dict[tuple[int], tuple[int, ...]]:
    """Append a repeated label, 0, or omega to a singleton base."""

    # Map every ordered pair to its rooted-state representative.
    rep_by_state = {state_key(relation_state(labels, 2)): labels for labels in pairs}
    answer: dict[tuple[int], tuple[int, ...]] = {}
    singleton_reps: dict[str, int] = {}
    for label in range(4):
        key = state_key(relation_state((label,), 2))
        singleton_reps[key] = min(singleton_reps.get(key, label), label)
    for label in sorted(singleton_reps.values()):
        signature: list[int] = []
        for appended in sorted({0, 3, label}):
            rep = rep_by_state[state_key(relation_state((label, appended), 2))]
            signature.extend(pairs[rep]["graphs"])
            signature.extend(pairs[rep]["weighted_bridges"])
            signature.extend(pairs[rep]["adaptive_omega_roots"])
        answer[(label,)] = tuple(signature)
    return answer


def summarize(
    signatures: dict[tuple[int, ...], tuple[int, ...]],
    m: int,
) -> dict[str, object]:
    by_signature: defaultdict[tuple[int, ...], list[tuple[int, ...]]] = defaultdict(list)
    state_of: dict[tuple[int, ...], str] = {}
    for labels, signature in signatures.items():
        by_signature[signature].append(labels)
        state_of[labels] = state_key(relation_state(labels, m))
    collisions = [values for values in by_signature.values() if len(values) > 1]
    return {
        "orbit_states": len({state_of[x] for x in signatures}),
        "semantic_classes": len(by_signature),
        "exposed_log2_classes": math.log2(len(by_signature)),
        "largest_collision": max((len(x) for x in by_signature.values()), default=0),
        "collision_classes": [[list(labels) for labels in values] for values in collisions[:20]],
        "all_collisions_listed": len(collisions) <= 20,
    }


def nested_report(raw: dict[tuple[int, ...], dict[str, tuple[int, ...]]], m: int) -> dict[str, object]:
    keys = list(next(iter(raw.values())).keys())
    result: dict[str, object] = {}
    cumulative: dict[tuple[int, ...], tuple[int, ...]] = {labels: () for labels in raw}
    for key in keys:
        for labels in raw:
            cumulative[labels] += raw[labels][key]
        result[key] = summarize(cumulative, m)
    return result


def digest_payload(payload: object) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--m2-chunk", type=int, default=128)
    parser.add_argument("--skip-m2", action="store_true")
    parser.add_argument("--skip-m2-roots", action="store_true")
    parser.add_argument("--include-raw", action="store_true")
    args = parser.parse_args()
    started = time.time()

    m1_k4 = m1_signatures(4)
    m1_k3 = m1_signatures(3)
    ext1 = extension_signatures_m1(3, m1_k4)

    payload: dict[str, object] = {
        "schema": "linear-label-walsh-semantic-wind-tunnel-v1",
        "classification": "exact finite exhaustive computation; not asymptotic evidence",
        "query_families_frozen_before_search": {
            "graphs": "all labelled unweighted simple graphs on existing blocks",
            "signed_weights": "all existing edge weights in {-1,0,1}",
            "omega_roots": "all unweighted graphs and all subsets of fields q*s_omega",
            "extensions": "append one block labelled 0, omega, or a repeated base label",
        },
        "software": {"python": platform.python_version(), "numpy": np.__version__},
        "m1_k4_nested": nested_report(m1_k4, 1),
        "m1_k3_nested": nested_report(m1_k3, 1),
        "m1_k3_extensions": summarize(ext1, 1),
    }

    m1_k4_raw = {str(labels): {k: list(v) for k, v in sig.items()} for labels, sig in m1_k4.items()}
    m1_k3_extension_raw = {str(labels): list(sig) for labels, sig in ext1.items()}
    payload["raw_signature_sha256"] = {
        "m1_k4": digest_payload(m1_k4_raw),
        "m1_k3_extensions": digest_payload(m1_k3_extension_raw),
    }
    if args.include_raw:
        payload["m1_k4_raw"] = m1_k4_raw
        payload["m1_k3_extension_raw"] = m1_k3_extension_raw

    if not args.skip_m2:
        pairs = m2_pair_signatures(args.m2_chunk, evaluate_roots=not args.skip_m2_roots)
        payload["m2_roots_evaluated"] = not args.skip_m2_roots
        ext2 = extension_signatures_m2(pairs)
        payload["m2_k2_nested"] = nested_report(pairs, 2)
        payload["m2_k1_extensions"] = summarize(ext2, 2)
        m2_k2_raw = {str(labels): {k: list(v) for k, v in sig.items()} for labels, sig in pairs.items()}
        m2_k1_extension_raw = {str(labels): list(sig) for labels, sig in ext2.items()}
        payload["raw_signature_sha256"].update(
            {
                "m2_k2": digest_payload(m2_k2_raw),
                "m2_k1_extensions": digest_payload(m2_k1_extension_raw),
            }
        )
        # The m=2 table is small and is the useful collision witness.
        payload["m2_k2_raw"] = m2_k2_raw
        payload["m2_k1_extension_raw"] = m2_k1_extension_raw

    payload["elapsed_seconds"] = time.time() - started
    payload["content_sha256_without_digest"] = digest_payload(payload)
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
        print(f"wrote {args.output}")
    else:
        print(text)
    print(f"elapsed={payload['elapsed_seconds']:.3f}s digest={payload['content_sha256_without_digest']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
