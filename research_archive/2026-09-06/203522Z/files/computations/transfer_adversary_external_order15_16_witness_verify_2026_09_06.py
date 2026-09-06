"""Independent integer replay of attributed external witnesses; no lower claim.

Default --replay-preserved reads only our preserved mathematical witness data.
Optional --extract-upstream parses pinned source as literal syntax, never
imports/executes it. The Gray-cube scan is independent of external evaluators.
"""
import ast
import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXTERNAL = ROOT / "external_research/2026-09-06/robby955_mo_413935_61ac2268b6f9"
COMMIT = "61ac2268b6f9234c3ea268e7a7072ac72141f36b"


def literal(path, name):
    for statement in ast.parse(path.read_text()).body:
        if isinstance(statement, ast.Assign):
            targets = statement.targets
        elif isinstance(statement, ast.AnnAssign):
            targets = [statement.target]
        else:
            continue
        if any(isinstance(t, ast.Name) and t.id == name for t in targets):
            return ast.literal_eval(statement.value)
    raise KeyError(name)


def direct_energy(a, x):
    return sum(a[i][j] * x[i] * x[j] for i in range(len(a))
               for j in range(i + 1, len(a)))


def scan(a):
    n = len(a)
    assert all(len(row) == n for row in a)
    assert all(a[i][i] == 0 for i in range(n))
    assert all(a[i][j] == a[j][i] and abs(a[i][j]) == 1
               for i in range(n) for j in range(i + 1, n))
    x = [1] * n
    fields = [sum(row) for row in a]
    energy = direct_energy(a, x)
    histogram = Counter()
    extreme = {}
    for step in range(1 << (n - 1)):
        if step:
            k = (step & -step).bit_length()  # coordinate zero remains fixed
            old_spin = x[k]
            energy -= 2 * old_spin * fields[k]
            for j in range(n):
                fields[j] -= 2 * a[j][k] * old_spin
            x[k] = -old_spin
        if step % 257 == 0:
            assert energy == direct_energy(a, x)
            assert fields == [sum(a[i][j] * x[j] for j in range(n))
                              for i in range(n)]
        histogram[energy] += 1
        extreme.setdefault(energy, x[:])
    assert sum(histogram.values()) == 1 << (n - 1)
    # Full-cube orthogonality; the global-sign quotient has identical moments.
    assert sum(e * c for e, c in histogram.items()) == 0
    assert sum(e * e * c for e, c in histogram.items()) == (
        (1 << (n - 1)) * n * (n - 1) // 2)
    cap = max(abs(e) for e in histogram)
    return {
        "order": n, "cap": cap, "projective_states_checked": 1 << (n - 1),
        "projective_energy_histogram": dict(sorted(histogram.items())),
        "positive_ground_count": histogram[cap],
        "negative_ground_count": histogram[-cap],
        "positive_maximizer": extreme.get(cap),
        "negative_maximizer": extreme.get(-cap),
    }


def extract_upstream():
    file15 = EXTERNAL / "verification/research_order15_certify.py"
    file16 = EXTERNAL / "verification/verify_framed_hadamard_lift_30.py"
    manifest = json.loads((ROOT / "research_archive/reviewed_external_dependencies.json").read_text())
    for source in (file15, file16):
        assert hashlib.sha256(source.read_bytes()).hexdigest() == manifest[str(source.relative_to(ROOT))]["sha256"]
    rows = literal(file15, "WITNESS_ROWS")
    a15 = [[{"0": 0, "+": 1, "-": -1}[c] for c in row] for row in rows]
    hadamard = literal(file16, "ORIENTED_H")
    base = literal(file16, "BASE_SIGNS")
    p = literal(file16, "P_EDGES")
    r = literal(file16, "R_EDGES")
    pairs = list(itertools.combinations(range(4), 2))
    a16 = [[0] * 16 for _ in range(16)]
    for cloud, edges in enumerate((p, r, p, r)):
        for (i, j), sign in zip(pairs, edges):
            a16[4 * cloud + i][4 * cloud + j] = sign
            a16[4 * cloud + j][4 * cloud + i] = sign
    for (left, right), sign in zip(pairs, base):
        for i in range(4):
            for j in range(4):
                a16[4 * left + i][4 * right + j] = sign * hadamard[i][j]
                a16[4 * right + j][4 * left + i] = sign * hadamard[i][j]
    result = {
        "external_repository": "https://github.com/Robby955/mo-413935-research",
        "external_author": "Rob Sneiderman / Robby955",
        "external_commit": COMMIT,
        "external_license": "No LICENSE file in pinned tree; GitHub API license=null. No reuse grant inferred.",
        "scope": "Independently checked witness caps; no independently certified global lower bound.",
        "cases": [],
    }
    for a, source, expected in ((a15, file15, 27), (a16, file16, 30)):
        case = scan(a)
        assert case["cap"] == expected
        case["matrix"] = a
        case["external_source"] = str(source.relative_to(ROOT))
        case["external_source_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
        case["matrix_repr_sha256"] = hashlib.sha256(repr(tuple(map(tuple, a))).encode()).hexdigest()
        deletions = []
        for omitted in range(len(a)):
            vertices = [i for i in range(len(a)) if i != omitted]
            deletions.append(scan([[a[i][j] for j in vertices] for i in vertices])["cap"])
        case["single_vertex_deletion_caps"] = deletions
        case["deletion_cap_histogram"] = dict(sorted(Counter(deletions).items()))
        result["cases"].append(case)
        print(json.dumps({key: case[key] for key in (
            "order", "cap", "positive_ground_count", "negative_ground_count",
            "deletion_cap_histogram", "matrix_repr_sha256")}), flush=True)
    assert result["cases"][1]["matrix_repr_sha256"] == literal(file16, "EXPECTED_MATRIX_SHA256")
    expected_full_hist = dict(literal(file16, "EXPECTED_FULL_HISTOGRAM"))
    assert {e: 2 * c for e, c in result["cases"][1]["projective_energy_histogram"].items()} == expected_full_hist
    output = ROOT / "computations/results/transfer_adversary_external_order15_16_witness_verify_2026_09_06.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS: independent integer witnesses, full histogram match at n=16, and all single deletions")


def replay_preserved():
    """Use only our preserved mathematical data; never read external cache."""
    source = ROOT / "computations/results/transfer_adversary_external_order15_16_witness_verify_2026_09_06.json"
    saved = json.loads(source.read_text())
    assert [case["order"] for case in saved["cases"]] == [15, 16]
    for case, expected_cap in zip(saved["cases"], (27, 30)):
        a = case["matrix"]
        fresh = scan(a)
        assert fresh["cap"] == expected_cap == case["cap"]
        for key in ("order", "projective_states_checked", "positive_ground_count",
                    "negative_ground_count", "positive_maximizer", "negative_maximizer"):
            assert fresh[key] == case[key], (case["order"], key)
        assert fresh["projective_energy_histogram"] == {
            int(e): count for e, count in case["projective_energy_histogram"].items()}
        digest = hashlib.sha256(repr(tuple(map(tuple, a))).encode()).hexdigest()
        assert digest == case["matrix_repr_sha256"]
        deletions = []
        for omitted in range(len(a)):
            vertices = [i for i in range(len(a)) if i != omitted]
            deletions.append(scan([[a[i][j] for j in vertices] for i in vertices])["cap"])
        assert deletions == case["single_vertex_deletion_caps"]
        assert dict(Counter(deletions)) == {int(q): count for q, count in case["deletion_cap_histogram"].items()}
        print(json.dumps({"order": len(a), "cap": fresh["cap"],
                          "full_projective_histogram_replayed": True,
                          "all_single_deletions_replayed": len(deletions),
                          "external_cache_read": False}), flush=True)
    print("PASS: standalone preserved-data replay; no external source dependency")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--extract-upstream", action="store_true",
                       help="Parse pinned external source and regenerate preserved witness data")
    modes.add_argument("--replay-preserved", action="store_true",
                       help="Re-enumerate preserved matrices without external cache (default)")
    args = parser.parse_args()
    if args.extract_upstream:
        extract_upstream()
    else:
        replay_preserved()


if __name__ == "__main__":
    main()
