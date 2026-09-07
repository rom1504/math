"""Replay every anti-weave interval cover, preserving all informative output.

This is a replay harness, not an independent implementation of interval math.
The separate construct-role chord implementation is the independent checker.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--precision", type=int, default=60)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    stem = root / "computations/results/principle_director_2026_09_07_anti_certificate_replay"
    command = [sys.executable, str(root / "computations/principle_synthesis_2026_09_07_anti_weave_interval.py"),
               "--precision", str(args.precision)]
    started = time.monotonic()
    result = subprocess.run(command, cwd=str(root), text=True, capture_output=True)
    stem.with_suffix(".txt").write_text(result.stdout + "\nSTDERR\n" + result.stderr)
    if result.returncode:
        raise RuntimeError("Certificate failed; inspect preserved text output")
    final = [line for line in result.stdout.splitlines() if line.startswith("FINAL_RESULT ")]
    if len(final) != 1:
        raise AssertionError("Missing unique full certificate")
    payload = json.loads(final[0][len("FINAL_RESULT "):])
    points = payload["point_results"]
    assert [p["index"] for p in points] == list(range(1, 21))
    assert all(p["status"] == "DIRECTED_INTERVAL_CERTIFICATE" for p in points)
    assert len(payload["chords"]["intervals"]) == 10
    report = {"status": "ALL_TWENTY_COVERS_REPLAYED", "precision": args.precision,
              "elapsed_seconds": time.monotonic() - started,
              "checked_boxes": sum(p["checked"] for p in points),
              "accepted_leaves": sum(p["accepted"] for p in points),
              "certificate": payload}
    stem.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key != "certificate"}))


if __name__ == "__main__":
    main()
