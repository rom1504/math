#!/usr/bin/env python3
"""Run/resume the fixed-child bridge test grid through total order 12.

For every exact-value pair m<=n with m+n<=12 and both relative child signs,
the driver asks whether the saved exact representatives can attain the known
optimal parent cap.  A feasible result is therefore an optimal fixed-child
composition.  Infeasibility only obstructs those representatives and that
relative orientation, not all minimizers.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


EXACT_M = {3: 3, 4: 4, 5: 4, 6: 5, 7: 9, 8: 10, 9: 12, 10: 13, 11: 17, 12: 18}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--time-limit", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("computations/results/bridge_grid_through_12.json"),
    )
    args = parser.parse_args()
    result_dir = args.output.parent
    log_dir = Path("computations/logs")
    result_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for m in range(3, 10):
        for n in range(m, 10):
            total = m + n
            if total > 12:
                continue
            for sign_b in (1, -1):
                stem = f"bridge_grid_{m}_{n}_sign{sign_b}_cap{EXACT_M[total]}"
                output_path = result_dir / f"{stem}.json"
                log_path = log_dir / f"{stem}.log"
                if args.force or not output_path.exists():
                    command = [
                        sys.executable,
                        "computations/bridge_block_cpsat.py",
                        f"computations/results/exact_m{m}.json",
                        f"computations/results/exact_m{n}.json",
                        "--sign-b",
                        str(sign_b),
                        "--decision-cap",
                        str(EXACT_M[total]),
                        "--time-limit",
                        str(args.time_limit),
                        "--workers",
                        str(args.workers),
                        "--output",
                        str(output_path),
                    ]
                    with log_path.open("w") as log:
                        completed = subprocess.run(
                            command,
                            stdout=log,
                            stderr=subprocess.STDOUT,
                            check=False,
                        )
                    if completed.returncode not in (0, 2):
                        raise RuntimeError((command, completed.returncode))
                if not output_path.exists():
                    rows.append(
                        {
                            "m": m,
                            "n": n,
                            "sign_b": sign_b,
                            "target_cap": EXACT_M[total],
                            "status": "NO_RESULT",
                            "result": str(output_path),
                            "log": str(log_path),
                        }
                    )
                    continue
                payload = json.loads(output_path.read_text())
                status = payload["solver"]["status"]
                row = {
                    "m": m,
                    "n": n,
                    "sign_b": sign_b,
                    "target_cap": EXACT_M[total],
                    "status": status,
                    "wall_time_seconds": payload["solver"]["wall_time_seconds"],
                    "result": str(output_path),
                    "log": str(log_path),
                }
                if "parent_profile" in payload:
                    row["verified_parent_M"] = payload["parent_profile"]["M"]
                    row["matrix_sha256"] = payload["parent_matrix_sha256"]
                rows.append(row)
                print(
                    f"{m}+{n} sign={sign_b:+d} target={EXACT_M[total]} "
                    f"status={status} wall={row['wall_time_seconds']:.3f}s",
                    flush=True,
                )
    summary = {
        "schema": "quadratic-signing-fixed-child-bridge-grid-v1",
        "classification": (
            "solver-certified finite tests of saved exact representatives; "
            "not an enumeration of all minimizers"
        ),
        "known_exact_M": {str(n): value for n, value in sorted(EXACT_M.items())},
        "rows": rows,
    }
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
