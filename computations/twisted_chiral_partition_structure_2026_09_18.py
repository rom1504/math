#!/usr/bin/env python3
"""Audit the exact induced-interval structure of a Boolean bilinear maximizer."""

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def endpoints(a):
    n = len(a)
    ids = np.arange(1 << n, dtype=np.int64)
    x = 1-2*((ids[:, None] >> np.arange(n)) & 1)
    e = np.einsum("bi,ij,bj->b", x, a, x)//2
    return int(e.max()), int(-e.min())


def main():
    source = ROOT / "computations/results/twisted_chiral_bilinear_audit_2026_09_18.json"
    data = json.loads(source.read_text())
    records = []
    for item in data["profiles"]:
        a = np.asarray(item["matrix"], dtype=np.int64)
        x = np.asarray(item["beta_witness_x"], dtype=np.int64)
        y = np.asarray(item["beta_witness_y"], dtype=np.int64)
        switched = a*x[:, None]*x[None, :]
        t = x*y
        inside, outside = t > 0, t < 0
        aa, bb = switched[np.ix_(inside, inside)], switched[np.ix_(outside, outside)]
        ui, vi = endpoints(aa)
        uj, vj = endpoints(bb)
        internal = (switched*(t[:, None] == t[None, :])).sum(axis=1)*t
        cross = (switched*(t[:, None] != t[None, :])).sum(axis=1)
        assert np.all(internal >= abs(cross))
        assert 2*(ui+vj) == item["beta"]
        assert int(aa.sum()//2) == ui
        assert int(-bb.sum()//2) == vj
        midpoint_difference_twice = ui-vi-uj+vj
        width_sum_twice = ui+vi+uj+vj
        assert width_sum_twice + midpoint_difference_twice == item["beta"]
        assert max(ui+uj, vi+vj) <= item["Q"]
        record = {"n": len(a), "Q": item["Q"], "beta": item["beta"],
                  "source_sha256_int8": item["sha256_int8"],
                  "partition_sizes": [int(inside.sum()), int(outside.sum())],
                  "I": np.flatnonzero(inside).tolist(), "J": np.flatnonzero(outside).tolist(),
                  "U_I": ui, "V_I": vi, "U_J": uj, "V_J": vj,
                  "midpoint_difference": midpoint_difference_twice/2,
                  "sum_induced_half_ranges": width_sum_twice/2,
                  "internal_signed_fields": internal.tolist(), "cross_row_sums": cross.tolist()}
        records.append(record)
        print(json.dumps({k: record[k] for k in ["n", "Q", "beta", "partition_sizes", "U_I", "V_I", "U_J", "V_J"]}))
    destination = ROOT / "computations/results/twisted_chiral_partition_structure_2026_09_18.json"
    destination.write_text(json.dumps({"source": str(source.relative_to(ROOT)), "records": records}, indent=2)+"\n")


if __name__ == "__main__":
    main()
