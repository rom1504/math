# Finite regression of the multicolor hard compiler

2026-09-07. Status: **Numerical / exact finite assertions PASS**, not an
asymptotic proof. Reproduction:

```bash
.venv/bin/python computations/principle_invent_2026_09_07_multicolor_compiler_check.py
```

The tested magnitude probabilities were `(1/4,1/4,1/2)` with magnitudes
`(0,1,2)`. The script separately verifies every color-degree correction,
connectivity of all repaired color graphs, exact final signed row types,
and the squared defect from the selected parity edges. Random seeds and
the complete printed summary are preserved below.

```json
{
  "status": "finite regression PASS",
  "reports": [
    {
      "m": 257,
      "seed": 20260907,
      "targets": [64, 64, 128],
      "initial_max_degree_error": 24,
      "edge_recolorings": 2931,
      "max_total_incident_changes": 65,
      "max_internal_incident_changes": 34,
      "internal_load_cap": 32.09370429797532,
      "transfer_two_edge_steps": 1317,
      "positive_three_edge_steps": 99,
      "negative_three_edge_steps": 0,
      "parity_edits": [[2, 0, 7]],
      "undirected_squared_defect": 16,
      "all_exact_row_types": true,
      "all_color_graphs_connected": true
    },
    {
      "m": 513,
      "seed": 20260908,
      "targets": [128, 128, 256],
      "initial_max_degree_error": 40,
      "edge_recolorings": 8538,
      "max_total_incident_changes": 107,
      "max_internal_incident_changes": 54,
      "internal_load_cap": 53.89618538422196,
      "transfer_two_edge_steps": 3651,
      "positive_three_edge_steps": 78,
      "negative_three_edge_steps": 334,
      "parity_edits": [],
      "undirected_squared_defect": 0,
      "all_exact_row_types": true,
      "all_color_graphs_connected": true
    }
  ]
}
```

The internal load can exceed its admission threshold by at most two,
because a chosen internal vertex acquires two incidence changes in that
step. This agrees with the proof's threshold-plus-two accounting.
