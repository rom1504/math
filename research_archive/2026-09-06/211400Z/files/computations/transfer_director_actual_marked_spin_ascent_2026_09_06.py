"""Exact finite endpoint/ascent check, not asymptotic feedback evidence.

Source f(g,y)=sign(y) 1{|y|>1} is frozen. Since
Y=A[S*((AS)^2-(n-1))]/((n-1)*sqrt(2*(n-1))), its source threshold
is decided by integer squares. No Gaussian surrogate is sampled.
"""
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def check_case(item):
    source = item["source"]
    source_path = ROOT / source["path"]
    data = json.loads(source_path.read_bytes())
    matrix = data
    for key in source["json_pointer"].split("/")[1:]:
        matrix = matrix[int(key)] if isinstance(matrix, list) else matrix[key]
    a = np.asarray(matrix, dtype=np.int64)
    n = len(a)
    m = n-1
    spins = np.asarray(list(itertools.product((-1, 1), repeat=n)), dtype=np.int64)
    energies = np.einsum("bi,ij,bj->b", spins, a, spins)//2
    q = int(np.max(np.abs(energies)))
    assert q == item["cap_verified_exact"]
    linear = spins @ a
    marked = spins * (linear*linear-m)
    y_numerator = marked @ a
    f = np.sign(y_numerator) * (y_numerator*y_numerator > 2*m**3)
    holes = 1-f*f
    af = f @ a
    c = holes*np.where(af < 0, -1, 1)
    j_values = np.einsum("bi,bi->b", f, c @ a)
    assert np.all(j_values >= 0)
    trials = len(spins)
    records = []
    oriented_sums = []
    for endpoint, orientation in ((1, 1), (-1, -1)):
        u = endpoint*f+c
        assert np.all(np.abs(u) == 1)
        au = u @ a
        local = orientation*u*au
        bad = local < 0
        negative_sum = int(np.sum(np.where(bad, -local, 0)))
        mean_negative_per_vertex = F(negative_sum, trials*n)
        # Gershgorin gives ||A||op<=m, with no numerical norm assumption.
        alpha = mean_negative_per_vertex/m
        assert 0 <= alpha <= 1
        masked = u*bad
        quadratic_sum = int(np.einsum("bi,ij,bj->", masked, a, masked))
        baseline_sum = int(orientation*np.einsum("bi,bi->", u, au))
        improvement = alpha*F(negative_sum, trials) + orientation*alpha*alpha*F(quadratic_sum, 2*trials)
        guaranteed = n*mean_negative_per_vertex**2/(2*m)
        assert improvement >= guaranteed >= 0
        baseline = F(baseline_sum, 2*trials)
        assert baseline+improvement <= q
        oriented_sums.append(baseline)
        records.append({
            "endpoint": endpoint, "energy_orientation": orientation,
            "full_cube_seed_count": trials,
            "unstable_seed_coordinate_pairs": int(bad.sum()),
            "mean_negative_field_per_vertex_raw_A": str(mean_negative_per_vertex),
            "alpha": str(alpha), "expected_oriented_energy_before": str(baseline),
            "exact_expected_fractional_ascent": str(improvement),
            "rigorous_Gershgorin_gain_lower": str(guaranteed),
            "expected_oriented_energy_after": str(baseline+improvement),
        })
    first_response = F(int(j_values.sum()), trials)
    assert (oriented_sums[0]+oriented_sums[1])/2 == first_response
    return {
        "case": item["case"], "order": n, "actual_cap_replayed": q,
        "source": source, "source_sha256_replayed": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "frozen_source": "sign(Y) 1{|Y|>1}; tie sign(AF)=+1",
        "source_nonzero_seed_coordinates": int(np.count_nonzero(f)),
        "mean_first_response_raw_energy": str(first_response),
        "normalization": "divide raw energies by n*sqrt(n-1) for this finite B convention",
        "endpoints": records,
    }


def main():
    source = ROOT / "computations/results/transfer_adversary_minimizer_isotropy_2026_09_06.json"
    cases = json.loads(source.read_bytes())["cases"]
    results = [check_case(item) for item in cases]
    payload = {
        "status": "EXACT FINITE INTEGER/RATIONAL REGRESSION; no asymptotic or new optimality claim",
        "actual_signing_cases": len(results),
        "full_seed_configurations": sum(2**r["order"] for r in results),
        "strict_both_orientations": sum(all(F(e["exact_expected_fractional_ascent"]) > 0 for e in r["endpoints"]) for r in results),
        "cases": results,
    }
    destination = ROOT / "computations/results/transfer_director_actual_marked_spin_ascent_2026_09_06.json"
    destination.write_text(json.dumps(payload, indent=2)+"\n")
    print(json.dumps({k: v for k, v in payload.items() if k != "cases"}))
    print("PASS: all seed spins, both energy orientations, exact fractional-ascent inequality")


if __name__ == "__main__":
    main()
