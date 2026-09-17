"""Solver-free exact geometry of complete ground codes of stored minimizers.

Global optimality labels are imported; caps and every identity below are
independently enumerated.  No frozen LP witness is modified.
"""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import importlib.util
import itertools
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "stored", ROOT / "computations/transfer_adversary_minimizer_isotropy_2026_09_06.py")
stored = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stored)
reports = []
for name, matrix, provenance in stored.cases():
    n = len(matrix)
    if n not in (4, 6, 8, 12, 14):
        continue
    a = np.asarray(matrix, dtype=object)
    words = np.asarray([(1,) + s for s in itertools.product((-1, 1), repeat=n-1)], dtype=object)
    edges = list(itertools.combinations(range(n), 2))
    energy = np.asarray([sum(a[i,j]*x[i]*x[j] for i,j in edges) for x in words], dtype=object)
    cap = int(max(map(abs, energy)))
    selected = np.asarray([abs(v) == cap for v in energy])
    code = words[selected]
    signs = np.asarray([1 if v > 0 else -1 for v in energy[selected]], dtype=object)
    size = len(code)
    cov_num = code.T @ code
    signed_num = code.T @ (signs[:,None]*code)
    identity = np.eye(n, dtype=object)
    gram = code @ code.T
    responses_num = np.abs(words @ code.T).sum(axis=1)
    row_hist = [dict(sorted(Counter(map(int, abs(row))).items())) for row in gram]
    report = {"case":name, "n":n, "cap":cap, "projective_ground_size":size,
              "uniform_covariance_is_identity":bool(np.array_equal(cov_num,size*identity)),
              "uniform_min_all_columns":str(F(int(min(responses_num)),size)),
              "uniform_own_response_range":[str(F(int(min(abs(gram).sum(axis=1))),size)),
                                            str(F(int(max(abs(gram).sum(axis=1))),size))],
              "gram_row_histograms_constant":all(row == row_hist[0] for row in row_hist),
              "first_gram_row_histogram":row_hist[0]}
    if n == 8:
        # Sigma=(1/2)I+bP; exact low-rank projection P.
        b = 2 if name == "n8_class0" else 1
        proj_num = 2*cov_num-size*identity
        proj_den = 2*size*b
        assert np.array_equal(proj_num@proj_num,proj_den*proj_num)
        rank = F(int(np.trace(proj_num)),proj_den)
        projected = [F(int(x@proj_num@x),proj_den) for x in code]
        assert len(set(projected)) == 1
        lower = projected[0]/rank
        assert lower == F(1,2)+b
        report["exact_every_query_law_covariance_certificate"] = {
            "projection_rank":str(rank), "constant_xPx":str(projected[0]),
            "lower_operator_norm":str(lower),
            "projector":[[str(F(int(v),proj_den)) for v in line] for line in proj_num]}
    if n == 12:
        assert np.array_equal(cov_num,size*identity)
        assert sum(signs)==0
        assert np.array_equal(signed_num@signed_num,size*size*identity)
        assert np.array_equal(code@signed_num,size*signs[:,None]*code)
        assert len(set(row_hist[0])) == 3 and row_hist[0] == {0:13,4:6,12:1}
        for sector in (-1,1):
            sector_code = code[signs == sector]
            assert np.array_equal(2*(sector_code.T@sector_code),size*identity+sector*signed_num)
            assert len(sector_code)==10
        assert min(responses_num)==36
        report["exact_signed_covariance_involution"] = True
        report["exact_sector_eigenword_identity"] = True
        report["opposite_sector_overlap_zero"] = bool(np.all((code[signs==1]@code[signs==-1].T)==0))
        report["signed_covariance"]= [[str(F(int(v),size)) for v in line] for line in signed_num]
        assert not all(signed_num[i,j]*a[0,1] == signed_num[0,1]*a[i,j] for i,j in edges)
        report["signed_covariance_radial"] = False
    if n == 14:
        assert np.array_equal(a@a,13*identity)
        assert np.array_equal(cov_num,size*identity)
        assert np.array_equal(13*signed_num,3*size*a)
        assert np.array_equal(signed_num@signed_num, F(9,13)*size*size*identity)
        assert min(responses_num)==392
        minimizers = responses_num == min(responses_num)
        assert np.array_equal(minimizers,selected)
        report["conference_identity_A_squared_13I"] = True
        report["signed_covariance_radial_coefficient"] = "3/13"
        report["all_minimizing_columns_exactly_ground_code"] = True
    reports.append(report)
print(json.dumps({"status":"PASS exact full-cube identities, no optimization", "cases":reports},indent=2))
