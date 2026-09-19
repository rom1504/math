#!/usr/bin/env python3
"""Exact diagonal-Hadamard-completion audit of saved joint-search witnesses."""
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def main():
    source_rows=[];records=[];seen=set()
    for name in ['joint_target38','joint_target38_child23']:
        source=ROOT/f'computations/results/twisted_chiral_2026_09_19_{name}.json'
        if source.exists():source_rows.extend((source,row) for row in json.loads(source.read_text())['records'])
    for source,row in source_rows:
        if row['matrix_sha256'] in seen:continue
        seen.add(row['matrix_sha256'])
        a=np.asarray(row['child_matrix'],dtype=np.int64)
        d=np.asarray(row['parent_matrix'],dtype=np.int64);n=len(d);square=d@d
        identity=np.eye(n,dtype=np.int64)
        difference=square-(n-1)*identity
        special_polynomial=(square-25*identity)@(square@square-30*square+205*identity)
        candidates=[]
        for h0 in [-1,1]:
            # Every off-diagonal equation is D^2_ij+(h_i+h_j)D_ij=0.
            # Row0 forces all other h_j once the single sign h0 is chosen.
            h=np.asarray([h0]+[int(-square[0,j]*d[0,j]-h0) for j in range(1,n)])
            bad=np.flatnonzero(np.abs(h)!=1)
            candidate=dict(h0=h0,forced_h=h.tolist(),all_entries_signs=not len(bad))
            if len(bad):
                candidate.update(excluded=True,reason='row0 forces a nonsign',first_bad_vertex=int(bad[0]))
            else:
                error=(d+np.diag(h))@(d+np.diag(h))-n*np.eye(n,dtype=np.int64)
                bad_pairs=np.argwhere(error!=0)
                candidate['excluded']=bool(len(bad_pairs))
                if len(bad_pairs):
                    i,j=map(int,bad_pairs[0]);candidate.update(reason='forced diagonal violates another square entry',first_bad_pair=[i,j],square_error=int(error[i,j]))
                else:candidate['reason']='valid Hadamard diagonal'
            candidates.append(candidate)
        records.append(dict(source=str(source.relative_to(ROOT)),matrix_sha256=row['matrix_sha256'],
            child_cap=row['child_cap'],parent_cap=row['cap'],
            child_is_conference=bool(np.array_equal(a@a,(len(a)-1)*np.eye(len(a),dtype=np.int64))),
            square_error_row_supports=np.count_nonzero(difference,axis=1).tolist(),
            satisfies_D2_minus19I_squared_equals36I=bool(np.array_equal(difference@difference,36*identity)),
            satisfies_polynomial_D6_minus55D4_plus955D2_minus5125I=bool(not np.any(special_polynomial)),
            projective_extreme_multiplicities={k:row['independent_verification'][k] for k in ['min_count','max_count']},
            diagonal_completion_candidates=candidates,
            admits_sign_diagonal_Hadamard_completion=any(not c['excluded'] for c in candidates),
            proof='Both possible h0 signs exhausted; row0 uniquely determines all remaining h_j, then every square equation is checked with exact integers'))
    output=ROOT/'computations/results/twisted_chiral_joint_witness_structure_2026_09_19.json'
    output.write_text(json.dumps({'records':records},indent=2)+'\n')
    for r in records:print(r['matrix_sha256'],r['projective_extreme_multiplicities'],r['admits_sign_diagonal_Hadamard_completion'])


if __name__=='__main__':main()
