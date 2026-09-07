"""Exact finite shell-completion test; no asymptotic assertion."""
import argparse
import itertools
import json
import time
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model


def frames(X):
    n = X.shape[1]
    one = len(X) - 1
    candidates = [i for i in range(len(X)) if X[i] @ X[one] == 0]
    normalized = []

    def visit(chosen, remaining):
        if len(chosen) == n - 1:
            normalized.append([one] + chosen)
            return
        if len(chosen) + len(remaining) < n - 1:
            return
        for j, i in enumerate(remaining):
            visit(chosen + [i], [k for k in remaining[j+1:] if X[i] @ X[k] == 0])

    visit([], candidates)
    lookup = {tuple(x): i for i, x in enumerate(X)}
    answer = set()
    for f in normalized:
        for x in X:
            answer.add(tuple(sorted(lookup[tuple(v * x)] for v in X[f])))
    return normalized, sorted(answer)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cap', type=int, default=28)
    ap.add_argument('--seconds', type=float, default=3)
    ap.add_argument('--output', default='computations/results/flatify_independent_2026_09_07_hadamard_free_child.json')
    args = ap.parse_args()
    begin = time.time()
    source = Path('computations/results/flatify_adversary_2026_09_07_hadamard_bridge_phase_search.json')
    representatives = json.loads(source.read_text())['records'][1]['child_representatives']
    X = np.array([[1] + [1 if mask >> i & 1 else -1 for i in range(7)] for mask in range(128)], dtype=np.int64)
    normalized, allframes = frames(X)
    edges = list(itertools.combinations(range(8), 2))
    T = np.array([[int(x[i]*x[j]) for i,j in edges] for x in X], dtype=np.int64)
    records = []
    result = {'status': 'RUNNING', 'normalized_frame_count': len(normalized), 'projective_frame_count': len(allframes), 'cap_target': args.cap, 'records': records}
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    for representative in representatives:
        A = np.array(representative['matrix'], dtype=np.int64)
        hA = np.einsum('bi,ij,bj->b', X,A,X)//2
        assert int(np.max(np.abs(hA))) == 10
        unique_intervals = set()
        for fi, frame in enumerate(allframes):
            B = X[list(frame)].T.copy()
            assert np.array_equal(B @ B.T, 8*np.eye(8,dtype=int))
            bridge = np.abs(X @ B @ X.T)
            lower = np.maximum(np.max(-args.cap+bridge-hA[:,None],axis=0),-10)
            upper = np.minimum(np.min(args.cap-bridge-hA[:,None],axis=0),10)
            key = (tuple(lower),tuple(upper))
            if key in unique_intervals:
                continue
            unique_intervals.add(key)
            record = {'atlas_index': representative['atlas_index'], 'frame_index':fi,'frame':list(frame)}
            if np.any(lower > upper):
                record['status'] = 'EMPTY_INTERVAL'
                records.append(record)
                continue
            model = cp_model.CpModel()
            z = [model.new_bool_var(f'z{i}') for i in range(28)]
            for coefficients, lo, hi in zip(T,lower,upper):
                expression = sum(int(a)*(2*b-1) for a,b in zip(coefficients,z))
                model.add(expression >= int(lo))
                model.add(expression <= int(hi))
            solver = cp_model.CpSolver()
            solver.parameters.max_time_in_seconds = args.seconds
            solver.parameters.num_search_workers = 1
            solver.parameters.random_seed = 0
            code = solver.solve(model)
            record['status'] = solver.status_name(code)
            record['solver_seconds'] = solver.wall_time
            if code in (cp_model.OPTIMAL,cp_model.FEASIBLE):
                D = np.zeros((8,8),dtype=np.int64)
                for variable,(i,j) in zip(z,edges):
                    D[i,j]=D[j,i]=2*solver.value(variable)-1
                hD = np.einsum('bi,ij,bj->b',X,D,X)//2
                value = int(np.max(np.abs(hA[:,None]+hD[None,:])+bridge))
                assert value <= args.cap and int(np.max(np.abs(hD))) == 10
                record.update(A=A.tolist(),B=B.tolist(),D=D.tolist(),parent_cap=value,child_caps=[10,10])
                records.append(record)
                result.update(status='EXACT WITNESS FOUND',elapsed_seconds=time.time()-begin)
                out.write_text(json.dumps(result,indent=2)+'\n')
                print(json.dumps({'status':result['status'],'atlas_index':representative['atlas_index'],'frame':fi,'cap':value}),flush=True)
                return
            records.append(record)
            if len(records)%20 == 0:
                result['elapsed_seconds']=time.time()-begin
                out.write_text(json.dumps(result,indent=2)+'\n')
                print(json.dumps({'tested':len(records),'atlas_index':representative['atlas_index'],'frame':fi,'status':record['status']}),flush=True)
        print(json.dumps({'completed_atlas':representative['atlas_index'],'unique_intervals':len(unique_intervals)}),flush=True)
    result.update(status='EXHAUSTED FINITE CLASS' if all(r['status'] in ('INFEASIBLE','EMPTY_INTERVAL') for r in records) else 'EXHAUSTED WITH UNKNOWN SOLVER CASES',elapsed_seconds=time.time()-begin)
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'records':len(records),'seconds':result['elapsed_seconds']}),flush=True)


if __name__ == '__main__':
    main()
