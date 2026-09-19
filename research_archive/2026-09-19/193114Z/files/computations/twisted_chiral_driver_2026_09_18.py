#!/usr/bin/env python3
"""Run exact-cap twisted chiral searches and independently verify all bests."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SCRATCH=ROOT/'tmp/twisted_chiral_2026_09_18'
SEARCH=SCRATCH/'search_frozen'
VERIFY=SCRATCH/'verify'
SEARCH_SOURCE=ROOT/'computations/twisted_chiral_search_frozen_2026_09_18.cpp'


def ensure_binaries():
    """Tracked C++ sources are sufficient; scratch binaries are rebuildable."""
    SCRATCH.mkdir(parents=True,exist_ok=True)
    for source,binary in [(SEARCH_SOURCE,SEARCH),
                          (ROOT/'computations/exact_fixed_signing_gray.cpp',VERIFY)]:
        if not binary.exists() or source.stat().st_mtime>binary.stat().st_mtime:
            subprocess.run(['g++','-O3','-march=native','-std=c++17',str(source),'-o',str(binary)],check=True)


def matrix_input(a):
    return str(len(a))+'\n'+'\n'.join(' '.join(map(str,row)) for row in a)+'\n'


def load_seeds(r):
    source=ROOT/f'computations/results/m{r}_minimizer_orbits.json'
    if source.exists():
        data=json.loads(source.read_text())
        return [(f'm{r}_class{c["class"]}',c['representative_matrix'],str(source.relative_to(ROOT))) for c in data['classes']]
    source=ROOT/f'computations/results/exact_m{r}.json'
    data=json.loads(source.read_text())
    return [(f'm{r}_exact',data['matrix'],str(source.relative_to(ROOT)))]


def verify(a,record):
    a=np.asarray(a,dtype=np.int64)
    p=np.asarray(record['p']);s=np.asarray(record['s']);d=np.asarray(record['d'])
    b=s[:,None]*a[np.ix_(p,p)]*s[None,:]
    c=b+np.diag(d)
    parent=np.block([[a,c],[c,-a]])
    assert np.all(np.diag(parent)==0)
    assert np.all(np.abs(parent+np.eye(len(parent),dtype=np.int64))==1)
    result=json.loads(subprocess.check_output([str(VERIFY)],input=matrix_input(parent),text=True))
    child=json.loads(subprocess.check_output([str(VERIFY)],input=matrix_input(a),text=True))
    assert result['cap']==record['cap'],(record,result)
    assert result['min_energy']==-result['max_energy']
    record.update(child_matrix=a.tolist(),bridge_hollow_matrix=b.tolist(),parent_matrix=parent.tolist(),
        child_cap=child['cap'],independent_verification=result,
        matrix_sha256=hashlib.sha256(parent.astype(np.int8).tobytes()).hexdigest(),
        ratio_to_child=record['cap']/(2**1.5*child['cap']),normalized_cap=record['cap']/len(parent)**1.5)
    return record


def main():
    global SEARCH_SOURCE, SEARCH
    parser=argparse.ArgumentParser()
    parser.add_argument('--orders',type=int,nargs='+',default=list(range(3,11)))
    parser.add_argument('--mode',choices=['search','exhaustive','quotient','bound','fastbound','matching','width'],default='search')
    parser.add_argument('--seconds',type=float,default=10)
    parser.add_argument('--seed',type=int,default=20260918)
    parser.add_argument('--tag',default='initial')
    parser.add_argument('--seed-file',type=Path,help='JSON list of {label,matrix}; overrides --orders')
    args=parser.parse_args()
    if args.mode=='matching':
        SEARCH_SOURCE=ROOT/'computations/twisted_chiral_structured_matching_frozen_2026_09_18.cpp'
        SEARCH=SCRATCH/'search_matching_frozen'
    elif args.mode=='fastbound':
        SEARCH_SOURCE=ROOT/'computations/twisted_chiral_search_extreme_order_frozen_2026_09_18.cpp'
        SEARCH=SCRATCH/'search_extreme_order_frozen'
    elif args.mode=='width':
        SEARCH_SOURCE=ROOT/'computations/twisted_chiral_profile_width_frozen_2026_09_18.cpp'
        SEARCH=SCRATCH/'search_width_frozen'
    ensure_binaries()
    records=[]
    seed_groups=([[(item['label'],item['matrix'],str(args.seed_file)) for item in json.loads(args.seed_file.read_text())]]
                 if args.seed_file else [load_seeds(r) for r in args.orders])
    for seed_group in seed_groups:
        for label,a,source in seed_group:
            print(f'START {label} {args.mode}',flush=True)
            log=ROOT/f'computations/results/twisted_chiral_2026_09_18_{args.tag}_{label}.jsonl'
            cmd=[str(SEARCH),'bound' if args.mode=='fastbound' else args.mode,str(args.seconds),str(args.seed)]
            proc=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
            proc.stdin.write(matrix_input(a));proc.stdin.close()
            lines=[]
            with log.open('w') as log_handle:
                for line in proc.stdout:
                    print(label,line.strip(),flush=True);item=json.loads(line);lines.append(item)
                    log_handle.write(json.dumps(item)+'\n');log_handle.flush()
                    if item['kind']=='best':
                        checkpoint=verify(a,dict(item))
                        checkpoint.update(label=label,source=source,command=cmd,random_seed=args.seed,
                            search_source=str(SEARCH_SOURCE.relative_to(ROOT)),search_source_sha256=hashlib.sha256(SEARCH_SOURCE.read_bytes()).hexdigest(),
                            evidence='exact fixed-witness cap; ongoing family enumeration or search',
                            sign_convention='B[i,j]=s[i]*s[j]*A[p[i],p[j]]; p maps destination to source')
                        (ROOT/f'computations/results/twisted_chiral_2026_09_18_{args.tag}_{label}_checkpoint.json').write_text(json.dumps(checkpoint,indent=2)+'\n')
            assert proc.wait()==0
            # These are generated computational outputs, not hand-edited source files.
            best=verify(a,lines[-1])
            best.update(label=label,source=source,command=cmd,random_seed=args.seed,
                search_source=str(SEARCH_SOURCE.relative_to(ROOT)),search_source_sha256=hashlib.sha256(SEARCH_SOURCE.read_bytes()).hexdigest(),
                evidence=('exhaustive optimum over every signed permutation B and matching d for this fixed child'
                    if best['kind']=='complete_family' else ('exhaustive optimum in signed skew-matching G subclass; not a full-family optimum'
                    if best['kind']=='complete_structured_subfamily' else 'exact fixed-witness cap; family search not exhaustive')),
                sign_convention='B[i,j]=s[i]*s[j]*A[p[i],p[j]]; p maps destination indices to source indices')
            records.append(best)
            if best['kind']=='complete_width_family':
                best['evidence']='exhaustive minimum conditional profile half-width over B; parent cap is a witness, not a family optimum'
            output=ROOT/f'computations/results/twisted_chiral_2026_09_18_{args.tag}.json'
            output.write_text(json.dumps({'schema':'twisted-chiral-exact-witness-v1','records':records},indent=2)+'\n')
            print(f'VERIFIED {label}: child={best["child_cap"]}, parent={best["cap"]}, ratio={best["ratio_to_child"]:.9f}',flush=True)


if __name__=='__main__':main()
