#!/usr/bin/env python3
"""One explicitly authorized1200-second continuation with child cap at most23."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import twisted_chiral_driver_2026_09_18 as base


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--seconds',type=float,default=1200)
    parser.add_argument('--seed',type=int,default=2026091903)
    parser.add_argument('--child-limit',type=int,default=23)
    parser.add_argument('--tag',default='joint_target38_child23')
    parser.add_argument('--seed-result',type=Path,default=base.ROOT/'computations/results/twisted_chiral_2026_09_19_joint_target38.json')
    args=parser.parse_args()
    base.ensure_binaries()
    source=base.ROOT/'computations/twisted_chiral_joint_extended_2026_09_19.cpp'
    dependency=base.ROOT/'computations/twisted_chiral_search_frozen_2026_09_18.cpp'
    binary=base.SCRATCH/'joint_extended_2026_09_19'
    if not binary.exists() or max(source.stat().st_mtime,dependency.stat().st_mtime)>binary.stat().st_mtime:
        subprocess.run(['g++','-O3','-march=native','-std=c++17',str(source),'-o',str(binary)],check=True)
    seeds=json.loads(args.seed_result.read_text())['records']
    seed=min(seeds,key=lambda r:(r['cap'],r['score']))
    command=[str(binary),str(args.seconds),str(args.seed),str(args.child_limit),'38']
    provenance=dict(command=command,source=str(source.relative_to(base.ROOT)),
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependency=str(dependency.relative_to(base.ROOT)),dependency_sha256=hashlib.sha256(dependency.read_bytes()).hexdigest(),
        seed_source=str(args.seed_result),seed_parent_sha256=seed['matrix_sha256'],
        evidence='single authorized bounded continuation allowing Q(A)<=23; failure is not an exclusion certificate')
    prefix=base.ROOT/f'computations/results/twisted_chiral_2026_09_19_{args.tag}'
    records=[];statistics=None
    proc=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
    text=base.matrix_input(seed['child_matrix'])+' '.join(map(str,seed['p']))+'\n'+' '.join(map(str,seed['s']))+'\n'
    proc.stdin.write(text);proc.stdin.close()
    with prefix.with_suffix('.jsonl').open('w') as log:
        for line in proc.stdout:
            log.write(line);log.flush();print(line.strip(),flush=True)
            record=json.loads(line)
            if record['kind'] in ['best','joint_search_finished']:
                record=base.verify(record['searched_child_matrix'],record)
                assert record['child_cap']<=args.child_limit
                record.update(provenance);records.append(record)
                prefix.with_suffix('.json').write_text(json.dumps({'records':records},indent=2)+'\n')
            elif record['kind']=='search_statistics':statistics=record
    assert proc.wait()==0 and statistics is not None
    prefix.with_suffix('.json').write_text(json.dumps({'records':records,'statistics':statistics},indent=2)+'\n')


if __name__=='__main__':main()
