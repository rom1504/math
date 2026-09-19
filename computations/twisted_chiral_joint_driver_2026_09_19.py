#!/usr/bin/env python3
"""Rebuild and log the bounded joint search; independently check every best."""
import argparse
import hashlib
import json
import subprocess

import twisted_chiral_driver_2026_09_18 as base


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--seconds',type=float,default=3600)
    parser.add_argument('--seed',type=int,default=20260919)
    parser.add_argument('--tag',default='joint_target38')
    args=parser.parse_args()
    base.ensure_binaries()
    source=base.ROOT/'computations/twisted_chiral_joint_search_2026_09_19.cpp'
    dependency=base.ROOT/'computations/twisted_chiral_search_frozen_2026_09_18.cpp'
    binary=base.SCRATCH/'joint_search_2026_09_19'
    if not binary.exists() or max(source.stat().st_mtime,dependency.stat().st_mtime)>binary.stat().st_mtime:
        subprocess.run(['g++','-O3','-march=native','-std=c++17',str(source),'-o',str(binary)],check=True)
    seed_source=base.ROOT/'computations/results/conference_order10_gf9.json'
    a=json.loads(seed_source.read_text())['conference_matrix']
    command=[str(binary),str(args.seconds),str(args.seed),'17','38']
    provenance={'command':command,'source':str(source.relative_to(base.ROOT)),
        'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'dependency':str(dependency.relative_to(base.ROOT)),
        'dependency_sha256':hashlib.sha256(dependency.read_bytes()).hexdigest(),
        'seed_source':str(seed_source.relative_to(base.ROOT)),
        'evidence':'bounded heuristic over child signings with Q(A)<=17 and twists; all d optimized exactly; not a lower-bound certificate'}
    prefix=base.ROOT/f'computations/results/twisted_chiral_2026_09_19_{args.tag}'
    records=[]
    proc=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
    proc.stdin.write(base.matrix_input(a));proc.stdin.close()
    with prefix.with_suffix('.jsonl').open('w') as log:
        for line in proc.stdout:
            log.write(line);log.flush();print(line.strip(),flush=True)
            record=json.loads(line)
            if record['kind'] in ['best','joint_search_finished']:
                record=base.verify(record['searched_child_matrix'],record)
                assert record['child_cap']<=17
                record.update(provenance)
                records.append(record)
                prefix.with_suffix('.json').write_text(json.dumps({'records':records},indent=2)+'\n')
            elif record['kind']=='search_statistics':
                statistics=record
    assert proc.wait()==0
    prefix.with_suffix('.json').write_text(json.dumps({'records':records,'statistics':statistics},indent=2)+'\n')


if __name__=='__main__':main()
