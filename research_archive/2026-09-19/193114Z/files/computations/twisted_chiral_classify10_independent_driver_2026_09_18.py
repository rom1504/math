#!/usr/bin/env python3
"""Bounded two-shard complete root-signing census at child order ten."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'computations/twisted_chiral_classify10_independent_2026_09_18.cpp'
BINARY=ROOT/'tmp/twisted_chiral_2026_09_18/classify10_independent'


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--shards',type=int,default=2)
    parser.add_argument('--seconds',type=float,default=2700);args=parser.parse_args()
    subprocess.run(['g++','-O3','-march=native','-std=c++17',str(SOURCE),'-o',str(BINARY)],check=True)
    seeds=json.loads((ROOT/'computations/results/twisted_chiral_2026_09_18_sampled_m10_seeds.json').read_text())
    inp=str(len(seeds))+'\n'+'\n'.join(' '.join(map(str,row)) for a in seeds for row in a['matrix'])+'\n'
    def run(shard):
        proc=subprocess.Popen([str(BINARY),str(shard),str(args.shards),str(args.seconds)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
        proc.stdin.write(inp);proc.stdin.close();records=[]
        with (ROOT/f'computations/results/twisted_chiral_2026_09_18_classify10_independent_shard{shard}.jsonl').open('w') as handle:
            for line in proc.stdout:
                print(f'shard{shard} {line.strip()}',flush=True);handle.write(line);handle.flush();records.append(json.loads(line))
        assert proc.wait()==0
        return records
    with ThreadPoolExecutor(max_workers=args.shards) as pool:shards=list(pool.map(run,range(args.shards)))
    complete=all(rows[-1]['kind']=='complete_shard' for rows in shards)
    classes={};counts={}
    for rows in shards:
        local=[r for r in rows if r['kind']=='class']
        for cls,count in zip(local,rows[-1]['class_counts']):
            key=cls['canonical_code'];classes.setdefault(key,cls);counts[key]=counts.get(key,0)+count
    if complete:
        assert sum(rows[-1]['root_signings_checked'] for rows in shards)==1<<36
        assert all(counts[key]==classes[key]['root_orbit_size'] for key in classes)
    output=dict(schema='twisted-chiral-independent-order10-census-v1',complete=complete,
        classification=('complete integer enumeration of all root-gauged order10 signings' if complete else 'bounded incomplete enumeration; not a class census or lower certificate'),
        classes=[dict(classes[key],enumerated_count=counts[key]) for key in sorted(classes)],
        shard_summaries=[rows[-1] for rows in shards],source=str(SOURCE.relative_to(ROOT)),method='signed half-edge energy tables; no popcount identity in evaluation')
    if complete:
        prior=json.loads((ROOT/'computations/results/twisted_chiral_2026_09_18_classify10.json').read_text())
        assert prior['complete']
        assert {c['canonical_code']:c['enumerated_count'] for c in prior['classes']}==counts
        output['comparison']='Every class count agrees with separate popcount census'
    (ROOT/'computations/results/twisted_chiral_2026_09_18_classify10_independent.json').write_text(json.dumps(output,indent=2)+'\n')
    print('COMPLETE',complete,'CLASSES',len(classes),flush=True)


if __name__=='__main__':main()
