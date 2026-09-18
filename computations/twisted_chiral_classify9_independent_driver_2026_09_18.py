#!/usr/bin/env python3
"""Independent all-root-signing census using half-edge signed-sum tables."""
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'computations/twisted_chiral_classify9_independent_2026_09_18.cpp'
BINARY=ROOT/'tmp/twisted_chiral_2026_09_18/classify9_independent'


def main():
    subprocess.run(['g++','-O3','-march=native','-std=c++17',str(SOURCE),'-o',str(BINARY)],check=True)
    prior=json.loads((ROOT/'computations/results/twisted_chiral_2026_09_18_classify9.json').read_text())
    seeds=prior['classes']
    inp=str(len(seeds))+'\n'+'\n'.join(' '.join(map(str,row)) for a in seeds for row in a['matrix'])+'\n'
    proc=subprocess.Popen([str(BINARY)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
    proc.stdin.write(inp);proc.stdin.close();records=[]
    with (ROOT/'computations/results/twisted_chiral_2026_09_18_classify9_independent.jsonl').open('w') as handle:
        for line in proc.stdout:
            print(line.strip(),flush=True);handle.write(line);handle.flush();records.append(json.loads(line))
    assert proc.wait()==0
    for key in ['root_signings_checked','global_minimum','minimizer_count','class_counts']:
        assert records[-1][key]==prior['census'][key],key
    output=dict(schema='twisted-chiral-independent-order9-half-edge-census-v1',source=str(SOURCE.relative_to(ROOT)),
        method='Every root-gauged signing; direct addition of independently precomputed signed half-edge energy tables; not popcount energy identity',
        comparison='All counts and minimum exactly match popcount census',census=records[-1])
    (ROOT/'computations/results/twisted_chiral_2026_09_18_classify9_independent.json').write_text(json.dumps(output,indent=2)+'\n')


if __name__=='__main__':main()
