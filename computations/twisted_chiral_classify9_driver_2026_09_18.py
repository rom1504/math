#!/usr/bin/env python3
"""Build/run the independent exhaustive order-nine class census."""
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'computations/twisted_chiral_classify9_2026_09_18.cpp'
BINARY=ROOT/'tmp/twisted_chiral_2026_09_18/classify9'


def main():
    subprocess.run(['g++','-O3','-march=native','-std=c++17',str(SOURCE),'-o',str(BINARY)],check=True)
    seeds=json.loads((ROOT/'computations/results/twisted_chiral_2026_09_18_sampled_m9_seeds.json').read_text())
    inp=str(len(seeds))+'\n'+'\n'.join(' '.join(map(str,row)) for a in seeds for row in a['matrix'])+'\n'
    proc=subprocess.Popen([str(BINARY)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
    proc.stdin.write(inp);proc.stdin.close();records=[]
    with (ROOT/'computations/results/twisted_chiral_2026_09_18_classify9.jsonl').open('w') as handle:
        for line in proc.stdout:
            print(line.strip(),flush=True);handle.write(line);handle.flush();records.append(json.loads(line))
    assert proc.wait()==0
    output=dict(schema='twisted-chiral-independent-order9-census-v1',source=str(SOURCE.relative_to(ROOT)),
        classification='complete integer enumeration of every root-gauged signing and all Boolean spins as needed to reject cap>12',
        classes=[x for x in records if x['kind']=='class'],census=records[-1])
    (ROOT/'computations/results/twisted_chiral_2026_09_18_classify9.json').write_text(json.dumps(output,indent=2)+'\n')


if __name__=='__main__':main()
