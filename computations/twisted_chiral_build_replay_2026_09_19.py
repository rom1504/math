#!/usr/bin/env python3
"""Recompile retained C++ sources and compare executable code/data sections."""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import re
import subprocess

import twisted_chiral_build_provenance_2026_09_18 as provenance


def sections(path):
    data=path.read_bytes()
    listing=subprocess.check_output(['readelf','-SW',str(path)],text=True)
    result={}
    for name in ['.text','.rodata']:
        match=re.search(r'\]\s+'+re.escape(name)+r'\s+\S+\s+\S+\s+(\S+)\s+(\S+)',listing)
        offset,size=(int(value,16) for value in match.groups())
        result[name]=hashlib.sha256(data[offset:offset+size]).hexdigest()
    return result


def rebuild(item):
    source,names=item
    output=provenance.SCRATCH/f'provenance_rebuild_{names[0]}'
    command=['g++','-O3','-march=native','-std=c++17',str(provenance.ROOT/source),'-o',str(output)]
    subprocess.run(command,check=True,capture_output=True,text=True)
    rebuilt_sections=sections(output)
    comparisons=[]
    for name in names:
        original=provenance.SCRATCH/name
        original_sections=sections(original)
        comparisons.append(dict(original=str(original.relative_to(provenance.ROOT)),
            sha256=provenance.digest(original),sections=original_sections,
            executable_code_and_readonly_data_match=original_sections==rebuilt_sections))
    record=dict(source=source,source_sha256=provenance.digest(provenance.ROOT/source),
        rebuilt_path=str(output.relative_to(provenance.ROOT)),rebuilt_sha256=provenance.digest(output),
        build_command=command,rebuilt_sections=rebuilt_sections,comparisons=comparisons)
    print(source,all(x['executable_code_and_readonly_data_match'] for x in comparisons),flush=True)
    return record


def main():
    by_source={}
    for name,source in provenance.PAIRS:
        if (provenance.SCRATCH/name).exists():by_source.setdefault(source,[]).append(name)
    with ThreadPoolExecutor(max_workers=2) as pool:records=list(pool.map(rebuild,by_source.items()))
    complete=all(c['executable_code_and_readonly_data_match'] for r in records for c in r['comparisons'])
    output=provenance.ROOT/'computations/results/twisted_chiral_build_replay_2026_09_19.json'
    output.write_text(json.dumps(dict(all_executable_sections_match=complete,records=records),indent=2)+'\n')
    print(output,complete)


if __name__=='__main__':main()
