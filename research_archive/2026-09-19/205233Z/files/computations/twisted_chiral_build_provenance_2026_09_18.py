#!/usr/bin/env python3
"""Generate candidate build-product provenance for the root archive manifest."""
import hashlib
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
SCRATCH=ROOT/'tmp/twisted_chiral_2026_09_18'
PAIRS=[
    ('search','computations/twisted_chiral_search_2026_09_18.cpp'),
    ('search_frozen','computations/twisted_chiral_search_frozen_2026_09_18.cpp'),
    ('search_bound_initial_frozen','computations/twisted_chiral_search_bound_initial_frozen_2026_09_18.cpp'),
    ('verify','computations/exact_fixed_signing_gray.cpp'),
    ('classify9','computations/twisted_chiral_classify9_2026_09_18.cpp'),
    ('classify9_independent','computations/twisted_chiral_classify9_independent_2026_09_18.cpp'),
    ('classify10','computations/twisted_chiral_classify10_2026_09_18.cpp'),
    ('classify10_independent','computations/twisted_chiral_classify10_independent_2026_09_18.cpp'),
    ('search_matching_frozen','computations/twisted_chiral_structured_matching_frozen_2026_09_18.cpp'),
    ('search_extreme_order_frozen','computations/twisted_chiral_search_extreme_order_frozen_2026_09_18.cpp'),
    ('search_width_frozen','computations/twisted_chiral_profile_width_frozen_2026_09_18.cpp'),
    ('joint_search_2026_09_19','computations/twisted_chiral_joint_search_2026_09_19.cpp'),
    ('switching_only_2026_09_19','computations/twisted_chiral_switching_only_2026_09_19.cpp'),
    ('joint_extended_2026_09_19','computations/twisted_chiral_joint_extended_2026_09_19.cpp'),
    ('stability_histogram_2026_09_19','computations/twisted_chiral_stability_histogram_2026_09_19.cpp'),
    ('running_initial_bound_elf','computations/twisted_chiral_search_bound_initial_frozen_2026_09_18.cpp'),
    ('rebuilt_initial_bound_elf','computations/twisted_chiral_search_bound_initial_frozen_2026_09_18.cpp'),
]


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    records=[]
    for name,source in PAIRS:
        binary=SCRATCH/name;src=ROOT/source
        if not binary.exists():continue
        records.append(dict(path=str(binary.relative_to(ROOT)),sha256=digest(binary),
            source=source,source_sha256=digest(src),
            build_command=f'g++ -O3 -march=native -std=c++17 {source} -o {binary.relative_to(ROOT)}',
            classification='rebuildable generated ELF; source and build command retained',mtime_utc_epoch=binary.stat().st_mtime))
        if name in ['joint_search_2026_09_19','switching_only_2026_09_19','joint_extended_2026_09_19']:
            dependency='computations/twisted_chiral_search_frozen_2026_09_18.cpp'
            records[-1]['included_source_dependencies']=[dict(path=dependency,sha256=digest(ROOT/dependency))]
    replay_path=ROOT/'computations/results/twisted_chiral_build_replay_2026_09_19.json'
    if replay_path.exists():
        replay=json.loads(replay_path.read_text())
        for row in replay['records']:
            records.append(dict(path=row['rebuilt_path'],sha256=row['rebuilt_sha256'],
                source=row['source'],source_sha256=row['source_sha256'],
                build_command=' '.join(row['build_command']),
                classification='independent rebuild for executable-section provenance comparison'))
            if row['source'].endswith(('twisted_chiral_joint_search_2026_09_19.cpp','twisted_chiral_switching_only_2026_09_19.cpp','twisted_chiral_joint_extended_2026_09_19.cpp')):
                dependency='computations/twisted_chiral_search_frozen_2026_09_18.cpp'
                records[-1]['included_source_dependencies']=[dict(path=dependency,sha256=digest(ROOT/dependency))]
    result=dict(compiler=subprocess.check_output(['g++','--version'],text=True).splitlines()[0],records=records,
        initial_bound_section_replay=dict(
            note='Initial bound source was reconstructed by removing later width diagnostics. Entire .text and .rodata match the captured running ELF byte for byte; ELF filenames/symbol labels differ.',
            text_sha256=digest(SCRATCH/'running_bound_text.bin'),
            rebuilt_text_sha256=digest(SCRATCH/'rebuilt_bound_text.bin'),
            rodata_sha256=digest(SCRATCH/'running_bound_rodata.bin'),
            rebuilt_rodata_sha256=digest(SCRATCH/'rebuilt_bound_rodata.bin'),
            section_dump_recipe='objcopy --dump-section .text=OUTPUT_TEXT --dump-section .rodata=OUTPUT_RODATA INPUT_ELF OUTPUT_COPY_ELF'),
        executable_section_replay=(str(replay_path.relative_to(ROOT)) if replay_path.exists() else None),
        note='Candidate for root-owned reviewed_build_products manifest; does not mutate that shared manifest.')
    output=ROOT/'computations/results/twisted_chiral_2026_09_18_build_products_candidate.json'
    output.write_text(json.dumps(result,indent=2)+'\n');print(output)


if __name__=='__main__':main()
