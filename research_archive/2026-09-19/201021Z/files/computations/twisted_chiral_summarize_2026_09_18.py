#!/usr/bin/env python3
"""Consolidate only completed family certificates, retaining source pointers."""
from datetime import datetime,timezone
import json
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/'computations/results'
PREFIX='twisted_chiral_2026_09_18_'
FAMILY_TAGS=['exhaustive_histograms','bound_r9','bound_r10','census9_group1',
             'census9_group2','census9_quick','census10_new','census9_fastgroup1','census10_fastnew']
SKEW_TAGS=['skew_matching','skew_matching_census9','skew_matching_census10_new']
WIDTH_TAGS=['width_min_small','width_bad_census9']
CHILD_CAP={3:3,4:4,5:4,6:5,7:9,8:10,9:12,10:13}
CLASS_COUNTS={3:1,4:1,5:1,6:1,7:3,8:2,9:9,10:2}
KNOWN_GLOBAL_PARENT={6:5,8:10,10:13,12:18,14:21}


def key(label):
    r=int(re.match(r'm(\d+)',label).group(1))
    if 'exact' in label:cls=0
    else:cls=int(re.search(r'(?:class|census|sample)(\d+)',label).group(1))
    return r,cls


def collect(tags,kind):
    records={}
    for tag in tags:
        path=RESULTS/f'{PREFIX}{tag}.json'
        if not path.exists():continue
        data=json.loads(path.read_text())
        for row in data['records']:
            if row['kind']!=kind:continue
            k=key(row['label']);assert row['child_cap']==CHILD_CAP[k[0]]
            assert row['independent_verification']['cap']==row['cap']
            if k in records:assert records[k]['cap']==row['cap'],(k,records[k],row)
            records[k]=dict(row,result_source=str(path.relative_to(ROOT)))
    return records


def main():
    families=collect(FAMILY_TAGS,'complete_family')
    skew=collect(SKEW_TAGS,'complete_structured_subfamily')
    width=collect(WIDTH_TAGS,'complete_width_family')
    untwisted_path=RESULTS/'twisted_chiral_untwisted_audit_2026_09_19.json'
    untwisted={row['label']:row for row in json.loads(untwisted_path.read_text())['records']}
    switching_path=RESULTS/'twisted_chiral_switching_audit_2026_09_19.json'
    switching={row['label']:row for row in json.loads(switching_path.read_text())['records']}
    rows=[];orders=[]
    for r in CLASS_COUNTS:
        available=[]
        for cls in range(CLASS_COUNTS[r]):
            f=families.get((r,cls));g=skew.get((r,cls))
            w=width.get((r,cls))
            u=untwisted[f'order{r}_class{cls}']
            so=switching[f'order{r}_class{cls}']
            row=dict(child_order=r,child_class=cls,child_cap=CHILD_CAP[r],parent_order=2*r,
                full_family_complete=f is not None,full_family_minimum=f['cap'] if f else None,
                skew_subfamily_minimum=g['cap'] if g else None,
                untwisted_best_matching_cap=u['cap'],
                untwisted_source=str(untwisted_path.relative_to(ROOT)),
                switching_only_minimum=so['cap'],
                switching_only_source=str(switching_path.relative_to(ROOT)),
                permutations_indispensable_for_full_minimum=so['permutations_indispensable_for_full_minimum'])
            if w:
                row.update(minimum_conditional_profile_radius=w['profile_radius'],
                    profile_radius_source=w['result_source'],
                    profile_radius_interpretation='minimum over every B of max_t (hi[t]-lo[t])/2; arbitrary independent sector recentering, not a matching-family cap')
            if f:
                available.append(f['cap']);row.update(result_source=f['result_source'],
                    witness_sha256=f['matrix_sha256'],distinct_B_count=f.get('full_B_count',f['count']),
                    evaluated_bridge_representatives=f['count'],
                    quotient_representatives=(f['count'] if f.get('signed_aut_antiaut_group_size',0)>0 else None),
                    enumeration_mode=('signed symmetry quotient' if f.get('signed_aut_antiaut_group_size',0)>0 else 'all distinct B directly'),
                    elapsed_seconds=f['elapsed'])
            rows.append(row)
        orders.append(dict(child_order=r,parent_order=2*r,total_optimal_child_classes=CLASS_COUNTS[r],
            completed_child_classes=len(available),all_optimal_child_classes_covered=len(available)==CLASS_COUNTS[r],
            best_parent_cap_witnessed=min(available) if available else None,
            minimum_over_all_optimal_children=min(available) if len(available)==CLASS_COUNTS[r] else None,
            separately_known_global_parent_minimum=KNOWN_GLOBAL_PARENT.get(2*r)))
    replay=json.loads((RESULTS/f'{PREFIX}classify10_independent.json').read_text())
    assert replay['complete']
    alternative_path=RESULTS/f'{PREFIX}alternative_conference.json'
    alternative=json.loads(alternative_path.read_text())['records'][0]
    assert alternative['kind']=='complete_family' and alternative['cap']==40 and alternative['child_cap']==15
    result=dict(schema='twisted-chiral-finite-family-summary-v1',
        generated_utc=datetime.now(timezone.utc).isoformat(),normalization='Q=max |sum_(i<j) Aij xi xj|',
        evidence='Only complete finite family runs enter minima. Every witness independently replayed over all projective spins. Global parent optima beyond14 are not asserted.',
        class_coverage_sources=dict(order9=f'computations/results/{PREFIX}classify9.json',order10=f'computations/results/{PREFIX}classify10.json'),
        orders=orders,classes=rows,
        nonoptimal_conference_child=dict(child_order=10,child_cap=15,parent_order=20,parent_cap=40,
            fixed_child_family_minimum=40,source=str(alternative_path.relative_to(ROOT)),
            witness_sha256=alternative['matrix_sha256']),
        decisive_finite_statement='The minimum over all globally optimal order10 children and every allowed twist/matching is44, while an explicit nonoptimal conference child gives40. Thus selectable-optimal-child sufficiency fails at parent20, regardless of the unknown global M20.')
    (RESULTS/f'{PREFIX}summary.json').write_text(json.dumps(result,indent=2)+'\n')
    lines=['# Exact twisted-double finite-family table','',
        'Normalization: Q=max absolute upper-triangle quadratic energy. A class-family minimum is not a global parent certificate.','',
        '| Child order | Optimal-child class | Child Q | Parent order | Untwisted minimum | Full family minimum | Skew matching subclass minimum |',
        '|---|---|---|---|---|---|---|']
    for row in rows:
        lines.append('| {child_order} | {child_class} | {child_cap} | {parent_order} | {untwisted_best_matching_cap} | {full} | {skew} |'.format(
            **row,full=row['full_family_minimum'] if row['full_family_complete'] else 'pending',
            skew=row['skew_subfamily_minimum'] if row['skew_subfamily_minimum'] is not None else 'pending'))
    lines+=['','Optimal-child class lists are complete through10. Both order9 and order10 censuses were replayed with independent signed-half-edge energy arithmetic.','',
        'Untwisted means B=A, still optimizing every matching d. It was independently checked by direct integer matrix products over all projective parent states. Twists improve25 to21 for two order7 classes and32 to30 for both order8 classes.','',
        'A nonoptimal conference child of order10 has child Q=15 and exact twisted-family minimum40, attained by B=A and d=all+. In contrast every optimal order10 child (Q=13) has family minimum44. This disproves selectable-optimal-child sufficiency at parent20 without claiming a global value for M20.','',
        'See `computations/results/twisted_chiral_2026_09_18_summary.json` for source pointers, exact state counts, and completeness flags.','']
    (ROOT/'artifacts/twisted_chiral_computation_table_2026_09_18.md').write_text('\n'.join(lines))
    print(json.dumps(orders,indent=2))


if __name__=='__main__':main()
