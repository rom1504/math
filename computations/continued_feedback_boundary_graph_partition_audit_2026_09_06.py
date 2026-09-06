"""Exact combinatorial stress test for marked boundary coverage.

This is not a proof by finite enumeration: graph sizes are bounded here.
It independently verifies every vertex-equality pattern for named source
monomials, including root identifications, loops, and retained multiedges.
"""
import argparse
import itertools
import json
from pathlib import Path


def partitions(n):
    labels=[0]*n
    def visit(position,maximum):
        if position==n:
            yield tuple(labels)
            return
        for label in range(maximum+2):
            labels[position]=label
            yield from visit(position+1,max(maximum,label))
    yield from visit(1,0)


def source_diagram(word):
    marks=[0];edges=[];stars=[]
    for letter in word:
        if letter=='S':marks[0]+=1
        elif letter=='G':
            leaf=len(marks);marks.append(1);edges.append((0,leaf))
        elif letter=='Y':
            center=len(marks);marks.extend([1,1,1])
            edges.extend([(0,center),(center,center+1),(center,center+2)])
            stars.append((center,center+1,center+2))
        else:raise ValueError(letter)
    return marks,edges,stars


def components(vertices,edges,skip=None):
    parent=list(range(vertices))
    def find(a):
        while parent[a]!=a:
            parent[a]=parent[parent[a]];a=parent[a]
        return a
    for index,(a,b) in enumerate(edges):
        if index!=skip:parent[find(a)]=find(b)
    return [find(a) for a in range(vertices)]


def audit(word):
    original_marks,original_edges,original_stars=source_diagram(word)
    records={'word':word,'original_vertices':len(original_marks),'partitions':0,
        'high_degree_patterns':0,'degree_two_patterns':0,'positive_global_patterns':0,
        'exceptional_cubic_patterns':0,
        'failures':[]}
    for labels in partitions(len(original_marks)):
        records['partitions']+=1
        size=max(labels)+1;marks=[0]*size
        for index,count in enumerate(original_marks):marks[labels[index]]+=count
        edges=[(labels[a],labels[b]) for a,b in original_edges]
        degree=[0]*size
        for a,b in edges:degree[a]+=1;degree[b]+=1
        root=labels[0];free={i for i,v in enumerate(marks) if v%2};q=len(free)
        if any((degree[i]-marks[i])%2 for i in range(size) if i!=root):
            raise AssertionError(('nonroot_parity',word,labels))
        if (degree[root]-marks[root]-len(word))%2:
            raise AssertionError(('root_parity',word,labels))
        bridges=[]
        for index,(a,b) in enumerate(edges):
            comp=components(size,edges,index)
            if comp[a]!=comp[b]:bridges.append(index)
        nonbridges=[edge for index,edge in enumerate(edges) if index not in bridges]
        comp=components(size,nonbridges);forest_degree={c:0 for c in comp}
        for index in bridges:
            a,b=edges[index];forest_degree[comp[a]]+=1;forest_degree[comp[b]]+=1
            full=components(size,edges,index)
            away={i for i in range(size) if full[i]!=full[root]}
            if len(away&free)>3:
                raise AssertionError(('bridge_three_mark_bound',word,labels,index))
        leaves={c for c,d in forest_degree.items() if d<=1}
        covered={comp[i] for i in free}
        if q==3 and not leaves<=covered:
            records['exceptional_cubic_patterns']+=1
            assert leaves-covered=={comp[root]}
            assert forest_degree[comp[root]]==1 and not free&{i for i in range(size) if comp[i]==comp[root]}
            separating=[index for index in bridges if comp[root] in (comp[edges[index][0]],comp[edges[index][1]])]
            assert len(separating)==1
            side_components=components(size,edges,separating[0])
            away={i for i in range(size) if side_components[i]!=side_components[root]}
            assert away==free and len(away)==3
            matching=[star for star in original_stars if {labels[i] for i in star}==away]
            assert len(matching)==1
        if q>3:
            records['high_degree_patterns']+=1
            if not leaves<=covered:records['failures'].append(['high_degree',labels])
        if q==2:
            records['degree_two_patterns']+=1
            if not leaves<=covered:records['failures'].append(['degree_two',labels])
        if q>0:
            records['positive_global_patterns']+=1
            if not leaves<=covered|{comp[root]}:records['failures'].append(['global',labels])
    assert not records['failures'],records
    return records


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--words',default='YYY,YYG,YGG,GGG,SYY,SSY,SSS,YY,SGY')
    args=parser.parse_args()
    results=[]
    for word in args.words.split(','):
        record=audit(word);results.append(record);print(json.dumps(record),flush=True)
    result={'records':results,'status':'exact finite equality-pattern audit PASS; general proof remains graph argument'}
    args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
