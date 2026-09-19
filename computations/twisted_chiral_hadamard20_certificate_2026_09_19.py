#!/usr/bin/env python3
"""Standalone exact cap-40 certificate, using only Python's standard library.

The upper bound is algebraic, not an enumeration: for symmetric H^2=20I,
put r_i=z_i(Hz)_i. Orthogonal +/-1 rows imply one common residue 0 or 2
modulo4, and sum r_i^2=400. In residue2, (r-2)(r-6)>=0 gives sum r<=80.
In residue0, (r-4)(r-8)>=0 gives sum r<=86+2/3, hence at most84.
At84, w=(Hz-4z)/4 is integral, ||w||_2^2=3, z.w=1, and Hw=z-4w.
Thus w^THw=-11, but |w^THw|<=||w||_1^2=9, a contradiction. Applying
the same argument to -H gives |z^THz|<=80. Here traceH=0, so deleting
its diagonal leaves quadratic cap<=40, attained by the displayed spin.

Canonical proof: artifacts/twisted_chiral_hadamard20_2026_09_19.md.
No external input, packages, solver or exhaustive parent scan is required.
"""
import json


A=[
    [0,1,1,1,1,1,1,1,1,1],
    [1,0,1,1,1,-1,-1,1,-1,-1],
    [1,1,0,1,-1,1,-1,-1,1,-1],
    [1,1,1,0,-1,-1,1,-1,-1,1],
    [1,1,-1,-1,0,1,1,1,-1,-1],
    [1,-1,1,-1,1,0,1,-1,1,-1],
    [1,-1,-1,1,1,1,0,-1,-1,1],
    [1,1,-1,-1,1,-1,-1,0,1,1],
    [1,-1,1,-1,-1,1,-1,1,0,1],
    [1,-1,-1,1,-1,-1,1,1,1,0],
]
D_MATCHING=[1]*10
Z=[1,-1,-1,-1,1,1,1,1,1,1,-1,-1,-1,1,1,1,1,1,1,1]


def main():
    assert all(A[i][j]==A[j][i] for i in range(10) for j in range(10))
    assert all(A[i][i]==0 for i in range(10))
    assert all(sum(A[i][k]*A[k][j] for k in range(10))==9*(i==j)
               for i in range(10) for j in range(10))
    d=[[0]*20 for _ in range(20)]
    for i in range(10):
        for j in range(10):
            d[i][j]=A[i][j];d[i+10][j+10]=-A[i][j]
            d[i][j+10]=d[j+10][i]=A[i][j]+D_MATCHING[i]*(i==j)
    h=[row[:] for row in d]
    for i in range(20):h[i][i]=-1 if i<10 else 1
    assert all(h[i][j] in [-1,1] and h[i][j]==h[j][i]
               for i in range(20) for j in range(20))
    assert all(sum(h[i][k]*h[k][j] for k in range(20))==20*(i==j)
               for i in range(20) for j in range(20))
    assert sum(h[i][i] for i in range(20))==0
    row_sums=list(map(sum,h))
    switched_rows=[Z[i]*sum(h[i][j]*Z[j] for j in range(20)) for i in range(20)]
    assert row_sums==[18]+[2]*19
    assert sorted(switched_rows)==[2]*10+[6]*10
    energy=sum(d[i][j]*Z[i]*Z[j] for i in range(20) for j in range(i+1,20))
    assert energy==40 and sum(switched_rows)==80
    # Exact scalar inequalities used in the uniform algebraic upper proof.
    assert all((r-2)*(r-6)>=0 for r in range(-18,19,4))
    assert all((r-4)*(r-8)>=0 for r in range(-20,21,4))
    assert (400+12*20)//8==80
    assert 4*((400+32*20)//48)==84
    assert (400-8*84+16*20)//16==3
    assert (84-4*20)//4==1
    assert abs(1-4*3)>3**2
    result=dict(conference_child=A,matching_d=D_MATCHING,maximizing_spin=Z,
        H_definition='D+diag(-I10,+I10), where D=[[A,A+I],[A+I,-A]]',
        H_squared_equals_20I=True,H_trace=0,H_row_sums=row_sums,
        switched_H_row_sums=switched_rows,attained_energy=energy,
        exact_parent_cap=40,upper_bound_method='algebraic symmetric-Hadamard20 bound, not parent spin enumeration',
        maximizing_spin_source='stored independent Gray argmax3591; bit i-1 denotes negative spin i, spin0=+1',
        canonical_proof='artifacts/twisted_chiral_hadamard20_2026_09_19.md')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
