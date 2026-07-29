# Scale transfer by insertion profiles: an exact closure theorem and a finite no-go

## Setup and normalization

For a symmetric zero-diagonal signing \(A=(a_{ij})_{i,j\le n}\), write

\[
H_A(x)=\sum_{1\le i<j\le n}a_{ij}x_ix_j,
\qquad
M(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|,
\qquad
M_n=\min_A M(A).
\]

If a new vertex has signed row \(b\in\{\pm1\}^n\), the extended
Hamiltonian is

\[
H_C(x,y)=H_A(x)+y\,b\cdot x.
\]

For fixed \(x\),

\[
\max_{y=\pm1}|H_A(x)+y\,b\cdot x|
=|H_A(x)|+|b\cdot x|.
\]

Switching the old vertices by \(b\) gauges the new row to \(\mathbf1\).
As the switched old signing still ranges over all signings,

\[
\boxed{
M_{n+1}
=
\min_A\max_x
\left(
|H_A(x)|+\left|\sum_i x_i\right|
\right).
}
\]

This verifies the affine recurrence with the convention that \(H_A\)
contains each edge once. If one instead writes \(x^\top Ax=2H_A(x)\),
both quadratic terms must be doubled.

## The magnetization-extrema profile is not closed

Define

\[
U_A(m)=\max_{\sum_i x_i=m}H_A(x),
\qquad
L_A(m)=\min_{\sum_i x_i=m}H_A(x).
\]

This pair computes the canonical all-\(+1\) extension, but it does not
compute the best signed-row extension

\[
E(A)=
\min_{b\in\{\pm1\}^n}
\max_x\bigl(|H_A(x)|+|b\cdot x|\bigr).
\]

Here is an exact counterexample at \(n=7\):

\[
A_1=
\begin{pmatrix}
0&1&1&1&1&1&1\\
1&0&-1&1&-1&1&1\\
1&-1&0&1&-1&1&-1\\
1&1&1&0&1&-1&-1\\
1&-1&-1&1&0&-1&-1\\
1&1&1&-1&-1&0&-1\\
1&1&-1&-1&-1&-1&0
\end{pmatrix},
\]

\[
A_2=
\begin{pmatrix}
0&1&1&1&1&1&1\\
1&0&-1&-1&1&1&1\\
1&-1&0&1&-1&1&-1\\
1&-1&1&0&1&-1&-1\\
1&1&-1&1&0&-1&-1\\
1&1&1&-1&-1&0&-1\\
1&1&-1&-1&-1&-1&0
\end{pmatrix}.
\]

Their complete signed extrema tables agree:

| \(m\) | \(L_{A_1}(m)=L_{A_2}(m)\) | \(U_{A_1}(m)=U_{A_2}(m)\) | number of spins |
|---:|---:|---:|---:|
| \(-7\) | \(3\) | \(3\) | \(1\) |
| \(-5\) | \(-9\) | \(7\) | \(7\) |
| \(-3\) | \(-9\) | \(7\) | \(21\) |
| \(-1\) | \(-9\) | \(7\) | \(35\) |
| \(1\) | \(-9\) | \(7\) | \(35\) |
| \(3\) | \(-9\) | \(7\) | \(21\) |
| \(5\) | \(-9\) | \(7\) | \(7\) |
| \(7\) | \(3\) | \(3\) | \(1\) |

Nevertheless,

\[
\boxed{E(A_1)=12,\qquad E(A_2)=10.}
\]

More precisely, over all \(128\) rows \(b\), the histograms of

\[
R_A(b)=\max_x\bigl(|H_A(x)|+|b\cdot x|\bigr)
\]

are

\[
\begin{array}{c|rrrr}
 &10&12&14&16\\ \hline
A_1&0&42&72&14\\
A_2&2&62&56&8.
\end{array}
\]

For \(A_2\), an optimal row is

\[
b=(-1,1,-1,-1,1,1,-1),
\]

or its negative. Direct evaluation gives \(R_{A_2}(b)=10\).
Conversely \(M(A_2)=9\), and \(b\cdot x\) is always odd, so every row
has \(R_{A_2}(b)\ge10\). The \(A_1\) lower bound follows from the finite
enumeration reproduced below.

Thus two order-\(7\) minimizers with the same full signed
magnetization-extrema profile behave differently under optimal
insertion: one reaches the known \(M_8=10\), while the other cannot do
better than \(12\).

## The exact closed profile

For \(\sigma\in\{\pm1\}\) and \(h\in\mathbb R^n\), define the
external-field support functions

\[
F_A^\sigma(h)
=
\max_{x\in\{\pm1\}^n}
\bigl(\sigma H_A(x)+h\cdot x\bigr).
\]

They obey exact recursions.

### Switching

For \(s\in\{\pm1\}^n\), let \(A^s=\operatorname{diag}(s)A
\operatorname{diag}(s)\). Then

\[
\boxed{F_{A^s}^\sigma(h)=F_A^\sigma(s\odot h).}
\]

### Vertex insertion

For

\[
C=\begin{pmatrix}A&b\\b^\top&0\end{pmatrix},
\]

an external field on the new system is written \((h,t)\). Maximizing
first over its last spin \(y\) gives

\[
\boxed{
F_C^\sigma(h,t)
=
\max_{y=\pm1}
\left[
ty+F_A^\sigma(h+\sigma yb)
\right].
}
\]

Hence the pair \((F_A^+,F_A^-)\) is closed under arbitrary switching and
vertex insertion.

## Closure is full information

The closed support profile is not a lower-dimensional compression. For
every cube vertex \(x_0\),

\[
\boxed{
H_A(x_0)
=
\lim_{\lambda\to\infty}
\left(F_A^+(\lambda x_0)-\lambda n\right).
}
\]

Indeed, the field term \(\lambda x_0\cdot x\) loses
\(2\lambda d_H(x,x_0)\) away from \(x_0\). Once

\[
\lambda>
\max_{x\ne x_0}
\frac{H_A(x)-H_A(x_0)}{2d_H(x,x_0)},
\]

\(x_0\) is the unique maximizer.

Equivalently, for the Fenchel conjugate

\[
(F_A^+)^*(z)=\sup_h\bigl(h\cdot z-F_A^+(h)\bigr),
\]

one has

\[
\boxed{(F_A^+)^*(x)=-H_A(x)\quad\text{for every }x\in\{\pm1\}^n.}
\]

Therefore the exact convex/Legendre profile contains the entire
\(2^{n-1}\)-entry switching-energy word. The one-dimensional radial
restriction \(h=t\mathbf1\) is too small: \(A_1,A_2\) above have
identical \(F_A^\pm(t\mathbf1)\) for every real \(t\), but different
optimal extensions.

This stops the proposed exact profile-dynamic-programming route. A
scale-transfer proof needs a genuinely lossy asymptotic theorem, such
as a structural discrepancy bound for the high-energy layers, rather
than an exact finite profile recursion.

## Compact independent enumeration

```python
from itertools import product

A1 = [
 [0,1,1,1,1,1,1], [1,0,-1,1,-1,1,1],
 [1,-1,0,1,-1,1,-1], [1,1,1,0,1,-1,-1],
 [1,-1,-1,1,0,-1,-1], [1,1,1,-1,-1,0,-1],
 [1,1,-1,-1,-1,-1,0],
]
A2 = [
 [0,1,1,1,1,1,1], [1,0,-1,-1,1,1,1],
 [1,-1,0,1,-1,1,-1], [1,-1,1,0,1,-1,-1],
 [1,1,-1,1,0,-1,-1], [1,1,1,-1,-1,0,-1],
 [1,1,-1,-1,-1,-1,0],
]

X = list(product((-1, 1), repeat=7))

def energy(A, x):
    return sum(A[i][j]*x[i]*x[j]
               for i in range(7) for j in range(i+1, 7))

def audit(A):
    H = [energy(A, x) for x in X]
    profile = {
        m: (min(h for h,x in zip(H,X) if sum(x) == m),
            max(h for h,x in zip(H,X) if sum(x) == m))
        for m in range(-7, 8, 2)
    }
    values = []
    for b in X:
        values.append(max(abs(h) + abs(sum(bi*xi for bi,xi in zip(b,x)))
                          for h,x in zip(H,X)))
    histogram = {r: values.count(r) for r in sorted(set(values))}
    return profile, min(values), histogram

print(audit(A1))
print(audit(A2))
```
