# Exact order-25 universal-row replenishment tower

Checkpoint date: 2026-07-25.

This file records a completely explicit finite tower for the max-plus
universal-row insertion mechanism. It is a certificate for the statement
that the hard suffix follow-the-leader replenishment gap can already exceed
the terminal Boolean quadratic norm by a factor \(3/2\) at order \(21\).
It is **not** an asymptotic construction.

## Definitions

For a symmetric zero-diagonal sign matrix \(A_m\), write

\[
E_m(z)=z^\top A_m z,\qquad
q_m=Q(A_m)=\max_{z\in\{\pm1\}^m}|E_m(z)|.
\]

Given \(x_m\in\{\pm1\}^m\), append a new vertex with incident row \(x_m\):

\[
A_{m+1}=\begin{pmatrix}A_m&x_m\\x_m^\top&0\end{pmatrix}.
\]

The listed row \(x_m\) satisfies the exact max-plus fixed-point condition

\[
E_m(x_m)=
\max_z\bigl(|E_m(z)|-4\delta_H(x_m,z)\bigr),
\qquad
\delta_H(x,z)=\min(d_H(x,z),m-d_H(x,z)).
\]

Consequently \((x_m,1)\) is a ground state of \(A_{m+1}\) and

\[
q_{m+1}=E_m(x_m)+2m,\quad
g_m=q_m-E_m(x_m),\quad
d_m=q_{m+1}-q_m=2m-g_m.
\]

Deleting vertices in reverse insertion order, \(25,24,\ldots,10\), realizes
the listed replenishment gaps \(g_{24},g_{23},\ldots,g_9\).

## Exact history

All energies use the doubled convention \(z^\top A z\), so they are twice
the original half-sum objective.

| \(m\) | \(q_m\) | \(E_m(x_m)\) | \(g_m\) | \(d_m\) | \(q_{m+1}\) | \(\sum_{k=9}^m g_k\) | \(q_m/m^{3/2}\) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 9 | 32 | 16 | 16 | 2 | 34 | 16 | 1.185185185 |
| 10 | 34 | 18 | 16 | 4 | 38 | 32 | 1.075174404 |
| 11 | 38 | 30 | 8 | 14 | 52 | 40 | 1.041584645 |
| 12 | 52 | 32 | 20 | 4 | 56 | 60 | 1.250925583 |
| 13 | 56 | 40 | 16 | 10 | 66 | 76 | 1.194738884 |
| 14 | 66 | 50 | 16 | 12 | 78 | 92 | 1.259945855 |
| 15 | 78 | 62 | 16 | 14 | 92 | 108 | 1.342634227 |
| 16 | 92 | 68 | 24 | 8 | 100 | 132 | 1.437500000 |
| 17 | 100 | 76 | 24 | 10 | 110 | 156 | 1.426680147 |
| 18 | 110 | 86 | 24 | 12 | 122 | 180 | 1.440402702 |
| 19 | 122 | 106 | 16 | 22 | 144 | 196 | 1.473090502 |
| 20 | 144 | 112 | 32 | 8 | 152 | 228 | 1.609968944 |
| 21 | 152 | 136 | 16 | 26 | 178 | 244 | 1.579481872 |
| 22 | 178 | 146 | 32 | 12 | 190 | 276 | 1.724987614 |
| 23 | 190 | 166 | 24 | 22 | 212 | 300 | 1.722510377 |
| 24 | 212 | 180 | 32 | 16 | 228 | 332 | 1.803096616 |

Thus

\[
q_{25}=228,\qquad
\sum_{m=9}^{24}g_m=332,\qquad
\frac{\sum g_m}{q_{25}}=\frac{83}{57}=1.456140\ldots .
\]

The maximum prefix ratio occurs after the \(m=20\) insertion:

\[
\frac{\sum_{m=9}^{20}g_m}{q_{21}}
=\frac{228}{152}=\frac32.
\]

## Appended rows

Row \(m\) below is \(x_m\), in the coordinate order \(1,\ldots,m\).

```text
m=9:  1 1 -1 1 -1 -1 -1 -1 -1
m=10: -1 -1 -1 1 1 -1 1 -1 1 -1
m=11: 1 -1 1 -1 -1 1 -1 1 -1 -1 -1
m=12: -1 1 1 -1 1 1 1 1 1 -1 1 -1
m=13: 1 -1 -1 1 -1 1 1 1 1 1 -1 1 -1
m=14: 1 1 1 1 -1 -1 -1 -1 -1 1 -1 -1 -1 -1
m=15: -1 1 1 -1 1 1 -1 -1 -1 -1 1 1 1 -1 -1
m=16: -1 -1 -1 1 1 -1 1 -1 1 -1 1 -1 1 -1 -1 -1
m=17: -1 1 1 1 1 -1 -1 -1 -1 1 1 -1 -1 -1 1 1 -1
m=18: 1 1 1 -1 -1 1 -1 1 -1 -1 -1 1 1 -1 1 1 -1 -1
m=19: 1 -1 -1 1 -1 -1 1 1 1 1 -1 1 -1 1 1 -1 -1 -1 -1
m=20: 1 1 1 1 -1 -1 -1 -1 -1 1 -1 -1 -1 -1 1 -1 -1 1 1 -1
m=21: -1 -1 -1 1 1 -1 1 -1 1 1 1 -1 -1 1 -1 -1 1 1 -1 1 -1
m=22: 1 -1 -1 -1 -1 1 1 1 1 -1 -1 1 -1 1 -1 -1 -1 -1 1 1 -1 -1
m=23: 1 1 -1 1 -1 -1 -1 -1 -1 1 -1 -1 -1 1 1 -1 -1 1 -1 1 1 1 -1
m=24: -1 1 -1 -1 1 1 1 -1 1 -1 1 -1 1 -1 -1 1 1 1 -1 -1 -1 1 -1 -1
```

## Final \(25\times25\) matrix

The seed occupies vertices \(1,\ldots,9\); vertices \(10,\ldots,25\) are
the appended rows above.

```text
 0 -1  1  1 -1 -1 -1  1  1  1 -1  1 -1  1  1 -1 -1 -1  1  1  1 -1  1  1 -1
-1  0  1  1  1  1 -1 -1 -1  1 -1 -1  1 -1  1  1 -1  1  1 -1  1 -1 -1  1  1
 1  1  0 -1  1  1 -1  1 -1 -1 -1  1  1 -1  1  1 -1  1  1 -1  1 -1 -1 -1 -1
 1  1 -1  0 -1 -1  1 -1  1  1  1 -1 -1  1  1 -1  1  1 -1  1  1  1 -1  1 -1
-1  1  1 -1  0 -1 -1 -1  1 -1  1 -1  1 -1 -1  1  1  1 -1 -1 -1  1 -1 -1  1
-1  1  1 -1 -1  0  1  1  1 -1 -1  1  1  1 -1  1 -1 -1  1 -1 -1 -1  1 -1  1
-1 -1 -1  1 -1  1  0  1  1 -1  1 -1  1  1 -1 -1  1 -1 -1  1 -1  1  1 -1  1
 1 -1  1 -1 -1  1  1  0  1 -1 -1  1  1  1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1
 1 -1 -1  1  1  1  1  1  0 -1  1 -1  1  1 -1 -1  1 -1 -1  1 -1  1  1 -1  1
 1  1 -1  1 -1 -1 -1 -1 -1  0 -1 -1 -1  1  1 -1 -1  1 -1  1  1  1 -1  1 -1
-1 -1 -1  1  1 -1  1 -1  1 -1  0 -1  1 -1 -1  1  1  1 -1 -1 -1  1 -1 -1  1
 1 -1  1 -1 -1  1 -1  1 -1 -1 -1  0 -1  1 -1  1 -1 -1  1  1 -1 -1  1 -1 -1
-1  1  1 -1  1  1  1  1  1 -1  1 -1  0 -1 -1  1  1 -1  1 -1 -1 -1 -1 -1  1
 1 -1 -1  1 -1  1  1  1  1  1 -1  1 -1  0 -1 -1 -1 -1 -1  1 -1  1  1  1 -1
 1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  0 -1 -1  1  1  1  1 -1 -1  1 -1
-1  1  1 -1  1  1 -1 -1 -1 -1  1  1  1 -1 -1  0 -1  1  1 -1 -1 -1 -1 -1  1
-1 -1 -1  1  1 -1  1 -1  1 -1  1 -1  1 -1 -1 -1  0 -1 -1 -1 -1  1 -1 -1  1
-1  1  1  1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1  1 -1  0 -1 -1  1  1 -1  1  1
 1  1  1 -1 -1  1 -1  1 -1 -1 -1  1  1 -1  1  1 -1 -1  0 -1  1 -1  1 -1 -1
 1 -1 -1  1 -1 -1  1  1  1  1 -1  1 -1  1  1 -1 -1 -1 -1  0 -1  1  1  1 -1
 1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1  1 -1  0 -1 -1  1 -1
-1 -1 -1  1  1 -1  1 -1  1  1  1 -1 -1  1 -1 -1  1  1 -1  1 -1  0 -1  1  1
 1 -1 -1 -1 -1  1  1  1  1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1  0 -1 -1
 1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1  1 -1  1  1  1 -1  0 -1
-1  1 -1 -1  1  1  1 -1  1 -1  1 -1  1 -1 -1  1  1  1 -1 -1 -1  1 -1 -1  0
```

## Verification algorithm

For each prefix \(A_m\):

1. enumerate all \(2^m\) sign vectors and compute \(E_m(z)\);
2. initialize \(F_0(z)=|E_m(z)|\);
3. for each hypercube coordinate \(b\), apply the exact distance transform
   \[
   F_{b+1}(z)=\max(F_b(z),F_b(z\oplus e_b)-4);
   \]
4. after all coordinates, check \(F_m(x_m)=E_m(x_m)\);
5. independently enumerate \(A_{m+1}\) and check
   \(Q(A_{m+1})=E_m(x_m)+2m\).

Antipodal symmetry makes ordinary Hamming distance in step 3 equivalent to
\(\delta_H\), because \(|E_m(z)|=|E_m(-z)|\). The search and verification
used exact integer arithmetic throughout.

## Interpretation and limitation

The identity

\[
\sum_{m=9}^{N-1}g_m=(N-1)N-q_N-\bigl(72-q_9\bigr)
\]

specializes here to \(332\) at \(N=25\). Therefore an **infinite** tower of
this type with \(q_N=o(N^2)\) would make \(\sum g_m/q_N\to\infty\), disproving
any universal \(O(Q)\) bound on hard replenishment.

This finite tower does not establish that. Its normalized norm rises to
\(q_{25}/25^{3/2}=1.824\), already well above the known asymptotic
conference benchmark \(Q/n^{3/2}\le1+o(1)\). The unresolved task is to make
the insertion rule extend indefinitely while keeping \(q_N=O(N^{3/2})\).
