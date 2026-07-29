# Regular peeling and replenishment gaps

## Setup

For a symmetric zero-diagonal sign matrix \(A\), write

\[
E_A(x)=x^{T}Ax,\qquad Q(A)=\max_{x\in\{\pm1\}^n}|E_A(x)|.
\]

Fix the vertex order \(1,\ldots,n\), and put

\[
V_i=\{i,i+1,\ldots,n\},\qquad A_i=A[V_i],\qquad q_i=Q(A_i).
\]

For every \(i<n\), choose an absolute maximizer \(x^{(i)}\) of \(A_i\).
Let

\[
\sigma_i=\operatorname{sgn}E_{A_i}(x^{(i)})
\]

(either sign may be used if the energy is zero), so that

\[
\sigma_iE_{A_i}(x^{(i)})=q_i.
\]

The oriented energy left after deleting vertex \(i\) is

\[
e_i=
\sigma_i E_{A_{i+1}}\!\left(x^{(i)}|_{V_{i+1}}\right),
\]

and the replenishment gap is

\[
g_i=q_{i+1}-e_i\ge0.
\]

When a suffix has several absolute maximizers, the certificate below uses
the one maximizing \(g_i\).

## FTL/mosaic identity

Encode the \(i\)-th suffix leader by its augmented cut

\[
Z^{(i)}_{jk}=\sigma_i x^{(i)}_j x^{(i)}_k,
\qquad j,k\in V_i.
\]

For an augmented cut \(Z\), define the \(i\)-th star reward

\[
f_i(Z)=2\sum_{j>i}a_{ij}Z_{ij}.
\]

Then

\[
q_i=\sum_{k\ge i}f_k(Z^{(i)}),
\qquad
e_i=\sum_{k\ge i+1}f_k(Z^{(i)}),
\]

and consequently

\[
f_i(Z^{(i)})=q_i-e_i=(q_i-q_{i+1})+g_i.
\]

Summing and using \(q_n=0\) gives the exact identity

\[
\boxed{
\sum_{i=1}^{n-1}g_i
=
\sum_{i=1}^{n-1}f_i(Z^{(i)})-q_1.
}
\]

Equivalently, define the triangular FTL mosaic

\[
W_{ij}=Z^{(i)}_{ij}\qquad(i<j).
\]

Then

\[
\boxed{
\sum_i g_i=2\sum_{i<j}a_{ij}W_{ij}-Q(A).
}
\]

Thus bounding cumulative replenishment is exactly an adaptivity-gap
problem: compare the payoff of suffix follow-the-leader, whose \(i\)-th
row may use a different augmented cut, with the best single augmented cut.
Without the suffix-optimality condition the mosaic \(W\) can be an
arbitrary edge signing, so triangle identities alone cannot bound it.

## Exact counterexample to \(\sum g_i\le2Q(A)\)

The following \(15\times15\) sign matrix uses the natural deletion order.

```text
 0 -1  1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1  1
-1  0 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1  1  1
 1 -1  0 -1  1  1 -1  1  1 -1 -1 -1 -1 -1  1
 1 -1 -1  0 -1 -1 -1 -1  1 -1  1 -1  1  1 -1
-1 -1  1 -1  0  1  1  1  1  1  1 -1  1 -1 -1
 1  1  1 -1  1  0 -1 -1 -1  1  1 -1 -1  1 -1
 1 -1 -1 -1  1 -1  0 -1  1  1 -1  1  1  1 -1
-1 -1  1 -1  1 -1 -1  0 -1 -1 -1 -1  1  1  1
-1 -1  1  1  1 -1  1 -1  0 -1 -1  1 -1 -1  1
-1 -1 -1 -1  1  1  1 -1 -1  0 -1  1  1  1 -1
-1 -1 -1  1  1  1 -1 -1 -1 -1  0 -1 -1 -1 -1
-1  1 -1 -1 -1 -1  1 -1  1  1 -1  0 -1 -1 -1
 1  1 -1  1  1 -1  1  1 -1  1 -1 -1  0  1  1
-1  1 -1  1 -1  1  1  1 -1  1 -1 -1  1  0 -1
 1  1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1  0
```

Exact suffix enumeration gives:

| suffix order \(m\) | \(q_m\) | next \(q_{m-1}\) | maximum \(g\) over leader ties | old restricted energy \(e\) |
|---:|---:|---:|---:|---:|
| 15 | 62 | 58 | 20 | 38 |
| 14 | 58 | 52 | 16 | 36 |
| 13 | 52 | 44 | 12 | 32 |
| 12 | 44 | 42 | 16 | 26 |
| 11 | 42 | 30 | 8 | 22 |
| 10 | 30 | 28 | 16 | 12 |
| 9 | 28 | 28 | 16 | 12 |
| 8 | 28 | 18 | 0 | 18 |
| 7 | 18 | 14 | 8 | 6 |
| 6 | 14 | 12 | 8 | 4 |
| 5 | 12 | 8 | 4 | 4 |
| 4 | 8 | 6 | 4 | 2 |
| 3 | 6 | 2 | 0 | 2 |
| 2 | 2 | 0 | 0 | 0 |

Therefore

\[
Q(A)=62,\qquad
\sum_i g_i=128,
\]

and

\[
\boxed{
\sum_i g_i=128>124=2Q(A).
}
\]

So the conjectured universal inequality

\[
\sum_i g_i\le2Q(A)
\]

is false, even for a dense sign matrix and singleton suffix deletion.

## Verification method

For each suffix of order \(m\):

1. Enumerate all \(2^m\) Boolean vectors.
2. Compute \(x^TA_ix\) in exact integer arithmetic and set
   \(q_i=\max_x|x^TA_ix|\).
3. Enumerate every tie with \(|x^TA_ix|=q_i\).
4. Orient it by \(\sigma=\operatorname{sgn}(x^TA_ix)\).
5. Delete its first coordinate, compute the exact oriented core energy
   \(e_i\), and maximize \(q_{i+1}-e_i\) over ties.

The table and total were reproduced independently by:

- a C++ bit-enumeration evaluator; and
- a Python/NumPy exhaustive evaluator.

All energies have the required parity.

## Surviving replacement question

The false factor \(2\) statement should be replaced by:

> Is the suffix-FTL adaptivity gap
> \[
> \frac{\sum_i g_i}{Q(A)}
> \]
> bounded by an absolute constant over all dense sign matrices, orders,
> and choices among suffix ground-state ties?

A positive answer with any absolute constant would close the cumulative-gap
part of the peeling estimate. A scalable family with a diverging ratio
would show that hard ground-state peeling must be replaced by an
entropy-regularized or otherwise averaged procedure.

## Cumulative visibility theorem

There is nevertheless a cumulative theorem for the part of each deleted
block visible to the *successor* ground-state layer.

Write a peeling step as

\[
A_t=
\begin{pmatrix}
D_t&B_t\\
B_t^T&C_t
\end{pmatrix},
\qquad C_t=A_{t+1},
\qquad d_t=Q(A_t)-Q(C_t).
\]

Define

\[
V_t=
\max_{y\in\operatorname{GS}(C_t)}
\|B_ty\|_1,
\]

where \(\operatorname{GS}(C_t)\) contains the Boolean vectors attaining
either sign of \(Q(C_t)\). Then

\[
\boxed{
2V_t\le d_t+Q(D_t).
}
\]

Indeed, orient \(C_t\) so a selected \(y\) has energy \(Q(C_t)\), and
choose \(z=\operatorname{sgn}(B_ty)\). Evaluating the two parent states
\((z,y)\) and \((-z,y)\) gives

\[
Q(A_t)
\ge
Q(C_t)-Q(D_t)+2\|B_ty\|_1.
\]

For disjoint peeled blocks,

\[
\sum_tQ(D_t)\le2Q(A_0),
\qquad
\sum_td_t\le Q(A_0),
\]

and hence

\[
\boxed{
\sum_tV_t\le\frac32Q(A_0).
}
\]

The same proof gives an energy-layer version. If

\[
s\,y^TC_ty\ge Q(C_t)-\delta
\]

for some orientation \(s\), then

\[
\boxed{
2\|B_ty\|_1
\le d_t+Q(D_t)+\delta.
}
\]

Therefore any tower with unbounded cumulative replenishment must obtain
its large fields from predecessor states that lie well below every
successor ground-state frame. Pointwise inversion is impossible: there
are regular finite examples with \(d_t=0\), \(V_t=0\), but \(g_t>0\).
The remaining obstruction is genuinely a succession of different,
ground-state-invisible energy layers.

## Universal-row insertion criterion

There is a concrete mechanism that could make the hard-FTL ratio
unbounded.

Let \(C\) have order \(m\), let

\[
E(z)=z^TCz,
\]

and adjoin a vertex whose row is a Boolean vector \(x\):

\[
A_x=
\begin{pmatrix}
0&x^T\\
x&C
\end{pmatrix}.
\]

The exact affine identity is

\[
Q(A_x)
=
\max_z\left(|E(z)|+2|x^Tz|\right).
\]

Writing

\[
\delta_H(x,z)=\min\{d_H(x,z),m-d_H(x,z)\},
\]

this becomes

\[
Q(A_x)
=
2m+\max_z\left(|E(z)|-4\delta_H(x,z)\right).
\]

Consequently, if

\[
\boxed{
E(x)=
\max_z\left(|E(z)|-4\delta_H(x,z)\right),
}
\]

then \((1,x)\) is a parent ground state, its new-vertex field is exactly
\(m\), and

\[
Q(A_x)=E(x)+2m.
\]

Let \(q_m=Q(C)\), \(d_m=Q(A_x)-q_m\), and let \(g_m\) be the
replenishment gap exposed by deleting the new vertex. Then

\[
\boxed{
d_m+g_m=2m.
}
\]

Thus an infinite recursive insertion tower satisfying the boxed
max-plus condition and \(q_m=o(m^2)\) would give

\[
\sum_{k<m}g_k
=
m(m-1)-q_m+O(1),
\]

so its hard-FTL adaptivity ratio would diverge. The max-plus condition can
be tested for all \(x\) in \(O(m2^m)\) operations by the hypercube
distance transform

\[
F(x)=\max_z\left(|E(z)|-4\delta_H(x,z)\right).
\]

Finite searches find many non-ground-state fixed points \(E(x)=F(x)\)
and recursive towers through order \(25\), but no proof currently shows
that a subquadratic tower continues indefinitely. This is the sharp
construction target for disproving every constant cumulative-gap bound.

## Exact cut and layer criteria for extending the tower

Switch \(C\) by a proposed row \(x\), so that the proposed row becomes
all \(+1\), and call the switched matrix \(B=D_xCD_x\). Put

\[
e=E_C(x)=\mathbf1^\top B\mathbf1
\]

and, for \(S\subset[m]\), let

\[
c_B(S)=\sum_{i\in S,\;j\notin S}b_{ij}.
\]

Flipping the coordinates in \(S\) changes the energy by

\[
E_B(\mathbf1^S)=e-4c_B(S).
\]

It follows, without loss restricting to \(|S|\le m/2\), that the
max-plus fixed-point condition is equivalent to the pair of cut
inequalities

\[
\boxed{
-|S|\le c_B(S)\le |S|+\frac e2
\qquad (|S|\le m/2).
}
\]

This is a concrete probabilistic construction target: find a switching
whose total excess is \(q_m-2m+O(\sqrt m)\), while every signed cut lies
inside the displayed asymmetric corridor.

There is also an exact near-ground-layer covering criterion. All energies
are congruent modulo \(4\). Define

\[
\mathcal P_r^+
=\{x:E_C(x)=q_m-4r\ge0\},
\qquad
\mathcal Z_s
=\{z:|E_C(z)|=q_m-4s\}.
\]

A vector \(x\in\mathcal P_r^+\) is a max-plus fixed point if and only if

\[
\boxed{
\delta_H(x,\mathcal Z_s)\ge r-s
\quad\text{for every }0\le s<r.
}
\]

Thus a sufficient (union-bound) extension criterion is

\[
|\mathcal P_r^+|
>
\sum_{s=0}^{r-1}
|\mathcal Z_s|\,
V_m(r-s-1),
\]

where

\[
V_m(R)=2\sum_{j=0}^R\binom mj
\]

is the size of a projective Hamming ball in the full cube when
\(R<m/2\). A sharper, still exact, version replaces \(V_m(R)\) with the
largest intersection of such a ball with \(\mathcal P_r^+\). For a
scale-preserving extension one needs

\[
r=\frac m2-O(\sqrt m),
\]

so this criterion exposes the real difficulty: a state only \(2m+O(\sqrt
m)\) below the optimum must lie simultaneously in equatorial bands around
every higher absolute-energy layer.

More generally, if the insertion has increment \(d_m\) and
\(|E_C(z)|\ge q_m-s\), fixed-point maximality gives

\[
\boxed{
|\langle x,z\rangle|
\le \frac{d_m+s}{2}.
}
\]

In particular, the ground state created by the preceding insertion is a
ground state of the current core, so, writing it \(u_m\),

\[
d_m\ge2|\langle x_m,u_m\rangle|.
\]

Telescoping yields the monotone obstruction

\[
\boxed{
\sum_{m=m_0}^{N-1}|\langle x_m,u_m\rangle|
\le \frac{q_N-q_{m_0}}2.
}
\]

Hence a tower with \(q_N=O(N^{3/2})\) must generate a new row nearly
orthogonal to the preceding ground state on all but a sparse set of
orders. Any rule having \(|\langle x_m,u_m\rangle|\ge\varepsilon m\) on
a positive density of steps necessarily pays \(q_N=\Omega(N^2)\).

This rules out the most literal periodic extension rules. More generally,
suppose an eventual construction has finitely many vertex types,
\(c_i\in[r]\), and

\[
a_{ij}=K(c_i,c_j)
\]

for a fixed symmetric sign kernel \(K\). On a subsequence, let \(p_a\) be
the limiting type densities. Spins constant on each type give

\[
\frac{Q(A_N)}{N^2}
\ \ge\
\max_{y\in\{\pm1\}^r}
\left|
\sum_{a,b}p_ap_bK_{ab}y_ay_b
\right|+o(1).
\]

The maximum is strictly positive: the square average over \(y\) contains
the nonzero Fourier mass

\[
\left(\sum_a p_a^2K_{aa}\right)^2
+4\sum_{a<b}p_a^2p_b^2.
\]

Therefore every fixed finite-type pair rule has \(Q(A_N)=\Omega(N^2)\).
A successful infinite tower must keep producing genuinely new,
asymptotically orthogonal types; bounded-period or bounded-type
regularity cannot provide it.
