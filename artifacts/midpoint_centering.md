# Midpoint centering of optimal signings

## 1. Definitions

For a signing \(A\), write

\[
P(A)=\max_xH_A(x),\qquad
Q(A)=-\min_xH_A(x),
\]

\[
M(A)=\max(P,Q),\qquad
W(A)=\frac{P+Q}{2},\qquad
d(A)=\frac{P-Q}{2}.
\]

Then

\[
M(A)=W(A)+|d(A)|.
\tag{1}
\]

The hoped-for theorem is that an exact, or at least asymptotic, minimizer of
\(M\) can be chosen with \(|d(A)|=O(1)\) or \(o(n^{3/2})\).

## 2. Exact coding translation

Let \(E=\binom n2\), encode the negative edges by
\(a\in\mathbb F_2^E\), and let \(C\) be the cut code of \(K_n\). Then

\[
P(A)=E-2d(a,C),
\]

\[
Q(A)=E-2d(a,\mathbf1+C).
\tag{2}
\]

The augmented code

\[
D=C\cup(\mathbf1+C)
\]

is antipodal, and minimizing \(M(A)\) is equivalent to choosing a deepest
hole of \(D\). Midpoint centering is exactly the following coding lemma:

> There is a deepest hole \(a\) of \(D\) whose distances to the two
> \(C\)-cosets differ by at most one.

Such a hole would satisfy

\[
|P(A)-Q(A)|\le2.
\tag{3}
\]

This lemma survived exhaustive tests for all antipodal subsets through
ambient length \(5\), thousands of random linear index-two examples through
length \(9\), and every cut-code minimizer checked below. A proof was not
found. The missing step is an antipodal plateau/connectivity statement for
the deepest-hole set; it is not true for an arbitrary antipodal vertex set
without the distance-to-code structure.

## 3. Exact small-order evidence

Exhaustive switching-class enumeration gives the following profiles among
optimal signings:

| \(n\) | \(M_n\) | optimal \((P,Q)\) profiles |
|---:|---:|---|
| 2 | 1 | \((1,1)\) |
| 3 | 3 | \((3,1),(1,3)\) |
| 4 | 4 | \((4,4)\) |
| 5 | 4 | \((4,4)\) |
| 6 | 5 | \((5,5)\) |
| 7 | 9 | \((9,7),(7,9)\) |

Thus every exact minimizer through \(n=7\) has

\[
|P-Q|\le2.
\tag{4}
\]

Orders \(3\) and \(7\) show that exact equality \(P=Q\) need not be possible;
the energy lattice makes (3) the natural sharp statement.

## 4. What edge-flip optimality proves

Suppose \(P\ge Q+4\) and a single edge flip cannot reduce \(M=P\). For every
edge \(e\), flipping \(e\) must leave a positive energy at least \(P\).
Therefore there is a configuration \(x_e\) satisfying

\[
H_A(x_e)\ge P-2,\qquad
a_e x_{e_1}x_{e_2}=-1.
\tag{5}
\]

So every coefficient disagreement occurs in the \(P-2\) near-ground layer.
This is a strong local certificate, but the universal Hamming noise cloud
shows that a large near-ground layer is not contradictory. Edge-local
optimality alone does not force \(P-Q\le2\).

In the coding picture, the obstruction is a potentially isolated local
maximum of the distance to one \(C\)-coset inside the deepest-hole set of
\(D\).

## 5. Third moments and endpoint bounds

For a uniform spin \(X\),

\[
\mathbb EH_A(X)=0,\qquad
\mathbb EH_A(X)^2=E,
\]

\[
\mathbb EH_A(X)^3
=6\sum_{i<j<k}a_{ij}a_{jk}a_{ki}.
\tag{6}
\]

The endpoint inequalities

\[
-Q\,E\le\mathbb EH_A^3\le P\,E
\]

and the Bhatia--Davis bound

\[
PQ\ge E
\tag{7}
\]

do not control \(P-Q\) on the \(n^{3/2}\) scale. Orientation-odd triangle
information can also vanish identically on self-complementary signings.

## 6. What centering would and would not imply

If (3), or its \(o(n^{3/2})\) version, holds for optimal signings, then

\[
M_n=W(A_n)+o(n^{3/2}),
\]

and the problem becomes asymptotically a centered signed-cut-norm problem.
However, the exact partition inequality

\[
\sum_rW(A[S_r])\le W(A)
\]

scales as \(\alpha\), not \(\alpha^{3/2}\), so centering alone does not prove
convergence.

Centering does remove the orientation loss in the max-extension coupling.
For a partition

\[
A=\begin{pmatrix}B&C\\C^\top&D\end{pmatrix},
\]

let \(u\) maximize \(B\), \(v\) minimize \(B\), and set

\[
R_C(u,v)=
\max\{\|C^\top(u-v)\|_1,\|C^\top(u+v)\|_1\}.
\]

The exact coupled extension inequality gives

\[
P(B)+Q(B)\le2M(A)-R_C(u,v).
\tag{8}
\]

If \(B\) is centered up to \(o(n^{3/2})\), then

\[
M(B)\le M(A)-\frac12R_C(u,v)+o(n^{3/2}).
\tag{9}
\]

Consequently, proportional scale transfer would follow from the quantitative
rectangular-separation statement

\[
R_C(u,v)
\ge2(1-\alpha^{3/2})M(A)-o(n^{3/2}).
\tag{10}
\]

Thus midpoint balance is useful but not terminal: it converts the remaining
problem into the explicit cross-separation inequality (10).

## 7. Verdict

All exact evidence supports the sharp bound \(|P-Q|\le2\) for a suitably
chosen minimizer, but the deepest-hole balancing lemma remains unproved.
Even if granted, the existing hereditary \(W\)-inequality loses
\(\sqrt\alpha\); the next required theorem is the rectangular separation
(10), not centering alone.

