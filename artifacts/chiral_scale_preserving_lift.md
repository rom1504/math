# Chiral scale-preserving compressed lifts

## Status

Let \(A\) be a centered/chiral signing:

\[
S^2=-I,\qquad S^\top AS=-A.
\]

For an exact compressed four-lift \(B\), the desired equality is

\[
Q(B)=8Q(A).
\tag{0.1}
\]

The current audit gives three rigorous conclusions:

1. chiral symmetry alone cannot guarantee an equality lift;
2. equality has strong exact centering, local-field, and cut
   conditions;
3. the natural inherited-chiral Clifford block family is completely
   ruled out for the order-\(12\), \(Q=40\) witness.

For that specific witness, existence in the unrestricted compressed
class remains unresolved. A finite \(1128\)-variable feasibility
system has been constructed. Its LP relaxation remains feasible even
after all one-fibre cut inequalities are imposed; an integer solve did
not finish at this checkpoint. Thus there is not yet an all-class
UNSAT certificate.

## 1. Equality in the range theorem

Let

\[
U(A)=\max_xx^\top Ax,\qquad
L(A)=\min_xx^\top Ax.
\]

For a compressed \(s\)-lift, every constant-fibre spin has energy

\[
E_B(x\otimes\mathbf1_s)=s^{3/2}x^\top Ax+d_B,
\tag{1.1}
\]

where

\[
d_B=\sum_i\mathbf1^\top D_i\mathbf1
\]

is the total diagonal-fibre shift.

If \(A\) is chiral, then

\[
U(A)=-L(A)=Q(A).
\]

Therefore the two repeated seed extrema have energies

\[
s^{3/2}Q(A)+d_B,\qquad
-s^{3/2}Q(A)+d_B.
\]

It follows immediately that

\[
Q(B)\ge s^{3/2}Q(A)+|d_B|.
\tag{1.2}
\]

Hence equality forces

\[
\boxed{d_B=0.}
\tag{1.3}
\]

It also forces every repeated positive seed maximizer to be a global
maximizer of \(B\), and every repeated negative seed maximizer to be a
global minimizer.

These are only the first equality conditions. The nonconstant fibre
states must all remain inside the same interval:

\[
\boxed{
|y^\top By|\le s^{3/2}Q(A)
\quad\text{for every }y\in\{\pm1\}^{ns}.
}
\tag{1.4}
\]

## 2. Exact local-field and cut conditions

Let \(X=x\otimes\mathbf1_s\), where \(x^\top Ax=Q(A)\). Since \(X\)
must globally maximize the lifted energy, switching by \(X\) makes
every cut nonnegative:

\[
\boxed{
\sum_{\substack{u\in T\\v\notin T}}
b_{uv}X_uX_v\ge0
\quad\text{for every }T\subset[ns].
}
\tag{2.1}
\]

For a repeated negative seed minimizer, every corresponding switched
cut must be nonpositive:

\[
\boxed{
\sum_{\substack{u\in T\\v\notin T}}
b_{uv}X_uX_v\le0.
}
\tag{2.2}
\]

Singleton cuts give the local-field conditions

\[
X_u(BX)_u\ge1
\]

at a repeated maximum, and

\[
X_u(BX)_u\le-1
\]

at a repeated minimum. The \(\pm1\) rather than weak zero bounds use
the fact that every microvertex has odd degree \(ns-1\).

For a fixed macro vertex \(i\) and microcoordinate \(\alpha\), define
the row sum into macro block \(j\) by

\[
r_{i\alpha,j}
=\sum_{\beta=1}^s
b_{(i,\alpha),(j,\beta)}
\]

and the internal row sum by

\[
d_{i\alpha}
=\sum_{\beta\ne\alpha}
b_{(i,\alpha),(i,\beta)}.
\]

Then a repeated macro spin has exact local field

\[
\boxed{
X_{i\alpha}(BX)_{i\alpha}
=x_i\sum_{j\ne i}r_{i\alpha,j}x_j+d_{i\alpha}.
}
\tag{2.3}
\]

For \(s=4\),

\[
r_{i\alpha,j}\in\{-4,-2,0,2,4\},\qquad
d_{i\alpha}\in\{-3,-1,1,3\}.
\]

Thus the repeated-extremizer constraints form a small-integer linear
system on block row sums before the individual block entries are even
considered.

## 3. Chiral symmetry is not sufficient

Take the order-two signing

\[
A_2=
\begin{pmatrix}
0&1\\
1&0
\end{pmatrix},
\qquad
S_2=
\begin{pmatrix}
0&-1\\
1&0
\end{pmatrix}.
\]

Then

\[
S_2^2=-I,\qquad S_2A_2=-A_2S_2,\qquad Q(A_2)=2.
\]

An equality four-lift would be an order-eight signing with

\[
Q=8Q(A_2)=16.
\]

But the certified exact value is

\[
\min_{\text{order }8}Q=2F(8)=20.
\]

Therefore:

\[
\boxed{
\text{There is no equality four-lift theorem valid for every chiral
signing.}
}
\tag{3.1}
\]

Any positive construction must use the full action profile of the
particular chiral seed, not merely the existence of a signed complex
structure.

## 4. The inherited-chiral Clifford ansatz

For the order-\(12\) witness, write \(J\) for its signed complex
structure and set

\[
C=JA.
\]

Then \(C\) is symmetric. The natural complementary-support ansatz is

\[
\widehat A
=A\otimes R
+C\otimes K
+J\otimes E
+\operatorname{diag}(D_1,\ldots,D_{12}),
\tag{4.1}
\]

where:

- \(R,K\) are symmetric partial-sign \(4\times4\) matrices with
  complementary support;
- \(R\) contains the microdiagonal and has total sum \(8\);
- \(K\) has total sum zero;
- \(E\) is skew-symmetric and supported on the off-diagonal support
  of \(K\);
- the \(D_i\) fill the remaining diagonal-fibre entries compatibly
  with inherited chirality.

These support conditions make \(\widehat A\) a valid symmetric
zero-diagonal signing and an exact compressed four-lift. A suitable
micro involution \(P\) makes \(J\otimes P\) an inherited signed complex
structure.

The complete exact enumeration in
`verify_chiral_clifford_no_go.py` gives:

- \(1008\) compatible \((P,R,K,E)\) masks;
- \(912\) already have a uniform-fibre profile value above \(320\);
- the remaining \(96\) form four microcoordinate-permutation orbits;
- for each orbit, one two-pattern Boolean fibre witness has energy at
  least \(392\) or \(480\), for every one of the \(2^{24}\) compatible
  diagonal fillings.

The last statement is certified by exact subset-sum dynamic
programming, not heuristic search. Consequently,

\[
\boxed{
\text{No lift in the complementary-support Clifford family attains
the equality target }320.
}
\tag{4.2}
\]

## 5. Unrestricted compressed-lift feasibility system

For the specific order-\(12\) witness and \(s=4\), an unrestricted
lift has

\[
\binom{48}{2}=1128
\]

binary edge variables. The script
`search_chiral_equality_constraints.py` imposes the following exact
necessary conditions:

1. all \(66\) cross-block compression equalities;
2. the zero diagonal-fibre shift (1.3);
3. all \(48\) singleton local-field inequalities for each of the nine
   positive and nine negative antipodal extremizer classes;
4. every cut inequality (2.1)--(2.2) whose flipped set lies entirely
   in one macro fibre.

For the first witness this gives:

\[
1128\ \text{binary variables},\qquad
3307\ \text{linear constraints}.
\]

For the second witness it gives:

\[
1128\ \text{binary variables},\qquad
2227\ \text{linear constraints}.
\]

The LP relaxations of both systems are feasible. An integer
feasibility run on the first system did not settle within the
checkpoint window. Hence these necessary conditions do not presently
yield either:

- a fractional/Farkas impossibility certificate; or
- a Boolean candidate equality lift.

The exact full problem is the same binary system augmented by the cut
constraints for every subset of all \(48\) microvertices. A practical
exact route is cutting-plane SAT/MILP:

1. solve the current master system;
2. maximize/minimize the candidate quadratic form over the Boolean
   cube;
3. if an energy exceeds \(320\), add that spin inequality;
4. repeat until a \(Q\le320\) witness or an UNSAT certificate appears.

The separation step is itself a \(48\)-vertex signed MaxCut problem,
so a certified run requires a serious SAT/branch-and-bound backend.

## 6. Hamming-slice second-moment conditions

There is a second exact family of necessary constraints. Fix arbitrary
grid magnetizations in every macro fibre and sample each fibre
uniformly on its Hamming slice. Equality implies

\[
\boxed{
\mathbb E\left[(X^\top BX)^2\right]\le(8Q(A))^2
}
\tag{6.1}
\]

for every such product distribution.

The left side has an exact Hoeffding decomposition into:

- the squared conditional mean;
- nonnegative one-fibre field variances;
- nonnegative two-fibre interaction variances;
- the corresponding within-fibre Hoeffding components.

Unlike the linear cut system, (6.1) is quadratic in the block entries.
It can be linearized with pair-product variables and semidefinite/SAT
consistency constraints. This is the next principled strengthening if
the pure cutting-plane system remains feasible.

## 7. Verdict

The equality route has narrowed substantially:

\[
\boxed{
\begin{array}{c}
\text{chiral symmetry alone: insufficient;}\\
\text{common tensor continuation: fails;}\\
\text{natural Clifford continuation: exactly impossible;}\\
\text{all dependent compressed lifts of the }Q=40\text{ seed:
unresolved.}
\end{array}
}
\]

Any surviving equality lift must be nonlocal across macro edges and
lie outside complementary-support Clifford formulas. The unrestricted
question is now a finite, reproducible \(1128\)-binary feasibility
problem with an exact separation oracle.

## 8. Bounded nonlocal cutting-plane pass

A final bounded pass was run after the analytic audit.

### Available solver

No CP-SAT, Z3, PySAT, or standalone SAT solver is installed in the
workspace. The pass therefore used the HiGHS mixed-integer solver
through SciPy and was deliberately time-bounded.

### Symmetry-restricted lazy master

The output was required to inherit the seed complex structure as

\[
\widehat S=S\otimes I_4,\qquad
\widehat S^\top B\widehat S=-B.
\]

The initial master contained:

- all compression equalities;
- exact zero centering;
- every repeated-extremizer local-field inequality;
- all inherited-chiral edge-orbit equalities.

This is a \(1128\)-binary system with \(1507\) linear constraints.
HiGHS found an integer solution. Exact separation then found a
two-pattern/local-search spin with

\[
\boxed{x^\top Bx=440>320.}
\]

The rejected matrix and its exact separating spin are stored in
`chiral_equality_rejected_candidate.npz`. The verifier
`verify_chiral_equality_bounded_checkpoint.py` confirms:

1. the matrix is a valid order-\(48\) signing;
2. all \(1507\) initial constraints hold;
3. inherited chirality holds exactly;
4. the separating energy is exactly \(440\).

After adding that one spin inequality, the integer master had \(1508\)
constraints. HiGHS returned no primal solution or UNSAT proof within
the \(45\)-second bounded window. The precise sparse instance is stored
in `chiral_equality_lazy_instance.npz`.

### Unrestricted packaged master

The full nonlocal necessary system, without inherited-chiral
restriction, has also been packaged:

\[
\boxed{
1128\text{ binary variables},\quad
3307\text{ constraints},\quad
315624\text{ nonzero coefficients}.
}
\]

It contains:

- all \(66\) compression equations;
- exact zero centering;
- all repeated positive/negative local fields;
- every extremal cut supported inside one macro fibre.

The file `chiral_equality_unrestricted_master.npz` includes the sparse
matrix, bounds, edge-variable ordering, labels, seed, all seed
extremizers, and target. It is produced reproducibly by
`package_chiral_equality_instances.py`.

### Bounded-pass conclusion

The pass produced neither a target lift nor a verifiable unrestricted
UNSAT certificate. It did produce:

- one exact rejected candidate and separator;
- one symmetry-restricted lazy master;
- one complete unrestricted necessary-condition master.

This route is stopped here rather than continuing generic MILP search.
