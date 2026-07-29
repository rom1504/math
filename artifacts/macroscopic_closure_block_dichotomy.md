# Macroscopic closure blocks after one-sided discrepancy balance

Checkpoint date: 2026-07-26.

## 1. Verdict

The adaptive closure identity

\[
g_{j+1}\le \frac12g_j
\]

already removes the cumulative replenishment obstruction.  The
remaining macroscopic-block problem cannot be closed from range
superadditivity alone: the bent/Walsh cross-block construction in the
master note suppresses the required ground-layer interaction excess.

There is, however, a new strict principal-decomposition theorem once
the one-sided graph-discrepancy product is added.  In the one-copy
normalization

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
q(A)=\max_x|H_A(x)|,
\]

every fixed-ratio principal split of a competitive signing satisfies

\[
\boxed{
q(A)\ge q(B)+
\frac{|D|^3}{12800\,q(D)}(1-o(1)),
\qquad
q(A)\ge q(D)+
\frac{|B|^3}{12800\,q(B)}(1-o(1)).
}
\tag{1.1}
\]

Consequently, for the minima \(q_n=M_n\), uniformly when
\(\delta n\le k\le(1-\delta)n\),

\[
\boxed{
q_{n-k}
\frac{k^3}{12800\,q_n}(1-o_\delta(1))
\le q_n
}
\tag{1.2}
\]

and symmetrically with \(k,n-k\) interchanged.  Equivalently,

\[
\boxed{
q_n^2
\ge
q_{n-k}^2+\frac{k^3}{12800}(1-o_\delta(1)).
}
\tag{1.3}
\]

This is a genuine scale-transfer inequality, but its gain is cubic in
the deleted proportion.  The \(3/2\)-homogeneous decrement needed for
convergence is linear at a small deleted proportion.  A slowly
oscillating scalar sequence satisfies (1.2), monotonicity, local
continuity, the rigorous cage, and one-sided block superadditivity.
Thus (1.1)--(1.3), even together with the closure half-contraction, do
not force convergence.

The second conclusion concerns the feature-kernel route.  At exact
cap, the asserted conditional block degeneracy is algebraically
tautological:

\[
h_T\cdot\phi_T=d-h_{T^c}\cdot\phi_{T^c}.
\tag{1.4}
\]

It therefore localizes no vertex block by itself.  A useful inverse
theorem must produce a **sparse or boundary-avoiding predictor**, not
merely a zero Schur complement.

## 2. One-sided notation

For a signing \(A\), put

\[
p(A)=\max_xH_A(x),\qquad
n(A)=-\min_xH_A(x),
\]

\[
q(A)=\max\{p(A),n(A)\},\qquad
\ell(A)=\min\{p(A),n(A)\}.
\tag{2.1}
\]

Thus

\[
w(A)=\frac{p(A)+n(A)}2,\qquad
d(A)=\frac{p(A)-n(A)}2,
\]

and

\[
q(A)=w(A)+|d(A)|,\qquad
\ell(A)=w(A)-|d(A)|.
\tag{2.2}
\]

For a vertex partition, write

\[
A=
\begin{pmatrix}
B&C\\
C^\mathsf T&D
\end{pmatrix}.
\tag{2.3}
\]

The relative global sign of the two block spin vectors changes the
cross term and leaves both internal energies fixed.  Hence

\[
p(A)\ge p(B)+p(D),
\qquad
n(A)\ge n(B)+n(D).
\tag{2.4}
\]

If the dominant side of \(B\) is positive, the first inequality gives
\(q(A)\ge q(B)+p(D)\); if it is negative, the second gives
\(q(A)\ge q(B)+n(D)\).  Therefore, exactly,

\[
\boxed{
q(A)\ge q(B)+\ell(D),
\qquad
q(A)\ge q(D)+\ell(B).
}
\tag{2.5}
\]

This is the orientation-balanced form of block superadditivity.

## 3. Bollobás--Scott input in the present normalization

The one-sided graph-discrepancy product, applied after switching a
positive ground state to \(\mathbf1\), gives the following asymptotic
statement for every sequence with \(q(A_m)=O(m^{3/2})\):

\[
\boxed{
\ell(A_m)\ge
\frac{m^3}{12800\,q(A_m)}(1-o(1)).
}
\tag{3.1}
\]

Here is the normalization check.  In doubled notation let

\[
P=2p,\qquad N=2n,\qquad Q=2q,\qquad R=P+N.
\]

The graph discrepancy theorem and the switching reduction give

\[
PR\ge\frac{m^3}{1600}(1-o(1)),
\qquad
NR\ge\frac{m^3}{1600}(1-o(1)).
\tag{3.2}
\]

Since \(R\le2Q\),

\[
\min(P,N)\ge
\frac{m^3}{3200Q}(1-o(1)).
\]

Dividing by two and using \(Q=2q\) yields (3.1).

For a fixed-ratio principal block of a signing with
\(q(A_n)=O(n^{3/2})\), principal monotonicity gives
\[
q(A[S])\le q(A_n)=O_\delta(|S|^{3/2}),
\]
so (3.1) applies uniformly to every block of order at least
\(\delta n\).

## 4. Strict principal decomposition

Substitute (3.1) into (2.5).  If both blocks have macroscopic order,

\[
q(A)
\ge
q(B)+
\frac{|D|^3}{12800q(D)}(1-o(1)),
\]

and the symmetric inequality follows identically.  This proves
(1.1).

There are two useful corollaries.

First, because \(q(B),q(D)\le q(A)\),

\[
\boxed{
q(B)\le q(A)-
\frac{|D|^3}{12800q(A)}(1-o(1)),
}
\tag{4.1}
\]

and likewise for \(D\).  Thus a macroscopic closure block really does
force a strict drop in the absolute norm of the surviving core,
independently of the sign of the block's dominant extremum.

Second, multiplying the two inequalities in (1.1) and using
\[
q(B)\ell(B)\ge\frac{|B|^3}{12800}(1-o(1)),
\qquad
q(D)\ell(D)\ge\frac{|D|^3}{12800}(1-o(1)),
\]
gives

\[
\boxed{
q(A)^2
\ge
q(B)q(D)
+
\frac{|B|^3+|D|^3}{12800}(1-o(1)).
}
\tag{4.2}
\]

Indeed,
\[
q(A)^2
\ge
(q(B)+\ell(D))(q(D)+\ell(B)),
\]
and the two diagonal products already give the displayed cubic
terms.

Now take \(A\) to be a global order-\(n\) minimizer and choose any
partition with \(|D|=k\).  Since
\[
q_{n-k}\le q(B),\qquad q(D)\le q(A)=q_n,
\]
(1.1) gives (1.2).  Multiplying its decrement by \(q_n\) gives
(1.3), because
\[
q_n^2-q_{n-k}^2
\ge q_n(q_n-q_{n-k}).
\]

## 5. Exact orientation-aware replacement

The same notation produces a refinement of the induced-block
replacement theorem.  Let \(A\) be a global minimizer of \(q\), and
let \(\widetilde B\) be any signing on the \(B\)-block.  Replacing
\(\widetilde B\) by its edgewise negative swaps
\(\widetilde p,\widetilde n\), hence changes
\(\widetilde d\) to \(-\widetilde d\) and leaves \(\widetilde w\)
fixed.  If
\[
R(C)=\|C\|_{\infty\to1},
\]
then
\[
q\!\begin{pmatrix}
\widetilde B&C\\ C^\mathsf T&D
\end{pmatrix}
\le
\widetilde w+w(D)
+|\widetilde d+d(D)|+R(C).
\]

Choose the better of \(\widetilde B,-\widetilde B\), and compare with
the internal-block lower bound for the original \(A\).  Global
minimality yields

\[
\boxed{
w(B)+|d(B)+d(D)|
\le
\widetilde w+
\bigl||\widetilde d|-|d(D)|\bigr|
+R(C).
}
\tag{5.1}
\]

Taking \(\widetilde B\) to be an order-\(|B|\) absolute-norm
minimizer gives the simpler consequence

\[
\boxed{
q(B)\le q_{|B|}+2|d(D)|+R(C)
=q_{|B|}+q(D)-\ell(D)+R(C).
}
\tag{5.2}
\]

Combining (5.2) with (3.1) gives the product-balanced exchange bound

\[
\boxed{
q(B)
\le q_{|B|}+q(D)+R(C)
-
\frac{|D|^3}{12800q(D)}(1-o(1)).
}
\tag{5.3}
\]

The new term is a strict macroscopic discount.  It is the strongest
consequence of one-sided balance available without controlling the
boundary norm or the complement midpoint.

## 6. Why the cubic decrement does not prove convergence

Put
\[
a_n=\frac{q_n}{n^{3/2}}.
\]
For \(k=\alpha n\), (1.2) reads
\[
a_n
\ge
(1-\alpha)^{3/2}a_{(1-\alpha)n}
+
\frac{\alpha^3}{12800a_n}
+o_\alpha(1).
\tag{6.1}
\]

The homogeneity deficit that must be paid to compare the two
normalized values is
\[
1-(1-\alpha)^{3/2}
=\frac32\alpha+O(\alpha^2),
\]
whereas the new payment in (6.1) is \(O(\alpha^3)\).  The mismatch is
quadratic in the small deleted proportion.

This is not bookkeeping.  Fix, for example,
\[
a(t)=0.42+0.02\sin(\log\log(t+e^e)),
\qquad
f(t)=t^{3/2}a(t).
\tag{6.2}
\]
After changing finitely many initial values, \(f(n)\) is
nondecreasing and
\[
f(n+1)-f(n)=O(\sqrt n).
\]
For every fixed \(0<\alpha<1\),
\[
\frac{f(n)-f((1-\alpha)n)}{n^{3/2}}
=a(n)\bigl[1-(1-\alpha)^{3/2}\bigr]+o(1).
\]
Since
\[
\inf_n a(n)^2
>
\frac1{12800}
\sup_{0\le\alpha\le1}
\frac{\alpha^3}{1-(1-\alpha)^{3/2}},
\]
the sequence satisfies the asymptotic form of (6.1), uniformly on
every fixed-ratio compact interval.  Its normalized values nevertheless
oscillate forever.

The same scalar model can be made compatible with one-sided
superadditivity by setting
\[
p(n)=n(n)=f(n).
\]
The margin between \(n^{3/2}\) and
\(k^{3/2}+(n-k)^{3/2}\), together with the small amplitude in (6.2),
makes both one-sided block inequalities hold asymptotically; finite
values can again be adjusted.  It also satisfies the universal range
lower bound and the discrepancy-product lower bound by a wide
constant margin.

Therefore the present scalar information permits slow oscillation.
To obtain convergence, (6.1) must be upgraded to a payment comparable
to
\[
q_n\left[1-(1-\alpha)^{3/2}\right],
\tag{6.3}
\]
or an equivalent profile/chaining theorem must control the cross
boundary.

## 7. Conditional feature kernels do not localize by themselves

At an exact cross cap, the internal feature vector
\(\phi(z)=(z_iz_j)_{ij\in E_h}\) obeys
\[
h\cdot\phi(z)=d
\qquad\text{on the whole cap family}.
\tag{7.1}
\]
For every coordinate block \(T\subseteq E_h\), this is exactly
\[
h_T\cdot\phi_T(z)
=d-h_{T^c}\cdot\phi_{T^c}(z).
\tag{7.2}
\]
Thus \(h_T\cdot\phi_T\) is linearly predicted by the outside features
for **every** \(T\).  The zero conditional Schur-complement identity
is simply (7.2) in covariance language.  It does not become stronger
when \(T\) is the edge set of an induced vertex block.

At positive cap slack the same issue persists: the global tube
\[
|h\cdot\phi-d|\le O(rn^{3/2})
\]
immediately gives the same error for (7.2).  No spatial localization
has occurred.

A viable inverse theorem therefore needs an additional conclusion,
for example:

* the predictor in (7.2) can be replaced by one supported only on a
  small vertex boundary;
* the outside coefficient mass can be sparsified while preserving the
  cap tube;
* or one finds an induced block \(T\) for which the predictor avoids
  the edges altered by the global replacement.

Without such a sparse-predictor statement, (10.136)--(10.137) cannot
be converted into (10.128)--(10.130): the outside term in (7.2)
records precisely the obstruction that a block replacement leaves
untouched.

## 8. Surviving target

The sharp next theorem suggested by this audit is a tangent-strength
principal decrement:

\[
q(B)
\le
q(A)
-
c\,|D|\sqrt{|A|}
+o(|A|^{3/2})
\tag{8.1}
\]

for the disagreement block \(D\) selected by adaptive closure, at
least after summing over a closure chain.  The product theorem proves
only
\[
q(B)
\le
q(A)
-
c\,\frac{|D|^3}{q(A)}.
\tag{8.2}
\]

For \(|D|=\alpha|A|\), (8.1) is linear in \(\alpha\), while (8.2) is
cubic.  Any successful use of the flat feature kernel must recover
the missing factor \(\alpha^{-2}\) through cross-boundary traffic or a
sparse conditional predictor.
