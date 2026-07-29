# Exact low-traffic child grounds: complementarity, uncrossing, and inverse reductions

Checkpoint date: 2026-07-26.

## 1. Scope and status

Let \(A=(a_{ij})\) be a symmetric zero-diagonal signing and use the
one-copy energy

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
M(A)=\max_x|H_A(x)|.
\]

Fix a full spin \(z\), switch by \(z\), and call an \(m\)-set \(S\)
**positive good** when

\[
\mathbf 1_S\text{ is a positive absolute ground of }A[S],
\qquad
A_{S^c,S}\mathbf 1_S=0.                         \tag{1.1}
\]

An arbitrary joint ground/traffic family can be split into its two
orientations; this loses only a factor two in cardinality and no
exponential rate.

This note proves:

1. exact traffic plus the ground condition is a strict
   complementarity system;
2. two crossing good sets obey an exact signed uncrossing law;
3. adjacent good sets force a literal partial twin-row relation;
4. \(o(n)\) traffic differs from exact traffic on only \(o(n)\)
   outside coordinates;
5. incidence saturation gives an exact multivariate
   Littlewood--Offord problem for a macroscopic collection of flat
   row vectors;
6. the two canonical maximum-entropy pair systems, full pair-union
   slices and full transversals, are incompatible with the child
   ground condition at large order;
7. any constant-weight affine family at the incidence threshold is an
   \(O(\log n)\)-defect paired quotient, and exact boundary Fourier
   separation rules it out completely using the universal child lower
   bound.

The missing step remains an inverse theorem for a code-separated
family: a speed-\(n\) family can have no adjacent pairs, so the local
twin lemma by itself does not classify the family.

There is a useful exact optimization/network reformulation.  The
supports in question are strictly complementary fixed points of a
symmetric zero-threshold Boolean network, and equivalently binary
solutions of one symmetric linear complementarity problem.  This
reduction is recorded in Proposition 2.2 below.  It also exposes another
necessary warning: even strict complementarity, without the absolute
child-ground condition, admits \(2^{n-o(n)}\) middle-layer supports.

## 2. Strict complementarity

### Proposition 2.1

If \(S\) is positive good, then \(m\) is even and

\[
\boxed{
A\mathbf 1_S\ge \mathbf 1_S
\quad\text{coordinatewise},\qquad
\operatorname{supp}(A\mathbf 1_S)=S.
}                                                   \tag{2.1}
\]

Equivalently,

\[
(A\mathbf 1_S)_j=0\quad(j\notin S),\qquad
(A\mathbf 1_S)_i\in\{1,3,\ldots,m-1\}\quad(i\in S). \tag{2.2}
\]

#### Proof

For \(j\notin S\), the first assertion in (1.1) is not needed:
the second assertion gives

\[
(A\mathbf 1_S)_j=\sum_{i\in S}a_{ji}=0.
\]

This is a sum of \(m\) signs, so \(m\) is even.

For \(i\in S\), flipping only \(i\) in the positive ground
\(\mathbf 1_S\) cannot increase the induced energy.  Hence

\[
(A\mathbf 1_S)_i
=\sum_{j\in S\setminus\{i\}}a_{ij}\ge0.
\]

This is a sum of the odd number \(m-1\) of signs, so it is in fact
at least one. \(\square\)

The negative orientation gives the same statement with
\(A\mathbf 1_S\le-\mathbf 1_S\) on \(S\).

### Proposition 2.2 (fixed-point and LCP formulations)

Put \(\eta=\mathbf 1_S\).  The strict-complementarity part of being
positive good is exactly

\[
\boxed{
\eta_i=\mathbf 1\{(A\eta)_i>0\}
\quad(i\in[n]),
\qquad
(A\eta)_i\in\{0,1,3,\ldots,m-1\}.
}                                                   \tag{2.3}
\]

Equivalently, if

\[
x=\mathbf 1-\eta,\qquad
w=A\mathbf 1-Ax=A\eta,
\]

then

\[
\boxed{
x\ge0,\quad w=A\mathbf 1-Ax\ge0,\quad x^\mathsf Tw=0,
\quad x\in\{0,1\}^n,
}                                                   \tag{2.4}
\]

with strict complementarity \(w_i\ge1\) whenever \(x_i=0\).
Thus the family is a fixed-weight binary subfamily of the solution
set of the symmetric LCP

\[
\operatorname{LCP}(q=A\mathbf 1,M=-A).
\]

For every \(0<\varepsilon<1\), it is also a family of strict local
maxima of

\[
f_\varepsilon(\eta)
=\sum_{i<j}a_{ij}\eta_i\eta_j-\varepsilon|\eta|.
\tag{2.5}
\]

Indeed, deleting an active coordinate changes \(f_\varepsilon\) by
at most \(-1+\varepsilon<0\), while adding an inactive coordinate
changes it by exactly \(-\varepsilon<0\).

Hunter Spink's sharp theorem on local maxima of quadratic Boolean
functions therefore gives the universal bound

\[
|\mathcal F|\le \binom n{\lfloor n/2\rfloor}.
\tag{2.6}
\]

That theorem alone is exponentially too weak for the present
\(2^{n/2+o(n)}\) target.  More decisively, Proposition 5.3 below
shows that both (2.3) and the bound (2.6) can be saturated at
exponential scale while the absolute child-ground condition fails.
Any useful refinement of the local-maxima/LCP route must retain the
global cut inequalities inside every selected support.

### Proposition 2.3 (exact complement decoupling)

If \(S\) is positive good and \(D=A[S^c]\), then for every
\(y\in\{\pm1\}^{S^c}\),

\[
\boxed{
H_A(\mathbf 1_S,y)
=M(A[S])+H_D(y).
}                                                   \tag{2.7}
\]

Consequently,

\[
\boxed{
M(A)\ge
\max\{M(A[S])+p(D),\ \nu(D)-M(A[S])\},
}                                                   \tag{2.8}
\]

where \(p(D)=\max H_D\) and \(\nu(D)=-\min H_D\).

Also, if \(\eta=\mathbf 1_S\) and
\(u=2\eta-\mathbf 1\), then

\[
\boxed{H_A(u)=H_A(\mathbf 1).}                      \tag{2.9}
\]

Thus an exact-good family gives a balanced equal-energy level set of
the full quadratic form, together with the much stronger rowwise
complementarity (2.1).

#### Proof

The cross term in the first displayed energy is

\[
(A_{S^c,S}\mathbf 1_S)\cdot y=0.
\]

The internal \(S\)-energy is \(M(A[S])\), proving (2.4).  Substitute
positive and negative grounds of \(D\) to get (2.5).

Flipping every spin in \(S\) changes the full all-one energy by twice
the total signed cut \(e(S,S^c)\).  The rowwise boundary equation
implies that cut is zero, proving (2.6). \(\square\)

### Corollary 2.4 (pairwise positivity)

If \(S,T\) are positive good, then

\[
\boxed{
\mathbf 1_T^\mathsf TA\mathbf 1_S
\ge |S\cap T|.
}                                                   \tag{2.10}
\]

Indeed, \(A\mathbf 1_S\) is nonnegative, supported on \(S\), and at
least one on that support.

Thus, if \(X\) is the matrix whose columns are the indicators of a
positive good family, then

\[
X^\mathsf T(A-I)X
\]

is entrywise nonnegative.  This is stronger than boundary cancellation
alone, although it has not yet yielded a sufficiently rigid global
rank bound.

## 3. Exact conditional-covariance kernel

Let \(\mathcal F\) be a family of positive good \(m\)-sets and let
\(Y=\mathbf 1_S\), where \(S\) is uniform on \(\mathcal F\).  For a
vertex \(j\), put

\[
q_j=\Pr(Y_j=0),\qquad
\mu_j=\mathbb E(Y\mid Y_j=0),\qquad
\Sigma_j=\operatorname{Cov}(Y\mid Y_j=0),
\]

and let \(a_j\) denote row \(j\) of \(A\).

### Proposition 3.1

For every \(j\) with \(q_j>0\),

\[
\boxed{
a_j\cdot\mu_j=0,\qquad
\Sigma_j a_j=0.
}                                                   \tag{3.1}
\]

#### Proof

Conditional on \(Y_j=0\), exact boundary traffic says
\(a_j\cdot Y=0\) almost surely.  Its conditional mean and variance
therefore vanish. \(\square\)

This is the clean covariance formulation of the inverse problem:
there are \(n\) compatible flat kernel vectors, but the kernel law
depends on which coordinate is conditioned to be absent.  In the
paired model these kernels are generated by pair differences.

There is a quantitative version.

### Proposition 3.2 (\(o(n)\)-traffic purification)

Suppose instead that every \(S\in\mathcal F\) satisfies

\[
L(S):=\|A_{S^c,S}\mathbf 1_S\|_1\le L_n.
\]

Then

\[
\boxed{
\sum_j q_j\,|a_j\cdot\mu_j|\le L_n,
\qquad
\sum_j q_j\,a_j^\mathsf T\Sigma_j a_j\le mL_n.
}                                                   \tag{3.2}
\]

If \(m\) is even, all nonzero outside fields have absolute value at
least two.  Hence every \(S\) has at most \(L_n/2\) exceptional
outside vertices; on all the remaining outside vertices the boundary
equation is exact.  If \(m\) is odd, every outside field is odd, and

\[
L(S)\ge n-m.                                       \tag{3.3}
\]

Consequently, for \(m/n\to\alpha\in(0,1)\), the regime
\(L_n=o(n)\) is possible only at even \(m\), and then differs from
the exact-zero system on only \(o(n)\) outside rows per member.

#### Proof

Averaging the definition of \(L(S)\) gives

\[
\sum_j q_j\,
\mathbb E\bigl(|a_j\cdot Y|\mid Y_j=0\bigr)
\le L_n.
\]

The first inequality follows by Jensen.  Since
\(|a_j\cdot Y|\le m\),

\[
\operatorname{Var}(a_j\cdot Y\mid Y_j=0)
\le
\mathbb E((a_j\cdot Y)^2\mid Y_j=0)
\le
m\,\mathbb E(|a_j\cdot Y|\mid Y_j=0),
\]

which proves the second inequality.  The parity assertions follow
because an outside field is a sum of \(m\) signs. \(\square\)

## 4. Signed uncrossing

For disjoint vertex sets \(U,V\), write

\[
e(U,V)=\sum_{i\in U,\ j\in V}a_{ij}.
\]

### Theorem 4.1 (exact crossing law)

Let \(S,T\) be positive good and partition the vertices as

\[
X=S\cap T,\qquad
Y=S\setminus T,\qquad
Z=T\setminus S,\qquad
W=(S\cup T)^c.
\]

Then

\[
\boxed{
e(X,Y)=e(X,Z)=-e(Y,Z)=:t\ge0.
}                                                   \tag{4.1}
\]

Moreover, row by row on \(W\),

\[
\boxed{
e(\{w\},Y)=e(\{w\},Z)=-e(\{w\},X)
\quad(w\in W).
}                                                   \tag{4.2}
\]

In particular,

\[
\boxed{M(A)\ge M(A[Y\cup Z])\ge t.}                 \tag{4.3}
\]

If \(t=0\), then the cuts \(X\mid Y\) in \(A[S]\) and
\(X\mid Z\) in \(A[T]\) have zero value.  Flipping either side
therefore creates another exact positive child ground.  Thus the
zero branch produces additional ground-state directions rather than
an unstructured cancellation.

#### Proof

Since \(Z\subset S^c\), summing the rowwise boundary equations of
\(S\) over \(Z\) gives

\[
e(X,Z)+e(Y,Z)=0.
\]

Similarly, \(Y\subset T^c\) gives

\[
e(X,Y)+e(Y,Z)=0.
\]

This proves the equalities in (4.1).  Since the all-one state is a
positive ground of \(A[S]\), its cut \(X\mid Y\) is nonnegative,
so \(t=e(X,Y)\ge0\).  The same follows from \(A[T]\).

For \(w\in W\), the two rowwise boundary equations are

\[
e(\{w\},X)+e(\{w\},Y)=0,\qquad
e(\{w\},X)+e(\{w\},Z)=0,
\]

which proves (4.2).

Finally, on \(Y\cup Z\), compare the spin which is all one with the
spin obtained by flipping all of \(Z\).  Their energies differ by
\(2e(Y,Z)=-2t\), so one has absolute value at least \(t\).
A principal energy is the expectation of the full energy over
uniform outside spins, proving (4.3).  If \(t=0\), flipping \(X\)
inside \(S\) changes the child energy by \(-2e(X,Y)=0\); the
statement for \(T\) is identical. \(\square\)

The theorem is a precise uncrossing dichotomy:

* \(t>0\) exposes a negative interaction between the two symmetric
  difference cells;
* \(t=0\) enlarges the child ground cloud.

Competitive scale alone only gives \(t=O(n^{3/2})\), so a further
entropy-to-structure argument is still needed.

### Corollary 4.2 (one-swap twin rigidity)

Suppose

\[
T=S\setminus\{i\}\cup\{j\}
\]

and both \(S,T\) are positive good.  Then

\[
\boxed{
a_{ij}=-1,\qquad
\sum_{v\in S\setminus\{i\}}a_{iv}
=
\sum_{v\in T\setminus\{j\}}a_{jv}
=1,
}                                                   \tag{4.4}
\]

and

\[
\boxed{
a_{ki}=a_{kj}
\quad\text{for every }k\notin S\cup\{j\}.
}                                                   \tag{4.5}
\]

#### Proof

In Theorem 4.1, \(Y=\{i\}\), \(Z=\{j\}\).  Thus

\[
t=-a_{ij}\in\{-1,1\}.
\]

Since \(t\ge0\), \(a_{ij}=-1\) and \(t=1\).  The first two identities
are \(e(X,Y)=e(X,Z)=1\).  For \(k\in W\), subtracting the two boundary
equations gives \(a_{ki}=a_{kj}\). \(\square\)

This is a literal repeated-row conclusion on the common outside
coordinates.  It cannot finish the inverse theorem by itself:
constant-weight codes of the required exponential rate can avoid all
one-swap pairs.

### Theorem 4.3 (disjoint uncrossing deficits add exactly)

Let

\[
(S_\ell,T_\ell),\qquad \ell=1,\ldots,L,
\]

be pairs of positive good sets, and suppose their symmetric
difference supports

\[
U_\ell=S_\ell\triangle T_\ell
\]

are pairwise disjoint.  Then

\[
\boxed{
M\!\left(A\left[\bigcup_{\ell=1}^L U_\ell\right]\right)
\ge
\sum_{\ell=1}^L t(S_\ell,T_\ell).
}                                                   \tag{4.6}
\]

#### Proof

Write

\[
Y_\ell=S_\ell\setminus T_\ell,\qquad
Z_\ell=T_\ell\setminus S_\ell,\qquad
t_\ell=t(S_\ell,T_\ell).
\]

The crossing law gives

\[
e(Y_\ell,Z_\ell)=-t_\ell.
\]

Put

\[
b_\ell=H_{A[Y_\ell]}(\mathbf1)
+H_{A[Z_\ell]}(\mathbf1).
\]

On \(U_\ell\), the configuration which is uniform on both cells has
energy \(b_\ell-t_\ell\), while the configuration which has opposite
signs on \(Y_\ell,Z_\ell\) has energy \(b_\ell+t_\ell\).  Choose one
of these two relative configurations for each \(\ell\).  Since

\[
\max_{\delta_\ell\in\{\pm1\}}
\left|\sum_\ell(b_\ell+\delta_\ell t_\ell)\right|
=\left|\sum_\ell b_\ell\right|+\sum_\ell t_\ell,
\tag{4.7}
\]

they can be chosen so that the absolute sum of their internal
energies is at least \(\sum_\ell t_\ell\).

Now multiply every spin in block \(U_\ell\) by an independent global
sign \(\varepsilon_\ell\).  This leaves its internal energy unchanged,
while every interaction between two distinct blocks has mean zero
over the \(\varepsilon_\ell\)'s.  The expected full energy is
therefore the chosen sum of internal energies.  Some realization has
absolute energy at least the absolute value of that expectation,
which proves (4.6). \(\square\)

Thus a weighted matching in the hypergraph

\[
\{S\triangle T:S,T\in\mathcal F,\ t(S,T)>0\},
\]

with edge weight \(t(S,T)\), is already a principal Boolean
quadratic witness of the same total weight.  No cancellation between
different packed deficits is lost.

### Proposition 4.4 (exact mean-deficit identity)

Let \(S,T\) be independent uniform members of \(\mathcal F\), and put

\[
p_i=\Pr(i\in S),\qquad
p_{ij}=\Pr(i,j\in S).
\]

Then

\[
\boxed{
\mathbb E\,t(S,T)
=
\sum_{i<j}a_{ij}p_{ij}
\bigl(p_i+p_j-2p_{ij}\bigr).
}                                                   \tag{4.8}
\]

#### Proof

For a fixed edge \(ij\), it contributes to \(e(S\cap T,S\setminus
T)\) in the two disjoint cases

\[
i\in S\cap T,\ j\in S\setminus T,
\quad\text{or}\quad
j\in S\cap T,\ i\in S\setminus T.
\]

Independence of \(S,T\) makes their probabilities

\[
p_{ij}(p_i-p_{ij})
\quad\text{and}\quad
p_{ij}(p_j-p_{ij}),
\]

respectively.  Summing the signed edge contributions proves
(4.8). \(\square\)

Equation (4.8) is an exact energy/entropy interface, but its right
side is a signed sum and currently has no useful lower bound from
the one- and two-coordinate entropy data alone.

Theorem 4.3 also identifies the precise packing obstruction.  A
rate-\(1/2\) constant-weight code may have all pairwise symmetric
differences of linear size, so the associated hypergraph can have
matching number one even when every pair has positive deficit.
Conversely, a positive density of zero-deficit pairs does not by
itself give additive quadruples in \(\mathbb F_2^n\): \(t(S,T)=0\)
is a scalar cut equality, not the group relation
\(S_1\triangle S_2=S_3\triangle S_4\) required by
Balog--Szemerédi--Gowers.  An additional implication from zero
deficit to XOR closure is still needed before BSG/Freiman theory can
feed Theorem 7.3.

## 5. Exact incidence-to-vector-Littlewood--Offord reduction

The zero-boundary family gives a multivariate concentration problem
without any selector-dependent spin.

### Proposition 5.1

Let \(\mathcal F\subseteq\binom{[n]}m\) satisfy the zero-boundary
part of (1.1), put \(k=n-m\), and fix \(0\le t\le k\).  Then some
\(t\)-set \(J\) satisfies

\[
\boxed{
\Pr_{\xi\in\{0,1\}^{J^c}}
\left(A_{J,J^c}\xi=0\right)
\ge
2^{-(n-t)}|\mathcal F|
\frac{\binom kt}{\binom nt}.
}                                                   \tag{5.1}
\]

If

\[
|\mathcal F|\ge2^{-m}\binom nm,
\]

then the right side is at least

\[
\boxed{
2^{-(n-t+m)}
\binom nm\frac{\binom kt}{\binom nt}.
}                                                   \tag{5.2}
\]

At \(m=k=n/2\) and \(t=o(n)\), this is

\[
\boxed{
2^{-n/2}
\exp\!\left(
-O\!\left(\frac{t^2}{n}+\log n\right)
\right).
}                                                   \tag{5.3}
\]

#### Proof

For \(J\in\binom{[n]}t\), let

\[
\mathcal F_J=\{S\in\mathcal F:J\subseteq S^c\}.
\]

Double counting gives

\[
\sum_{|J|=t}|\mathcal F_J|
=|\mathcal F|\binom kt.
\]

Choose \(J\) for which \(|\mathcal F_J|\) is at least its average.
Every \(S\in\mathcal F_J\) is one of the \(2^{n-t}\) Boolean vectors
on \(J^c\), and its boundary equations on the rows \(J\) say

\[
A_{J,J^c}\mathbf 1_S=0.
\]

This proves (5.1), and (5.2) follows from the incidence threshold.
For the balanced asymptotic, use

\[
\frac{\binom{n/2}{t}}{\binom nt}
=2^{-t}\exp\!\left(-O(t^2/n)\right)
\]

and the central binomial estimate. \(\square\)

For columns \(v_i\) of \(A_{J,J^c}\), (5.1) is the atom

\[
\Pr\left(\sum_{i\in J^c}\xi_i v_i=0\right).
\]

The elementary rank bound

\[
\sup_u\Pr\left(\sum_i\xi_i v_i=u\right)
\le2^{-\operatorname{rank}\{v_i\}}                 \tag{5.4}
\]

follows by conditioning on all but a maximal independent set of
columns.  At the incidence scale this rank bound is too weak.  What
is needed is a genuinely exponential inverse Littlewood--Offord
theorem: an atom as large as (5.2), simultaneously with the positive
child-ground inequalities, should force most flat column patterns
into repeated/opposite types or a low-complexity affine quotient.

Proposition 5.1 gives the exact constants and removes an ambiguity in
that target.

### Proposition 5.2 (the boundary-only \(2^{n/2}\) conjecture is false)

Let \(n\) be divisible by four, choose
\(s\in\{\pm1\}^n\) with \(n/2\) signs of each kind, and put

\[
A=ss^\mathsf T-I.
\]

Then the middle-layer boundary system

\[
(1-\eta_j)(A\eta)_j=0\quad(j\in[n]),
\qquad |\eta|=n/2,                                 \tag{5.5}
\]

has

\[
\boxed{
\binom{n/2}{n/4}^2
=2^{n-o(n)}
}                                                   \tag{5.6}
\]

solutions.

#### Proof

Take every \(S\) containing \(n/4\) coordinates of each \(s\)-type.
Then

\[
\sum_{i\in S}s_i=0.
\]

For \(j\notin S\),

\[
(A\mathbf 1_S)_j
=s_j\sum_{i\in S}s_i=0.
\]

The number of such sets is (5.6). \(\square\)

This example does not satisfy the ground condition.  On every such
\(S\), with \(m=n/2\),

\[
H_{A[S]}(\mathbf 1)
=\frac{(\sum_{i\in S}s_i)^2-m}{2}
=-\frac m2,
\]

whereas the spin \(x=s_S\) has

\[
H_{A[S]}(x)=\binom m2.
\]

Thus the proposed boundary-only estimate
\(2^{n/2+o(n)}\) is false by an exponential factor.  The corrected
extremal target must include strict complementarity and the
**absolute** child-ground condition.  This rank-one construction also
explains why a generic row-kernel or covariance theorem cannot be
enough: all row constraints collapse onto the single type-balance
equation \(\sum_{i\in S}s_i=0\).

### Proposition 5.3 (strict complementarity alone is still false)

Keep the same balanced type vector \(s\), but reverse the signing:

\[
A=I-ss^\mathsf T,
\qquad a_{ij}=-s_is_j\quad(i\ne j).
\tag{5.7}
\]

For every middle set \(S\) containing \(n/4\) vertices of each type,

\[
\boxed{
A\mathbf 1_S=\mathbf 1_S.
}                                                   \tag{5.8}
\]

Thus all

\[
\binom{n/2}{n/4}^2=2^{n-o(n)}
\tag{5.9}
\]

sets satisfy the full strict-complementarity/fixed-point system
(2.1), not merely the outside boundary equations.

#### Proof

For \(j\notin S\),

\[
(A\mathbf 1_S)_j
=-s_j\sum_{i\in S}s_i=0.
\]

For \(j\in S\), omission of the diagonal term gives

\[
(A\mathbf 1_S)_j
=-s_j\sum_{i\in S\setminus\{j\}}s_i
=-s_j(-s_j)=1.
\]

This proves (5.8). \(\square\)

Again the absolute child-ground condition fails maximally.  The
all-one child energy is \(m/2\), while \(x=s_S\) gives

\[
H_{A[S]}(s_S)=-\binom m2.
\tag{5.10}
\]

In fact, as \(S\) ranges over all type-balanced sets of all even
sizes, this example has

\[
\sum_r\binom{n/2}{r}^2=\binom n{n/2}
\]

strict local maxima, attaining Spink's general bound (2.6).
Therefore a fixed-point, LCP, or strict-local-maximum inverse theorem
which omits the absolute child-ground inequalities is impossible even
at the exact extremal count.

### Computational audit 5.4 (the first competitive balanced case)

For \(n=8\), an exhaustive enumeration of all \(2^{21}\) switching
classes (fixing the seven edges incident to vertex \(0\) positive)
finds \(4{,}200\) classes with \(M(A)=F(8)=10\).  Independently
re-running the enumeration and evaluating the reported witness by a
separate Python implementation gives

\[
\boxed{
\max_{\substack{M(A)=10\\z\in\{\pm1\}^8}}
\#\{S\in\tbinom{[8]}4:z_S\text{ is positive good}\}=4.
}                                                   \tag{5.11}
\]

The incidence threshold is

\[
2^{-4}\binom84=4.375,
\]

so an integer family would need at least five members.  One maximizing
witness has canonical edge bits \(53014\), switching mask \(24\), and
the four support masks

\[
71,\ 114,\ 141,\ 184.
\]

They form two complementary pairs and a two-dimensional affine
\(\mathbb F_2\)-flat.  This is only finite evidence, not an
asymptotic theorem, but it is consistent with the proposed
affine/type alternative and shows that the exact optimum already lies
strictly below the incidence threshold at the first nontrivial
balanced order.

## 6. The two canonical pair saturators fail the ground condition

Traffic alone has natural pair systems at exactly the incidence
entropy.  The next two propositions show that the complete versions
of both systems are killed by the child-ground condition.

### Proposition 6.1 (complete pair-union slices are impossible)

Partition \(2d\) vertices into pairs \(V_1,\ldots,V_d\).  Fix
\(2\le\ell\le d-2\), and suppose every union

\[
S_Q=\bigcup_{r\in Q}V_r,
\qquad Q\in\binom{[d]}\ell,
\]

is positive good.  Then no such signing exists.

#### Proof

Fix \(r\) and a vertex \(v\in V_r\).  When \(r\notin Q\), its
boundary equation is

\[
\sum_{s\in Q}\sum_{u\in V_s}a_{vu}=0.
\]

All \(\ell\)-subset sums of the \(d-1\) displayed block-row sums
vanish.  Taking differences of two subsets which differ by one
element shows that the block-row sums are all equal; their
\(\ell\)-fold sum is zero, so every one is zero.  Symmetry gives zero
row and column sums in every \(2\times2\) interpair block.  Hence

\[
A[V_r,V_s]
=c_{rs}
\begin{pmatrix}1&-1\\-1&1\end{pmatrix},
\qquad c_{rs}\in\{\pm1\}.                           \tag{6.1}
\]

For a selected pair \(V_r\), its cross-pair contribution to either
local field is zero.  Positive singleton stability therefore forces
its internal edge to be \(+1\).

Fix \(Q\).  The all-one child energy is \(\ell\), so the absolute
ground hypothesis says \(M(A[S_Q])=\ell\).  On every pair use instead
an antiuniform spin \(t_r(1,-1)\).  Its energy is

\[
-\ell+4H_{C[Q]}(t),
\qquad C=(c_{rs}).
\]

The bound

\[
\left|-\ell+4H_{C[Q]}(t)\right|\le\ell
\]

would force

\[
0\le H_{C[Q]}(t)\le\ell/2
\quad\text{for every }t.
\]

But a nonzero Boolean quadratic polynomial has mean zero and cannot
be nonnegative everywhere.  This is a contradiction. \(\square\)

### Proposition 6.2 (complete transversal cubes are impossible
asymptotically)

Again partition \(2d\) vertices into pairs
\(\{p_r,q_r\}\).  Suppose every transversal

\[
S_u=\{p_r:u_r=1\}\cup\{q_r:u_r=0\},
\qquad u\in\{0,1\}^d,
\]

is positive good.  Then the induced quotient signing \(C\) of order
\(d\) satisfies

\[
\boxed{M(C)=d/2.}                                   \tag{6.2}
\]

Consequently such a family cannot exist for all sufficiently large
\(d\), by the universal \(\Omega(d^{3/2})\) lower bound.

#### Proof

Hold \(u_r\) so that one vertex of pair \(r\) is outside, and vary
\(u_s\), \(s\ne r\), in its boundary equation.  The outside vertex
must have equal signs to \(p_s\) and \(q_s\).  Doing this from both
sides and using symmetry shows that every interpair block is
constant:

\[
A[\{p_r,q_r\},\{p_s,q_s\}]
=c_{rs}J_2,\qquad c_{rs}\in\{\pm1\}.                \tag{6.3}
\]

Let \(h_r=a_{p_rq_r}\).  The boundary equation of the unselected
vertex in pair \(r\) becomes

\[
h_r+\sum_{s\ne r}c_{rs}=0.                          \tag{6.4}
\]

The selected vertex has child local field
\(\sum_{s\ne r}c_{rs}=-h_r\).  Positive singleton stability forces
\(h_r=-1\), and hence

\[
\sum_{s\ne r}c_{rs}=1.                              \tag{6.5}
\]

Every transversal induces the same signing \(C=(c_{rs})\), and its
all-one vector is a positive absolute ground.  Equation (6.5) gives

\[
H_C(\mathbf 1)=\frac d2.
\]

Thus \(M(C)=d/2\), proving (6.2). \(\square\)

These propositions explain exactly where the ground condition enters:
the maximum-entropy pair families which saturate traffic counting are
not merely atypical; their complete forms are impossible.  A surviving
pair exception must therefore encode a strict subfamily in a quotient.
For the exact paired lift \(\mathcal L(C)\), that subfamily descends
by

\[
M(\mathcal L(C))=d+4M(C),\qquad
c(\mathcal L(C))=\sqrt2\,c(C)+O(d^{-1/2}).           \tag{6.6}
\]

Since \(2c_*>1/2\), \(c_*=0.336493364431\ldots\), two consecutive
exact pair lifts cannot lie on a minimizing sequence.  Thus any pair
exception produced by the missing inverse theorem has depth at most
one.

## 7. Constant-weight affine families are paired quotients

The affine branch of the inverse target admits an exact
classification independent of \(A\).

### Theorem 7.1 (Fourier classification and stability)

Let

\[
\mathcal V=v_0+W\subseteq\mathbb F_2^n
\]

be an affine subspace of dimension \(d\), and suppose every
\(v\in\mathcal V\) has the same Hamming weight.  Then

\[
\boxed{d\le n/2.}                                  \tag{7.1}
\]

More precisely, after parametrizing \(\mathcal V\) by
\(w\in\mathbb F_2^d\), every nonconstant coordinate has the form

\[
v_i(w)=\frac{1-\sigma_i(-1)^{\lambda_i\cdot w}}2,
\qquad
\lambda_i\in\mathbb F_2^d\setminus\{0\},
\quad \sigma_i\in\{\pm1\}.                         \tag{7.2}
\]

For each represented nonzero character \(\lambda\), its fibre

\[
V_\lambda=\{i:\lambda_i=\lambda\}
\]

has even size and contains equally many phases
\(\sigma_i=1\) and \(\sigma_i=-1\).  If \(q\) is the number of
represented nonzero characters and \(r\) the number of nonconstant
coordinates, then

\[
\boxed{
d\le q\le r/2\le n/2.
}                                                   \tag{7.3}
\]

Consequently:

1. if \(d=n/2\), there are no constant coordinates, there are exactly
   \(d\) character fibres, every fibre is one opposite-phase pair,
   and the represented characters form a basis.  Hence
   \(\mathcal V\) is, after an affine reparametrization, precisely the
   full transversal cube of \(n/2\) hidden coordinate pairs;
2. if \(d\ge n/2-r_0\), all but at most \(2r_0\) coordinates lie in
   opposite-phase pairs, and all but at most \(r_0\) represented
   character fibres can be chosen as independent basis fibres.

#### Proof

The coordinate functions on an affine \(\mathbb F_2\)-space are
constant or affine linear, giving (7.2).  Since the weight is
constant,

\[
\sum_i v_i(w)
=c+\frac r2
-\frac12
\sum_{\lambda\ne0}
\left(\sum_{i\in V_\lambda}\sigma_i\right)
(-1)^{\lambda\cdot w}
\tag{7.4}
\]

is independent of \(w\).  Independence of the nontrivial characters
forces

\[
\sum_{i\in V_\lambda}\sigma_i=0
\quad(\lambda\ne0).                                \tag{7.5}
\]

Every represented fibre therefore has at least two coordinates, with
equal phase counts, so \(q\le r/2\).  The coordinate functions
separate the \(2^d\) points of \(\mathcal V\), hence the represented
characters span \((\mathbb F_2^d)^*\), giving \(q\ge d\).  This proves
(7.3).

If \(d=n/2\), equality holds throughout (7.3).  Thus \(r=n\), every
fibre has size two, and the \(q=d\) spanning characters are a basis.
The two coordinates in each fibre have opposite phases, which is
exactly the transversal description.

Finally suppose \(d\ge n/2-r_0\).  The number of constant coordinates
plus the total excess

\[
\sum_\lambda(|V_\lambda|-2)
=r-2q
\]

is

\[
(n-r)+(r-2q)=n-2q\le n-2d\le2r_0.                 \tag{7.6}
\]

Thus at most \(2r_0\) coordinates are constant or belong to the
excess part of a fibre.  Since the \(q\) represented characters span
a \(d\)-space, a basis leaves only \(q-d\le n/2-d\le r_0\)
nonbasis fibres. \(\square\)

### Corollary 7.2 (the incidence-scale affine branch)

At the balanced layer \(m=n/2\), if an affine family satisfies

\[
|\mathcal V|
\ge2^{-m}\binom nm,
\]

then

\[
d=\log_2|\mathcal V|
\ge\frac n2-\frac12\log_2 n-O(1).                 \tag{7.7}
\]

Theorem 7.1 therefore makes it an
\(O(\log n)\)-defect paired transversal quotient: all but
\(O(\log n)\) coordinates occur in opposite-phase pairs, and only
\(O(\log n)\) pair characters are dependent on a chosen character
basis.

The apparent \(O(\log n)\) boundary margin can in fact be removed by
Fourier separation.

### Theorem 7.3 (incidence-scale affine good families are impossible)

For all sufficiently large even \(n\), there is no affine family

\[
\mathcal V\subseteq\binom{[n]}{n/2}
\]

of positive good sets satisfying

\[
|\mathcal V|\ge2^{-n/2}\binom n{n/2}.              \tag{7.8}
\]

#### Proof

Write \(d=\dim\mathcal V\) and choose, using Theorem 7.1, one
opposite-phase coordinate pair

\[
P_r=\{p_r,q_r\},
\qquad r=1,\ldots,d,
\]

whose characters form a basis.  Relabel the affine parameter
\(w\in\mathbb F_2^d\) so that

\[
\mathbf 1_{\{p_r\in S(w)\}}
=\frac{1+\chi_r(w)}2,\qquad
\mathbf 1_{\{q_r\in S(w)\}}
=\frac{1-\chi_r(w)}2,
\quad
\chi_r(w)=(-1)^{w_r}.                              \tag{7.9}
\]

Put all remaining coordinates in an exceptional set \(E\).  Its size
is

\[
h=n-2d=O(\log n)                                   \tag{7.10}
\]

by (7.7).  Every \(e\in E\) has a constant or affine membership
function

\[
\mathbf 1_{\{e\in S(w)\}}
=\frac{1+\sigma_e\chi_{\lambda_e}(w)}2
\]

when nonconstant.

Call a basis index \(r\) dirty if it belongs to
\(\operatorname{supp}\lambda_e\) for some exceptional coordinate
whose character has weight one or two.  There are at most

\[
b\le2h                                               \tag{7.11}
\]

dirty indices.  Delete their core pairs conceptually, and call the
remaining \(D=d-b\) pair indices clean.

We claim that for two distinct clean indices \(r,s\), the entire
interpair block is constant:

\[
\boxed{
A[P_r,P_s]=c_{rs}J_2,\qquad c_{rs}\in\{\pm1\}.
}                                                   \tag{7.12}
\]

Consider first the outside-field equation for row \(p_r\), restricted
to the hyperplane \(w_r=1\), on which \(p_r\notin S(w)\).  The
contribution of core pair \(s\) is

\[
\frac{a_{p_rp_s}+a_{p_rq_s}}2
+
\frac{a_{p_rp_s}-a_{p_rq_s}}2\chi_s(w).            \tag{7.13}
\]

On this hyperplane, an exceptional character
\(\chi_{\lambda_e}\) restricts to \(\pm\chi_{\lambda_e\setminus
\{r\}}\).  It can contribute to the singleton Fourier coefficient
\(\chi_s\) only when

\[
\lambda_e=e_s
\quad\text{or}\quad
\lambda_e=e_r+e_s.                                 \tag{7.14}
\]

Either possibility would make \(s\), or both \(r,s\), dirty.  Since
\(r,s\) are clean, Fourier uniqueness in the identity

\[
\sum_{i\in S(w)}a_{p_ri}=0
\quad(w_r=1)
\]

forces the singleton coefficient in (7.13) to vanish:

\[
a_{p_rp_s}=a_{p_rq_s}.                             \tag{7.15}
\]

The row \(q_r\) is outside on \(w_r=0\), and the identical argument
gives

\[
a_{q_rp_s}=a_{q_rq_s}.                             \tag{7.16}
\]

Applying (7.15)--(7.16) with \(r,s\) interchanged and using symmetry
of \(A\) makes the two row constants equal, proving (7.12).
Constant exceptional coordinates only affect the constant Fourier
coefficient and cause no exception to this argument.

Fix a clean \(r\) and use either vertex of \(P_r\) on the hyperplane
where it is outside.  Its boundary equation contains:

* its selected partner in \(P_r\);
* exactly one selected vertex from each other clean pair, contributing
  \(c_{rs}\);
* at most one selected vertex from each of the \(b\) dirty pairs;
* at most \(h\) selected exceptional coordinates.

Therefore

\[
\boxed{
\left|\sum_{\substack{s\ {\rm clean}\\s\ne r}}c_{rs}\right|
\le 1+b+h\le1+3h.
}                                                   \tag{7.17}
\]

Now fix any child \(S(w)\).  It contains one vertex from every core
pair and, because its weight is \(n/2\), exactly

\[
n/2-d=h/2
\]

exceptional vertices.  Separate its selected vertices into the
\(D\) clean representatives and a remainder \(R\) of size

\[
|R|\le b+h/2=O(h).
\]

By (7.12), the clean--clean all-one energy is

\[
\frac12\sum_{r\ {\rm clean}}
\sum_{\substack{s\ {\rm clean}\\s\ne r}}c_{rs}
=O(Dh)                                             \tag{7.18}
\]

using (7.17).  All edges incident to \(R\) contribute at most

\[
D|R|+\binom{|R|}{2}=O(Dh+h^2).                    \tag{7.19}
\]

Hence the all-one child energy obeys

\[
H_{A[S(w)]}(\mathbf 1)=O(n\log n).                 \tag{7.20}
\]

But \(\mathbf 1\) is, by hypothesis, a positive absolute ground of
the child.  The universal bound therefore gives

\[
H_{A[S(w)]}(\mathbf 1)
=M(A[S(w)])
\ge(c_*+o(1))(n/2)^{3/2},
\]

contradicting (7.20). \(\square\)

Thus the affine alternative in the inverse theorem is not merely a
bounded-depth exception: at the incidence threshold it is excluded
outright by the absolute child-ground inequality.  A surviving
inverse-Littlewood--Offord output must either have genuinely
non-affine code-separated entropy or produce a smaller quotient
before it reaches an affine constant-weight family.

## 8. Current exact inverse target

The remaining theorem can now be stated without generic
concentration language.

Let \(\mathcal F\) be a positive-good family of cardinality at least
\(2^{-m}\binom nm\).  Equivalently, every member satisfies

\[
A\mathbf 1_S\ge\mathbf 1_S,\qquad
\operatorname{supp}(A\mathbf 1_S)=S,               \tag{8.1}
\]

together with the full absolute-ground cut inequalities

\[
0\le e(R,S\setminus R)\le M(A[S])
\quad(R\subseteq S).                               \tag{8.2}
\]
\]

Prove one of:

1. a macroscopic negative interaction \(t\) from Theorem 4.1 yields
   a tangent-scale principal witness;
2. the zero/low-\(t\) crossing relations generate enough exact child
   grounds to form an affine family, which Theorem 7.3 now rules out;
3. the multivariate atom (5.2) forces repeated/opposite flat column
   types, after which the exact paired descent applies.

The principal unresolved obstruction is a code-separated family:
the incidence entropy does not force one-swap neighbors, and the
rank-only form (5.4) does not distinguish a generic high-rank atom
from a structured pair quotient.  Any next proof must use the
positive inequalities in (7.1) together with inverse
Littlewood--Offord structure; traffic equations alone are already
known to be false at the target threshold.
