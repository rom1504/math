# Action convergence / graphop audit for the Boolean quadratic norm

## 1. Normalization

Put \([n]\) under the uniform probability measure and let
\[
T_A=\frac{A}{\sqrt n}
\]
act by ordinary matrix multiplication.  With the normalized inner product,
\[
\langle f,T_Af\rangle
=\frac1n f^\top\frac A{\sqrt n}f
=\frac{f^\top Af}{n^{3/2}}.
\]
Hence the normalized objective is
\[
\Phi(T_A)
=\max_{f\in\{\pm1\}^n}|\langle f,T_Af\rangle|.
\tag{1}
\]

For a near-competitive signing, cube polarization gives
\[
\|A\|_{\infty\to1}^{\rm matrix}
\le2\max_x|x^\top Ax|=O(n^{3/2}),
\]
so the \(P\)-operator norm of \(T_A\) is uniformly bounded:
\[
\|T_A\|_{\infty\to1}
=\frac{\|A\|_{\infty\to1}^{\rm matrix}}{n^{3/2}}
=O(1).
\tag{2}
\]

## 2. What the primary action-convergence theorem gives

Backhausz and Szegedy define the \(k\)-profile of a \(P\)-operator from the
joint laws
\[
\mathcal D(f_1,\ldots,f_k,Tf_1,\ldots,Tf_k),
\qquad |f_i|\le1,
\]
and action convergence as Hausdorff convergence of all these profile sets.
They prove:

- every uniformly \(\|\cdot\|_{\infty\to1}\)-bounded sequence has an
  action-Cauchy subsequence;
- representation by an actual limit \(P\)-operator, with compactness of
  weak-equivalence classes, is proved under a stronger uniform
  \(p\to q\) bound with \(p<\infty\) and \(q>1\);
- action profiles imply quotient convergence under the same stronger
  bounds.

Primary source:
[Backhausz--Szegedy, *Action convergence of operators and graphs*,
Canadian Journal of Mathematics 74 (2022), 72--121](https://doi.org/10.4153/S0008414X2000070X),
especially Definitions 2.1--2.5 and Theorems 2.9--2.10.

The distinction at the endpoint \((\infty,1)\) is substantive.
Hrušková constructs action-convergent graph sequences having no other
uniform \(p\to q\) bound and shows that their limit operators can be
nonunique and fail the usual graphop properties:
[Hrušková, *Limits of action convergent graph sequences with unbounded
\((p,q)\)-norms*](https://arxiv.org/abs/2210.10720).

## 3. The Boolean objective is not action-continuous

The failure can be shown inside the class of sign matrices whose objective
is \(O(n^{3/2})\).

Let
\[
s=\lfloor n^{3/4}\rfloor,\qquad t=n-s.
\]
Choose:

1. a \(t\times t\) conference-type signing \(D\) such that
   \[
   M(D)\le(1/2+o(1))t^{3/2};
   \]
   replace \(D\) by \(-D\), if necessary, so that
   \(P(D)=M(D)\).  The universal lower bound gives
   \[
   P(D)\ge c_*t^{3/2}
   \]
   for an absolute \(c_*>0\);
2. any \(s\times s\) signing \(E\) with
   \[
   M(E)=O(s^{3/2});
   \]
3. an \(s\times t\) signing \(C\) with
   \[
   \|C\|_{\infty\to1}
   =O\!\left(\sqrt{st(s+t)}\right)
   =O(n^{11/8}).
   \tag{3}
   \]

Such a \(C\) exists by taking independent signs and applying Hoeffding's
inequality plus a union bound over the \(2^{s+t}\) row/column sign pairs.

Form
\[
A_n=\begin{pmatrix}E&C\\C^\top&D\end{pmatrix}
\]
and let \(B_n\) be obtained by replacing \(E\) with the all-positive
off-diagonal signing \(J_s-I_s\).

The triangle inequality gives
\[
M(A_n)
\le M(E)+M(D)+\|C\|_{\infty\to1}
\le(1/2+o(1))n^{3/2}.
\tag{4}
\]

For \(B_n\), take the all-one spin on the \(s\)-block and a positive endpoint
witness on the \(t\)-block.  Flipping every spin in the \(s\)-block preserves
the clique energy and reverses the cross term, so its sign can be made
nonnegative.  Therefore
\[
\begin{aligned}
M(B_n)
&\ge\binom s2+P(D)\\
&\ge(1/2+c_*-o(1))n^{3/2}.
\end{aligned}
\tag{5}
\]
Thus their normalized Boolean objectives stay a positive distance apart.

On the other hand,
\[
T_{A_n}-T_{B_n}
\]
is supported on the \(s\times s\) principal block.  For every fixed
\(k\)-tuple of bounded input functions, use the same functions for both
operators.  Their output coordinates agree outside this \(s\)-set.  Under
the uniform coordinate coupling, the two joint empirical laws therefore
differ with probability at most
\[
\frac sn=n^{-1/4}+o(1).
\]
Consequently, for every \(k\),
\[
d_H\bigl(S_k(T_{A_n}),S_k(T_{B_n})\bigr)
\le\frac sn,
\]
in both directions, and hence
\[
\boxed{d_M(T_{A_n},T_{B_n})\longrightarrow0.}
\tag{6}
\]

The objective gap is carried by output values of size
\[
s/\sqrt n=n^{1/4}
\]
on a set of probability
\[
s/n=n^{-1/4}.
\]
Their \(L^1\) contribution is order one, but their mass vanishes weakly.
This is exactly a failure of uniform integrability.  Therefore
\[
\boxed{\Phi\text{ is not continuous for action convergence under a
pure }\infty\to1\text{ bound}.}
\tag{7}
\]

## 4. The sign constraint and realization problem

The same pair \(A_n,B_n\) consists entirely of off-diagonal
\(\{\pm1\}\)-matrices, yet action convergence identifies their normalized
limits while (1) separates them.  Thus the microscopic sign constraint is
not retained in a way strong enough to determine the extremal objective.

There is a second, independent gap.  A subsequential compact limit does not
show that the same limit object is realizable by sign matrices at **every**
large order.  Backhausz--Szegedy explicitly note that not every graphop is
an action limit of finite graphs; the local-global subtheory already
contains non-realizable examples (Remark 3.4 of the cited paper).  No
sampling theorem in the action-convergence theory preserves simultaneously:

- exact off-diagonal \(\{\pm1\}\) entries;
- the \(1/\sqrt n\) normalization;
- arbitrary prescribed orders;
- the Boolean extremum (1).

## 5. Verdict

Action convergence is compact enough to extract profile subsequences, but
at the available \(\infty\to1\) endpoint it is not fine enough to retain
the normalized Boolean quadratic maximum.  The explicit \(n^{3/4}\)-vertex
planted block proves discontinuity inside the competitive sign-matrix
class.  Stronger \(2\to2\) bounds would restore uniform integrability, but
they are not uniform for all competitive signings, and even then the
all-orders sign-matrix realization theorem needed for convergence is
absent.
