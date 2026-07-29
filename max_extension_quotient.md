# Max-extension quotient for proportional restriction

## 1. Exact block identities

Write a signing in block form

\[
A=
\begin{pmatrix}
B&C\\
C^\top&D
\end{pmatrix},
\]

where \(B=A[S]\) and \(D=A[S^c]\). For \(y\in\{\pm1\}^S\), put

\[
h_y=C^\top y.
\]

Then

\[
H_A(y,z)=H_B(y)+h_y\cdot z+H_D(z).
\]

The positive and negative max-extension profiles are

\[
G_S^+(y)=\max_zH_A(y,z),\qquad
G_S^-(y)=-\min_zH_A(y,z),
\]

and the extension gaps are therefore

\[
\boxed{
D_S^+(y)
=G_S^+(y)-H_B(y)
=\max_z\{H_D(z)+h_y\cdot z\},}
\tag{1}
\]

\[
\boxed{
D_S^-(y)
=G_S^-(y)+H_B(y)
=\max_z\{-H_D(z)-h_y\cdot z\}.}
\tag{2}
\]

Since \(H_D(-z)=H_D(z)\), these can be written

\[
D_S^+(y)=\max_z\{H_D(z)+|h_y\cdot z|\},
\]

\[
D_S^-(y)=\max_z\{-H_D(z)+|h_y\cdot z|\}.
\tag{3}
\]

Thus the max-extension quotient is exactly the external-field support
profile already encountered in the rooted transfer recursion.

## 2. Exact lower bounds and centered/oriented split

Let

\[
P_D=\max_zH_D(z),\qquad Q_D=-\min_zH_D(z).
\]

From (3),

\[
\boxed{
D_S^+(y)\ge
\max\{P_D,\ \|C^\top y\|_1-Q_D\},}
\tag{4}
\]

\[
\boxed{
D_S^-(y)\ge
\max\{Q_D,\ \|C^\top y\|_1-P_D\}.}
\tag{5}
\]

Adding the two separate maxima after evaluating them on the same \(z\) gives

\[
\boxed{
D_S^+(y)+D_S^-(y)\ge2\|C^\top y\|_1.}
\tag{6}
\]

Equivalently, the centered extension range is at least the rectangular
field norm. Moreover, viewing \(H_D\) as a perturbation of the symmetric
linear form \(h_y\cdot z\) gives

\[
|D_S^+(y)-D_S^-(y)|\le2M(D),
\tag{7}
\]

and hence again

\[
D_S^\pm(y)\ge\|C^\top y\|_1-M(D).
\tag{8}
\]

The scale-transfer target for a positive restricted state is the *oriented*
bound

\[
D_S^+(y)\ge
(1-\alpha^{3/2})M(A)-o(n^{3/2}),
\tag{9}
\]

not merely a bound on the sum in (6). The outside quadratic form can shift
the entire rectangular gain into the wrong orientation.

## 3. Ground-state extensions

If \(y\) is the restriction of a positive global ground state \(x^\star\),
then its matching outside spin is an admissible extension attaining
\(M(A)\). Consequently,

\[
\boxed{
D_S^+(y)=M(A)-H_B(y).}
\tag{10}
\]

There is an analogous identity on the negative side. Thus ground-state
restrictions have the perfect extension gap automatically.

The unresolved issue has a sharp interpretation:

> Can one choose \(S\) so that every high restricted state has a near-ground
> full extension?

This is more precise than controlling all complement extensions, but it is a
nonlocal ground-state compatibility statement.

## 4. Variance bound is too small

For fixed \(y\), the outside affine chaos

\[
W_y(z)=H_D(z)+h_y\cdot z
\]

has

\[
\mathbb EW_y=0,\qquad
\mathbb EW_y^2=\binom{|S^c|}{2}+\|h_y\|_2^2.
\tag{11}
\]

Degree-two hypercontractivity implies that both its positive and negative
maxima are at least a universal constant times the square root of (11).
For a typical flat cross block, however,
\(\|h_y\|_2=\Theta(n)\), so this gives only \(O(n)\), whereas the required
gap is \(\Theta(n^{3/2})\). The needed gain comes from the
\(\ell_1\) field norm in (4)--(6), not from variance.

## 5. Sharp orientation/kernel obstruction for a fixed block

Fix a high state \(y\) of \(B\). If \(|S|\) is even, the cross columns can be
chosen so that each is balanced against \(y\):

\[
C^\top y=0.
\tag{12}
\]

Then

\[
\boxed{
D_S^+(y)=P_D,\qquad D_S^-(y)=Q_D.}
\tag{13}
\]

This construction is compatible with all three blocks having
\(O(n^{3/2})\) Boolean norm and \(O(\sqrt n)\) operator norm: take
conference/random-like \(B,D\), and choose balanced pseudorandom columns for
\(C\). Thus neither spectral regularity nor a large *global* rectangular
discrepancy forces the needed gap at a specified high state.

More generally, (8) shows that a lower bound through the cross block must
prove a restricted singular-value statement on the high-energy layer:

\[
\inf_{\{y:H_B(y)\text{ high}\}}\|C^\top y\|_1
\quad\text{is large.}
\tag{14}
\]

Ordinary Grothendieck--Pietsch factorization is an upper-norm theorem and
cannot provide (14); a rectangular signing can have a large
\(\infty\to1\) norm while annihilating a prescribed Boolean vector.

## 6. Coupling the two restricted extrema

Let \(u\) maximize \(H_B\), and let \(v\) minimize it. Evaluating the two
correctly oriented extension maxima on the same outside spin \(z\) cancels
\(H_D(z)\), while evaluating the second on \(-z\) gives the other parity.
Thus

\[
\boxed{
D_S^+(u)+D_S^-(v)
\ge
R_C(u,v),}
\tag{15}
\]

where

\[
R_C(u,v)=
\max\left\{
\|C^\top(u-v)\|_1,\,
\|C^\top(u+v)\|_1
\right\}.
\tag{16}
\]

Since every extension is bounded by \(K=M(A)\),

\[
P_B+D_S^+(u)\le K,\qquad
Q_B+D_S^-(v)\le K.
\]

Consequently,

\[
\boxed{
P_B+Q_B\le2K-R_C(u,v),\qquad
W(B)\le K-\frac12R_C(u,v).}
\tag{17}
\]

This is an exact centered-range restriction theorem. It does not bound
\(M(B)=\max(P_B,Q_B)\): the entire cross gain can be absorbed by the gap
belonging to the non-dominant orientation.

There is a sharp kernel obstruction even after coupling. Choose every cross
column \(c_j\) to satisfy

\[
c_j\cdot u=c_j\cdot v=0.
\tag{18}
\]

Such sign columns exist whenever the joint sign cells of \(u,v\) have the
appropriate even parities; random columns conditioned on these two linear
constraints can still have regular spectral and rectangular behavior. Then

\[
C^\top u=C^\top v=0,\qquad R_C(u,v)=0,
\]

and the correctly oriented gaps reduce to \(P_D,Q_D\).

The smallest concrete seed is the all-negative signing \(B_4\). It has a
positive maximizer

\[
u=(-1,-1,1,1)
\]

and a negative maximizer

\[
v=(-1,-1,-1,-1).
\]

The sign column

\[
c=(1,-1,1,-1)
\]

is orthogonal to both. Thus even coupling the two exact restricted extrema
does not force a cross extension gain.

## 7. Convex dual form

Define the finite lifted set

\[
\mathcal V_{C,D}
=\{(Cz,H_D(z)):z\in\{\pm1\}^{S^c}\}.
\]

Then

\[
D_S^+(y)=h_{\operatorname{conv}\mathcal V_{C,D}}(y,1),
\qquad
D_S^-(y)=h_{\operatorname{conv}\mathcal V_{C,D}}(-y,-1),
\tag{19}
\]

where \(h_K\) is the support function. This is an exact convex dual, but it
also shows why a scalar recursion does not close: the whole lifted polytope,
equivalently the full external-field profile, is required.

## 8. Verdict

The max-extension quotient removes the artificial extension redundancy from
the selector process and gives the right exact state. It does not yet yield a
scale-transfer inequality. The centered rectangular gain is automatic, but
the needed oriented gain can vanish on a high state lying in the cross
kernel.

A successful theorem must select \(S\) so that its high positive and negative
energy layers simultaneously avoid the corresponding oriented
cross-kernel obstructions. Coupling the two extrema yields the exact centered
bound (17), but the codimension-two construction (18) shows that it does not
repair orientation. This is not implied by operator norm, cut norm, or the
marginal energy profile, and is the exact residual subset-selection problem
left by this route.
