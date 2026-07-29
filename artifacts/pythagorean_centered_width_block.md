# Nonlinear block lower bounds for centered width

Let
\[
A=\begin{pmatrix}B&C\\ C^\top&D\end{pmatrix},
\qquad
H_A(y,z)=H_B(y)+H_D(z)+y^\top Cz.
\]
Write
\[
P(E)=\max H_E,\qquad Q(E)=-\min H_E,\qquad
W(E)=\frac{P(E)+Q(E)}2.
\]

## 1. Exact block endpoint formula

For fixed \(y,z\), replacing \(y\) by \(-y\) preserves \(H_B(y)\) and
\(H_D(z)\) and reverses \(y^\top Cz\).  Therefore
\[
\boxed{
P(A)=\max_{y,z}
\left(H_B(y)+H_D(z)+|y^\top Cz|\right),}
\tag{1}
\]
and
\[
\boxed{
Q(A)=\max_{y,z}
\left(-H_B(y)-H_D(z)+|y^\top Cz|\right).}
\tag{2}
\]

Put
\[
h(y,z)=H_B(y)+H_D(z),
\quad
P_0=P(B)+P(D),
\quad
Q_0=Q(B)+Q(D),
\]
\[
W_0=W(B)+W(D)=\frac{P_0+Q_0}{2},
\]
and define the two endpoint slacks
\[
s_+(y,z)=P_0-h(y,z),\qquad
s_-(y,z)=Q_0+h(y,z).
\]
Equations (1)--(2) give the exact identity
\[
\boxed{
W(A)=W_0+\frac12\left[
\max_{y,z}\bigl(|y^\top Cz|-s_+(y,z)\bigr)
+
\max_{y,z}\bigl(|y^\top Cz|-s_-(y,z)\bigr)
\right].}
\tag{3}
\]
Both maxima in (3) are nonnegative: at a positive or negative endpoint of
the block-diagonal Hamiltonian the corresponding slack is zero.

Thus cross energy increases centered width only when it exceeds the
appropriate internal energy-layer slack.

## 2. Exact counterexample to every Pythagorean gain

Consider the signing
\[
A=
\begin{pmatrix}
0& 1&-1& 1&-1\\
1& 0& 1&-1&-1\\
-1&1& 0&-1&-1\\
1&-1&-1& 0&-1\\
-1&-1&-1&-1&0
\end{pmatrix}
\]
and split it after the first vertex.  Thus \(B=(0)\), \(C\) is the row
\[
C=(1,-1,1,-1),
\]
and \(D\) is the lower-right \(4\times4\) block.

Direct enumeration gives
\[
\operatorname{range}H_D=[-4,4],\qquad W(D)=4,
\]
\[
\|C\|_{\infty\to1}=4,
\]
and nevertheless
\[
\operatorname{range}H_A=[-4,4],\qquad W(A)=4.
\]
In fact, for every \(y\in\{\pm1\}^4\),
\[
\boxed{
|C y|
\le
\min\{4-H_D(y),\,4+H_D(y)\}.}
\tag{4}
\]
The cross field fits entirely inside the positive and negative endpoint
slacks, so both correction terms in (3) vanish.

Consequently, for every \(c>0\), the proposed inequality
\[
W(A)^2\ge
\bigl(W(B)+W(D)\bigr)^2
+c\|C\|_{\infty\to1}^2
\]
is false.  The strongest uniform coefficient is
\[
\boxed{c_{\rm uniform}=0.}
\]

## 3. Strongest scalar consequence

The usual two lower bounds remain
\[
W(A)\ge W(B)+W(D)
\]
and
\[
W(A)\ge\|C\|_{\infty\to1}.
\]
Equation (3) shows why no strictly stronger function of the three scalar
quantities can be obtained without controlling the joint layer profile
\[
\bigl(h(y,z),|y^\top Cz|\bigr).
\]
Large cross discrepancy can be supported entirely on internal configurations
with enough endpoint slack, as (4) demonstrates.

## 4. Verdict

The desired uniformly convex/Pythagorean correction is blocked by exact
energy-layer anti-alignment, already at order five.  Any nonlinear block
lower bound capable of ruling out slow oscillation must retain at least the
two functions
\[
s_\pm(y,z)
\]
or an entropy/distributional summary strong enough to prevent
\(|y^\top Cz|\) from hiding under both of them.  The scalars
\(W(B),W(D)\), and \(\|C\|_{\infty\to1}\) are insufficient.

