# Averaged anisotropic selector criterion and hub obstruction

## 1. Exact averaged joint-tail criterion

Let

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad K=M(A),
\]

and fix \(m,L,t\) with \(0\le t<L<K\). Put

\[
\delta=\frac{L-t}{K-t}.
\]

For a uniform \(m\)-subset \(S\) and an independent uniform full spin \(X\),
define

\[
p_+(t,L)=
\mathbb P\{H_A(X)\ge t,\ H_{A[S]}(X_S)\ge L\},
\]

\[
p_-(t,L)=
\mathbb P\{H_A(X)\le-t,\ H_{A[S]}(X_S)\le-L\}.
\]

If \(b_+(L)\) and \(b_-(L)\) are the fractions of \(m\)-supports whose
positive or negative restricted maxima are at least \(L\), then

\[
\boxed{
b_+(L)\le\delta^{-1}2^m p_+(t,L),\qquad
b_-(L)\le\delta^{-1}2^m p_-(t,L).}
\tag{1}
\]

Consequently,

\[
\boxed{
p_+(t,L)+p_-(t,L)<\delta2^{-m}
\quad\Longrightarrow\quad
\exists\,|S|=m:\ M(A[S])<L.}
\tag{2}
\]

To prove (1), select one restricted witness \(y_S\) for each positive-bad
support. Uniformly extending \(y_S\) has mean energy at least \(L\) and energy
at most \(K\), so at least a \(\delta\)-fraction of its extensions have full
energy at least \(t\). Sum over bad supports, then enlarge the selected
witness event to the joint event in the definition of \(p_+\). The negative
case is identical.

The universal joint second moments are

\[
\mathbb E H_A(X)^2=\binom n2,\qquad
\mathbb E H_{A[S]}(X_S)^2=\binom m2,
\]

\[
\mathbb E\!\left[H_A(X)H_{A[S]}(X_S)\right]=\binom m2.
\tag{3}
\]

Thus the two coordinates have asymptotic correlation \(\alpha=m/n\), but
these second moments do not determine the required exponential joint tail.

## 2. A \(\sqrt n\)-hub obstruction to any uniform joint exponent

For every fixed \(\alpha\in(0,1)\), there are complete signings with
\(M(A)=O(n^{3/2})\) for which the joint probability in (2), at
\(L=\alpha^{3/2}M(A)\), is exponentially larger than \(2^{-\alpha n}\).

Choose a set \(R\) of

\[
r=\lfloor\gamma\sqrt n\rfloor
\]

hub vertices, and put \(T=[n]\setminus R\). Set every edge incident to a hub
to \(+1\), including edges internal to \(R\). On \(T\), choose any signing
\(B\) with

\[
M(B)\le Cn^{3/2};
\]

such signings exist with an absolute \(C\) by the elementary random-signing
union bound.

For \(s_R=\sum_{i\in R}x_i\) and \(s_T=\sum_{j\in T}x_j\),

\[
H_A(x)
=\frac{s_R^2-r}{2}+s_Rs_T+H_B(x_T).
\tag{4}
\]

It follows that

\[
K=M(A)\le(\gamma+C+o(1))n^{3/2}.
\tag{5}
\]

Take \(m=\alpha n+o(n)\) and condition a uniform \(m\)-subset \(S\) to contain
all hubs. This costs only

\[
\frac{\binom{n-r}{m-r}}{\binom nm}
=2^{-O(\sqrt n)}.
\tag{6}
\]

Set all hub spins to \(+1\), require the selected non-hub spins to have bias
at least

\[
v_\gamma
=\frac{\alpha^{3/2}(\gamma+C)+C}{\alpha\gamma}+o(1),
\tag{7}
\]

and require the unselected non-hub spins to have nonnegative sum. If

\[
\gamma>
\frac{C(1+\alpha^{3/2})}{\alpha-\alpha^{3/2}},
\tag{8}
\]

then \(v_\gamma<1\). The hub cross term, even after a worst-case contribution
\(-Cn^{3/2}\) from \(B\), guarantees simultaneously

\[
H_{A[S]}(X_S)\ge\alpha^{3/2}K,\qquad
H_A(X)\ge\alpha^{3/2}K.
\tag{9}
\]

The hub-spin condition costs \(2^{-O(\sqrt n)}\); the outside nonnegative-sum
condition costs only a constant. The selected-spin bias has probability

\[
2^{-\left(
\alpha\left[
1-h_2\!\left(\frac{1-v_\gamma}{2}\right)
\right]+o(1)\right)n}.
\tag{10}
\]

Letting \(\gamma\to\infty\), its rate tends to

\[
\boxed{
\alpha\left[
1-h_2\!\left(\frac{1-\sqrt\alpha}{2}\right)
\right]<\alpha.}
\tag{11}
\]

Thus the joint probability is much larger than \(2^{-\alpha n}\). Since the
event in (9) also implies \(H_A(X)\ge t\) for every \(t<L\), optimizing the
lower threshold in (2) does not remove the obstruction.

## 3. Interpretation

The counterexample does not refute proportional restriction: all hubs occupy
only \(O(\sqrt n)=o(n)\) vertices and can be avoided by the restricted set.
It refutes a *uniform* averaged joint-tail theorem that ignores where the
large selector leverage is located.

The correct next structural branch is therefore:

> remove \(o(n)\) high-leverage vertices, then prove the exponential averaged
> selector bound on the remaining principal submatrix.

This is the averaged analogue of the operator-norm obstruction in the
max-multiplicity argument.

