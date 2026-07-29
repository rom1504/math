# Sharp vertex insertion as endpoint-cut discrepancy

## Normalization

Let

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
M=M(A)=\max_x|H_A(x)|
\]

and

\[
E(A)=\min_{b\in\{\pm1\}^n}\max_x
\bigl(|H_A(x)|+|b\cdot x|\bigr).
\]

Put

\[
\Delta_A(x)=M-|H_A(x)|.
\]

Then the insertion overhead is exactly the nonuniform discrepancy

\[
\boxed{
E(A)-M=
\min_{b\in\{\pm1\}^n}\max_x
\bigl(|b\cdot x|-\Delta_A(x)\bigr).
}
\tag{1}
\]

Thus a sharp insertion theorem is not a cardinality statement about
the extremal layer.  It is a simultaneous discrepancy theorem with one
slab for every energy state and with the energy deficit as the slab's
additional width.

## Endpoint-cut form

Switch a positive absolute maximizer to \(\mathbf1\), and call the
switched signing \(D\).  Thus

\[
H_D(\mathbf1)=M.
\]

For \(S\subset[n]\), let

\[
C_D(S)=\sum_{i\in S,\ j\notin S}d_{ij}
\]

and let \(y_S\) be \(1\) off \(S\) and \(-1\) on \(S\).  Then

\[
H_D(y_S)=M-2C_D(S).
\]

Absolute maximality is equivalent to

\[
0\le C_D(S)\le M\qquad(S\subset[n]).
\]

Consequently (1) becomes

\[
\boxed{
E(D)-M=
\min_b\max_{S\subset[n]}
\left(
|b\cdot y_S|
-2\min\{C_D(S),M-C_D(S)\}
\right).
}
\tag{2}
\]

The only relevant constraints are therefore the two endpoint layers of
a signed cut function which takes all its values in \([0,M]\).  This is
the structural form that a weighted discrepancy argument must exploit.

## Exact ground-state obstruction

Let

\[
\mathcal G(A)=\{x:|H_A(x)|=M\}.
\]

Equation (1) immediately gives

\[
\boxed{
E(A)-M\ge
\operatorname{disc}\mathcal G(A):=
\min_b\max_{x\in\mathcal G(A)}|b\cdot x|.
}
\tag{3}
\]

In particular, if \(\mathcal G(A)\) contains a full orthogonal Boolean
basis \(x^{(1)},\ldots,x^{(n)}\), then Parseval gives

\[
\sum_{k=1}^n\bigl(b\cdot x^{(k)}\bigr)^2=n^2
\]

for every Boolean \(b\), and hence

\[
E(A)-M\ge\sqrt n.
\tag{4}
\]

This shows that an \(O(\sqrt n)\) theorem, if true, has the correct
natural scale.  It also shows that the sharper derivative coefficient
cannot be obtained from a discrepancy theorem that ignores the signs
and relative deficits of the extremizers.

This obstruction does **not** automatically apply to conference
matrices.  In fact, a conference matrix of order \(n>2\) cannot have a
full orthogonal Boolean eigenbasis.  If it did, writing the Boolean
eigenbasis as a Hadamard matrix \(K\) would give

\[
C=\frac1nK^\top\operatorname{diag}(\lambda_1,\ldots,\lambda_n)K,
\qquad
\lambda_i\in\{\pm\sqrt{n-1}\}.
\]

If \(n-1\) is not a square, the entries cannot be integral.  If
\(n-1=r^2\), every off-diagonal entry would have the form

\[
C_{ij}=\frac r n\,z_{ij},\qquad z_{ij}\in\mathbb Z,
\]

so \(C_{ij}=\pm1\) would require \(z_{ij}=\pm n/r
=\pm(r+1/r)\), again impossible for \(r>1\).

## A pointwise mixed insertion inequality

There is a useful exact mixed-strategy statement.  For a spin \(x\)
with \(H_A(x)\ne0\), orient its signed local fields by

\[
\ell_j(x)=
\operatorname{sgn}(H_A(x))\,x_j(Ax)_j
\]

and put \(d=M-|H_A(x)|\).  Flipping coordinate \(j\) and using the
global bound \(M\) gives

\[
-\frac d2\le\ell_j(x)\le M-\frac d2,
\qquad
\sum_j\ell_j(x)=2(M-d).
\]

It follows that

\[
\sum_j|\ell_j(x)|
\le 2M+(n-2)d
\]

and hence

\[
\boxed{
\frac1n\sum_{j=1}^n
\left(
|H_A(x)|+|\ell_j(x)|
\right)
\le
M+\frac{2M}{n}-\frac{2d}{n}.
}
\tag{5}
\]

Take as a candidate new row the \(j\)-th row of \(A\), with its zero
diagonal entry filled by either sign.  Its field at \(x\) has absolute
value at most \(|\ell_j(x)|+1\).  Therefore the uniform distribution
over these \(n\) duplicate-row insertions has, against every fixed
state \(x\), expected payoff at most

\[
M+\frac{2M}{n}+1-\frac{2\Delta_A(x)}n.
\tag{6}
\]

The quantifier order is the sole gap: (6) is a mixed insertion theorem,
whereas \(E(A)\) requires one pure row that works for every state.
Purifying (6) is a concrete joint-level-set problem.  A scalar entropy
profile does not encode enough information to do it.

## Exact low-temperature certificate

Define

\[
W_A(\beta)=\sum_{x,\sigma=\pm1}e^{\beta\sigma H_A(x)}
=2\sum_x\cosh(\beta H_A(x)).
\]

For the extension by \(b\),

\[
W_{A,b}^{\rm ext}(\beta)
=\sum_{x,y,\sigma}
e^{\beta\sigma(H_A(x)+y\,b\cdot x)}
=4\sum_x\cosh(\beta H_A(x))\cosh(\beta b\cdot x).
\]

Averaging over a uniform Boolean row \(b\) gives

\[
\mathbb E_b W_{A,b}^{\rm ext}(\beta)
=2(\cosh\beta)^nW_A(\beta).
\]

Thus some row satisfies

\[
\boxed{
E(A)
\le
\frac1\beta
\left(
\log W_A(\beta)+n\log\cosh\beta+\log2
\right).
}
\tag{7}
\]

At \(\beta=a/\sqrt n\), (7) yields an \(O(\sqrt n)\) insertion
overhead whenever

\[
\log W_A(a/\sqrt n)-\frac{aM}{\sqrt n}=O(1).
\]

This is only a sufficient condition, not a proved property of
minimizers.  It identifies the needed statement as a frozen
low-temperature bound at the \(\sqrt n\) energy-gap scale.

## Current frontier

The most direct missing lemma is the following weighted endpoint-cut
statement: for every asymptotically minimizing cut-positive signing
\(D\), there is a Boolean \(b\) such that

\[
|b\cdot y_S|
\le
2\min\{C_D(S),M-C_D(S)\}+O(\sqrt n)
\quad\text{for every }S.
\]

Equation (2) shows that this is exactly an \(O(\sqrt n)\) insertion
theorem.  The mixed bound (5) suggests trying to purify duplicate-row
insertions using the joint geometry of their bad-coordinate sets.
The orthogonal-basis calculation (3)--(4) shows why no argument can
hope for a scale below \(\sqrt n\) in general.
