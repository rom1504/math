# Correlated chaining and block amplification

Checkpoint date: 2026-07-25.

## 1. Setup and verdict

Use the one-copy normalization

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
M(A)=\max_x|H_A(x)|.
\]

For

\[
G=\begin{pmatrix}A&B\\B^\top&D\end{pmatrix},
\]

the block objective is exactly

\[
M(G)
=\max_{x,y}
\left(
|H_A(x)+H_D(y)|+|x^\top By|
\right). \tag{1.1}
\]

When \(B\) is an independent random sign matrix, the exponential union
bound can be replaced by an exact **polynomial energy-layer
reduction**.  There are only \(O(n^2m^2)\) pairs of possible internal
energy values, and the correlated bilinear supremum concentrates on
each pair.  Consequently

\[
\mathbb E_BM(G)
=
\max_{p,q}
\left(
|p+q|+W_B(\mathcal X_A(p),\mathcal X_D(q))
\right)
+o((n+m)^{3/2}), \tag{1.2}
\]

where

\[
\mathcal X_A(p)=\{x:H_A(x)=p\}
\]

and

\[
W_B(X,Y)=
\mathbb E_B\max_{x\in X,y\in Y}|x^\top By|.
\]

This is a proof-level sufficient width condition for amplification,
given in Section 4.

There is also a decisive limitation.  Every complete-graph quadratic
Hamiltonian has an exact energy level \(p=O(n)\) containing a
\(1/\operatorname{poly}(n)\) fraction of the cube.  Such polynomially
dense levels have asymptotically full Gaussian width.  More strongly,
for an independent random sign cross block, restricting both sides to
these bulk levels changes the full rectangular Boolean norm by only
\(o((n+m)^{3/2})\).  Therefore

\[
\boxed{
\mathbb E_BM(G)
\ge
\mathbb E_B\|B\|_{\infty\to1}
-o((n+m)^{3/2}).
}
\tag{1.3}
\]

No chaining argument can use the child energies to lower this bulk
floor.  Independent random cross blocks can prove same-constant
amplification only if their own bipartite ground-energy constant is at
most the target constant.  The remaining viable construction is a
flat structured block whose rectangular norm is low and whose
high-energy layer widths satisfy (4.3); independence by itself gives
no purification.

## 2. Exact polynomial layer reduction

The possible values of \(H_A\) lie in

\[
\left[-\binom n2,\binom n2\right]
\]

and have fixed parity.  Hence there are at most \(O(n^2)\) nonempty
sets \(\mathcal X_A(p)\), and similarly \(O(m^2)\) for \(D\).

For a fixed nonempty pair \((p,q)\), put

\[
S_{p,q}(B)
=\max_{\substack{x\in\mathcal X_A(p)\\
                 y\in\mathcal X_D(q)}}
|x^\top By|.
\]

Since each energy layer is invariant under \(x\mapsto-x\), the
absolute value can equivalently be removed from this supremum.
Equation (1.1) becomes the exact finite identity

\[
M(G)=\max_{p,q}\bigl(|p+q|+S_{p,q}(B)\bigr). \tag{2.1}
\]

Changing one entry of a Rademacher \(B\) changes every \(S_{p,q}\) by
at most \(2\).  Bounded-difference concentration and a union bound over
only \(O(n^2m^2)\) layers give

\[
\mathbb E\max_{p,q}
\left(S_{p,q}-\mathbb ES_{p,q}\right)
=O\!\left(\sqrt{nm\log(n+m)}\right). \tag{2.2}
\]

The same statement follows from Gaussian concentration when \(B\) is
Gaussian.  Since

\[
\mathbb E\max_{p,q}Z_{p,q}\ge\max_{p,q}\mathbb EZ_{p,q},
\]

(2.1)--(2.2) prove the two-sided asymptotic formula (1.2), uniformly
when \(n,m\) are comparable.

This is the correct replacement for the failed state-by-state union
bound: correlations are retained exactly inside each energy layer.

## 3. Canonical metric and a marginal-width upper bound

First take \(B\) Gaussian.  For

\[
Z_{x,y}=x^\top By,
\]

the canonical metric is

\[
d\bigl((x,y),(x',y')\bigr)^2
=2nm-2(x\cdot x')(y\cdot y'). \tag{3.1}
\]

For \(X\subseteq\{\pm1\}^n\), define its Gaussian width

\[
w(X)=\mathbb E_g\max_{x\in X}g\cdot x.
\]

Compare \(Z\) with

\[
U_{x,y}=\sqrt m\,g\cdot x+\sqrt n\,h\cdot y.
\]

Writing \(a=x\cdot x'\), \(b=y\cdot y'\), the increment variances obey

\[
\mathbb E(U_{x,y}-U_{x',y'})^2
-\mathbb E(Z_{x,y}-Z_{x',y'})^2
=2(n-a)(m-b)\ge0. \tag{3.2}
\]

Sudakov--Fernique therefore gives the constant-one comparison

\[
\boxed{
W_{\rm Gauss}(X,Y)
\le\sqrt m\,w(X)+\sqrt n\,w(Y).
}
\tag{3.3}
\]

Equivalently, one may use the full generic-chaining functional

\[
W_{\rm Gauss}(X,Y)\asymp
\gamma_2(X\otimes Y,d)
\]

or Dudley's metric-entropy integral in the metric (3.1).  Equation
(3.3) is a simple explicit sufficient bound which already exploits the
strong tensor correlations.

For a Rademacher cross block, conditioning a Gaussian coefficient into
its sign and magnitude gives the general comparison

\[
W_{\rm Rad}(X,Y)
\le\sqrt{\frac\pi2}\,
W_{\rm Gauss}(X,Y), \tag{3.4}
\]

and hence

\[
W_{\rm Rad}(X,Y)
\le\sqrt{\frac\pi2}
\left(\sqrt m\,w(X)+\sqrt n\,w(Y)\right). \tag{3.5}
\]

The factor in (3.4) is too expensive for a sharp constant.  A
constant-sharp Rademacher proof must estimate the actual correlated
width \(W_B(X,Y)\), or prove a low-influence invariance principle for
the relevant energy layers.

## 4. Proof-level sufficient profile condition

For an independent Rademacher block, (1.2) gives the exact sufficient
condition

\[
\boxed{
\max_{p,q}
\left[
|p+q|+W_B(\mathcal X_A(p),\mathcal X_D(q))
\right]
\le c(n+m)^{3/2}+o((n+m)^{3/2}).
}
\tag{4.1}
\]

Then

\[
\mathbb E_BM(G)\le c(n+m)^{3/2}+o((n+m)^{3/2}).
\]

For Gaussian \(B\), the more checkable sufficient condition is

\[
\boxed{
\max_{p,q}
\left[
|p+q|
+\sqrt m\,w(\mathcal X_A(p))
+\sqrt n\,w(\mathcal X_D(q))
\right]
\le c(n+m)^{3/2}+o((n+m)^{3/2}).
}
\tag{4.2}
\]

Let \(N=n+m\), \(n/N\to\alpha\), \(m/N\to\beta\).  If

\[
p=u\,n^{3/2}+o(n^{3/2}),\qquad
w(\mathcal X_A(p))=\omega_A(u)n+o(n),
\]

and similarly for \(D\), then (4.2) becomes the dimensionless profile
condition

\[
\boxed{
\left|\alpha^{3/2}u+\beta^{3/2}v\right|
+\alpha\sqrt\beta\,\omega_A(u)
+\beta\sqrt\alpha\,\omega_D(v)
\le c
}
\tag{4.3}
\]

for every pair of nonempty limiting energy layers.

This is a genuinely lossy asymptotic profile: unlike the exact
external-field support function, it records only the Gaussian widths
of equal-energy layers.

## 5. Every quadratic Hamiltonian has a maximal-width bulk level

For uniform \(x\),

\[
\mathbb EH_A(x)^2=\binom n2. \tag{5.1}
\]

Taking, for example, \(K=2\), Chebyshev gives

\[
\left|
\left\{x:|H_A(x)|\le Kn\right\}
\right|
\ge c_K2^n \tag{5.2}
\]

with an absolute \(c_K>0\).  The interval \([-Kn,Kn]\) contains only
\(O(n)\) possible energy values.  Hence some exact value \(p_A\) with

\[
|p_A|\le Kn
\]

satisfies

\[
\frac{|\mathcal X_A(p_A)|}{2^n}
\ge\frac{c_K}{Cn}. \tag{5.3}
\]

The same holds for \(D\).

There is a useful width lemma.  If \(X\subseteq\{\pm1\}^n\) has cube
density \(\pi\), then

\[
\boxed{
w(X)
\ge
\sqrt{\frac2\pi}
\left(
n-2\sqrt{\frac n2\log\frac1\pi}
\right).
}
\tag{5.4}
\]

To prove it, let \(U=\operatorname{sign}(g)\), choose an unweighted
nearest point \(x(U)\in X\), and use independence of \(U\) and the
magnitudes \((|g_i|)\):

\[
\mathbb E[g\cdot x(U)]
=\sqrt{\frac2\pi}
\left(n-2\mathbb E d_H(U,X)\right).
\]

The cube transport/isoperimetric inequality gives

\[
\mathbb E d_H(U,X)
\le\sqrt{\frac n2\log\frac1\pi}.
\]

Applying (5.4) to (5.3),

\[
\boxed{
w(\mathcal X_A(p_A))
=\left(\sqrt{\frac2\pi}+o(1)\right)n.
}
\tag{5.5}
\]

Thus complete-graph quadratic structure does not force a narrow
zero-energy layer.  It forces the opposite: at least one
\(o(n^{3/2})\)-energy layer has asymptotically the full cube width.

At a balanced split, inserting (5.5) into the right side of (4.2)
gives the bulk marginal-width scale

\[
\frac1{\sqrt\pi}=0.564189\ldots
\]

after normalization by \((2n)^{3/2}\).  This is an upper comparison,
not a lower bound on the true bipartite width, but it shows that the
simple marginal-width sufficient criterion cannot reach any target
below \(1/\sqrt\pi\).

## 6. Bulk preservation under an independent random cross block

The stronger obstruction uses the true correlated width.

Let \(X\subseteq\{\pm1\}^n\), \(Y\subseteq\{\pm1\}^m\) have densities
at least inverse-polynomial.  Let

\[
R(B)=\|B\|_{\infty\to1}
=\max_{x,y}x^\top By.
\]

Choose a maximizer \((X_B,Y_B)\) uniformly among all maximizing pairs.
Row and column switching invariance of an iid Rademacher matrix implies
that each marginal \(X_B,Y_B\) is uniform on its cube.

Let \(x'\in X,y'\in Y\) be nearest Hamming points to these marginals.
Then

\[
R(B)-x'^\top By'
\le
2\|B\|_{\mathrm{op}}
\left(
\sqrt{m\,d_H(X_B,X)}
+\sqrt{n\,d_H(Y_B,Y)}
\right). \tag{6.1}
\]

The cube transport bound gives

\[
\mathbb E d_H(X_B,X)=O(\sqrt{n\log n}),
\qquad
\mathbb E d_H(Y_B,Y)=O(\sqrt{m\log m}). \tag{6.2}
\]

The standard random-matrix estimate

\[
\mathbb E\|B\|_{\mathrm{op}}^2
=O(n+m)
\]

and Cauchy--Schwarz in (6.1) imply, for comparable \(n,m\),

\[
\boxed{
\mathbb E\max_{x\in X,y\in Y}x^\top By
\ge
\mathbb E R(B)-o((n+m)^{3/2}).
}
\tag{6.3}
\]

Take \(X=\mathcal X_A(p_A)\), \(Y=\mathcal X_D(p_D)\) from Section 5.
Their internal energies are only \(O(n+m)\), so (1.1) and (6.3) prove
(1.3).

The argument also applies to a fixed flat block randomized by
independent row and column sign switches, provided its operator norm
is \(O(\sqrt n+\sqrt m)\).  Its rectangular norm is switching
invariant, and the random switched maximizer has uniform marginals.
Thus a randomized Hadamard block retains its full rectangular norm on
the polynomially dense bulk energy layers, up to \(o(N^{3/2})\).

## 7. Consequences for amplification

Define the balanced iid bipartite sequence by

\[
\gamma_{{\rm bip},n}
=\frac{\mathbb E\|B_n\|_{\infty\to1}}{(2n)^{3/2}}.
\]

There is a rigorous explicit lower bound obtained by two greedy
updates.  Take odd \(n\), start with \(y^{(0)}=\mathbf1\), and set

\[
x_i=\operatorname{sign}\left(\sum_kB_{ik}\right),
\qquad
y_j=\operatorname{sign}\left(\sum_iB_{ij}x_i\right).
\]

The resulting energy is

\[
x^\top By=\sum_j\left|\sum_iB_{ij}x_i\right|.
\]

For a fixed \(j\), the summands

\[
W_i=B_{ij}\operatorname{sign}\left(\sum_kB_{ik}\right)
\]

are iid over \(i\), and

\[
\mathbb EW_i
=\frac1n\mathbb E|S_n|
=\left(\sqrt{\frac2\pi}+o(1)\right)n^{-1/2}.
\]

The central limit theorem and uniform integrability give

\[
\frac1{\sqrt n}\mathbb E\left|\sum_iW_i\right|
\longrightarrow
\kappa
:=\mathbb E\left|Z+\sqrt{\frac2\pi}\right|
=1.0391966601\ldots ,
\]

where \(Z\) is standard Gaussian.  Therefore

\[
\boxed{
\mathbb E\|B_n\|_{\infty\to1}
\ge(\kappa-o(1))n^{3/2},
}
\tag{7.1}
\]

and

\[
\boxed{
\liminf_n\gamma_{{\rm bip},n}
\ge\frac{\kappa}{2\sqrt2}
=0.3674115027\ldots .
}
\tag{7.2}
\]

Equation (1.3) gives the necessary condition

\[
c\ge\limsup_n\gamma_{{\rm bip},n}
\]

for any same-constant amplification theorem using an independent
random cross block.  The child energy profiles cannot reduce this
threshold, because bulk configurations of polynomial density recover
the full cross optimum at leading order.

In particular, an iid random cross block cannot preserve any
one-copy normalized constant below \(0.3674115\ldots\).  It is already
incompatible with the current universal lower constant

\[
\frac{c_*}{2}=0.3364933644\ldots .
\]

For a randomized flat block \(B\), the analogous necessary threshold
is

\[
\frac{\|B\|_{\infty\to1}}{(n+m)^{3/2}}.
\]

This can be smaller than the iid threshold, so the flat-block route is
not killed.  But the remaining burden is now exact: besides its low
global rectangular norm, the block must satisfy the high-energy
layer-width inequalities in (4.1).  Random switching does not itself
create the required anti-alignment.

## 8. Surviving target and no-go boundary

The correlated-chaining campaign leaves one viable statement:

> Construct flat cross blocks \(B\) with low
> \(\|B\|_{\infty\to1}\) and prove that their restrictions to every
> pair of high child-energy layers obey (4.1), while the unavoidable
> bulk floor remains below the target constant.

For iid random \(B\), (1.2) is the exact reduction and (1.3) is the
decisive no-go boundary.  Generic chaining cannot improve it: the
full cross ground state can be moved into exact bulk layers by changing
only \(O(\sqrt{n\log n})\) spins per side, at an
\(o(N^{3/2})\) operator-norm cost.

For flat/Hadamard blocks, concentration over energy layers is weaker
because the randomness may consist of only \(n+m\) row/column switches.
A proof there needs either a richer flat ensemble or a deterministic
restriction theorem.  The exact object to control is no longer layer
cardinality but the correlated restricted norm

\[
\max_{x\in\mathcal X_A(p),\,y\in\mathcal X_D(q)}
|x^\top By|.
\]
