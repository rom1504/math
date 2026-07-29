# Sparse repair versus block purification

## 1. Normalizations

For a symmetric zero-diagonal signing \(A\) of order \(n\), use the
one-copy energy

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
M(A)=\max_x|H_A(x)|.
\]

Put

\[
d_A(x)=M(A)-|H_A(x)|
\]

and

\[
Z_A(\lambda)=\sum_x e^{-\lambda d_A(x)}.
\]

For the augmented-cut representation, write

\[
h_v=a\cdot v,\qquad
g_v=\frac{M(A)-h_v}{2}.
\]

The two orientations of one projective spin have gaps \(g\) and
\(M-g\), while its absolute gap is

\[
d_A(x)=2\min\{g,M-g\}.
\tag{1.1}
\]

All factors below use this one-copy normalization.

## 2. The two exact randomized criteria

### 2.1 Fixed-cardinality repair

Let

\[
L=\binom n2
\]

and flip a uniformly random \(K\)-subset of the \(L\) edges.  Put
\(p=K/L<1/2\).  For an exact minimizer, no realization lowers
\(M(A)\).  The exact fixed-cardinality repair converse therefore
implies

\[
\boxed{
\sum_v e^{-I_p(g_v)}\ge1,
}
\tag{2.1}
\]

where terms with \(g_v>K\) may be omitted and

\[
I_p(g)=
2pL\left(
\frac{M}{2L}+
\frac{(1-2p)g}{2pL}
\right)^2.
\tag{2.2}
\]

The exact hypergeometric relative-entropy exponent can replace
\(I_p\); (2.2) is the Pinsker form.

### 2.2 Balanced duplicate-row purification

Let \(D\) be an order-\(kn\) signing, with \(k\) a positive integer.
Use every row of \(A\) exactly \(k\) times in the old--new block and
randomize the global column signs and diagonal fills.  The proved
balanced purification criterion has logarithmic prefactor

\[
\boxed{
2Z_D(\lambda)(\cosh\lambda)^{kn}S_{A,k}(\lambda),
}
\tag{2.3}
\]

where

\[
S_{A,k}(\lambda)=
\sum_x e^{-\lambda d_A(x)}
\prod_{j=1}^n
\cosh\!\left(\lambda(Ax)_j\right)^k.
\tag{2.4}
\]

If the logarithm of (2.3) is at most \(\lambda T\), the constructed
block signing has norm at most

\[
M(A)+M(D)+T.
\tag{2.5}
\]

Thus repair and purification depend on the same absolute-gap layers.
Repair tries to prove that a weighted layer sum is *small*.  Exact
optimality says that it is *large*.  Purification also needs a
weighted layer sum to be small.

## 3. Exact variational bridge

Assume \(K\ge M\), as is automatically true when \(p>0\) is fixed and
\(M=O(n^{3/2})\).  Since \(I_p(g)\) is increasing and the two
orientations have gaps \(g,M-g\), (2.1) gives

\[
\sum_{\substack{v\\0\le g_v\le M/2}}
e^{-I_p(g_v)}
\ge\frac12.
\tag{3.1}
\]

For a state represented by the smaller orientation gap \(g\),

\[
|H_A(x)|=M-2g.
\]

Moreover,

\[
\sum_j |(Ax)_j|\ge 2|H_A(x)|.
\]

Convexity of \(\log\cosh\) consequently gives

\[
\prod_j\cosh(\lambda(Ax)_j)^k
\ge
\exp\left\{
kn\log\cosh\left(
\frac{2\lambda(M-2g)}n
\right)
\right\}.
\tag{3.2}
\]

Comparing (3.1) term by term with (2.4) proves the exact bridge

\[
\boxed{
\begin{aligned}
\log S_{A,k}(\lambda)
\ge -\log2+
\inf_{0\le g\le M/2}\Bigg[
&I_p(g)-2\lambda g\\
&+kn\log\cosh\left(
\frac{2\lambda(M-2g)}n
\right)
\Bigg].
\end{aligned}
}
\tag{3.3}
\]

In particular, dropping the nonnegative local-field term gives

\[
\log Z_A(\lambda)
\ge-\log2+
\sup_{0<p<1/2}
\inf_{0\le g\le M/2}
\{I_p(g)-2\lambda g\}.
\tag{3.4}
\]

This makes the direction of the duality explicit: failure of repair
is a **lower bound** on the very partition function which purification
needs to upper-bound.

### 3.1 Closed asymptotic form

Suppose

\[
M(A)=c\,n^{3/2}+o(n^{3/2}),\qquad
\lambda=\frac b{\sqrt n}.
\]

Set \(g=rn^{3/2}\).  Uniformly on compact parameter sets,

\[
\frac1n\{I_p(g)-2\lambda g\}
=
pc^2+
\bigl(2c(1-2p)-2b\bigr)r
+\frac{(1-2p)^2}{p}r^2+o(1).
\tag{3.5}
\]

For \(0<b<2c\), maximize the minimum in (3.5).  The optimizer is

\[
1-2p_*=\frac{2b}{2c+b},
\qquad
r_*=\frac{2c-b}{8}.
\tag{3.6}
\]

Substitution gives

\[
\boxed{
\liminf_{n\to\infty}\frac1n
\log Z_A(b/\sqrt n)
\ge
\psi_c(b):=\frac{(2c-b)_+^2}{8}.
}
\tag{3.7}
\]

At \(b=0\) the statement follows by continuity.  For \(b\ge2c\) the
right side is zero, which is already supplied by a ground state.

Equation (3.7) is exactly the free-energy shadow of fixed-cardinality
repair optimality.  It is not an independent entropy assumption.

## 4. Universal no-overlap for independent purification

The repair floor (3.7) already has the wrong sign for purification.
In fact two signing-independent floors give a complete no-go.

Let

\[
\ell(t)=\mathbb E_{G\sim N(0,1)}\log\cosh(tG).
\]

For a uniform Boolean spin, every row field \((Ax)_j\) is a sum of
\(n-1\) independent signs, for every signing \(A\).  Since
\(d_A(x)\le M(A)\), Jensen's inequality and the central limit theorem
give

\[
\liminf\frac1n\log S_{A,k}(b/\sqrt n)
\ge
\log2-bc+k\ell(b).
\tag{4.1}
\]

A ground state gives the independent floor

\[
\liminf\frac1n\log S_{A,k}(b/\sqrt n)
\ge k\log\cosh(2bc).
\tag{4.2}
\]

If \(D\) has order \(kn\) and

\[
M(D)\le c(kn)^{3/2}+o(n^{3/2}),
\]

then \(0\le d_D(y)\le M(D)\), so

\[
\liminf\frac1n\log Z_D(b/\sqrt n)
\ge
\bigl(k\log2-bc\,k^{3/2}\bigr)_+.
\tag{4.3}
\]

Suppose one tries to preserve the normalized value \(c\).  The
available cross-block allowance is

\[
T=c\delta_k n^{3/2}+o(n^{3/2}),
\qquad
\delta_k=(1+k)^{3/2}-1-k^{3/2}.
\tag{4.4}
\]

Put \(s=bc\) and \(L_0=\log2\).  Since every competitive sequence has
\(c\le1/2\), monotonicity of \(\ell\) gives
\(\ell(b)\ge\ell(2s)\).  Equations (4.1)--(4.4) show that the logarithm
of the purification prefactor minus \(\lambda T\) is at least

\[
\boxed{
\begin{aligned}
R_k(s)={}&(kL_0-sk^{3/2})_+\\
&+\max\{0,L_0-s+k\ell(2s),
k\log\cosh(2s)\}
-s\delta_k.
\end{aligned}
}
\tag{4.5}
\]

This lower bound is strictly positive for every integer \(k\ge1\)
and every \(s\ge0\).

Here is a short verification.  The two terms meeting at

\[
s_0=\frac{L_0}{\sqrt k}
\]

make (4.5) decrease to the cusp and increase after it.  For
\(1\le k\le4\), the ground-state branch suffices and the values at the
cusp are respectively

\[
0.1795498765,\quad
0.1672944948,\quad
0.1511332305,\quad
0.1369259824.
\tag{4.6}
\]

For \(k\ge5\), use

\[
t-t^3\le\ell'(t)\le t,
\qquad
\ell(t)\ge\frac{t^2}{2}-\frac{t^4}{4}.
\tag{4.7}
\]

These inequalities again put the minimum at \(s_0\).  Since

\[
\frac{(k+1)^{3/2}-k^{3/2}}{\sqrt k}
\le\frac32+\frac{3}{8k},
\]

the value there is at least

\[
2L_0^2-\frac{L_0}{2}
-\frac{4L_0^4+3L_0/8}{k}>0
\qquad(k\ge5).
\tag{4.8}
\]

Therefore the balanced duplicate-row union bound cannot certify a
scale-preserving block, for any block ratio or inverse temperature.
The obstruction is exponential in \(n\), not a missing polynomial
prefactor.

### 4.1 Exact logical barrier

Let \(N(g)\) be the population of an absolute-gap layer.  Both

\[
\sum_gN(g)e^{-I_p(g)}
\quad\text{and}\quad
\sum_gN(g)e^{-2\lambda g}
\]

are coordinatewise increasing in every \(N(g)\).  The local-field
factor in (2.4) is also positive.  Thus:

* low layer entropy can make sparse repair succeed;
* if repair fails, exact optimality forces a large weighted layer
  population;
* that same branch makes purification harder, never easier.

Consequently no dichotomy which remembers only the one-replica gap
profile can have a “repair fails, therefore purification succeeds”
branch.  Equations (3.3) and (4.5) quantify the non-overlap.

This is a no-go for the **independent exponential/union-bound method**,
not for all correlated block constructions.

## 5. The natural correlated escape and its exact limit

The most direct correlated replacement is to center a current
internal signing \(h\in\{\pm1\}^m\) by choosing
\(\beta\in\{\pm1\}^m\) uniformly subject to

\[
\langle\beta,h\rangle=0.
\tag{5.1}
\]

Equivalently, choose exactly \(m/2\) coordinates and flip \(h\) there.
This is the \(p=1/2\) endpoint of fixed-cardinality repair.

Fix a feature vector \(\phi\in\{\pm1\}^m\) and put

\[
u=\langle h,\phi\rangle.
\]

If \(K\) is the number of selected coordinates on which
\(h_e\phi_e=+1\), then

\[
K\sim\operatorname{Hyp}
\left(m,\frac{m+u}{2},\frac m2\right),
\qquad
\langle\beta,\phi\rangle=m+u-4K.
\tag{5.2}
\]

Thus the replacement is exactly centered for every feature:

\[
\mathbb E\langle\beta,\phi\rangle=0,
\]

and

\[
\boxed{
\operatorname{Var}\langle\beta,\phi\rangle
=\frac{m^2-u^2}{m-1}.
}
\tag{5.3}
\]

For a family of feature/slack pairs \((\phi_z,r_z)\), an exact
cap-conditioned fill exists whenever

\[
\boxed{
\sum_z
\Pr\left[
\left|m+u_z-4K_z\right|>r_z
\right]<1,
}
\tag{5.4}
\]

where \(K_z\) has the hypergeometric law in (5.2).  This is the
fixed-cardinality correlated analogue of independent internal refill.

### 5.1 The rank-one gain is only subexponential

Put

\[
q=\frac um,\qquad \tau=\frac tm.
\]

For a boundary value
\(\langle\beta,\phi\rangle=t\), Stirling's formula gives the rate

\[
\begin{aligned}
J(q,\tau)=\log2
&-\frac{1+q}{2}
H\left(\frac{1+q-\tau}{2(1+q)}\right)\\
&-\frac{1-q}{2}
H\left(\frac{1-q+\tau}{2(1-q)}\right),
\end{aligned}
\tag{5.5}
\]

where \(H\) is binary entropy.  The two tails have the same rate.
Uniformly for \(q,\tau=O(n^{-1/2})\),

\[
\boxed{
mJ(q,\tau)
=\frac{t^2}{2m}+O(1),
}
\tag{5.6}
\]

and the logarithm of the tail differs from \(-mJ\) by
\(O(\log m)\).

In the capped-cut application,

\[
m=\Theta(n^2),\qquad
u=O(n^{3/2}),\qquad
t=O(n^{3/2}).
\]

Therefore (5.6) is the same speed-\(n\) exponent as independent sign
refill, up to \(o(n)\).  Conditioning out the current kernel direction
does not change the leading traffic rate.

The covariance makes the reason transparent.  The law (5.1) has

\[
\operatorname{Cov}(\beta)
=\frac m{m-1}
\left(I-\frac{hh^\top}{m}\right).
\tag{5.7}
\]

It removes only the rank-one variance

\[
\frac{\langle h,\phi\rangle^2}{m}=O(n)
\]

from a total feature variance \(m=\Theta(n^2)\).  At a threshold
\(t=\Theta(n^{3/2})\), this changes the logarithmic tail exponent by
only \(O(1)\), while the known traffic deficit is \(\Theta(n)\).

More generally, if every available kernel direction \(h_r\) satisfies

\[
|\langle h_r,\phi\rangle|=O(n^{3/2})
\]

on the cap family, then each normalized direction can remove only
\(O(n)\) variance.  A constant-factor reduction of the
\(\Theta(n^2)\) feature variance therefore needs

\[
\boxed{\Omega(n)\text{ essentially independent kernel directions}.}
\tag{5.8}
\]

The single conditional kernel supplied by the current
feature-covariance theorem is insufficient, even when imposed by an
exact fixed-cardinality constraint.

## 6. Verdict and surviving target

The attempted repair--purification synthesis closes as a no-go:

1. the exact repair converse lower-bounds the same gap partition
   function used by purification;
2. universal state-entropy and local-field floors make every
   independent duplicate-row union bound miss the scale-preserving
   allowance by an exponential factor;
3. the natural correlated half-flip removes only one covariance
   direction and has the same leading large-deviation rate as
   independent refill.

The viable replacement must therefore be genuinely chained.  It has
to turn cap geometry into \(\Theta(n)\) low-variance feature
directions, or directly control the supremum without a union bound.
A concrete next lemma would be:

> From the conditional feature-kernel identity on every cap
> restriction, construct \(\Theta(n)\) approximately orthogonal sign
> directions whose span captures a fixed fraction of every thick-cap
> feature vector.

If such a lemma held, a multi-constraint slice replacement could
improve the exponent in (5.6) at speed \(n\).  Without a rank-growth
statement of this kind, sparse repair and block purification do not
overlap in a way that can prove scale transfer.

## 7. Audit of the multi-block dual-measure escape

For an endpoint pair \(p=(x,y)\), define its width gradient and
midpoint feature by

\[
s_e(p)=\frac{x_ix_j-y_iy_j}{2},
\qquad
r_e(p)=\frac{x_ix_j+y_iy_j}{2}.
\tag{7.1}
\]

Coordinatewise,

\[
s_e(p)r_e(p)=0,\qquad s_e(p)^2+r_e(p)^2=1.
\tag{7.2}
\]

The global block minimax theorem supplies, separately for each edge
block \(T\), a law \(\mu_T\) with

\[
\mathbb E_{\mu_T}\operatorname{score}_A
\ge W(A)-\eta_T,
\qquad
2\sum_{e\in T}
\bigl(a_e\mathbb E_{\mu_T}s_e\bigr)_+
\le\eta_T,
\tag{7.3}
\]

where

\[
\eta_T=\sqrt{2|T|(2v(T)\log2+1)}.
\]

This does **not** yet produce many midpoint-kernel directions.  There
are two exact losses.

First, (7.3) is a first-moment statement about \(s_T\), while the
desired directions are the complementary vectors \(r(p)\).  This
separation is graph-realizable.  If \(T\) is bipartite with shores
\(U,V\), take endpoint pairs satisfying

\[
y_i=-x_i\ (i\in U),\qquad y_i=x_i\ (i\in V).
\]

Then on every edge of \(T\),

\[
r_e(p)=0,\qquad s_e(p)=x_ix_j.
\tag{7.4}
\]

A sign-symmetric distribution of \(x\)'s can cancel all first moments
of \(s_T\), while producing no midpoint direction at all on the
controlled coordinates.  Near-activity may rule out a particular
such law for a particular minimizer, but that conclusion uses
information absent from (7.3).  Thus anti-aligned gradients alone do
not imply midpoint-rank growth.

Second, the laws in (7.3) depend on \(T\).  Combining blocks by
replacing them with their union pays the rounding error for the
union.  For a \(d\)-regular graphing on all \(n\) vertices,

\[
|T|=\frac{dn}{2},\qquad v(T)=\Theta(n),
\]

and hence

\[
\eta_T=\Theta(n\sqrt d).
\tag{7.5}
\]

If \(d=o(n)\), this gives a genuinely near-active common law but
controls only an \(o(1)\) fraction of the \(\Theta(n^2)\) edge
features.  If \(d=\Theta(n)\), it controls a constant fraction, but
\(\eta_T=\Theta(n^{3/2})\), the full width scale.  At that point the
near-activity statement can be quantitatively vacuous.

Equivalently, one may obtain strong laws on many sparse graphings, but
they can be mutually singular.  The present minimax theorem has no
consistency clause that places their midpoint directions in one
covariance matrix.  Passing to one common law recreates the critical
\(\sqrt{|T|v(T)}\) loss in (7.5).

Therefore the current global block duality plus the one-direction
tangent kernel does not prove the \(\Omega(n)\)-direction statement
in (5.8).  A sufficient new statement would need one of:

1. a common near-active law valid on \(\Theta(n)\) overlapping sparse
   graphings;
2. a rank-sensitive separation theorem showing that low-dimensional
   midpoint span lets one choose a single block \(T\) violating
   (7.3);
3. a direct entropy/rank lower bound on midpoint features inside one
   near-active dual law.

None follows from the current first-moment inequalities, and (7.4)
is the precise local obstruction.
