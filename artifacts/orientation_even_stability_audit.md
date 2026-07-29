# Orientation-even stability: audited theorem and scale obstruction

Let \(A\) be a symmetric zero-diagonal \(n\times n\) signing, put

\[
m=n-1,\qquad
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
M(A)=\max_x|H_A(x)|.
\]

## 1. The finite arcsine remainder is valid, with a stronger constant

For \(s=\pm1\), let

\[
X^{(s)}
=\operatorname{sign}\bigl((sA+\sqrt m\,I)g\bigr),
\qquad g\sim N(0,I).
\]

Every pre-sign coordinate has variance \(2m\).  For \(i\ne j\), after
gauging by \(a_{ij}\), its correlation is

\[
s\,a+b_{ij},\qquad
a=\frac1{\sqrt m},\qquad
b_{ij}=\frac{a_{ij}(A^2)_{ij}}{2m}.
\]

Since

\[
(A^2)_{ij}
=\sum_{k\ne i,j}a_{ik}a_{kj},
\]

we have

\[
|b_{ij}|\le\frac{m-1}{2m}.
\]

Moreover

\[
\frac1{\sqrt m}+\frac{m-1}{2m}\le1,
\]

so \(a\pm b_{ij}\in[-1,1]\) throughout the full admissible domain.

Put

\[
F_a(b)
=\arcsin(a+b)+\arcsin(a-b)-2\arcsin a.
\]

For

\[
h(u)=\frac{u}{(1-u^2)^{3/2}},
\]

\[
F_a''(b)=h(a+b)+h(a-b).
\]

Now

\[
h'(u)=\frac{1+2u^2}{(1-u^2)^{5/2}}
\]

is even and increasing as a function of \(|u|\).  Therefore, for
\(b\ge0\),

\[
h(a+b)-h(a)
=\int_0^b h'(a+t)\,dt
\ge
\int_0^b h'(a-t)\,dt
=h(a)-h(a-b).
\]

Thus

\[
F_a''(b)\ge2h(a)
=\frac{2a}{(1-a^2)^{3/2}}.
\]

Since \(F_a(0)=F_a'(0)=0\), evenness gives

\[
\boxed{
F_a(b)\ge
\frac{a}{(1-a^2)^{3/2}}\,b^2.
}
\tag{1}
\]

The opposite-orientation Gaussian argument then yields, for \(m>1\),

\[
\boxed{
M(A)\ge
\frac{nm}{\pi}\arcsin\frac1{\sqrt m}
+
\frac{\|A^2-mI\|_F^2}
{8\pi m^{5/2}(1-1/m)^{3/2}}.
}
\tag{2}
\]

The previously recorded denominator \(8\pi m^{5/2}\) is therefore
valid but slightly nonsharp.  The \(n=2\) case has \(A^2=I\) and is
separate.

## 2. Why the correction cannot force conference structure

The coefficient in (2) is on the scale \(n^{-5/2}\).  A Wigner-like
signing has

\[
\|A^2-mI\|_F^2=\Theta(n^3),
\]

so (2) distinguishes it from a conference matrix by only
\(\Theta(\sqrt n)\), not by a leading \(\Theta(n^{3/2})\) term.

This obstruction persists inside the exact self-complementary class.
Here is a concrete ensemble.  Let \(n=4r\), index vertices by
\((a,u)\in\mathbb Z_4\times[r]\), and choose independently:

- a symmetric zero-diagonal signing \(C_0\);
- an arbitrary \(r\times r\) signing \(C_1\);
- a symmetric \(r\times r\) signing \(C_2\), including its diagonal.

Put

\[
C_3=-C_1^\top
\]

and define the \(r\times r\) block from layer \(a\) to layer \(a+d\) by

\[
A_{a,a+d}=(-1)^a C_d.
\]

This is symmetric, has zero diagonal and signs everywhere else.  If
\(P(a,u)=(a+1,u)\), then

\[
\boxed{PAP^\top=-A.}
\]

A direct orbit-monomial count gives

\[
\boxed{
\mathbb E\|A^2-(n-1)I\|_F^2
=n^3-2n^2-4n.
}
\tag{3}
\]

Indeed, after expanding each off-diagonal entry of \(A^2\), independent
orbit signs kill unequal monomials.  Fix \(i=(a,u)\), put
\(j=(a+d,v)\), and write \(d\) modulo \(4\).  The surviving second
moments are

\[
\begin{array}{c|cc}
d &u\ne v&u=v\\ \hline
0&4r-2&\text{diagonal}\\
1,3&4r-2&4r-4\\
2&4r-2&8r-4 .
\end{array}
\]

Summing these over \(j\ne i\), and then over the \(4r\) choices of
\(i\), gives

\[
64r^3-32r^2-16r.
\]

Consequently some self-complementary members satisfy

\[
\|A^2-mI\|_F^2=\Theta(n^3).
\]

To justify that the expectation is not carried by rare high-norm
members, use the standard subgaussian norm bounds for the two symmetric
Wigner blocks \(C_0,C_2\) and the iid block \(C_1\).  With probability
\(1-e^{-\Omega(r)}\), all three have operator norm \(O(\sqrt r)\), and
the \(4\times4\) block form then gives
\(\|A\|_{\rm op}=O(\sqrt n)\).  On that event the defect is \(O(n^3)\).
Off the event it is at most \(O(n^4)\), whose contribution to the
expectation is \(o(n^3)\).  Hence the conditional expectation on the
good event is \((1-o(1))n^3\), giving a realization with both
\(\|A\|_{\rm op}=O(\sqrt n)\) and defect \(\Theta(n^3)\).

For such a member,

\[
\operatorname{tr}A^4
=nm^2+\|A^2-mI\|_F^2
\ge(2-o(1))n^3.
\]

Since

\[
\operatorname{tr}A^4
\le \|A\|_{\rm op}^2\operatorname{tr}A^2,
\]

it also has

\[
\|A\|_{\rm op}\ge(\sqrt2-o(1))\sqrt n.
\]

Nevertheless (2) adds only \(O(\sqrt n)\).  Hence the finite arcsine
stability theorem, even with its optimal quadratic Taylor coefficient,
cannot imply

\[
\|A\|_{\rm op}=(1+o(1))\sqrt n.
\]

## 3. A broader pair-response ceiling

The scale loss is not peculiar to the sign/arcsine kernel.  If a
coordinatewise Gaussian Boolean rule has noise-stability kernel

\[
K(\rho)=\sum_{k\ge0}w_k\rho^k,
\qquad w_k\ge0,\quad \sum_k w_k=1,
\]

then, for \(|\rho|\le1/2\),

\[
|K''(\rho)|
\le
\sup_{k\ge2}k(k-1)|\rho|^{k-2}
\le3.
\tag{4}
\]

If the rule is odd, only odd \(k\) occur and

\[
|K''(\rho)|\le6|\rho|.
\tag{5}
\]

In all row-field Gaussian constructions the baseline edge correlation
is \(O(n^{-1/2})\), while the row-correlation perturbation is
\((A^2)_{ij}/O(n)\).  Therefore:

- a general coordinatewise pair kernel can contribute at most
  \(O(\|A^2-mI\|_F^2/n^2)\);
- an odd kernel contributes at most
  \(O(\|A^2-mI\|_F^2/n^{5/2})\).

On the self-complementary Wigner-scale ensemble (3), these are at most
\(O(n)\) and \(O(\sqrt n)\), respectively.  Neither changes the leading
\(n^{3/2}\) constant.  A successful orientation-even theorem must use
a genuinely global statistic, an \(n\)-dependent singular response, or
dependent multichannel selection; a fixed smooth pair response cannot
do it.

## 4. Exact spectral-tail localization lemmas

Although the orientation-even correction does not control the top
eigenvalue, two exact reductions isolate the remaining anomalous
branch.

If \(Av=\lambda v\), \(\|v\|_2=1\), independently choose Boolean \(X_i\)
with

\[
\mathbb E X_i=\frac{v_i}{\|v\|_\infty}.
\]

The zero diagonal gives

\[
\mathbb E X^\top AX
=\frac{\lambda}{\|v\|_\infty^2}.
\]

Therefore

\[
\boxed{
Q(A):=\max_x|x^\top Ax|
\ge\frac{|\lambda|}{\|v\|_\infty^2}.
}
\tag{6}
\]

Thus any eigenvalue \(|\lambda|\gg\sqrt n\) in a sequence with
\(Q(A)=O(n^{3/2})\) has a localized eigenvector:

\[
\frac1{\|v\|_\infty^2}
\le\frac{Q(A)}{|\lambda|}
=o(n).
\]

There is also a tail-dimension bound.  Let \(E_\lambda\) be the span of
the eigenvectors with eigenvalues at least \(\lambda>0\), and put
\(k=\dim E_\lambda\).  An extreme point of

\[
E_\lambda\cap[-1,1]^n
\]

has at least \(k\) active coordinate constraints, hence squared
Euclidean norm at least \(k\).  Since a zero-diagonal quadratic form
can be rounded coordinate by coordinate from the cube to a vertex
without decreasing it,

\[
\max_xx^\top Ax\ge\lambda k.
\]

Applying the same statement to \(-A\) gives

\[
\boxed{
\#\{i:\lambda_i(A)\ge\lambda\}
\le\frac{Q(A)}{\lambda},
\qquad
\#\{i:\lambda_i(A)\le-\lambda\}
\le\frac{Q(A)}{\lambda}.
}
\tag{7}
\]

Hence the super-\(\sqrt n\) spectral branch is necessarily low
dimensional, and its top eigenvectors are localized.  This is the
branch a vertex-peeling argument must remove.  Broad Wigner-like
constant-factor spectral spread lies in the \(O(\sqrt n)\) regular
branch but cannot be distinguished at leading order by (2).

## 5. A leading-scale orientation-even \(A^2\)-energy theorem

The lower-order obstruction above applies to the **mean quadratic
energy**, but not to the squared local-field energy of the same
witnesses.  Put

\[
X^\sigma=\operatorname{sign}\bigl((\sigma A+\sqrt m I)g\bigr).
\]

For \(i\ne j\), abbreviate

\[
c=(A^2)_{ij},\qquad
u=\frac{c}{2m},\qquad
v=\frac{a_{ij}}{\sqrt m}.
\]

The correlation averaged over the two orientations is

\[
C_{ij}
=\frac1\pi\left[
\arcsin(u+v)+\arcsin(u-v)
\right].
\]

Let

\[
G_v(u)=\arcsin(u+v)+\arcsin(u-v).
\]

Then \(G_v(0)=0\), while

\[
G_v'(u)
=g(u+v)+g(u-v),
\qquad
g(z)=\frac1{\sqrt{1-z^2}}.
\]

The function \(g\) is even and convex.  Hence the even convex function
\(G_v'\) has its minimum at zero:

\[
G_v'(u)\ge2g(v)
=\frac2{\sqrt{1-1/m}}.
\]

It follows that

\[
cC_{ij}
\ge
\frac{c^2}{\pi m\sqrt{1-1/m}}.
\]

Expanding \(X^\top A^2X\) and summing gives the exact finite theorem

\[
\boxed{
\frac12\sum_{\sigma=\pm1}
\mathbb E (X^\sigma)^\top A^2X^\sigma
\ge
nm+
\frac{\|A^2-mI\|_F^2}
{\pi m\sqrt{1-1/m}}.
}
\tag{8}
\]

Unlike (2), this detects a defect \(\Theta(n^3)\) by a leading
\(\Theta(n^2)\) increase in squared local-field energy.

## 6. Exact capped-field conversion

There is an exact way to combine oriented energy and squared fields
without first selecting one sample having both.

For any Boolean witness \(x\), orient \(A\) so

\[
q=x^\top Ax,\qquad
r_i=x_i(Ax)_i,\qquad
S_K=\sum_i r_i\,\operatorname{clip}
\bigl(r_i,[-K\sqrt m,K\sqrt m]\bigr).
\]

In the switched gauge put

\[
u_i=\operatorname{clip}
\bigl(r_i,[-K\sqrt m,K\sqrt m]\bigr),
\qquad
\mu=(1-\alpha)\mathbf1+
\frac{\alpha u}{K\sqrt m},
\quad 0\le\alpha\le1.
\]

The vector \(\mu\) lies in the cube.  Since \(A\) has zero diagonal,
coordinatewise randomized rounding shows that its quadratic value is
at most \(Q(A)\).  Also

\[
\mathbf1^\top Du=r^\top u=S_K,
\qquad
|u^\top Du|\le K^2mQ(A).
\]

Expanding \(\mu^\top D\mu\) therefore proves

\[
\boxed{
(1+\alpha^2)Q(A)
\ge
(1-\alpha)^2q+
\frac{2\alpha(1-\alpha)}{K\sqrt m}S_K.
}
\tag{9}
\]

This may be averaged over any witness distribution.  Set

\[
C=\frac{Q(A)}{n\sqrt m},\quad
c=\frac{\mathbb E q}{n\sqrt m},\quad
z=\frac{\mathbb E S_K}{Knm}.
\]

Optimizing (9) in \(\alpha\) gives

\[
\boxed{
C\ge {\cal F}(c,z),
\qquad
{\cal F}(c,z)=
\begin{cases}
c,&z\le c,\\[2mm]
c-z+\sqrt{z^2+(z-c)^2},&z>c.
\end{cases}
}
\tag{10}
\]

For the Gaussian witnesses in (8),

\[
c\ge
\frac{2\sqrt m}{\pi}
\arcsin\frac1{\sqrt m}
=\frac2\pi+o(1).
\]

Write

\[
\delta_n=
\frac{\|A^2-mI\|_F^2}{nm^2}
\]

and define the normalized field-square tail

\[
\Psi_n(K)
=\frac1{nm}\,
\mathbb E\sum_i r_i^2
\mathbf1_{\{|r_i|>K\sqrt m\}}.
\]

Since \(S-S_K\) is at most the displayed tail square,
(8)--(10) yield

\[
\boxed{
C\ge
{\cal F}\left(
\frac2\pi+o(1),
\frac{
1+\delta_n/\bigl(\pi\sqrt{1-1/m}\bigr)-\Psi_n(K)}
K
\right).
}
\tag{11}
\]

For reference,

\[
{\cal F}(2/\pi,z)>0.672986728863\ldots
\quad\Longleftrightarrow\quad
z>0.8942308433\ldots .
\]

Thus (11) is a quantitative defect-versus-\(Q\) theorem whenever a
fixed cap retains enough of the squared local fields.  The sole loss is
the heavy-field tail \(\Psi_n(K)\).  The universal-vertex construction
shows that this tail cannot be bounded from the first two moments or
cut positivity alone.

The negative-field flip correction gives the complementary exact fact

\[
L_-^2\le Q(A)(Q(A)-q),
\]

so a witness whose oriented energy nearly saturates \(Q(A)\) cannot
hide (8) in heavy **negative** fields.  For the Gaussian witnesses the
baseline is only \(2/\pi\), leaving a constant energy gap to the current
\(0.672986\ldots\) bound; consequently (11), rather than this
near-saturation observation, is the effective unconditional
combination.

## 7. Exact finite defect gain for smoothed field-plus-spin

Fix \(t\in\mathbb R\) and \(\tau>0\), and let

\[
X_i^\sigma
=\operatorname{sign}\left(
\sigma\frac{(A\xi)_i}{\sqrt m}
+t\xi_i+\tau Z_i
\right),
\qquad \sigma\in\{\pm1\},
\]

where the \(\xi_i\)'s are independent Rademachers and the \(Z_i\)'s
are independent standard Gaussians.  Write

\[
\psi_\tau(s)=2\Phi(s/\tau)-1.
\]

There is first a completely finite, though numerically crude, defect
theorem for these exact Rademacher-field witnesses.  Put

\[
R_t=2\sqrt{1+t^2},\qquad
d_{t,\tau}
=\min_{|z|\le R_t+1}\psi_\tau'(z)>0,
\qquad
\kappa_{t,\tau}=\frac{d_{t,\tau}^2}{2}.
\]

Then

\[
\boxed{
\frac12\sum_{\sigma=\pm1}
\mathbb E (X^\sigma)^\top A^2X^\sigma
\ge
nm+\kappa_{t,\tau}\frac{\|A^2-mI\|_F^2}{m}.
}
\tag{12}
\]

Here is the exact discrete proof.  Fix \(i\ne j\), switch so that
\(a_{ij}=1\), and change common Rademachers so that the row-\(i\)
coefficients are all \(+1\).  If the row-\(j\) common coefficients
are \(c_k\in\{\pm1\}\), put

\[
d=\sum_{k\ne i,j}c_k=(A^2)_{ij}.
\]

By exchangeability, the orientation-averaged pair correlation

\[
C_m(d)=\frac12\sum_{\sigma=\pm1}
\mathbb E X_i^\sigma X_j^\sigma
\]

depends only on \(d\), and a measure-preserving sign change gives
\(C_m(-d)=-C_m(d)\).

Explicitly, in this gauge the two pre-threshold fields are

\[
\sigma m^{-1/2}\left(\xi_j+\sum_k\eta_k\right)+t\xi_i+\tau Z_i,
\qquad
\sigma m^{-1/2}\left(\xi_i+\sum_kc_k\eta_k\right)+t\xi_j+\tau Z_j.
\]

In the law with coefficients \(-c_k\), make the substitution
\[
(\eta_k,\sigma,\xi_j,Z_j)
\longmapsto(-\eta_k,-\sigma,-\xi_j,-Z_j).
\]
The first field is unchanged and the second is negated.  This proves
the claimed oddness without an approximation.

Change one common coefficient from \(-1\) to \(+1\), so that
\(d\) increases by \(2\), and condition on all variables except the
corresponding common Rademacher.  If \(s=m^{-1/2}\), the exact increase
is

\[
\frac14\sum_{\sigma=\pm1}
\mathbb E\left[
\Delta_s(W_i^\sigma)\Delta_s(W_j^\sigma)
\right],
\qquad
\Delta_s(w)=\psi_\tau(w+s)-\psi_\tau(w-s).
\tag{13}
\]

The residual fields satisfy
\(\mathbb E(W_i^\sigma)^2,\mathbb E(W_j^\sigma)^2\le1+t^2\).
Thus, by a union bound and Chebyshev,

\[
\mathbb P\{|W_i^\sigma|\le R_t,\ |W_j^\sigma|\le R_t\}\ge\frac12.
\]

On this event both increments in (13) are at least
\(2d_{t,\tau}/\sqrt m\).  Consequently each one-coordinate increase
of \(C_m\) is at least \(d_{t,\tau}^2/m\).  Anchoring at \(d=0\) when
the parity is even, or using
\(2C_m(1)=C_m(1)-C_m(-1)\) when it is odd, gives

\[
d\,C_m(d)\ge\frac{d_{t,\tau}^2}{2m}d^2.
\]

Summing over ordered pairs proves (12).  Explicitly,

\[
\kappa_{t,\tau}
=\frac1{\pi\tau^2}
\exp\left(
-\frac{(1+2\sqrt{1+t^2})^2}{\tau^2}
\right),
\]

so this finite coefficient degenerates badly as \(\tau\downarrow0\).

## 8. Relative Lindeberg refinement

**Status.**  The derivation below has passed an independent proof audit,
including the derivative-product replacement, all orientation factors,
both parity anchors, and exact small-\(m\) enumeration.

The preceding proof has considerably more quantitative content if the
smooth replacement is applied to the **factored increment** (13),
rather than to the original pair correlation.

Define

\[
f_{t,\tau}(u)
=\frac{\psi_\tau(u+t)+\psi_\tau(u-t)}2,
\qquad
K_{t,\tau}(q)
=\mathbb E f_{t,\tau}(G)f_{t,\tau}(H),
\]

where \(G,H\) are standard Gaussians with correlation \(q\).
For every fixed \((t,\tau)\), the following uniform estimate holds:

\[
\boxed{
\left|C_m(d)-K_{t,\tau}(d/m)\right|
\le
C_{t,\tau}\frac{|d|}{m^{3/2}}.
}
\tag{14}
\]

The constant is independent of the coefficient pattern, of \(d\),
and of \(m\).

For completeness, the normalization behind (14) is as follows.  In a
step from \(d=r-1\) to \(d=r+1\), remove the flipped common
Rademacher.  The remaining common field is a sum of \(m-2\)
independent two-dimensional vectors, each of Euclidean size
\(O(m^{-1/2})\).  In (13), \(\Delta_s\), together with each of its
first three derivatives, is \(O_{t,\tau}(s)\).  Thus every derivative
of total order at most three of the product test function is
\(O_{t,\tau}(s^2)\).  The two-dimensional Lindeberg replacement error
is consequently

\[
O_{t,\tau}\left(
(m-2)m^{-3/2}s^2
\right)
=O_{t,\tau}(m^{-3/2}).
\]

This estimate is uniform even when the limiting covariance is
singular, because it uses only bounded derivatives of the smoothed
test function.

After replacement, average over the two endpoint spins and the two
orientations.  Taylor expansion of each factored increment gives

\[
C_m(r+1)-C_m(r-1)
=\frac2m K_{t,\tau}'(r/m)
+O_{t,\tau}(m^{-3/2}).
\tag{15}
\]

The direct endpoint terms are \(O(m^{-1/2})\) inside the smooth
arguments and hence contribute only \(O(m^{-3/2})\) to (15).
The Gaussian common fields have variance \(1-2/m\); changing this to
one has the same smaller-order cost.  Finally telescope (15), using
the parity anchor \(C_m(0)=0\), or
\(2C_m(1)=C_m(1)-C_m(-1)\).  The midpoint quadrature error is
\(O_{t,\tau}(|d|/m^2)\), while the accumulated replacement error is
\(O_{t,\tau}(|d|/m^{3/2})\).  This proves (14), including when
\(|d|=\Theta(m)\).

The function \(f_{t,\tau}\) is odd.  Its Gaussian noise stability has
the Hermite expansion

\[
K_{t,\tau}(q)
=\sum_{k\ {\rm odd}}w_kq^k,
\qquad w_k\ge0.
\]

The first coefficient is

\[
w_1
=\left(\mathbb E f_{t,\tau}'(G)\right)^2
=4\phi_{1+\tau^2}(t)^2.
\]

Therefore

\[
qK_{t,\tau}(q)
\ge
4\phi_{1+\tau^2}(t)^2q^2.
\tag{16}
\]

Combining (14) and (16), pair by pair, gives

\[
d\,C_m(d)
\ge
\left(
w_1-O_{t,\tau}(m^{-1/2})
\right)\frac{d^2}{m}.
\]

Thus the Gaussian coefficient transfers to the original Rademacher
row fields:

\[
\boxed{
\frac12\sum_{\sigma=\pm1}
\mathbb E (X^\sigma)^\top A^2X^\sigma
\ge
nm+
\left(
4\phi_{1+\tau^2}(t)^2-O_{t,\tau}(m^{-1/2})
\right)
\frac{\|A^2-mI\|_F^2}{m}.
}
\tag{17}
\]

At the optimal zero-dither threshold \(t_*=0.8769009856\ldots\), the
coefficient is

\[
4\phi(t_*)^2=0.2950713629\ldots .
\]

Hence, taking \(n\to\infty\) first and then \(\tau\downarrow0\), the
exact field-plus-spin witness has

\[
\frac{\mathbb ES}{nm}
\ge1+
0.2950713629\ldots\,\delta_n.
\tag{18}
\]

The order of limits is necessary because the constants controlling
the smoothed derivatives diverge as \(\tau\downarrow0\).

## 9. The resulting scalar inequality and its exact bottleneck

The energy theorem for these same witnesses gives, along every
competing sequence,

\[
\frac{\mathbb E R}{n\sqrt m}
\ge e(t,\tau)-o(1),
\]

where \(R=\sigma(X^\sigma)^\top AX^\sigma\) and

\[
e(t,\tau)
=4\phi_{1+\tau^2}(t)
\left[
2\Phi\left(\frac{t}{\sqrt{1+\tau^2}}\right)-1
\right].
\]

Put

\[
w(t,\tau)=4\phi_{1+\tau^2}(t)^2,\qquad
C_n=\frac{Q(A)}{n\sqrt m},\qquad
\delta_n=\frac{\|A^2-mI\|_F^2}{nm^2}.
\]

For a field cap \(H>0\), define the tail of these same witnesses by

\[
\Psi_{t,\tau,n}(H)
=\frac1{nm}\,
\mathbb E\sum_i r_i^2
\mathbf1_{\{|r_i|>H\sqrt m\}},
\qquad
r_i=\sigma X_i^\sigma(AX^\sigma)_i.
\]

Substitution of (17) into the exact capped conversion (9)--(10)
produces the scalar inequality

\[
\boxed{
C_n\ge
{\cal F}\left(
e(t,\tau)-o(1),
\frac{
1+\bigl(w(t,\tau)-o(1)\bigr)\delta_n
-\Psi_{t,\tau,n}(H)}
H
\right).
}
\tag{19}
\]

This is the strongest conclusion available from energy plus
orientation-even defect without a separate tail theorem.  In
particular, if the rounded local fields are uniformly capped by
\(H\sqrt m\), then \(\Psi=0\).

Let

\[
c_*=\max_t e(t,0)=0.672986728863\ldots .
\]

For a target \(C_0\ge e\), inversion of
\({\cal F}(e,z)=C_0\) gives

\[
z_{\rm req}(e,C_0)
=C_0+\sqrt{2C_0(C_0-e)}.
\tag{20}
\]

Consequently (19) strictly improves \(c_*\) whenever, for some
\((t,\tau,H)\),

\[
1+w(t,\tau)\delta
-\Psi_{t,\tau}(H)
>
H\left[
c_*+\sqrt{2c_*(c_*-e(t,\tau))}
\right].
\tag{21}
\]

At the energy optimizer \(t=t_*\), after \(\tau\downarrow0\), this
reduces to

\[
\boxed{
\Psi_{t_*,0}(H)
<
1+0.2950713629\ldots\,\delta
-0.672986728863\ldots\,H.
}
\tag{22}
\]

With no tail, saturation at \(c_*\) therefore requires

\[
H\ge
\frac{1+0.2950713629\ldots\,\delta}
{0.672986728863\ldots}.
\tag{23}
\]

For \(\delta=1\), the right side is

\[
1.92436\ldots .
\]

Numerical optimization of (21) over \(t\) confirms that the threshold
in (23) is attained at the cusp \(t=t_*\): moving away from the energy
optimizer loses more baseline energy than it gains in the first
Hermite coefficient.  For example, at
\(H=1.941916296\ldots\), a tail-free defect must satisfy

\[
\delta>1.04004\ldots
\]

before (19) improves \(c_*\).

The negative-field correction

\[
L_-^2\le Q(A)(Q(A)-R)
\]

shows that samples whose oriented energy nearly equals \(Q(A)\) cannot
place the tail in negative fields.  It does **not** control sparse
heavy positive fields.  A single field of order \(n\) can carry order
\(n^2\) squared-field energy while contributing only order \(n\) to
the oriented energy.  This is precisely the universal-vertex
obstruction.  Therefore the exact missing input for (19) is a
positive-heavy-field tail theorem, or a peeling theorem which deletes
the vertices carrying that tail and preserves the normalized
min--max value.

## 10. Conditional-independence and \(2\to4\) tail reduction

The dither gives a rigorous diffuse-versus-localized formulation of
that missing theorem.  Condition on \((\xi,\sigma)\), and put

\[
\mu_j
=\mathbb E_Z[X_j^\sigma\mid\xi,\sigma]
=\psi_\tau\left(
\sigma\frac{(A\xi)_j}{\sqrt m}+t\xi_j
\right),
\qquad
b=A\mu.
\]

Conditional on \((\xi,\sigma)\), the coordinates \(X_j^\sigma\) are
independent.  In particular, since the diagonal of \(A\) is zero,

\[
(AX^\sigma)_i=b_i+\eta_i,
\]

where \(\eta_i\) is a centered sum of \(m\) independent variables of
range length \(2\).  Suppose additionally that

\[
\|A\|_{\rm op}\le K\sqrt m.
\]

Hoeffding's inequality and the split

\[
\{|b_i+\eta_i|>H\sqrt m\}
\subset
\{|b_i|>H\sqrt m/2\}
\cup
\{|\eta_i|>H\sqrt m/2\}
\]

give, for \(H\ge1\),

\[
\boxed{
\Psi_{t,\tau,n}(H)
\le
\frac{20}{H^2nm^2}\,
\mathbb E_\xi\|A\mu(\xi)\|_4^4
+
\left(4K^2+H^2+8\right)e^{-H^2/8}.
}
\tag{24}
\]

Thus the only non-product obstruction is the normalized \(2\to4\)
mass

\[
\Theta_4
=\frac{\mathbb E_\xi\|A\mu(\xi)\|_4^4}{nm^2}.
\]

For fixed \((K,t,\tau)\), this quantity is bounded uniformly in \(n\).
Here is a direct influence proof.  Let

\[
L_1=\|\psi_\tau'\|_\infty,\qquad
L_2=\|\psi_\tau''\|_\infty.
\]

For fixed \(i\), write \(b_i(\xi)=(A\mu(\xi))_i\), and let
\(\Delta_kb_i\) be its change when \(\xi_k\) is flipped.  Taylor
expansion in the \(m-1\) small off-diagonal influences gives

\[
\Delta_kb_i
=-\frac{2\sigma\xi_k}{\sqrt m}(ADA)_{ik}+e_{ik},
\qquad
D=\operatorname{diag}\psi_\tau'\left(
\sigma A\xi/\sqrt m+t\xi
\right),
\]

with

\[
|e_{ik}|\le2+2L_2.
\]

The term \(2\) is the single self-spin response \(j=k\); every other
Taylor remainder is at most \(2L_2/m\).  Therefore, pointwise in
\(\xi\),

\[
\begin{aligned}
\sum_k(\Delta_kb_i)^2
&\le
\frac8m\|e_i^\top ADA\|_2^2
+2n(2+2L_2)^2\\
&\le
\left[
8K^2L_1^2+2(2+2L_2)^2
\right]n
=:\Gamma_{K,t,\tau}n.
\end{aligned}
\tag{25}
\]

Also \(b_i(-\xi)=-b_i(\xi)\), so \(\mathbb E b_i=0\).
The hypercube log-Sobolev inequality applied to the pointwise
carré-du-champ bound (25) implies

\[
\mathbb E b_i^4
\le C\Gamma_{K,t,\tau}^2n^2.
\]

After summing over \(i\),

\[
\boxed{
\Theta_4\le C\Gamma_{K,t,\tau}^2.
}
\tag{26}
\]

Equations (24)--(26) prove uniform integrability of the local-field
squares on every fixed spectrally regular, fixed-dither branch.
They also show why this fact alone does not satisfy (21).  Explicitly,

\[
L_1=\sqrt{\frac2\pi}\,\tau^{-1},
\qquad
L_2=\sqrt{\frac2{\pi e}}\,\tau^{-2},
\]

so the displayed bound deteriorates like a power of \(\tau^{-1}\).
Moreover the product-noise term in (24) is already much too large at
the decisive cap \(H\approx1.9\).  The theorem is therefore a genuine
tail reduction and a peeling/localization criterion, but not yet the
sharp fixed-cap estimate required to beat \(c_*\).

## Verdict

The claimed finite stability refinement is proved, with a stronger
coefficient.  It is an orientation-even theorem and survives exact
self-complementarity.  It is nevertheless a factor \(n\) too weak to
force conference structure, and the self-complementary orbit ensemble
shows that this is a genuine scale obstruction rather than a missing
Taylor constant.

The usable dichotomy is:

1. diffuse Frobenius defect, detected at leading \(A^2\)-energy scale
   by (8), and for the exact field-plus-spin witnesses by (17);
2. heavy local-field tails or sparse macroscopic row correlations,
   isolated quantitatively by (19)--(23) and requiring peeling; or
3. a low-dimensional, localized super-\(\sqrt n\) spectral tail,
   quantified by (6)--(7), again requiring peeling.
