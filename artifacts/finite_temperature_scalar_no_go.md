# Scalar finite-temperature composition does not control the SK-scale diagonal

Status: **proved abstract no-go theorem and exact finite-order falsifier**.
This note does not disprove convergence of the actual signing pressure. It
proves that the proposed annealed bridge composition, even after exact
centering and with the usual scalar regularity properties added, cannot by
itself prove that convergence.

## 1. Exact normalization and bridge centering

For an order-\(n\) signing \(A\), put

$$
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad E_n=\binom n2,
$$

and use the normalized pressure

$$
\overline Z_n(A,t)=2^{-n}\sum_{x\in\{\pm1\}^n}\cosh(tH_A(x)),
\qquad
F_n(t)=\min_A\log\overline Z_n(A,t).
\tag{1.1}
$$

The pressure in the proposed formulation is therefore

$$
\Phi_n(\beta)
=\left(1+\frac1n\right)\log2+
 \frac1nF_n\left(\frac{\beta}{\sqrt n}\right).
\tag{1.2}
$$

Indeed, its unnormalized partition function is
\(2^{n+1}\overline Z_n(A,\beta/\sqrt n)\).

Split \(N=m+n\), choose internal signings \(A,D\), a bridge
\(B\in\{\pm1\}^{m\times n}\), and an orientation
\(\epsilon\in\{\pm1\}\) for \(D\). For fixed \(x,y\), uniform independent
bridge signs give

$$
\mathbb E_B\cosh\left(t\{H_A(x)+\epsilon H_D(y)+x^{\mathsf T}By\}\right)
=(\cosh t)^{mn}
\cosh\left(t\{H_A(x)+\epsilon H_D(y)\}\right).
$$

Averaging also over \(\epsilon\) and using

$$
\frac{\cosh(a+b)+\cosh(a-b)}2=\cosh a\cosh b
$$

gives the exact factorization

$$
\mathbb E_{\epsilon,B}\overline Z_N
=(\cosh t)^{mn}\overline Z_m(A,t)\overline Z_n(D,t).
\tag{1.3}
$$

Some deterministic pair \((\epsilon,B)\) is no larger than its average, so

$$
F_{m+n}(t)\le F_m(t)+F_n(t)+mn\log\cosh t.
\tag{1.4}
$$

Thus the exact edge-annealed centering

$$
R_n(t)=F_n(t)-E_n\log\cosh t
\tag{1.5}
$$

is genuinely subadditive at every fixed physical parameter:

$$
\boxed{R_{m+n}(t)\le R_m(t)+R_n(t).}
\tag{1.6}
$$

This incorporates the full \(\log\cosh t\) term, not merely its quadratic
approximation.

Define the required diagonal quantity

$$
r_n(\beta)=\frac1nR_n\left(\frac{\beta}{\sqrt n}\right).
\tag{1.7}
$$

Then

$$
\Phi_n(\beta)
=\left(1+\frac1n\right)\log2
+\frac{E_n}{n}\log\cosh\frac{\beta}{\sqrt n}
+r_n(\beta),
\tag{1.8}
$$

where the first two terms converge to
\(\log2+\beta^2/4\). Consequently, convergence of \(\Phi_n(\beta)\) is
exactly the unresolved convergence of \(r_n(\beta)\).

For \(\theta=m/(m+n)\), (1.6) gives only

$$
r_{m+n}(\beta)
\le
\theta r_m(\beta\sqrt\theta)
+(1-\theta)r_n(\beta\sqrt{1-\theta}).
\tag{1.9}
$$

In particular,

$$
r_{2n}(\beta)\le r_n(\beta/\sqrt2).
\tag{1.10}
$$

Exact centering removes the extensive bridge defect, but it does not remove
the contraction of the scaled inverse temperature.

There is also a universal reverse comparison, but it is much weaker. For
every fixed block signing, condition on \(x\) and average over \(y\). Both
\(H_D(y)\) and \(x^{\mathsf T}By\) have uniform \(y\)-mean zero, so Jensen's
inequality gives

$$
\overline Z_{m+n}(A,B,D;t)\ge\overline Z_m(A;t).
$$

Interchanging the blocks gives the analogous inequality for \(D\), and
hence

$$
\boxed{F_{m+n}(t)\ge\max\{F_m(t),F_n(t)\}.}
\tag{1.11}
$$

After subtracting the edge-annealed term, (1.11) has a leading-order gap at
\(t=\beta/\sqrt{m+n}\), so it does not match (1.6).

## 2. An analytic oscillating countermodel

The following theorem isolates what scalar composition cannot accomplish.

**Theorem 2.1 (scalar-pressure no-go).** There are abstract pressures
\(F_n^*:\mathbb R\to[0,\infty)\), \(n\ge2\), with all of the following
properties:

1. \(F_n^*\) is even, real analytic, and convex, with
   \(F_n^*(0)=0\) and \((F_n^*)''(0)=E_n\).
2. Its zero-temperature slope \(a_n=\lim_{t\to\infty}F_n^*(t)/t\)
   is positive, nondecreasing, and has the form
   \(a_n=c_n n^{3/2}\), where \(c_n\) does not converge.
3. It obeys the same spin-entropy squeeze as a normalized partition
   function:

   $$
   a_n|t|-n\log2\le F_n^*(t)\le a_n|t|.
   \tag{2.1}
   $$

4. \(F_n^*(t)\) is nondecreasing in \(n\) for each fixed \(t\ge0\), as in
   the restriction/Jensen comparison (1.11).
5. The exactly centered functions

   $$
   R_n^*(t)=F_n^*(t)-E_n\log\cosh t
   $$

   obey \(R_{m+n}^*(t)\le R_m^*(t)+R_n^*(t)\).
6. The diagonal functions

   $$
   \beta\longmapsto
   \frac1nF_n^*(\beta/\sqrt n)
   $$

   are uniformly Lipschitz, on all of \([0,\infty)\), but for every
   fixed \(\beta>0\) their values fail to converge as \(n\to\infty\).

**Construction.** Fix

$$
0<c_-<c_+<\frac1{2\sqrt2}.
$$

Choose \(c_0,\varepsilon>0\) with
\(c_0-\varepsilon=c_-\) and \(c_0+\varepsilon=c_+\), and, for \(n\ge3\),
set

$$
c_n=c_0+\varepsilon\sin(\log\log(n+n_0)).
\tag{2.2}
$$

Taking \(n_0\) sufficiently large makes the oscillation slow enough that

$$
a_n=c_nn^{3/2}\ \text{increases},\qquad
\vartheta_n=\frac{a_n}{E_n}\ \text{decreases},\qquad
L_n=\frac{a_n^2}{E_n}\ \text{increases}.
\tag{2.3}
$$

for every \(n\ge3\). Define the exceptional initial value explicitly by

$$
\vartheta_2=\vartheta_3,\qquad
a_2=E_2\vartheta_2=\vartheta_3,\qquad
c_2=\frac{a_2}{2^{3/2}},\qquad
L_2=E_2\vartheta_2^2=\vartheta_3^2.
\tag{2.4}
$$

Then \(a_n\) and \(L_n\) are nondecreasing and \(\vartheta_n\) is
nonincreasing for every \(n\ge2\). Moreover, \(\vartheta_n<1\) and
\(L_n\le n\): for \(n\ge3\),

$$
\vartheta_n=\frac{2c_n\sqrt n}{n-1}
\le c_+\sqrt3<1,
\qquad
\frac{L_n}{n}=\frac{2c_n^2n}{n-1}
\le3c_+^2<1,
\tag{2.5}
$$

and the order-two claims follow from \(\vartheta_2=\vartheta_3\).
Now put

$$
\boxed{
F_n^*(t)=
\frac{a_n^2}{E_n}
\log\cosh\left(\frac{E_n}{a_n}t\right)
=L_n\log\cosh(t/\vartheta_n).
}
\tag{2.6}
$$

For completeness, the monotonicities asserted in (2.3) follow
by differentiating the continuous interpolation of (2.2). Since

$$
\left|\frac{c'(x)}{c(x)}\right|
\le
\frac{\varepsilon}
{c_-(x+n_0)\log(x+n_0)},
$$

\(n_0\) can be chosen so that this is at most \(1/(8x)\) for all
\(x\ge2\). The logarithmic derivatives are

$$
\frac{a'}a=\frac{c'}c+\frac{3}{2x},
$$

$$
\frac{\vartheta'}{\vartheta}
=\frac{c'}c+\frac1{2x}-\frac1{x-1},
$$

and

$$
\frac{L'}L
=2\frac{c'}c+\frac2x-\frac1{x-1}.
$$

The first is positive and the second is negative for \(x\ge3\). The
third is positive there as well: its last two terms have value at least
\(1/(2x)\), whereas \(2c'/c\ge-1/(4x)\). The explicit order-two
definition then completes all three monotonicity claims. Notice also that
\(E_n/a_n=1/\vartheta_n\) increases.

**Proof of the analytic properties and squeeze.** Formula (2.6) is even,
analytic, convex, and zero at the origin. Direct differentiation gives

$$
(F_n^*)''(0)
=L_n\vartheta_n^{-2}
=E_n,
$$

while

$$
\lim_{t\to\infty}\frac{F_n^*(t)}t
=\frac{L_n}{\vartheta_n}
=a_n.
$$

The elementary inequality

$$
|u|-\log2\le\log\cosh u\le|u|
$$

gives

$$
a_n|t|-L_n\log2\le F_n^*(t)\le a_n|t|.
$$

Since \(L_n\le n\), this proves (2.1). Because both \(L_n\) and
\(1/\vartheta_n\) are nondecreasing, (2.6) also proves that
\(F_n^*(t)\) is nondecreasing in \(n\) for every fixed \(t\ge0\).

**Proof of exact centered subadditivity.** For \(0<\vartheta\le1\), define

$$
h_\vartheta(t)
=\vartheta^2\log\cosh(t/\vartheta)-\log\cosh t.
\tag{2.7}
$$

For \(u=t/\vartheta\),

$$
\frac{\partial}{\partial\vartheta}
\left\{\vartheta^2\log\cosh(t/\vartheta)\right\}
=\vartheta\{2\log\cosh u-u\tanh u\}.
\tag{2.8}
$$

The expression in braces is nonnegative. Indeed, if

$$
k(u)=2\log\cosh u-u\tanh u,
$$

then

$$
k'(u)=\tanh u-u\,\operatorname{sech}^2u,
$$

and the derivative of the right side is
\(2u\operatorname{sech}^2u\tanh u\ge0\) for \(u\ge0\). It vanishes at
zero, so \(k'(u)\ge0\) and then \(k(u)\ge0\). Evenness covers \(u<0\).
Thus \(h_\vartheta(t)\) is nondecreasing in \(\vartheta\). Since
\(h_1(t)=0\),

$$
h_\vartheta(t)\le0\qquad(0<\vartheta\le1).
\tag{2.9}
$$

Now

$$
R_n^*(t)=E_nh_{\vartheta_n}(t).
$$

For \(N=m+n\), monotonicity of \(\vartheta_n\) gives
\(\vartheta_N\le\vartheta_m,\vartheta_n\), while

$$
E_N=E_m+E_n+mn.
$$

Using first (2.9) and then monotonicity in \(\vartheta\),

$$
\begin{aligned}
R_N^*(t)
&=(E_m+E_n+mn)h_{\vartheta_N}(t)\\
&\le(E_m+E_n)h_{\vartheta_N}(t)\\
&\le E_mh_{\vartheta_m}(t)+E_nh_{\vartheta_n}(t)\\
&=R_m^*(t)+R_n^*(t).
\end{aligned}
\tag{2.10}
$$

This proves the exact counterpart of (1.6).

**Proof of diagonal oscillation and equicontinuity.** From (2.6),

$$
\frac1nF_n^*\left(\frac{\beta}{\sqrt n}\right)
=
\frac{2c_n^2n}{n-1}
\log\cosh\left(
\frac{\beta(n-1)}{2c_nn}
\right).
\tag{2.11}
$$

Uniformly for \(c_n\in[c_-,c_+]\), this differs by \(o(1)\) from

$$
J_\beta(c_n)
=2c_n^2\log\cosh\left(\frac{\beta}{2c_n}\right).
\tag{2.12}
$$

For \(u=\beta/(2c)>0\),

$$
J_\beta'(c)
=2c\{2\log\cosh u-u\tanh u\}>0,
\tag{2.13}
$$

where strict positivity follows from the argument following (2.8). The slowly
oscillating sequence (2.2) has subsequences tending to \(c_-\) and \(c_+\):
one may take the nearest integers to the real solutions of
\(\log\log(x+n_0)=2\pi j+3\pi/2\) and
\(2\pi j+\pi/2\), since the derivative of \(\log\log(x+n_0)\) tends to
zero. Hence (2.11) has two distinct subsequential limits for every
\(\beta>0\).

Finally,

$$
\frac{\partial}{\partial\beta}
\left[
\frac1nF_n^*\left(\frac{\beta}{\sqrt n}\right)
\right]
=c_n\tanh\left(\frac{\beta(n-1)}{2c_nn}\right),
\tag{2.14}
$$

whose absolute value is at most \(c_+\). Thus the entire diagonal family
is uniformly Lipschitz. Subtracting

$$
\frac{E_n}{n}\log\cosh(\beta/\sqrt n)\longrightarrow\frac{\beta^2}{4}
$$

shows that the centered diagonal \(R_n^*(\beta/\sqrt n)/n\) oscillates as
well. This completes the proof of Theorem 2.1.

The construction is abstract: it is not asserted to arise from complete
quadratic signings. Its point is precisely that any successful proof must
use a property that distinguishes complete quadratic signing pressures
from scalar functions satisfying all properties above.

## 3. Exact order-four obstruction to quadratic scale transport

One tempting way to repair the contracted parameter in (1.9) is a
quadratic scaling inequality for the centered pressure. It already fails
for actual complete signings at order four.

Write \(\rho=\tanh t\). The high-temperature expansion gives

$$
\frac{\overline Z_4(A,t)}{(\cosh t)^6}
=
\sum_{\substack{F\subseteq E(K_4)\\
                 \deg_F(v)\ {\rm even\ for\ all}\ v\\
                 |F|\ {\rm even}}}
\left(\prod_{e\in F}a_e\right)\rho^{|F|}.
\tag{3.1}
$$

Besides the empty graph, the only contributing subgraphs are the three
Hamilton four-cycles. Let their sign products be \(p_1,p_2,p_3\). Every
edge of \(K_4\) lies in exactly two of these cycles, so

$$
p_1p_2p_3=1.
$$

Their sum is consequently either \(3\) or \(-1\). The latter is attained,
for example, by making exactly one edge negative. Therefore

$$
\boxed{
R_4(t)=\log(1-\rho^4).
}
\tag{3.2}
$$

For \(0<\lambda,\rho<1\), convexity of
\(L(s)=-\log(1-s)\), together with \(L(0)=0\), gives

$$
L(\lambda^4\rho^4)
\le\lambda^4L(\rho^4)
<\lambda^2L(\rho^4).
$$

Equivalently,

$$
\boxed{
R_4(\operatorname{arctanh}(\lambda\rho))
>
\lambda^2R_4(\operatorname{arctanh}\rho).
}
\tag{3.3}
$$

Thus the scale-transport inequality with the opposite sign,

$$
R_n(\operatorname{arctanh}(\lambda\rho))
\le\lambda^2R_n(\operatorname{arctanh}\rho),
$$

is false already at \(n=4\). That inequality is the natural direction for
turning fixed-\(\rho\) centered subadditivity into upper control on a
shrinking-\(\rho\) diagonal.

## 4. What remains a genuinely new theorem

The proposed finite-temperature criterion remains correct: convergence of
\(\Phi_n(\beta)\) for every fixed \(\beta>0\) would imply convergence of
\(M_n/n^{3/2}\). What fails is the claim that the annealed bridge law puts
that criterion within reach of scalar Fekete theory.

The exact bridge calculation has now been pushed to its strongest scalar
form:

- the favorable child orientation removes the harmless factor of two;
- exact edge annealing removes the entire bridge defect;
- centered subadditivity still follows the characteristic
  \(\beta\mapsto\beta\sqrt{m/(m+n)}\);
- restriction Jensen, convexity, the ground-state squeeze, and uniform
  equicontinuity do not prevent diagonal oscillation;
- the simplest quadratic scale-transport inequality is false for a real
  signing.

A successful finite-temperature proof therefore needs an additional
signing-specific theorem not present in these scalar data. Examples of
sufficiently new information would be a matching cross-block lower law
with a subextensive error, a coefficient-level Eulerian constraint forcing
uniqueness of the diagonal pressure, or an overlap/correlation state whose
interpolation has a definite sign. Proving convergence of the scalar
\(R_n(\beta/\sqrt n)/n\) directly would of course suffice, but merely
restating that convergence is not a reduction.
