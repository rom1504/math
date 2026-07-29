# Positive heavy fields at a global Boolean maximizer

Let \(D\) be a symmetric zero-diagonal sign matrix of order \(m\), and
suppose (after switching and, if needed, replacing \(D\) by \(-D\)) that

\[
R=\mathbf 1^{\mathsf T}D\mathbf 1=Q(D)
 :=\max_{x\in\{\pm1\}^m}|x^{\mathsf T}Dx|.
\]

Write \(r=D\mathbf 1\). Then \(r_i\ge 0\), and every signed cut has weight
between \(0\) and \(R/2\).

## 1. Exact universal-vertex obstruction

For \(k\ge 1\), set

\[
\widetilde D_k=
\begin{pmatrix}
J_k-I_k & J_{k,m}\\
J_{m,k} & D
\end{pmatrix}.
\]

Then

\[
\boxed{
Q(\widetilde D_k)
=R+2km+k(k-1).
}
\]

Indeed, if \(z\in\{\pm1\}^k\), \(y\in\{\pm1\}^m\),
\(s=\sum z_i\), and \(t=\sum y_i\), then

\[
\begin{aligned}
\left|
\binom zy^{\mathsf T}\widetilde D_k\binom zy
\right|
&=\left|s^2-k+2st+y^{\mathsf T}Dy\right|\\
&\le k(k-1)+2km+R,
\end{aligned}
\]

and equality holds at \(z=\mathbf1,y=\mathbf1\).

Consequently, \(\widetilde D_k\) is again globally cut-stable at
\(\mathbf1\). Its row sums are

\[
\widetilde r_i=m+k-1 \quad (i\le k),
\qquad
\widetilde r_{k+j}=r_j+k.
\]

In particular, for one added vertex,

\[
\sum_i\widetilde r_i^2
=m^2+\sum_j r_j^2+2R+m.
\]

Thus one can add a row of size \(m\), and add \(m^2+o(m^2)\) to the
row-square mass, while changing \(Q\) by only \(2m=o(m^{3/2})\).
More generally, any \(k=o(\sqrt m)\) adds arbitrarily many heavy positive
rows without changing the leading normalized energy:

\[
\frac{Q(\widetilde D_k)}{(m+k)^{3/2}}
-\frac{Q(D)}{m^{3/2}}\longrightarrow0.
\]

This is an actual signed-complete-graph construction, not a moment
relaxation. Therefore no argument whose only extra inputs are

- all cuts are nonnegative,
- the mean row sum is \(c\sqrt n\), and
- the row-square mean is at least \(n\),

can turn the mere existence of a heavy positive row into a strict
leading-order improvement. Exceptional heavy rows can be manufactured
at \(o(n^{3/2})\) energy cost.

## 2. Exact exhaustion of row-threshold information for a conference
matrix

Now additionally suppose \(D^2=(n-1)I\), put

\[
U=\frac D{\sqrt{n-1}},\qquad
u=U\mathbf1,\qquad
\rho=\frac{\mathbf1^{\mathsf T}D\mathbf1}
{n\sqrt{n-1}}.
\]

With empirical inner product \(n^{-1}\sum_i f_i g_i\),

\[
\mathbb Eu=\rho,\qquad \mathbb Eu^2=1,\qquad
U\mathbf1=u,\qquad Uu=\mathbf1.
\]

For any Boolean \(x\), define

\[
\alpha=\mathbb Ex,\qquad \beta=\mathbb Exu.
\]

The vector

\[
w_-=\frac{\mathbf1-u}{\sqrt{2(1-\rho)}}
\]

is a unit \((-1)\)-eigenvector of \(U\). Hence

\[
\frac1n x^{\mathsf T}Ux
\le
1-\frac{(\alpha-\beta)^2}{1-\rho}.
\]

Since global absolute maximality gives
\(x^{\mathsf T}Ux/n\ge-\rho\), every Boolean \(x\) satisfies

\[
|\alpha-\beta|\le\sqrt{1-\rho^2}.
\]

Optimizing over all Boolean functions of the row sum is attained by the
single threshold \(x_i=\operatorname{sign}(1-u_i)\), and gives the exact
necessary condition

\[
\boxed{
\mathbb E|1-u|\le\sqrt{1-\rho^2}.
}
\]

This completely exhausts what the two-dimensional
\(\operatorname{span}\{\mathbf1,u\}\) conference geometry can obtain
from arbitrary row thresholds.

It still cannot control rare heavy rows. If \(0\le u\le L\),
\(\mathbb Eu=\rho\), and \(\mathbb Eu^2=1\), then (when \(L\rho\ge1\))

\[
\inf \mathbb E|1-u|
=(1-\rho)\left(1+\frac2L\right).
\]

The extremizer is supported on \(\{0,1,L\}\). The proof is the identity

\[
\mathbb E|1-u|=1-\rho+2\mathbb E(u-1)_+
\]

together with

\[
1-\rho=\mathbb E(u^2-u)
\le L\,\mathbb E(u-1)_+,
\]

with equality for that three-point law. Since a conference row can have
\(L=\sqrt{n-1}\), this lower bound tends only to \(1-\rho\), and the
threshold constraint becomes asymptotically vacuous.

## 3. The statistic that remains

The correct unresolved datum is not the first two row moments, but
uniform integrability of the row-square profile. In normalized variables
\(u_i=r_i/\sqrt n\), define

\[
\Psi_n(K)=\frac1n\sum_{i=1}^n
u_i^2\mathbf1_{\{u_i>K\}}.
\]

The universal-vertex extension has a nonvanishing \(\Psi_n(K)\) for
every fixed \(K\), while preserving the leading \(n^{3/2}\) objective.
A viable heavy-field proof therefore needs a genuine two-case theorem:

1. a regular branch controlling configurations for which
   \(\lim_{K\to\infty}\limsup_n\Psi_n(K)=0\);
2. a peeling/structural branch showing that concentration of
   \(\Psi_n(K)\) can be removed or converted using information beyond
   cut positivity and the first two moments.

Threshold cuts, independent correlated cuts based only on \(r_i\), and
the two-dimensional conference involution do not supply the second
branch.
