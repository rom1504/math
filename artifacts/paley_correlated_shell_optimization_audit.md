# Correlated-shell audit of sparse Paley perturbation

## Status

This note independently checks Sections 7--9 of
`perturbed_conference_stratified_entropy.md` and corrects the interpretation
of the Hamming-cloud obstruction in
`paley_energy_layer_cloud_obstruction.md`.

The conclusions are:

1. the absolute per-configuration entropy union bound, even after inserting
   the exact Hamming-cloud entropy, can **never** certify a constant below
   \(1/2\);
2. this is only a failure of that scalar proof architecture;
3. the exact correlated shell process has a much smaller width;
4. conditional on the ideal deterministic shell profile, the rigorous
   Gaussian-width upper bound would yield \(0.454258819884\ldots\);
5. the sole missing input is now a deterministic upper bound on the maximum
   Paley energy within each shell around every resonant center.

## 1. The absolute cloud union bound has exact minimum \(1/2\)

For a cap seed, the independent-sample cloud expression is

\[
B_{\rm abs}(\delta)
=\max_{0\le\theta\le1/2}
\left[
\frac{(1-2\delta)(1-2\theta)^2}{2}
+2\sqrt{\delta(1-\delta)h(\theta)}
\right]. \tag{1}
\]

Numerically, \(B_{\rm abs}(\delta)>1/2\) for every
\(\delta\in(0,1/2]\).  This has a short analytic proof.

Take \(\theta=\delta=d\), put \(\mu=1-2d\), and define

\[
P(d)=3-6d+4d^2.
\]

Since

\[
1-\mu^3=2dP(d),
\]

the inequality

\[
\frac{\mu^3}{2}+2\sqrt{d(1-d)h(d)}\ge\frac12 \tag{2}
\]

is equivalent, after squaring nonnegative quantities, to

\[
h(d)\ge q(d):=\frac{dP(d)^2}{4(1-d)}. \tag{3}
\]

Let \(\psi=h-q\).  Polynomial division gives

\[
q(d)=-4d^4+8d^3-7d^2+2d-\frac14+\frac1{4(1-d)}.
\]

Writing \(t=d(1-d)\),

\[
\psi''(d)
=14-48t-\frac1t-\frac1{2(1-d)^3}. \tag{4}
\]

By AM--GM,

\[
48t+\frac1t\ge2\sqrt{48}=8\sqrt3,
\]

and \(1/[2(1-d)^3]\ge1/2\).  Hence

\[
\psi''(d)\le\frac{27}{2}-8\sqrt3<0
\qquad(0<d\le1/2). \tag{5}
\]

The function \(\psi\) is concave, with

\[
\psi(0)=0,\qquad
\psi(1/2)=\log2-\frac14>0.
\]

It lies above its endpoint chord, so \(\psi(d)>0\) for \(d>0\).
Therefore

\[
\boxed{\displaystyle
\inf_{0\le\delta\le1/2}B_{\rm abs}(\delta)=\frac12,} \tag{6}
\]

with equality only at \(\delta=0\).

This is a rigorous stopping theorem for every proof that treats the
Hamming-cloud perturbations as independent absolute variables.

## 2. Exact covariance of the actual perturbation process

Let

\[
Z_x=\sum_{i<j}c_{ij}x_ix_j(\eta_{ij}-\mu),
\qquad
\operatorname{Var}(\eta_{ij}-\mu)=4D,
\quad D=\delta(1-\delta).
\]

For \(x,y\in\{\pm1\}^n\), put \(z_i=x_iy_i\).  Then

\[
\begin{aligned}
\operatorname{Cov}(Z_x,Z_y)
&=4D\sum_{i<j}z_iz_j\\
&=2D\left[\left(\sum_i z_i\right)^2-n\right]. \tag{7}
\end{aligned}
\]

Since

\[
\operatorname{Var}Z_x=2Dn(n-1),
\]

configurations at Hamming distance \(\alpha n\) have correlation

\[
\operatorname{Corr}(Z_x,Z_y)
=\frac{(n-2\alpha n)^2-n}{n(n-1)}
=(1-2\alpha)^2+o(1). \tag{8}
\]

At the radius

\[
\alpha_*=\frac{1-\sqrt{r_*}}2
=0.00800258218364\ldots ,
\]

this correlation is approximately \(r_*=0.9682458\).  The cloud variables
are therefore very far from independent.

For a center \(x\), write \(y=x\odot s\), with
\(S=\{i:s_i=-1\}\).  The common noise cancels:

\[
Z_y-Z_x
=-2\sum_{i\in S,\ j\notin S}
c_{ij}x_ix_j(\eta_{ij}-\mu). \tag{9}
\]

For two \(k\)-subsets \(S,T\), \(d=|S\triangle T|\), the Gaussianized
version of (9) has the exact increment metric

\[
\boxed{\displaystyle
\mathbb E(Y_S-Y_T)^2=16D\,d(n-d).} \tag{10}
\]

This verifies the canonical metric stated in the perturbation note.

## 3. Independent verification of the Gaussian-width sandwich

Let \(g_1,\ldots,g_n\) be independent standard Gaussians and compare \(Y\)
with

\[
L_S^+=4\sqrt{Dn}\sum_{i\in S}g_i,\qquad
L_S^-=4\sqrt{D(n-2k)}\sum_{i\in S}g_i.
\]

If \(d=|S\triangle T|\le2k\), then

\[
16D(n-2k)d
\le16Dd(n-d)
\le16Dnd.
\]

Sudakov--Fernique therefore gives, for \(k/n\to\alpha\le1/2\),

\[
\boxed{\displaystyle
4\sqrt{D(1-2\alpha)}\,\phi(z_\alpha)
\le
\liminf\frac{\mathbb E\max_{|S|=k}Y_S}{n^{3/2}}
\le
\limsup\frac{\mathbb E\max_{|S|=k}Y_S}{n^{3/2}}
\le
4\sqrt D\,\phi(z_\alpha),} \tag{11}
\]

where

\[
z_\alpha=\Phi^{-1}(1-\alpha).
\]

Indeed, the sum of the largest \(\alpha n\) among \(n\) standard normals,
divided by \(n\), converges to

\[
\int_{z_\alpha}^{\infty}t\phi(t)\,dt=\phi(z_\alpha).
\]

The replacement of Gaussian edges by centered Bernoulli edge flips is
valid for fixed \(\delta>0\).  For the soft maximum

\[
\beta^{-1}\log\sum_{|S|=k}e^{\beta Y_S},
\]

the entropy error is \(O(n/\beta)\) and the summed third-order Lindeberg
error is \(O(n^2\beta^2)\).  The choice
\(\beta=n^{-3/8}\) makes both \(o(n^{3/2})\).

At \(\alpha_*\), direct recomputation gives

\[
\begin{array}{c|c}
\text{penalty divided by }\sqrt D & \text{coefficient}\\ \hline
\text{absolute scalar union} & 2\sqrt{h(\alpha_*)}
  =0.431772489417\ldots\\
\text{increment union} &
\sqrt{32\alpha_*(1-\alpha_*)h(\alpha_*)}
  =0.108810374085\ldots\\
\text{Gaussian lower comparison} &
4\sqrt{1-2\alpha_*}\,\phi(z_{\alpha_*})
  =0.086998464392\ldots\\
\text{Gaussian upper comparison} &
4\phi(z_{\alpha_*})
  =0.087703147047\ldots
\end{array} \tag{12}
\]

Thus retaining correlation removes roughly eighty percent of the original
scalar entropy charge.

## 4. Conditional optimization with the ideal shell profile

Suppose, only for this diagnostic calculation, that every cap center obeyed
the ideal deterministic shell upper profile

\[
u(\alpha)=\frac12(1-2\alpha)^2. \tag{13}
\]

Combining (13) with the rigorous upper side of (11) gives

\[
\mathcal G_{\rm up}(\delta)
=\max_{0\le\alpha\le1/2}
\left[
\frac{(1-2\delta)(1-2\alpha)^2}{2}
+4\sqrt{\delta(1-\delta)}\,\phi(z_\alpha)
\right]. \tag{14}
\]

Independent numerical optimization gives

\[
\boxed{
\begin{aligned}
\delta_{\rm opt}&=0.06297045\ldots,\\
\alpha_{\rm worst}&=0.05443909\ldots,\\
\min_\delta\mathcal G_{\rm up}(\delta)
&=0.454258819884\ldots .
\end{aligned}} \tag{15}
\]

For comparison, replacing the true correlated upper width by the increment
union bound gives

\[
\begin{aligned}
\delta_{\rm opt}^{\rm inc}&=0.02971507\ldots,\\
\alpha_{\rm worst}^{\rm inc}&=0.02700412\ldots,\\
\min_\delta\mathcal G_{\rm inc}(\delta)
&=0.475722453966\ldots .
\end{aligned} \tag{16}
\]

The lower Gaussian comparison in (11) gives the diagnostic value
\(0.445679375462\ldots\), but it is not an upper bound and therefore cannot
certify a construction.

The constants in (15)--(16) are conditional diagnostics, not theorems about
Paley signings.  Equation (13) is the shell **mean** around an exact cap
center; the required input is an upper bound on the shell maximum.

## 5. The exact remaining deterministic lemma

For every relevant resonant center \(x\), define

\[
u_{C,x}(\alpha)
=\frac1{n^{3/2}}
\max_{d_H(x,y)=\alpha n}H_C(y).
\]

The sparse-perturbation construction would give a strict constant below
\(1/2\) if one could prove a uniform profile \(\bar u\) such that

\[
u_{C,x}(\alpha)\le\bar u(\alpha)+o(1)
\]

for all cap-near centers and

\[
\inf_{\delta>0}\sup_{0\le\alpha\le1/2}
\left[
(1-2\delta)\bar u(\alpha)
+4\sqrt{\delta(1-\delta)}\,\phi(z_\alpha)
\right]
<\frac12. \tag{17}
\]

Pure spectral information gives only \(\bar u(\alpha)\le1/2\), which is too
weak.  The average shell identity gives (13), which is strong enough
numerically but does not control the maximum.

Thus the Paley inverse problem has been reduced from counting exceptional
vectors to a concrete local rigidity statement:

> **Shell-max target.**  Prove that a cap-near Paley vector cannot have
> another cap-near vector at every prescribed positive Hamming radius, with
> a quantitative profile strong enough for (17).

## 6. Verdict

- The naive full-cloud entropy branch has exact optimum \(1/2\), so that
  route is rigorously stopped.
- The exact covariance, metric, Gaussian comparison bounds, universality
  scaling, and coefficients in the updated perturbation note all check.
- Correlation reopens sparse perturbation; it is not killed by the
  exponential Hamming cloud.
- The next proof target is entirely deterministic: a Paley
  shell-maximum/inverse theorem around resonant centers.

