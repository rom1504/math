# Bounded-degree compatible graphing

## Status

This file is the durable checkpoint for the bounded-degree route.
There are three logically separate results:

1. the Pietsch/threshold reduction is rigorous;
2. the small-edge perturbative theorem is rigorous after smoothing
   both sign steps;
3. global monotonicity along the compatible homotopy is currently a
   conjecture under audit.  The numerical data below are explicitly
   not a proof.

The doubled normalization is
\[
Q(A)=\max_{x\in\{\pm1\}^n}|x^\top Ax|,
\qquad c_2=0.783387533648\ldots.
\]

## 1. Rigorous bounded-degree reduction

For a signing \(A\), delete at most \(\varepsilon n\) vertices by
Grothendieck--Pietsch and call the remaining order \(r\).  With
\[
B=\frac{A_R}{\sqrt{r-1}},\qquad C=B^2,\qquad
C_{ii}=1,\qquad q_{ij}=C_{ij},
\]
one has
\[
\|C\|_{\rm op}\le
L_\varepsilon+o(1),\qquad
L_\varepsilon=
\frac{(4K_Gc_2)^2}{\varepsilon^2(1-\varepsilon)}.
\tag{1}
\]
Consequently
\[
\sum_jq_{ij}^2=(C^2)_{ii}\le L_\varepsilon+o(1).
\tag{2}
\]
At threshold
\[
\zeta_{\varepsilon,\eta}
=\sqrt{\eta/L_\varepsilon},
\]
the graph of edges \(|q_{ij}|\ge\zeta_{\varepsilon,\eta}\) has
\[
\Delta_{\varepsilon,\eta}
\le L_\varepsilon^2/\eta,
\tag{3}
\]
while the deleted correlations satisfy
\[
\frac1r\sum_{|q_{ij}|<\zeta_{\varepsilon,\eta}}q_{ij}^4
\le\eta+o(1).
\tag{4}
\]
Also \(Q(A)\ge Q(A_R)\), so this costs only
\((1-\varepsilon)^{3/2}\) in the final normalization.

Thus all non-negligible rank-three residual correlations live on a
uniformly bounded-degree, finitely labelled graph after the order of
limits
\[
r\to\infty,\quad \eta\downarrow0,\quad\varepsilon\downarrow0.
\]

## 2. Compatible transport

The residual covariance is
\[
K=s^2I+\kappa[Q],\qquad
\kappa(q)=\sum_{\substack{\ell\ge3\\\ell\ {\rm odd}}}
u_\ell^2q^\ell.
\tag{5}
\]
The fresh field covariance is
\[
D=BKB.
\tag{6}
\]
Although \(B\) is dense, it is not an uncontrolled extra graph.  For
the row cuts
\[
z^{(i)}_j=\sqrt{r-1}\,B_{ij}\in\{\pm1\},
\]
\[
q_{jk}=\frac1r\sum_i z^{(i)}_jz^{(i)}_k+o(1)
\tag{7}
\]
and
\[
D_{ii}
=s^2+\frac2r\sum_{j<k}\kappa(q_{jk})
 z^{(i)}_jz^{(i)}_k+o(1).
\tag{8}
\]
In particular
\[
\frac1r\operatorname{tr}D
=s^2+\frac2r\sum_{j<k}q_{jk}\kappa(q_{jk})+o(1)
\ge s^2.
\tag{9}
\]

## 3. Audited positive edge Hessian

Smooth the final sign by
\[
\psi_\tau(u)=2\Phi(u/\tau)-1.
\]
Take a nonzero stationary scalar response
\[
y=\psi_\tau(U),\qquad
U=hS+aG+jR+mZ,
\]
where
\[
a=\mathbb ESy,\quad h=\mathbb EGy,\quad
m=\mathbb ERy,\quad j=\mathbb EZy.
\]
Define
\[
\ell_\tau=\mathbb E\psi_\tau'(U)>0,\qquad
k_\tau=-\mathbb E[U\psi_\tau''(U)]
=\frac2{\tau^3}\mathbb E[U^2\phi(U/\tau)]>0.
\tag{10}
\]

For \(C=I+Q\), use the neighbor-spin perturbation
\[
y_i=\psi_\tau\left(
U_i+\lambda h\sum_{j\ne i}q_{ij}S_j
\right).
\tag{11}
\]
A line-by-line response-matrix expansion gives
\[
\mathcal E(Q)
=c_{2,\tau}
+h^2(2\lambda\ell_\tau-\lambda^2k_\tau)
\frac{\operatorname{tr}Q^2}{r}
+O(Q^3).
\tag{12}
\]
No residual covariance derivative is hidden in (12):

- \(\mathcal A_{ij}=\lambda h\ell_\tau q_{ij}+O(Q^3)\);
- \(\mathcal H,\mathcal J\) are diagonal through degree two by
  Gaussian integration by parts;
- \(\mathcal M_{ij}=O(Q^3)\) off the diagonal because the residual
  has cross-site Hermite rank three;
- the four diagonal response shifts combine, by stationarity, into
  \(-\lambda^2h^2k_\tau\operatorname{tr}Q^2/r\).

Taking
\[
\lambda_\tau=\min\{1,\ell_\tau/k_\tau\}
\]
gives the explicit safe coefficient
\[
a_{\tau}^{\rm safe}
=h^2\lambda_\tau\ell_\tau>0.
\tag{13}
\]

The earlier perturbation by
\(\sum_jq_{ij}(\beta S_j+\delta R_j)\) is withdrawn; its claimed
coefficient omitted residual-channel covariance terms.

## 4. Uniform small-edge theorem

Smooth the first sign as well, with parameter \(\sigma>0\), and
regress its \(S,G\) components exactly.  Then
\[
|\kappa_\sigma(q)|\le s_\sigma^2|q|^3.
\tag{14}
\]
Let
\[
\widehat Q=(|q_{ij}|),\qquad
\varrho=\|\widehat Q\|_{\rm op}
\le\Delta q_{\max}.
\tag{15}
\]
Interpolate by
\[
C_u=I+uQ,\qquad
B_u=O\,C_u^{1/2},\qquad O=\operatorname{sgn}(B).
\tag{16}
\]
For \(\varrho\le1/8\), all \(C_u\) are uniformly nonsingular and
\[
\|\kappa_\sigma[uQ]\|_{\rm op}
\le s_\sigma^2q_{\max}^2\varrho,
\tag{17}
\]
\[
\|B_u\kappa_\sigma[uQ]B_u\|_{\rm op}
\le(1+\varrho)s_\sigma^2\varrho^3.
\tag{18}
\]
Equation (18) is the operator-norm control of the dense transport
that was missing from the first cluster argument.

Three Gaussian covariance differentiations, with two edge insertions
in Hilbert--Schmidt norm and every remaining insertion in
\(\widehat Q\)-operator norm, give a dimension-free constant
\(\mathfrak R_{\sigma,\tau,\Delta}<\infty\) such that
\[
\left|
\mathcal E(Q)-c_{2,\sigma,\tau}
-h^2(2\lambda\ell_\tau-\lambda^2k_\tau)
 \frac{\operatorname{tr}Q^2}{r}
\right|
\le
\mathfrak R_{\sigma,\tau,\Delta}\,
\varrho\frac{\operatorname{tr}Q^2}{r}.
\tag{19}
\]
The constant depends only on the two fixed dithers, the scalar
coefficients, and \(\Delta\), not on \(r\) or \(B\).  It is obtained
from derivatives through order six of the two smoothed one-site
maps.  Star contractions from (11) cost a factor depending on
\(\Delta\); the dense \(B\)'s cost only their operator norms.

Therefore, with
\[
\rho_{\sigma,\tau,\Delta}
=\min\left\{
\frac18,\,
\frac{a_{\tau}^{\rm safe}}
{2\mathfrak R_{\sigma,\tau,\Delta}}
\right\},
\tag{20}
\]
the concrete condition
\[
\boxed{\Delta q_{\max}\le\rho_{\sigma,\tau,\Delta}}
\tag{21}
\]
implies
\[
\boxed{
\mathcal E(Q)\ge
c_{2,\sigma,\tau}
+\frac{a_{\tau}^{\rm safe}}2
\frac{\operatorname{tr}Q^2}{r}.
}
\tag{22}
\]
This is uniform in the matrix order.

## 5. Global compatible homotopy: current audit

The proposed globalization is
\[
C_u=I+u(C-I),\qquad
B_u=\operatorname{sgn}(B)\,C_u^{1/2},
\qquad0\le u\le1.
\tag{23}
\]
Let \(V(u)\) be the optimized, dithered finite-channel paired value.
If
\[
V'(u)\ge0
\tag{24}
\]
for every compatible path, then \(V(1)\ge V(0)\ge c_{2,\sigma,\tau}\)
and the bounded-degree graphing gap closes at once.

At a stationary response the envelope theorem removes the derivative
of the response itself.  However, a channelwise I--MMSE proof is not
available: each paired term is a cross bilinear form,
\[
\operatorname{tr}(B\mathcal A B\mathcal H^\top),
\qquad
\operatorname{tr}(B\mathcal M B\mathcal J^\top),
\tag{25}
\]
whose polarization is a difference of two squares.  Positivity, if
true, must come from cancellation between the explicit transport
derivative, the residual regression derivative, and the fresh-field
Onsager derivative.  Those terms are still being simplified.

### Computation, not proof

`_homotopy_finite_channel.py` implements the exact finite-channel
functional and repeated conditional-score optimization along (23).
With common scrambled-Sobol samples:

- every nonsingular symmetric sign fibre of sizes \(2\) and \(3\)
  showed no decrease;
- for
  \[
  R=J_3-2I_3,
  \]
  the computed values at
  \(u=0,.1,.25,.5,.75,1\) were
  \[
  0.783401,\ 0.783953,\ 0.786690,\ 0.795803,\
  0.810877,\ 0.830730.
  \]

A coarse size-four sweep produced one apparent initial decrease, but
the computation used insufficient quadrature/optimization accuracy
and is not yet classified as a counterexample.  The candidate fibre
is
\[
\begin{pmatrix}
-1&1&-1&-1\\
1&1&-1&1\\
-1&-1&1&1\\
-1&1&1&-1
\end{pmatrix}.
\tag{26}
\]
It must be rerun with a smoothed stationary response, common samples,
and interval/error control before drawing any conclusion.

## 6. Remaining target

The immediate target is one of:

1. an exact envelope identity expressing \(V'(u)\) as a positive
   square or Schur complement after all regression/Onsager terms
   cancel; or
2. a certified finite compatible fibre for which \(V'(u)<0\).

Failure of monotonicity would not kill the lower-bound program: (22)
already handles the perturbative background, and the remaining
large-edge graph has bounded degree and edge weights bounded below.
It would then be attacked through finite components, matchings, and
local weak approximation instead of a global homotopy.

## 7. Exact envelope derivative

The full derivative can be written cleanly after whitening.  This
also shows why the hoped-for I--MMSE square has not appeared.

Fix one \(u\) and use a common probability space
\[
S\in\{\pm1\}^s,\qquad Z,V\sim N(0,I_s)
\]
with \(S,Z,V\) independent.  Put
\[
G_u=C_u^{1/2}Z,\qquad
R_u=R_\sigma(S,G_u),\qquad
K_u=\mathbb E R_uR_u^\top,
\]
\[
D_u=T_uK_uT_u,\qquad W_u=D_u^{1/2}V,
\qquad T_u=B_u.
\tag{27}
\]
For a response \(Y=Y(S,Z,V)\), define the moments on this fixed
space
\[
A=\mathbb EYS^\top,\qquad
H_0=\mathbb EYZ^\top,\qquad
M_u=\mathbb EYR_u^\top,\qquad
J_0=\mathbb EYV^\top.
\tag{28}
\]
When the envelope theorem is applied, \(Y\), and hence
\(A,H_0,J_0\), are held fixed; only \(M_u\) retains explicit
\(u\)-dependence.

Since
\[
T_uC_u^{-1/2}=O:=\operatorname{sgn}(B),
\tag{29}
\]
the first paired channel is exactly
\[
\operatorname{tr}(T_u A O H_0^\top).
\tag{30}
\]
Likewise,
\[
\mathcal J=J_0D_u^{-1/2},
\]
so the second paired channel is
\[
\operatorname{tr}
\left(T_uM_uT_uD_u^{-1/2}J_0^\top\right).
\tag{31}
\]
Let \(V_{\sigma,\tau}(u)\) include the one-site regularizer whose
Euler equation is the dithered conditional-score rule.  The
regularizer has no explicit \(u\)-dependence on the fixed whitened
space.  At any differentiable stationary maximizer, the exact
envelope identity is therefore
\[
\boxed{
\begin{aligned}
V_{\sigma,\tau}'(u)=\frac2s\{&
\operatorname{tr}(T_u' A O H_0^\top)\\
&+\operatorname{tr}(T_u'M_uT_uD_u^{-1/2}J_0^\top)\\
&+\operatorname{tr}(T_uM_u'T_uD_u^{-1/2}J_0^\top)\\
&+\operatorname{tr}(T_uM_uT_u'D_u^{-1/2}J_0^\top)\\
&+\operatorname{tr}(T_uM_uT_u
 (D_u^{-1/2})'J_0^\top)\}.
\end{aligned}
}
\tag{32}
\]
Here
\[
T_u'=\frac12\,O\,C_u^{-1/2}Q
\tag{33}
\]
and
\[
M_u'=\mathbb E Y(R_u')^\top
\tag{34}
\]
contains the exact derivative of the regressed first residual.  The
last line of (32), together with the two adjacent \(T_u'\)-terms, is
the fresh-field Onsager derivative.  Thus (32) includes all terms
that were absent from the withdrawn neighbor-\((S,R)\) argument.

An equivalent polar form is useful.  Normalize the residual and put
\[
\widetilde R_u=K_u^{-1/2}R_u,\qquad
\widetilde M_u=\mathbb EY\widetilde R_u^\top,
\]
\[
O_{2,u}=K_u^{1/2}T_uD_u^{-1/2}.
\tag{35}
\]
Then \(O_{2,u}\) is orthogonal and the second channel is
\[
\operatorname{tr}(T_u\widetilde M_uO_{2,u}J_0^\top).
\tag{36}
\]
Its derivative contains
\[
\operatorname{tr}
\left(T_u\widetilde M_u O_{2,u}'J_0^\top\right).
\tag{37}
\]
Writing
\[
\Omega_u=O_{2,u}'O_{2,u}^\top,
\]
one has \(\Omega_u^\top=-\Omega_u\).  The term (37) is therefore a
contraction with a skew-symmetric polar velocity.  It has no fixed
sign.  Nor does the direct term (30), because \(Q\), and hence
\(T_u'\), is generally indefinite.

Stationarity removes \(Y'\), but it does **not** separately remove
either the residual term (34) or the skew-polar term (37).  The exact
unresolved residual in a monotonicity proof is
\[
\boxed{
\mathscr I_u=
\operatorname{tr}(T_u\widetilde M_u' O_{2,u}J_0^\top)
+\operatorname{tr}(T_u\widetilde M_u
 \Omega_uO_{2,u}J_0^\top).
}
\tag{38}
\]
No Schur-complement representation with a definite sign is presently
known for \(\mathscr I_u\).  Formula (38) is an explicit obstruction,
not an omission hidden under the word ``Onsager''.

The finite-fibre search was repeated after correcting a normalization
bug in the initial response: the stored finite-channel residuals are
unnormalized, so their scalar coefficients must be divided by their
standard deviation.  With the correction and a damped conditional
score iteration, the entire nonsingular symmetric sign-fibre census
for \(s=4\) showed no decrease on the coarse grid
\[
u=0,.2,.4,.6,.8,1.
\]
The apparent candidate (26) is therefore withdrawn as a numerical
counterexample; its small negative increments were quadrature and
iteration error.  This still does not prove (24).
