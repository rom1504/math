# Wave 30 route 7: matched parent common mode without a Taylor remainder

## Status

- **Verified, all deletion sizes:** the matched parent likelihood is the
  normalized exponential tilt by one explicit harmonic extension cost
  \(\Phi_\beta\). Its entropy is exactly an integral of tilted variances of
  \(\Phi_\beta\). This is a uniform sufficient reduction, not a fixed-\(A\)
  expansion.
- **Verified, one deletion:** for every signing,
  \[
  \operatorname{Ent}_{\nu_\beta}(f_\beta)
  \le 2\beta\sqrt{\mathbb E_{\nu_\beta}R_2/n}
  \le 2\beta\|A\|_{\rm op}.
  \]
  For an exact minimizer and
  \(\beta=\Theta(n^{-1/2+c})\), this is
  \(O(n^{1/4+c})\), and therefore proves the first estimate in (10.800)
  at the one-deletion landing when \(c\le1/12\).
- **Verified obstruction:** centered \(R_2\) is only the first common mode.
  The exact minimizer \(A_6\) has \(A_6^2=5I\), so \(R_2\equiv30\), but
  its matched parent likelihood is nonconstant and its parent entropy is
  positive.
- **Open:** the fixed-density case \(n-m=\Theta(n)\). The one-deletion
  estimate cannot be iterated to prove (10.800) on an active ratio window.

All finite normalizations and audits are checked by
`tmp/parent_common_mode_r30.py`.

## 1. Exact all-\(k\) harmonic-likelihood reduction

Let \(k=n-m\), let \(y=d[S]\), and write the energy of a full extension
as

\[
\langle A,d\rangle=c_S(y)+H_{S,y}(z),
\qquad
K_{\beta,S}(y)=\sum_{z\in\{\pm1\}^k}e^{\beta H_{S,y}(z)}.
\tag{R30.1}
\]

The uniform outside completion has
\(\mathbb E_zH_{S,y}(z)=0\): every external edge monomial contains at
least one centered outside spin. Hence Jensen gives

\[
F_{\beta,S}(y):=\log\frac{K_{\beta,S}(y)}{2^k}\ge0,\qquad
u_{\beta,S}(y):=e^{-F_{\beta,S}(y)}=\frac{2^k}{K_{\beta,S}(y)}\le1.
\tag{R30.2}
\]

Define the common harmonic score and its cost

\[
U_\beta(d)=\mathbb E_{S\sim U_m}u_{\beta,S}(d[S]),
\qquad
\Phi_\beta(d)=-\log U_\beta(d)\ge0.
\tag{R30.3}
\]

At \(\gamma=\beta\), (10.884) gives

\[
f_\beta(d)=\frac{C_\beta}{2^k}U_\beta(d).
\]

Since \(\mathbb E_{\nu_\beta}f_\beta=1\), all constants cancel exactly:

\[
\boxed{
f_\beta(d)=
\frac{e^{-\Phi_\beta(d)}}
{\mathbb E_{\nu_\beta}e^{-\Phi_\beta}},
\qquad
\mathbb E_{\nu_\beta}U_\beta
=\frac{2^kZ_0(\beta)}{Z_A(\beta)}.
}
\tag{R30.4}
\]

This was checked directly also at \(k=3\), not inferred from the
one-deletion formula.

Put

\[
\psi(t)=\log\mathbb E_{\nu_\beta}e^{-t\Phi_\beta},\qquad
\nu_{\beta,t}(d)=e^{-t\Phi_\beta(d)-\psi(t)}\nu_\beta(d).
\]

Then \(\psi''(t)=\operatorname{Var}_{\nu_{\beta,t}}(\Phi_\beta)\),
\(\nu_{\beta,1}=f_\beta\nu_\beta\), and two integrations give the exact,
nonperturbative identity

\[
\boxed{
\operatorname{Ent}_{\nu_\beta}(f_\beta)
=\psi'(1)-\psi(1)
=\int_0^1t\,
\operatorname{Var}_{\nu_{\beta,t}}(\Phi_\beta)\,dt.
}
\tag{R30.5}
\]

Thus a genuine all-\(k\) sufficient theorem is

\[
\boxed{
\int_0^1t\,
\operatorname{Var}_{\nu_{\beta,t}}(\Phi_\beta)\,dt
=O(n^{1/2-2c}),
}
\tag{R30.6}
\]

uniformly for target-specific exact minimizers and every target order in
one fixed ratio window. A stronger but sometimes simpler sufficient input
is

\[
\log\frac{Z_A(\beta)}{2^kZ_0(\beta)}
=-\log\mathbb E_{\nu_\beta}U_\beta
=O(n^{1/2-2c}),
\tag{R30.7}
\]

because

\[
\operatorname{Ent}_{\nu_\beta}(f_\beta)
=-\mathbb E_{\nu_{\beta,1}}\Phi_\beta
-\log\mathbb E_{\nu_\beta}e^{-\Phi_\beta}
\le-\log\mathbb E_{\nu_\beta}U_\beta.
\tag{R30.8}
\]

The Rényi bound

\[
\operatorname{Ent}_{\nu_\beta}(f_\beta)
\le
\log\frac{\mathbb E_{\nu_\beta}U_\beta^2}
{(\mathbb E_{\nu_\beta}U_\beta)^2}
\tag{R30.9}
\]

is another exact sufficient reduction. No claim is made that the stronger
free-energy gap (R30.7) already has the target scale.

## 2. Exact one-deletion theorem

For \(m=n-1\), put

\[
h_i(d)=\sum_{j\ne i}A_{ij}d_{ij},\qquad
R_2(d)=\sum_i h_i(d)^2.
\]

There are two extensions of the deleted child and their external energies
are \(\pm2h_i\). Therefore

\[
\boxed{
U_\beta(d)=\frac1n\sum_{i=1}^n\operatorname{sech}(2\beta h_i(d)).
}
\tag{R30.10}
\]

The elementary inequalities

\[
\operatorname{sech}t\ge e^{-|t|},
\qquad
\operatorname{sech}t\ge e^{-t^2/2}
\tag{R30.11}
\]

and Jensen give, pointwise,

\[
U_\beta(d)\ge
\exp\!\left[-2\beta\sqrt{\frac{R_2(d)}n}\right],
\qquad
U_\beta(d)\ge
\exp\!\left[-2\beta^2\frac{R_2(d)}n\right].
\tag{R30.12}
\]

Apply Jensen once more under \(\nu_\beta\) and use (R30.8):

\[
\boxed{
\operatorname{Ent}_{\nu_\beta}(f_\beta)
\le
\min\left\{
2\beta\sqrt{\frac{\mathbb E_{\nu_\beta}R_2}{n}},
\frac{2\beta^2}{n}\mathbb E_{\nu_\beta}R_2
\right\}.
}
\tag{R30.13}
\]

Since \(R_2(d)=x^{\mathsf T}A^2x\le n\|A\|_{\rm op}^2\),

\[
\boxed{
\operatorname{Ent}_{\nu_\beta}(f_\beta)
\le2\beta\|A\|_{\rm op}.
}
\tag{R30.14}
\]

For an exact minimizer, the proved bound
\(\|A\|_{\rm op}^2\le2q_n=O(n^{3/2})\) yields

\[
\operatorname{Ent}_{\nu_\beta}(f_\beta)
=O(n^{1/4+c})
=O(n^{1/2-2c})\quad\text{when }c\le1/12.
\tag{R30.15}
\]

This is uniform in \(n\) and uses no Taylor remainder. It only handles
\(m=n-1\). In particular, it does not establish the fixed-density
requirement in (10.800), and summing this cost over \(\Theta(n)\)
one-vertex deletions is far too expensive.

A sharper data-dependent successor is explicit: the first bound in
(R30.13) reaches the target whenever
\(\mathbb E_{\nu_\beta}R_2=O(n^{3-6c})\). Proving such a Gibbs row-square
bound over a fixed-density endpoint is still open.

## 3. Exact row-square falsifier and the higher-mode wall

Let \(P_\beta=f_\beta\nu_\beta\), and let \(R=R_2(D)\). The KL chain rule
gives

\[
\boxed{
D(P_\beta\Vert\nu_\beta)
=D((P_\beta)_R\Vert(\nu_\beta)_R)
+\mathbb E_{r\sim(P_\beta)_R}
D(P_\beta(\cdot\mid R=r)\Vert\nu_\beta(\cdot\mid R=r)).
}
\tag{R30.16}
\]

Thus an exact non-Taylor row-square falsification criterion is available.
For any row event \(G\), put \(p=P_\beta(G)\), \(q=\nu_\beta(G)\). Data
processing gives

\[
\boxed{
\operatorname{Ent}_{\nu_\beta}(f_\beta)
\ge {\rm kl}(p\Vert q)
\ge p\log(1/q)-\log2.
}
\tag{R30.17}
\]

Consequently, if along an asymptotic exact-minimizer family a row event has
\(P_\beta(G)\ge\delta>0\) but
\(-\log\nu_\beta(G)=\omega(n^{1/2-2c})\), then the first estimate in
(10.800) is falsified. Equivalently, for
\(\Delta_R=\mathbb E_{P_\beta}R-\mathbb E_{\nu_\beta}R\),

\[
\boxed{
\operatorname{Ent}_{\nu_\beta}(f_\beta)
\ge\sup_{\lambda\in\mathbb R}
\left\{\lambda\Delta_R-
\log\mathbb E_{\nu_\beta}
e^{\lambda(R-\mathbb E_{\nu_\beta}R)}\right\}.
}
\tag{R30.18}
\]

These are genuine sufficient falsifiers, not assertions that such a
minimizer family exists.

Centered \(R_2\) cannot by itself prove the upper bound. The exact
order-six minimizer satisfies

\[
A_6^2=5I,\qquad R_2(d)=30\quad\text{for every cut }d.
\tag{R30.19}
\]

Hence the first term of (R30.16) is identically zero. Nevertheless, at one
deletion its harmonic score has two values,

\[
U_{\rm ground}
=\frac{5\operatorname{sech}(2\beta)+\operatorname{sech}(10\beta)}6,
\qquad
U_{\rm excited}
=\frac{\operatorname{sech}(2\beta)+\operatorname{sech}(6\beta)}2.
\tag{R30.20}
\]

They are unequal (for example \(0.5422911084\ldots\) and
\(0.3736911005\ldots\) at \(\beta=1/2\)). Thus \(f_\beta\) is
nonconstant although centered \(R_2\) vanishes exactly. Direct enumeration
gives

\[
\operatorname{Ent}_{\nu_{1/2}}(f_{1/2})
=0.008845026487227\ldots>0.
\tag{R30.21}
\]

The high-temperature \(R_2\) term is therefore only the first common mode;
higher field-profile modes survive inside a fixed row-square level. This is
a finite mechanism wall, not an asymptotic falsifier of (10.800).

## 4. Audits of the named walls

- **\(A_4\).** At the matched point \(\beta=\gamma=1/2\),
  \(\|A_4\|_{\rm op}=\sqrt5\), the exact parent entropy is
  \(0.012957183462626\ldots\), and (R30.14) gives the valid upper bound
  \(\sqrt5\). The earlier \(A_4\) coefficient failure used mismatched
  \((\beta,\gamma)=(1/2,1/10)\); it is not being silently imported into
  the matched theorem.

- **\(J-I\).** This actual complete signing is not an exact minimizer.
  Its norm is \(n-1\), so (R30.14) gives only
  \(O(\beta n)=O(n^{1/2+c})\), not the target scale. Exact enumeration of
  the fixed-\(n\) leading coefficients gives
  \[
  4(n-1)(n-2)^2/n,\qquad4(n-2),
  \]
  whose ratio is \((n-1)(n-2)/n\), reproducing (10.888). This audit uses
  the exact coefficients only; no fixed-\(A\) Taylor remainder is treated
  as uniform at growing \(n\).

## 5. Frontier

The all-\(k\) object to control is now the harmonic extension cost
\(\Phi_\beta\), not merely its quadratic tangent. Proving (R30.6), or a
target-scale Rényi bound (R30.9), for selected exact minimizers at fixed
density would establish the parent half of (10.800). The one-deletion
theorem supplies a real endpoint bound and exposes the precise
\(c\le1/12\) range, but does not close a convergence window. The
\(A_6\) wall requires any fixed-density successor to control higher
within-row modes as well as centered \(R_2\).
