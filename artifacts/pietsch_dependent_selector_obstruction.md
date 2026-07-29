# Pietsch-weighted vertex selectors: exact variance tradeoff and a
# first-order obstruction

## 1. Factorized selector calculus

Let

\[
B=\|A\|_{\infty\to1}
\]

and use Grothendieck--Pietsch factorization in the form

\[
A=D_\mu^{1/2}TD_\nu^{1/2},
\qquad
\|T\|_{\mathrm{op}}\le K_GB,
\tag{1}
\]

where \(\mu,\nu\) are probability vectors.

Let \(\delta_i\) be independent Bernoulli selectors with means \(p_i\),
put \(\xi_i=\delta_i-p_i\), and let
\[
r_i=p_i(1-p_i).
\]
For a fixed Boolean vector \(x\), the selected quadratic form has the
exact Hoeffding/ANOVA expansion

\[
q_\delta(x)-q_p(x)
=
\sum_i \xi_i h_i(x)
+
\sum_{i<j}a_{ij}x_ix_j\xi_i\xi_j,
\tag{2}
\]

where

\[
q_\delta(x)=\sum_{i<j}a_{ij}x_ix_j\delta_i\delta_j,\qquad
q_p(x)=\sum_{i<j}a_{ij}x_ix_jp_ip_j,
\]

and

\[
h_i(x)=x_i\sum_{j\ne i}a_{ij}x_jp_j.
\]

The two terms in (2) are orthogonal in \(L_2\).  Factorization gives

\[
\boxed{
\sum_i r_i h_i(x)^2
\le
(K_GB)^2
\left(\max_i r_i\mu_i\right)
\left(\sum_j\nu_jp_j^2\right).
}
\tag{3}
\]

Moreover,

\[
\boxed{
\|D_{\sqrt r}AD_{\sqrt r}\|_{\mathrm{op}}
\le
K_GB
\sqrt{\left(\max_i r_i\mu_i\right)
      \left(\max_j r_j\nu_j\right)}
}
\tag{4}
\]

and, because \(|a_{ij}|=1\),

\[
\|D_{\sqrt r}AD_{\sqrt r}\|_F^2
=\sum_{i\ne j}r_ir_j
\le\left(\sum_i r_i\right)^2.
\tag{5}
\]

Equations (3)--(5) are the exact bias/variance inputs for a
Pietsch-adapted Hanson--Wright argument.  Taking \(p_i\) near zero or one
on vertices with large \(\mu_i,\nu_i\) suppresses their variance.

## 2. Bias from nonuniform probabilities

Write \(p=\alpha{\bf1}+e\), with \(\|e\|_\infty\le\eta\).  Box
polarization gives

\[
\boxed{
Q(D_pAD_p)
\le
\alpha^2Q(A)+B(2\alpha\eta+\eta^2).
}
\tag{6}
\]

Thus eliminating a high-weight Pietsch hub by moving \(p_i\) from
\(\alpha\) to zero uses \(\eta\asymp\alpha\) and can already cost a
leading-order fraction of \(Q(A)\).  In the regular worst case
\(\mu_i=\nu_i=1/n\), uniform \(p_i=\alpha\) is the natural bias minimizer,
and (3) has scale

\[
(K_GB)^2\frac{\alpha^3(1-\alpha)}n.
\tag{7}
\]

After maximizing over exponentially many spin vectors this becomes a
leading \(n^{3/2}\) term.  Its square-root dependence on \(1-\alpha\)
does not fit the scale-transfer margin

\[
(\alpha^{3/2}-\alpha^2)Q(A),
\]

which is linear in \(1-\alpha\) near one.

## 3. No-go theorem for cancelling the first-order field

The proposed stronger strategy was to use dependent fixed-cardinality
rounding so that the first-order matrix

\[
L_e=D_eAD_p+D_pAD_e
\tag{8}
\]

has Boolean quadratic norm \(o(n^{3/2})\).

This is impossible already at density \(\alpha=1/2\), for every signing
and every selector.

Let \(S\) be any half-set, \(T=S^c\), \(p_i=1/2\), and
\(e_i=1/2\) on \(S\), \(e_i=-1/2\) on \(T\).  Directly from the entries,

\[
\boxed{
L_e=\frac12\left(A_{S,S}\oplus(-A_{T,T})\right).
}
\tag{9}
\]

Therefore

\[
\boxed{
Q(L_e)
\ge
\frac12\max\{Q(A_{S,S}),Q(A_{T,T})\}
\ge
\frac12\,Q_{\min}(n/2)
=\Theta(n^{3/2}).
}
\tag{10}
\]

The last step uses any universal \(c\,m^{3/2}\) lower bound for signings
of order \(m\).  This obstruction is independent of the Pietsch weights
and survives perfect spectral regularization.

Consequently no vector discrepancy theorem, Gram--Schmidt walk,
interlacing-polynomial method, swap rounding, or matrix potential can
satisfy the requested \(o(n^{3/2})\) first-order target.  The obstruction
is in the target itself, not in the rounding algorithm.

## 4. Why first and second order must be controlled jointly

At \(\alpha=1/2\), put \(\sigma=2\delta-1\).  Then

\[
D_\delta AD_\delta
=\frac14
\left(A+D_\sigma A+AD_\sigma+D_\sigma AD_\sigma\right).
\tag{11}
\]

The second-order term cancels the entire leading block from (9) on
\(T\), and reinforces it on \(S\).  Indeed the right side of (11) is
exactly \(A_{S,S}\oplus0_T\).

Thus separating the ANOVA orders by triangle inequalities necessarily
loses a leading term.  A viable dependent-rounding theorem must prove
cancellation between first and second order in the Boolean supremum.
At half density the exact sufficient condition is simply

\[
Q\!\left(
A+D_\sigma A+AD_\sigma+D_\sigma AD_\sigma
\right)
\le\sqrt2\,Q(A),
\tag{12}
\]

for some balanced \(\sigma\).  But the matrix on the left is exactly
\(4(A_{S,S}\oplus0)\), so (12) is equivalent to the desired proportional
restriction inequality itself.  The ANOVA decomposition does not make
that core problem easier unless a new joint-cancellation invariant is
introduced.
