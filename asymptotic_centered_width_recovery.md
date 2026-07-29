# Asymptotic centered-width recovery

## Status

Fix a square fibre size \(s\ge4\).  For a signing \(A\) of order
\(n\), let \(\mathcal L_s(A)\) be the exact compressed-lift class
\[
\mathbf1^\top B_{ij}\mathbf1
=a_{ij}s^{3/2}\qquad(i\ne j).
\tag{1}
\]
The range theorem gives
\[
W(B)\ge s^{3/2}W(A)
\qquad(B\in\mathcal L_s(A)).
\tag{2}
\]
The objective here is the reverse inequality with error
\[
W(B)\le s^{3/2}W(A)+e_s(n),
\qquad e_s(n)=o(n^{3/2}).
\tag{3}
\]

This note gives an exact excess identity and a rigorous probabilistic
recovery criterion.  It reduces (3) to a quantitative two-sided
low-slack traffic bound.  The standard consequences of width
minimality—edge-flip certificates, bounded operator norm, and
mean-square row pseudorthogonality—do not yet prove that traffic bound.

## 1. Exact upper/lower slack identity

Write
\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
P(A)=\max_xH_A(x),\qquad
m(A)=\min_xH_A(x),
\]
\[
W(A)=\frac{P(A)-m(A)}2.
\]
For a microstate
\[
X=(x_1,\ldots,x_n),\qquad x_i\in\{\pm1\}^s,
\]
put
\[
\mu_i(X)=\frac1s\mathbf1^\top x_i.
\tag{4}
\]
Since \(H_A(\mu)\) is multi-affine, its extrema on
\([-1,1]^n\) occur at Boolean corners.  Therefore the two coarse
slacks
\[
u_A(X)=s^{3/2}\bigl(P(A)-H_A(\mu(X))\bigr),
\tag{5}
\]
\[
\ell_A(X)=s^{3/2}\bigl(H_A(\mu(X))-m(A)\bigr)
\tag{6}
\]
are nonnegative.

Choose arbitrary diagonal-fibre signings \(D_i\), and write
\[
D(X)=\sum_iH_{D_i}(x_i),\qquad
d=\sum_iH_{D_i}(\mathbf1).
\tag{7}
\]
For \(B\in\mathcal L_s(A)\), define the exact cross-block residual
\[
R_B(X)=
\sum_{i<j}
x_i^\top
\left(B_{ij}-\frac{a_{ij}}{\sqrt s}J_s\right)x_j.
\tag{8}
\]
Then
\[
H_B(X)=s^{3/2}H_A(\mu(X))+D(X)+R_B(X).
\tag{9}
\]

Define
\[
E_+(A,B)=
\max_X\bigl(R_B(X)-u_A(X)+D(X)-d\bigr),
\tag{10}
\]
\[
E_-(A,B)=
\max_X\bigl(-R_B(X)-\ell_A(X)+d-D(X)\bigr).
\tag{11}
\]
Constant-fibre lifts of a seed maximum and minimum make both quantities
nonnegative.  Equations (9)--(11) give the exact endpoint formulas
\[
P(B)=s^{3/2}P(A)+d+E_+(A,B),
\tag{12}
\]
\[
m(B)=s^{3/2}m(A)+d-E_-(A,B).
\tag{13}
\]
Consequently
\[
\boxed{
W(B)-s^{3/2}W(A)
=\frac{E_+(A,B)+E_-(A,B)}2.
}
\tag{14}
\]

Thus recovery is exactly a two-sided residual-versus-slack problem.
Neither a scalar norm of \(A\) nor the covariance of its rows records
the information in (5)--(6).

## 2. Independent exact-sum microblocks

For every macro edge \(ij\), choose \(B_{ij}\) independently and
uniformly among the \(s\times s\) sign matrices with total sum
\[
a_{ij}s^{3/2}.
\tag{15}
\]
Such matrices exist for square \(s\).  Exchangeability gives
\[
\mathbb E(B_{ij})=\frac{a_{ij}}{\sqrt s}J_s,
\tag{16}
\]
and hence
\[
\mathbb E H_B(X)
=s^{3/2}H_A(\mu(X))+D(X).
\tag{17}
\]

There is also an exact statewise variance formula.  Put \(M=s^2\) and
let \(K\) be the number of \(+1\) entries prescribed by (15).  For
fixed \(x,y\in\{\pm1\}^s\), set
\[
\mu=\frac1s\mathbf1^\top x,\qquad
\nu=\frac1s\mathbf1^\top y.
\]
Sampling the \(K\) positive positions without replacement gives
\[
\boxed{
\operatorname{Var}(x^\top B_{ij}y)
=\kappa_s(1-\mu^2\nu^2),\qquad
\kappa_s=\frac{s^3}{s+1}.
}
\tag{18}
\]
Indeed, the population weights \(w_{\alpha\beta}=x_\alpha y_\beta\)
have mean \(\mu\nu\), and converting the sampled indicators to signs
gives
\[
4\frac{K(M-K)}{M-1}
=\frac{s^4-s^3}{s^2-1}
=\frac{s^3}{s+1}.
\]
The centered block contribution is also bounded in absolute value by
\(2s^2\).

For a complete microstate \(X\), write \(z_i=\mu_i(X)^2\).  Independence
of the macro-edge blocks makes its exact total variance
\[
\boxed{
V_s(X)
=\kappa_s\sum_{i<j}(1-z_iz_j)
=\frac{\kappa_s}{2}
\left[
n(n-1)-\left(\sum_i z_i\right)^2+\sum_i z_i^2
\right].
}
\tag{18a}
\]

Let
\[
N_e=\binom n2
\]
and define the uniform one-sided Bernstein envelope
\[
\Phi_{n,s}(t)=
\exp\left[
-\frac{t^2}
{2\left(N_e\kappa_s+\frac{2s^2t}{3}\right)}
\right],
\qquad t\ge0.
\tag{19}
\]
Then, for every fixed microstate \(X\),
\[
\Pr\{R_B(X)\ge t\},
\quad
\Pr\{-R_B(X)\ge t\}
\le\Phi_{n,s}(t).
\tag{20}
\]
The sharper state-dependent version replaces \(N_e\kappa_s\) in
(19) by \(V_s(X)\).

## 3. Exact traffic criterion

For arbitrary diagonal fibres,
\[
|D(X)-d|\le\delta_s n,\qquad
\delta_s=s(s-1).
\tag{21}
\]
Fix \(e>\delta_sn\), and define
\[
\boxed{
\begin{aligned}
\mathcal T_{A,s}(e)
=\sum_{X\in\{\pm1\}^{ns}}\bigg[
&\Phi_{n,s}\bigl(e-\delta_sn+u_A(X)\bigr)\\
+&\Phi_{n,s}\bigl(e-\delta_sn+\ell_A(X)\bigr)
\bigg].
\end{aligned}}
\tag{22}
\]

### Recovery theorem

If
\[
\boxed{\mathcal T_{A,s}(e)<1,}
\tag{23}
\]
then there is an exact compressed lift \(B\in\mathcal L_s(A)\) such
that
\[
\boxed{
W(B)\le s^{3/2}W(A)+e.
}
\tag{24}
\]

### Proof

If \(E_+(A,B)>e\), some \(X\) satisfies
\[
R_B(X)>
e+u_A(X)+d-D(X)
\ge e-\delta_sn+u_A(X).
\]
Equation (20) bounds this event by the first summand in (22).
The same argument applied to \(-R_B\) and (11) gives the second
summand.  The union bound and (23) produce a realization with
\[
E_+(A,B)\le e,\qquad E_-(A,B)\le e.
\]
Equation (14) proves (24).  Because every block ensemble is finite,
the realization can in principle be found by revealing the blocks one
at a time and applying conditional expectation to the number of
violated endpoint slabs.

The theorem is nonlocal: the chosen value of a block may depend on the
complete two-sided slack profile of the seed.

## 4. Explicit asymptotic rate threshold

For a state whose normalized slack satisfies
\[
u_A(X)=r n^{3/2},
\]
equations (18)--(19) give
\[
\Phi_{n,s}(u_A(X)+o(n^{3/2}))
=
\exp\left[-\left(\frac{r^2}{\kappa_s}+o(1)\right)n\right].
\tag{25}
\]
The same holds for \(\ell_A\).  Thus the explicit Bernstein traffic
threshold is
\[
\boxed{\text{low-slack entropy rate at }r
<\frac{r^2}{\kappa_s}
=\frac{s+1}{s^3}r^2.}
\tag{26}
\]

The exact magnetization-dependent threshold is stronger.  If
\[
\rho(X)=\frac1n\sum_i\mu_i(X)^2,
\]
then (18a) gives, uniformly away from the deterministic case
\(\rho=1\),
\[
\Pr\{R_B(X)\ge rn^{3/2}\}
\le
\exp\left[
-\left(
\frac{r^2}{\kappa_s(1-\rho(X)^2)}+o(1)
\right)n
\right].
\tag{26a}
\]
Thus a shell with normalized slack \(r\) and squared-magnetization
mean \(\rho\) must have entropy rate strictly below
\[
\boxed{
\frac{r^2}{\kappa_s(1-\rho^2)}.
}
\tag{26b}
\]
When \(\rho=1\), every fibre is constant and the residual is
identically zero.

This is a genuine quantitative threshold, not merely a requirement
for some positive endpoint rate.  Among nontrivial exact square fibre
sizes \(s\ge4\),
\[
\max_s\frac1{\kappa_s}
=\frac1{\kappa_4}
=\frac5{64}.
\tag{26c}
\]
Larger \(s\) makes the coefficient smaller.  Choosing the fibre size
therefore cannot turn an arbitrarily weak positive LDP rate into the
traffic estimate.

A convenient sufficient LDP formulation is the following.  Suppose
there is a sequence \(\varepsilon_n\downarrow0\) and
\[
e_n=\varepsilon_n n^{3/2},\qquad
e_n\gg n,
\tag{27}
\]
such that the upper- and lower-slack shell counts, uniformly over a
mesh \(r\in\varepsilon_n\mathbb Z_{\ge0}\), satisfy
\[
\log\#\left\{
X:\ rn^{3/2}\le u_A(X)<(r+\varepsilon_n)n^{3/2}
\right\}
\le
\left(
\frac{(r+\varepsilon_n/2)^2}{\kappa_s}
-\omega_n
\right)n,
\tag{28}
\]
and analogously for \(\ell_A\), where
\[
\omega_n n\to\infty
\]
uniformly over occupied shells.  With the harmless adjustment of the
mesh constants needed near \(r=0\), (28) makes (22) tend to zero and
therefore yields
\[
W(B)\le s^{3/2}W(A)+o(n^{3/2}).
\tag{29}
\]

In particular, the zero-slack layers must be subexponential strongly
enough that
\[
\log\#\{X:u_A(X)=0\},
\quad
\log\#\{X:\ell_A(X)=0\}
=o(e_n^2/n^2).
\tag{30}
\]
Equations (26)--(30) are sufficient conditions, not facts currently
known for centered-width minimizers.

## 5. What width minimality currently forces

Let \(A\) be a centered-width minimizer.  Flipping one edge cannot
decrease its range.  The exact two-sided edge certificate is:
for every edge \(ij\), at least one of the following holds:

1. there is an upper state \(x\) with
   \[
   H_A(x)\ge P(A)-2,\qquad
   a_{ij}x_ix_j=-1;
   \]
2. there is a lower state \(y\) with
   \[
   H_A(y)\le m(A)+2,\qquad
   a_{ij}y_iy_j=+1.
   \]

Otherwise flipping \(a_{ij}\) lowers the upper endpoint by two and
raises the lower endpoint by two.

There are also exact spectral consequences.  For Boolean \(u,v\),
split the coordinates according to \(u_i=v_i\) or \(u_i=-v_i\).
The cross terms cancel, and principal-submatrix monotonicity gives
\[
\|A\|_{\infty\to1}\le4W(A).
\tag{31}
\]
Together with
\[
\|A\|_{\mathrm{op}}^2\le\|A\|_{\infty\to1},
\]
this yields, on every \(W(A)=O(n^{3/2})\) sequence,
\[
\|A\|_{\mathrm{op}}=O(n^{3/4}).
\tag{32}
\]
For normalized row correlations
\[
q_{ij}=\frac{(A^2)_{ij}}{n-1},
\]
\[
\frac1{n^2}\sum_{i\ne j}q_{ij}^2
=O(n^{-1/2}).
\tag{33}
\]

The gap is now precise:

- the edge certificate supplies *coverage* by near-endpoint witnesses,
  not an upper bound on their number;
- (32)--(33) control second-order geometry, while (22) depends on the
  complete large-deviation profile of the multi-affine slacks;
- universal Hamming clouds show that raw near-endpoint entropy is
  inevitably exponential at proportional slack.

Thus none of (31)--(33), separately or together, proves (26).

## 6. Generic residuals and the fixed-gap warning

Without a slack LDP, the independent-block construction gives only
\[
R_B(X)=O_s(n^{3/2})
\]
uniformly over all \(2^{sn}\) microstates, hence only
\[
W(B)\le s^{3/2}W(A)+O_s(n^{3/2}).
\tag{34}
\]
This is a fixed normalized loss, not amplification.

The obstruction is not merely the number of states.  Every sign block
with the prescribed coarse mean has compulsory microscopic
Frobenius/ANOVA mass.  Independent block choices turn that mass into a
finite-alphabet vector-spin process at scale \(n^{3/2}\).  To reduce
(34) to \(o(n^{3/2})\), a globally correlated choice must make this
process fit simultaneously under both slack roofs in (10)--(11).
The exact order-five Pythagorean counterexample shows that such
anti-alignment is possible in finite systems, so a generic residual
lower bound alone cannot rule it out.

## 7. Surviving global replacement target

The useful next statement is an entropy-versus-replacement dichotomy.
For a width minimizer \(A\), prove one of:

1. its two-sided slack shells obey the traffic threshold (26), giving
   (29) by the recovery theorem; or
2. excessive low-slack traffic has enough coherent covariance or
   majority structure to choose a global family of edge-block
   replacements that lowers \(W(A)\), contradicting minimality.

The replacement must be global.  Local edge-flip certificates permit
large near-endpoint layers and cannot distinguish a coherent
low-slack cloud from unrelated witnesses.  A successful argument
needs a quantitative statement about the signed covariance
\[
\mathbb E_{\mathcal L_+}(x_ix_j)
-\mathbb E_{\mathcal L_-}(x_ix_j)
\tag{35}
\]
of upper and lower low-slack layers, or an equivalent separating
functional on their cut vectors.

At present, (35) is the exact missing inequality.  The recovery side
is closed by (22)--(24); the global-minimality side is not.

## 8. Exact global-replacement audit

The replacement side has an exact finite formulation.  For a pair of
Boolean states \(p=(x,y)\), put
\[
s_e(p)=\frac{x_ix_j-y_iy_j}{2}\in\{-1,0,1\},
\qquad
\operatorname{score}_A(p)=\sum_e a_es_e(p),
\]
and
\[
\delta_A(p)=W(A)-\operatorname{score}_A(p).
\tag{36}
\]
Let \(T\) be any set of edges and replace \(a_T\) by
\(\beta\in\{\pm1\}^T\).  Direct substitution gives
\[
\boxed{
W(A^{T\to\beta})-W(A)
=
\max_p\left\{
(\beta-a_T)\mathbin{\cdot}s_T(p)-\delta_A(p)
\right\}.
}
\tag{37}
\]
Thus only profiles with slack comparable to the possible gain on
\(T\) can obstruct a replacement.

Equation (37) also gives a randomized sufficient criterion.  If a
distribution \(\Pi\) on replacements and some \(\gamma>0\) satisfy
\[
\sum_p
\Pr_{\beta\sim\Pi}
\left\{
(\beta-a_T)\mathbin{\cdot}s_T(p)>
\delta_A(p)-\gamma
\right\}<1,
\tag{38}
\]
then some replacement lowers the width by at least \(\gamma\).
Conversely, at a global width minimizer, for every \(\beta\) there is
a profile \(p\) with
\[
\delta_A(p)\le(\beta-a_T)\mathbin{\cdot}s_T(p)\le2|T|.
\tag{39}
\]
Hence the projected profiles with slack at most \(2|T|\) cover every
replacement.  This is a precise warning: global minimality naturally
forces *lower* coverage by near-active profiles, whereas recovery
needs an *upper* traffic bound.

There is a useful approximate dual form of the same obstruction.
Suppose \(T\) has \(m\) edges incident to \(v\) vertices, and set
\[
\eta_T=\sqrt{2m(2v\log2+1)}.
\tag{40}
\]
Continuous relaxation of \(\beta\), followed by independent sign
rounding, yields a probability measure \(\mu\) on pair profiles such
that
\[
\mathbb E_\mu\operatorname{score}_A(p)\ge W(A)-\eta_T
\tag{41}
\]
and
\[
\boxed{
2\sum_{e\in T}
\left(a_e\mathbb E_\mu s_e(p)\right)_+
\le\eta_T.
}
\tag{42}
\]
For completeness, minimize the convex function
\[
f(z)=\max_p\{(z-a_T)\cdot s_T(p)-\delta_A(p)\}
\]
over \(z\in[-1,1]^T\).  Rounding a minimizer \(z_*\) and union-bounding
over the at most \(2^{2v}\) pair profiles shows \(f(z_*)\ge-\eta_T\).
A subgradient measure supported on the active profiles at \(z_*\)
then gives (41).  The box optimality conditions and comparison with
\(z=a_T\) give (42).

For an induced \(k\)-vertex block, (40) is \(O(k^{3/2})\); for a
\(k\times n\) cross block it is \(O(n\sqrt k)\).  Any successful
replacement theorem must therefore beat a near-active distribution
whose mean width gradient is almost coordinatewise anti-aligned with
the original signing.

There is also an exact cut-factorization of the slack traffic.  Switch
a maximizing state to the all-plus state and define
\[
C(S)=\sum_{ij\in\partial S}a_{ij}.
\]
Then \(0\le C(S)\le W(A)\), and
\[
H_A(x^S)=P(A)-2C(S).
\]
For \(p=(x^S,x^T)\),
\[
\operatorname{score}_A(p)=C(T)-C(S),
\qquad
\delta_A(p)=C(S)+[W(A)-C(T)].
\tag{43}
\]
Consequently the pair-profile Laplace sum factors exactly, up to the
harmless two-to-one spin/cut multiplicities:
\[
\boxed{
\sum_{S,T}e^{-\lambda\delta_A(S,T)}
=
\left(\sum_Se^{-\lambda C(S)}\right)
\left(\sum_Te^{-\lambda(W(A)-C(T))}\right).
}
\tag{44}
\]
Writing \(d=P(A)-W(A)\), uniform random cuts satisfy
\[
\mathbb E C=\frac{W+d}{2},\qquad
\mathbb E(W-C)=\frac{W-d}{2},\qquad
\operatorname{Var}(C)=\frac14\binom n2
\tag{45}
\]
exactly.  Positive midpoint \(d\) helps the low-cut factor, but
variance alone gives only polynomial endpoint rarity.  The missing
input is a cone-specific exponential tail, or an equivalent
entropy-versus-replacement theorem strong enough to beat (26b).

## 9. Verdict

The fixed-\(s\) centered-width problem has been reduced to an explicit
large-deviation statement:
\[
\boxed{
\text{quadratic slack traffic below }r^2/\kappa_s
\Longrightarrow
\text{scale-preserving recovery}.
}
\]
No unconditional \(e_s(n)=o(n^{3/2})\) bound is proved.  Conversely,
no fixed-gap obstruction is proved for globally correlated
microblocks.  The remaining route is the global
entropy-versus-replacement dichotomy in Sections 7--8, not another
local block gadget.  Equations (39)--(42) show why minimality alone
does not automatically provide it.

Related sidecars:

- `centered_width_amplification_reboot.md`
- `spectral_filter_microblock_recovery.md`
- `pythagorean_centered_width_block.md`
- `entropy_energy_dichotomy.md`
