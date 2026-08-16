# Referee report: total-variance bilinear/profile rounding

Date: 2026-08-16.

## Verdict

The rounding claims (1.4)--(1.7), the normalization in (1.12), and the
uniform one-profile consequence are correct.  In fact, (1.12) gives more:
at the same order it makes every fixed (k)-profile, and hence the standard
full action distance (d_M), tend to zero between (T_A) and (T_W).  No
maximum-row variance condition, exceptional-row deletion, or operator-norm
bound on the rounded (A) is needed for this comparison.

Conditional on Lemma L, the convergence implication (1.11) is also correct.
It uses the operator bound only for the weighted-to-target continuity step,
then uses the direct scalar defect for (A-W).  It does not obtain an
operator bound on (T_A), nor does it need one.  The final transfer from a
ratio-dense recovery order (m) to an arbitrary order (N) still uses a
principal submatrix; this transfers the objective, not action profiles.

The profile conclusion is not a genuinely new consequence relative to the
proof in `sign_near_weighted_recovery.md`: its one-spin Bernstein event,
together with real polarization, already implies the same
(L^\infty\!\to L^1) error.  The direct two-spin proof makes this consequence
explicit and gives a better direct estimate for (B(E)), but when reused for
the scalar objective it has worse finite constants than the one-spin proof.
There is no loss at the normalized (m^{3/2}) leading scale under
(V=o(m^2)), because either error is (o(m^{3/2})).

## 1. Check of (1.4)--(1.7)

For fixed (x,y\in\{\pm1\}^m), write

\[
x^{\mathsf T}Ey
=\sum_{i<j}c_{ij}(A_{ij}-w_{ij}),
\qquad c_{ij}=x_i y_j+x_j y_i\in\{-2,0,2\}.
\]

Each summand is centered, is bounded in absolute value by (4), and has
variance

\[
c_{ij}^2(1-w_{ij}^2)\le 4(1-w_{ij}^2).
\]

Thus its total variance is at most (4V(W)).  Bernstein with increment bound
(4), followed by the union bound over (4^m) ordered sign pairs, gives
exactly

\[
\Pr\{B(E)\ge t\}
\le 2\,4^m\exp\!\left[-\frac{t^2}{2(4V+4t/3)}\right].
\]

There is no missing factor from the symmetric double occurrence of an edge:
that occurrence is precisely the coefficient (c_{ij}) and accounts for the
factor (4) in both the variance bound and the increment bound.

Taking natural logarithms and

\[
a_m=2m\log 2+\log4,
\qquad t_m=\sqrt{8Va_m}+\frac83a_m
\]

is a valid Bernstein inversion with variance proxy (4V) and increment
bound (4).  Indeed

\[
\frac{t_m^2}{2(4V+4t_m/3)}\ge a_m,
\]

so the displayed probability is at most
(2\,4^m e^{-a_m}=1/2).  Hence a supported outcome with (B(E)\le t_m)
exists.

Since (H_E(x)=\tfrac12x^{\mathsf T}Ex),

\[
|Q(A)-Q(W)|\le Q(E)
=\frac12\max_x|x^{\mathsf T}Ex|
\le\frac12B(E).
\]

Thus (1.6) is correct (and has the two-sided strengthening displayed above).
Finally (a_m=\Theta(m)), so (V=o(m^2)) gives

\[
t_m=O(\sqrt{mV}+m)=o(m^{3/2}),
\]

which proves every assertion in (1.7).

## 2. Norm and action normalization

The unnormalized identity is

\[
\sup_{\|f\|_\infty\le1}\|Ef\|_{\ell^1}=B(E).
\]

Indeed,

\[
\|Ef\|_{\ell^1}
=\max_{y\in\{\pm1\}^m}y^{\mathsf T}Ef,
\]

and maximizing the resulting linear form over the (f)-cube puts (f) at
a sign vertex.  The absolute value in the definition of (B) is immaterial,
because one sign vector can be flipped.  The note should preferably write
(\ell^1) here to distinguish counting norm from normalized (L^1).

With uniform vertex measure, the conventions are

\[
(T_Ef)_i=\frac1{\sqrt m}(Ef)_i,
\qquad
\langle f,T_Ef\rangle
=\frac1{m^{3/2}}f^{\mathsf T}Ef.
\]

Consequently there is no missing (m) or (2):

\[
\|T_E\|_{L^\infty\to L^1}
=\sup_{\|f\|_\infty\le1}
  \frac1m\sum_i\left|\frac{(Ef)_i}{\sqrt m}\right|
=\frac{B(E)}{m^{3/2}}=:r_m.
\]

Likewise, hollow symmetry and separate affinity give

\[
\Phi(T_W)=\frac{2Q(W)}{m^{3/2}},
\]

so (1.9) is normalized correctly.  Directly,

\[
|\Phi(T_A)-\Phi(T_W)|\le r_m,
\]

because

\(
|\langle f,T_Ef\rangle|
\le\|f\|_\infty\|T_Ef\|_{L^1}.
\)

## 3. Profile consequence, including all fixed (k)

For one test function, couple the two laws at the same uniform vertex.  Their
first coordinates agree and the expected absolute output discrepancy is at
most (r_m).  Markov at threshold (sqrt{r_m}), followed by the coupling
characterization of Levy--Prokhorov distance, gives

\[
d_{LP}\bigl(\mathcal L(f,T_Af),\mathcal L(f,T_Wf)\bigr)
\le \min\{1,\sqrt{r_m}\}.
\]

This is uniform in (f), proves (1.12), and pairs the same (f) in both
directions, so it bounds the Hausdorff as well as either directed one-profile
distance.

The same argument applies to a (k)-tuple (f_1,\ldots,f_k).  Under the
identity coupling the input coordinates again agree, while

\[
\mathbb E\sum_{j=1}^k
 |(T_A-T_W)f_j|\le kr_m.
\]

For any of the usual truncated product metrics on \(\mathbb R^{2k}\), Markov
therefore yields

\[
d_H\bigl(\mathcal S_k(T_A),\mathcal S_k(T_W)\bigr)
\le \min\{1,\sqrt{kr_m}\}.                              \tag{R.1}
\]

In the standard action metric
(d_M=\sum_{k\ge1}2^{-k}d_H(\mathcal S_k(\cdot),
\mathcal S_k(\cdot))), (R.1) gives

\[
d_M(T_A,T_W)
\le \left(\sum_{k\ge1}2^{-k}\sqrt{k}\right)\sqrt{r_m}
=o(1).
\]

No output moment or operator-norm bound is used: rare large outputs are
handled by the displayed Markov tail estimate.  This upgrades only the
rounded-versus-weighted comparison.  Lemma L assumes merely directed
one-profile convergence of (T_W) to (T_\ell), so it implies the same
directed one-profile convergence for (T_A), but not full action convergence
to (T_\ell).  Full convergence to the target would follow if the weighted
input itself satisfied (d_M(T_W,T_\ell)\to0).

## 4. Exact convergence implication and its limits

For fixed (ell), Lemma L and directed continuity first give

\[
\Phi(T_{W_m})
\le\Phi(T_\ell)+5D_m\sqrt{\delta_m}+\delta_m
=\Phi(T_\ell)+o(1).
\]

The total-variance rounding then gives, at the same orders,

\[
\frac{Q(A_m)}{m^{3/2}}
\le \frac12\Phi(T_\ell)+o(1)
\le L+\frac12\eta_\ell+o(1).
\]

Upward ratio density and objective monotonicity under a principal submatrix
therefore yield

\[
\limsup_N\frac{M_N}{N^{3/2}}
\le L+\frac12\eta_\ell,
\]

and then convergence after (ell\to\infty).  This proof needs neither
exceptional-row pruning nor an operator bound for (T_A).  It does still use
the (D_m\sqrt{\delta_m}\to0) hypothesis for (T_W\to T_\ell), and it uses
principal deletion solely to cover all orders.  Without further control,
that last deletion should not be claimed to preserve profiles.

Also, weak profile closeness by itself is not a substitute for operator
control when passing an unbounded quadratic observable.  Here the objective
passes because the stronger (L^\infty\!\to L^1) estimate directly bounds
the quadratic form.  Total variance alone only makes (A) close to its given
(W); Lemma L's externally certified weighted recovery remains the
substantive missing input.

## 5. Comparison with the existing one-spin rounding

The direct two-spin argument is valid but is not needed to deduce profile
preservation.  For symmetric hollow (E), real polarization gives

\[
B(E)\le2q(E),
\qquad q(E)=\max_{x\in\{\pm1\}^m}|x^{\mathsf T}Ex|.
\]

Indeed, for sign vectors (x,y), put
(u=(x+y)/2) and (v=(x-y)/2).  Then

\[
x^{\mathsf T}Ey=u^{\mathsf T}Eu-v^{\mathsf T}Ev,
\]

and (u,v\in[-1,1]^m); hollowness makes the absolute quadratic maximum on
that cube equal to (q(E)).  The one-spin event proved in
`sign_near_weighted_recovery.md` controls
(q(E)=2\max_x|H_E(x)|), and hence already gives
(B(E)=o(m^{3/2})), (1.12), (R.1), and (d_M(T_A,T_W)=o(1)).

For the objective, the two-spin union bound is quantitatively worse.  Its
error (t_m/2) is

\[
\sqrt{2Va_m}+\frac43a_m,
\qquad a_m=2m\log2+\log4.
\]

The direct one-spin union bound may instead use

\[
\sqrt{2Va'_m}+\frac43a'_m,
\qquad a'_m=m\log2+\log4.
\]

Thus the two-spin route costs asymptotically a factor (sqrt2) in the
Gaussian term and (2) in the linear Bernstein term when it is used for the
scalar objective.  These are finite/critical-scale constant losses, not a
surviving (m^{3/2})-scale loss under (V=o(m^2)).  The right assessment is
therefore: correct and useful as an explicit action-norm corollary, stronger
than the existing artifact's stated profile route, but already implicit in
its scalar rounding proof and polarization.
