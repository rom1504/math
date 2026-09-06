# Wave 28 Route 8: completion entropy and collisions

## Status

- **Proved:** a low-information, low-row completion channel contains a single
  low-row cut with exponentially large selector degree. Randomization cannot
  evade the fixed-cut collision problem.
- **Proved:** a microcanonical child-state entropy bound controls both mutual
  information and row square. The dependence on
  \(\lVert A[:,S]y\rVert_2^2\) is retained and bounded by an
  entropy/Hanson--Wright argument.
- **Proved:** conditioning a row-Gibbs prior on the admissible incidence
  relation gives an exact pressure and derivative criterion.
- **Numerical/exhaustive:** all incidence, entropy, row, and pressure
  identities pass on \(A_6,A_8,A_9\).
- **Open:** no asymptotic exact-minimizer theorem supplies the required
  microcanonical entropy or incidence pressure.

## 1. The correct microcanonical relation

Let \(\mathcal S=\binom{[n]}m\), \(N=|\mathcal S|\), and let
\(\mathcal D\) be the \(2^n\) oriented projective full cuts. For
\(d=(\sigma,x)\), define the child deficit

\[
\delta_S(d)=Q(A[S])-\sigma x_S^{\mathsf T}A[S]x_S\ge0.
\tag{R28.1}
\]

The definition (10.792) gives the exact identity

\[
\widehat\ell(S,d)
=\delta_S(d)+p_2\sigma x^{\mathsf T}Ax-p^{3/2}q_n.
\]

Since every oriented parent payoff is at most \(q_n\),

\[
\boxed{\widehat\ell(S,d)\le\delta_S(d)-B_{n,m}.}
\tag{R28.2}
\]

Thus for any threshold \(t\), the admissible incidence relation is

\[
S\sim_t d
\quad\Longleftrightarrow\quad
\delta_S(d)\le B_{n,m}+t,
\tag{R28.3}
\]

and every incident pair has \(\widehat\ell(S,d)\le t\). Exact child
grounds are the stronger subrelation \(\delta_S=0\).

## 2. Low information forces a single collision

For either the exact or microcanonical relation, put

\[
\mathcal C_t(d)=\{S:S\sim_t d\},\qquad r_t(d)=|\mathcal C_t(d)|.
\tag{R28.4}
\]

Consider any channel with \(S\sim U_m\) supported on \(S\sim_tD\).
Conditional on \(D=d\), the selector posterior has support at most
\(r_t(d)\). Hence

\[
\boxed{
I(S;D)\ge
\mathbb E\log\frac{N}{r_t(D)}.
}
\tag{R28.5}
\]

Suppose \(I(S;D)\le K\) and \(\mathbb ER_2(D)\le C\), with
\(K,C>0\). For every \(\alpha>1\), some output \(d\) satisfies

\[
\boxed{
\frac{r_t(d)}N\ge e^{-\alpha K},
\qquad
R_2(d)\le\frac{\alpha}{\alpha-1}C.
}
\tag{R28.6}
\]

Proof: set \(L(d)=\log(N/r_t(d))\) and
\(a=(\alpha-1)K/C\). Equations (R28.5) and the cost assumption give
\(\mathbb E[L+aR_2]\le K+aC=\alpha K\). Some output is no larger
than this mean, and its two nonnegative summands give (R28.6). For
\(\alpha=2\), the constants are exactly \(e^{-2K}\) and \(2C\).
When \(K=0\), (R28.5) gives \(r_t(D)=N\) almost surely.

At the target scales

\[
K=O(n^{3/4-c}),\qquad C=O(n^{9/4-c}),
\tag{R28.7}
\]

(R28.6) is already the single-cut lemma (10.795). Therefore any successful
randomized completion rule automatically contains a successful fixed cut.
For exact child completions this is stronger than (10.795), since (R28.2)
gives threshold \(-B_{n,m}\). Per-selector completion bounds cannot bypass
the missing collision theorem.

There is also an exact falsification test. If, for some \(\alpha>1\), every
cut of row square at most \(\alpha C/(\alpha-1)\) has coverage below
\(e^{-\alpha K}\), no channel supported on the relation can have both
information at most \(K\) and average row square at most \(C\).

## 3. Microcanonical entropy gives a complete channel criterion

Let \(\mathcal F_S(t)\) be the oriented projective child states satisfying
\(\delta_S\le B_{n,m}+t\), and write \(f_S(t)=|\mathcal F_S(t)|\).
There are \(2^m\) oriented projective child states in total. Define

\[
\boxed{
\mathcal K_t
=\mathbb E_{S\sim U_m}\log\frac{2^m}{f_S(t)}.
}
\tag{R28.8}
\]

Choose an element of \(\mathcal F_S(t)\) uniformly and fill every outside
relative spin uniformly. Equivalently, condition the uniform law on full
oriented cuts on \(S\sim_tD\). The conditional support has size
\(f_S(t)2^{n-m}\), so information radius relative to the uniform full-cut
law gives

\[
\boxed{I(S;D)\le\mathcal K_t.}
\tag{R28.9}
\]

The output entropy is generally linear because the common outside
randomness already contributes \((n-m)\log2\). Mutual information is the
right quantity: this irrelevant conditional entropy cancels.

The row term requires a separate audit. Let \(Y_S\) denote the selected
relative child spin (forgetting orientation, while retaining its possible
multiplicity). Uniform outside completion gives the exact identity

\[
\boxed{
\mathbb E[R_2(D)\mid S]
=(n-m)(n-1)+
\mathbb E\lVert A[:,S]Y_S\rVert_2^2.
}
\tag{R28.10}
\]

The second term cannot simply be replaced by its uniform-cube mean. Let
\(P_S\) be the law of a signed representative of \(Y_S\), let \(U\) be
uniform on \(\{\pm1\}^m\), and put \(B_S=A[:,S]\),
\(M_S=B_S^{\mathsf T}B_S\). Data processing from the uniform law on
oriented projective states gives

\[
D(P_S\Vert U)\le K_S,\qquad
K_S=\log\frac{2^m}{f_S(t)}.
\tag{R28.11}
\]

The Rademacher Hanson--Wright mgf and entropy duality imply, for an absolute
constant \(C_{\rm HW}\),

\[
\mathbb E_{P_S}Y^{\mathsf T}M_SY
\le \operatorname{tr}M_S+
C_{\rm HW}\left(
\lVert(M_S)_{\rm off}\rVert_F\sqrt{K_S}
+\lVert(M_S)_{\rm off}\rVert_{\rm op}K_S
\right).
\tag{R28.12}
\]

For completeness, the input is the Bernstein-form mgf

\[
\log\mathbb E_U
\exp\{u(Y^{\mathsf T}M_SY-\operatorname{tr}M_S)\}
\le
\frac{C u^2\lVert(M_S)_{\rm off}\rVert_F^2}
{1-Cu\lVert(M_S)_{\rm off}\rVert_{\rm op}}.
\]

Substituting this in the entropy variational inequality and optimizing \(u\)
proves (R28.12). The matrix norms are

\[
\begin{aligned}
\operatorname{tr}M_S&=m(n-1),\\
\lVert(M_S)_{\rm off}\rVert_F
&\le\lVert A\rVert_{\rm op}\sqrt{m(n-1)},\\
\lVert(M_S)_{\rm off}\rVert_{\rm op}
&\le\lVert A\rVert_{\rm op}^2+n-1.
\end{aligned}
\tag{R28.13}
\]

For an exact minimizer,
\(\lVert A\rVert_{\rm op}^2\le2q_n=O(n^{3/2})\).
Averaging (R28.10)--(R28.13) and using Jensen yields the audited row bound

\[
\boxed{
\mathbb ER_2(D)
\le n(n-1)
+O\left(n^{7/4}\sqrt{\mathcal K_t}
+n^{3/2}\mathcal K_t\right).
}
\tag{R28.14}
\]

Consequently, for fixed \(0<c<1/4\),

\[
\boxed{
t=O(n^{3/2-c}),\qquad
\mathcal K_t=O(n^{3/4-c})
}
\tag{R28.15}
\]

is a complete sufficient lemma: (R28.9) gives the information budget,
(R28.14) gives \(\mathbb ER_2=O(n^{9/4-c})\), and (R28.2)
gives the loss budget. The first error in (R28.14) has exponent
\(17/8-c/2<9/4-c\), while the second has exactly the target exponent.
This proves the power-saving restriction edge through (10.794).

For exact grounds, if \(g_S\) is the number of oriented projective exact
grounds, then

\[
\mathcal K_{\rm exact}
=\mathbb E\log\frac{2^m}{g_S}
=\log2+\mathbb E\log\frac{2^{m-1}}{g_S}.
\tag{R28.16}
\]

Thus high exact-ground entropy is a special, much stronger sufficient
condition. The microcanonical formulation (R28.8) is the correct target.

## 4. The macroscopic-threshold \(n^{3/4}\) wall

The microcanonical entropy criterion has an exact conditional obstruction.
Write

\[
h_S(t)=Q(A[S])-B_{n,m}-t.
\tag{R28.17}
\]

Under a uniformly random oriented projective child state, membership in
\(\mathcal F_S(t)\) is the upper-tail event
\(\sigma y^{\mathsf T}A[S]y\ge h_S(t)\). Rademacher Hanson--Wright gives

\[
\frac{f_S(t)}{2^m}
\le
\exp\left[
-c_0\min\left\{
\frac{(h_S(t)_+)^2}{m(m-1)},
\frac{h_S(t)_+}{\lVert A[S]\rVert_{\rm op}}
\right\}\right]
\tag{R28.18}
\]

up to an inessential absolute prefactor, which contributes only \(O(1)\)
to its negative logarithm. Since
\(\lVert A[S]\rVert_{\rm op}\le\lVert A\rVert_{\rm op}=O(n^{3/4})\),
if \(h_S(t)\ge\gamma n^{3/2}\) on a selector set of fixed positive
measure, then

\[
\boxed{\mathcal K_t=\Omega(n^{3/4}).}
\tag{R28.19}
\]

This misses the desired \(O(n^{3/4-c})\) by a power. Thus the elementary
uniform microcanonical channel cannot close in any regime with a
macroscopic residual energy threshold on a positive fraction of selectors.
Escaping this wall requires at least one of:

1. \(h_S(t)=o(n^{3/2})\) for almost all selectors, which is already a
   strong restriction statement;
2. exceptional spectral/large-deviation structure giving far more upper-tail
   mass than generic Hanson--Wright permits; or
3. a nonuniform prior whose incidences concentrate on shared low-row cuts.

Statement (R28.19) is conditional on the displayed macroscopic lower bound;
no such lower bound is asserted uniformly without an additional theorem.

## 5. Row-Gibbs incidence pressure

The most general canonical version of option 3 is exact. For the chosen
relation \(S\sim_td\), define

\[
Z(\lambda)=\sum_de^{-\lambda R_2(d)},\qquad
Z_S(\lambda)=\sum_{d:S\sim_td}e^{-\lambda R_2(d)},
\]

\[
\pi_\lambda(d)=\frac{e^{-\lambda R_2(d)}}{Z(\lambda)},\qquad
a_S(\lambda)=\frac{Z_S(\lambda)}{Z(\lambda)},\qquad
\mathcal P_t(\lambda)=\mathbb E_S[-\log a_S(\lambda)].
\tag{R28.20}
\]

Condition \(\pi_\lambda\) on \(S\sim_tD\). Information radius and exact
differentiation give

\[
\boxed{
\begin{aligned}
I(S;D)&\le\mathcal P_t(\lambda),\\
\mathbb ER_2(D)
&=\mathbb E_{\pi_\lambda}R_2(D)+\mathcal P_t'(\lambda)\\
&=\mathbb E_S\mathbb E_{\pi_\lambda}
[R_2(D)\mid S\sim_tD].
\end{aligned}
}
\tag{R28.21}
\]

Indeed

\[
\frac d{d\lambda}[-\log a_S]
=\mathbb E_{\pi_\lambda}[R_2\mid S\sim_tD]
-\mathbb E_{\pi_\lambda}R_2.
\]

There is no general sign for \(\mathcal P_t'\). The exact surviving pressure
lemma is the existence of a target-specific \(\lambda\ge0\) for which

\[
\boxed{
\mathcal P_t(\lambda)=O(n^{3/4-c}),\qquad
\mathbb E_{\pi_\lambda}R_2+\mathcal P_t'(\lambda)
=O(n^{9/4-c}).
}
\tag{R28.22}
\]

By (R28.6), even this criterion necessarily produces a single low-row cut
with the required microcanonical coverage. It is therefore a reformulation
of the collision problem, not a way around it.

## 6. A scoped combinatorial wall

Large completion fibers alone do not force overlap. Assign each \(m\)-set
an independent uniform pattern \(y_S\), and let its fiber be
\(F_S=\{x:x_S=y_S\}\). Every fiber has \(2^{n-m}\) elements. For fixed
density \(p=m/n\), however, there is a constant \(\kappa(p)>0\) and a
deterministic labeling for which

\[
\max_x|\{S:x\in F_S\}|\le Ne^{-\kappa n}.
\tag{R28.23}
\]

For fixed \(x\), its degree is \(\operatorname{Bin}(N,2^{-m})\).
Taking \(0<\kappa<\min\{p\log2,H(p)\}\), Chernoff at
\(Ne^{-\kappa n}\), followed by a union bound over \(2^n\) points, proves
(R28.23). One can assign zero abstract cost to every edge. Thus
per-selector low cost and full completion freedom do not imply collisions.

This is deliberately abstract: arbitrary labels need not be realizable as
child optima of one signing matrix. A proof must use shared exact-minimizer
structure absent from (10.826)--(10.830).

## 7. Finite audit and frontier

The checker completion_collision_r28.py exhaustively verifies (R28.5),
(R28.6), (R28.10), and (R28.21) on all slices
\(m\ge\lceil n/2\rceil\) of \(A_6,A_8,A_9\). For these small matrices
\(B_{n,m}\) lies below the first positive deficit lattice step, so the
zero-threshold microcanonical relation coincides with exact-ground incidence.
Actual mutual information ranges from \(0.4377\) to \(2.3271\) nats, and
the largest exact selector coverage ranges from \(1/4\) to \(5/6\).
These are exhaustive finite observations, not asymptotic evidence.

The best new exact target is (R28.22), with (R28.15) as a simple
microcanonical-entropy sufficient condition and (R28.19) as its conditional
power-scale wall. Because (R28.6) extracts fixed-cut coverage from every
successful channel, this route is no easier than approximate alignment
(10.837). It should remain a diagnostic or backup unless a new
exact-minimizer exchange theorem controls the incidence pressure or directly
forces a high-degree low-row cut.
