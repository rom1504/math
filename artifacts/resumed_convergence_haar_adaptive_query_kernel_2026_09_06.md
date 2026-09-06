# Exact conditional kernel for finitely many adaptive Haar-involution queries

Date: 2026-09-06. This independently audits the paired-query recursion
in `resumed_response_mixed_charge_involution_feedback_2026_09_06.md`.
The theorem is finite-depth and nondegenerate. It does not describe the
unrestricted Boolean optimum or a number of queries growing with size.

## 1. Exact posterior, including the longitudinal term

Let `B=O diag(I_r,-I_s) O^T`, where `O` is Haar orthogonal,
`n=r+s`, and `r/n,s/n -> 1/2`. All initial row seeds are independent
of `B`. Query vectors are measurable functions of those seeds and the
previous vectors and responses. After finitely many queries, write

    K=span{v_1,Bv_1,...,v_t,Bv_t}.

This subspace is invariant under `B`, since `B^2=I`. The action on it is
known: it interchanges each recorded query and response. Conditional on
the full history and seeds, the restriction to `C=K^perp` is a Haar
involution with the remaining positive and negative multiplicities.

Indeed, the fiber of matrices compatible with the history consists
exactly of the known restriction on `K` and an arbitrary involution of
these multiplicities on `C`. The original Haar-orbit measure and the
conditioning event are invariant under the orthogonal stabilizer of
`K`. That stabilizer acts transitively on the fiber. The unique invariant
conditional measure is consequently its Haar-orbit measure. Adaptivity
adds no constraints: once previous replies and seeds are fixed, the
subsequent query vectors are already fixed.

For the next query remove its known projection and put
`u=(I-Pi_K)v`, `L=||u||`, `e=u/L`, assuming `L>0`. Let `d=dim C`
and `r',s'` be its remaining multiplicities. Exactly,

\[
z=\langle e,B e\rangle
\ \stackrel{d}=\ 2\operatorname{Beta}(r'/2,s'/2)-1.
\tag{1}
\]

Writing `mu=(r'-s')/d`, its mean is `mu` and its variance is
`2(1-mu^2)/(d+2)`. Conditional on `z`,

\[
B u=L\bigl(z e+\sqrt{1-z^2}\,\xi\bigr),            \tag{2}
\]

where `xi` is uniform on the unit sphere in `C intersect e^perp`,
independent of `z`. Formula (1) follows by expressing `e` in a Haar
eigenbasis and summing the first `r'` squared coordinates. Formula (2)
then follows from the stabilizer of `e` in `C`.

When both multiplicities are positive, `|z|<1` almost surely. In the
basis `(e,xi)`, the restriction to the newly revealed plane is exactly

\[
\begin{pmatrix}z&\sqrt{1-z^2}\\
\sqrt{1-z^2}&-z\end{pmatrix}.                    \tag{3}
\]

Thus each nondegenerate query removes precisely one positive and one
negative eigendirection. The remaining complement again has the exact
conditional law above. In particular, the transverse sphere is exact,
but orthogonality of `u` and `Bu` is only asymptotic, not exact.

## 2. Gaussian replacement at fixed depth

For a fresh independent standard Gaussian vector `g`, let `P` project
onto `C intersect e^perp`. Couple `xi=Pg/||Pg||`. At any fixed query
depth, `P` has fixed codimension, so conditionally

    E||Pg-g||^2=codim(P)=O(1),   ||Pg||/sqrt(n) -> 1.

Also `z -> 0` by (1). If `L^2/n -> v>0`, (2) consequently implies

\[
\frac1n\|Bu-\sqrt v\,g\|^2\longrightarrow0
\tag{4}
\]

in probability. With uniformly integrable normalized query moments it
also holds in mean. No coordinate incoherence of the exposed span is
needed for this normalized-L2 statement: its finite rank alone controls
the projection error.

Induction now gives an iid-row Gaussian representation, deterministic
empirical second moments, and empirical Wasserstein-2 convergence for
any fixed number of globally Lipschitz coordinate queries with suitable
seed moments and strictly positive limiting residual variances. The
Gram matrices of newly appended residual pairs are asymptotically
orthogonal by (4), so inverse-Gram coefficients converge. Bounded
queries automatically supply the needed normalized-moment bounds.

Hard coordinate thresholds can be passed by globally Lipschitz
approximation once their limiting scalar score has no atom at the
threshold. This extension is not asserted without the no-atom check.
Neither growing depth nor a degenerate inverse-Gram limit is covered.

## 3. Scalar recursion and energy orientation

Write the scalar history as the column vector
`H=(V_1,W_1,...,V_t,W_t)^T`, and let `J` interchange each query-response
pair. Put `Gamma=E[H H^T]`. Orthogonality and symmetry of `B` give
`J Gamma J=Gamma` and symmetry of `Gamma J`. For a new coordinate query
`V`, set

\[
b=\mathbb E[HV],\qquad c=\Gamma^{-1}b,
\qquad v=\mathbb E V^2-b^T\Gamma^{-1}b>0.
\]

The next response has the exact limiting representation

\[
W=H^T Jc+\sqrt v\,Z,\qquad
Z\sim N(0,1)\text{ independent of prior row history and seeds}.
\tag{5}
\]

The normalized half-energy is therefore

\[
\frac12\mathbb E[VW]=\frac12 c^T\Gamma Jc.        \tag{6}
\]

Equivalently, orthonormalize the residual pairs. For an orthonormal
history `(q_j,p_j)` whose limiting involution interchanges its entries,

    a_j=E[V q_j], b_j=E[V p_j],
    residual=V-sum_j(a_j q_j+b_j p_j),
    W=sum_j(a_j p_j+b_j q_j)+sigma Z,

and the half-energy is `sum_j a_j b_j`. This is precisely the response
agent's recursion. For the first nontrivial example, start with a
Rademacher `S` and independent `Z_1`; if `V=f(S,Z_1)`, then

    BV = a Z_1+b S+sqrt(EV^2-a^2-b^2) Z_2,
    a=E[SV], b=E[Z_1 V],   half-energy=ab.

The coefficients above are population-analysis coefficients, not
additional empirical operations in the original coordinate algorithm.

## 4. Inertial signs: nondegeneracy and the response audit

For the fixed rule `u_(t+1)=sign(Bu_t+alpha_t u_t)`, every fixed step
has positive residual variance in (5). Inductively, condition on the
row history before the last fresh Gaussian `Z`. The next spin is a
nonconstant sign of an affine function of `Z` with nonzero slope.
Every element of the already exposed linear span is, under the same
conditioning, affine in `Z`: its only new Gaussian basis coordinate is
that last response, and the corresponding residual query predates it.
A nonconstant sign of a nondegenerate Gaussian cannot agree almost
surely with an affine function. Hence the new spin is outside the
finite exposed span in L2. This also proves the threshold's no-atom
condition and licenses soft-sign approximation at every fixed depth.

The strict damped improvement in the response artifact also passes.
For a Boolean `u`, `y=Bu`, `v=sign y`, and `||B||op<=1`,

\[
\frac{H_B((1-\eta)u+\eta v)-H_B(u)}n
\ge\eta g_n-2\eta^2,\qquad
g_n=\frac1n\sum_i(|y_i|-u_i y_i),\quad H_B(x)=x^TBx/2.
\tag{7}
\]

If the scalar response is `y=k+sigma Z` with `sigma>0`, then
`g=E[E_Z|k+sigma Z|-u k]>0`. The strict pointwise inequality
`E_Z|k+sigma Z|>|k|>=u k` suffices. Also `g<=2`, since
`E y^2=1`. Taking the fixed `eta=g/4` gives a limiting gain at least
`g^2/8`. This does not give a uniform-in-depth gain. The displayed
Jensen refinement in that artifact is valid: for
`Gamma_sigma(z)=E|z+sigma Z|-z`, `z>=0`, the second derivative of
`Gamma_sigma(sqrt x)` has numerator
`z Gamma_sigma''(z)-Gamma_sigma'(z)>0`.

The independent strengthening in
`resumed_bound_audit_mixed_charge_feedback_2026_09_06.md`, Sections 4--5,
has also been reconstructed and passes: a conditional best-affine
Gaussian projection gives an explicit positive innovation lower bound
uniform over bounded depth and bounded inertia. Conversely, random-gated
damped feedback produces fixed finite-depth states with arbitrarily
small positive stability gap. In view of the separate Haar performance
ceiling, a positive gain uniform solely in the deficit from `1/2` is
therefore false for this larger feedback class. This does not assert
that the particular synchronous inertial schedule converges.

## 5. Scope of deterministic-matrix transfer and primary corroboration

Transfer to deterministic matrices uses the separately audited WZF
fixed-rule universality theorem, with signed-permutation-equivariant
globally Lipschitz approximants and the stated matrix hypotheses.
The safe classes here are exact flat symmetric involutions (including
normalized symmetric Hadamards), or operator-close perturbations of an
explicit comparator satisfying those hypotheses. Merely assuming
`||B^2-I||op=o(1)` does not by itself verify the required flat power
conditions. This note does not strengthen that scope.

The posterior for an involution in Section 1 is proved here, not
imported from a theorem about an unconstrained Haar orthogonal matrix.
For closely related primary conditioning and finite-rank sphere facts,
see Zhou Fan, *Approximate Message Passing algorithms for rotationally
invariant matrices*, arXiv:2008.11892v5, Appendix F, Propositions F.1 and
F.2, pp.84--85:

<https://arxiv.org/pdf/2008.11892>

Proposition F.1 states the Haar orthogonal conditional representation;
Proposition F.2 and its proof state the Gaussian empirical limit for
a uniform sphere embedded in a fixed-codimension complement. Neither
is being misquoted as an involution-specific posterior theorem.
