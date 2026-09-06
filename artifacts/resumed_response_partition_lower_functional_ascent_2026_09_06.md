# Feasible ascent of a conditional-partition response lower functional

Date: 2026-09-06. Status: exact derivative reconstructed locally; numerical
test not yet performed. This is a lower-certificate optimization mechanism,
not a proof of globally convergent response iteration.

Let (F,q) be an odd/even conditional feasible pair, |F|<=q<=1, and let
K=U*F, t=Eq-||K||^2>0. Fix a finite even partition C_i of the old Gaussian
space (in the application, symmetric bins of V). Put

    M_i=E[1_(C_i)(1-q)], I_i=E[1_(C_i)(1-q)K], k_i=I_i/M_i,
    L(F,q)=sum_i M_i Psi(k_i,t).

Terms with M_i=0 are zero. Weighted conditional Jensen and actual even-gate
purification make L an attainable ORIGINAL signing lower functional.
Its status does not depend on the conditional law of a matrix response
being identified completely.

Suppose M_i>0 on the included cells. Define

    s_i=2Phi(k_i/sqrt(t))-1,
    D_i=Psi(k_i,t)-k_i*s_i=2sqrt(t)phi(k_i/sqrt(t)),
    B=sum_i M_i*phi(k_i/sqrt(t))/sqrt(t).

Differentiating the perspective M Psi(I/M,t), including the dependence of
K and t on F, gives exactly

    dL = <U[(1-q)s_bin-2B K],delta F>
          +E[(B-D_bin-s_bin K)delta q].              (1)

Indeed dI_i=E[1_i(1-q)U*delta F]-E[1_i K delta q],
dM_i=-E[1_i delta q], and dt=E delta q-2<K,U*delta F>.
Thus (1) follows without replacing K by a conditional expectation inside
its norm. Partitioning lowers the functional, but not by a hidden change
in first-chaos geometry.

For variations restricted to a finite standard independent Gaussian
frame X=(X_1,...,X_d), conditional expectation of the first-chaos score
in (1) is the linear Gaussian field

    R(X)=sum_j r_j X_j,
    r_j=E[(1-q)s_bin h_j]-2B<K,h_j>, h_j=U*X_j.

The q coefficient is conditioned on X, giving the pointwise linear
oracle

    F_new=sign R(X)*q_new,
    q_new=1{|R(X)|>D_bin+s_bin E[K|X]-B}.            (2)

Ties may be chosen arbitrarily. Formula (2) maximizes the directional
derivative over the entire conditional feasible domain measurable in X.
It need not be the full unrestricted gradient update. A line search
along the conditional convex segment preserves feasibility and must use
t(theta)=Eq(theta)-||K(theta)||^2, with even-gate realization before the
matrix limit. Positive derivative ensures improvement for a sufficiently
small step, but does not license a full step.

In the rich-core application X=(V,Z), h_V=g and h_Z=(f-A g)/nu.
The only non-scalar object needed by (2) is the SAME bivariate polynomial
E[g|V,Z] already certified for the first response birth. The exact
principal-angle tail estimate controls truncation. Rectangle versions of
(2), including arbitrary signs on finitely many Z intervals, have all
their sufficient statistics given by Gaussian Hermite endpoint sums.
Thus further certified ascent within this fixed response frame does not
automatically require a third correlated Gaussian or a refreshed-spin
assumption.

This theorem does not ensure a positive restricted gradient everywhere,
give a frame-independent gain, or identify the global maximum of L.
It offers a precise, reproducible finite-data test of whether the first
rich-core response has exhausted its own actual two-Gaussian frame.
