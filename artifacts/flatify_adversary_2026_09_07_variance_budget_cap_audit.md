# Independent audit of the direct variance-budget cap operation

Date: 2026-09-07. Verdict: PASS for
`flatify_director_variance_budget_cap_comparison_2026_09_07.md`.
The direct N^(5/4) count-cleanup bound is proved independently of pressure
or temperature smoothing. It is not a flatification endpoint theorem.

## Rounding and contraction

Let v,w be profiles on d<=N vertices and w_ij>=s_i^2 s_j^2 v_ij, with
0<=s_i<=1. For an optimal signing A of v choose independent C_ij in
{+1,-1} with mean s_i s_j sqrt(v_ij/w_ij) A_ij whenever w_ij>0.
When w_ij=0 domination forces the desired mean coefficient to be zero,
and its arbitrary sign has no effect. The centered coefficient has
variance w_ij-s_i^2 s_j^2 v_ij and absolute value at most b=2K.

For each spin vector the variance sum is T. Bernstein's bound with
L=(N+2)log2 and t=sqrt(2TL)+2bL/3 gives failure probability at most
2exp(-L). This value of t is safe: t^2>=2L(T+bt/3). A union bound over
at most 2^N vectors gives total failure probability at most 1/2.

The deterministic mean polynomial at x equals the original polynomial at
(s_i x_i). Separate affinity bounds its absolute value throughout the
cube by the original vertex cap. Therefore an actual full choice of signs
attains the stated q(v)+sqrt(2TL)+4KL/3 bound. No covariance surrogate is
substituted for the output signing.

## Simultaneous deletion and count comparison

Restriction cannot increase the minimum achievable cap: average over the
deleted spins for each fixed signing, then minimize. Conversely retain an
optimal restricted signing and round all newly incident edges independently
with mean zero. Their total variance is at most r(N-1), since summing the
r removed row budgets counts internal removed edges twice. This proves
the simultaneous insertion bound, with only one Bernstein linear term.

For two posterior counts differing by D, delete D known vertices so the
remaining known labels coincide. The common known-known weights agree;
nonzero known-unknown weights have relative difference O_delta(D/u), as
do unknown-unknown weights (the latter have a uniform positive lower
bound for u>=4). Scaling only the unknown vertices therefore gives
bidirectional coordinate domination for D/u sufficiently small.

The exact edge-mass bookkeeping is useful here. If E_removed denotes
the internal mass among the D deleted vertices, then

    remainder mass = N(N-1)/2 - D(N-1) + E_removed.

Thus the two remainder masses differ by O_delta(D^2), not by O(ND).
The contraction loses at most O_delta(ND), using only the unknown row
budgets; hence the actual introduced variance T is O_delta(ND). Its
nonnegativity follows from coordinate domination, regardless of the sign
of the remainder mass difference. Since D<=N, the O(D^2) discrepancy
fits this budget. Rounding and simultaneous restoration cost
O_delta(N(sqrt(D)+1)) in either direction.

If D/u is not sufficiently small, split the count interval into a bounded
number of shorter intervals, depending only on delta. This is valid for
u sufficiently large depending on delta, which is all the bulk argument
needs. The finitely small-u cases are covered directly by the terminal
deletion bound; one need not force the subdivision through indivisible
integer count steps. For D=0 the comparison is exact, though the stated
upper bound harmlessly includes an additive N term.

## Uniform expected error and exact scope

With u>=sqrt(N), the central event has E sqrt(D)=O(u^(1/4)) and complement
probability O_delta(1/u). The central contribution is O_delta(N^(5/4));
the complement is O_delta(N^(3/2)/u), using unbiased full-sign rounding
for a universal row-regular cap bound. With u<sqrt(N), each profile and
any consistent terminal completion agree on the known submatrix.
Simultaneous insertion gives distance O_delta(Nsqrt(u)+N) from the same
terminal optimum. This is again O_delta(N^(5/4)).

All comparisons concern actual weighted full signings. Nevertheless,
direct terminal-to-uniform flatification by these contractions introduces
Theta(N^2) variance and therefore still pays Theta(N^(3/2)) in the
rounding estimate. The operation rigorously handles small variance-budget
changes; it does not supply the missing favorable endpoint drift.
