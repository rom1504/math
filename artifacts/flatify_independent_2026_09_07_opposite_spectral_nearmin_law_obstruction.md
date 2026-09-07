# Opposite-spectral Gaussian-sign law: an actual near-minimizer obstruction

2026-09-07. Director's proposed falsifier, derived using the new audited
quenched universality theorem. This is a failure of a specified RANDOM
BRIDGE LAW on actual additive near-minimizers. It is not a counterexample
for every bridge selection, nor a claim about every exact minimizing child.

## 1. Weak spectral normalization makes the Gaussian bridge asymptotically iid

Let A_n,D_n be hollow symmetric sign matrices with Q(A_n),Q(D_n)=O(n^(3/2)),
and put p_n=||A_n||op||D_n||op. Assume p_n/n->infinity. The actual bridge
law is C=sign G with

    Cov(vec G)=I-rho A_n tensor D_n/p_n.

Here rho<1 can be any fixed value, or the explicitly paid sequence
rho_n=1-n^(-1/7) from the quenched universality theorem. Its matched
Gaussian bridge Y has covariance

    S=I-kappa A_n tensor D_n,
    kappa=(2/pi)arcsin(rho/p_n)=O(1/p_n).

For any two Boolean rank-one features v=x tensor y and w=x' tensor y',
the covariance difference from an iid standard Gaussian bridge Z obeys

 |v^T(S-I)w|
    <=kappa beta(A_n)beta(D_n)=O(n^3/p_n).             (1)

The same estimate holds after adjoining the polarity used for an
absolute maximum. Gaussian maximum comparison yields

    |E beta(Y)-E beta(Z)|<=C n^2/sqrt(p_n)
                                  =o(n^(3/2)).        (2)

To check (2), a covariance-entry bound delta bounds increment-variance
differences by4delta. Add independent Gaussian noises of variance2delta
at every index to the smaller-increment process and use Sudakov--Fernique.
Their expected maximum is at most C sqrt(delta log M), with
M<=2^(2n+1). Reverse the roles for the other direction. This proves the
displayed estimate without a fictitious entrywise coupling.

By the already banked iid Gaussian bipartite bridge/SK comparison,

    E beta(Z)>=[2 P_SK-o(1)]n^(3/2),
    P_SK>3/4.                                        (3)

The control proof and its certified stronger value are in
`flatify_independent_2026_09_07_heat_martingale_sk_bound.md` and
`flatify_independent_2026_09_07_balanced_iid_bridge_exact_limit.md`.
Only the conservative P_SK>3/4 is needed here.

Apply the new quenched Gaussian-sign universality theorem with zero
child offsets. Equations (2)--(3) give for the ACTUAL sign bridge

    E beta(C)>=[2P_SK-o(1)]n^(3/2).                   (4)

## 2. Concentration of the actual bridge cap

More generally let C=sign G with epsilon I<=Cov(G)<=KI. Write
G=W+sqrt(epsilon)Z as in the universality proof. Conditional on W,
the signs are independent. Changing one sign entry changes beta(C)
by at most2, so conditional Efron--Stein gives

    E Var(beta(C)|W)<=2n^2.

Let m(W)=E[beta(C)|W]. The i-th conditional success probability is
Phi_normal(W_i/sqrt(epsilon)); differentiating this product expectation
and using the same coordinate-change bound gives

    |partial_i m(W)|<=2/sqrt(2pi epsilon).

Gaussian Poincare, with ||Cov W||op<=K, therefore gives

    Var m(W)<=2K n^2/(pi epsilon).

Thus

    Var beta(C)<=2n^2+2K n^2/(pi epsilon).             (5)

For a fixed gap, or epsilon_n=n^(-1/7), this is o(n^3). Chebyshev
and (4) show

    beta(C)>=[2P_SK-o(1)]n^(3/2)

with probability tending to one. The o(1) threshold may be chosen
slowly enough to dominate the variance bound.

For ANY internal children, flipping one entire block reverses the
bridge and preserves both internal energies. Hence the full absolute
parent cap is at least beta(C), deterministically. With N=2n,

    Q(parent)/N^(3/2)>=P_SK/sqrt(2)-o(1)
                      >3/(4sqrt(2))-o(1)>0.53-o(1)   (6)

with high probability. This is strictly above the banked original
upper constant below0.494.

## 3. Actual near-minimizing children satisfying the hypothesis

Start from an exact minimizing full signing A_n. Choose

    s_n=floor(n^(3/4)/log n)

vertices and set every internal edge of this set to +1. Call the
result B_n. At most binom(s_n,2) coefficients change, each by at most2,
so

    M_n<=Q(B_n)<=M_n+s_n(s_n-1)=M_n+o(n^(3/2)).       (7)

The unit vector constant on the planted set and zero elsewhere has
Rayleigh quotient s_n-1; therefore ||B_n||op>=s_n-1. Take both actual
children to be B_n. Then

    p_n=||B_n||op^2 >=(s_n-1)^2,
    p_n/n >=(1-o(1))n^(1/2)/(log n)^2 ->infinity.

The law obstruction (6) therefore applies to a sequence of genuine
full-sign, additive near-minimizers. Neither (7) nor this proof assumes
that M_n/n^(3/2) converges.

## Exact scope

The specified p-normalized opposite-spectral random law does not provide
a near-minimizer-uniform favorable composition: for the children above,
almost every sampled bridge has the wrong leading constant. Its rare
outputs have not been excluded, and an independently selected favorable
bridge for these fixed children is not ruled out. The planting step does
not prove the same obstruction for every exact minimizer. The result
also does not reject other covariance kernels, global internal-edge
changes, or the original convergence conjecture.
