# Low-rank orthogonality does not remove macroscopic energy width

Status: proved; independently reconstructed by the constructive and independent
researchers. This is a theorem about actual sign matrices, not a surrogate
ensemble. It restricts a scalar certificate, not the actual cap of a bridge.

## Theorem

Let A be a symmetric hollow sign matrix of order n with Q(A)<=C n^(3/2),
where C is fixed. Let X be ANY prescribed subspace of rank r=o(sqrt(n)).
There is a sequence tau_n=o(1) such that the set

    S_X={x in {+1,-1}^n: ||P_X x||^2 <= tau_n n}

has energy half-width at least

    [4/(3*pi*sqrt(3))-o(1)] n^(3/2).                 (1)

The result is uniform over A and X with the indicated cap/rank bounds.
For 1<=r<=n^(1/2-2*delta), one can take tau_n=n^(-delta), with the
o(1) in (1) replaced by O_C(n^(-delta)+n^(-1)).
Half-width means (max H_A-min H_A)/2 on the declared set; its midpoint
need not be zero. No spectral flatness, minimizer hypothesis, or convergence
assumption is used.

## Elementary spectral consequence of bounded Boolean cap

Write beta(A)=max_{u,v Boolean}|u^T A v|. Polarization on the cube gives
beta(A)<=4Q(A). For completeness, u=(x+y)/2 and v=(x-y)/2 give
x^T A y=u^T A u-v^T A v, and each quadratic cube value has magnitude
at most 2Q(A), since A is hollow and its quadratic is multiaffine.

Every row a_i of A belongs to the cube. Thus

    sum_j |(A^2)_ij| = max_{s Boolean} a_i^T A s <= beta(A).

Since A^2 is symmetric, its operator norm is at most its maximum absolute
row sum. Consequently ||A||op^2<=beta(A)<=4C n^(3/2). This proof avoids
any complex interpolation normalization.

## Gaussian witness with bounded covariance

Fix a vertex split I,J, with k=|I|=floor(n/3), l=n-k, and B=A[J,I].
For a standard Gaussian g in R^k define

    z_I=sign(g),    z_J=sign(Bg),
    z^+=(z_I,z_J),  z^-=(z_I,-z_J).

All Gaussian coordinates have positive variance, so ties have probability
zero. The cross energy L=z_J^T B z_I satisfies exactly

    E L = k*l*(2/pi)*arcsin(1/sqrt(k))
        = [4/(3*pi*sqrt(3))+O(1/n)] n^(3/2).        (2)

Indeed each entry B_ji is a sign and the correlation between g_i and
(Bg)_j/sqrt(k) is B_ji/sqrt(k). The bivariate Gaussian sign identity
then gives each summand (2/pi)*arcsin(1/sqrt(k)).

The standardized pre-sign Gaussian vector has correlation matrix

    R = [Id; B/sqrt(k)] [Id; B/sqrt(k)]^T,
    ||R||op = 1+||B||op^2/k <= 1+4C n^(3/2)/k.

Its sign covariance is K=(2/pi)*arcsin[R], entrywise. Its operator norm
is at most ||R||op. To verify this without a hidden entrywise norm claim,
expand arcsin into its positive odd-power series; after multiplying by
2/pi the coefficients sum to one. Schur multiplication by a PSD matrix
with diagonal one is a positive unital map and hence contracts the
operator norm of symmetric matrices: apply it to -t Id<=D<=t Id.
Induction gives ||R^{circ j}||op<=||R||op, proving the claim by the
convergent series. The covariance of z^- is a diagonal-sign conjugate
of K and has the same bound.

It follows for either polarity that

    E ||P_X z^+/-||^2 <= r*(1+4C n^(3/2)/k).       (3)

The probability that either projection exceeds tau*n is therefore at
most 2r*(1+4C n^(3/2)/k)/(tau*n).

## Keeping the width witness after conditioning

The internal energies are unchanged by flipping J, so

    H_A(z^+)-H_A(z^-)=2L,    |L|<=Q(A).

The cross term L need not be nonnegative pointwise; the proof only uses
its positive expectation and the uniform absolute bound. On the event
G that BOTH projections obey the threshold,

    E[L 1_G] >= E L-Q(A)*Pr(G^c).

Hence some good Gaussian sample has L at least the right side, whenever
it is positive (dividing by Pr(G)<=1 only improves this bound).
Its two spins prove (1). For r=o(sqrt(n)), choose tau tending to zero
more slowly than r/sqrt(n). For the stated power bound take tau=n^-delta.
All constants are uniform in the deterministic split and prescribed X.

## What this rules out, and what it does not

Combine two children of order n and two arbitrary subspaces of total rank
o(sqrt(n)). The newly constructed sign bridge has the upper certificate

    |x^T C y| <= n^(3/2)*sqrt((1-p_x)*(1-p_y))+o(n^(3/2)),
    p_x=||P_X x||^2/n, p_y=||P_Y y||^2/n.

Even if the child's signed energy envelopes are known EXACTLY, maximizing
this scalar upper certificate costs at least

    [1+8/(3*pi*sqrt(3))-o(1)] n^(3/2) > 1.4900 n^(3/2). (4)

Proof: on the two near-orthogonal sets, the sum of the energy intervals
has width at least four times the constant in (1); its largest absolute
value is at least half that width, irrespective of interval midpoints or
child polarity. The residual factor is at least 1-tau. The target for
equal optimal children is 2*sqrt(2)*m_n n^(3/2), whose limsup is below
2*sqrt(2)*.493608094 < 1.397. Thus this scalar certificate cannot close
favorable flatification anywhere in the entire subleading-rank range.

Equation (4) is NOT a lower bound on the actual parent cap. It maximizes
an upper envelope for the bridge; the bridge need not attain its residual
norm on the high-energy pairs. Joint response information or global changes
of old edges may still make the actual construction favorable. Neither the
new full-sign subspace construction nor this obstruction proves convergence
or nonconvergence of the original optimum.
