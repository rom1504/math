# Recursive full Hadamard frames: exact joint law and the norm-loss boundary

Date: 2026-09-07. This follows the director's concrete proposal to replace
fixed Walsh frames by independently randomized recursive Hadamard frames.
The ensemble and recursive law below are exact. A strict-subhalf bound on
its full cube has NOT been proved here.

## 1. Actual full sign recursion

Given independent full Hadamards F1,F2 of order k, form

    F_(2k) = [F1  F2; F1 -F2].

This is a full Hadamard: the two column groups are orthogonal and each
has Gram2k I. Fresh signed row/column permutations may be inserted at
every node. Thus it supplies actual rank-two cross sign tiles at every
admissible order, with independent frame draws in distinct fibres.

Put U_i=F_i/sqrt(k). For an input row x=(xL,xR), the exact transform is

    F_(2k)^T x/sqrt(2k)
       = (U1^T z+, U2^T z-),
    z+=(xL+xR)/sqrt(2), z-=(xL-xR)/sqrt(2).               (1)

For Boolean x, z+ and z- are complementary ternary words with nonzero
magnitude sqrt(2). Under a fresh parent input gauge their support sizes
are J and k-J with J binomial(k,1/2). Independent input gauges and frame
draws in the two children make their OUTPUT profile laws independent
conditional on J. Complementary support energies are retained in (1);
the two channels have not been separately charged as full Boolean cubes.

For a general finite input word, uniformly signed-permute its coordinates
before forming pairs. Conditional on the resulting signed pair table pi,
the two child absolute empirical laws are those of |A+B|/sqrt(2) and
|A-B|/sqrt(2). Their output laws are independent conditional on this table,
and their parent output histogram is the half-half mixture. Therefore the
FULL histogram transition law obeys an exact finite recursion:

    P_(r,2k)(nu -> rho)
      = sum_pi Pr_nu(pi)
        sum_(rho+,rho-): (rho++rho-)/2=rho
          P_(r-1,k)(nu+ -> rho+)
          P_(r-1,k)(nu- -> rho-).                         (2)

Here Pr_nu is the literal finite random signed-pairing probability, not a
large-deviation replacement. Terminal frames may be fixed Hadamards with
fresh input gauges. This is a computable actual ensemble law covering ALL
input spectra, rather than an assumed Gaussian limit.

## 2. Annealed counts legitimately enter the typed certificate

For a random frame F let c_F(rho) be the number of Boolean inputs with
output magnitude type rho. Set cbar_k(rho)=E_F c_F(rho). Equivalently,
cbar_k(rho)=2^k times the output-type probability for a uniform Boolean
input and independent F; equation(2) computes this law.

In the director's finite typed certificate replace c_k by cbar_k.
This is VALID: conditional on the row types, the fresh output column
permutations/signs have the same kernel regardless of the original frame.
Choose potentials by type tuple, apply the exact conditional bound, and
average over the independent frames. Products of row counts average to
products of cbar_k. Markov/union then produces some actual choice of both
frames and output gauges. No single frame is required to realize every
averaged count simultaneously. The subexponential type overhead remains
valid because every output frame is an integer sign Hadamard.

## 3. What the existing Gaussian-Fock norm argument actually gives

Let R=H2/sqrt(2). Gaussian edge defects factor as

    exp[-t||a-Rb||^2]=<phi_t(a), Gamma(R)phi_t(b)>,

where Gamma(R) is unitary on Gaussian Fock space. After averaging each
row's gauges, the full graph is a tensor contraction with these unitary
edge maps. Absorb each unitary into one endpoint leg and apply the usual
product-of-vertex-Hilbert-norm bound. Every vertex norm is unchanged.
The result is the same scalar orbit quantity L_t(w_i) as in the earlier
rank-one proof; the H2 mixing has been erased exactly.

Let T be the cross feature involution, zero on unused loop coordinates.
With ||w||^2=N, use the nonnegative defect2N-2sigma<w,Tw>. Its cross-edge
part is twice sum||a-sigma Rb||^2; unused coordinates contribute a further
nonnegative term. The scalar row-norm/recursive estimate gives the cap
coefficient expression

    1/2 + [log2+E_t(nu_1)]/(2t),                          (3)

where nu_1 is the symmetric Boolean input law and E_t is the existing
conditional-variance envelope. Choosing the label to reveal the input
shows E_t(nu_1)>=-log2. Thus (3) is NEVER below1/2. Independent recursive
children do not repair the loss after this norm bound.

This is the precise specialization of the preserved unitary-edge audit,
`decisive_bridge_vector_seed_weave_unitary_obstruction_2026_09_06.md`.
It does not claim the actual construction has a half-floor.

## 4. Even local sign/permutation projection has no strict operator gap

One might hope that coordinate-sign and pair-swap averaging compresses
Gamma(R) to a strict contraction. It does not. The invariant homogeneous
degree4 polynomial

    p(x,y)=x^4-6x^2 y^2+y^4

is unchanged by either coordinate sign and by swapping x,y, but
p(R(x,y))=-p(x,y). Radial polynomials are fixed. More generally angular
modes cos(4j theta) acquire (-1)^j. These are nonzero invariant Fock
directions with eigenvalues+1 and-1, so the compressed norm is still1.

Any improvement must retain GLOBAL cancellation between these signed
angular sectors, or control the actual annealed histogram law in (2)
inside the joint typed certificate. The old scalar Gaussian-boundary
argument controls row norms but not that cancellation. No valid closed
joint inequality retaining it has yet been derived here.

The exact-flat sector theorem remains available for every realized
recursive frame, since its proof used only power-of-two Hadamard order.
Potapov's Walsh-specific counting theorem is not transferred to these
general random Hadamard frames.

## 5. Direct substitution of the joint terminal potential is false

For a symmetric scalar source nu of variance1, consider the genuine
four-marginal transport potential

 Psi_t(nu)=-t+(1/4) sup_pi
      [2t E_pi(a^T R b)-D(pi||nu^4)],                     (4)

where all four scalar marginals of pi are nu. This is the natural
rank-two analogue of the old self-transport terminal potential. It cannot
replace that terminal potential in its uniform orbital-comparison lemma.

For nu1=(delta_-1+delta_1)/2, a^T R b is always +/-sqrt(2). The Gibbs tilt
preserves every fair scalar marginal by the kernel symmetries. Therefore

    Psi_t(nu1)=-t+(1/4)log cosh(2sqrt(2)t).               (5)

For a Gaussian source, choose jointly Gaussian a,b with cross covariance
rR. Then E a^T R b=2r and D=-log(1-r^2), giving

    Psi_t(gamma)>=-t+tr+(1/4)log(1-r^2).                 (6)

At t=4,r=7/8 the latter is -1/2+(1/4)log(15/64)>-1,
while (5) is less than -4+2sqrt(2)<-1. Thus a strict gap is proved without
relying on floating-point optimization.

Now take any sequence of normalized sign Hadamards U_k and a uniformly
random Boolean word epsilon. The symmetrized empirical coordinate law of
U_k epsilon converges in probability to the standard Gaussian. Indeed
each coordinate satisfies the ordinary bounded-summand CLT; each distinct
pair has a joint independent-Gaussian limit, uniformly since two Hadamard
rows agree in exactly half their positions after a sign gauge. Bounded
test-function variances tend to zero. The empirical second moment is
exactly1, so this is also convergence in quadratic Wasserstein distance.

To transfer the lower bound (6) to empirical laws, push its Gaussian
coupling through the same coordinatewise quantile map to each empirical
marginal. Relative entropy cannot increase under this map, while the
bilinear energy converges by L2 convergence and Cauchy--Schwarz. Hence
with probability tending to1 the empirical Psi_t is at least the bound
in (6) minus any fixed epsilon. Consequently

    E_g exp[k Psi_4(empirical(U_k g 1))]

has an exponential rate strictly larger than Psi_4(nu1). A proposed
uniform comparison by exp(o(k)) exp[k Psi_4(nu1)] is therefore FALSE,
even when U_k itself is a normalized actual Hadamard. This rules out one
specific direct reuse of the old orbital lemma, not the recursive ensemble
or the finite annealed-count certificate in Sections1--2.

## 6. An exact independent two-label Bellman obstruction

Take independent fair signs X,Y. The normalized sum and difference each
have law nu0=(1/2)delta0+(1/4)(delta_sqrt2+delta_-sqrt2), while I(X;Y)=0.
Thus a scalar joint-potential Bellman inequality would require
Psi_t(nu0)<=Psi_t(nu1). This already fails at t=4.

Here is an explicit four-marginal coupling for nu0. With probability2/3,
choose uniformly the eight positive-maximal configurations in which one
endpoint is one-hot with nonzero magnitude sqrt2 and the other endpoint
has both coordinates of magnitude sqrt2; their signs make a^T R b=2sqrt2.
With probability1/3 choose the all-zero configuration. Each scalar is
nonzero with probability(2/3)(3/4)=1/2 and has balanced signs. Hence all
four marginals are nu0. Every active atom has product-law probability1/128
and coupling probability1/12; the zero atom has product probability1/16
and coupling probability1/3. Therefore

 D=(2/3)log(32/3)+(1/3)log(16/3)=4log2-h(2/3),
 E a^T R b=4sqrt2/3,
 Psi_4(nu0)>=-4+8sqrt2/3-log2+h(2/3)/4 > -0.763.

In contrast Psi_4(nu1)=-4+(1/4)logcosh(8sqrt2)<-1.344.
This is a finite, exact, zero-information counterexample to using Psi
itself as the Bellman supersolution. It does not rule out a larger
conditional-mean-and-variance envelope.

The finite-quadrature diagnostic
`computations/flatify_construct_2026_09_07_mean_variance_probe.py` tests
the symmetric binary posterior channel E[X|L]=+/-a. Its Gaussianized
posterior mixture is nu_a=(N(a,1-a^2)+N(-a,1-a^2))/2 and its information
cost is log2-h((1+a)/2). At t=4 and quadrature order16, the largest
sampled value of Psi_4(nu_a)-I is about -0.777690 at a=0; the tested
nonzero values a=.25,.5,.75,.9,1 are all lower. This is only a diagnostic,
not an upper bound for this channel family or for arbitrary channels.
In particular it cannot establish that a conditional envelope fails:
unsampled, asymmetric, or multi-label channels may raise the envelope.

## 7. The retained Gaussian means need joint residual covariance

For any fixed pair of labels, let A=m1+sqrt(v1)G1 and
B=m2+sqrt(v2)G2 with independent standard Gaussian residuals. Applying R
gives conditional means (m1+m2,m1-m2)/sqrt2 and residual covariance

 [[(v1+v2)/2,(v1-v2)/2],[(v1-v2)/2,(v1+v2)/2]].

Thus replacing the two outputs by independent Gaussians of their
individual conditional variances loses covariance whenever v1!=v2.
This is an exact identity, not an inequality. Equal-variance binary
posterior labels avoid this particular loss, but already produce the
mean pattern {0,+/-sqrt2 a} under independent label pairing. A proposed
closure must control that label transformation as well as its information
cost. No valid inequality doing so for all labels is asserted here.
