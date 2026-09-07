# Fresh reconstruction of the <0.494515125 all-order upper chain

Date: 2026-09-07. Status: PASS. This reconstruction uses the mathematical
proofs, not their earlier audit verdicts. It needs neither H=E nor a lower
conditional-copy construction. The original convergence question stays open.

## 1. The minimal analytic closure

Fix t>0. Write F for entropic Gaussian self-transport, Phi=-F/2,
G(nu)=g_t(m_2(nu)), and E(nu)=sup_L[g_t(E Var(X|L))-I(X;L)].
The precise source proofs are the standalone reconstruction Sections 1--3,
terminal-gap reduction Sections 1--2, and precision/Schur Sections 1--2.

The scalar precision proof gives BE<=E without optimizer attainment:
positive child precisions lambda_1,lambda_2 yield arithmetic/harmonic
parent precisions a,h in (0,t], ah=lambda_1 lambda_2. Sequential parent
labels (B,M,N) and (M,N) have total information at most the child total
plus I(A;B). Completing the weighted square bounds their distortion.
Concavity of c_t(exp(s)) then compares c_t(a)+c_t(h) with the child
sum. Approximate child channels and independent child precision suprema
complete the proof. No singular covariance inversion occurs here.

The bounds G<=E<=Phi follow from the constant label, Gaussian extremality
of F, and conditional self-coupling. Gaussian extremality for an asymmetric
conditional source follows by centering and symmetrizing first. It is not
restricted to symmetric posteriors. Gaussian sources satisfy G=E=Phi.

The unbounded strictness proof is substantive: finite scaling measures
mu_j have Gaussian energy one and uniformly bounded mass on every unit
cube. Vague compactness and uniformly summable Gaussian tails pass
alpha_j=(k mu_j)mu_j to alpha=(k mu)mu. This also prevents a zero limit.
Kernel Cauchy--Schwarz gives 0<k mu<=1; retaining one positive-mass ball
bounds -log(k mu) by a quadratic. The resulting potential is integrable.
Entropy comparison therefore identifies the exact self-coupling. Strict
positivity of Gaussian energy for nonzero signed differences follows from
Gaussian convolution and its nonvanishing Fourier multiplier, applied to
tempered distributions. Auxiliary measures need not have finite total mass.
The log-sum inequality plus strict energy convexity proves strict source
concavity of F on finite-second-moment laws.

Self-transport tensorization and strict concavity make product sources
the UNIQUE equality cases in F(X,Y)<=F(X)+F(Y). For a safe Bellman pair,
the information identity in the standalone proof gives actual Phi drift
at least I(A;B)/4. An attained zero-drift pair must consequently be iid,
and its sum and difference independent. The characteristic-function
identity and the finite-variance expansion at zero identify a centered
Gaussian, including variance zero.

Moment escape is controlled through costs, not through continuity of the
second moment. Both F and each rate-distortion J_lambda obey the binary
mixture sandwich, while J_lambda<=F_t<=2K_t(m_2), K_t=-g_t=O(log(1+v)).
A tail of mass q in a moment-C ball thus changes the cost by at most
h(q)+2q[K_t(C/q)+K_t(C/(1-q))], uniformly in lambda<=t. Compact-source
Wasserstein continuity then gives weak continuity of Phi and E on the
whole moment ball. Compactness of safe coupling sets and lower
semicontinuity of mutual information give lower semicontinuity of
Gamma=Phi-BPhi and attained maximizing pairs.

Let f_r=B^r Phi, D=Phi-E, E0=m_2(nu0), and B0=Phi(nu0)-G(nu0).
For epsilon>0 and C>E0, compactness and Gaussian rigidity give positive
kappa=inf{Gamma: m_2<=C, D>=epsilon}, unless that set is empty.
For an approximately optimal FINITE depth-r policy, its total expected
actual Phi drift is at most B0+zeta, because f_r>=G. Follow a uniform
branch and stop at D<=epsilon, moment>C, or depth r. The moment is an
exact nonnegative martingale, and the unexited probability is at most
(B0+zeta)/(r kappa). Iterating BE<=E only to that bounded stopping time
gives, after zeta decreases to zero,

    f_r(nu0)-E(nu0)
      <=epsilon+E0 sup_{s>C} K_t(s)/s+K_t(C)B0/(r kappa).

Choose C, then epsilon, then finite r. This proves limsup_r f_r<=E.
The budget is Phi-G, NOT Phi-E. No lower comparison E<=f_r was used.

## 2. Fixed-depth row bound, including arbitrarily large coordinates

For the signed permutation orbit define P_t(v)=E_g exp(-t||v-gv||^2)
and L_t=sqrt(P_t). Its Gaussian Fock feature orbit covariance has norm
exactly P_t(v): its finite orbit Gram matrix is nonnegative with that
constant row sum. At degrees <=D, the invariant rank is bounded by
sum_{j<=D/2}p(j)<=exp(pi sqrt(D/3)). The remaining squared feature norm
is a Poisson tail. If ||v||^2<=Cm and
D=ceil(e^2(2tC+1)m), that tail is at most exp(-(4tC+1)m)<=P_t(v).
Hence for EVERY orthogonal U, uniformly in the vector and in U,

    E_g L_t(Ugv)<=exp(O_{t,C}(sqrt(m))) L_t(v).

This directly covers a coordinate carrying order-m energy; no coordinate
truncation or delocalization assumption is being substituted.

The full signed-permutation invariant space lies inside the block-product
invariant space. Gaussian Fock factorization therefore gives
L(v_+,v_-)<=L(v_+)L(v_-). Also the permanent terms fixing coordinate i
give P(v)>=K_t(v_i,v_i)P(v without i)/m, and K_t(a,a)>=1/2. Thus

    max_i L(v without i)<=sqrt(2m) L(v).

This diagonal deletion costs only a polynomial factor, even for a large
deleted coordinate or repeated spectrum values.

At depth r, exact signed-word counting gives a pair-table probability
exp[-m D(pi||nu tensor nu)/2+O(log m)]. The constraint is the average
ABSOLUTE marginal. Symmetric reference probabilities make the displayed
relative entropy valid even if the pair marginals themselves are unequal.
Simultaneous reversal and swap then make both signed marginals nu without
changing either symmetrized child law and only decreasing information.

At FIXED r, all alphabets and numbers of node types are finite/polynomial.
Every terminal energy is bounded by the root energy; its moment per
terminal dimension is at most 2^r C. There are only 2^r orbital estimates.
Their total log loss and all table errors are o(m), uniformly over the
terminal Hadamards. Consequently

    limsup_m m^-1 log E L_t(U_m v_m)<=B^r Phi_t(nu).

Changing k/m to its limit p causes no alphabet issue in the application:
the root nonzero amplitude is 1/sqrt(k/m), and every fixed-depth amplitude
is a fixed finite linear combination of this common scale.

## 3. Actual signing, all spins, both polarities

Use independent recursive bases H_i=sqrt(m)U_i^T and fresh output-column
permutations in each fibre, and independent fair off-diagonal S_ij=S_ji.
For k retained rows per fibre the full sign matrix W_T has order N=mk.
With h_i=H_i[T_i,:]^T x_i, the exact identities are

    sum_i ||h_i||^2=m^2 k,
    D_sigma=2(m^2 k-sigma x^T W_T x).

Each undirected edge appears twice in D_sigma. Averaging its sign in
exp[-tD_sigma/(2k)] therefore gives EXACTLY the folded Gaussian kernel
at h_i(j)/sqrt(k),h_j(i)/sqrt(k), with no extra factor two. Drop only
nonnegative diagonal defects. Conditional on the diagonal absolute value,
the remaining row entries are uniformly permuted multisets. PSD Gram
factorization and graph Cauchy--Schwarz bound the moment by the product
of their L_t values. Repeated values retain the same factorial quotient.

Define Z_i as the sum over ALL 2^k row spins of the maximum deleted-row
L_t. Markov and the union over all full spins AND sigma=+/-1 give

    P(exists x,sigma: D_sigma<=2 gamma m^2 k)
      <=2 exp(t gamma m^2) product_i Z_i.

The normalized input vector is (1_T x)/sqrt(k/m), with second moment
exactly one. Thus limsup m^-1 log E Z_i<=p log2+B^r Phi_t(nu_p).
Independent fibres give (E Z_i)^m, not E[Z_i^m]. Any fixed selectors of
size k suffice; selectors are not optimized after seeing the bases.

For t=4,p=31/32, the certified phase E_4(nu_p)=g_4(1) and Section 1
allow any fixed strict exponent margin at some finite depth. The resulting
hollow cap, paying at most N/2 to delete the matrix diagonal, satisfies

    limsup Q/N^(3/2)<=(4+p log2+g_4(1))/(8 sqrt(p)).

## 4. All orders and the numerical endpoint

The explicit H_12=I+R from the quadratic character of F_11 has RR^T=11I
and R^T=-R, so it is Hadamard. H_2 and H_12 supply terminal orders
2^a 12^b. Since log12/log2 is irrational, finitely many nonnegative b
give arbitrarily small circular gaps of their residues. For every large
target logarithm, a nonnegative a then places an available order just
above it with relative error tending to zero. Multiplying by fixed 2^r
preserves this property, as does N=m floor(pm).

The terminal matrices may themselves have many tensor factors: this does
NOT introduce growing recursive depth into the row theorem, because its
orbital estimate is uniform in the arbitrary terminal orthogonal matrix.
Principal monotonicity Q(A_I)<=Q(A), proved by averaging omitted unbiased
spins, fills every remaining order with normalized multiplicative loss
(N/n)^(3/2)=1+o(1).

The complete limit order is: strict cap margin; finite r; all sufficiently
large admissible m and then every n; finally margin to zero. There is no
growing-depth type theorem, spectral-tail substitution, H=E premise, or
assumption about actual optimizing seeds.

The phase support and weight reductions were freshly reconstructed:
the KKT Gaussian mixture has at most one positive local maximum, source
endpoints cannot be maxima, and the weight transformation covers [0,1]
bijectively. Low precisions lambda<=1/2 are handled by unit subgaussianity;
high precisions lie in the closed certified rectangle. Zero amplitude
and singular optimizer-weight formulas are handled by their endpoint
branches. The interval exponential/logarithm implementations were inspected
in the installed mpmath1.3.0 and use downward/upward rounding separately.

The exact rational endpoint conversion was rerun and gives

    (4+p log2+g_4(1))/(8 sqrt(p))
      <=7787631971809/15748015748016<0.494515125.

The denominator uses a LOWER square-root endpoint because the numerator
is positive. The standalone exact checker was also rerun: 240 signed words,
50 pair tables, 24 permanent projections, the repeated-value graph test,
all 4096 spins of its finite weave, and 138 interval primitives passed.
The full high-precision rectangle verifier was freshly rerun with `--verify`:
425009 tree nodes, 212505 accepted rectangles, strict threshold -0.777671,
and Gaussian lower endpoint -622136276211/800000000000; it completed with
`status="directed interval exclusion"` in 93.43 seconds. The floating
diagnostic maximum is not used as a certified numerical endpoint.

No material mathematical correction was found. Older dependency text that
requires H=E for this upper bound can be shortened to the direct-E closure
in Section 1; no claimed lower construction is needed.
