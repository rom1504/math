# Independent audit of uniform almost-stability entropy loss

2026-09-07. **PASS, including total violation, arbitrary fields, and removal
of the operator bound.** Audited source:
`principle_director_universal_stability_entropy_2026_09_07.md`.

The theorem counts, for every full signing with Q(A)<=Cn^(3/2), all
spins whose total unfavorable local-field magnitude is at most
n^(3/2)/1024. It is uniform over every deterministic external field.
I reconstructed the proof independently rather than only checking its
exact-stability specialization.

In the normalized fixed-op core, the many-large-field case is correct:
every mismatched coordinate either consumes at least T/2 of violation
or at least T^2/4 of squared internal field. With T=128L these give
k/4096+k/(32768L)<k/1024 total mismatches. The Hamming-ball count on
at least k/4 prescribed signs gives a fixed exponential loss.

In the few-large-field case, conditioning on those spins leaves a field
of l2 norm at most 129L sqrt(k), independently of their values. The
convex l1 field norm has Lipschitz norm at most L sqrt(k), and symmetry
plus the sharp-enough elementary Khintchine lower bound gives mean at
least k/8. Applying the previously checked convex-zero inequality to
(F-k/16)_+ gives the stated lower tail with exponent k/(4096L^2).

The exact total-violation identity implies, outside that lower tail,
that the quadratic plus linear terms exceed k/32. The linear MGF uses
the full conditioned field norm, so it does not accidentally omit the
large coordinates' contribution. The quadratic form is centered because
the core is hollow; its Frobenius square is at most k and its operator
norm at most L. I checked the convention in
[Rudelson--Vershynin, Theorem 1.1](https://arxiv.org/pdf/1306.2872),
which gives precisely the two required quadratic tail scales, up to its
universal constant.

The simultaneous two-sided diagonal majorant normalization also checks:
the symmetric block SDP dual has objective (Tr D_u+Tr D_v)/2; swapping
and averaging gives D=(D_u+D_v)/2 with D>=+/-B and TrD<=8Cn. Thus a
core of at least 3n/4 coordinates has norm at most 32max(C,1), and its
conditional exterior field is arbitrary but covered by the proved core
lemma. Dividing the original total-violation threshold by sqrt(n) gives
n/1024<=|S|/512 as required. No full-operator or exterior-field bound
is used at this last step.

Consequently the output-entropy obstruction is valid for ACTUAL low-cap
full signings, not just the earlier weighted bipartite involution. The
scope is one bound for each arbitrary field, not a small common stable
set over all possible fields. It is a structural entropy theorem, not
a convergence proof or an entropy-preserving cleanup kernel.

## Sharp threshold and quantitative extension: PASS

I subsequently read the completed Sections 5--6 of the director's source.
The arbitrary-field lower-tail theorem holds for every fixed
v<1/sqrt(2pi). Its limit order is correct: choose the discarded and heavy
fractions first so that sqrt(2/pi) rho^(3/2)>2v, then a sufficiently large
fixed heavy-field cutoff, and finally n. The many-heavy mismatch count
uses the full total-violation budget vn. In the remaining case the exact
Rademacher mean, rather than a joint Gaussian approximation, leaves a
strict positive gap. The convex, linear, and hollow-quadratic estimates
are uniform under both conditionings. The exact zero-field mean proves
the claimed endpoint optimality, but does not assert an exponential
bound at the endpoint itself.

The sharper zero-original-field argument also passes. From
B=D^(1/2) T D^(1/2), ||T||op<=1, one has
B_S D_S^(-1) B_S<=D_S. Hence for every l1 subgradient s in [-1,1]^S,

    ||B_S s||^2 <= dmax s^T B_S D_S^(-1) B_S s
                 <= dmax Tr D_S.

The same contraction gives the conditional exterior-field norm bound
||B_(S,S^c) x_(S^c)||^2<=dmax Tr D_(S^c). Neither estimate needs a
global operator bound or bounds on the discarded diagonal entries.
Taking epsilon proportional to eta gives the exponent eta^3 n/C0^2.
The O(1/n) exact-mean error matters only when n eta is bounded, where
the enlarged prefactor makes the asserted tail estimate trivial.
Integrating the negative tail and using the exact universal mean proves
the stated L1 rate. Arbitrary added external fields are correctly
excluded from this quantitative norm shortcut.

## Independent fourth-moment variance complement: PASS

I also read and independently reconstructed
`principle_construct_2026_09_07_instability_variance.md`. Convexity gives
Delta_j F<=2x_j(A^T sign(Ax))_j for F=||Ax||_1, and flip symmetry gives
E Delta_j^2=2 E(Delta_j)_+^2. Thus the cube Poincare constant yields
Var F<=2 E||A^T sign(Ax)||^2 exactly; nonsmooth boundary errors are not
being summed in this step.

For distinct hollow rows, removing their two unmatched endpoint signs
costs O(n^(-1/2)) in sign covariance by the central-binomial bound.
The agreement/disagreement decomposition leaves two independent
one-dimensional sign sums. Scalar Berry--Esseen, the absolute-value
CDF replacement, and an O(n^(-1/2)) tie estimate give
|K_ij|<=C(|(A^2)_ij|/n+n^(-1/2)). Frobenius Cauchy--Schwarz then yields

    Var I_A <= C[n^2+Tr(A^4)/n].

The previously audited actual-sign Schatten-four bound consequently
gives the stronger zero-field L2 rate O_C(n^(-1/2)) for bounded-cap
signings, and convergence in L2 whenever Q(A)=o(n^2). Its tails are
only polynomial; it does not supersede the director's exponential
lower-deviation theorem.
