# Gaussian regularization: exact flip-stability identity and endpoint control

This is an original-limit reduction, not a convergence proof. All sign matrices below are deterministic and are minimized **outside** Gaussian expectation. Write N=n(n−1)/2 and Q(A)=max_x |Σ_{i<j} Aij xi xj|.

## Uniform reduction

For f_n(t)=min_A E Q(A+tG)/n^(3/2), convexity and the triangle inequality give

    m_n ≤ f_n(t) ≤ m_n+t sqrt(log 2).

Indeed the signed cut family has at most 2^n distinct linear forms, each of Gaussian variance N, and E Q(G)≤sqrt(2N n log 2). Thus convergence for every fixed t>0 implies convergence of m_n by the uniform Cauchy criterion.

For β>0 define χ_e(σ,x)=σ xi xj and

    P_n(A;a,s)=(βn)^−1 E log Σ_{σ=±1,x∈{±1}^n}
                    exp[(β/√n) Σ_e(a A_e+s G_e)χ_e].

At a=1,s=t, the minimum of P differs from f_n(t) by a number in [0,(1+1/n)log(2)/β]. Convergence for all fixed β,t therefore also suffices.

## Fixed total variance and exact nonnegative local cost

Fix β,t>0, u∈[0,1], a=√u, s=√(1+t²−u), and

    F_n(u)=min_A P_n(A;a,s).

Choose a minimizer A(u) independent of the realized G. The minimum of finitely many smooth functions is locally Lipschitz; the envelope derivative exists almost everywhere and equals an active branch derivative. For A^e obtained by flipping edge e put

    Δ_e(u)=P_n(A(u)^e;a,s)−P_n(A(u);a,s) ≥0,
    S_n(u)=(4u)^−1 Σ_e Δ_e(u),       u>0.

Then, almost everywhere,

    F_n'(u)=−S_n(u)+R_n(u),
    |R_n(u)|≤β²√u/(3√n).

Here is a direct check including constants. Conditional on all Gaussian couplings, the partition sum as a function of the deterministic contribution on edge e is a positive constant times cosh(h+v), where v=β√u A_e/√n. Set m=tanh(h+v), m'=sech²(h+v). Differentiating the mean and integrating by parts in the Gaussian variance gives the contribution

    A_e E m/(2√u n√n) − β E m'/(2n²).

The contribution of −Δ_e/(4u) is

    E[log cosh(h+v)−log cosh(h−v)]/(4uβn).

Their difference has numerator 2vm−2v²m'−[log cosh(h+v)−log cosh(h−v)]. Taylor expansion at h+v bounds its absolute value by (8/3)|v|³, since |(log cosh)'''|≤2. Thus the individual error is at most (2/3)β²√u/n^(5/2); summing N≤n²/2 proves the claim. The cavity variable may depend on the Gaussian realization; the Taylor bound is pointwise, so this causes no issue.

## Uniform endpoint bound: no concentration of stability mass at zero

For every fixed A,s, a↦P_n(A;a,s) is even and convex. Evenness follows by G↦−G and σ↦−σ; convexity is log-sum-exp convexity. Its derivative for a≥0 is nonnegative. Consequently the negative Gaussian-variance derivative bounds the active branch derivative below:

    F_n'(u) ≥ −βN/(2n²) ≥ −β/4.

Together with the exact identity this gives the pointwise, dimension-uniform bound

    0 ≤ S_n(u) ≤ β/4+β²√u/(3√n).

There can therefore be no escaping delta mass near u=0, nor at any other point, for fixed β. In particular

    F_n(1)=F_n(0)−∫_0^1 S_n(u)du+ε_n,
    |ε_n|≤2β²/(9√n).

Also, directly by even convexity and the Gaussian variance derivative,

    F_n(0)−F_n(u) ≤ βu/4.

At ground state write c_n=E Q(G)/n^(3/2)≤√log 2. Jensen in the deterministic mean yields the sharper scale-sensitive endpoint bound

    c_n√(1+t²)−min_A E Q(√u A+√(1+t²−u)G)/n^(3/2)
       ≤ c_n[√(1+t²)−√(1+t²−u)]
       = c_n u/[√(1+t²)+√(1+t²−u)].

This is an upper bound on the Gaussian-reference advantage, not an absolute-value comparison.

The same bounded-density theorem also holds for the one-sided partition sum (omit σ). Evenness of its radial mean function is not generally true and is not needed: convexity remains true, and the derivative at a=0 is zero because Gaussian switching symmetry gives E〈xi xj〉=0 for every edge. Its radial derivative is therefore nonnegative for a≥0. The cavity identity only uses χ_e²=1 and is unchanged. In particular the one-sided fixed-total-variance optimized pressure is asymptotically nonincreasing with the same uniformly bounded density.

## Exact remaining lemma

The Gaussian reference F_n(0) has the usual SK pressure limit. Using the absolute partition only changes the one-sided quenched pressure by o(1): log(Z(G)+Z(−G)) lies between log Z(G) and max(log Z(G),log Z(−G))+log 2; Gaussian concentration bounds the expected excess maximum by O(√n) for fixed β,t. Hence convergence of the optimized pressure is equivalent, up to the displayed O(n^−1/2) error, to convergence of

    ∫_0^1 S_n(u)du.

The densities are uniformly bounded and nonnegative, so subsequential weak-* limits exist. What is missing is uniqueness of their total mass across subsequences. No comparison between different n follows merely from boundedness or stationarity. This isolates an exact quantitative hard lemma rather than resolving the original problem.

## Why ordinary block interpolation is still blocked

With general deterministic edge weights w_e(r) and Gaussian variances v_e(r), differentiation gives

    ∂_r P = Σ_e w'_e A_e E〈χ_e〉/(n√n)
             + β/(2n²) Σ_e v'_e E[1−〈χ_e〉²].

At an optimizer and for w_e>0, flipping an edge and Taylor expanding implies

    A_e E〈χ_e〉 ≤ (βw_e/√n)E[1−〈χ_e〉²]
                       +(4/3)β²w_e²/n.

This is an upper bound only. A child-to-parent interpolation decreases internal normalization weights √(n/n_i)→1 and increases cross weights 0→1. Reversing the path increases internal weights but decreases cross weights. Either direction contains a class with w'_e<0, for which the inequality reverses and supplies no desired upper control. Optimizing at every intermediate point does not alter this one-sided fact. Removing cross deterministic means needs a separate argument.

At w_e=0 the edge flip is identically ineffective and does not constrain the arbitrary unused sign A_e. A statement at that endpoint requires either a positive-weight limiting optimizer or an independent zero-correlation argument; the inequality above is not asserted there.

For the absolute objective there is an additional polarity issue: a parent has one σ, while the sum of absolute child caps permits independent child σ_i. Gaussian averaged one-sided pressures for A and −A need not coincide, so this cannot silently be replaced by one-sided SK subadditivity.

The applicable classical interpolation theorem is Guerra–Toninelli, *The infinite volume limit in generalized mean field disordered models*, https://arxiv.org/pdf/cond-mat/0208579, Theorem 1. Its deterministic mean is required to be a fixed finite-dimensional function of bounded additive order parameters, and its covariance a convex function of bounded additive overlaps, with uniform finite-size errors. An arbitrary optimized sign matrix does not provide the required additive mean representation. Gaussian smoothing alone does not verify that hypothesis.

## A genuine block lower comparison, and its two remaining losses

There is a separate way to remove cross means for a LOWER bound, which does not require edge stationarity. Partition the vertices into two blocks. For fixed internal means and Gaussian noise, the expected log partition is even convex in the common scalar multiplying all deterministic cross edges. Flip every spin in the first block and simultaneously flip all Gaussian cross entries to prove evenness. Therefore setting the deterministic cross means to zero decreases pressure. This works for both one-sided and absolute partition functions.

Now restrict attention to the one-sided pressure P_n^+(A;a,t), with only x in its partition sum. Once deterministic cross means are zero, ordinary Gaussian SK interpolation keeps all internal deterministic fields fixed, and gives

    P_n^+(A;1,t) ≥ Σ_i r_i P_{n_i}^+(A_i;√r_i,t),
    r_i=n_i/n.

The Gaussian covariance difference is n t²/2 times q²−Σ_i r_i q_i²≤0; its diagonal difference is zero up to the harmless common off-diagonal convention. The standard derivative consequently has the lower-bound direction displayed. Minimizing A gives the same inequality for optimized one-sided pressures, since the right side depends only on independently choosable internal matrices.

This is not almost-superadditivity of the desired normalized sequence: each child mean amplitude is √r_i, not 1. The fixed-total-variance monotonicity established above would compare the desired child endpoint if its noise were √(t²+1−r_i), but the available noise is only t. The missing variance is order one for comparable blocks.

For the absolute objective, one may retain the two one-sided expected pressures as a vector. Gaussian concentration shows that its absolute pressure differs from max(P_n^+(A),P_n^+(−A)) by O_{β,t}(n^−1/2), uniformly in deterministic A. The preceding componentwise comparison then yields a lower bound involving

    max_{σ=±1} Σ_i r_i P_{n_i}^+(σ A_i;√r_i,t).

It does not yield Σ_i r_i max_σ P_{n_i}^+(σ A_i;√r_i,t): the polarity must be common to the children. Thus even if the missing variance were supplied, the finite two-state optimization still needs a composition theorem. No equality of the one-sided pressures for A and −A was used or asserted.

## The one-sided optimized problem is actually trivial

This shows that the polarity issue is essential, not a technicality. For every a≥0,s≥0,

    P_n^+(0;0,s) ≤ min_A P_n^+(A;a,s)
                       ≤ P_n^+(0;0,s)+a/(2√n).

The lower bound is convexity and zero derivative at the zero mean, as proved above. For the upper bound take every Aij=−1. Then

    (a/√n) Σ_{i<j} Aij xi xj
       = a[n−(Σ_i xi)²]/(2√n) ≤ a√n/2.

This is a pointwise upper bound on the deterministic Hamiltonian, hence gives the displayed pressure upper bound for each Gaussian realization. It also proves the analogous ground-state bound. Thus the optimized one-sided regularized pressure always converges to the pure Gaussian pressure; it retains none of the original absolute-cap optimization difficulty. Along the fixed-total-variance path its limiting density is simply the Gaussian variance derivative. A successful original-limit interpolation must preserve simultaneous control of both polarities throughout.
