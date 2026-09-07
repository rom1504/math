# Uniform quadratic fluctuations and an extensive actual bridge penalty

Status: complete elementary proof. This is a genuine logarithmic-pressure comparison, not an inference from variance at the zero-bridge endpoint. It applies to every bounded-operator interaction, without assuming global optimality or mixing. It does not settle the contracted-child-temperature or common-polarity issues.

## 1. Product-mixture quadratic covariance

Let J be hollow symmetric, ||J||op≤L, and consider the joint law

    μ(s,x) ∝ exp[s Σ_{i<j} Jij xi xj + b s],

where b is arbitrary. A quenched scalar b=hg can subsequently be averaged; every bound below is uniform in b. Put K_s=L I+sJ≥0. Given (s,x), sample

    Y=K_s x+K_s^(1/2)G,       G standard Gaussian.

Conditional on (s,Y), the coordinates X_i are independent with means tanh Y_i. To verify this, multiply the Ising density exp(x^T K_s x/2), whose added diagonal is constant on the cube, by the Gaussian channel density. The quadratic terms cancel and leave exp(Y^T x). Singular K_s follows by adding εI and taking a limit.

For any deterministic hollow quadratic coefficients b_ij, write F=Σ_{i<j}b_ij XiXj. Under a product law, expand X_i=m_i+Z_i into centered independent coordinates. The degree-two terms b_ij Z_iZ_j are pairwise orthogonal and orthogonal to the constant and linear terms. Therefore

    Var(F | s,Y) ≥ Σ_{i<j} b_ij² sech²Y_i sech²Y_j.

The same inequality holds for sF. The law of total variance and sech²y≥exp(−2|y|) give

    Var(F), Var(sF) ≥ Σ_{i<j} b_ij² a_i a_j,
    a_i=exp[−2 E|Y_i|].                              (1)

Indeed Jensen applied to exp[−2(|Y_i|+|Y_j|)] yields a_i a_j; no independence of the latent fields is assumed.

This product-mixture mechanism is standard; a recent primary use for Ising anti-concentration is Daskalakis–Kandiros–Yao, *Estimating Ising Models in Total Variation Distance*, https://arxiv.org/pdf/2511.21008, §6.2. Their statistical theorems impose additional hypotheses; none of those theorems is invoked here. The identity and the present bound have been proved directly above.

## 2. A uniform weight budget on every macroscopic subset

Since diag(K_s)=L I,

    E|Y_i| ≤ L+E|(JX)_i|+√(2L/π).

For U⊂[n] of size rn, Cauchy–Schwarz and the operator bound give, deterministically for every spin vector,

    Σ_{i∈U}|(Jx)_i| ≤√|U| ||Jx||₂≤L√(n|U|).

Consequently Jensen supplies the dimension-uniform subset bound

    Σ_{i∈U} a_i ≥ |U| exp{−2[L+L/√r+√(2L/π)]}.       (2)

This does not assert a uniform lower bound on every a_i. It is enough for a flat macroscopic rectangular coefficient pattern.

For a bipartition U,V of proportions r,1−r and coefficients b_ij of absolute value β/√n on every cross edge, (1)–(2) imply

    Var(sF) ≥ κ(L,r) β²r(1−r)n,                      (3)

where

    κ(L,r)=exp{−2[2L+L(r^(−1/2)+(1−r)^(−1/2))
                         +2√(2L/π)]}>0.

No asymptotic approximation appears in (3).

## 3. The lower bound holds throughout the actual bridge tilt

Now suppose J=βA/√n for an actual full sign matrix A and ||J||op≤L. Let J_0 retain only the two diagonal blocks, and let B=J−J_0 be the cross bridge. For 0≤α≤1 put

    J_α=J_0+αB=(1−α)J_0+αJ.

Each principal block of J has norm at most L, hence ||J_0||op≤L and ||J_α||op≤L throughout. The *same* lower bound (3) therefore applies under every tilted Gibbs law along this path, with F=Σ_cross Jij xi xj.

Define ψ_b(α)=log Σ_{s,x}exp[s Σ_{i<j}(J_α)ij xi xj+b s]. Then exactly

    ψ_b''(α)=Var_α(sF)≥κ(L,r)β²r(1−r)n.

Switching all spins in U shows that ψ_b is even in α, so ψ_b'(0)=0. Twice integrating gives the extensive comparison

    log Z_b(J)−log Z_b(J_0)
       ≥ [κ(L,r)/2] β²r(1−r)n.                     (4)

The estimate holds for every b, and therefore after averaging any quenched orientation field. It is uniform over every actual parent signing satisfying the operator bound. There is no sign-reoptimization switch error: each parent is held fixed along its bridge tilt, and the comparison is established for all such parents separately before taking any minimum.

## 4. Exact scope of the gain

Equation (4) upgrades a zero-bridge variance lower bound to a genuine extensive pressure penalty by proving the fluctuation estimate under all intermediate tilts. Covariance of X alone would not justify this step; the latent product decomposition and a pathwise operator bound are essential to the argument given here.

The endpoint J_0 still has internal coefficients β/√n, so child temperatures are contracted by √r and √(1−r). Its two-sided partition has one shared energy orientation, not independently chosen orientations for the children. Moreover κ(L,r) is a possibly very small positive constant; it has not been shown to match the exact normalization loss. Thus (4) is a new quantitative bridge inequality, not an almost-superadditive optimized-pressure theorem.

The fixed-temperature spectral-refill theorem in decisive_audit_regularization_flip_stability_2026_09_07.md, §5, separately permits approximation by bounded-operator actual signings. Its O(1/L) loss cannot silently be absorbed into the exponentially small κ(L,r) in (4). No such absorption or optimized recurrence is asserted.
