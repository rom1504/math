# Independent audit of the diffuse-variance stability penalty

2026-09-07. Full-file independent audit of `principle_director_stability_orthant_penalty_2026_09_07.md`: **PASS**. No numerical optimization or unproved high-dimensional orthant approximation is used. The bad-profile counting step remains open, exactly as stated in the director's file.

Later same-session update: the selector net argument closes that bad-profile step for the existing restricted-weave certificate. Its separate complete audit is `principle_construct_2026_09_07_selector_stability_gap_audit.md`. The discussion below records exactly what the row lemma itself proves, not the later global conclusion.

## 1. Imported product-space theorem checked

I read Talagrand's Theorem 4.1.1, Lemma 4.1.2, and the coordinate-induction setup directly in the [primary paper](https://arxiv.org/pdf/math/9406212), printed pages 50-52. The stated exponential moment is E exp(d_T(X,A)^2/4)<=1/Pr(A). I also read the arbitrary-product statement and convex-hull exponential inequality in Ledoux's [primary exposition](https://www.numdam.org/item/SPS_1999__33__120_0.pdf), Section 3.2, printed page 163. There is no requirement that the independent coordinate laws be identical or unbiased.

For a nonnegative convex K-Lipschitz violation L and a point z, a subgradient g with norm at most K gives, for every y in its zero set,

    L(z)<=g.(z-y)<=2 sum_j |g_j| 1_(z_j!=y_j).

The same nonnegative vector |g|/K works for every y; the supremum/infimum order in convex distance is therefore correct. It follows that d_T(z,{L=0})>=L(z)/(2K). The imported exponential moment, followed by Jensen twice, yields precisely

    Pr{L(X)=0}<=exp[-(E L(X))^2/(16K^2)].

The proof never applies an upper-tail convex-function theorem to a concave negative. Empty zero sets and constant functions are handled directly.

## 2. Heavy conditioning and the scalar approximation

After conditioning on heavy signs, let w_j=v_j S_j for heavy j and w_j=v_j mu_j for light j. Then ||w||_2<=sqrt(Cm). The signed row transform H/sqrt(m) is a contraction, so the conditional mean vector has norm at most sqrt(Cm)+D sqrt(k). With k>=p0 m this is at most (sqrt(C/p0)+D)sqrt(k). Hence at least k/2 coordinates have absolute mean at most the director's M.

All physical coordinates have the same light variance s^2=m^(-1)sum_J v_j^2(1-mu_j^2), because every Hadamard entry has magnitude one. It lies in [kappa,C]. Each centered summand has magnitude at most 2V/sqrt(m), so the sum of third absolute moments is at most 2VC/sqrt(m). The matched Gaussian summands have a third-moment sum of the same order, since each standard deviation is at most V/sqrt(m).

The smoothed hinge (sqrt(z^2+delta^2)-z)/2 has uniform error at most delta/2 and third derivative bounded by a universal multiple of delta^(-2). Thus independent Lindeberg replacement gives an additive expectation error tending to zero uniformly in the conditional means, once delta is fixed. This is a one-coordinate Lipschitz expectation approximation, not an estimate of an exponentially small joint probability.

For N(b,s^2), E(-N)_+=s phi(b/s)-b Phi(-b/s). It decreases with b and increases with s. Hence b<=M and s>=sqrt(kappa) give the positive lower bound g in the director's theorem. At least k/2 coordinates therefore contribute at least g/2 each for large m, even after arbitrary heavy conditioning.

## 3. Convex violation and exact constants

For L=sum_a(-Z_a)_+, the outer hinge sum is sqrt(k)-Lipschitz. The light-coordinate linear map has norm at most V, since H/sqrt(m) is a contraction. Thus L has Lipschitz constant at most V sqrt(k), while its conditional mean is at least gk/4. Substitution in Section 1 gives

    Pr(all Z_a>=0 | heavy signs)
       <=exp[-g^2 k/(256V^2)].

Every bound was uniform in the conditioned signs, so averaging them is legitimate.

## 4. Actual weave, diagonal, and pressure payment

Writing C=W-diag(W), direct differentiation of the hollow quadratic energy gives

    sigma x_i(a)(Cx)_i(a)/sqrt(mk)
      =sigma x_i(a)sum_(j!=i) H_i(a,j)u_j(i)S_ij/sqrt(m)
        +sigma S_ii x_i(a)H_i(a,i)u_i(i)/sqrt(m)
        -sigma S_ii/sqrt(mk).

This is the stated field, including the physical diagonal subtraction. Its deterministic term has norm at most sqrt(k)|u_i(i)|/sqrt(m)+1/sqrt(m)<=2sqrt(k), since |u_i(i)|<=sqrt(k) and k<=m.

Under the fixed-spin tilt the independent edge means are tanh(2t sigma u_i(j)u_j(i)). Each independent edge enters exactly two stability events, so graph Holder assigns exponent 1/2 to each row probability. For theta m covered rows the added factor is exp[-theta g^2 mk/(512V^2)]. Dividing by m^2 gives the penalty theta p g^2/(512V^2).

The Markov threshold factor for Q(C)>=c(mk)^(3/2) is exp[-2tc sqrt(p)m^2]. Therefore the director's final conditional cap coefficient (K-Delta)/(2t sqrt(p)) has the correct normalization. This consequence still REQUIRES the explicitly stated strict bound on the spin weight of bad profiles. The row lemma alone does not establish such a bound.

Minor notational suggestion sent to the director: use C, rather than the pre-hollowing W, in the Section 3 tilt. Keeping the two matrices distinct makes the already correct physical diagonal accounting unambiguous.
