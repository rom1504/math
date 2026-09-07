# Quadratic fluctuations and pressure curvature from an actual cap bound

Status: proved; the Gaussian product-mixture argument was independently
reconstructed by both `decisive_bridge` and `transfer_seeds`. This strengthens
the preceding optimizer covariance theorem: optimality is not needed for
the aggregate assertion below. It is not a convergence theorem.

## 1. General finite theorem

Let J be ANY real symmetric hollow N by N matrix. Put

    H_J(x)=sum_(i<j) Jij xi xj,       Q(J)=max_x |H_J(x)|.

Consider the Gibbs law proportional to exp[s H_J(x)+hgs], where s=+-1 and
g is standard Gaussian, and average over g after normalization. The same
conclusion holds with a fixed s or with any externally randomized weights
of the two branches. There are no external fields acting directly on x.

Choose D=diag(d_i)>=+-J, and write

    q=Q(J)/N,   d=tr(D)/N,
    K=min{4q,2q+1}+d+sqrt(2d/pi).

For arbitrary real hollow quadratic coefficients b=(bij), let

    V_b=sum_(i<j) bij^2,
    eta_b=N max_i sum_(j!=i) bij^2 / V_b,       if V_b>0.

Then, with variance taken in the joint averaged law,

    Var(H_b(x)) >= V_b exp(-2 eta_b K),
    Var(s H_b(x)) >= V_b exp(-2 eta_b K).                    (1)

In fact (1) holds for E_g Var_(mu_g), not merely for the variance after
mixing g. This distinction is needed for the quenched pressure consequence.

A valid simultaneous diagonal majorant satisfies d<=4 K_G q, by the
Grothendieck/SDP proof reconstructed in
`transfer_director_sparse_restriction_audit_2026_09_06.md`, Section 1.
Thus bounded normalized actual cap and bounded eta_b give a positive
dimension-free multiple of the full squared coefficient norm. No bound on
||J||op, high-temperature condition, or sign assumption on b is required.

## 2. Exact Gaussian product mixture

For fixed s the matrix K_s=D+sJ is PSD. The identity

    exp(x^T K_s x/2)=E_(G standard Gaussian) exp(x^T K_s^(1/2)G)

represents the Ising law as a mixture of product laws. Diagonal addition
contributes only tr(D)/2, independent of x and s, so it does not change the
original branch probabilities. A concrete joint realization is:

1. draw (g,s,x) from the stated original law;
2. draw Y conditional on (g,s,x) as N(K_s x,K_s).

Conditional on (g,s,Y), the coordinates xi are independent with means
tanh Yi and variances sech^2 Yi. This also holds when K_s is singular, by
using its PSD square root or by a vanishing diagonal regularization.

Let y_i=E|Yi|, averaged over the whole joint law. Since diag K_s=d_i,

    sum_i y_i <= tr(D)+E||Jx||_1+sqrt(2N tr(D)/pi).           (2)

There are two independent bounds on its middle term. Pointwise cube
polarization gives ||Jx||_1<=4Q(J). Alternatively, writing h_i=(Jx)_i and
conditioning on xi gives

    E[s xi h_i]=E[h_i tanh h_i],
    sum_i E[h_i tanh h_i]=2E[s H_J(x)]<=2Q(J).

As |u|<=u tanh u+1, this yields E||Jx||_1<=2Q(J)+N. Therefore

    sum_i y_i<=NK.                                         (3)

Both bounds remain true for arbitrary branch weights because the conditional
law of xi, given s and the other spins, is unchanged.

## 3. Product chaos, Jensen, and coefficient delocalization

In the conditional product law expand xi=mi+zi, with centered independent zi.
The degree-two part of H_b is sum bij zi zj. It is orthogonal to all linear
and constant terms, and its summands are pairwise orthogonal. Consequently

    Var(H_b |g,s,Y)>=sum_(i<j) bij^2 sech^2 Yi sech^2 Yj.

The same formula holds for s H_b because s is conditioned. Variance
decomposition, sech^2 u>=exp(-2|u|), and Jensen give

    E_g Var_(mu_g)(s H_b)
      >=sum_(i<j) bij^2 exp[-2(y_i+y_j)]
      >=V_b exp[-2 sum_i (sum_j bij^2)y_i/V_b]
      >=V_b exp(-2 eta_b K).

This proves every version of (1). Correlation among the latent fields Yi is
allowed; Jensen, not independence of those fields, is used.

For a flat sign bridge between blocks m and n, N=m+n, V_b=mn and
eta_b=N/min(m,n). At comparable splits it is bounded. Thus every actual
bounded-cap full Gibbs law has flat-bridge variance Omega(mn), including
after the bridge has been turned on. This improves the preceding result
which only treated a product of two child laws.

## 4. Genuine pressure consequence and the optimizer-switch limitation

Fix matrices J_0,B and an interval [a,b]. Suppose Q(J_0+tB)<=CN throughout
the interval and eta_B<=eta. Define the QUENCHED pressure

    p(t)=E_g log sum_(s,x) exp[s H_(J_0+tB)(x)+hgs].

Differentiating this fixed finite sum gives

    p''(t)=E_g Var_(mu_g,t)(s H_B)>=c(C,eta) V_B,            (4)

where c(C,eta)=exp[-2eta(min{4C,2C+1}+4K_G C+
sqrt(8K_G C/pi))]. Therefore

    p(b)>=p(a)+(b-a)p'(a)+c(C,eta)V_B(b-a)^2/2.             (5)

For a cross-block B, a block-diagonal J_0, and a=0, switching one block
makes p(t) even. Thus p'(0)=0. When B is a flat bridge scaled as beta/sqrt N,
V_B=beta^2 mn/N, and (5) is a TRUE order-N pressure increase for comparable
blocks, uniformly over bounded-cap paths. No variance-to-exponential-tail
leap is being made: the variance lower bound is valid along the entire path.

This still does not prove the desired normalized cross-order recurrence.
First, child temperatures are contracted at the block-diagonal endpoint;
the gained coefficient is not shown to pay their temperature change.
Second, if signs are reoptimized with t, the minimum is an envelope of
smooth convex branches. Downward jumps of its derivative at optimizer
switches can offset the branchwise curvature. Equation (4) is a statement
on every fixed branch, not a lower distributional-curvature bound for that
minimum. Any such envelope claim needs a separate proof.

## 5. Interpretation

The proof combines classical Gaussian augmentation and orthogonal product
chaos with an actual-cap-controlled diagonal majorant. The resulting
delocalized quadratic anti-concentration and interacting pressure bound do
not assume the answer as a Gibbs-state hypothesis. They apply to dense
quadratic systems beyond exact sign minimizers. They do not identify the
original limiting constant or establish its existence.
