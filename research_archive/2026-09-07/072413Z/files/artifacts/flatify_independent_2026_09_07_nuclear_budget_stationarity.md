# Convex nuclear-budget surgery and a single global near-ground law

Date: 2026-09-07. Status: extension of the audited actual finite-rank surgery
construction. This strengthens the quantifiers in its stationarity conclusion;
it does not yet give a cross-order comparison.

## 1. Rank is not needed in the implementation

In `flatify_construct_2026_09_07_finite_rank_surgery_stationarity.md`, replace
G=UTU^T of rank at most r by ANY real symmetric matrix satisfying

```math
||G||op<=1,  ||G||*<=r.                              (1)
```

Full rank is allowed. The same actual-sign surgery and the same uniform
error bound hold. Here are all places in the proof that used rank:

* Replace the leverage projector UU^T by |G|. Its diagonal entries are
  nonnegative, at most 1, and sum to at most r. Mask coordinates where
  |G|_ii>mu; at most r/mu coordinates are removed.
* The mask error |x^TGx-(Zx)^TG(Zx)|<=2sqrt(n b) uses only ||G||op<=1.
* Diagonal errors are bounded by sum_i |G_ii|<=Tr|G|<=r, both in masking
  and when replacing an off-diagonal form by the full form.
* With G=sum_l lambda_l v_l v_l^T, the nonnegative PSD envelope
  K=sum_l |lambda_l| |Zv_l||Zv_l|^T has diagonal <=diag|G|<=mu on retained
  vertices. Thus the same sign-compatible flip probabilities are valid,
  and the same Schur-bias bound applies.
* The full-cube rounding variance uses only
  sum_l |lambda_l| ||Zv_l||_1^2<=n sum_l |lambda_l|<=nr.

There is no dependence on the number of nonzero eigenvalues. Consequently
for every G in (1) an actual full signing A'(G) exists with

```math
Q(A'(G)) <= Q(A-theta sqrt(n) offdiag G)+e_n n^(3/2),
```

with exactly the e_n of the cited theorem. In particular the construction
is uniform on the COMPACT CONVEX set (1), not only on a union of fixed-rank
operator balls. With alpha=r/sqrt(n), mu=alpha^(1/3)/sqrt(n),

```math
e_n <= theta(1+2 K_G C) alpha^(1/3)+theta r/(2n)
       +2sqrt(theta[(n+2)log2/n]alpha)
       +4(n+2)log2/(3n^(3/2)),                      (2)
```

provided theta alpha^(1/3)<=1 and Q(A)<=C n^(3/2).

## 2. A global law via the Ky Fan support function

Let A be a normalized delta-near-minimizer, set q=Q(A)/n^(3/2), and write
h(x)=H_A(x)/n^(3/2). Global optimality and the actual construction give

```math
min_{G satisfying (1)} max_{sigma,x}
 [sigma h(x)-(theta/2) Tr(G sigma xx^T/n)] >= q-epsilon,
epsilon=delta+e_n+theta r/(2n).                     (3)
```

For an integer 1<=r<=n, the support function of (1), on a symmetric C,
is sum_{j=1}^r s_j(C), the Ky Fan r-norm (the r largest singular values).
Indeed diagonalization and the trace inequality reduce the maximum to
choosing r eigen-directions of largest absolute eigenvalue, with matching
signs. The noninteger version linearly interpolates the last singular value.

Finite-dimensional minimax applied to (3) therefore gives ONE probability
law mu on signed spins such that, with C_mu=E_mu[sigma xx^T/n],

```math
E_mu[q-sigma h(x)]+(theta/2)||C_mu||_(r)<=epsilon.     (4)
```

Both terms are nonnegative. This law works simultaneously for ALL feature
subspaces: every rank-r orthogonal projection P obeys

```math
||P C_mu P||*<=||C_mu||_(r)<=2epsilon/theta.
```

The earlier fixed-feature theorem selected a potentially different law for
each subspace. Equation (4) removes that dependence, because the actual
sign operation extends to the convex nuclear-budget domain.

At r=1, fixed theta and exact global optimality, one common near-ground law
has ||C_mu||op=O(n^(-1/6)). Also ||C_mu||*<=1, hence
||C_mu||F<=sqrt(||C_mu||op)=O(n^(-1/12)). This is SIGNED covariance balance,
not unsigned isotropy and not necessarily equal orientation probabilities.

## 3. One-sided gap forces globally spread near-ground states

Suppose the negative cap is at most (q-gamma)n^(3/2), with gamma>0.
Equation (4) implies negative-polarity mass at most epsilon/gamma. Let B_+
be the unnormalized positive-polarity covariance E[1_{sigma=+}xx^T/n].
Since C_mu=B_+-B_- and Tr B_-<=epsilon/gamma, Ky Fan triangle inequality
gives

```math
||B_+||_(r)<=2epsilon/theta+epsilon/gamma.
```

Condition additionally on q-h(x)<=eta. Positive matrices only decrease
before renormalization; the retained probability is at least
1-epsilon/gamma-epsilon/eta. Thus the conditional positive near-ground
law has covariance B of trace 1 with

```math
||B||_(r) <= (2epsilon/theta+epsilon/gamma)
              /(1-epsilon/gamma-epsilon/eta).        (5)
```

This is a single near-ground law for which NO r-dimensional subspace
captures more than the displayed trace fraction. If that fraction is kappa
and rank(B)=d>=r, then r/d<=||B||_(r)<=kappa, so d>=r/kappa. Taking r just
below sqrt(n) makes this a quantitative near-square-root-dimensional
spreading statement, rather than escape separately from selected subspaces.

It does not exclude a one-sided gap supported on a high-dimensional family.
Nor does it bound max_x(|H_A(x)|+|v dot x|) for an inserted row. Those are
genuine remaining steps, not consequences of covariance stationarity.
