# Codewise field universality after exact new-spin maximization

2026-09-17. **Verified:** director derivation and independent discrepancy
and Bernoulli reconstructions of the finite argument. This
improves the preceding fourth-mass frame theorem on low-entropy old
codes. It uses scalar Stein replacement only AFTER retaining the full
new-child maximum, plus one empirical-frame event. No external novelty
claim is established.

## 1. Finite theorem

Let n=kp+ell, 0<=ell<k, p>=1. Let A_n,B_q be arbitrary full
sign children, and let C be any nonempty old-spin code with H=log|C|.
Suppose q>=3k log(4k), and put delta=sqrt(3k log(4k)/q)<=1.
There exist deterministic labels h_1,...,h_q in {+-1}^k such that

```math
 \sum_jh_jh_j^T\preceq q(1+\delta)I_k,                 \tag{1}
 R:=\sup_{x\in C}\sum_j\max\{1,\max_{t\le p}
                  |h_j^Tx_{:,t}|\}
 \le q[1+\sqrt{2k\log(2p)}]
           +\sqrt{2qk[H+\log4]}.                     \tag{2}
```

Construct an exact sign bridge with core columns h_j tensor g_j,
where the g_(j,t) are independent fair signs, and independent fair
signs on all leftover edges. Let Q_C^R be the parent absolute cap
restricted ONLY by x in C. Let Q_C^G be the same maximum after all
scalar bridge drivers are replaced by independent standard Gaussians.
The internal children are unchanged. Then

```math
 |E Q_C^R-E Q_C^G|
 \le R+\sqrt{2q(1+\delta)n\log(2|C|)}.               \tag{3}
```

All new-spin assignments and both polarities occur in BOTH maxima.
The labels may depend on C but not on the subsequent bridge drivers.
No favorable value of the Gaussian maximum is asserted.

## 2. One empirical list works for the whole code

Sample h_j independently uniformly from {+-1}^k. The elementary
rank-one matrix Chernoff proof in
[critical realization, Section 1](paper_localization_critical_realization_2026_09_17.md)
gives (1) with failure probability at most1/4.

For each fixed x and j, each h_j dot x_col is a k-sign sum. The
usual exponential-moment maximum bound gives
E max_t|h_j dot x_col|<=sqrt(2k log(2p)). Thus the mean of the
sum in (2) is at most its stated first term. Changing one of its qk
independent mode bits changes the sum by at most2. Bounded differences
gives centered MGF proxy qk and hence deviation
sqrt(2qk(H+log4)) with failure at most exp(-H)/4. Union over C
proves (2) with failure at most1/4. The two events intersect with
probability at least1/2. Select one list in their intersection.

## 3. Replace scalar fields, not the selected old word

Fix x and polarity s. The exact inner maximum is

```math
 F_{x,s}=sH_A(x)+\max_{y\in\{\pm1\}^q}
       \{sH_B(y)+s\sum_jv_j(x)y_j\},
 \qquad v_j(x)=\sum_t(h_j^Tx_{:,t})\epsilon_{j,t}
                         +\sum_{i\le\ell}x_i\xi_{j,i}.
```

The fields v_j are independent across j for this FIXED x. Conditional
on all other fields, the inner maximum is 1-Lipschitz in v_j; it
retains the complete new-child optimization and its offsets.
The reconstructed scalar Stein bound says, for real a_l and v=sum a_l^2,

```math
 d_W(\sum_l a_l\epsilon_l,\sqrt v\,G)
    \le\frac{\sum_l|a_l|^3}{v}\le\max_l|a_l|.       \tag{4}
```

At v=0 the field is identically zero and the replacement costs zero.
Applying (4) sequentially to the independent fields gives
|E F_(x,s)^R-E F_(x,s)^G|<=R, uniformly in x,s. This comparison
is NOT applied conditionally to an optimizer chosen using the fields.

## 4. Pay only the OLD code for selecting its best word

For fixed x,s, changing one scalar sign driver changes F by at most
twice its coefficient magnitude. Its centered MGF proxy is consequently

```math
 V(x)=\sum_{j,t}|h_j^Tx_{:,t}|^2+q\ell
       \le q(1+\delta)n=:V_* .                      \tag{5}
```

The Gaussian F is sqrt(V(x))-Lipschitz in its independent standard
Gaussian drivers, so it has the SAME centered MGF bound. These facts
hold even at ties, because a maximum of linear functions with the same
coefficient norm is Lipschitz. No independence across different x,s is
needed. For EITHER law,

```math
 0\le E\max_{x\in C,s}F_{x,s}-\max_{x\in C,s}E F_{x,s}
       \le\sqrt{2V_*\log(2|C|)}.                    \tag{6}
```

Combine (4) and (6) in each direction. Only ONE copy of the latter
error is needed in each direction, since the other expected maximum
is at least its largest mean. This proves (3).

## 5. Critical-scale consequence for actual regularized minimizers

For q=Theta(n), k=Theta(sqrt n), (3) is

```math
 |E Q_C^R-E Q_C^G|
  =O(n^{5/4}\sqrt{\log n}+n\sqrt{H+1}).             \tag{7}
```

Thus EVERY exp(o(n)) old code permits critical-scale universality,
without assuming its members have diffuse fourth moments in advance.
The cost is smaller than the previous bound derived from W4. The
new-spin entropy has disappeared from the scalar CLT error because
it was maximized exactly before replacing fields; it has NOT been
discarded from the physical parent or its Gaussian comparison.

Use the independently audited exchangeable-sign regularizer with
q0=floor(n^(2/3)(log n)^(2/3)). Its actual microscopic code at
T=n/sqrt(q0) has

```math
 Q(W)\le M_n+O(n^{4/3}(\log n)^{1/3}),\qquad
 H=O(n^{2/3}(\log n)^{2/3}).                          \tag{8}
```

Consequently (7) on that code costs O(n^(4/3)(log n)^(1/3)),
the SAME order as the regularization loss. This is a construction
step for selectable asymptotically minimizing full signings at all
orders, not a condition imposed on unknown exact minimizers.

The exact two remaining gaps are unchanged: the Gaussian parent value
is uncontrolled, and a parent optimum outside this shrinking window
has not been excluded. The restricted comparison error is much larger
than the protected window itself. No original cap bound or recurrence
is inferred from (8).

## 6. Classical mechanisms and what the combination adds

The scalar Stein factor, matrix Chernoff inequality, and subGaussian
maximum estimate are classical and are reconstructed in the linked
campaign artifacts. The contribution here is their ordered use:
choose one physical empirical frame; fix the old word; maximize the
new child exactly; replace its independent scalar fields; only THEN
pay for choosing the old word. The exchangeable regularizer supplies
the small old code on actual full signings without an optimizer
structure assumption. No external novelty claim is made without
further literature checking.
