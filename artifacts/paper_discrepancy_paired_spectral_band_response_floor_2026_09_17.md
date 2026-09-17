# Sharp response floor for the paired Gaussian-sign spectral band

2026-09-17. This is a limitation of a specified physical-law class,
not a lower bound for arbitrary physical sign laws or all deterministic
parent constructions. It explains why improving constants in the
[low-cap response theorem](paper_director_low_cap_uniform_response_2026_09_17.md)
cannot by itself reach the original near-order slope.

## 1. Exact covariance bound and uniform asymptotic response floor

Fix 0<a<1. For any symmetric hollow matrix T with ||T||_op<=a, let
nu_T be the equal mixture of sign N(0,I+T) and sign N(0,I-T). Let
L_(n,a) be the convex hull of all these physical sign laws. Its members
are centered, exactly isotropic, and have full physical support.

Put kappa=sqrt(2/pi) and

```
s_a=(2/pi)asin(a),
f(s)=(sqrt(1+s)+sqrt(1-s))/2,
beta(a)=kappa f(s_a).
```

Uniformly over EVERY Boolean query x and EVERY nu in L_(n,a),

```
E_nu |h.x|/sqrt(n) >= beta(a)-epsilon_n(a),
epsilon_n(a)=O_a(n^(-1/6)sqrt(log(en))) ->0.          (1)
```

The dependence on a is fixed before n tends to infinity. No claim of
uniformity as a approaches one with n is made.

Here is the proof. For each positive integer j, the entrywise power
T^(circ j) is a compression of T^(tensor j): use the isometry sending
e_i to e_i^(tensor j). Therefore

```
||T^(circ j)||_op<=||T||_op^j.
```

All coefficients in the odd arcsine series are positive. Operator-norm
convergence of that series for ||T||_op<=a<1 gives

```
||asin(T)||_op<=asin(a).                              (2)
```

This is entrywise arcsine, not spectral functional calculus. The exact
sign covariance matrices are I+V and I-V, where V=kappa^2 asin(T).
Consequently the paired query variance contrast obeys

```
v_x=x^T Vx/n,              |v_x|<=s_a.               (3)
```

The uniform bounded-spectrum scalar Gaussian-sign theorem, reconstructed
in [the quenched proof](flatify_independent_2026_09_07_gaussian_sign_quenched_universality.md),
applies with the fixed latent interval [1-a,1+a]. It gives the actual
physical response kappa f(v_x) with error at most epsilon_n(a).
Since f is even and decreases with |v|, (3) proves (1) for every pair.
The same pointwise lower bound is preserved by arbitrary mixtures.

In particular the exact spectral band used by the current positive
theorem has a=1/2 and s_a=1/3. Its floor is

```
beta(1/2)=(sqrt(2)+1)/sqrt(3pi)
          =.7863938738970605... .                    (4)
```

This strengthens the weaker floor kappa f(1/2) that follows only from
the Gaussian Holder subGaussian proxy. Neither value is asserted as an
exact finite-n lower bound without the error term in (1).

## 2. Sharpness: independent two-coordinate blocks

The constant in (1) is sharp for every fixed a. Take n=2m and let P
be the permutation matrix of a perfect matching, with zero diagonal.
Set T=aP and query x=(1,...,1). The matrices I plus or minus T split
into independent two-coordinate blocks. Gaussian signs in one block
have correlation plus or minus s_a, respectively. Therefore under the
two branches,

```
Var(h.x)/n=1+s_a,       respectively 1-s_a.
```

The bounded independent-block central limit theorem, followed by
uniform integrability from the bounded second moments, gives

```
E_(nu_T) |h.x|/sqrt(n) -> beta(a).                    (5)
```

Thus there is no uncharged slack in the Schur estimate relevant to
this class-wide floor. The example is a law/query example, not an
assertion that this query is an optimizer of a low-cap signing.

At a=1/2, the half-sum of a sign pair is a lazy walk step in{-1,0,1}.
In the positive branch its probabilities are (1,1,1)/3, and in the
negative branch (1,4,1)/6. Exact integer convolution verifies every
finite response and variance in the accompanying replay. For example,
the normalized paired responses at n=128,512,1024 are approximately
.784669408, .785963174, .786178558, approaching (4). These decimals
are diagnostics; the proof of convergence is the elementary block CLT.

## 3. Why the current certificate cannot meet the parent slope

The near-order scalar target at normalized old cap c is response below
3c/2. Throughout the reported regime c<=.493608094 this is at most
.740412141, strictly below (4). Even the weaker target at c=1/2 is
3/4, and

```
beta(1/2)>3/4.                                      (6)
```

An entirely rational check of (6) uses sqrt(2)>7/5 and pi<22/7:
beta(1/2)^2=(3+2sqrt(2))/(3pi)>203/330>9/16.

Consequently NO improved covariance localization, cap-budget estimate,
choice of dual law, or optimization of the same minimax proof can
produce a mean-response certificate below 3c/2 while its output stays
inside L_(n,1/2). This is stronger than observing that the displayed
constant 2^(-67) is small. It is a sharp class obstruction independent
of the energy geometry.

It does not rule out non-Gaussian physical laws, unequal or nonpaired
Gaussian mixtures, a larger fixed spectral band, structured dependent
columns, rare selected realizations, or a different parent construction.

## 4. A typical independent-column extension also has the wrong first slope

There is a direct actual-sign consequence at the usual independent-
column deployment. Fix A_n, a signed ground word x_0, and q=floor(epsilon n).
Independently sample physical bridge columns h_j, each from L_(n,a).
Their laws may differ. Let the new child D be any fixed full signing,
or one chosen independently of these columns. Complete zero field ties
by independent fair signs and choose

```
y_j=sign(H_A(x_0))*sign(h_j.x_0).
```

Every column law is globally sign symmetric, so the y_j are independent
fair signs. Hence the child energy has mean zero and variance binom(q,2),
regardless of its cap. Also exact isotropy gives
Var(|h_j.x_0|)<=n. Thus for each fixed epsilon>0, as n tends to infinity,

```
sum_j |h_j.x_0| >= q sqrt(n)[beta(a)-o(1)]
                         -o_probability(n^(3/2)),
H_D(y)=o_probability(n^(3/2)).
```

Evaluating the full parent at this ONE actual Boolean word gives

```
Q(parent)/n^(3/2)
 >= Q(A_n)/n^(3/2)+epsilon beta(a)-o_probability(1). (7)
```

On a subsequence with Q(A_n)/n^(3/2)->c, the lower normalized parent
curve is (c+epsilon beta(a))/(1+epsilon)^(3/2). Its derivative at
epsilon=0 is beta(a)-3c/2. Therefore typical independent deployment
of the current half-band laws cannot have the desired improving local
slope in the reported cap regime, even after adding an arbitrary
independent full-sign child. The order of limits is n first at fixed
epsilon, then epsilon down to zero.

This is not an impossibility theorem for exceptional deterministic
bridge choices from the full support, or for a child chosen adaptively
after the bridge. Those cases are not covered by the independent fair
new-spin argument.

## 5. Necessary enlargement of the latent band

For a prospective target r=3c/(2kappa) between 1/sqrt(2) and one,
the necessary condition beta(a)<3c/2 is exactly

```
a > sin(pi*r*sqrt(1-r^2)).                           (8)
```

For orientation only, the right side is approximately

| c | Required band radius a, strictly larger than |
|---|---:|
| .4333221116640807 | .996261746 |
| .45 | .988214969 |
| .493608094 | .884960093 |
| .5 | .845565195 |

These are necessary conditions for this particular paired class,
NOT sufficient energy-response constructions. Near the lower endpoint,
even its sharp class floor requires a latent minimum eigenvalue below
about .00374. Any theorem using such a band must keep that positive
floor fixed before its scalar-comparison limit. Letting it collapse
with n is a separate unproved uniformity problem here.

## 6. Replay

```
.venv/bin/python computations/paper_discrepancy_2026_09_17_paired_band_floor.py
```

The script checks exact lazy-walk distributions, symmetry and variances
through n=1024, and the rational strict inequality in (6). It separately
labels the general-band threshold decimals as diagnostic. The sharp
operator bound, asymptotic comparison, and deployment argument are
proved above rather than inferred from those finite computations.
