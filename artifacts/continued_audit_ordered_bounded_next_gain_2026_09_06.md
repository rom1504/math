# Ordered bounded next-feedback gain after exact coherent regression

Date: 2026-09-06. Consequence of the audited first marked-history theorem,
boundary graph estimate, and exact Boolean stable-noise module. This is
an ordered averaged GAIN statement, not a cutoff-free state evolution,
an innovation floor, or an arbitrary-depth closure.

## 1. Actual bounded fields and homogeneous coefficients

Assume the first marked-history bounded smooth hypotheses, including
feasibility `0<=H<=1-|f|` and `|psi|<=1`. Retain literally

```math
W=(S,G,Y,QS,QD),\quad V=b_0QS+b_1QD,\quad Z=Br(G,Y),
\quad C=H(G,Y)\psi(V+Z).
```

Use that theorem's exact deterministic conditional mean c0 and first
coefficient a, and define

```math
L=B(c_0+D_aZ),\quad \eta=BC-L,
\quad u=f+C,\quad J=1-|u|,\quad
j_{2,n}=\frac1n\mathbb E\sum_i J_i|(Bu)_i|.
```

Thus `0<=J<=1`, `J` is a function of the OLD `(W,Z)`, and exactly
`Bu=V+Z+L+eta`. No Gaussian substitution is made in these definitions.

Let eta_p be the exact original Boolean homogeneous projection of eta
at odd degree p. Put, for each root,

```math
v_{i,p}=\mathbb E\eta_{i,p}^2,\qquad
c_{i,p}=\mathbb E[\eta_{i,p}L_i].
```

Different homogeneous degrees are exactly orthogonal. The fields are
globally odd and hence centered. For fixed P>=5 and delta>0 set

```math
I_i(P,\delta)=
 \sum_{\substack{5\le p\le P,\ p\ odd\\v_{i,p}\ge\delta}}
       \left(\sqrt{v_{i,p}}+
                          \frac{c_{i,p}}{\sqrt{v_{i,p}}}\right)^2.
                                                               (1)
```

Only the declared variance-truncated coefficients enter this expression.
Define the uniform degree-tail quantity

```math
\theta_P^2=\limsup_{n\to\infty}\frac1n
             \mathbb E\|\eta_{>P}\|_2^2.
```

The established ordered polynomial approximation implies
`theta_P -> 0` and `n^-1 E||eta_1+eta_3||^2 -> 0`.

## 2. Precise gain statement

For every fixed P>=5 and delta>0,

```math
\boxed{\quad
\liminf_{n\to\infty}
 \left[j_{2,n}-\sqrt{\frac2\pi}\frac1n
       \sum_i\mathbb E[J_i]\sqrt{I_i(P,\delta)}\right]
 \ge-\theta_P-\sqrt{P\delta}.
\quad}                                                     (2)
```

The conservative factor P can be replaced by the number of retained
odd degrees. An equivalent useful ordered consequence is

```math
\liminf_n j_{2,n}\ge\sqrt{\frac2\pi}
 \lim_{P\to\infty}\lim_{\delta\downarrow0}\liminf_n
       \frac1n\sum_i\mathbb E[J_i]\sqrt{I_i(P,\delta)}.    (3)
```

The quantities on the right increase when degrees are added or delta
is lowered and are uniformly bounded. The limits are taken in the
displayed order, with n FIRST. Equation (3) does not allow the variance
cutoff to depend on n or replace it by an untruncated coefficient.

## 3. Uniform L2 approximation and why degree tails are available

The bounded first-stage proof supplies polynomial approximants to the
actual C, c0, and Z with arbitrarily small limiting averaged L2 error.
Its ideal-channel coefficient control supplies bounded deterministic
approximants to a in the needed weighted L2 norm. Bounded-op transport
then gives the same approximation property for L and eta. At every
fixed stage all source coefficients come from a finite catalog, with
uniformly bounded row coefficients.

Every such approximant has a fixed finite original Boolean degree.
Orthogonal projection onto higher degrees is L2-contractive; hence the
approximation property implies uniform limsup degree tails for L, eta,
and Z. At each finite stage the main eta terms have degree at least five,
and the collision/variance-normal-ordering errors are o in averaged L2.
The same projection argument removes eta_1 and eta_3 in the bounded
limit. This is not an assumption that the bounded response itself has
finite degree.

There is also a uniform averaged second-moment bound for all the actual
fields here. C and c0 are bounded, a is bounded, Z has bounded averaged
variance, and all transports use bounded B. These bounds will justify
the root variance cutoff and the final absolute-value tests.

## 4. Fixed degree, variance floor, and root cap first

For the proof fix P,delta and a root cap R. Keep only roots where the
total variances of L, eta, and Z are at most R. On such a root retain
the degrees in (1), define `beta_{i,p}=c_{i,p}/v_{i,p}`, and keep the
ACTUAL regression residual

```math
R_i^{\mathrm{coh}}=L_i-
       \sum_{p\ retained}\beta_{i,p}\eta_{i,p}.          (4)
```

The coefficients are bounded by sqrt(R/delta). Regression is exact at
each retained degree, and

```math
\mathbb E[\eta_{i,p}R_i^{\mathrm{coh}}]=0,
\qquad \mathbb E|R_i^{\mathrm{coh}}|^2\le\mathbb E L_i^2.
```

The finite-polynomial result in
`continued_audit_boundary_graph_and_next_return_2026_09_06.md` controls
all full-noise contractions into higher degrees of L. Proper noise
contractions are already small, equal-degree covariance is removed by
(4), and old W has degree at most three. The audited zero covariance
of eta with the old Z components is retained as a separate input; it
is not inferred from unequal local Hermite labels.

To apply this to actual homogeneous projections, approximate in averaged
L2 first, keeping all target homogeneous degrees fixed. On the capped
roots discard also those where the local approximation error exceeds
a chosen small threshold. Their fraction is controlled by the averaged
error. On the remaining roots tensor-contraction differences are bounded
by local Frobenius errors times capped norms. Fixed-degree Boolean
hypercontractivity controls the associated finite-degree Stein moments.
Send n to infinity at the fixed approximation stage, and then remove
the approximation error and the additional discarded roots. This
transfers the mixed derivative conditions to the exact retained
projections; it does not assume a uniform operator bound for arbitrary
bounded-function approximants.

First truncate L and Z at a fixed original degree Q as well. The exact
Boolean conditional Stein/characteristic-function argument compares the
retained eta vector with independent centered Gaussian coordinates of
variances v_i,p, relative to the ACTUAL variables
`(W,Z_{<=Q},R_coh,<=Q)`. Send Q to infinity afterwards using the uniform
averaged L2 tails. This retains the distribution of `(W,Z,R_coh)` and
does not Gaussianize the regression residual or either marked return.

## 5. From the comparison to absolute gain

On the retained roots, omit the unretained eta terms from the next
field. The averaged L1 error is at most

```math
\theta_P+\sqrt{P\delta}+o(1).
```

Low degrees contribute o(1), the high-degree tail contributes theta_P,
and each of at most P omitted retained-range components has variance
below delta. Orthogonality gives the displayed bound.

The comparison in Section 4 turns the retained sum in the next field
into a centered Gaussian of variance I_i(P,delta), independent of the
literal `(W,Z,R_coh)`. Its coherent shift is still
`V_i+Z_i+R_i^coh`. For every real shift b and s>=0,

```math
\mathbb E_N|b+sN|\ge\sqrt{2/\pi}\,s.
```

Since J depends only on old `(W,Z)` and is nonnegative, this yields
the right-hand main term of (2) on the capped roots. Bounded continuous
approximations to J are justified by the first-stage L2 approximations
and the Lipschitz operation `J=1-|u|`; Gaussian-a.e. regularity of the
old-frame functions is the same hypothesis used there. Absolute-value
tests are obtained from bounded tests by the uniform averaged second
moments. No new discontinuous sign comparison is required for this gain.

Finally remove the root cap. Bessel's inequality for the orthogonal
eta_p gives

```math
\sum_p c_{i,p}^2/v_{i,p}\le\mathbb E L_i^2,
\qquad I_i(P,\delta)\le
     2\mathbb E\eta_i^2+2\mathbb E L_i^2.
```

Zero-variance terms are omitted in this inequality. The fraction of
discarded roots is O(1/R); Cauchy--Schwarz shows that their contribution
to the averaged square-root variance is O(R^-1/2), uniformly in P,delta.
The actual gain on those roots is nonnegative and may be discarded.
Taking R to infinity after the prior limits proves (2), and the ordered
degree/variance limits give (3).

## 6. Scope and the exact endpoint certificate

The two feasible Boolean-rounding mean endpoints
`+u+J sign(Bu)` and `-u+J sign(Bu)` give the usual exact finite certificate

```math
j_{2,n}+|e(u)+e(J\,\operatorname{sign}(Bu))|.
```

The present result controls its next cross gain through (2)--(3), while
retaining both endpoint self-energies. It does not prove that any term
in (1) is bounded below away from zero. A covariance c=-v can cancel a
retained innovation component exactly.

Very small eta variances can have very large regression coefficients.
They are left inside the actual coherent residual until the ordered
cutoff limits license otherwise. The result is therefore deliberately
not a cutoff-free Gaussian state law, nor an assertion that every later
feedback depth inherits the same primitive graph structure.

## 7. Two-second-moment scalar corollary

There is a simpler, potentially weaker consequence which has no degree
cutoff in its final observable. Define the ACTUAL row second moments

```math
\sigma_i^2=\mathbb E\eta_i^2,\qquad
\rho_i=\mathbb E[\eta_iL_i].
```

For zero variance set the expression below to zero; at every positive
cutoff it is already omitted. Then

```math
\boxed{\quad
\liminf_n j_{2,n}\ge\sqrt{\frac2\pi}
 \lim_{\varepsilon\downarrow0}\liminf_n
 \frac1n\sum_i\mathbb E[J_i]
 \left|\sigma_i+\frac{\rho_i}{\sigma_i}\right|
 \mathbf1_{\{\sigma_i^2\ge\varepsilon\}}.
\quad}                                                       (5)
```

Only the two actual moments of the returned nonlinear residual and its
coherent partner appear. In particular the mask remains
`J=1-|f+C|`, not `1-|C|`.

Here is a derivation directly from (2), which avoids assuming any new
cutoff-free Gaussian law. For a set A of retained homogeneous degrees,
write `v_A=sum_A v_p` and `c_A=sum_A c_p`. Cauchy--Schwarz gives

```math
\sqrt{\sum_{p\in A}(\sqrt{v_p}+c_p/\sqrt{v_p})^2}
  \ge \frac{|v_A+c_A|}{\sqrt{v_A}}.                         (6)
```

Work first on roots where `sigma^2>=epsilon` and
`sigma^2+E L^2<=K`. The omitted variance
`d=sigma^2-v_A` is nonnegative, and orthogonality gives
`|rho-c_A|<=sqrt(d E L^2)`. On `d<=epsilon/2`, the scalar function in
(6) therefore differs from `|sigma^2+rho|/sigma` by at most
`C(epsilon,K)(sqrt(d)+d)`. The other roots have proportion at most
`2 average(d)/epsilon`; on the capped set their target contribution is
bounded by `2sqrt(K)` times that proportion. Averaged d is controlled
by the high-degree tail, the vanishing degrees one and three, and
`P delta`. Thus first n, then the retained-degree variance cutoff
delta, and then P tending to infinity make this comparison error zero.
The errors on the right of (2) vanish in the same order.

Finally remove the cap K. Pointwise

```math
\left|\sigma+\rho/\sigma\right|^2
 \le 2\sigma^2+2\mathbb E L^2,
```

by Cauchy--Schwarz for the covariance. Uniform averaged second moments
make the cap's root fraction O(1/K), and another Cauchy--Schwarz bound
makes its contribution O(K^-1/2). This proves (5) for each fixed
epsilon, after which its displayed monotone cutoff limit is allowed.

At a fixed finite polynomial stage there is also an equivalent direct
Stein explanation. The inverse-generator derivative of
`eta=sum_p eta_p` is `sum_p D eta_p/p`; its mixed carré-du-champ with L
has expectation rho, while the noise self term has expectation sigma^2.
The audited proper and unequal-degree full-contraction estimates make
their nonconstant parts negligible. Thus the actual scalar regression
`R=L-rho eta/sigma^2` has vanishing mixed Stein term. On fixed
`sigma^2>=epsilon` and a row cap, eta is asymptotically Gaussian relative
to the retained actual `(R,W,Z)`. This explanation is consistent with
(5), but (6) and the controlled tail argument already prove the stated
bounded gain corollary.

No positive innovation floor follows. The cancellation `rho=-sigma^2`
can make this lower contribution zero. The small-variance cutoff is
taken after n and cannot be silently replaced by an n-dependent cutoff
or by an unconditional state-evolution assertion.

## 8. Independent audit of the stronger fixed-cutoff scalar formula

The full companion
`continued_feedback_second_query_scalar_gain_2026_09_06.md` was read
after the preceding corollary was banked. Its stronger fixed-cutoff
comparison (4) and low-variance completion (5) pass independently.

Use the actual sigma, rho and L above, and write T=V+Z. On
`G={sigma^2>=epsilon, E L^2+sigma^2<=K}`, set
`R=L-rho eta/sigma^2` and `t=sigma+rho/sigma`. For each fixed epsilon,K,
the difference between the two averaged quantities

```math
\frac1n\sum_{i\in G}\mathbb E[J_i|(Bu)_i|],\qquad
\frac1n\sum_{i\in G}\mathbb E[J_i\mathbb E_N|T_i+R_i+t_iN|]
```

tends to zero. This keeps the actual joint law of R and the old fields;
it is not an independent replacement of L. The Gaussian is independent
of that retained joint law only within the stated averaged comparison.

At a fixed polynomial stage the inverse-generator derivative gives the
correct covariance constants even when eta contains multiple original
degrees. The new full-contraction criterion handles smaller noise degrees
against larger coherent degrees; proper noise cuts handle the opposite
inequality; scalar regression subtracts the full equal-degree constant.
The corresponding old-Z covariance is separately negligible in averaged
absolute value and is not discarded by a false degree-alias argument.

For the bounded passage, averaged L2 approximation gives averaged L1
error of sigma^2 and rho. Delete roots where these statistic errors are
large. On the remaining roots the approximating variance stays above
epsilon/2 and all regression coefficients lie in a compact set.
Coefficient error times local second moments is consequently controlled,
as is the L2 difference of the actual and approximating R. For the
mask use `(1-|u_polynomial|)_+`, a bounded Lipschitz approximation to J.
The finite-degree stable comparison first applies to bounded tests;
uniform averaged second moments then give the required linear-growth
tests. T need not have a pointwise variance cap: its bounded AVERAGED
second moment suffices both here and on deleted rows.

On the complementary low-noise rows with the same root cap, omitting
eta changes the averaged gain by at most sqrt(epsilon), leaving the
actual `J|T+L|`. The high-cap rows have fraction O(1/K), while bounded
u and bounded-op B give averaged second moment of Bu at most ||B||^2.
Their gain is therefore O(K^-1/2). This reconstructs the source's
fixed-cutoff completion formula with errors
`O(sqrt(epsilon)+K^-1/2)+o_n(1)`.

The finite endpoint normalization is
`max_x |x^T B x|/(2n)=Qabs(A)/(n sqrt(n-1))`.
For endpoints `u +/- J sign(Bu)`, half the difference of these normalized
energies is exactly j_(2,n). Thus the completed formula is an actual
second-query certificate for the original signing. It still supplies
no lower bound on sigma, no protection against rho=-sigma^2, and no
automatic extension beyond the audited first marked history.
