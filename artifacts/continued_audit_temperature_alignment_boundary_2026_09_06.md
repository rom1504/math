# Independent audit: information slack and separated Gaussian boundaries

Date: 2026-09-06. Independent reconstruction by `continued_audit`.

Sources read in full:

- `continued_convergence_temperature_alignment_next_target_2026_09_06.md`;
- `continued_feedback_two_gaussian_correlation_boundary_2026_09_06.md`,
  including its appended singular-correlation Section 6.

**Verdict: PASS at the stated scopes.** The first file isolates a sufficient
optimal-channel information inequality, not a proof of it. The second
rigorously excludes two particular separated-Gaussian candidate
counterexamples. Neither proves `B E_t <= E_t`, evaluates the limiting
Bellman fixed point, or changes the original minimax bound.

## 1. Exact information slack and its scope

For an admissible pair with equal signed marginals, choose local channels
`M|U`, `N|V` independently given `(U,V)`, and put `L=(M,N)`. Conditional
independence gives

`I(U,V;M,N)=I(U;M)+I(V;N)-I(M;N)`.

The interaction-information identity, followed by the orthogonal change
of variables, gives

`I(A;L)+I(B;L)`
`= I(U;M)+I(V;N)+I(A;B)-I(M;N)-I(A;B|L)`.

The fair selector between A and B is independent of its selected source
value because their marginals agree. Therefore dividing by two is
exactly the information of the proposed parent label, including the
selector. There is no lost `I(M;N)` term. The variance comparison follows
from conditioning and the pointwise orthogonal identity. As `g_t` is
decreasing, replacing the actual refined residual variance by the average
of the two child distortions makes the sufficient inequality stronger,
not weaker.

The sparse-ternary diagonal-pair example correctly refutes the proposed
inequality for arbitrary channels: its information slack tends to zero,
whereas strict convexity leaves a fixed positive variance Jensen gap.
This example does not use optimal child channels.

For a finite source, the reproduction-distribution variational formula
and directional derivative imply `Z(y)<=1` everywhere and equality at
active reproduction points. Projection to the source convex hull
justifies the assertion also outside that hull. At an active point y,
the ratio of kernels at y+h and y is exactly
`exp(2 lambda h (X-y)-lambda h^2)`. Thus the posterior moment bound and
posterior mean identity follow without any Gaussian replacement. The
variance proxy is `1/(2 lambda)`; it is not the posterior's actual
variance. The warning against upgrading it is essential.

For independent Gaussian residuals of fixed variances a,b, rotation
gives conditional mutual information
`(1/2) log((a+b)^2/(4ab))`. Direct differentiation of the Gaussian
potential gives `v^2 g_t''(v)=rho^2/[2(1+rho^2)]<=1/4`. Concavity of
`g_t(v)+(1/4)log v` therefore pays the requisite Jensen gap with the
stated factor-of-two margin. This calculation does not generalize to
arbitrary subGaussian residual laws merely from their variances.

## 2. Uniform rate-distortion separation and temperature localization

Let `nu_V=w N(0,v)+q N(0,V)`, where `w=1-q`, `0<q<1/2`, v,t are fixed.
For temperatures in any compact subset of `(0,t]`, revealing the
component and then using conditional Gaussian channels gives an
`O(log V)` upper bound, uniformly in temperature. This is a legitimate
channel of X: one may first sample its component posterior given X.

For the lower bound, extend an arbitrary X-channel by the canonical
Markov law `J-X-L`. Its conditional-mean reproduction has squared error
`O(log V)` because lambda has a positive lower bound. For any fixed
`0<beta<1/2`, classifying by whether the reproduction exceeds `V^beta`
in magnitude makes error at most the sum of a low-Gaussian tail,
`O(V^(beta-1/2))`, and `O(log V/V^(2 beta))`. All tend to zero uniformly.
Binary Fano gives `H(J|L)=o(1)`. The exact Markov information identity is

`I(X;L)=I(X;L|J)+h(q)-H(J|L)`.

Adding J can only decrease reproduction error. Conditional Gaussian
rate-distortion lower bounds yield the asserted uniform expansion of
`J_lambda(nu_V)`. Finite reproduction quantization preserves the
finite-latent convention to arbitrary accuracy.

The optimization really is localized away from lambda=0. If
`lambda V>=1/2`, source concavity and the high-component Gaussian lower
bound make the normalized objective at most
`C+(1-2q)log(lambda)/4`. If `lambda V<1/2`, the elementary nonnegative
cost bound gives `C-(1-2q)log V/4`. Both tend to minus infinity along
any vanishing-temperature sequence. A fixed positive temperature has a
finite normalized limit, so neither regime can maximize. The hypothesis
`q<1/2` is used precisely here.

## 3. Scalar optimizer, curvature, and fixed correlations

With `alpha=1-2q`, the localized objective has derivative

`alpha/(4 lambda)-1/[4(2t-lambda)]-wv`.

It is strictly decreasing, positive at zero and negative at t. If the
low Gaussian were also in its logarithmic rate-distortion branch, the
derivative would instead be
`-1/(4 lambda)-1/[4(2t-lambda)]<0`. Thus the maximizer is indeed in the
linear low-component branch. Multiplication of the derivative equation
gives exactly

`2wv lambda^2-w(1+4tv)lambda+t alpha=0`.

The smaller root is the unique optimizer. Differentiation gives
`A_v=-w lambda`. The identity

`v lambda=(alpha-wz)/[2w(2-z)]`, where `z=lambda/t`,

has strictly negative derivative with respect to z, while z decreases
with v. Hence `v lambda` increases and is bounded by `alpha/(4w)`.
It follows that `0<=A_vv<=alpha/(4v^2)`. This justifies the concavity of
`A(v)+(alpha/4)log v` and the exact logarithmic midpoint bound.

Component entropy decomposition for the Gaussian pair has errors
controlled by its one-coordinate component-classification entropy.
Therefore the information expansion remains valid for every fixed
nonsingular pair of component correlations. Substitution of the two
child scale shifts cancels the high-component correlation exactly.
The remaining coefficient is
`w/4-alpha/8=1/8`, giving

`limsup gap <= -h(q)/2+(1/8)log(1-r_low^2)<0`.

## 4. Singular matching-variance policy

For `r_low=1-c/v`, `r_high=1-c/V`, with fixed `0<c<2v`, the minus child
is exactly `N(0,c)`. It is independent of the plus child: conditional
independence and the identical minus law in both component classes
factor the joint mixture. The plus child has variances `2v-c,2V-c`.

The pair-information expansion is still justified by the recoverability
of the component from one input coordinate, independent of correlation.
Its growing `q log V/2` term cancels the remaining source/child scale
terms and gives exactly equation (8) of the source file.

Set `a=(1-2q)/4`, `B(v)=A(v)+a log v`, and
`G(v)=g_t(v)+(1/4)log v`. At the Gaussian optimizer, write
`y=lambda_G v=rho/[2(1+rho)]<1/4`. Evaluation at this allowed temperature
gives

`B(v)-G(v) >= q[y-1/2-(1/2)log(2y)]`
`              >= q[(log 2)/2-1/4]`.

The scalar expression decreases on `(0,1/2)`, so the second inequality
has the correct direction. Expanding all logarithms in the policy gap
gives coefficients `1/8,1/8,-1/4` on `log(2v-c),log c,log v`, respectively.
Thus the exact rewritten gap is

`B(2v-c)/2+G(c)/2-B(v)+(1/8)log[c(2v-c)/v^2]`.

Replace G(c) by its strict B(c) upper bound and use concavity of B. The
result is

`gap <= -q(2 log 2-1)/8+(1/8)log[c(2v-c)/v^2]`
`    <= -q(2 log 2-1)/8 < 0`.

This is an analytic exclusion for the specified singular policy, not
evidence that every correlation depending on V, every Gaussian mixture,
or every reachable Bellman state obeys temperature alignment.

## 5. New conditional-variance envelope supersolution: finite-source proof

The following is a separate successful supersolution, not the unresolved
original envelope inequality. Define

`Ebar(nu)=sup_L [E g_t(Var(X|L))-I(X;L)]`.

For finite source alphabets and finite child labels, the complete proof
below has no measure-theoretic extension requirement. These are all the
states and labels needed in a fixed finite Bellman tree started from the
ternary source.

Choose the child channels and L as in Section 1. Use the parent label
`(L,B)` for A, and L for B, then use the independent fair selector. All
these labels are finite. The exact information sum is

`I(A;L,B)+I(B;L)=I(A;B)+I(A,B;L)`
` <= I(A;B)+I(U;M)+I(V;N)`.

For each fixed L, write the conditional covariance of (A,B) as Sigma,
with diagonal entries a,b and off-diagonal c. If b>0, conditional-mean
prediction does at least as well as linear prediction, so

`E[Var(A|L,B)|L] <= a-c^2/b = det(Sigma)/b`.

Convexity and monotone decrease of g give

`E[g(Var(A|L,B))|L]+g(b)`
` >= g(det(Sigma)/b)+g(b)`.

Let the eigenvalues be `lambda_min,lambda_max`. Both b and det/b lie
between the eigenvalues and have the same product. Their logarithms
therefore form a less spread pair with the same sum. Directly,

`d/d log v [g_t(v)] = -rho/[2(1+rho)]`,

which decreases with v; hence `g_t(exp s)` is concave in s. It follows
that the Schur-pair sum is at least the eigenvalue-pair sum. Convexity of
g, together with diagonal majorization by eigenvalues, then bounds this
from below by

`g(Var(U|L))+g(Var(V|L))`.

Additional conditioning increases expected g of conditional variance:
apply total variance, monotone decrease, then convexity. Consequently
the last expression, after averaging L, is at least

`E g(Var(U|M))+E g(Var(V|N))`.

Combining reward and information inequalities proves
`B Ebar <= Ebar` on finite sources. All inequality directions are as
displayed. If b=0, c=0 and the corresponding assertion is immediate.
If det=0, approximate by positive definite covariance or take the
continuous limit `g(0)=0`; no covariance inverse is needed in the final
statement.

The director's comparison `E<=Ebar<=Phi` also checks. For the upper
comparison, mix optimal conditional self-couplings: the chain bound is
`I(X;Y)<=I(X;L)+I(X;Y|L)`, so
`F(nu)<=I(X;L)+E F(nu_L)`. Together with Gaussian extremality this is
stronger than the required `Ebar<=Phi` comparison. Thus the already
proved Gaussian-boundary iteration lies below Ebar at every finite
depth. A certified negative root offset for Ebar would be a genuine
upper-construction certificate; the earlier certificate for E alone
cannot be substituted for it.

## 6. A finite-mesh upper-error bound for the ternary Ebar evaluation

For posterior parameters `(z,s)` in `[0,1]^2`, let

`q(z,s)=h(z)+z h((1+s)/2)+g_t((z-z^2 s^2)/p)`.

On a mesh of width delta, round z randomly to its adjacent mesh points,
preserving its mean, and round s to its nearest mesh point. The posterior
ternary laws change in total variation by at most `5 delta/4`. Entropy
continuity therefore costs at most

`h(5 delta/4)+(5 delta/4)log 2`.

The variance derivatives satisfy the sharper uniform bounds
`|partial_z v|<=1/p`, `|partial_s v|<=2/p`. Since `|g_t'|<=t`, its
change costs at most `2t delta/p`. Thus the full continuous concave-hull
value is at most its discrete mesh concave hull plus

`h(5 delta/4)+(5 delta/4)log 2+2t delta/p`.

The randomized z rounding preserves the required barycentre exactly.
This bound can turn an outward-certified mesh evaluation into a global
upper bound; an unverified floating grid remains only diagnostic.
