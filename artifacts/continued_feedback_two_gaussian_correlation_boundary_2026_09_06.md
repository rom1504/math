# Two separated Gaussian components with label-dependent correlations

Date: 2026-09-06. A rigorous non-falsification result for one asymptotic
temperature-segregation family. It does not prove general `B E<=E`.
The complete proof, including the singular case in Section 6, passed
independent audit-agent reconstruction.

Fix t>0, v>0, and 0<q<1/2. Put w=1-q and alpha=1-2q. Consider the
symmetric source

`nu_V=w N(0,v)+q N(0,V)`, with V tending to infinity.

Given the component label, pair two Gaussian coordinates of the same
component with fixed correlations r_low,r_high in (-1,1), respectively.
This is a valid Bellman pair with both input marginals nu_V. The two
children remain two-Gaussian mixtures, and their component variances are
`v(1+/-r_low)` and `V(1+/-r_high)`.

The result is that the selected Bellman value, with the FULL latent
envelope optimized separately in both children, satisfies

`limsup_(V->infinity) [selected B E(nu_V)-E(nu_V)]`
` <= -h(q)/2+(1/8)log(1-r_low^2) < 0`.                      (1)

In particular fixed label-dependent correlations do not falsify
temperature alignment in this widely separated-variance family. The
high-component correlation cancels from the exact limiting gap.

## 1. Uniform fixed-temperature rate-distortion asymptotic

Let J denote the binary component label. For lambda in a fixed compact
subset of (0,t],

`J_lambda(nu_V)=h(q)+w J_lambda(N(0,v))`
`                    +q J_lambda(N(0,V))+o(1)`.             (2)

The upper bound samples the component label from its posterior given X
and gives this label to the reproduction channel, then uses the Gaussian
channels conditionally. Thus it is a valid channel from X, not an
assumption that the channel is given an inaccessible external label.
For the lower bound, extend any nearly optimal X-channel L by the
canonical component label J, conditionally independent of L given X.
Then J--X--L is a Markov chain and the conditional source laws remain
the two specified Gaussians. Replace the reproduction by E[X|L].
The upper construction bounds its objective by O(log V), uniformly over
the compact temperature set, so its squared error is O(log V).

Choose any beta in (0,1/2). Classifying the component according to whether
the reproduction magnitude exceeds V^beta has error tending to zero:
the low Gaussian exceeds this threshold with vanishing probability,
the high Gaussian falls below a fixed multiple of it with probability
O(V^(beta-1/2)), and reproduction errors of size V^beta cost at most
O(log V/V^(2beta)). Thus H(J|L)=o(1). The Markov information identity
`I(X;L)=I(J;L)+I(X;L|J)` yields

`I(X;L)>=h(q)-o(1)+I(X;L|J)`.

Conditioning the prediction on J only decreases its squared error.
The conditional Gaussian lower bounds prove (2), uniformly on the
declared compact temperature set. Finite quantizations of the Gaussian
reproduction channels justify the finite-latent-alphabet convention.

## 2. The source envelope has an explicit scalar limit

Write `c_t(lambda)=log[lambda(2t-lambda)/t^2]/4`. The Gaussian Lagrangian
equals lambda v below lambda v=1/2 and `(1+log(2lambda v))/2` above it.
Equation (2) gives

`E_t(nu_V)=-q log(V)/2-h(q)+A(q,v,t)+o(1)`,                 (3)

where

`A(q,v,t)=sup_(0<lambda<=t)`
` {c_t(lambda)-w lambda v-q[1+log(2lambda)]/2}`.              (4)

Two details justify this formula. First, temperatures approaching zero
cannot maximize the normalized expression in (3). Concavity in the
source gives `J_lambda(nu_V)>=q J_lambda(N(0,V))`. If lambda V>=1/2,
the normalized envelope is at most a constant plus
`(1-2q)log(lambda)/4`, tending to minus infinity as lambda tends to zero.
If lambda V<1/2, the bound `J>=0` instead gives a constant plus
`(q/2-1/4)log V`, also tending to minus infinity. This localizes the
optimization in a compact temperature set where (2) is uniform.

Second, the maximizing temperature has lambda v<1/2. Above that point,
using the logarithmic formula for BOTH Gaussian components gives
derivative `c_t'(lambda)-1/(2lambda)<0`. At the matching boundary the
derivative from below is already negative. Thus the lower-variance
component is optimally left without a reproduction channel, and its
cost in (4) is indeed lambda v.

The objective in (4) has derivative

`alpha/(4lambda)-1/[4(2t-lambda)]-wv`,

which is strictly decreasing from infinity to a negative value. The
unique optimizer solves

`2wv lambda^2-w(1+4tv)lambda+t alpha=0`,                    (5)

using the smaller root. Every scalar quantity in the limit is explicit;
there is no growing latent-state optimization in this family.

## 3. A curvature inequality for A

Let lambda(v) denote the optimizer in (5) and put a=alpha/4. The envelope
derivative is `A_v=-w lambda(v)`. Its stationarity relation gives, with
z=lambda/t,

`v lambda=[alpha-wz]/[2w(2-z)]`.

The right side is strictly decreasing in z, while lambda decreases as v
increases. Hence v lambda increases with v and is bounded above by
alpha/(4w). Consequently

`0<=A_vv=-w lambda'(v)<=alpha/(4v^2)=a/v^2`.

Thus `A(q,v,t)+a log v` is concave as a function of v. Applying its
midpoint inequality gives

`[A(v(1+r))+A(v(1-r))]/2-A(v)`
` <= -(a/2)log(1-r^2)`.                                    (6)

The omitted q,t arguments are fixed throughout this display.

## 4. Exact pair information and cancellation of the high correlation

For fixed nonsingular correlations, a single input coordinate identifies
the component asymptotically. Therefore its component information tends
to h(q), and so does the information from the pair. Gaussian conditional
entropy gives

`I(A;B)=h(q)-(w/2)log(1-r_low^2)`
`              -(q/2)log(1-r_high^2)+o(1)`.                 (7)

The two child expansions in (3) contain high-variance scale shifts
`-q log(1+/-r_high)/2`. Averaging them and subtracting half of (7)
cancels their entire r_high dependence. The exact limiting difference
from the source envelope is

`-h(q)/2+[A(v(1+r_low))+A(v(1-r_low))]/2-A(v)`
`                      +(w/4)log(1-r_low^2)`.

Insert (6). Since `w/4-alpha/8=1/8`, this proves (1).

## 5. Scope and the remaining singular boundary

The proof allows different temperatures in the two child envelopes and
different fixed correlations in the two component classes. It specifically
rules out this proposed asymptotic obstruction; it is not a numerical
non-falsification statement.

It does not cover correlations approaching +/-1 as V grows. In particular,
`r_high=1-c/V` leaves a bounded-variance high component in the minus child,
so its two-scale expansion changes. A further special choice
`r_low=1-c/v`, with 0<c<2v, makes that entire minus child exactly N(0,c).
The next section resolves this particular singular matching-variance
regime as well. General component-dependent singular correlations remain
outside the scope; no general temperature alignment is inferred.

## 6. Singular correlation that makes one child exactly Gaussian

Fix 0<c<2v and now set

`r_low=1-c/v`, `r_high=1-c/V`.

The minus child is exactly N(0,c), independent of the plus child: the
conditional sum and difference are independent Gaussian coordinates,
and the difference has the SAME law in both component classes. The plus
child has component variances 2v-c and 2V-c. Entropy decomposition, with
the asymptotically recoverable label, gives

`I(A;B)=h(q)-(w/2)log[c(2v-c)/v^2]`
`                   +(q/2)log V-(q/2)log(2c)+o(1)`.

Using the source and plus-child expansions from Section 2, the exact
limiting policy-minus-source gap is

`Delta=A(q,2v-c,t)/2+g_t(c)/2-A(q,v,t)`
`              +(w/4)log[c(2v-c)/v^2]+(q/4)log c`.           (8)

This quantity is strictly negative for ALL the parameters in the stated
scope, not merely at large t. To see it, define

`B(v)=A(q,v,t)+a log v`, `G(v)=g_t(v)+(1/4)log v`,

where a=(1-2q)/4. Section 3 proves concavity of B. Moreover,

`B(v)-G(v)>=q[(log2)/2-1/4]>0`.                            (9)

Indeed evaluate the supremum defining A at the Gaussian optimizer
lambda_G for g_t(v). Write y=lambda_G v. The difference obtained is

`q[y-1/2-(1/2)log(2y)]`.

The Gaussian optimizer obeys
`y=r/[2(1+r)]<1/4`, where r is its positive correlation parameter.
The displayed scalar expression decreases on 0<y<1/2, so evaluating
it at y=1/4 proves (9).

Rewriting (8) in these variables gives exactly

`Delta=B(2v-c)/2+G(c)/2-B(v)`
`                           +(1/8)log[c(2v-c)/v^2]`.

Use (9), then the midpoint concavity of B. The resulting quantitative
upper bound is

`Delta<=-q(2log2-1)/8+(1/8)log[c(2v-c)/v^2]`
`      <=-q(2log2-1)/8<0`.                                 (10)

Thus this natural singular correlation policy, which deliberately makes
one output Gaussian and lets the other retain the large-variance class,
also cannot violate `B E<=E` in the separated-variance limit. Equations
(1) and (10) are scoped analytic exclusions of candidate counterexample
families, not a proof covering all pairings or all Gaussian mixtures.
