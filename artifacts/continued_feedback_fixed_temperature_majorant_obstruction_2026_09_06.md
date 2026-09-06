# A uniform additive fixed-temperature terminal majorant cannot close the live point

Date: 2026-09-06. Exact obstruction to one specific upper-certificate
strategy at p=31/32,t=4. It does not falsify the recursive Bellman candidate,
the latent lower theorem, or a state-dependent/adaptive upper construction.
The full proof and exact checker passed an independent audit and replay;
the replay returned the same rational lower endpoint as Section 4.

## 1. The precise strategy under test

Use the notation of
`continued_convergence_deep_latent_envelope_2026_09_06.md`, fully read
including its fixed-temperature supersolution proof. For 0<lambda<=t,

`J_lambda(nu)=inf_L [I(X;L)+lambda E Var(X|L)]`,
`c_t(lambda)=(1/4)log[lambda(2t-lambda)/t^2]`,
`e_lambda(nu)=c_t(lambda)-J_lambda(nu)`.

The verified inequality is `B e_lambda<=e_lambda`. Adding a scalar
constant d preserves it. If a uniform additive correction d makes
`Phi_t<=e_lambda+d` on the relevant terminal states, monotonicity and
this supersolution inequality give the propagated root certificate

`B^r Phi_t(nu_p) <= e_lambda(nu_p)+d`.                       (1)

The target needs the right side, plus `p log2+t(1-sqrt(p))`, to be
negative. We show that the zero state and the Gaussian limit ALONE
force this propagated exponent to be at least

`19688465594483/320000000000000 > .06152645498`               (2)

when p=31/32,t=4, for EVERY lambda in (0,4]. This is an obstruction to
closing the calculation by the fixed-temperature supersolution inequality
with a uniform additive terminal error. Computing further Bellman iterates
of that majorant explicitly could in principle improve on (1); that is
not excluded by this result.

## 2. The two unavoidable terminal tests

At the zero law, Phi_t=J_lambda=0. Hence terminal domination requires

`d+c_t(lambda)>=0`.                                        (3)

At a unit-variance Gaussian, put `g=Phi_t(N(0,1))`. The rate-distortion
Lagrangian is exactly

`J_lambda(N(0,1))=lambda` for lambda<=1/2,
`J_lambda(N(0,1))=(1+log(2lambda))/2` for lambda>=1/2.         (4)

For completeness, Gaussian maximum entropy conditional on the latent
variable gives `I(X;L)>=max(0,(1/2)log(1/D))`, where
`D=E Var(X|L)<=1`. Minimizing this bound plus lambda D gives (4).
Gaussian reproduction channels attain every D in (0,1], and finite
quantizations approximate their information and squared prediction error.
Thus restricting the defining infimum to finite latent alphabets does
not change this value.

Terminal domination at the Gaussian limit also requires

`d+c_t(lambda)>=g+J_lambda(N(0,1))`.                        (5)

Consequently the proposed root upper bound in (1) cannot be smaller than

`max(0,g+J_lambda(N(0,1)))-J_lambda(nu_p)`.                  (6)

The auxiliary constant c_t CANCELS from this obstruction. In particular,
making its temperature-dependent offset very negative does not help:
the required terminal correction increases by exactly the same amount.

## 3. One explicit root channel bounds every temperature

Let s=9/10. Take a binary latent sign which, at the positive/negative
nonzero source, agrees with its sign with probability (1+s)/2. At the
zero source it is fair. Its marginal is fair and its exact information
and residual variance are

`I_s=p[log2-h((1+s)/2)]`,
`V_s=1-p s^2=689/3200` at p=31/32.

Therefore, for every lambda,

`J_lambda(nu_p)<=I_s+lambda V_s`.                           (7)

This is a single explicit admissible channel, not a variational upper
estimate based on a numerical optimizer.

Write `b=p log2+t(1-sqrt(p))`. For 0<lambda<=9/10, equations (3) and
(7) give the following lower bound on the propagated exponent:

`b+e_lambda(nu_p)+d`
` >= b-I_s-lambda V_s`
` >= p h(19/20)+4(1-sqrt(31/32))-(9/10)V_s`.                 (8)

For 9/10<=lambda<=4, use (4)--(7):

`b+e_lambda(nu_p)+d`
` >= p h(19/20)+4(1-sqrt(31/32))+g`
`       +(1+log(2lambda))/2-lambda V_s`.                     (9)

The last expression is concave in lambda. Its minimum on the declared
closed interval is therefore attained at one of its TWO endpoints.
Equations (8)--(9) reduce the entire continuum of temperatures to three
explicit scalar inequalities.

## 4. Exact rational enclosure of the three tests

The verifier
`computations/continued_feedback_fixed_temperature_majorant_obstruction_2026_09_06.py`
imports only the rational logarithm/entropy/square-root enclosure routines
from the independently replayed latent-envelope certificate. It uses no
optimizer or floating transcendental values. Its saved output is
`computations/results/continued_feedback_fixed_temperature_majorant_obstruction_2026_09_06.json`.

From the repository root, the exact replay command is:

```sh
.venv/bin/python computations/continued_feedback_fixed_temperature_majorant_obstruction_2026_09_06.py --output computations/results/continued_feedback_fixed_temperature_majorant_obstruction_2026_09_06.json
```

It uses only the Python standard library plus the neighboring certificate
module; Python automatically puts that script directory on its import
path. The output is a reproducible generated result, not an optimization
checkpoint. All three inequalities and their strict positivity are
asserted in exact rational arithmetic before any display conversion.

The three exact lower endpoints are:

| Region/test | Rational lower endpoint | Decimal display |
|---|---|---|
| lambda<=9/10 | 19688465594483/320000000000000 | .061526454982759375 |
| high interval, lambda=9/10 | 24879821494403/320000000000000 | .07774944217000937 |
| high interval, lambda=4 | 49954601778723/320000000000000 | .15610813055850936 |

The Gaussian potential is evaluated from

`r=(sqrt(1+16t^2)-1)/(4t)`,
`g=-t(1-r)+(1/4)log(1-r^2)`.

For its lower endpoint, use the lower endpoint of r in the increasing
linear term, and the upper endpoint of r in the decreasing logarithmic
term. This gives the rational lower value
`-622136276211/800000000000` for g. Entropy is enclosed from below and
sqrt(p) from above in (8)--(9). The minimum of the three resulting
rational values is (2).

## 5. The obstruction is relevant to finite-alphabet terminal trees

No actual Gaussian alphabet is needed. Let nu_r be the normalized sum
of 2^r independent nu_p variables. This is an admissible finite terminal
law along the iid policy, with variance exactly one. Its W2 distance
omega_r from the unit Gaussian tends to zero: the central limit theorem
gives weak convergence and the uniformly bounded fourth moments give
uniform integrability of squared values.

The rate-distortion Lagrangian is uniformly continuous along this
approximation for 0<lambda<=t. Couple X~nu_r with X' Gaussian at their
W2 distance, and, given X, pass it through a nearly optimal finite latent
channel. The chain L--X--X' contracts mutual information. Using the same
reproduction for X' increases squared error by at most
`2 omega_r sqrt(D)+omega_r^2`, with D<=1. Consequently

`J_lambda(N(0,1)) <= J_lambda(nu_r)
                         +t(2 omega_r+omega_r^2)`.          (10)

Also `Phi_t(nu_r)>=g` by the already proved Gaussian lower bound for
every finite symmetric unit-variance source. Thus even if the additive
correction d_r is required to dominate the terminal ONLY at the zero
state and this iid terminal state, its propagated root exponent obeys

`b+e_(lambda_r)(nu_p)+d_r`
` >= .061526454982759375-t(2 omega_r+omega_r^2)`               (11)

for every possibly r-dependent choice 0<lambda_r<=t. The zero terminal
is attainable at every positive depth by a diagonal pairing and retaining
its minus child. The iid terminal is attainable at that same depth by
all iid pairings. A uniform terminal error must cover both possible laws,
even though they arise from different policies.

The error in (11) tends to zero uniformly in temperature. Hence this is
a scalable large-depth obstruction to the uniform-additive-error method,
not merely a test on an inadmissible limiting state. It does not rule out
state-dependent errors, paying different errors on different branches,
an adaptive temperature with controlled information cost, or a genuinely
tighter full Bellman upper certificate.
