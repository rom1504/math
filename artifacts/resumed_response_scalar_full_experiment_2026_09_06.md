# Scalar full-center optimization: reproducible diagnostic and one exact birth

Date: 2026-09-06. Numerical status: FLOATING DIAGNOSTIC ONLY. No displayed
decimal in this file is an exact signing lower certificate or an upper
bound on the scalar optimum. The actual variational and realization
theorems are in
`resumed_response_scalar_sobolev_full_center_variation_2026_09_06.md`.

## 1. Constrained full-functional optimization

The script
`computations/resumed_response_scalar_full_optimizer_2026_09_06.py`
optimizes the exact one-dimensional full-center functional, not merely
the correlation c. Nonconstant coefficients are parameterized as
x_r=sqrt(r)g_r, so D(g)=sum x_r^2<=1 and
g_0=sqrt(1-sum x_r^2/r). SLSQP uses the analytically derived coefficient
gradient; alpha is varied as well. Central Gauss--Legendre integration
is independently replayed with more nodes.

The best values from the initial degree progression were:

| Highest even degree | Full functional | Threshold alpha |
| --- | --- | --- |
| 20 | 0.4247419614 | 0.73791672 |
| 40 | 0.4260517763 | 0.74882579 |
| 80 | 0.4266940442 | 0.74558998 |
| 160 | 0.4269472615 | 0.74637089 |
| 320 | 0.4270169877 | 0.74681210 |

Multiple nearby starts agreed at each degree; separate broad random
starts at degree 80 were also run. The derivative constraint was active
to floating tolerance, consistent with the proved criticality theorem.
Some raw optimizer outputs have D slightly ABOVE one by about 1e-12,
which is another reason not to treat them as exact feasible certificates.
The one-birth diagnostic below explicitly shrinks the derivative-energy
parameter norm by 1e-8 before evaluation, but still requires exact audit
before it could count as a certificate.

Saved result files:

* `computations/results/resumed_response_scalar_full_optimizer_2026_09_06.json`
* `computations/results/resumed_response_scalar_full_optimizer_high_degree_2026_09_06.json`
* `computations/results/resumed_response_scalar_full_optimizer_broad_starts_2026_09_06.json`

The first two runs' independent quadrature differences were about 1e-15.
Random directional finite differences of the analytic coefficient
gradient are included in the broad-start result file. These are checks
of the numerical implementation, not global optimization certificates.

The tested branch is well below the existing 0.4314603928237005 frontier.
Increasing polynomial degree or finding the same local optimum from
many starts does NOT prove that the scalar class has no better branch.
No exact candidate certification is warranted solely by these values.

## 2. A first unrestricted birth still admits a one-dimensional line search

There is a useful exact reduction, independent of the numerical outcome.
At the original scalar response, K(v)=lambda g(v)+gamma H(v) depends
only on V. Write

    B=E[H phi(K/sqrt(t))]/sqrt(t),
    w(v)=H(v)s(K(v))-2BK(v),
    A_1=<w,g>, nu=||w-A_1g||_2,
    h_U=(w-A_1g)/nu.

Assume nu>0. Then U=U h_U is a standard Gaussian independent of V and
the full gradient is A_1 V+nu U. This is a legitimate causal extension
of the scalar core. The old response need not be measurable in this
two-coordinate pair because its W coordinate can have another Gaussian
component; the union frame is used for the exact feasible convex line.

Conditional on V=v, let m=A_1v and b=max(Psi(K(v),t)-B,0). The full
linearized maximizer has conditional support and sign mean

    q_new(v)=Phi((m-b)/nu)+Phi((-m-b)/nu),
    f_mean(v)=Phi((m-b)/nu)-Phi((-m-b)/nu).

Its first-chaos moments are the one-dimensional integrals

    a_V=E[V f_mean(V)],
    a_U=E[phi((b-m)/nu)+phi((b+m)/nu)],
    p_new=E q_new(V).

No other first-chaos direction occurs. Therefore its inverse is again
a function of V alone:

    K_new(v)=a_V g(v)+a_U h_U(v).

Along the feasible convex line of conditional pairs, put

    K_theta=(1-theta)K+theta K_new,
    p_theta=(1-theta)(1-p)+theta p_new,
    t_theta=p_theta-||K_theta||_2^2.

The exact full-center objective is

    J_theta=E_V[(1-(1-theta)(1-H(V))-theta q_new(V))
                                  Psi(K_theta(V),t_theta)].

All norms and cross moments of K and K_new are known from g,h_U
orthonormality and the one-dimensional integrals above. Thus ONE
unrestricted birth and its full line search are still exactly reducible
to one outer Gaussian integration. No conditional law of a matrix output
has been assumed.

The next birth does not automatically enjoy the same reduction: the
new center mask depends on U, so its inverse gradient need no longer
be a function of V alone. The one-step collapse must not be iterated
without accounting for that new dependence.

## 3. Reproducible first-birth diagnostic

The script
`computations/resumed_response_scalar_first_birth_2026_09_06.py`
implements Section 2. It integrates positive V separately on [0,alpha]
and [alpha,12], exploiting parity; it uses exact inverse-feature moments
for squared norms rather than integrating a high-degree polynomial
squared on a finite interval. It replays with twice as many nodes.
Results are saved in
`computations/results/resumed_response_scalar_first_birth_2026_09_06.json`.

Starting from the slightly shrunk degree-320 inverse, the diagnostic gives:

    old full value      0.4270169872828883
    full gradient gap   0.0075727184518151
    best theta          0.5764057585
    new full value      0.4291454376404552
    gain                0.0021284503575684
    full-step value     0.4281112724690156

The difference between the 384- and 768-node best values is about
1.13e-12. The finite-difference derivative of the line at theta=0 agrees
with the independently evaluated full gap to the expected first-order
finite-difference error. Neither the finite cutoff nor the quadrature
error has been enclosed by exact arithmetic, so these remain diagnostics.

The unrestricted birth gives a substantial increase within this simpler
class, and the optimal line step materially outperforms taking theta=1.
It still does not beat the banked 21-anchor certificate. This is useful
negative evidence for prioritizing that richer geometry, not a proof of
a scalar or one-birth ceiling.
