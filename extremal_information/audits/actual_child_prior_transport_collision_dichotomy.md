# A child-prior transport gap converts overlap into collision tails

**Status.** Rigorous conditional theorem for every centrally symmetric
rank-one child prior, hence for the actual optimizing-child prior.  It
reduces the strong-channel overlap obstruction to one child-only
transport/subgaussian modulus.  A uniform positive transport gap forces an
extensive posterior-collision tail on positive actual path mass.  Failure
of the gap certifies a one-query MGF-violating direction, but does not by
itself make that direction macroscopic or operationally concise.  The
theorem does not decide which branch actual minimizing children occupy.

## 1. The child-only coherence modulus

Let `mu` be a centrally symmetric law on rank-one sign matrices
`Q=XY^T`, put `d=mn`, and let

```math
\mathcal C_\mu=\operatorname {conv}(\operatorname {supp}\mu).
```

For `0<c<1/2`, say that `mu` has rank-one transport gap `c` if, for every
nonzero `A in C_mu`,

```math
\boxed{
\log E_\mu\exp\left\{{\langle A,Q\rangle\over2\sqrt d}\right\}
\le\left({1\over2}-c\right){\|A\|_F^2\over\sqrt d}.}         \tag{PT.1}
```

This is a property of the two child Gibbs laws before a bridge is sampled.
It asks a family of one-dimensional exponential-moment queries and contains
no inverse escort, posterior table, or target-order optimum.  Whether the
modulus or an optimizing query can be found without full convex-hull
information is not asserted.

The two extremal calibration examples are exact.  The uniform projective
rank-one prior satisfies (PT.1) with `c=1/3` by the determinant MGF bound
(37.294).  A two-antipodal-atom prior fails every fixed positive gap for all
large `d`, already along arbitrarily small multiples of its atom direction.

### Lemma PT.1 (transport gap gives entropy curvature)

If (PT.1) holds and `nu<<mu` has barycenter `M=E_nu Q`, then

```math
\boxed{D(\nu\Vert\mu)\ge c{\|M\|_F^2\over\sqrt d}.}          \tag{PT.2}
```

*Proof.*  The barycenter belongs to `C_mu`.  Entropy duality with the test
`<M,Q>/(2sqrt(d))` gives

```math
D(\nu\Vert\mu)
\ge {\|M\|_F^2\over2\sqrt d}
 -\log E_\mu e^{\langle M,Q\rangle/(2\sqrt d)}.
```

Apply (PT.1).  The case `M=0` is immediate.  `square`

### Lemma PT.1a (a checkable necessary susceptibility bound)

Suppose `mu` has full rank-one support, as every finite-temperature actual
child Gibbs law does.  Let `Sigma_Q` be the covariance operator of
`vec(Q)`.  A transport gap `c` necessarily implies

```math
\boxed{\|\Sigma_Q\|_{op}\le(4-8c)\sqrt d.}       \tag{PT.2a}
```

For a pure product law, write

```math
C_X=E(XX^{\mathsf T}),\qquad C_Y=E(YY^{\mathsf T}).
```

Then

```math
\boxed{
\Sigma_Q=C_Y\otimes C_X,qquad
\|\Sigma_Q\|_{op}=\|C_X\|_{op}\|C_Y\|_{op}.}   \tag{PT.2b}
```

*Proof.*  Full rank-one support spans the matrix space, and its centrally
symmetric convex hull contains a neighborhood of zero.  Hence for every
matrix `H`, `epsilon H` belongs to `C_mu` for all sufficiently small
positive `epsilon`.  Expand (PT.1) at zero.  Since `E_muQ=0`,

```math
\log E_\mu e^{\epsilon\langle H,Q\rangle/(2\sqrt d)}
={\epsilon^2\over8d}
 \operatorname {Var}_\mu\langle H,Q\rangle+o(\epsilon^2).
```

Comparing quadratic coefficients gives
`Var <H,Q><=(4-8c)sqrt(d)||H||_F^2`, proving (PT.2a).
Independence of the two child factors gives
`E(Q_(ij)Q_(kl))=E(X_iX_k)E(Y_jY_l)`, which is (PT.2b).  `square`

Thus a product susceptibility exceeding `4sqrt(d)` rules out every fixed
positive transport gap and supplies a rank-one violating direction from
the two top child-correlation eigenvectors.  The converse is false without
a nonlinear concentration hypothesis: bounded covariance does not control
the finite exponential moment in (PT.1).

The exact actual zero-bridge law is a two-sector mixture, not generally one
pure product.  With the notation (LE.2), its covariance is

```math
\boxed{
\Sigma_Q=\sum_{s=\pm1}\pi_s^{(\epsilon)}
  \bigl(C_{D,\epsilon s}\otimes C_{A,s}\bigr).}  \tag{PT.2c}
```

Therefore `||Sigma_Q||_op<=(4-8c)sqrt(d)` is the correct necessary actual-
child test.  More simply, if either sector satisfies

```math
\pi_s^{(\epsilon)}
 \|C_{A,s}\|_{op}\|C_{D,\epsilon s}\|_{op}>4\sqrt d,
```

then the gap fails and the two sector top eigenvectors give a rank-one
violating direction.  No factorization of the complete sector mixture is
claimed.

### Lemma PT.1b (sector linear subgaussianity is sufficient)

For an even sign law `nu` on `{-1,1}^k`, define its linear subgaussian
proxy by

```math
\sigma^2(\nu)=\inf\left\{s:
 \log E_\nu e^{\langle u,X\rangle}
 \le {s\|u\|_2^2\over2}\quad\hbox{for every }u\right\}.
```

For the actual two-sector law, put

```math
\kappa_*=max_{s=\pm1}
 \sigma^2(\mu_{A,s})\sigma^2(\mu_{D,\epsilon s}).
```

If `kappa_*<4`, then for every `A in C_mu`,

```math
\boxed{
\log E_\mu e^{\langle A,Q\rangle/(2\sqrt d)}
\le {\kappa_*\|A\|_F^2
 \over8d(1-\kappa_*/4)}.}                       \tag{PT.2d}
```

In particular, (PT.1) holds with `c=1/4` whenever

```math
\sqrt d\ge {2\kappa_*\over4-\kappa_*}.          \tag{PT.2e}
```

*Proof.*  Condition first on one factor and apply its linear subgaussian
bound.  Gaussian linearization of the resulting quadratic exponential,
followed by the other factor's bound, gives in each sector

```math
\log Ee^{\theta X^{\mathsf T}AY}
\le-{1\over2}\log\det(I-\kappa_*\theta^2A^{\mathsf T}A).
```

Every sector obeys the same upper bound, so their mixture does too.  Bound
the log determinant by
`kappa_* theta^2||A||_F^2/[2(1-kappa_*theta^2||A||op^2)]`.
For `theta=1/(2sqrt(d))` and `A in C_mu`, one has
`||A||op^2<=d`; this proves (PT.2d).  Comparing its coefficient with
`1/(4sqrt(d))` proves (PT.2e).  `square`

Thus an actual child log-Sobolev or transportation inequality which yields
the displayed Herbst bounds with `kappa_*<4` would settle the extensive-
collision branch.  No current optimizer identity supplies such a uniform
strong-mixing estimate at the dense physical temperature.

The covariance test is sharply insufficient.  Let
`nu_m=(1-p_m)U_m+(p_m/2)(delta_1+delta_(-1))`, with
`p_m=exp(-km)`, and define `nu_n` similarly.  Its covariance operator norm
is `1+(m-1)p_m=1+o(1)`, and its off-diagonal correlations and largest atoms
are exponentially small.  Yet for the product rank-one law, `C=11^T`, and
`A=C/2`,

```math
\log Ee^{\langle A,Q\rangle/(2\sqrt d)}
\ge {\sqrt d\over4}-kN-O(1),
```

whereas the right side of (PT.1) even at `c=0` is `sqrt(d)/8`.  At
comparable splits this violates every positive gap when
`k<sqrt(gamma_0)/8`.  Rare coherent tails, invisible to covariance, are
exactly the nonlinear information retained by (PT.1).  This spike law is
not an actual quadratic Gibbs minimizer and need not satisfy the strongest
conditional-spread theorem.

## 2. Positive overlap forces a positive collision tail

Let `q` be any bridge law, `M(B)=E[Q|B]` the complete posterior mean, and
`r(B)` the deleted-edge cavity matrix evaluated at the complete bridge.
Write `rho=tanh(t)` and assume

```math
{1\over d}E_q\|r(B)\|_F^2\ge\eta,
\qquad a:={\eta\over2}-4\rho^2>0.                \tag{PT.3}
```

### Theorem PT.2 (transport--collision dichotomy)

If `mu` has transport gap `c`, then

```math
\boxed{
q\left\{\log K_0(B)\ge {ca\sqrt d\over2}\right\}
\ge {a\over2-a}.}                               \tag{PT.4}
```

For every deleted edge, on the same event,

```math
\boxed{\log K_e(B_{-e})\ge {ca\sqrt d\over2}-4t.}           \tag{PT.5}
```

*Proof.*  Exact insertion gives
`||M-r||_F^2<=4rho^2d`, and therefore

```math
E_q\|M\|_F^2\ge {1\over2}E_q\|r\|_F^2-4\rho^2d\ge ad.
```

Since every posterior barycenter lies in the convex hull of sign matrices,
`Z=||M||_F^2/d` lies in `[0,1]`.  If
`p=q{Z>=a/2}`, then

```math
a\le E_qZ\le p+(1-p)a/2,
```

so `p>=a/(2-a)`.  Lemma PT.1 and `D_2>=D` imply on this event

```math
\log K_0=D_2(\mu_B\Vert\mu)
\ge c\|M\|_F^2/\sqrt d\ge ca\sqrt d/2.
```

Finally GC.3 gives `K_e>=e^(-4t)K_0` pointwise for every edge.  `square`

Thus a positive transport gap does more than make the annealed collision
mean exponential: it violates every subexponential collision threshold on
a fixed positive amount of the declared bridge mass.

## 3. Application to the actual negative path

For comparable actual contracted-temperature minimizing children at any
fixed

```math
\beta>\beta_{\rm BG}(\gamma_0),
```

Theorem 37.56 supplies a constant `eta_*>0` such that the negative-path
mixture

```math
\bar q={1\over\lambda}\int_{-\lambda}^0q_s\,ds
```

obeys

```math
{1\over d}E_{\bar q}\|r\|_F^2\ge\eta_*+o(1).   \tag{PT.6}
```

Consequently, if the corresponding actual child priors have a uniform
transport gap `c>0`, Theorem PT.2 proves that for some constants
`c_1,c_2>0`,

```math
\boxed{
(\bar q\otimes U_{\rm edge})
 \{\log K_e\ge c_1N\}\ge c_2}                  \tag{PT.7}
```

for all large `N`.  This is a genuine actual-minimizer obstruction to the
tail condition (37.291), conditional only on the child-prior property
(PT.1).

If no uniform gap exists, then for every fixed `c>0` along a subsequence
there is one `A_N in C_(mu_N)` with

```math
\boxed{
\log E_{\mu_N}e^{\langle A_N,Q\rangle/(2\sqrt d)}
>\left({1\over2}-c\right){\|A_N\|_F^2\over\sqrt d}.}         \tag{PT.8}
```

This certifies the existence of one child-only MGF-violating query.  It is
absent from the uniform-factor model, while a two-antipodal-atom prior
strongly violates (PT.1).  The witness may have vanishing norm or excess,
require `Theta(N^2)` real coordinates, and be costly to discover from the
convex hull.  Thus (PT.8) is not yet a coherent low-dimensional phase, a
reusable response quotient, or a target-relevance theorem.

## 4. What the dichotomy changes

The raw-overlap decay target is false in the strong-channel regime.  The
new child-only decision problem is:

> prove a fixed transport gap (PT.1), which forces a positive-mass
> extensive collision obstruction; or strengthen a violating query
> (PT.8) to have quantitative norm/excess, a concise child-generated
> description, and target relevance.

This is semantically lower-output-information than the complete inverse-
escort bridge landscape: it is a modulus of the zero-bridge child prior.
No operational description-complexity separation is proved.  It is not yet
known to be mathematically easier for exact minimizing signings, and
neither branch alone supplies the target-reach inequality required for a
Level-6 recurrence.
