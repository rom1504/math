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
