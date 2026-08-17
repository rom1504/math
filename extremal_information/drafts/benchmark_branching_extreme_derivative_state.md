# Orthogonal benchmark: derivative-martingale state for branching extremes

Status: primary-source import with an explicit branching law.  Portfolio verdict:
**promote as a scoped rare-event benchmark, not as a finite-port theorem**.

## 1. Imported extremal-process theorem

Consider a branching random walk in the boundary case

```math
\mathbb E\sum_{|u|=1}1>1,\qquad
\mathbb E\sum_{|u|=1}e^{-V(u)}=1,\qquad
\mathbb E\sum_{|u|=1}V(u)e^{-V(u)}=0.                 \tag{BD.1}
```

Write `X=sum_(|u|=1)e^(-V(u))` and
`X_tilde=sum_(|u|=1)V(u)_+e^(-V(u))`.  Assume that the offspring displacement
process is non-lattice and

```math
\mathbb E\sum_{|u|=1}V(u)^2e^{-V(u)}<\infty,\qquad
\mathbb E[X(\log_+X)^2]<\infty,\qquad
\mathbb E[\widetilde X\log_+\widetilde X]<\infty.       \tag{BD.1a}
```

Put

```math
W_n=\sum_{|u|=n}e^{-V(u)},\qquad
Z_n=\sum_{|u|=n}V(u)e^{-V(u)},\qquad
a_n={3\over2}\log n.                                  \tag{BD.2}
```

Then `W_n -> 0` almost surely, while the derivative martingale converges to
`Z_infinity`, which is positive on survival.  Madaule's Theorem 1.1 states
that, conditional on survival,

```math
\left(\sum_{|u|=n}\delta_{V(u)-a_n+\log Z_\infty},Z_n\right)
\Longrightarrow (\mathcal E,Z_\infty),                         \tag{BD.3}
```

where the two limiting coordinates are independent.  The process `mathcal E`
is a decorated Poisson point process: its cluster leaders have intensity
`lambda e^x dx`, and every leader receives an independent copy of a fixed
decoration process `mathcal D` supported on `R_+`.

Equivalently, without the realization-dependent shift, the limiting process
`mathcal E_Z` is conditionally on `Z=Z_infinity` a Cox decorated Poisson
process whose leaders have intensity

```math
\lambda Z e^x\,dx.                                             \tag{BD.4}
```

For every nonnegative compactly supported test function `f`, its conditional
Laplace response is therefore exactly

```math
\mathbb E\left[e^{-\langle\mathcal E_Z,f\rangle}\mid Z\right]
=\exp\{-Z A_{\mathcal D}(f)\},                                \tag{BD.5}
```

where

```math
A_{\mathcal D}(f)
=\lambda\int_{\mathbb R}e^x
 \left(1-\mathbb E_{\mathcal D}
  e^{-\sum_{d\in\mathcal D}f(x+d)}\right)dx.                  \tag{BD.6}
```

Thus, once the model-dependent decoration law is fixed, the realized
environment enters **all Laplace-functional queries** through the single
positive scalar `Z`.  Aidekon's minimum theorem is the one-coordinate shadow:

```math
\mathbb P(M_n\ge a_n+x)\longrightarrow
\mathbb E\exp\{-C_*e^xZ_\infty\}.                             \tag{BD.7}
```

Primary sources are [Madaule, *Convergence in law for the branching random
walk seen from its tip*](https://arxiv.org/abs/1107.2543) and
[Aidekon, *Convergence in law of the minimum of a branching random
walk*](https://arxiv.org/abs/1101.1810).  The corresponding BBM decorated
process and its explicit decoration construction are proved by
[Aidekon--Berestycki--Brunet--Shi](https://arxiv.org/abs/1104.3738).

## 2. Branching composition and the prelimit qualification

Cut a finite-offspring branching random walk at generation `r`.  For a cut
particle `u`, let `W_m^(u),Z_m^(u)` be the additive and derivative martingales
of its descendant process, measured relative to `V(u)`.  The branching
property gives the exact finite-depth identity

```math
Z_{r+m}=\sum_{|u|=r}e^{-V(u)}
 \left(Z_m^{(u)}+V(u)W_m^{(u)}\right).                         \tag{BD.8}
```

The descendant processes are conditionally independent.  Sending `m` to
infinity and using `W_m^(u)->0` gives the smoothing transform

```math
\boxed{Z_\infty=\sum_{|u|=r}e^{-V(u)}Z_\infty^{(u)}}.          \tag{BD.9}
```

This law matches the Cox process exactly.  Translating the extremal process of
subtree `u` by `V(u)` multiplies its leader intensity by `e^{-V(u)}`;
independent decorated Poisson processes superpose by adding their intensities.
Consequently

```math
\biguplus_{|u|=r}
 \left(V(u)+\mathcal E_{Z_\infty^{(u)}}^{(u)}\right)
\ \stackrel d=\
\mathcal E_{\sum_{|u|=r}e^{-V(u)}Z_\infty^{(u)}}.             \tag{BD.10}
```

In logarithmic coordinates `S=-log Z`, the update is the soft minimum

```math
S=-\log\sum_{|u|=r}\exp\{-(V(u)+S_u)\}.                       \tag{BD.11}
```

This is a genuine asymptotic composition law.  It is not an exact scalar
finite-depth law: before taking the limit, (BD.8) requires the pair `(W,Z)`.
Madaule's equivalent law-level formulation is superposability of the
decorated exponential Poisson process under shifted independent unions.

## 3. Mapping to response state, and its boundary

Equation (BD.5) is a rank-one response family after taking negative logarithms:

```math
-\log L_Z(f)=Z A_{\mathcal D}(f).                              \tag{BD.12}
```

It is therefore a real compression phenomenon, but of a different kind from
the current finite-port theory:

- the port is an unbounded point configuration, not a fixed finite set;
- the statement is conditional and distributional, not equality of every
  realized deterministic future response;
- the scalar appears only after critical centering and martingale
  renormalization; and
- the decoration is a fixed infinite-dimensional law, not information encoded
  by the scalar state.

The current max-plus response-image theory would retain a full response table
or point configuration.  It does not by itself predict that a critical random
environment collapses, for this restricted query class, to a Cox intensity
given by a derivative martingale.

### Falsifier of an unrestricted scalar-state claim

Keep branch labels at the cut.  Two mass allocations `(z_1,z_2)` and
`(z'_1,z'_2)` can satisfy `z_1+z_2=z'_1+z'_2` while `z_1 ne z'_1`.  They give
the same unmarked scalar `Z` and hence the same unmarked extremal law.  A future
query that counts or forbids extremes descended from branch 1 has conditional
void probability depending on `z_1`, and distinguishes them.  Under such marked
futures the state must be the derivative-mass measure

```math
\sum_{|u|=r}e^{-V(u)}Z_\infty^{(u)}\delta_{\mathrm{mark}(u)},   \tag{BD.13}
```

not its total mass.  Likewise, if the model class is allowed to change its
decoration law, `Z` alone cannot answer cluster-multiplicity queries.

## 4. Portfolio verdict

**Promote, scoped.**  This supplies what a generic rare-event compactness
proposal lacks: a primary theorem with a unique limiting object, an exact
renormalized branching law, and a sharp recovery statement for a declared
query class.  Keep it orthogonal to the finite-port main theorem.  Promotion
does not license the claims that `Z` is an exact finite-time state, that it is
sufficient for genealogy-marked futures, or that derivative-martingale
compression holds outside the boundary-case universality class.
