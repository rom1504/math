# Wave 44: intrinsic-load increments and a phase-bottleneck wall

## Status

This note gives two exact increment calculi for the intrinsic harmonic load:
state flips are sums of two-coordinate curvature increments, and coordinate
deletion is governed by a symmetric relative-entropy Hessian.  Neither
calculus supplies the project bound.  In fact, the most direct
Efron--Stein/heat-bath implementation is falsified already by the exact
minimizer `A8,m=4`: its load has macroscopic variation between endpoint
phases while every local edge joining different phases crosses an
exponentially suppressed barrier.

The surviving concrete conjecture is a **minimizer-specific integrated
endpoint domination**

```math
\mathcal J_L^2\le C C_V(1),
```

not a local Poincare inequality.  It passes all current exact-minimizer
audits and the full low-temperature `A8` boundary-layer calculation, but it
is open and is false for the abstract Wave-36 matched endpoint.

Checks are in
`/home/math/quadra/tmp/harmonic_load_increments_r44_check.py`.

## 1. Exact state-flip increment

For a directed vertex flip from `d` to `d^i`, put

```math
\omega_i(d)=\log\frac{\nu(d^i)}{\nu(d)},
\qquad
\chi_i(d)=g(d^i)-g(d).
```

The conditional curvature and KL are

```math
v_i(t,d)
=\frac{\chi_i(d)^2}
{4\cosh^2((\omega_i(d)+t\chi_i(d))/2)},
\qquad
k_i(s,d_{-i})=\int_0^s t v_i(t,d)\,dt.
\tag{R44.1}
```

Both quantities are unchanged when `d` is replaced by `d^i`.  Therefore,
for the intrinsic vertex load `L_s(d)=sum_(i in V_*) k_i(s,d_-i)`, flipping
a different vertex `j` gives exactly

```math
\boxed{
L_s(d^j)-L_s(d)
=\sum_{i\in V_*\setminus\{j\}}
\int_0^s t\{v_i(t,d^j)-v_i(t,d)\}\,dt.
}
\tag{R44.2}
```

The `i=j` summand is exactly zero; it is not bounded and discarded.

For the quadratic Gibbs base,

```math
\boxed{
\omega_i(d^j)-\omega_i(d)=8\beta a_{ij}d_{ij},
\qquad
\chi_i(d^j)-\chi_i(d)=\Delta_j\Delta_i g(d).
}
\tag{R44.3}
```

The mixed score difference is symmetric in `i,j`.  Thus every local
increment of `L_s` splits into a known single-entry base perturbation and a
mixed discrete Hessian of the matched log likelihood.  This is the exact
interface at which minimality or completion structure would have to enter a
state-flip proof.

For reference, if

```math
v(\omega,\chi;t)=\frac{\chi^2}{4}
\operatorname{sech}^2\frac{\omega+t\chi}{2},
```

then

```math
\partial_\omega v=-v\tanh\frac{\omega+t\chi}{2},
```

```math
\partial_\chi v
=\frac\chi2\operatorname{sech}^2\frac{\omega+t\chi}{2}
-\frac{t\chi^2}{4}\operatorname{sech}^2
\frac{\omega+t\chi}{2}\tanh\frac{\omega+t\chi}{2}.
\tag{R44.4}
```

These formulas give pointwise Lipschitz bounds, but summing them loses both
resonant suppression and cross-coordinate cancellation.  The `A8` wall
below shows that inserting (R44.2) into a uniform local Poincare inequality
cannot be the missing proof.

## 2. Exact deletion Hessian

There is a separate, averaged deletion calculus.  For a chart-coordinate
set `B`, let

```math
F_s(B)=D((\mu_s)_B\Vert\nu_B).
```

For `i in B`, the chain rule gives the nonnegative deletion increment

```math
C_i^B(s)=F_s(B)-F_s(B\setminus\{i\})
=\mathbb E D(\mu_s(D_i\mid D_{B\setminus i})
\Vert\nu(D_i\mid D_{B\setminus i})).
\tag{R44.5}
```

At the full coordinate set this is `E_(mu_s) k_i(s)`.  The symmetric mixed
deletion Hessian is

```math
\boxed{
H_{ij}(s)
=F_s(V)-F_s(V\setminus i)-F_s(V\setminus j)
+F_s(V\setminus\{i,j\}).
}
\tag{R44.6}
```

Writing `R=V\setminus\{i,j\}`, direct expansion gives

```math
\boxed{
H_{ij}
=I_{\mu_s}(D_i;D_j\mid D_R)
-\mathbb E_{\mu_s}
\log\frac{\nu(D_i\mid D_j,D_R)}{\nu(D_i\mid D_R)}.
}
\tag{R44.7}
```

The first term is nonnegative; the base-interaction term has no sign.  If
`L_s^{(-j)}` denotes the local load formed after marginalizing coordinate
`j` from both laws, then

```math
\boxed{
\mathbb E_{\mu_s}L_s
-\mathbb E_{(\mu_s)_{-j}}L_s^{(-j)}
=C_j^V(s)+\sum_{i\in V_*\setminus\{j\}}H_{ij}(s).
}
\tag{R44.8}
```

Consequently deletion is not monotone.  At `s=1`, the exact-minimizer audit
finds both signs of `H_ij` on `A9,m=4` and `A9,m=7`; for the latter the range
is approximately `[-0.010894,0.0029085]`.  The total deletion increment also
has both signs, approximately `[-0.010610,0.0036307]`.  Therefore relative-
entropy submodularity, supermodularity, and monotone leave-one-out
self-bounding are all unavailable without an additional signed correction.

## 3. Exact `A8,m=4` tropical phase certificate

Let `a(d)` be the completion exponent from (10.1146), and let `c(d)` be the
sum of leading reciprocal completion multiplicities.  Exact enumeration for
`A8,m=4` gives endpoint line maximum

```math
\max_d\{E_d-a(d)\}=12.
```

There are twenty active states.  Their exact triples
`(E,a,c)` and tropical load slopes are

| count | `(E,a,c)` | `lim_(beta->infinity) L_1/beta` |
|---:|:---:|---:|
| 7 | `(20,8,1)` | `4/7` |
| 7 | `(16,4,1/6)` | `4/7` |
| 1 | `(20,8,1)` | `0` |
| 5 | `(16,4,1/6)` | `0` |

The fourteen positive-load states form seven active cube edges, each pairing
one state of each energy.  The other six active states are isolated in the
active-state graph.  Every cube edge whose two endpoint load slopes differ
has smaller endpoint line at least four below the maximum.  These are exact
integer/rational enumeration statements.

It follows directly that the endpoint active law gives mass `49/60` to the
positive-load phase and `11/60` to the zero-load phase, and hence

```math
\boxed{
\operatorname{Var}_{\mu_1}(L_1)
\sim\frac{11}{225}\beta^2.
}
\tag{R44.9}
```

No active cube edge sees the leading difference between these two phases.
Thus the heat-bath Dirichlet energy of `L_1` is `o(beta^2)`, while its
variance is asymptotic to (R44.9).  This rigorously falsifies every
temperature-uniform estimate

```math
\operatorname{Var}_{\mu_s}(L_s)
\le C\mathcal E_{\mu_s}(L_s)
```

even for exact minimizers and even after the orientation edges are included.
Exact quadrature indicates the sharper asymptotic ratio
`Theta(exp(4 beta))`; the checker values at `beta=2,4,8` are respectively
`279.6`, `8.32e5`, and `7.39e12`.  The exponential equivalent is labelled
**tropical/numerical**, since only the divergence, not the complete
subleading expansion, is needed for the rigorous wall.

The same phase description gives the complete leading boundary layer.  Put
`u=beta(1-s)`.  On the twenty active states the limiting partition,
bad-state probability inside a paired component, and total paired-phase
mass are

```math
Z(u)=8e^{8u}+2e^{4u},
\qquad p(u)=\frac1{1+6e^{4u}},
```

```math
P(u)=\frac{7e^{8u}+(7/6)e^{4u}}
{8e^{8u}+2e^{4u}}.
```

The load divided by `beta` tends to `4p(u)` on the paired phase and to zero
on the singleton phase.  Dominated finite-state Laplace asymptotics then
give the **Verified** limit

```math
\boxed{
\frac{\mathcal J_L^2}{\beta}
\longrightarrow
\int_0^\infty16P(u)(1-P(u))p(u)^2\,du
}
```

```math
\boxed{
=7\log\frac76+\frac{217}{36}\log\frac45+\frac{49}{180}
=0.0062172411468\ldots .
}
\tag{R44.10}
```

Also

```math
\frac{C_V(1)}\beta\longrightarrow\frac7{15},
\qquad
\frac{\mathcal J_L^2}{C_V(1)}
\longrightarrow0.0133226596003\ldots .
\tag{R44.11}
```

Thus `A8` destroys local Efron--Stein control but **does not** obstruct the
integrated project target or endpoint domination.

## 4. The narrow surviving theorem

The preceding wall says that any successful proof must charge phase-to-phase
load variation nonlocally.  The narrowest useful replacement found here is

```math
\boxed{
\mathcal J_L^2
=\int_0^1\frac{\operatorname{Var}_{\mu_s}(L_s)}s\,ds
\le C\,C_V(1),
}
\tag{R44.12}
```

uniformly over the project family of exact minimizers, fixed-density target
sizes, and prescribed temperature.  This is an **Open target**, not a
consequence of the flip or deletion identities.  It is materially stronger
than merely renaming `J_L=O(a_n)`: the endpoint vertex cost is an existing,
separate input, and (R44.12) would supply the migration norm for free.

Together with (10.1153), (R44.12) gives

```math
\mathcal A_V\le\sqrt{\mathscr H\,C C_V(1)},
```

so endpoint cost `C_V(1)=O(a_n)`, restoring, and a separate orientation
bound close the quadratic bootstrap (10.1156).  Adjacent-selector Hellinger
remains separate.

Finite evidence is consistent with (R44.12).  The largest displayed ratios
`J_L^2/C_V(1)` are about `0.01264` on `A8,m=4,beta=8`, `0.00962` on
`A9,m=7`, and `0.0140` on the tested `A9,m=4` grid.  Equation (R44.11)
checks the only growing low-temperature finite phase found here.

The claim is necessarily minimizer/signing-specific.  On the positive
canonical Wave-36 transient example, the ratio grows numerically as

```text
M=100:    0.0391
M=1000:   0.1068
M=10000:  0.3571,
```

while adverse migration divided by endpoint cost grows to `145.98`.
The flat endpoint causes no issue: both sides of (R44.12) vanish.  The `A9`
rare edge also causes no issue because its contribution is exponentially
small under the actual interpolation law.

## 5. Conclusion

- State flips reduce exactly to the cross-curvature formula (R44.2), with
  quadratic increment and mixed score Hessian (R44.3).
- Deletion reduces exactly to the sign-indefinite entropy Hessian
  (R44.6)--(R44.8); exact `A9` data rule out monotone deletion.
- Uniform local Poincare/Efron--Stein control of `L_s` is rigorously false on
  `A8`.  The failure is a phase bottleneck, not the orientation coordinate.
- The integrated norm itself remains well behaved on that wall and is
  asymptotically a small constant times endpoint vertex cost.
- The precise next harmonic lemma is the minimizer-specific endpoint
  domination (R44.12), or a phase-aware theorem of the same strength.  A
  proof must compare separated interpolation phases globally; local state
  gradients and raw reveal costs cannot do it.
