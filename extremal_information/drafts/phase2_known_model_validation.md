# Phase 2 known-model validation: three different compression mechanisms

**Status.** Specialist report.  The Gaussian and max-plus results cited below
are imported classical facts.  The response-law classification,
Shapley--Folkman response theorem, and vector-balancing corollaries are proved
here.  This file does not claim literature novelty for the
Shapley--Folkman specialization.

## Executive verdict

REM/GREM, finite-state trellises, and vector balancing validate three
different notions that should not be conflated.

| Model | Small state | Exact for what? | What makes it close? |
|---|---|---|---|
| REM/GREM | Gaussian covariance hierarchy (and, for one asymptotic query, its convex-hull quotient) | the **law** of every deterministic-field response, not a quenched realization | Gaussianity, independence, and hierarchical exchangeability |
| homogeneous finite-state trellis | boundary max-plus kernel, eventually a periodic orbit modulo linear drift | every deterministic endpoint response | a finite separator; max-plus cyclicity under repeated identical composition |
| fixed-dimensional vector balancing | zonotope / convex response body | exactly all linear support queries and approximately every Lipschitz aggregate query | Shapley--Folkman: only `d` summands need remain fractional |

The third row gives the useful new theorem-level conclusion.  Convexification
can preserve more than linear support queries after many additive
compositions: its error for a nonlinear extremal query is bounded by the
feature dimension rather than the number of components.  This produces an
exactly composable, strict sub-landscape state for fixed-dimensional vector
balancing.  A same-zonotope construction with a linear discrepancy gap when
the dimension grows proves that the dimension hypothesis is structural.

The REM/GREM conclusion is mainly diagnostic.  Its familiar low-dimensional
order parameter is not a deterministic Parisi-like quotient of a realized
landscape.  It is an exact parameter of an **ensemble response law**.  The
distinction explains how it can retain zero-entropy extremes without
contradicting the project's deterministic counterexamples.

## 1. Gaussian landscapes: exact ensemble sufficiency, not a quenched roof

Let `Omega` be finite and let `G=(G_x)_(x in Omega)` be a Gaussian random
landscape with mean vector `m` and covariance matrix `C`.  For a deterministic
field `u in R^Omega`, put

```math
V_G(u)=\max_{x\in\Omega}\{G_x+u_x\}.                 \tag{KM.1}
```

Here the observation is the random process `u -> V_G(u)`, not only its
expectation.

### Theorem KM.1 (Gaussian response-law sufficiency)

The pair `(m,C)` determines the complete joint law of

```math
(V_G(u))_{u\in\mathbb R^\Omega}.                    \tag{KM.2}
```

Conversely, the joint law in (KM.2) determines `(m,C)`.  Hence mean and
covariance are the coarsest exact parameter, up to one-to-one recoding, for
the full-field response experiment **within the class of Gaussian ensemble
laws**.

For independent Gaussian landscapes on the same state set,

```math
(m,C)_{G+G'}=(m+m',C+C').                           \tag{KM.3}
```

For the direct-product landscape

```math
\widetilde G(x,y)=G(x)+G'(y),
```

the covariance is

```math
\widetilde C((x,y),(x',y'))=C(x,x')+C'(y,y'),       \tag{KM.4}
```

and separable response values add.  Thus Gaussian response-law parameters
have an exact composition algebra.

#### Proof

A finite Gaussian vector is determined in law by its mean and covariance,
and (KM.1) is a deterministic measurable function of that vector.  This
proves sufficiency for every finite collection of fields.

For the converse, fix `x`.  Let `u^(x,M)_x=0` and
`u^(x,M)_y=-M` for `y!=x`.  Then

```math
V_G(u^{(x,M)})\longrightarrow G_x
```

almost surely as `M -> infinity`.  Taking several such pinning fields at
once recovers every finite-dimensional joint law of `G`, hence its mean and
covariance.  Equations (KM.3)--(KM.4) are the elementary covariance rules
for independent Gaussian sums. `square`

### REM: a one-parameter ensemble response

For the REM, `|Omega_n|=N_n` independent energies have common law
`N(0,a n)`.  For an arbitrary deterministic field,

```math
\mathbb P\{V_G(u)\le t\}
=\prod_{x\in\Omega_n}
 \Phi\!\left({t-u_x\over\sqrt{an}}\right).          \tag{KM.5}
```

Thus the response law depends on `u` only through its multiset of values,
and on the disorder only through `a`.  If

```math
\log N_n=sn+o(n),
```

standard Gaussian tail bounds give

```math
{1\over n}\max_xG_x\longrightarrow\sqrt{2as}
\quad\text{in probability}.                         \tag{KM.6}
```

Indeed, the union bound makes the probability of exceeding
`(sqrt(2as)+epsilon)n` tend to zero.  At
`(sqrt(2as)-epsilon)n`, the expected number of exceedances grows
exponentially and independence makes the probability of no exceedance tend
to zero.  For `Omega_n={-1,+1}^n`, `s=log 2`.

This is the cleanest example in which an exponentially rare maximum is
preserved by a scalar state.  The reason is not deterministic compression:
the scalar specifies an i.i.d. Gaussian **law** and the extreme-value
operation is performed after that law is generated.

### GREM: covariance hierarchy and a distributional composition law

On a rooted tree, attach independent increments `Y_e` to edges and let a
leaf energy be their path sum.  If every node at level `j-1` has `b_j`
children, the level-`j` increments have law `nu_j`, and `F_j` is the CDF of
the maximum below a level-`j` node, then

```math
F_{j-1}(t)
=\left(\int F_j(t-y)\,d\nu_j(y)\right)^{b_j},
\qquad F_k(t)=\mathbf 1_{[0,\infty)}(t).             \tag{KM.7}
```

This follows by conditioning on the independent edge increments and child
subtrees.  Equation (KM.7) is an exact distributional composition algebra:
one response distribution per remaining subtree type replaces the realized
table of all leaf energies.

For a Gaussian GREM with the rooted hierarchy and branching data declared as
side information, the covariance is determined by the level variance
increments.  Theorem KM.1 therefore says that this hierarchy is an exact
parameter for every deterministic-field **response law**.  For `k` homogeneous
levels it uses `O(k)` scalar variance parameters, even though a realization
has exponentially many leaf energies.

This elementary finite-law statement aligns with, but is weaker than, the
rigorous asymptotic GREM theory.  Bovier and Kurkova define their Gaussian
process by

```math
\mathbb E X_\sigma X_{\sigma'}=A(d_N(\sigma,\sigma'))
```

and prove in Theorem 3.1 that the leading expected ground-state energy is an
explicit functional of the derivative of the convex hull `bar A`.  Their
Theorem 3.3 gives the limiting free energy from the same hull data.  See
[Bovier--Kurkova, *Derrida's Generalized Random Energy Models 2*,
Ann. IHP 40 (2004), 481--495](https://doi.org/10.1016/j.anihpb.2003.09.003).
The paper explicitly warns that subleading corrections need not be functions
of this hull alone.  The original models are due to
[Derrida, REM](https://doi.org/10.1103/PhysRevLett.45.79) and
[Derrida, GREM](https://doi.org/10.1051/jphyslet:01985004609040100).

The resulting classification is exact:

1. the full covariance hierarchy is an exact ensemble response parameter;
2. `bar A` is a still smaller sufficient quotient for declared leading
   ground-state/free-energy queries;
3. neither is an exact summary of a fixed disorder realization.

For the last point, expose the one-hot feature `phi(x)=e_x`.  The deterministic
upper response roof at the vertex `e_x` equals the realized value `G_x`.
Full pinning responses therefore recover the complete quenched landscape.
Gaussian self-averaging cannot be silently substituted for the uniform
zero-temperature guarantees required in an adversarial landscape class.

## 2. Finite-state trellises: an exact deterministic roof

For a `Q`-state homogeneous transition landscape with score matrix
`K in R^(Q times Q)`, a length-`t` path has endpoint kernel

```math
K^{\odot t},
\qquad
(K\odot L)_{ac}=\max_b(K_{ab}+L_{bc}).               \tag{KM.8}
```

The feature-growth draft proves that this kernel is exactly the upper
response roof for endpoint fields, is the coarsest deterministic statistic
for that query interface, and is nonexpansive under gluing.  This is not an
ensemble statement.

There is an additional known-model conclusion for repeated **homogeneous**
composition.  Since a finite real matrix is irreducible as a max-plus matrix,
the max-plus cyclicity theorem gives a maximum cycle mean `lambda`, an
admissible ultimate period `gamma` determined by the critical graph, and a
transient `T` such that

```math
K^{\odot(t+\gamma)}
=K^{\odot t}+\gamma\lambda\mathbf 1
\qquad(t\ge T).                                     \tag{KM.9}
```

Consequently every endpoint-field response is eventually periodic modulo
the linear drift `lambda t`.  After the transient, the complete family of
all lengths is represented by `gamma` boundary kernels, rather than by a
new state at every scale.  This is a genuine saturation law for the
query-generated feature algebra.

Equation (KM.9) is imported, not reproved here.  A precise general treatment
is B. De Schutter,
[“On the ultimate behavior of the sequence of consecutive powers of a matrix
in the max-plus algebra,” *Linear Algebra Appl.* 307 (2000),
103--117](https://doi.org/10.1016/S0024-3795(00)00013-6).

The boundary is equally precise.  An inhomogeneous product
`K_1 odot ... odot K_t` still has the fixed `Q^2`-entry kernel state, but no
eventual periodicity follows.  A future query that touches an internal state
instead of an endpoint bypasses the separator and invalidates the quotient.

## 3. Vector balancing: when a convex roof controls a nonlinear minimum

Let `v_1,...,v_n in R^d` and let `||.||` be any norm.  Define the discrete
signed-sum set and its convex hull

```math
S(V)=\left\{\sum_{i=1}^n\epsilon_i v_i:
              \epsilon_i\in\{-1,+1\}\right\},
\qquad
Z(V)=\operatorname{conv}S(V)
     =\sum_{i=1}^n[-v_i,v_i].                        \tag{KM.10}
```

`Z(V)` is the usual zonotope.  Vector balancing asks for a point of `S(V)`
near the origin or, more generally, near a target.

### Theorem KM.2 (support roof versus target-distance state)

The following statements hold.

**(a) Linear support queries convexify exactly.**  For
`phi(epsilon)=sum_i epsilon_i v_i` and zero base energy, the upper response
roof is the constant zero roof on `Z(V)`, and

```math
\max_{\epsilon}\langle\theta,\phi(\epsilon)\rangle
=h_{Z(V)}(\theta)
=\sum_i|\langle\theta,v_i\rangle|.                  \tag{KM.11}
```

Thus the zonotope is the exact minimal closed convex state for all linear
support queries.

**(b) Exact target queries retain the holes.**  Put

```math
D_V(t)=\min_{s\in S(V)}\|t-s\|.                     \tag{KM.12}
```

The complete function `D_V` determines `S(V)` because its zero set is
exactly `S(V)`.  Hence `S(V)`, up to equality, is the coarsest exact state
for all target-distance queries.  Under concatenation,

```math
S(V\sqcup W)=S(V)+S(W),
\qquad Z(V\sqcup W)=Z(V)+Z(W).                      \tag{KM.13}
```

Both states compose by Minkowski addition, but the first may be exponential.

**(c) Fixed-dimensional convexification has bounded extremal distortion.**
Let `r=dim span{v_i}` and arrange the norms in decreasing order
`a_1>=...>=a_n`, where `a_i=||v_i||`.  Then

```math
d_H^{\|\cdot\|}(S(V),Z(V))
\le \sum_{i=1}^{\min(r,n)}a_i.                      \tag{KM.14}
```

Consequently, uniformly for every target `t`,

```math
0\le D_V(t)-\operatorname{dist}(t,Z(V))
\le \sum_{i=1}^{\min(r,n)}a_i.                     \tag{KM.15}
```

If `r` is fixed and `||v_i||<=R`, the error is at most `rR`, independent of
the number of composed vectors.  On an extensive `nR` scale it is `o(nR)`.

#### Proof

Parts (a) and (b) follow directly from the definitions.  For (c), fix
`z in Z(V)`.  There is a coefficient vector `t in [-1,1]^n` with

```math
z=\sum_i t_i v_i.
```

Choose an extreme point of the polytope

```math
P_z=\{u\in[-1,1]^n:\sum_i u_iv_i=z\}.
```

At most `r` coordinates of this extreme point lie strictly between `-1`
and `1`.  Otherwise the corresponding vectors are linearly dependent and a
small signed perturbation in that null direction stays inside the cube,
contradicting extremality.

Keep every integral coordinate and round each of the at most `r` fractional
coordinates to its nearer sign.  The resulting signed sum `s` obeys

```math
\|z-s\|
\le\sum_{i\text{ fractional}}|t_i-\epsilon_i|\|v_i\|
\le\sum_{i=1}^{\min(r,n)}a_i.
```

Since `S(V) subset Z(V)`, this proves (KM.14).  For a closest
`z in Z(V)` to a target, the triangle inequality gives (KM.15). `square`

This is the two-point specialization of the Shapley--Folkman principle, but
the floating-variable proof gives the sharper `sum a_i` constant instead of
a bound by the sum of exceptional diameters.

### Theorem KM.3 (response-body convexification under composition)

The same mechanism is not special to signs.  Let
`E_i subset R^p` be nonempty compact component response sets, put

```math
E=E_1+\cdots+E_n,
\qquad K=\operatorname{conv}E
       =\operatorname{conv}E_1+\cdots+\operatorname{conv}E_n,
```

and let `Delta_i=diam(E_i)` in a fixed norm, arranged decreasingly.  Put

```math
r=\dim\operatorname{span}\bigcup_i(E_i-E_i)\le p.  \tag{KM.16a}
```

Then

```math
d_H(E,K)\le\sum_{i=1}^{\min(r,n)}\Delta_i.           \tag{KM.16}
```

For every `L`-Lipschitz aggregate query `Psi:R^p -> R`,

```math
\left|\sup_{e\in E}\Psi(e)-\sup_{z\in K}\Psi(z)\right|
\le L\sum_{i=1}^{\min(r,n)}\Delta_i,                \tag{KM.17}
```

and the same holds for infima.

#### Proof

Translate one point of every `E_i` to the origin.  All remaining affine
differences lie in the `r`-dimensional space in (KM.16a).  The
Shapley--Folkman lemma therefore represents any `z in K` as a sum in which
all but at most `r` component terms lie in `E_i`, while every exceptional
term lies in `conv(E_i)`.  Replace each exceptional term by any point of its
component set.  The resulting point of `E` differs by at most the sum of the
corresponding diameters, bounded by the right side of (KM.16).  Since
`E subset K`, the Hausdorff estimate follows.  A Lipschitz function changes
by at most `L` times this distance, proving (KM.17). `square`

The lemma was first published in Appendix 2 of R. M. Starr,
[“Quasi-Equilibria in Markets with Non-Convex Preferences,” *Econometrica*
37 (1969), 25--38](https://doi.org/10.2307/1909201).  The response-body
consequence (KM.16)--(KM.17) and its query interpretation are the deductions
used here.

In the language of the project, each `conv(E_i)` is a finite-dimensional
response body and Minkowski addition is its exact composition algebra.
Theorem KM.3 says that after an arbitrary number of additive compositions,
this convex state also answers every Lipschitz nonlinear aggregate query up
to a defect charged to only `r` components.  It is an approximate
synchronization theorem between the discrete reachable set and its convex
roof.

For lifted points `E_i={(phi_i(x),H_i(x))}`, one has `r<=d+1`.  If the
queried objective is of the form `h+g(u)`, only the upper
response boundary is operationally relevant.  Thus fixed-dimensional roofs
can control a declared class of nonlinear composed extrema without storing
the complete product landscape.  This is automatically a succinct state for
zonotopes, generator descriptions, and other controlled response-body
families.  For arbitrary compact component sets, the convex bodies themselves
may have large description complexity, so (KM.16) alone is not a compression
theorem.

### A sharp dimensional obstruction

The fixed-dimensional conclusion cannot be made dimension free.  In
`R^d`, for each coordinate `j` compare the two vector pairs

```math
A_j=(2e_j,2e_j),
\qquad B_j=(3e_j,e_j),
```

and let `A` and `B` be the concatenations over all coordinates.  Both have
the identical response roof

```math
Z(A)=Z(B)=[-4,4]^d.                                 \tag{KM.18}
```

But for the `ell_1` target query at zero,

```math
\min_{s\in S(A)}\|s\|_1=0,
\qquad
\min_{s\in S(B)}\|s\|_1=2d.                       \tag{KM.19}
```

The first signing cancels each equal pair; every coordinate of a signed sum
from `B` has absolute value at least two.  There are `2d` input vectors, so
the gap remains a fixed positive fraction of the extensive scale.  Exact
linear support responses, the full zonotope, and every statistic derived
only from that roof all miss it.

This counterexample cleanly separates two regimes:

- fixed response rank: convexification loses at most `O(dR)=o(nR)`;
- rank proportional to the number of components: convexification can lose a
  leading amount.

The model is standard vector balancing in the sense of, for example,
[Banaszczyk's Gaussian-measure theorem](https://doi.org/10.1002/(SICI)1098-2418(199807)12:4%3C351::AID-RSA3%3E3.0.CO;2-S),
but Theorems KM.2--KM.3 concern information sufficiency rather than the best
universal discrepancy constant.

## 4. Classification against the current framework

### Which standard objects are roofs?

- A trellis boundary kernel is exactly an upper response roof for endpoint
  fields.
- A vector-balancing zonotope is the domain of the zero-energy upper response
  roof and is exactly the support-query quotient.
- A GREM covariance function is **not** a roof of a realized sample.  It is a
  parameter of a Gaussian ensemble law.  Its convex hull is a further
  asymptotic quotient for selected disorder-averaged queries.

### What information is forgotten?

- Gaussian covariance forgets the realization; Gaussianity regenerates its
  law, while a deterministic adversary cannot.
- A trellis kernel forgets internal paths, but the separator screens every
  permitted future context exactly.
- A zonotope forgets the holes in the signed-sum set.  Fixed-dimensional
  additive composition makes those holes uniformly shallow by
  Theorem KM.2; growing dimension does not.

### Does the state contain less than the landscape?

- A homogeneous `k`-level GREM law has `O(k)` variance parameters when the
  hierarchy is fixed, versus exponentially many realized leaf energies.
- A `Q`-state trellis segment has `Q^2` kernel entries, versus exponentially
  many paths.
- A vector list stores `O(nd)` generator coordinates and its zonotope is a
  Minkowski sum of `n` segments, versus `2^n` signed states.  For fixed `d`,
  it answers all target-distance queries within `O(dR)` without enumerating
  the signed-sum set.

## 5. Director checkpoint

### Did this explain something new, or merely rename old machinery?

The ingredients are classical: Gaussian sufficiency, hierarchical extreme
recursions, max-plus cyclicity, and Shapley--Folkman/floating variables.  It
would be wrong to present any of them as a newly discovered field.

The combined response classification is nevertheless generative at a
limited theorem level.  It isolates three mathematically distinct reasons a
small state can preserve a zero-entropy extreme:

1. **regeneration from a law** (REM/GREM);
2. **conditional screening by a separator** (trellises); and
3. **amortized nonconvexity in fixed feature dimension** (vector balancing).

The third mechanism was not explicit in the current surface framework.
Theorem KM.3 says that a convex response body can become approximately
sufficient for nonlinear extremal queries after composition, with a defect
that does not accumulate with the number of factors.  The construction
(KM.18)--(KM.19) supplies its sharp stopping condition.  This is more than a
vocabulary change, although it remains an application of a known theorem.

The max-plus example also answers the feature-growth question for a second
model outside fixed-rank mean field: a finite separator closes exactly, and
under homogeneous repetition its normalized state is eventually periodic.

### Single strongest next theorem

The most promising next statement is a **low effective-dimension
Shapley--Folkman theorem for exposed response faces**:

> Replace ambient feature dimension `p` in (KM.16) by a query-dependent
> dimension of the exposed normal cone or by the rank of the response-
> separation polytope, uniformly over a declared nonlinear query family.

If true, this would let a high-dimensional feature algebra have small
extremal distortion whenever only a low-dimensional exposed face is relevant.
It would apply simultaneously to vector balancing, finite-type CSP count
vectors, and lifted response bodies.  It is falsified if one can construct
bounded-diameter component sets whose queried exposed faces have bounded
local rank but whose discrete-versus-convex optimum gap grows linearly with
the number of components.

This target is genuinely different from asking for another exact roof: it
asks when composition itself regularizes the information that convexification
forgets.
