# Actual-child projective atoms and the flip-information ceiling

Status: **rigorous actual-minimizer theorem plus an abstract scope
ceiling**.  The actual statement below applies to every finite
contracted-temperature pressure minimizer.  The final half-atom construction
is only a law on the correct augmented projective state space satisfying all
of the *inequality signs* in the flip optimality system.  It is **not**
claimed to be the Gibbs law of a signing, let alone of an actual minimizing
child.

The purpose is to test whether the shared-latent common-sign mechanism in
Theorem 37.58 can already be excluded by (AC.32)/(FC.5).  It cannot.  Exact
minimality used only through the one-edge tangents gives a constant atom
bound.  The scalar pressure contraction plus the quadratic Hamming-sphere
identity gives the much stronger exponential spread theorem in
[`actual_child_sector_min_entropy.md`](actual_child_sector_min_entropy.md).
The values of an exponential-size, but only `m`-dimensional, flip subcube
recover the complete augmented child law.  The signs of those inequalities,
even the signs for all edge sets, still allow an atom arbitrarily close to
mass `1/2` in the abstract constraint cone.

## 1. Actual child and its projective atoms

Let `m>=3`, `E=binom([m],2)`, and let `A` minimize

```math
 Z_C(t)=E_{x,tau}\exp\{t\tau H_C(x)\},
 \qquad H_C(x)=\sum_{e=\{i,j\}}c_ex_ix_j                 \tag{PA.1}
```

over all signings `C`, at a fixed `t>0`.  Put `rho=tanh t`.  The augmented
Gibbs law and signed edge-cavity variables are

```math
 \nu_A(x,\tau)
 ={e^{t\tau H_A(x)}\over2^{m+1}Z_A(t)},
 \qquad Y_e=\tau a_ex_ix_j.                            \tag{PA.2}
```

Because the Hamiltonian is invariant under `x -> -x`, pass to

```math
 \Omega_m=\{+-1\}\mathbin\times
          (\{+-1\}^m/\{x\sim-x\}),
 \qquad |\Omega_m|=2^m.                                \tag{PA.3}
```

The projective atom at `z=(tau,[x])` has mass

```math
 \bar\nu_A(z)={e^{t\tau H_A(x)}\over2^mZ_A(t)}.         \tag{PA.4}
```

Exact flip minimality is

```math
 R_A(S):={Z_{A^S}(t)\over Z_A(t)}
 =E_{\bar\nu_A}\exp\{-2t\sum_{e\in S}Y_e\}\ge1
 \qquad(S\subseteq E).                                 \tag{PA.5}
```

Its one-edge tangent is

```math
 E_{\bar\nu_A}Y_e\le\rho\qquad(e\in E).               \tag{PA.6}
```

This is (FC.5)/(FC.15), with the signing absorbed into `Y_e`.

**Theorem PA.1 (actual-minimizer projective-atom ceiling).**  Every exact
pressure-minimizing child satisfies

```math
 \boxed{
 \|\bar\nu_A\|_\infty\le{1+\tanh t\over2}.}             \tag{PA.7}
```

More precisely, if `z=(tau,[x])` and `tau H_A(x)<=0`, then

```math
 \bar\nu_A(z)\le2^{-m}.                                \tag{PA.8}
```

If `tau H_A(x)>0`, then (PA.7) holds for that atom.

*Proof.*  Since `Z_A(t)=E_x cosh(tH_A(x))>=1`, (PA.4) immediately gives
(PA.8).  Otherwise

```math
 \tau H_A(x)=\sum_{e\in E}Y_e(z)>0,
```

so some edge `e` has `Y_e(z)=1`.  Write
`p=bar nu_A(z)`.  Since every other value of `Y_e` is at least `-1`,

```math
 E_{\bar\nu_A}Y_e\ge p-(1-p)=2p-1.
```

Combine this with (PA.6) to obtain `p<=(1+rho)/2`. `square`

At the contracted scale `t=beta/sqrt(N)`, Theorem PA.1 is only

```math
 \|\bar\nu_A\|_\infty\le{1\over2}+O_beta(N^{-1/2}).    \tag{PA.9}
```

Thus the direct optimizer tangent does not make the projective law diffuse.
In particular it does not rule out a fixed positive-mass latent atom.  This
is a scope statement about the tangent, not the final actual-child result:
the quadratic spread theorem quoted below does rule it out.

## 2. The exact low-information test for a two-word bridge atom

For two actual children `A,D` of orders `m,n`, condition their zero-bridge
augmented product law on relative orientation `epsilon=tau_1tau_2`.  Use

```math
 Z_A^s=E_xe^{stH_A(x)},
 \qquad
 \mathcal Z_\epsilon
 =Z_A^+Z_D^\epsilon+Z_A^-Z_D^{-\epsilon}.              \tag{PA.10}
```

Let `mu_epsilon` be the resulting law of the rank-one bridge word
`Q=tau_1XY^T`.  If `Q=xy^T`, direct summation over its four augmented-spin
preimages gives

```math
 \boxed{
 \mu_\epsilon(Q)
 ={2^{2-m-n}\over\mathcal Z_\epsilon}
   \cosh\{t(H_A(x)+\epsilon H_D(y))\}.}                \tag{PA.11}
```

Consequently

```math
 \boxed{
 \eta_\epsilon(A,D;t):=\|\mu_\epsilon\|_\infty
 ={2^{2-m-n}\over\mathcal Z_\epsilon}
   \cosh\!\left(t\max_{x,y}
      |H_A(x)+\epsilon H_D(y)|\right).}                \tag{PA.12}
```

The law is centrally symmetric, so the largest antipodal two-word block has
mass `2 eta_epsilon`.  Hence

```math
 \eta_\epsilon(A,D;t)=o(1)                             \tag{PA.13}
```

is the weakest direct scalar condition excluding a fixed positive-mass
copy of the exact two-word common-sign mechanism of Theorem 37.58 in that
orientation.  Formula (PA.12) is genuinely much smaller than the complete
bridge landscape: it uses the two sector partition functions and the four
signed child energy extrema.  It does **not** exclude a diffuse cluster of
many rank-one words producing the same row-retuning effect, and a lower
bound on `eta_epsilon` does not by itself prove linear retuning.

*Proof of (PA.11).*  Given `Q`, for each `s=tau_1` there are exactly two
spin pairs `(x,y),(-x,-y)` satisfying `sxy^T=Q`.  In the conditional law
from (LE.2), cancellation of the sector partition functions leaves

```math
 {2\over2^{m+n}\mathcal Z_\epsilon}
 \sum_{s=+-1}e^{st(H_A(x)+\epsilon H_D(y))},
```

which is (PA.11).  Maximization proves (PA.12). `square`

This identifies a strict, optimizer-specific observable for the literal
two-word obstruction.

There is a robust version which remains far smaller than the full child
table.  Let

```math
 \bar\mu_{A,s}([x])
 ={e^{stH_A(x)}\over2^{m-1}Z_A^s},
 \qquad
 \kappa_{A,s}(\delta)
 =\max_{u\in\{+-1\}^m}\bar\mu_{A,s}
   \left\{[x]:{|\langle x,u\rangle|\over m}
                    \ge1-\delta\right\}.              \tag{PA.13a}
```

Define the largest antipodal rank-one cap of the actual orientation law by

```math
 \Xi_\epsilon(\delta)
 =\max_{u,v}\mu_\epsilon\left\{Q:
   {|\langle Q,uv^{\mathsf T}\rangle|\over mn}
                         \ge1-\delta\right\}.           \tag{PA.13b}
```

With the exact sector weights `pi_s^(epsilon)` from (LE.2), one has the
two-sided deterministic comparison

```math
 \boxed{
 \max_s\pi_s^{(\epsilon)}
   \kappa_{A,s}(\delta/2)\kappa_{D,\epsilon s}(\delta/2)
 \le\Xi_\epsilon(\delta)
 \le\sum_s\pi_s^{(\epsilon)}
   \kappa_{A,s}(\delta)\kappa_{D,\epsilon s}(\delta).}  \tag{PA.13c}
```

Indeed, for `Q=sxy^T`,

```math
 {|\langle Q,uv^T\rangle|\over mn}
 ={|\langle x,u\rangle|\over m}
  {|\langle y,v\rangle|\over n}.                       \tag{PA.13d}
```

If the left side is at least `1-delta`, both factors are at least
`1-delta`, proving the upper bound after conditioning on `s`.  Conversely,
if both are at least `1-delta/2`, their product is at least
`(1-delta/2)^2>=1-delta`; choose the maximizing caps in one sector to prove
the lower bound.

Thus the four sector-resolved scalar functions `kappa_(A,s),kappa_(D,s)`
decide whether the actual latent law has a fixed positive-mass *geometric
common-sign cluster*.  They are genuinely coarser than the complete Gibbs
table.  They do not control coherent retuning carried by a diffuse family
with `Xi_epsilon(delta)=o(1)`, so they are not by themselves a closure of
`L_row-lifetime-closure`.

The new quadratic spread theorem makes this criterion effective.  Put

```math
 C_\beta=\log2+{\beta^2\over4},
 \qquad
 \underline\eta_\beta={e^{-1-\beta^2}\over16}.        \tag{PA.13e}
```

Equations (ME.2)--(ME.9) prove, uniformly over actual children and sectors,

```math
 \|\mu_{A,s}\|_\infty
 \le\exp\{-(\underline\eta_\beta-o(1))m\}.             \tag{PA.13f}
```

The input is genuinely coarser than the Gibbs table: the scalar optimizer
bound `Z_A(t)<=(cosh t)^{binom m2}` bounds `t max|H_A|/m`, and exact
averaging over a Hamming sphere of a quadratic form supplies exponentially
many comparison states.  No nonradial flip values are used.

Since an overlap cap in (PA.13a) contains at most
`2 exp{m h(delta/2)+o(m)}` full spin states, (PA.13f) gives

```math
 \boxed{
 \kappa_{A,s}(\delta)
 \le\exp\{-[\underline\eta_\beta-h(\delta/2)-o(1)]m\}.} \tag{PA.13g}
```

Together with (PA.13c), this proves exponential decay of
`Xi_epsilon(delta)` for every fixed `delta>0` small enough that
`h(delta/2)<underline eta_beta`.  It also proves

```math
 \eta_\epsilon(A,D;t)
 \le\exp\{-(\underline\eta_\beta-o(1))(m+n)\}          \tag{PA.13h}
```

at comparable contracted-temperature splits.  More generally, a union of
`exp(kappa N)` such rank-one caps has vanishing prior mass whenever

```math
 \kappa+h(\delta/2)<\underline\eta_\beta.              \tag{PA.13i}
```

If the negative-disorder posterior assigns fixed mass to such a union,
binary data processing forces a linear retuning KL.  Hence the exact
two-word obstruction and every sufficiently low-rate narrow-cluster variant
are incompatible with **sublinear** actual posterior retuning.  What remains
is diffuse retuning across an exponential-rate family; the spread theorem
does not bound its row-lifetime contribution.

## 3. A minimal edge set encodes the whole augmented law

The complete augmented child law can nevertheless be recovered from a much
smaller flip cube than the `2^{binom m2}` full signing cube.  Let

```math
 \mathcal B
 =\{\{1,i\}:2\le i\le m\}\mathbin\cup\{\{2,3\}\},
 \qquad |\mathcal B|=m.                               \tag{PA.14}
```

**Theorem PA.2 (cycle-basis flip inversion).**  The map

```math
 \Gamma_A:\Omega_m\longrightarrow\{+-1\}^{\mathcal B},
 \qquad z\longmapsto(Y_e(z))_{e\in\mathcal B}          \tag{PA.15}
```

is a bijection.  For every law `pi` on `Omega_m`, its `2^m` transform values

```math
 R_\pi(S)=E_\pi\exp\{-2t\sum_{e\in S}Y_e\},
 \qquad S\subseteq\mathcal B,                          \tag{PA.16}
```

determine `pi` exactly.  In particular, for an actual child the pressure
ratios `R_A(S)` on this `m`-edge flip subcube determine its complete
augmented projective Gibbs law.

*Proof.*  Fix the projective representative `x_1=1`.  Write
`z_i=x_i` for `i>=2`.  The star variables are

```math
 Y_{1i}=\tau a_{1i}z_i.
```

The one extra triangle edge recovers the sector bit:

```math
 Y_{12}Y_{13}Y_{23}=\tau a_{12}a_{13}a_{23}.           \tag{PA.17}
```

Thus `tau` and then every `z_i` are recovered from (PA.15).  Domain and
codomain both have size `2^m`, proving bijectivity.

Put `c=cosh(2t)`, `s=sinh(2t)`, and define the Walsh moments

```math
 M(T)=E_\pi\prod_{e\in T}Y_e.
```

Since `e^{-2tY_e}=c-sY_e`,

```math
 R_\pi(S)
 =\sum_{T\subseteq S}c^{|S|-|T|}(-s)^{|T|}M(T).       \tag{PA.18}
```

This subset-triangular transform has nonzero diagonal because `t>0`, so it
recovers every `M(T)`.  Finally, for `y in {+-1}^mathcal B`,

```math
 \pi\{Y=y\}
 =2^{-m}\sum_{T\subseteq\mathcal B}
     \left(\prod_{e\in T}y_e\right)M(T),              \tag{PA.19}
```

which is ordinary Walsh inversion. `square`

For arbitrary laws on `Omega_m`, exact recovery requires `2^m-1`
independent real parameters: the probability simplex has that dimension,
and (PA.18) is an invertible change of coordinates.  This dimension
statement is about the unrestricted law class.  It is not a lower bound on
the dimension of the special family of Gibbs laws coming from minimizing
signings.

## 4. The inequality signs alone allow a half atom

The preceding inversion uses the **values** of the flip ratios.  The signs
`R_A(S)>=1` are far weaker, even if retained for every edge set rather than
only the basis `mathcal B`.

**Proposition PA.3 (abstract flip-sign cone has fixed atoms).**  Fix any
`z_0=(tau_0,[x_0]) in Omega_m`, and let `bar z_0=(-tau_0,[x_0])`; then
`Y(bar z_0)=-Y(z_0)` on every edge.  For `0<epsilon<1`, define the law

```math
 \pi_\epsilon
 =(1-\epsilon){\delta_{z_0}+\delta_{\bar z_0}\over2}
   +\epsilon U_{\Omega_m}.                            \tag{PA.20}
```

For every `S subseteq E`,

```math
 \boxed{
 E_{\pi_\epsilon}e^{-2t\sum_{e\in S}Y_e}\ge1,}       \tag{PA.21}
```

while

```math
 \|\pi_\epsilon\|_\infty
 ={1-\epsilon\over2}+{\epsilon\over2^m}.             \tag{PA.22}
```

For every nonempty `S`, (PA.21) is strict.

*Proof.*  The antipodal part of (PA.20) contributes

```math
 \cosh\left(2t\sum_{e\in S}Y_e(z_0)\right)\ge1.
```

For the uniform part, averaging first over `tau` gives

```math
 E_{[x]}\cosh\left(2t\sum_{e\in S}a_ex_ix_j\right)
 \ge1.                                                \tag{PA.23}
```

If nonempty `S` made the polynomial inside (PA.23) zero at every cube
point, all of its distinct Walsh coefficients `a_e` would vanish, which is
impossible.  Hence the uniform contribution is strictly greater than one
for nonempty `S`.  Convexity proves (PA.21), and (PA.22) is immediate.
`square`

Proposition PA.3 is deliberately an **abstract sign-cone ceiling**.  The law
`pi_epsilon` need not have the Gibbs form (PA.2), and no actual minimizing
signing realizing it is asserted.  What it proves is the following precise
logical limitation:

```math
 \boxed{
 \text{the inequality directions in (AC.32)/(FC.5), by themselves, do not
 imply projective min-entropy }\longrightarrow\infty.} \tag{PA.24}
```

Any proof excluding common-sign coherent retuning must therefore use at
least one of:

1. quantitative values of nonradial flip responses (the basis cube in
   Theorem PA.2 is sufficient but exponentially large);
2. quadratic/Gibbs rigidity beyond the inequality signs, as in the sector
   spread theorem; or
3. a new theorem controlling diffuse exponential-rate coherent clusters.

This **strictly narrows but does not close** the row-lifetime SML.  Literal
finite-word and low-rate narrow-cluster retuning are now excluded at
sublinear posterior KL cost.  The surviving question is whether the actual
inverse-disorder posterior can retune diffusely over an exponential-rate
family and still create linear `J-I^leftarrow`, or whether child minimality
compresses that diffuse retuning to a controllable direction.  The note also
shows exactly why AC.32 used only as a family of lower inequalities cannot
answer that question.
