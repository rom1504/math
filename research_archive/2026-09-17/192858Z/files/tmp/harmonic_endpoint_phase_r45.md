# Wave 45B: an affine-phase endpoint theorem and its sharp degree wall

## Status

This note proves a **fixed finite-state low-temperature theorem** for the
integrated harmonic load.  It is not yet the project-uniform endpoint bound.
It identifies the missing geometric quantity precisely: endpoint domination
holds with the weighted overlap of simultaneous base-surprise crossings (and,
in particular, with their maximum active-state incidence degree), divided by
the surviving ground-phase mass.

The theorem passes the exact `A8,m=4` phase wall.  An abstract tropical star
then shows that an affine interpolation envelope and a ground-phase mass
bounded away from zero do **not** suffice by themselves: the ratio can grow
linearly with the crossing degree.  Thus a project proof must extract a
uniform weighted-overlap estimate from minimizer/completion structure, rather
than from generic posterior-mixture geometry.

Checks are in
`/home/math/quadra/tmp/harmonic_endpoint_phase_r45_check.py`.

## 1. Exact hypotheses

Let `Omega` be a fixed finite state space equipped with a fixed undirected
coordinate graph.  Consider a low-temperature family

```math
\nu_\beta(d)=Z_{0,\beta}^{-1}\exp\{\beta E_d+o(1)\},
\qquad
\mu_{1,\beta}(d)=Z_{1,\beta}^{-1}
c_d\exp\{\beta(E_d-a_d)\}(1+o(1)),
\tag{R45B.1}
```

where `E_d,a_d` are fixed real numbers, `c_d>0`, and all remainders are
uniform in `d`.  Interpolate exponentially,

```math
\mu_{s,\beta}(d)\propto
\nu_\beta(d)^{1-s}\mu_{1,\beta}(d)^s.
\tag{R45B.2}
```

For every coordinate edge `e={x,y}`, let `k_e(s)` be the conditional binary
KL from the `mu_s` conditional law to the `nu` conditional law.  It has the
same value at both endpoints.  Put

```math
L_s(d)=\sum_{e\ni d}k_e(s),
\qquad
\mathcal J_L^2
=\int_0^1\frac{\operatorname{Var}_{\mu_s}(L_s)}s\,ds,
\qquad
C_V(1)=\mathbb E_{\mu_1}L_1.
\tag{R45B.3}
```

The sum may be restricted to the vertex coordinates, as in the project.
Assume the tropical interpolation envelope is affine:

```math
\Lambda(s)=\max_d(E_d-sa_d)
=(1-s)\max_d E_d+s\max_d(E_d-a_d),
\qquad 0\le s\le1.
\tag{R45B.4}
```

For a finite collection of lines, (R45B.4) is equivalent to the existence of
a state maximizing at both endpoints.  Write

```math
\mathcal A=\operatorname*{argmax}_d(E_d-a_d),
\qquad E_*=\max_dE_d,
\qquad
\mathcal G=\{d\in\mathcal A:E_d=E_*\}.
\tag{R45B.5}
```

Thus `G` is nonempty.  Define

```math
Z(u)=\sum_{d\in\mathcal A}c_d e^{-u(E_*-E_d)},
\qquad
\pi_u(d)=\frac{c_de^{-u(E_*-E_d)}}{Z(u)},
\qquad
\eta=\frac{Z(\infty)}{Z(0)}
=\frac{\sum_{d\in\mathcal G}c_d}
{\sum_{d\in\mathcal A}c_d}.
\tag{R45B.6}
```

The positive leading edges incident to an active state are of exactly two
types.

1. If both endpoints are active and `E_h-E_l=h>0`, with leading coefficients
   `A,B` at the high and low endpoints, then

   ```math
   \kappa_e(u)=h\frac{Be^{-hu}}{A+Be^{-hu}}.
   \tag{R45B.7}
   ```

2. If only `x` is active and its inactive neighbor has base energy
   `E_x+h`, then

   ```math
   \kappa_e(u)=h.
   \tag{R45B.8}
   ```

All other edges have zero leading load.  Put

```math
\ell_u(d)=\sum_{e\ni d}\kappa_e(u),
\qquad
\Delta=\max_{d\in\mathcal A}
\#\{e\ni d:\kappa_e\not\equiv0\}.
\tag{R45B.9}
```

## 2. The affine-phase theorem

**Theorem (verified).** Under (R45B.1)--(R45B.4), for fixed finite `Omega`,

```math
\frac{\mathcal J_L^2}{\beta}
\longrightarrow
I:=\int_0^\infty
\operatorname{Var}_{\pi_u}(\ell_u)\,du,
\qquad
\frac{C_V(1)}\beta
\longrightarrow\mathbb E_{\pi_0}\ell_0,
\tag{R45B.10}
```

and

```math
\boxed{
I\le\frac{\Delta}{\eta}\,
\mathbb E_{\pi_0}\ell_0.}
\tag{R45B.11}
```

Consequently

```math
\boxed{
\mathcal J_L^2
\le\left(\frac{\Delta}{\eta}+o(1)\right)C_V(1).}
\tag{R45B.12}
```

The finite-state statement also admits the sharper replacement of `Delta`
by the least boundary-layer overlap constant `D` satisfying

```math
\mathbb E_{\pi_u}
\left(\sum_{e\ni d}\kappa_e(u)\right)^2
\le D\sum_e\sum_{d\in e\cap\mathcal A}
\pi_u(d)\kappa_e(u)^2
\quad(u\ge0).
\tag{R45B.13}
```

The elementary degree estimate gives `D<=Delta`.  A project theorem only
needs a uniform weighted `D/eta`, not a uniform pointwise degree.

### Proof

Put `u=beta(1-s)`.  Uniform finite-state Laplace asymptotics give

```math
\mu_{1-u/\beta,\beta}\longrightarrow\pi_u,
\qquad
\beta^{-1}L_{1-u/\beta}\longrightarrow\ell_u.
\tag{R45B.14}
```

Formulae (R45B.7) and (R45B.8) follow by taking the binary conditional KL.
For example, on an active-active edge the base-disfavored endpoint has
conditional probability `Be^(-hu)/(A+Be^(-hu))`, and its log likelihood
penalty is `beta h+O(1)`.  Affineness supplies a common base/endpoint ground
state, so all non-endpoint layers are exponentially localized and contribute
`o(beta)` to `J_L^2`; the `s=0` boundary layer contributes only `O(1)` because
`L_0=0`.  Dominated finite-state Laplace asymptotics give (R45B.10).

Use `Var(ell)<=E ell^2` and Cauchy--Schwarz over incident positive edges:

```math
\ell_u(d)^2
\le\Delta\sum_{e\ni d}\kappa_e(u)^2.
\tag{R45B.15}
```

It remains to integrate one edge.  On an active-active edge, let the higher
base energy be `E_*-r`, and use `A,B,h` as above.  Its contribution after
(R45B.15) is

```math
\int_0^\infty
\frac{h^2B^2e^{-(r+2h)u}}
{Z(u)(A+Be^{-hu})}\,du
\le\frac{hB}{Z(\infty)}.
\tag{R45B.16}
```

Indeed `Z(u)>=Z(infinity)`, discard `e^(-ru)<=1`, and substitute
`x=e^(-hu)`; the remaining integral is at most
`hB^2/Z(infinity) int_0^1 x/(A+Bx) dx <= hB/Z(infinity)`.
The same edge contributes exactly `hB/Z(0)` to
`E_(pi_0) ell_0`.

For an active-inactive reversal edge, let the active coefficient be `c` and
its base deficit be `r=E_*-E_x`.  Since its neighbor has energy `E_x+h` and
`E_*` is the global base maximum, `r>=h`.  Hence

```math
\int_0^\infty\frac{h^2ce^{-ru}}{Z(u)}\,du
\le\frac{hc}{Z(\infty)},
\tag{R45B.17}
```

whereas its endpoint contribution is `hc/Z(0)`.  Summing (R45B.16) and
(R45B.17) and using `Z(0)/Z(infinity)=1/eta` proves (R45B.11).

## 3. The exact `A8,m=4` wall passes

The Wave-44 exact tropical enumeration has twenty endpoint-active states.
Eight of their ten units of leading coefficient have maximal base energy, so

```math
\eta=\frac8{10}=\frac45.
```

Its seven active crossing edges are disjoint and every active state is
incident to at most one, hence `Delta=1`.  The theorem gives the valid (very
loose) coefficient `Delta/eta=5/4`.  The exact boundary calculation is much
smaller:

```math
\frac{\mathcal J_L^2}{C_V(1)}
\longrightarrow0.0133226596003\ldots .
```

The exact project examples `A6`, `A8`, and `A9` all satisfy the affine-envelope
test for every audited selector size.  Exhaustive signing enumeration gave the
same result for all exact minimizers through `n=6` (counts `48,192,384` for
`n=4,5,6`, respectively).  This is finite evidence, not a theorem.  A stored
non-minimizing `n=6,m=4` signing has no common endpoint maximizer and a genuine
interior envelope facet, so affineness is not generic signing geometry.

## 4. Sharp obstruction: an affine tropical star

For each `k`, take an abstract phase graph with a high-energy center, `k`
low-energy leaves joined only to the center, and one isolated high-energy
state.  Give the center and leaves coefficient one and the isolated state
coefficient `k+1`.  Let high energy be `h`, low energy zero, and set `a_d=E_d`
for every state.  Thus every endpoint line is tied, the envelope is affine,
and

```math
\eta=\frac{k+2}{2k+2}\longrightarrow\frac12,
\qquad \Delta=k.
\tag{R45B.18}
```

With `x=e^(-hu)`,

```math
Z(u)=k+2+kx,
\qquad p(u)=\frac{x}{1+x}.
\tag{R45B.19}
```

The center load is `khp(u)`, each leaf load is `hp(u)`, and the isolated load
is zero.  Direct integration gives

```math
\frac{C_V(1)}\beta\longrightarrow\frac{kh}{2(k+1)}
\longrightarrow\frac h2,
\qquad
I=\frac{kh}{8}+O(h),
\tag{R45B.20}
```

and therefore

```math
\boxed{
\frac{\mathcal J_L^2}{C_V(1)}
\sim\frac k4.}
\tag{R45B.21}
```

So affine phase geometry plus `eta>=1/2+o(1)` does not imply endpoint
domination.  Some crossing-overlap control is necessary, and the linear
dependence on `Delta` in (R45B.11) has the correct order for generic phase
data.  This is an abstract tropical graph, not an asserted realization by an
exact minimizer/completion posterior.

## 5. Exact remaining target

The phase theorem narrows, but does not prove, the project endpoint estimate.
A sufficient project result is:

```math
\text{(i) all interior tropical facets have total integrated load charge}
\le C C_V(1),
```

```math
\text{(ii) on every endpoint boundary layer, }D/\eta\le C,
\tag{R45B.22}
```

uniformly in the exact minimizer, selector size, and target temperature.
Condition (ii) may be attacked through the completion posterior: it asks for
an `L^2` bound on simultaneous base-surprise incidences under the actual
phase weights, which is strictly weaker than bounded maximum degree.  A
counterexample with unbounded `D/eta` inside the exact minimizer/completion
family would sharply falsify this phase route.

Even if (R45B.22) is proved, restoring, orientation control, and
adjacent-selector Hellinger remain separate obligations.  Local Poincare,
Efron--Stein, monotone deletion, and the raw-reveal baseline remain retired by
the earlier walls.
