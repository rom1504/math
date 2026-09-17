# Wave 36 harmonic attack: exclusion identities and a scalable migration wall

## Status and scope

- **Verified exact identity:** for every vertex coordinate, endpoint erasure
  KL is an exclusion-posterior logarithmic Jensen gap,
  `C_i=E_mu log(bar r_i(D_-i)/r_i(D))`; the global orientation is not part of
  this identity and remains separately charged.
- **Verified sufficient lemma:** for any finite endpoint likelihood
  `g=log f`, weighted context backtracking is at most
  `(exp(osc(g))-1) sum_j C_j`.  This closes migration only if signing
  structure gives a uniform `O(1)` likelihood oscillation, which is not
  known.
- **Scoped scalable falsifier:** fixed selector size, exact vertex omission,
  the cancellation `J_j=C_j+I(S;D_j|D_-j)`, full selector support, and a
  separate nonzero orientation coordinate do **not** imply either bounded
  backtracking or endpoint domination of actual adverse migration.  An
  explicit canonical mixture has `sum C_j=O(1/M)`, weighted backtracking
  `Theta(1)`, and adverse covariance `Omega(1/log^2 M)`.  Tensor powers give
  fixed selector density and separate the project scales.
- **Scope warning:** the falsifier is not asserted to be the matched endpoint
  of an exact quadratic-signing minimizer.  It rules out generic proofs from
  the named information/omission properties.  The minimizer-specific route
  survives, but must add a no-transient-mode-transfer hypothesis.

The finite calculations and all chain-rule checks are reproduced by
`tmp/harmonic_erasure_r36_check.py`.  Its migration integrals use adaptive
quadrature; the formulas and asymptotic bounds below are analytic.

## 1. Exact posterior-exclusion form of endpoint erasure

Use the canonical endpoint lift

```math
P(S,d)=q(S)\mu_S(d)=\mu(d)\pi_d(S),
\qquad \mu=\nu f,
```

and fix a projective vertex coordinate `i`.  Set

```math
Z_i=\mathbf 1_{\{i\notin S\}},\qquad
r_i(d)=P(Z_i=1\mid D=d),\qquad
\bar r_i(e)=P(Z_i=1\mid D_{-i}=e).
```

If `i notin S`, the component likelihood `h_S` is invariant under the
`i`-flip.  Consequently, context by context,

```math
P(D_i=b\mid D_{-i}=e,Z_i=1)=\nu(D_i=b\mid D_{-i}=e).
\tag{H36.1}
```

Bayes' rule applied to (H36.1) gives the pointwise identity

```math
\boxed{
\mu(D_i=b\mid e)r_i(e,b)
=\bar r_i(e)\nu(D_i=b\mid e).
}
\tag{H36.2}
```

Thus, with no inequality and with the endpoint context law,

```math
\boxed{
C_i
=\mathbb E_\mu\log\frac{\bar r_i(D_{-i})}{r_i(D)}.
}
\tag{H36.3}
```

The same calculation recovers the selector cancellation.  Conditional on a
context, the omitted component is exactly the base bit law, and applying the
KL chain rule before and after revealing the full selector gives

```math
J_i=C_i+I(S;D_i\mid D_{-i}).
\tag{H36.4}
```

For selectors of deterministic size `m`, one also has pointwise

```math
\boxed{\sum_i r_i(d)=n-m.}
\tag{H36.5}
```

Equations (H36.3)--(H36.5) expose the remaining gap.  Fixed size controls a
first moment of the exclusions, whereas endpoint cost contains an unweighted
negative logarithm of a possibly tiny exclusion posterior.  More
importantly, all three formulas are endpoint statements: they do not control
which edge contexts were heavy at intermediate `t` and then vanished before
`t=1`.  The global-orientation cost `C_0` has no `r_i` representation and
must still be added separately.

## 2. A general sufficient anti-evanescence lemma

Let `G=osc(g)=max g-min g`, let `mu_t` be proportional to `nu exp(tg)`, and
let `M_e(t)` be any coordinate-context mass.  Direct differentiation and
backward change of measure give

```math
|M'_e(t)|
\le G M_e(t),
\qquad
\frac{\mu_t(d)}{\mu_1(d)}
=\frac{e^{-(1-t)g(d)}}{\mathbb E_{\mu_1}e^{-(1-t)g}}
\le e^{(1-t)G}.
\tag{H36.6}
```

It follows that

```math
\boxed{
\int_0^1[-M'_e(t)]_+dt
\le(e^G-1)M_e(1).
}
\tag{H36.7}
```

Multiplying by the fixed endpoint binary costs and summing proves

```math
\boxed{
\sum_{j,e}K_e\int_0^1[-M'_e(t)]_+dt
\le(e^{\operatorname{osc}(\log f)}-1)\sum_jC_j.
}
\tag{H36.8}
```

This is a rigorous sufficient version of (10.1019), with an explicit
constant.  It is useful at the project exponent only under a uniform
`osc(log f)=O(1)` estimate (or a sharper cost-weighted replacement).
No such exact-minimizer estimate is presently known.  The construction below
has `osc(log f)=Theta(log M)` and shows that dependence on an
anti-evanescence parameter cannot simply be removed.

## 3. One-block canonical counterexample

Take a separate global-orientation bit `O` and vertex bits `u,v,y`.  Under
the base law, `O,u,y` are fair and independent, while `v=u` with probability
`1-delta` and `v\ne u` with probability `delta>0`.  Hence the base has full
support and `v` is a nearly redundant copy of `u`.

For `M>2` and `eta=M^{-2}`, put

```math
R_\eta(O)=1-\eta\quad(O=0),\qquad
R_\eta(O)=1+\eta\quad(O=1),
```

and define the normalized endpoint likelihood

```math
f_M(O,u,y)=R_\eta(O)\frac{2}{M+2}
\begin{cases}
1,&(u,y)=(0,0),\\
2,&(u,y)=(0,1),\\
M,&(u,y)=(1,0),\\
M+1,&(u,y)=(1,1).
\end{cases}
\tag{H36.9}
```

This is exactly a fixed-size singleton selector mixture on `{u,v,y}`.  Give
the three selectors weights

```math
q_u=\frac{M}{M+2},\qquad
q_v=\frac{1}{2(M+2)},\qquad
q_y=\frac{3}{2(M+2)},
\tag{H36.10}
```

and normalized component likelihoods

```math
\begin{aligned}
h_u&=R_\eta(O)(M^{-1},\,2-M^{-1})_u,\\
h_v&=R_\eta(O),\\
h_y&=R_\eta(O)(1/3,\,5/3)_y.
\end{aligned}
\tag{H36.11}
```

Then `f_M=sum_s q_s h_s`, every `h_s` is invariant under all omitted vertex
coordinates, all weights and likelihoods are strictly positive, and
`sum_i r_i=2` at every state.  The orientation factor is common to every
component but is nonconstant, so its KL is retained rather than quotiented
out.  This is also canonical in the algebraic sense of (10.938): choose a
small common `a>0`, set `exp(-F_s)=a q_s h_s`, and then the normalizations
recover exactly `q_s,h_s,f_M`, with every `F_s>=0`.

Let

```math
\kappa(a)=D(\operatorname{Ber}(e^a/(1+e^a))
\Vert\operatorname{Ber}(1/2)),
\qquad k_0=\kappa(\log2)>0.
```

For the `y` coordinate, the endpoint low-context mass and the two binary
costs are

```math
L_M(1)=\frac{3}{2M+4},\qquad
K_{\rm low}=k_0,\qquad
K_{\rm high}=\kappa\!\left(\log(1+M^{-1})\right).
```

Therefore exactly

```math
\boxed{
C_y=\frac{3k_0}{2M+4}
+\left(1-\frac{3}{2M+4}\right)
\kappa\!\left(\log(1+M^{-1})\right)
=O(M^{-1}).
}
\tag{H36.12}
```

Here `K_high<=1/(2M+1)^2` follows from `D<=chi^2`.  Also `C_v=0`,
`C_O<=eta^2=M^{-4}`, and, by continuity as `delta downarrow0`, one may
choose a strictly positive `delta=delta_M` so that `C_u<=M^{-2}`.  Hence

```math
\sum_jC_j=O(M^{-1}).
\tag{H36.13}
```

The cancellation and omission are nonvacuous.  Only the `y` selector
contributes to the lifted `y` cost, so

```math
J_y=\frac{3}{2(M+2)}
D(\operatorname{Ber}(5/6)\Vert\operatorname{Ber}(1/2)),
\qquad I(S;D_y\mid D_{-y})=J_y-C_y.
\tag{H36.14}
```

Thus the exact subtractive information is present, but it remains an
endpoint cancellation.

Along interpolation, the total low-`u` context mass for the `y` edges is

```math
L_M(t)=\frac{1+2^t}{1+2^t+M^t+(M+1)^t}.
\tag{H36.15}
```

It decreases strictly for `M>2`.  Since all low contexts have endpoint cost
`k_0`, their total drop proves

```math
\boxed{
\sum_{j,e}K_e\int_0^1[-M'_e(t)]_+dt
\ge k_0\left(\frac12-\frac{3}{2M+4}\right)=\Theta(1).
}
\tag{H36.16}
```

Together, (H36.13)--(H36.16) give a backtracking/endpoint ratio
`Omega(M)`; hence no uniform (10.1019) follows from the exact selector
identities, fixed size, and omission.

The failure is not merely an artifact of replacing `k_e(t)` by `K_e`.
Writing `k_low(t)=kappa(t log2)` and
`k_high(t)=kappa(t log(1+M^{-1}))`, the exact `y` migration covariance is

```math
\operatorname{Cov}_y(t)
=L'_M(t)\{k_{\rm low}(t)-k_{\rm high}(t)\}<0.
\tag{H36.17}
```

On `t in [1/log M,2/log M]`, the drop in `L_M` is bounded below by a
positive absolute constant, while
`k_low-k_high=Omega(1/log^2 M)`.  Consequently

```math
\boxed{
\mathcal A_y
=\int_0^1[-\operatorname{Cov}_y(t)]_+dt
=\Omega(1/\log^2 M),
\qquad
\frac{\sum_j\mathcal A_j}{\sum_jC_j}
=\Omega\!\left(\frac{M}{\log^2M}\right).
}
\tag{H36.18}
```

The constants can be made explicit.  For `M>=12`, the two endpoint values
in (H36.15) obey

```math
L_M(1/\log M)\ge\frac14,
\qquad L_M(2/\log M)\le\frac15.
```

For `0<=s<=1`, Pinsker, `tanh(s/2)>=s/4`, and the chi-square upper bound
give

```math
\frac{s^2}{32}\le\kappa(s)\le\frac{s^2}{4}.
```

Hence, after increasing the fixed threshold on `M` if necessary,

```math
\mathcal A_y\ge
\frac{(\log2)^2}{1280(\log M)^2}.
\tag{H36.19}
```

On the other hand, choosing `delta_M` as above gives
`sum_j C_j<=K_C/M` for the absolute constant
`K_C=3k_0/2+2`.  Thus (H36.18) holds with an explicit positive absolute
constant.  Similarly, for `M>=4`, (H36.16) is at least `k_0/4`.

For example, the checker obtains at `M=1000` a weighted-backtracking ratio
`332.51096` and an actual-adverse-migration ratio `25.43500`.

## 4. Fixed-density tensor scaling

Take `N` independent copies of the three vertex bits and one genuinely
separate global orientation bit.  Expand the product likelihood as a
selector mixture by choosing exactly one of `{u_b,v_b,y_b}` in each block.
Every supported selector therefore has size `N` among `3N` vertices,
the selector density is `1/3`, every component is invariant outside its
selector, and the pointwise exclusion sum is `2N`.  Multiplying every
component by the one orientation likelihood retains one, rather than `N`,
global-orientation coordinate.

The selector support can be made all of the fixed-size slice: give every
missing size-`N` selector a positive constant component of sufficiently
small total weight.  Finite-space continuity preserves all inequalities
within a factor two.  Thus zero selector weights are not essential to the
wall.

Before that harmless perturbation, tensor factorization gives exactly

```math
\sum_jC_j=O(N/M),\qquad
\sum_j\mathcal A_j=\Omega(N/\log^2M),\qquad
\sum_{j,e}K_e\int[-M'_e]_+=\Omega(N).
\tag{H36.20}
```

Fix `0<c<1/4`, take `M=N^p` with `p>1/2+2c`, and note that the ambient
dimension is `Theta(N)`.  Then the endpoint-erasure half has the desired
abstract scale,

```math
\sum_jC_j=o(N^{1/2-2c}),
```

while actual adverse migration is

```math
\sum_j\mathcal A_j
=\Omega(N/\log^2N)
=\omega(N^{1/2-2c}).
\tag{H36.21}
```

The positive marginal-information variation in (10.1018) also fails: its
sum is at least its net endpoint increase, and product redundancy makes it
macroscopic (indeed quadratic after summing erased coordinates).  Thus the
two clauses of (10.1018) are genuinely independent at the identity level.

## Frontier

The surviving harmonic target must be explicitly minimizer-specific.  It is
not enough to prove endpoint erasure from selector cancellation.  One must
also prove, uniformly at the prescribed signing temperature, either

```math
\sum_j\mathcal A_j=O(n^{1/2-2c})
```

directly, or a cost-weighted **no-transient-mode-transfer** estimate such as
(10.1019).  A uniform `osc(log f)=O(1)` would suffice by (H36.8), but is much
stronger than current information and is not proposed as the leading
lemma.  A sharper viable theorem may allow large likelihood oscillation but
must prevent a high-`K_e` context from carrying appreciable intermediate
mass and becoming endpoint-rare through a mode shift redundantly visible in
the other coordinates.  Exact quadratic minimality (and the prescribed
matched external-partition form of `F_S`) is the remaining structure capable
of doing this.  Restoring bounds, adjacent-selector Hellinger control, and
the separately charged orientation term remain independent inputs.
