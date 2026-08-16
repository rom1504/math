# Directed responses under an interacting min-plus continuation

**Status.** Working theorem draft awaiting an independent audit.  The general
exposure inequality and the binary Ising calculation are elementary exact
claims.  Their role is to test whether the directed table of Theorem 15.1 has
an algebra beyond direct products.

## 1. The witness that a continuation must preserve

Let `Y,X` be finite, let `K:X times Y -> R`, and define the common min-plus
continuation

```math
(T_Kf)(x)=\min_{y\in Y}\{K(x,y)+f(y)\}.             \tag{DC.1}
```

For functions on a common set write

```math
r(f,g)=\max_y\{f(y)-g(y)\}.                         \tag{DC.2}
```

This directed quantity may be negative and is not itself a metric.  The
two-sided projective separation is controlled only after both orientations
are retained.

For `eta>=0`, define the inputs exposed by `f` through the continuation:

```math
E_{K,\eta}(f)=\left\{y:\text{ for some }x,
 K(x,y)+f(y)\le (T_Kf)(x)+\eta\right\}.             \tag{DC.3}
```

Equivalently, the least exposure penalty of an input is

```math
e_K^f(y)=\min_x\{K(x,y)+f(y)-(T_Kf)(x)\}\ge0,       \tag{DC.3a}
```

and `E_(K,eta)(f)={y:e_K^f(y)<=eta}`.

### Theorem DC.1 (exposed-witness preservation)

For every `f,g` and `eta>=0`,

```math
r(f,g)-\Delta_{K,\eta}(f,g)-\eta
\le r(T_Kf,T_Kg)\le r(f,g),                         \tag{DC.4}
```

where

```math
\Delta_{K,\eta}(f,g)
=r(f,g)-\max_{y\in E_{K,\eta}(f)}(f(y)-g(y))\ge0.  \tag{DC.5}
```

The sharpest same-witness form, with no threshold chosen, is

```math
\max_y\{f(y)-g(y)-e_K^f(y)\}
\le r(T_Kf,T_Kg)\le r(f,g).                         \tag{DC.5a}
```

In particular, the directed response is preserved exactly whenever an exact
maximizer of `f-g` is exposed by `f`.  Both directed responses of a pair are
preserved if this holds also after exchanging `f,g`.

For a sequence of continuations `T_t` and evolved functions
`f_t=T_tf_(t-1)`, `g_t=T_tg_(t-1)`, arbitrary choices `eta_t>=0` give

```math
r(f_T,g_T)\ge r(f_0,g_0)
-\sum_{t=1}^T
 \{\Delta_{K_t,\eta_t}(f_{t-1},g_{t-1})+\eta_t\}.  \tag{DC.6}
```

Thus a sublinear total exposure defect, rather than strict contraction at
each step, is enough to preserve the leading directed scale.

#### Proof

For the upper bound in (DC.4), choose a minimizer `y_g` for `T_Kg(x)`:

```math
(T_Kf)(x)-(T_Kg)(x)
\le f(y_g)-g(y_g)\le r(f,g).
```

For the lower bound, take any `y` and an `x` attaining (DC.3a).  Then

```math
(T_Kf)(x)=K(x,y)+f(y)-e_K^f(y),
\qquad
(T_Kg)(x)\le K(x,y)+g(y).
```

Subtract and maximize over `y` to obtain (DC.5a).  Restricting to
`e_K^f(y)<=eta` gives (DC.4).  Iteration proves (DC.6). `square`

The theorem is not the assumption that the directed table is preserved.  Its
falsifiable content is that only the intersection of a directed maximizing
set with an exposed-input set matters.  It can be much smaller than the full
continuation table.  It also identifies why generic nonexpansiveness is too
weak: strict contraction requires the continuation to hide every exact
maximizing witness of the directed difference.  Hiding them is not by itself
sufficient, because different inputs may reproduce the same final gap.

## 2. A finite-alphabet synchronization criterion

For an alphabet `f_a`, `a in A`, define the directed maximizing sets

```math
M_{ab}=\operatorname*{argmax}_{y\in Y}(f_a(y)-f_b(y)). \tag{DC.7}
```

### Corollary DC.2 (exposure incidence closes the directed table)

If

```math
M_{ab}\cap E_{K,0}(f_a)\ne\varnothing
\quad\text{for every ordered pair }(a,b),           \tag{DC.8}
```

then the entire directed table is unchanged:

```math
r(T_Kf_a,T_Kf_b)=r(f_a,f_b).                         \tag{DC.9}
```

More generally, if for every ordered pair there is
`y in E_(K,eta)(f_a)` with

```math
f_a(y)-f_b(y)\ge r(f_a,f_b)-\delta,
```

then every table entry loses at most `delta+eta`.

Condition (DC.8) stores an incidence certificate—one exposed witness per
ordered alphabet pair—not all rows of `K`.  It is useful only when a model
supplies those witnesses structurally; checking it by reconstructing the full
continued functions would provide no compression.

## 3. Exact interacting validation: metric-isometry bottlenecks

Let `(Y,d)` be a finite metric space, let `g` be an isometry of `Y`, and for
`lambda>0` define

```math
K_{\lambda,g}(a,t)=\lambda d(t,g(a)).               \tag{DC.10}
```

### Theorem DC.3 (metric-isometry bottleneck algebra)

For isometries `g,h` and positive strengths `lambda,mu`,

```math
K_{\lambda,g}\star K_{\mu,h}
=K_{\min\{\lambda,\mu\},h\circ g}.                 \tag{DC.11}
```

Every ordered pair of rows satisfies

```math
r(K_{\lambda,g}(a,\cdot),K_{\lambda,g}(b,\cdot))
=\lambda d(a,b).                                    \tag{DC.12}
```

Consequently an arbitrary chain closes exactly on the bottleneck strength
and isometry holonomy:

```math
K_1\star\cdots\star K_T
=K_{\lambda_*,g_T\circ\cdots\circ g_1},
\qquad \lambda_*=\min_i\lambda_i.                 \tag{DC.13}
```

The directed table after the chain is simply `lambda_* d`; the holonomy is
needed only for labelled endpoint responses.

#### Proof

After moving the second endpoint by `h^(-1)`, the composite entry is

```math
\min_u\{\lambda d(u,g(a))+\mu d(h^{-1}(t),u)\}.
```

The triangle inequality makes this at least
`min(lambda,mu)d(h^(-1)(t),g(a))`.  Choosing `u=g(a)` or
`u=h^(-1)(t)` attains the smaller of the two endpoint costs, proving
(DC.11).  The reverse triangle inequality gives an upper bound
`lambda d(a,b)` for (DC.12), and evaluating at `t=g(b)` attains it.  Induct
for (DC.13). `square`

The algebra is robust.  If actual kernels satisfy

```math
\|\widetilde K_i-K_{\lambda_i,g_i}\|_\infty\le\eta_i,
```

min-plus nonexpansiveness gives

```math
\left\|\widetilde K_1\star\cdots\star\widetilde K_T
-K_{\lambda_*,g_T\circ\cdots\circ g_1}\right\|_\infty
\le\sum_i\eta_i.                                   \tag{DC.14}
```

Every directed row response is therefore within `2 sum_i eta_i` of
`lambda_* d(a,b)`.  In particular, total perturbation
`o(lambda_* d(a,b))` preserves the leading directed response for the declared
pair through arbitrarily many interacting steps.

This state is strictly smaller than a generic transfer kernel whenever the
isometry family has a compact description.  At matched scalar precision it
uses one scalar plus a holonomy label, rather than `|Y|^2` unrelated values.

## 4. Permutation-Potts and Ising specializations

Let the spin alphabet be `[q]`.  For a strength `J>0` and permutation
`pi in S_q`, define the min-plus edge kernel

```math
P_{J,\pi}(a,t)=-J{\bf1}\{t=\pi(a)\}.               \tag{DC.15}
```

This is a zero-temperature Potts constraint whose preferred relation is
labelled by `pi`.

### Corollary DC.4 (nonabelian permutation-Potts law)

For `J,L>0` and `pi,rho in S_q`,

```math
P_{J,\pi}\star P_{L,\rho}
=-\max\{J,L\}\boldsymbol1
 +P_{\min\{J,L\},\rho\circ\pi}.                    \tag{DC.16}
```

Consequently a chain of such edges has endpoint kernel

```math
-\left(\sum_{i=1}^T J_i-\mu\right)\boldsymbol1
 +P_{\mu,\Pi},
\qquad
\mu=\min_iJ_i,
\quad
\Pi=\pi_T\circ\cdots\circ\pi_1.                   \tag{DC.17}
```

For distinct input labels, both oriented directed row responses of
`P_(J,pi)` equal `J`; after the chain they all equal `mu`.  Thus the additive
baseline, the nonabelian permutation holonomy, and one bottleneck scalar form
an exact composable endpoint state.  The directed table itself closes on the
single bottleneck scalar.

#### Proof

Use the discrete metric `d(a,b)={\bf1}{a\ne b}`.  Then
`P_(J,pi)=K_(J,pi)-J boldsymbol1`, so (DC.16)--(DC.17) are Theorem DC.3 with
the scalar offsets restored.  Two distinct rows differ by `-J` at `pi(a)`,
by `+J` at `pi(b)`, and by zero elsewhere. `square`

Theorem DC.1 explains the threshold exactly.  If `L>=J`, every input label is
exposed through the second permutation edge for every first-edge row, so all
directed witnesses survive.  If `L<J`, only that row's preferred input is
exposed; every ordered-pair maximizing witness is hidden and the new link
clips all gaps to `L`.

This is a genuine interacting continuation.  A middle spin is optimized out,
the gauge state has the noncommutative law of `S_q`, and a weak edge causes
irreversible contraction.  Under labelled endpoint-field queries, fixed
baseline and `mu`, the `q!` permutation kernels are pairwise response-distance
`mu`: a row on which two permutations differ exposes both signs of the
kernel difference.  Thus accuracy below `mu/2` needs `log_2(q!)=Theta(q log
q)` discrete bits, while `(Pi,mu)` suffices, compared with `q^2` unrelated
real entries at matched precision.  The directed row table alone forgets
`Pi` completely.

### Binary Ising form

Let the spin alphabet be `{+1,-1}` and use the min-plus edge kernel

```math
K_J(s,t)=-Jst.                                      \tag{DC.18}
```

Its two rows are functions `f_a(s)=K_J(a,s)`.  Continue them through a second
edge `K_L`.  A direct minimization gives

```math
(T_{K_L}f_a)(t)
=-\left|Ja+Lt\right|
=-c-J' at,                                          \tag{DC.19}
```

where

```math
c=\max\{|J|,|L|\},
\qquad
J'=\operatorname {sgn}(JL)\min\{|J|,|L|\}.         \tag{DC.20}
```

### Corollary DC.5 (binary signed bottleneck law)

The directed row table after a genuine interacting continuation obeys

```math
r(T_{K_L}f_+,T_{K_L}f_-)
=r(T_{K_L}f_-,T_{K_L}f_+)
=2\min\{|J|,|L|\}.                                  \tag{DC.21}
```

Along a chain with nonzero couplings `J_1,...,J_T`, its projective endpoint
kernel therefore has

```math
\text{holonomy sign}=\operatorname {sgn}\prod_tJ_t,
\qquad
\text{directed amplitude}=2\min_t|J_t|.             \tag{DC.22}
```

More precisely, if `mu=min_t|J_t|` and
`s=sgn(prod_t J_t)`, the endpoint kernel is

```math
-\left(\sum_t|J_t|-\mu\right)-s\mu\,a t.           \tag{DC.23}
```

The additive baseline, the holonomy sign, and the bottleneck magnitude form
an exact three-scalar composable state for endpoint energies; the directed
distance table itself needs only the bottleneck magnitude.  Appending an edge
with `|L|>=|J|` preserves the old table exactly, while a weaker edge contracts
it to a completely quantified new table.  There is no cumulative additive
loss.

#### Proof

Equations (DC.19)--(DC.20) follow by evaluating the two possible intermediate
spins.  The difference of the two rows is `-2J't`, so both oriented maxima
are `2|J'|`.  For the full endpoint formula, suppose a prefix has kernel
`-b-s mu at`, where `b=sum_(i<T)|J_i|-mu`.  Composing with `J_T` adds
`-max(mu,|J_T|)` to the baseline, multiplies the sign by `sgn(J_T)`, and
replaces `mu` by `min(mu,|J_T|)`.  In either ordering of the magnitudes, the
new baseline is `sum_(i<=T)|J_i|-min_(i<=T)|J_i|`.  This proves
(DC.22)--(DC.23) by induction. `square`

This is not an identity continuation: the intermediate spin is optimized
out, the coupling sign composes as a nontrivial gauge cocycle, and the weaker
link irreversibly contracts the response.  Theorem DC.1 predicts the sharp
threshold.  When `|L|>=|J|`, both directed maximizing inputs remain exposed;
when `|L|<|J|`, both ordered-pair maximizing witnesses are hidden and the
transmitted amplitude is clipped to `2|L|`.

## 5. Limits of the result

1. Exposure is a sufficient certificate, not a necessary characterization:
   different minimizing inputs can sometimes reproduce the same final
   directed difference.
2. A common minimizing state, uniqueness, or a small minimizer fibre alone
   does not imply (DC.8); the fibre must meet the relevant *directed* exposed
   face.
3. Generic repeated continuations can incur linear total exposure defect.
   Equation (DC.6) states exactly what must be controlled and does not turn
   nonexpansiveness into convergence.
4. The Ising bottleneck state answers endpoint ground-state queries.  It does
   not preserve multiplicities, finite-temperature partition functions, or a
   future interaction that bypasses the exposed endpoint.
5. Even a unique or common minimizer does not force preservation.  With a
   singleton output, `K=0`, `f=(0,1)`, and `g=(1,0)`, both continued values
   are zero although the input oscillation is two: both directed maximizers
   are hidden.

The next theorem-level question is whether a higher-state CSP or code family
has a structural exposure incidence with sublinear cumulative defect.  The
binary chain proves that the directed algebra can survive a nonproduct
operation; it does not yet show that the mechanism scales to growing
interfaces.
