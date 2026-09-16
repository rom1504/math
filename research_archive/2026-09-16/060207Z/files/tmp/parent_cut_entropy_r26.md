### Wave 26, Route 2: parent-cut likelihood transport and its generic bottleneck

**Status.** The likelihood identities, heat-bath normalizations, exact
two-mode obstruction, low-temperature asymptotics, and finite enumerations
below are **Verified** by `tmp/parent_cut_entropy_r26.py`.  The
minimizer-specific target theorem in Section 8 is **Open**.  No bound of
`O(n^{1/2-2c})` is proved.

#### 1. The parent cut is exactly a barycenter likelihood ratio

Let `R=U_m\otimes\nu_\beta` on full lifted pairs `(S,d)` and put

```math
L(S,d)=C e^{-\ell_S(d[S])}
={d\widehat M_0\over dR}(S,d),
\qquad
f(d)=Ca(d).
\tag{R26.1}
```

The endpoint identities can then be written as conditional expectations:

```math
\boxed{
f(D)=\mathbb E_R[L\mid D],
\qquad
w_S=\mathbb E_R[L\mid S],
\qquad
\mathbb E_RL=1.
}
\tag{R26.2}
```

Moreover, with `\pi_0(S)=U_m(S)w_S`,

```math
\boxed{
{q_S(d)\over\nu_\beta(d)}={L(S,d)\over w_S}
={\mu_{\gamma,S}(d[S])\over\nu_{\beta,S}(d[S])},
\qquad
f\nu_\beta=\sum_S\pi_0(S)q_S.
}
\tag{R26.3}
```

Thus the missing parent law is the barycenter of the normalized lifted
channels.  Data processing/convexity gives

```math
D(f\nu_\beta\Vert\nu_\beta)
\le\mathbb E_{S\sim\pi_0}D(q_S\Vert\nu_\beta),
\tag{R26.4}
```

but there is an exact and informative gap.  Under `\widehat M_0`,

```math
\boxed{
\mathbb E_{\pi_0}D(q_S\Vert\nu_\beta)
=\operatorname{Ent}_{\nu_\beta}(f)+I_{\widehat M_0}(S;D).
}
\tag{R26.5}
```

Indeed, the two chain rules for `D(\widehat M_0\Vert U_m\otimes\nu_\beta)`
are

```math
\operatorname{Ent}_{\nu_\beta}(f)+I(S;D)+D(\pi_0\Vert U_m)
=D(\pi_0\Vert U_m)+\mathbb E_{\pi_0}D(q_S\Vert\nu_\beta).
```

Consequently

```math
\boxed{
\operatorname{Ent}_{\nu_\beta}(Ca)
=\mathbb E_{\pi_0}
D(\mu_{\gamma,S}\Vert\nu_{\beta,S})-I(S;D).
}
\tag{R26.6}
```

Equation (R26.6) is the sharp data-processing formulation.  Bounding the
average channel KL alone simply returns the original full endpoint target;
one needs either a direct barycenter theorem or a lower bound showing that
mutual information cancels nearly all of the average channel cost.

#### 2. Exact parent heat-bath forms

Encode the `2^n` oriented cuts by the global orientation and switches of
vertices `1,\ldots,n-1`.  Write `T_i d` for a coordinate flip.  The global
orientation flip has

```math
E(T_0d)-E(d)=-2E(d),
```

while a vertex switch has

```math
E(T_id)-E(d)=-4H_i(d),
\qquad H_i(d)=\sum_{j\ne i}a_{ij}d_{ij}.
\tag{R26.7}
```

For an unordered cube edge `e=\{d,d'\}`, define

```math
s_e=\nu(d)+\nu(d'),
\qquad
c_e={\nu(d)\nu(d')\over\nu(d)+\nu(d')},
\qquad
\bar f_e={\nu(d)f(d)+\nu(d')f(d')\over s_e}.
```

With every coordinate updated at rate one, the sum of one-site conditional
entropies and the heat-bath modified Dirichlet form are

```math
\boxed{
\begin{aligned}
\mathfrak A_\nu(f)
&=\sum_e\left[
\nu(d)f(d)\log{f(d)\over\bar f_e}
+\nu(d')f(d')\log{f(d')\over\bar f_e}
\right],\\
\mathfrak J_\nu(f)
&=\sum_e c_e(f(d)-f(d'))(\log f(d)-\log f(d')).
\end{aligned}}
\tag{R26.8}
```

If instead one coordinate is chosen uniformly at total rate one, both forms
in (R26.8) are divided by `n`; the corresponding log-Sobolev constants are
multiplied by `n`.

There is an exact two-point decomposition.  Let `p_e` be the conditional
Gibbs law on the pair and let `r_e` be its tilt by `f`.  Then

```math
\boxed{
c_e(f-f')(\log f-\log f')
=s_e\bar f_e
\left[D(r_e\Vert p_e)+D(p_e\Vert r_e)\right].
}
\tag{R26.9}
```

The first KL on the right is exactly the edge contribution to
`\mathfrak A_\nu(f)`.  Hence

```math
\boxed{\mathfrak A_\nu(f)\le\mathfrak J_\nu(f).}
\tag{R26.10}
```

The parent ordinary/Hellinger form is similarly

```math
\mathfrak H_\nu(f)
=\sum_e c_e(\sqrt{f(d)}-\sqrt{f(d')})^2.
\tag{R26.11}
```

Unlike the Johnson constant, no useful generic parent-cube constant is
available in this low-temperature regime: the conditional Gibbs
probabilities in (R26.7) can be exponentially unbalanced when `\beta n` is
large.

#### 3. The precise heat-bath sufficient theorem

For diagnostics, define the `f`-specific approximate-tensorization and
modified-LSI ratios

```math
K_{\rm AT}(f,\nu)
={\operatorname{Ent}_\nu(f)\over\mathfrak A_\nu(f)},
\qquad
K_{\rm HB}(f,\nu)
={\operatorname{Ent}_\nu(f)\over\mathfrak J_\nu(f)},
\tag{R26.12}
```

with the usual zero conventions.  A non-tautological sufficient theorem in
the sum-rate normalization (R26.8) is: there are constants `K_0,K_1<\infty`,
uniform over the active ratio window, target orders, and chosen endpoints,
such that

```math
\boxed{
\operatorname{Ent}_{\nu_\beta}(f)
\le K_0\mathfrak A_{\nu_\beta}(f),
\qquad
\mathfrak A_{\nu_\beta}(f)\le K_1n^{1/2-2c}.
}
\tag{R26.13}
```

Then the parent term is at most `K_0K_1n^{1/2-2c}`.  The modified-LSI
alternative is the equally explicit pair

```math
\boxed{
\operatorname{Ent}_{\nu_\beta}(f)
\le K_0\mathfrak J_{\nu_\beta}(f),
\qquad
\mathfrak J_{\nu_\beta}(f)\le K_1n^{1/2-2c}.
}
\tag{R26.14}
```

The product `K_0K_1`, rather than either constant separately, is the relevant
endpoint constant.  A global LSI for every function is far stronger than
necessary; these inequalities need hold only for the specific endpoint
likelihood `f`.

Conditional entropy is convex and homogeneous.  Since
`f=\mathbb E_UL_S`,

```math
\boxed{
\mathfrak A_{\nu_\beta}(f)
\le\mathbb E_{S\sim U_m}\mathfrak A_{\nu_\beta}(L_S).
}
\tag{R26.15}
```

After normalizing `L_S/w_S=q_S/\nu_\beta`, each summand is `w_S` times the
sum of conditional one-coordinate KLs between `q_S` and `\nu_\beta`.  A
vertex coordinate outside `S` contributes zero because `L_S` depends only
on `d[S]`.  The remaining terms, however, are under `q_S`, not the parent
marginal in (10.746).  Thus (10.746) does not bound (R26.15); this is the
same wrong-marginal obstruction in an exact heat-bath form.

#### 4. Exact signing obstruction to temperature-uniform selector absorption

The Wave 25 `A_4` value only defeated two displayed coefficients.  There is
a stronger statement if a comparison is required uniformly over
temperatures.

Set `\gamma=0` and keep one deletion.  Then every child partition function
is `2^m`, so `w_S=1` and `\pi_0=U_m`.  As `\beta\to\infty`, the endpoint law
chooses `S` and a child cut `y` uniformly and then chooses a maximum-energy
parent extension of `y`.  Put

```math
\overline E_\infty
=\mathbb E_{S\sim U_m,\,y\sim U(\mathcal D_S)}
\max_{d:d[S]=y}\langle A,d\rangle.
\tag{R26.16}
```

On the finite cut space, direct Gibbs asymptotics give

```math
\boxed{
\operatorname{Ent}_{\nu_\beta}(Ca)
=\beta\bigl(Q(A)-\overline E_\infty\bigr)+O_A(1).
}
\tag{R26.17}
```

Indeed, `\log Z_A(\beta)=\beta Q(A)+O_A(1)`, the endpoint law converges,
and its mean energy tends to `\overline E_\infty`.  Meanwhile the actual
conditional selector entropy is at most `\log n`, and, because `w_S=1`,
every mean adjacent `H^2(q_S,q_T)` is at most one.  Therefore any positive
deficit in (R26.17) rules out a constant-factor absorption into selector
entropy or selector Hellinger uniformly in `\beta`.

Exact enumeration for the four exact minimizers gives

| signing | `Q(A)` | `\overline E_\infty` | slope `Q-\overline E_\infty` |
|---|---:|---:|---:|
| `A_4` | `8` | `3` | `5` |
| `A_6` | `10` | `15/4` | `25/4` |
| `A_8` | `20` | `35/8` | `125/8` |
| `A_9` | `24` | `35/8` | `157/8` |

For example, on `A_9` at `\beta=8`, the cut entropy is
`154.0587248645`, the selector entropy is `0.4254410837`, and the mean
finite Hellinger edge is `0.3445333798`.  This is a genuine exact-minimizer
endpoint obstruction, stronger than the old coefficient-one `A_4` audit.

It is not a target-regime falsifier: `\gamma=0` supplies no child distortion
control, whereas the intended theorem may impose
`\beta\asymp\gamma\asymp n^{-1/2+c}`.  It proves that any viable comparison
must use that parameter relation and exact-minimizer structure, rather than
claim a temperature-uniform information inequality.

#### 5. Positive two-mode obstruction to generic parent-cube transport

The likelihood/conditional-expectation structure alone does not yield a
dimension-free heat-bath theorem.  Fix `0<\varepsilon<1/4` and put on the
two-bit cube

```math
\nu_\varepsilon(00)=\nu_\varepsilon(11)={1-2\varepsilon\over2},
\qquad
\nu_\varepsilon(01)=\nu_\varepsilon(10)=\varepsilon.
\tag{R26.18}
```

This is an exact positive ferromagnetic two-spin Gibbs law.  Fix
`0<\delta<1`, let `A_\delta=2-\delta`, `B_\delta=\delta`, and define two
strictly positive leave-one-coordinate likelihoods

```math
L_1(x)=
\begin{cases}A_\delta,&x_1=0,\\B_\delta,&x_1=1,\end{cases}
\qquad
L_2(x)=
\begin{cases}A_\delta,&x_2=0,\\B_\delta,&x_2=1.\end{cases}
\tag{R26.19}
```

Both satisfy `\mathbb E_{\nu_\varepsilon}L_i=1`.  Taking a uniform selector
and `f=(L_1+L_2)/2` therefore satisfies exactly the normalization and
conditional-expectation identities (R26.1)--(R26.3).  If

```math
K_\delta=A_\delta\log A_\delta+B_\delta\log B_\delta,
```

then direct calculation gives

```math
\boxed{
\operatorname{Ent}_{\nu_\varepsilon}(f)
={1-2\varepsilon\over2}K_\delta,
\qquad
I(S;X)=\varepsilon K_\delta,
}
\tag{R26.20}
```

and

```math
\boxed{
H^2(\nu L_1,\nu L_2)
=\varepsilon(\sqrt{A_\delta}-\sqrt{B_\delta})^2.
}
\tag{R26.21}
```

The parent heat-bath Jeffreys form is exactly

```math
\boxed{
\mathfrak J_{\nu_\varepsilon}(f)
=2\varepsilon(1-2\varepsilon)(1-\delta)
\log{A_\delta\over B_\delta}.
}
\tag{R26.22}
```

The conditional-entropy sum `\mathfrak A` also satisfies
`\mathfrak A=\Theta_\delta(\varepsilon)`.  More explicitly, with

```math
u_A=(1-2\varepsilon)A_\delta+2\varepsilon,
\qquad
u_B=(1-2\varepsilon)B_\delta+2\varepsilon,
```

it equals twice

```math
{1-2\varepsilon\over2}A_\delta\log{A_\delta\over u_A}
+\varepsilon\log{1\over u_A}
+\varepsilon\log{1\over u_B}
+{1-2\varepsilon\over2}B_\delta\log{B_\delta\over u_B}.
\tag{R26.23}
```

Hence the cut entropy tends to the positive constant `K_\delta/2`, while
the selector entropy, selector Hellinger edge, parent conditional-entropy
sum, parent heat-bath Jeffreys form, and parent heat-bath Hellinger form all
vanish linearly in `\varepsilon`.  In particular,

```math
{\operatorname{Ent}_{\nu_\varepsilon}(f)
 \over\mathfrak J_{\nu_\varepsilon}(f)}
={K_\delta
 \over4\varepsilon(1-\delta)\log(A_\delta/B_\delta)}
\longrightarrow\infty.
\tag{R26.24}
```

This is a sharp positive bottleneck, not a zero-support artifact.  It
falsifies every dimension-free parent approximate-tensorization, ordinary
LSI, modified-LSI, or selector-absorption bound purportedly deduced only
from the abstract endpoint likelihood structure; data processing by itself
cannot add the missing bound.  It is not generated by a complete-signing
endpoint, so it leaves a minimizer-specific theorem open.

#### 6. Named exact-minimizer audit

At the mismatched Wave 25 `A_4` parameters
`(\beta,\gamma)=(1/2,1/10)`, the parent forms are

| quantity | value |
|---|---:|
| cut entropy | `0.627697022731` |
| selector entropy | `0.147161088198` |
| `\mathfrak A_{\nu_\beta}(f)` | `0.689538798181` |
| `\mathfrak J_{\nu_\beta}(f)` | `1.097205021938` |
| parent heat-bath Hellinger form | `0.247248174955` |
| parent Metropolis Jeffreys form | `1.185870372681` |

Thus `A_4` breaks the two selector coefficients but does **not** break the
coefficient-one parent conditional-entropy or heat-bath Jeffreys bounds.

At matched temperature, the sampled ratios are:

| signing | `\beta` | cut entropy / `\mathfrak A` | cut entropy / `\mathfrak J` |
|---|---:|---:|---:|
| `A_4` | `0.5` | `1.395431` | `0.721944` |
| `A_6` | `0.5` | `0.373929` | `0.177738` |
| `A_8` | `0.5` | `1.027687` | `0.531332` |
| `A_9` | `0.5` | `0.833751` | `0.404226` |
| `A_9` | `4` | `1.142446` | `0.565912` |
| `A_9` | `8` | `1.142446` | `0.565912` |

For `A_4,A_6,A_8`, the matched cut entropy tends to zero as
`\beta\to\infty`.  For `A_9` it saturates at `0.0356868189`; the parent
conditional and Jeffreys forms also saturate, giving the last two stable
ratios.  The 25 `A_9` parent grounds have 26 cube-flip edges in components
of sizes `15,8,2`, and the endpoint density on the ground face ranges from
`0.44642857` to `1.33928571`.  Thus the named matched examples neither prove
nor falsify an `O(1)` minimizer-specific parent heat-bath comparison.

These are finite numerical diagnostics, not asymptotic evidence at the
target temperature scaling.

#### 7. Why exact minimality has not yet closed the form

The rates in (R26.8) contain logistic factors for the energy changes
(R26.7).  The cap `E(d)\le Q(A)` controls neither the conductance across a
specific low-mass cut of the parent cube nor the endpoint likelihood
variation across such a cut.  Exact minimality therefore does not by itself
give an `O(1)` parent LSI constant.

One could replace coordinate flips by nonlocal competitor moves connecting
the relevant ground components, but then a new theorem must show both:

1. sufficient conductance under `\nu_\beta`; and
2. an `O(n^{1/2-2c})` Dirichlet form for the specific barycenter density
   `f=Ca`.

Neither follows from the parent cap or (10.746).  The abstract obstruction
shows why both ingredients, not just a formal LSI invocation, are necessary.

#### 8. Minimal surviving parent theorem

The exact unresolved statement remains

```math
\boxed{
D\left(\sum_S\pi_0(S)q_S\middle\Vert\nu_\beta\right)
=O(n^{1/2-2c})
}
\tag{R26.25}
```

for target-specific exact minimizers, fixed retention density, and the
chosen target relation between `\beta` and
`\gamma=\Theta(n^{-1/2+c})`.  This is exactly
`\operatorname{Ent}_{\nu_\beta}(Ca)=O(n^{1/2-2c})` rewritten as a
barycenter theorem.

A non-tautological heat-bath version sufficient for (R26.25) is:

```math
\boxed{
\operatorname{Ent}_{\nu_\beta}(f)
\le K\,\mathfrak J_{\nu_\beta}(f),
\qquad
K\,\mathfrak J_{\nu_\beta}(f)=O(n^{1/2-2c}),
}
\tag{R26.26}
```

with constants uniform over the active ratio window and target orders.  It
is enough that (R26.26) hold for this specific `f`; a global Gibbs LSI is not
needed.  The analogous `\mathfrak A` or nonlocal-competitor statement is
equally sufficient.

Falsification requires an asymptotic exact-minimizer family in the actual
target parameter regime for which the left side of (R26.25) is
`\omega(n^{1/2-2c})`.  The `\gamma=0` signing examples and the abstract
two-mode Gibbs example sharply falsify generic shortcuts, but do not meet
that criterion.

#### 9. Disposition

The parent term is a barycenter bias, not selector variation.  Data
processing gives the exact identity (R26.5), but no saving.  Generic parent
cube LSI is blocked by the positive two-mode example, and uniform
constant-factor selector absorption is blocked within exact signing
endpoints by (R26.17).  The named matched-temperature examples leave only a
narrow, explicitly minimizer- and parameter-specific heat-bath or nonlocal
competitor theorem.  The adjacent Johnson Hellinger estimate remains a
separate open term throughout.
