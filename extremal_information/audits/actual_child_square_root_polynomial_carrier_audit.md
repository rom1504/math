# Square-root polynomial carriers: pressure recovery and the entropy blockade

Status: **rigorous pressure-side theorem; full variational recovery remains
open**.  Projecting the square root of an actual optimal row density and
squaring the projection produces an exactly nonnegative bounded-degree
polynomial density.  Its replacement cost against the actual bridge pressure
is `O(d^(-1/2))` per row, uniformly in the child orders.  This removes the
positive-part Fourier leakage for the pressure term.

The same argument does not recover the entropy term.  The available
Hellinger and square-root Dirichlet bounds do not imply entropy continuity;
an explicit polynomial spike shows that such an inference is false in
general.  No claim is made that this spike is itself the projection of an
actual optimal density.  A new projection-specific entropy lemma would be
needed to turn the construction into a feasible recovery of the full
row-product variational value.

## 1. Square-root Fourier truncation

Let `f` be a strictly positive density on the fair `n`-cube and suppose

```math
E_Uf=1,
\qquad \|f\|_2\le K,
\qquad
|\log f(x)-\log f(x^{(j)})|\le {A\over\sqrt n}.     \tag{SQ.1}
```

Put

```math
r=\sqrt f,
\qquad g_d=\Pi_{\le d}r,
\qquad \varepsilon_d=\|r-g_d\|_2,
\qquad z_d=E_Ug_d^2=1-\varepsilon_d^2,              \tag{SQ.2}
```

and, when `z_d>0`, define

```math
q_d={g_d^2\over z_d}.                               \tag{SQ.3}
```

Then `q_d` is a probability density and has exact Walsh degree at most `2d`.

**Lemma SQ.1 (uniform Hellinger truncation).**  One has

```math
\boxed{
\varepsilon_d^2\le {A^2\over16(d+1)},
\qquad
\|q_d-f\|_1\le2\varepsilon_d+\varepsilon_d^2,
\qquad
\|\sqrt{q_d}-\sqrt f\|_2
\le\varepsilon_d+\varepsilon_d^2.}                \tag{SQ.4}
```

Moreover, with the standard cube log-Sobolev normalization,

```math
\boxed{
D(q_dU\Vert U)\le {A^2\over8z_d},
\qquad
D(fU\Vert U)\le2\log K.}                          \tag{SQ.5}
```

*Proof.*  If adjacent values of `r` are `a,b>0`, (SQ.1) gives

```math
{|a-b|\over a+b}
\le\tau_n:=\tanh {A\over4\sqrt n}.
```

Consequently

```math
\sum_jE_U(r-r^{(j)})^2
\le4n\tau_n^2\le {A^2\over4}.                     \tag{SQ.6}
```

The Walsh Dirichlet identity and projection above degree `d` imply the first
bound in (SQ.4).  Orthogonality gives `z_d=1-epsilon_d^2`.  Also

```math
\|r^2-g_d^2\|_1
\le\|r-g_d\|_2\|r+g_d\|_2
\le\varepsilon_d(1+\sqrt {z_d}),                  \tag{SQ.7}
```

while normalization changes `g_d^2` in `L^1` by `1-z_d=epsilon_d^2`.
This proves the second bound.  Since `r>=0`, pointwise
`|r-|g_d||<=|r-g_d|`; rescaling `|g_d|` from norm `sqrt(z_d)` to one costs
`1-sqrt(z_d)<=epsilon_d^2`.  This proves the Hellinger bound.

Projection can only decrease the quadratic Dirichlet form.  Cube
log-Sobolev and (SQ.6) therefore give

```math
\operatorname {Ent}_U(g_d^2)
\le {1\over2}\sum_jE_U(g_d-g_d^{(j)})^2
\le {A^2\over8}.
```

Division by `z_d` proves the first part of (SQ.5).  The second is monotonicity
of Renyi divergences:
`D(fU||U)<=log E_Uf^2<=2log K`. `square`

## 2. Dimension-free pressure replacement

Let `G` be any row function whose oscillation under one bit flip is at most
`2u`, and put `X=G-E_UG`.  Bounded differences gives

```math
\log E_Ue^{sX}\le {s^2nu^2\over2}.                 \tag{SQ.8}
```

In particular, direct integration of the subgaussian tail gives

```math
E_U\exp\{X^2/(4nu^2)\}\le3.                       \tag{SQ.9}
```

Entropy duality consequently gives, for every density `p`,

```math
E_pX^2\le4nu^2\{D(pU\Vert U)+\log3\}.             \tag{SQ.10}
```

**Theorem SQ.2 (exact-polynomial pressure recovery).**  If `z_d>=1/2`, then

```math
\boxed{
|E_{q_d}G-E_fG|
\le2\sqrt2\,u\sqrt n
 (\varepsilon_d+\varepsilon_d^2)
 \left\{{A^2\over8z_d}+2\log K+2\log3\right\}^{1/2}.}         \tag{SQ.11}
```

*Proof.*  Constants cancel because both densities have mass one.  Thus, by
Cauchy--Schwarz and `(sqrt q+sqrt f)^2<=2(q+f)`,

```math
|E(q_d-f)X|
\le\|\sqrt{q_d}-\sqrt f\|_2
 \{2(E_{q_d}X^2+E_fX^2)\}^{1/2}.                  \tag{SQ.12}
```

Apply (SQ.4), (SQ.5), and (SQ.10). `square`

Apply this sequentially to a product of `m` row factors.  Averaging the
actual pressure `L` over every other current factor leaves a row function
with one-bit oscillation at most `2u`.  At physical scale
`u=beta/sqrt(N)`, (SQ.4) and (SQ.11) give

```math
\boxed{
|E_{\otimes_iq_{i,d}}L-E_{\otimes_if_i}L|
\le C_{A,K,\beta}{m\sqrt{n/N}\over\sqrt{d+1}}
\le C_{A,K,\beta}{N\over\sqrt{d+1}}.}             \tag{SQ.13}
```

This theorem concerns the exact actual-child pressure and uses no cumulant
or cluster expansion.

## 3. Why the entropy conclusion does not follow

The missing estimate is

```math
D(q_dU\Vert U)\le D(fU\Vert U)+o_d(1)              \tag{SQ.14}
```

uniformly in the cube dimension.  Lemma SQ.1 does not prove it.  Although
`q_d` is close to `f` in total variation and Hellinger distance, the only
immediate collision bound is the degree-dependent hypercontractive estimate

```math
\|q_d\|_2
={\|g_d\|_4^2\over z_d}
\le {3^d\over z_d}.                                \tag{SQ.15}
```

The dimension-free entropy modulus on a fixed `L^2` ball therefore has a
radius that grows exponentially with `d`; combining (SQ.15) with the
`O(d^(-1/2))` Hellinger error does not tend to zero.
In particular, `q_d` is not known to obey the original fixed carrier-norm
constraint `||q_d||_2<=K`; exact polynomial positivity alone does not make it
a feasible member of that class.

This is not merely a cosmetic weakness of one continuity inequality.
Hellinger closeness plus bounded square-root Dirichlet energy does not
control entropy uniformly.  On a `d`-cube let `C` be the all-plus vertex,
`p=2^(-d)`, `a=d^(-1/2)`, and set

```math
s=1+a p^{-1/2}\mathbf1_C,
\qquad q={s^2\over E_Us^2}.                        \tag{SQ.16}
```

Then `s` and `q` are polynomials of degree at most `d`, and direct
calculation gives

```math
\|\sqrt q-1\|_2=\Theta(d^{-1/2}),
\qquad
\sum_jE(\sqrt q-(\sqrt q)^{(j)})^2=O(1),           \tag{SQ.17}
```

but

```math
\boxed{D(qU\Vert U)\longrightarrow\log2.}          \tag{SQ.18}
```

Indeed `E s^2=1+2a\sqrt p+a^2`, the `q`-mass of `C` is
`a^2+o(a^2)`, and its log density there is
`d\log2-\log d+o(d)`, whose product tends to `log2`; the complement
contributes `o(1)`.

Example (SQ.16) is **not** asserted to equal
`(Pi_(<=d)sqrt f)^2/E(Pi_(<=d)sqrt f)^2` for an actual optimal `f`.
It proves only that the bounds established in Lemma SQ.1 cannot yield
(SQ.14): projection geometry or another actual-factor identity must be used.

Therefore the square-root construction supplies a genuine finite-row-Fourier
carrier for the pressure term, but not yet a feasible approximation to the
entropy-regularized product value.  The finite-degree child-closure SML
remains open at exactly this entropy/projection step, even before addressing
the exponentially many cross-row response coefficients.
