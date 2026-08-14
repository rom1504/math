# Microcanonical disorder-counting composition

**Status:** verified exact counting theorem; the theorem strengthens the
annealed bridge identity by giving many good parent signings, but it retains
the changing-temperature obstruction and does not control the endpoint needed
for convergence.

## 1. Definitions and exact bridge average

Let `\mathcal A_n={\{\pm1\}}^{E_n}` be the labeled edge signings of `K_n`,
where `E_n=\binom n2`, and put

```math
\overline Z_A(t)=2^{-n}\sum_{x\in\{\pm1\}^n}
 \cosh(tH_A(x)),\qquad
R_A(t)=\log\overline Z_A(t)-E_n\log\cosh t.       \tag{1}
```

Write

```math
N_n(t,u)=\#\{A\in\mathcal A_n:R_A(t)\le u\},
\qquad
f_n(t,u)=2^{-E_n}N_n(t,u).                         \tag{2}
```

For `A\in\mathcal A_m`, `B\in\mathcal A_n`, an `m\times n` bridge
`C`, and `\epsilon\in\{\pm1\}`, let

```math
S(A,B,\epsilon,C)=
\begin{pmatrix}A&C\\ C^{\mathsf T}&\epsilon B\end{pmatrix}.
```

The bridge-and-orientation average proved in the soft-cap audit becomes,
after the exact edge centering in (1),

```math
2^{-(mn+1)}\sum_{\epsilon,C}
  \exp R_{S(A,B,\epsilon,C)}(t)
=\exp\bigl(R_A(t)+R_B(t)\bigr).                    \tag{3}
```

Indeed, averaging the bridge signs contributes `(\cosh t)^{mn}` and
averaging `\epsilon` uses
`[\cosh(a+b)+\cosh(a-b)]/2=\cosh(a)\cosh(b)`.
The cross-edge factor is then cancelled exactly by the `E_{m+n}` centering.

## 2. Exact lower-tail product theorem

For all `m,n\ge2`, `t\in\mathbb R`, `u,v\in\mathbb R`, and
`\lambda>0`,

```math
\boxed{
N_{m+n}(t,u+v+\lambda)
\ge (1-e^{-\lambda})\,2^{mn}
       N_m(t,u)N_n(t,v).}                           \tag{4}
```

In particular, at `\lambda=\log2`,

```math
\boxed{
N_{m+n}(t,u+v+\log2)
\ge 2^{mn-1}N_m(t,u)N_n(t,v).}                      \tag{5}
```

Equivalently, the lower-tail fractions satisfy

```math
f_{m+n}(t,u+v+\lambda)
\ge(1-e^{-\lambda})f_m(t,u)f_n(t,v).                \tag{6}
```

### Proof, including the orientation multiplicity

The set `\{B:R_B(t)\le v\}` is closed under `B\mapsto-B`, because
`H_{-B}=-H_B` and the partition function uses `\cosh`.  Choose one
representative from every antipodal pair; there are `N_n(t,v)/2` such
representatives.

Fix a good `A` and one representative `B`.  Equation (3) and Markov's
inequality show that at least a fraction `1-e^{-\lambda}` of the
`2^{mn+1}` pairs `(\epsilon,C)` obey

```math
R_{S(A,B,\epsilon,C)}(t)\le u+v+\lambda.            \tag{7}
```

The map

```math
(A,B,\epsilon,C)\longmapsto S(A,B,\epsilon,C)       \tag{8}
```

is injective after this choice of antipodal representatives.  The first
principal block recovers `A`, the bridge block recovers `C`, and the second
principal block `D` uniquely recovers the representative and orientation
from `D=\epsilon B`.  Notice that replacing `A` by `-A` does not create a
collision, because the first principal block changes.  Consequently the
number of distinct parents in (7) is at least

```math
N_m(t,u)\,{N_n(t,v)\over2}
(1-e^{-\lambda})2^{mn+1},
```

which is (4).  Without choosing representatives, the same proof says that
the raw parameterization has multiplicity at most two, exactly from
`(B,\epsilon)` and `(-B,-\epsilon)`.

## 3. What ordinary subadditivity this does give

For fixed physical `t>0` and `a\in\mathbb R`, set

```math
g_n(t,a)=f_n(t,an-\log2),\qquad
J_n(t,a)=-\log g_n(t,a).
```

Equation (5) gives

```math
J_{m+n}(t,a)+\log2
\le[J_m(t,a)+\log2]+[J_n(t,a)+\log2].               \tag{9}
```

Thus, whenever the events are eventually nonempty, Fekete's lemma gives a
limit for `J_n(t,a)/n`.  This formal rate is nevertheless trivial at every
fixed `t>0`: it equals zero for every fixed `a`.

To see this, choose a signing uniformly.  For each fixed spin `x`, `H_A(x)`
is a sum of `E_n` independent signs, so for every fixed
`C>\sqrt{\log2}` a union bound gives

```math
\Pr\!\left(\max_x|H_A(x)|>Cn^{3/2}\right)
\le 2^{n+1}\exp\!\left(-{C^2n^3\over2E_n}\right)
=e^{-\Omega(n)}.                                   \tag{10}
```

On the complementary event,

```math
R_A(t)\le tCn^{3/2}-E_n\log\cosh t,                 \tag{11}
```

which is below `an-\log2` for all sufficiently large `n`.  Hence
`g_n(t,a)\to1` exponentially and the rate furnished by (9) is zero.  The
additive thresholds in (4) lie far above the fixed-`t` endpoint.

More precisely, the soft-cap sandwich gives

```math
\min_A R_A(t)
=-E_n\log\cosh t+tM_n+O(n).                         \tag{12}
```

Its order-`n^2` leading term is universal.  Combining two endpoint child
thresholds in (4) omits the parent cross term
`-mn\log\cosh t`; at fixed nonzero `t` this is order `n^2`.  Therefore
(4) is not an endpoint-aligned large-deviation composition.

## 4. The project diagonal retains parameter contraction

Define the scaled-parameter count

```math
N_n^{(\beta)}(u)=N_n(\beta/\sqrt n,u).
```

Applying (4) at the parent's physical temperature
`t=\beta/\sqrt{m+n}` gives exactly

```math
N_{m+n}^{(\beta)}(u+v+\lambda)
\ge(1-e^{-\lambda})2^{mn}
N_m^{(\beta\sqrt{m/(m+n)})}(u)
N_n^{(\beta\sqrt{n/(m+n)})}(v).                    \tag{13}
```

Thus a balanced split sends both child parameters to
`\beta/\sqrt2`.  Iteration drives the leaves to `\beta=0`, just as in the
scalar centered-pressure recurrence.  The counting strengthening supplies
many parents at the recursively produced threshold, but it supplies neither
a converse inequality nor a comparison that restores the child parameter to
`\beta`.  It therefore does **not** imply convergence of the minimized free
energy at fixed `\beta`.

At this diagonal temperature the omitted endpoint cross term is

```math
mn\log\cosh{\beta\over\sqrt{m+n}}
=\frac{\beta^2mn}{2(m+n)}+O(1),                    \tag{14}
```

which is order `m+n` for a balanced split: precisely the leading scale of
the desired pressure, not a summable error.

## 5. Is a speed-`n^2` disorder LDP the right object?

There are two distinct scales, and they should not be conflated.

First, `R_A` is invariant under vertex switching and under `A\mapsto-A`.
For `n\ge3` the augmented switching action is free and every orbit has
`2^n` elements.  Hence every level set is a union of equal-size orbits; the
raw fraction `f_n` is exactly the same as the fraction in the switching
quotient.  Switching removes only `n` bits from an ambient disorder space
with `E_n=\Theta(n^2)` bits.  It therefore does not, by itself, rule out a
speed-`n^2` rarity statement for exact minimizers.

Second, at the project temperature a single edge flip changes `R_A` by at
most `2\beta/\sqrt n`.  McDiarmid's inequality consequently gives, for
uniform disorder,

```math
\Pr\bigl(|R_A-\mathbb ER_A|\ge sn\bigr)
\le2\exp\!\left(-{s^2n^2\over\beta^2(n-1)}\right). \tag{15}
```

This makes speed `n`, rather than speed `n^2`, the natural first scale for
order-`n` free-energy deviations.  Consistently, repeated use of (6) loses
only a constant fraction at each merge and constructs any recursively
attainable threshold with relative abundance `e^{-O(n)}`; all such events
have zero speed-`n^2` rate.

An exact support edge could still be represented by only a bounded number of
switching orbits and hence be speed-`n^2` rare.  But a speed-`n^2` LDP would
then need sharp control of the boundary of its effective domain at order-`n`
energy resolution.  Neither switching multiplicity nor (4) provides that
control.  At fixed `t`, the speed-`n^2` endpoint is already the universal
value `-\tfrac12\log\cosh t` by (12); on the project diagonal, (13) misses
the endpoint by a leading-order cross term and contracts `\beta`.

**Conclusion.**  The microcanonical product theorem (4) is exact and
nontrivial as a disorder-counting statement.  Its endpoint consequence is a
no-go: counting many annealed bridges does not repair the changing-temperature
obstruction, and a bare speed-`n^2` disorder LDP is not a substitute for the
missing fixed-`\beta` support-edge theorem.
