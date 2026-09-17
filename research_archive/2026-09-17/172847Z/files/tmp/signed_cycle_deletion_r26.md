# Wave 26 Route 3: cancellation-preserving signed-cycle deletion

## Status

The multivariate cycle-polynomial identities, Boolean-noise representation,
star-coupling susceptibility formula, Jensen moment bound, sharp `9/16`
sufficient lemma, and deficit-partition recurrence below are **Verified**.
They are checked on `A_5,A_6,A_8,A_9` by
`tmp/signed_cycle_deletion_r26.py`.

The required minimizer-specific cavity second-moment lower bound is an
**Open target**.  `A_9` is a finite exact-minimizer obstruction to treating
that moment bound as automatic, not an asymptotic falsifier.  The final
deficit-ratio condition is a **Verified sufficient falsification criterion**;
no complete-signing family satisfying it asymptotically is presently known.
No convergence theorem or asymptotic signing counterexample is obtained.

## 1. A positive multivariate realization that retains the signed cycles

Fix an order-`r` complete signing `B`, put
`\rho=\tanh(2\beta)`, and let `\mathcal E_{\rm even}` be the even-cardinality
Eulerian edge sets from (10.778).  Introduce vertex variables by

```math
\boxed{
\mathscr P_B(z_1,\ldots,z_r)
=\sum_{F\in\mathcal E_{\rm even}}
b_F\rho^{|F|}\prod_{v=1}^r z_v^{\deg_F(v)}.}
\tag{R26.1}
```

This does not take absolute values of the cycle coefficients.  In fact it
has an exact positive primal realization on the whole cube
`0\le z_v\le1`.  For each edge set

```math
K_{uv}(z)=\operatorname{arctanh}(\rho z_uz_v).
```

Then, on the `2^r` oriented-projective states,

```math
\boxed{
\mathscr P_B(z)
=\frac{
2^{-r}\sum_{(\sigma,x)}
\exp\{\sum_{u<v}K_{uv}(z)\sigma b_{uv}x_ux_v\}}
{\prod_{u<v}\cosh K_{uv}(z)}>0.}
\tag{R26.2}
```

Expanding each edge factor proves (R26.2): the orientation average removes
odd edge cardinality and the projective spin average removes odd vertex
degrees, with the original signed coefficient `b_F` unchanged.

Set all variables except `z_i` equal to one and write

```math
\mathscr P_i(t)=\mathscr P_B(1,\ldots,1,t,1,\ldots,1).
```

Then

```math
\mathscr P_i(0)=\mathcal P_{B[-i]}(\rho),
\qquad
\mathscr P_i(1)=\mathcal P_B(\rho).
\tag{R26.3}
```

There is also an exact cancellation-preserving use of every star competitor.
If `S` contains each edge of `\delta(i)` independently with probability
`(1-t)/2`, Walsh character averaging gives

```math
\boxed{
\mathscr P_i(t)
=\mathbb E_S\mathcal P_{B^S}(\rho).}
\tag{R26.4}
```

Thus this interpolation is not a coefficientwise majorant: it is the
Boolean-noise average of the actual signed competitor partition functions.
For an exact minimizer, every competitor floor `D_S\ge1` gives a lower bound
on (R26.4).  It still gives no upper bound on the deletion ratio.  The
finite audit below finds strongly negative initial log curvature, including
on `A_9`, so exact minimality supplies no missing convexity sign.

## 2. Exact star-coupling and susceptibility identities

Let `C_i=B[-i]`, let `b_i` be the deleted row, and let
`\mu_{\beta,C_i}` be the oriented child Gibbs law.  Put

```math
L_i(y)=b_i^{\mathsf T}y,
\qquad
F_i(K)=\mathbb E_{\mu_{\beta,C_i}}\cosh(KL_i).
\tag{R26.5}
```

If

```math
K(t)=\operatorname{arctanh}(\rho t),
```

then direct integration of the deleted spin in (R26.2) gives

```math
\boxed{
\frac{\mathscr P_i(t)}{\mathcal P_{C_i}(\rho)}
=\frac{F_i(K(t))}{\cosh(K(t))^{r-1}}.}
\tag{R26.6}
```

At `t=1`, this is exactly (10.779), since

```math
e^{\beta\kappa_{\beta,i}(B)}=F_i(2\beta).
\tag{R26.7}
```

To expose the cancellation dynamically, extend the child law by a symmetric
spin and let `H_i` be its signed star field.  Tilt this extension by
`e^{KH_i}` and denote expectation by `\mathbb E_K`.  For
`f_i(K)=\log F_i(K)`, ordinary log-partition differentiation gives

```math
f_i'(K)=\mathbb E_KH_i,
\qquad
f_i''(K)=\operatorname{Var}_K(H_i),
\qquad
f_i(0)=f_i'(0)=0.
```

Consequently

```math
\boxed{
\begin{aligned}
\beta\kappa_{\beta,i}(B)
&=\int_0^{2\beta}\mathbb E_KH_i\,dK\\
&=\int_0^{2\beta}(2\beta-K)
\operatorname{Var}_K(H_i)\,dK.
\end{aligned}}
\tag{R26.8}
```

This is an exact star susceptibility, not a bound on unsigned cycle mass.
Differentiating (R26.6) at the deleted endpoint gives another useful exact
identity.  With

```math
v_i(\beta;B)
:=\mathbb E_{\mu_{\beta,C_i}}L_i^2
=f_i''(0),
```

one has

```math
\boxed{
\left.\frac{d^2}{dt^2}\log\mathscr P_i(t)\right|_{t=0}
=\rho^2\,[v_i(\beta;B)-(r-1)].}
\tag{R26.9}
```

Thus the cavity second moment is precisely the signed-cycle deletion
curvature relative to the uniform row value `r-1`.

## 3. The sharp cavity second-moment sufficient lemma

The function `x\mapsto\cosh(a\sqrt x)` is convex on
`[0,\infty)`: its power series is a positive constant plus a linear term and
positive multiples of `x^k`, `k\ge2`.  Jensen applied to (R26.5) therefore
gives

```math
\boxed{
\kappa_{\beta,i}(B)
\ge\frac1\beta
\log\cosh\left(2\beta\sqrt{v_i(\beta;B)}\right)
\ge2\sqrt{v_i(\beta;B)}-\frac{\log2}{\beta}.}
\tag{R26.10}
```

This yields a concrete theorem with the exact constant needed by the
constant-shortfall route.  Freeze a root scale

```math
\alpha=\frac{q_n}{n^{3/2}}
```

along an order window.  At a current order `r`, if some selectable vertex
satisfies

```math
\boxed{
v_i(\beta;B_r)\ge\frac9{16}\alpha^2r,}
\tag{R26.11}
```

then

```math
\boxed{
\kappa_{\beta,i}(B_r)
\ge
\alpha[r^{3/2}-(r-1)^{3/2}]
-\frac{\log2}{\beta}.}
\tag{R26.12}
```

Indeed, (R26.10) gives
`\kappa\ge(3/2)\alpha\sqrt r-(\log2)/\beta`, while convexity of
`x^{3/2}` gives

```math
r^{3/2}-(r-1)^{3/2}\le\frac32\sqrt r.
```

No positive margin beyond `9/16` is needed.  Therefore (R26.11) at every
selected step on the required paths proves (10.617) with the universal
shortfall

```math
K=\frac{\log2}{\beta}.
```

This is a genuine new sufficient lemma, but (R26.11) is still open for
exact minimizers.  It asks for one child Gibbs law whose deleted row has a
linear second moment.  It does not use the circular restriction-loss mean.

## 4. Exact finite obstruction and why competitor floors do not fill the gap

The moment premise is not automatic at finite order.  The relevant exact
values at `\beta=1` are:

| minimizer | `r` | `q` | range of `v_i` | range of actual `\kappa_i` | range of Jensen bound |
|:---|---:|---:|---:|---:|---:|
| `A_5` | 5 | 8 | `0.146476` | `0.880058` | `0.268146` |
| `A_6` | 6 | 10 | `1.004292` | `1.518841` | `1.329137` |
| `A_8` | 8 | 20 | `1.189376` | `2.261149` | `1.500691` |
| `A_9` | 9 | 24 | `0.149208--0.292952` | `1.094479--1.760116` | `0.272738--0.497983` |

For `A_9`, `\alpha=8/9`, so the right side of (R26.11) is exactly `4`.
Every coordinate misses it by more than an order of magnitude at
`\beta=1`; at `\beta=1/2` the range is still only
`1.226168--2.126418`.  This **falsifies** any claim that the sharp moment
premise follows pointwise from exact minimality at all finite orders.  It
does not rule out an asymptotic selected-root theorem, an additive
finite-order allowance, or another susceptibility functional.

The corresponding initial signed-cycle curvatures from (R26.9) are all
negative.  At `\beta=1` they range from `-3.5813` on `A_5` to
`[-7.2961,-7.1625]` on `A_9`.  Thus the positive multivariate polynomial is
not log-convex along the star-deletion interpolation.  The exact competitor
constraints in (R26.4) only put a floor under its noise averages, reproducing
the already known wrong-way direction rather than proving (R26.11).

## 5. A complete-signing falsification criterion

There is a second exact deletion recurrence which cleanly separates the
zero-temperature decrement and the finite-temperature degeneracy.  Let

```math
d_i=Q(B)-Q(C_i),
\qquad
D_\beta(B)=\sum_\omega e^{-\beta[Q(B)-e_B(\omega)]}.
```

Rearranging the exact cavity affinity (10.623) gives

```math
\boxed{
\kappa_{\beta,i}(B)
=d_i+rac1\beta
\log\frac{D_\beta(B)}{2D_\beta(C_i)}.}
\tag{R26.13}
```

This retains the full positive deficit partition functions, and hence all
signed-cycle cancellation through their primal values.  It neither takes
absolute cycle coefficients nor invokes a Parseval budget.

It gives the following testable minimizer-specific falsifier.  Suppose a
fixed `\beta>0` and an unbounded sequence of exact complete minimizers
`B_r` satisfy

```math
\boxed{
\max_i d_i=o(\sqrt r),
\qquad
\max_i\left[
\log\frac{D_\beta(B_r)}{2D_\beta(B_r[-i])}
\right]_+=o(\sqrt r).}
\tag{R26.14}
```

Then (R26.13) proves

```math
\max_i\kappa_{\beta,i}(B_r)=o(\sqrt r),
```

which is a literal fixed-temperature falsification of (10.617) on those
roots.  To falsify the convergence-weakened path version, the same condition
must hold for every eligible root or across every required deletion cutset,
as recorded in `STEERING.md`.

`A_9` is a finite template for (R26.14): all `d_i=0`, and at `\beta=1` the
second term ranges only from `1.094479` to `1.760116`.  No unbounded exact
minimizer family with (R26.14) is known, so this is a criterion rather than an
asymptotic counterexample.

## 6. Checker coverage

For each of `A_5,A_6,A_8,A_9` and
`\beta\in\{1/4,1/2,1\}`, the checker independently verifies:

1. the signed-cycle ratio (10.779) from positive primal partition sums;
2. the multivariate identity (R26.6) at four interpolation points;
3. the Boolean-noise competitor identity (R26.4) by enumerating every star
   perturbation;
4. both susceptibility integrals in (R26.8) by numerical quadrature;
5. the curvature identity (R26.9) and Jensen bound (R26.10); and
6. the exact deficit recurrence (R26.13).

The identities and inequalities are proved algebraically; quadrature and the
displayed decimals are **Numerical audits** of those proofs.

## Route conclusion

The signed-cycle deletion ratio admits a clean cancellation-preserving
interpolation.  Its full logarithmic derivative is a star magnetization and
its second derivative is a genuine susceptibility.  This converts the
square-root cavity problem into the concrete moment target (R26.11), with the
sharp sufficient coefficient `9/16` and constant shortfall
`(\log2)/\beta`.

Exact minimality does not by itself force that finite-order moment: `A_9`
and the negative deletion curvatures rule out the naive version.  The two
remaining precise outcomes are therefore:

- prove a minimizer-specific, path-selectable version of (R26.11), possibly
  using more than the pointwise competitor floors; or
- construct an unbounded exact-minimizer family satisfying (R26.14), which
  would falsify the literal fixed-temperature constant-shortfall route.
