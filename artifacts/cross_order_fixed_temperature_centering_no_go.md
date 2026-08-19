# Fixed-edge-temperature centering cannot imply a sublinear own-scale defect

Status: **proved scalar-method no-go**.  The theorem below does not construct
complete signings and does not lower-bound the true cross-order defect.  It
shows sharply that fixed-edge-temperature centered subadditivity, even
together with the usual scalar pressure regularity, cannot imply an
`o(N)` own-scale composition defect.  A signing-specific input is necessary.

## 1. The exact implication one would need

For a complete signing put

```math
F_n(t)=\min_A\log\mathbb E_x\cosh(tH_A(x)),
\qquad e_n={n\choose2},
```

and

```math
R_n(t)=F_n(t)-e_n\log\cosh t.
```

Independent bridge averaging gives the exact fixed-`t` inequality

```math
R_{m+n}(t)\le R_m(t)+R_n(t).                         \tag{1.1}
```

For the own-scale pressure

```math
P_n(\beta)=F_n(\beta/\sqrt n),
```

write `N=m+n`, `t_N=beta/sqrt(N)`, and

```math
\Delta_k(\beta;N)
=F_k(\beta/\sqrt k)-F_k(t_N),\qquad k\in\{m,n\}.
```

Then (1.1) has exactly the following quantitative consequence:

```math
\boxed{
P_N(\beta)-P_m(\beta)-P_n(\beta)
\le mn\log\cosh t_N-\Delta_m(\beta;N)-\Delta_n(\beta;N).
}                                                        \tag{1.2}
```

Thus a scalar fixed-temperature proof of a defect `E_{m,n}=o(N)` must
derive

```math
\Delta_m+\Delta_n
\ge mn\log\cosh(\beta/\sqrt N)-o(N)                    \tag{1.3}
```

from additional information.  Convexity alone only gives

```math
\Delta_k(\beta;N)
\ge (1-\sqrt{k/N})P_k(\beta),                           \tag{1.4}
```

which need not have the magnitude in (1.3).

The next theorem proves more than the failure of this particular convexity
estimate: all the standard scalar axioms are compatible with a *positive
linear* own-scale defect at infinitely many balanced pairs.

## 2. Linear-defect scalar countermodel

### Theorem 2.1

There are functions `F_n^*:R -> [0,infinity)`, `n>=2`, with all of the
following properties.

1. `F_n^*` is even, real analytic and convex,
   `F_n^*(0)=0`, and `(F_n^*)''(0)=e_n`.
2. If

   ```math
   a_n=\lim_{t\to\infty}{F_n^*(t)\over t},
   ```

   then `a_n>0` is nondecreasing and

   ```math
   a_n|t|-n\log2\le F_n^*(t)\le a_n|t|.                \tag{2.1}
   ```

3. For each fixed `t>=0`, `F_n^*(t)` is nondecreasing in `n`.
4. The exactly centered functions

   ```math
   R_n^*(t)=F_n^*(t)-e_n\log\cosh t
   ```

   satisfy

   ```math
   R_{m+n}^*(t)\le R_m^*(t)+R_n^*(t)                   \tag{2.2}
   ```

   for all `m,n>=2`.
5. The diagonal family

   ```math
   \beta\longmapsto n^{-1}F_n^*(\beta/\sqrt n)
   ```

   is uniformly Lipschitz on `[0,infinity)`.
6. Its zero-temperature coefficients obey the same known asymptotic window
   as complete-signing minimizers: for all sufficiently large `n`,

   ```math
   0.34<{a_n\over n^{3/2}}<0.495.                      \tag{2.3}
   ```

7. Nevertheless, for every fixed `beta>0` there are infinitely many `n`
   for which

   ```math
   F_{2n}^*(\beta/\sqrt{2n})
   -2F_n^*(\beta/\sqrt n)
   \ge c_\beta n,                                      \tag{2.4}
   ```

   where `c_beta>0` depends only on `beta`.

Consequently no argument whose hypotheses use only properties 1--6 can
prove an `o(N)` own-scale composition defect.  In the notation requested
for the cross-order campaign,

```math
\text{properties 1--6}\ \not\Longrightarrow\ E_{n,n}=o(n).             \tag{2.5}
```

indeed they permit `E_{n,n}>=c_beta n` infinitely often.

### Construction

Let `h:R -> [-3/16,3/16]` be the continuous, `3`-periodic triangular wave.
Writing `v=u mod 3` in `[0,3)`, set

```math
h(u)=
\begin{cases}
 {3\over16}-{v\over4},&0\le v\le3/2,\\
 -{3\over16}+{v-3/2\over4},&3/2\le v<3.
\end{cases}                                             \tag{2.6}
```

Take `kappa=69/100` and define

```math
\vartheta_n
=\kappa e_n^{-1/4}\exp[-h(\log e_n)],
\qquad
L_n=e_n\vartheta_n^2,
\qquad
F_n^*(t)=L_n\log\cosh(t/\vartheta_n).                  \tag{2.7}
```

On the first half of each period, `h'(u)=-1/4`; on the second,
`h'(u)=1/4`.  Hence, away from the harmless corners,

```math
{d\over du}\log\vartheta(e^u)
=-{1\over4}-h'(u)\in\{0,-1/2\},                        \tag{2.8}
```

and

```math
{d\over du}\log L(e^u)
={1\over2}-2h'(u)\in\{1,0\}.                          \tag{2.9}
```

Thus `vartheta_n` is nonincreasing and `L_n` is nondecreasing.  Also
`a_n=L_n/vartheta_n=e_n vartheta_n` is nondecreasing, because its two
logarithmic slopes in `u=log e` are `1` and `1/2`.

The numerical choice of `kappa` makes the elementary pressure constraints
automatic.  Since `e_n>=1` and `|h|<=3/16`,

```math
0<\vartheta_n\le {69\over100}e^{3/16}<1,                \tag{2.10}
```

and, using `sqrt(e_n)<=n/sqrt2`,

```math
{L_n\over n}
\le {1\over\sqrt2}\left({69\over100}\right)^2
e^{3/8}<0.49.                                           \tag{2.11}
```

Formula (2.7) now gives analyticity, convexity,
`(F_n^*)''(0)=L_n/vartheta_n^2=e_n`, and zero-temperature slope
`a_n`.  The inequalities

```math
|z|-\log2\le\log\cosh z\le|z|
```

and (2.11) prove (2.1).  Since both `L_n` and `1/vartheta_n` are
nondecreasing, `F_n^*(t)` is nondecreasing in `n` for every fixed `t>=0`.

Moreover,

```math
{a_n\over n^{3/2}}
={69\over100}{e_n^{3/4}\over n^{3/2}}e^{-h(\log e_n)}.
```

Because `e_n^{3/4}/n^{3/2}->2^{-3/4}`, its asymptotic range is the
interval with endpoints

```math
{69\over100}2^{-3/4}e^{-3/16}
=0.340131\ldots,
\qquad
{69\over100}2^{-3/4}e^{3/16}
=0.494887\ldots.                                       \tag{2.12}
```

This proves property 6.  In particular, the countermodel cannot be excluded
by either side of the current `0.336493...` to `1/2` ground-state window.

### Proof of centered subadditivity

For `0<vartheta<=1`, put

```math
q_\vartheta(t)
=\vartheta^2\log\cosh(t/\vartheta)-\log\cosh t.
```

The elementary derivative calculation

```math
{\partial\over\partial\vartheta}
\{\vartheta^2\log\cosh(t/\vartheta)\}
=\vartheta\{2\log\cosh z-z\tanh z\}\ge0,
\qquad z=t/\vartheta,                                  \tag{2.13}
```

shows that `q_vartheta(t)` is nondecreasing in `vartheta`, while
`q_1(t)=0`.  Therefore `q_vartheta<=0`.

For completeness, the sign in (2.13) is global, not a small-`t`
approximation.  For `z>=0`, if

```math
k(z)=2\log\cosh z-z\tanh z,
```

then `k(0)=k'(0)=0` and

```math
k''(z)=2z\,\operatorname{sech}^2z\tanh z\ge0.
```

Thus `k>=0`; evenness covers negative `z`.

Since

```math
R_n^*(t)=e_nq_{\vartheta_n}(t)
```

and `e_{m+n}=e_m+e_n+mn`, monotonicity of `vartheta_n` gives

```math
\begin{aligned}
R_{m+n}^*(t)
&=(e_m+e_n+mn)q_{\vartheta_{m+n}}(t)\\
&\le(e_m+e_n)q_{\vartheta_{m+n}}(t)\\
&\le e_mq_{\vartheta_m}(t)+e_nq_{\vartheta_n}(t),
\end{aligned}                                           \tag{2.14}
```

which is (2.2).

### Proof of the linear balanced defect

Put `delta=3/4-log 2>0` and choose integers `n_j` with

```math
\log e_{n_j}=3j+\delta+o(1).                            \tag{2.15}
```

Such integers exist because the mesh of `log e_n` tends to zero.  Moreover,

```math
\log e_{2n_j}-\log e_{n_j}\longrightarrow\log4.
```

Since `delta+log4=3/4+log2<3/2`, both exact points lie, for all large
`j`, in the first half `[3j,3j+3/2]` of the same period.  Equation (2.8)
does not merely give asymptotic closeness: throughout that whole interval,

```math
\log\vartheta(e^u)
=\log\kappa-u/4-\{3/16-(u-3j)/4\}
=\log\kappa-3/16-3j/4
```

is exactly constant.  Therefore the two integer orders satisfy the exact
identity

```math
\vartheta_{2n_j}=\vartheta_{n_j}.                       \tag{2.16}
```

Furthermore

```math
d_j=\sqrt{n_j}\,\vartheta_{n_j}
\longrightarrow
d={69\over100}>0.                                      \tag{2.17}
```

Here the simplification follows from
`h(delta)=3/16-delta/4=(log2)/4`.  Using `e_n/n^2 -> 1/2`, (2.7),
(2.16), and (2.17) yields

```math
\begin{aligned}
{1\over n_j}\bigg[
F_{2n_j}^*\!\left({\beta\over\sqrt{2n_j}}\right)
-2F_{n_j}^*\!\left({\beta\over\sqrt{n_j}}\right)
\bigg]
\longrightarrow
d^2\left[
2\log\cosh\!\left({\beta\over\sqrt2d}\right)
-\log\cosh\!\left({\beta\over d}\right)
\right].                                               \tag{2.18}
\end{aligned}
```

The bracket is strictly positive for every `beta>0`.  Indeed, with
`z=beta/d`,

```math
g(z)=2\log\cosh(z/\sqrt2)-\log\cosh z
```

satisfies `g(0)=0` and

```math
g'(z)=\sqrt2\tanh(z/\sqrt2)-\tanh z>0                 \tag{2.19}
```

because `tanh z/z` is strictly decreasing on `(0,infinity)`.  Taking one
half of the positive limit in (2.18) as `c_beta` proves (2.4).

Finally, the diagonal derivative is

```math
{d\over d\beta}\left[{1\over n}F_n^*(\beta/\sqrt n)\right]
={a_n\over n^{3/2}}
\tanh\!\left({\beta\over\sqrt n\,\vartheta_n}\right). \tag{2.20}
```

The prefactor is uniformly bounded because

```math
{a_n\over n^{3/2}}
=\kappa {e_n^{3/4}\over n^{3/2}}e^{-h(\log e_n)}
\le {69\over100}2^{-3/4}e^{3/16}.
```

This proves property 5 and completes the theorem.

## 3. Exact scope of the obstruction

The no-go covers deductions from:

- fixed-edge-temperature centered subadditivity;
- convexity, analyticity, evenness and the exact second derivative at zero;
- the normalized-partition entropy squeeze;
- restriction-type monotonicity in the order;
- monotone zero-temperature slope; and
- uniform diagonal Lipschitz regularity.

It does **not** cover a proof using coefficient-level Eulerian identities,
actual Gibbs overlaps, bridge-dependent cancellation, optimizer-specific
stationarity beyond scalar pressure, or another signing-specific constraint.
Those are exactly the kinds of data that could distinguish genuine complete
quadratic signing pressures from the countermodel.

The quantitative campaign conclusion is therefore:

```math
\boxed{
\text{fixed-}t\text{ centered subadditivity + universal scalar regularity}
\not\Longrightarrow E_{n,n}=o(n).
}                                                        \tag{3.1}
```

This is a method-class impossibility, not Level 6 and not a new bound on the
true defect for actual optimizing children.
