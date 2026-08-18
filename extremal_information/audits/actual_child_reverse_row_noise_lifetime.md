# Audit: row-noise lifetime representation of the actual-child negative tilt

Status: **rigorous exact identity and two sharp benchmark tests; no
actual-child bound**.  A raw coordinate Dirichlet form pays a Walsh
interaction once for every coordinate in its support.  The row-refresh noise
flow below instead pays that interaction over its lifetime: its instantaneous
multiplicity is cancelled by its proportionally shorter survival time.

Applied relative to the existing canonical row product, this gives an exact
upper certificate for the reverse row-product projection without invoking
the unknown optimal product.  It vanishes on the fixed-projective
row-factor sharpness channel and gives the correct `Theta(N^(-1))` cost for
the pure-parity ceiling.  Its limitation is equally precise: the lifetime
integral equals the canonical error `J`, so controlling it remains stronger
than controlling the best reverse projection and does not address coherent
factor retuning.

The evidentiary split in this note is:

1. **Exact:** (NR.5), (NR.6), and the specialization
   `I^leftarrow <= J` for the actual child-induced law.
2. **Exact benchmark calculations:** the fixed-projective and parity tests.
3. **Perturbative only:** the Hoeffding inverse-order expansion (NR.13)--
   (NR.14), for a fixed finite product space as `epsilon -> 0`.
4. **Not proved here:** that the lifetime representation is easier to bound
   for optimizing children, that `J=o(N)`, or any recurrence consequence.

## 1. Product-base row-refresh flow

Let

```math
\Omega=\prod_{i=1}^m\Omega_i,
\qquad r=\bigotimes_{i=1}^m r_i                       \tag{NR.1}
```

be any full-support row product on a finite product space.  Let `q` be a
full-support law, and write

```math
h={dq\over dr},\qquad E_rh=1.                        \tag{NR.2}
```

For each row, let `E_i` average that row against `r_i`, leaving the other
rows fixed.  Define

```math
\mathcal L_r=\sum_{i=1}^m(I-E_i),
\qquad
P_t^r=e^{-t\mathcal L_r},
\qquad
h_t=P_t^rh.                                          \tag{NR.3}
```

Operationally, `P_t^r` independently retains each row with probability
`e^(-t)` and otherwise refreshes it from `r_i`.

For a positive density `g`, define its reverse row-refresh Fisher
functional

```math
\mathfrak K_r(g)
=\sum_{i=1}^m E_{r_{-i}}
 \left[(E_ig)(E_i g^{-1})-1\right].                 \tag{NR.4}
```

Conditional Cauchy--Schwarz makes every summand nonnegative.  Equivalently,
if `X_i,X_i'` are conditionally independent samples from `r_i`, then

```math
\mathfrak K_r(g)
={1\over2}\sum_iE_rE_{X_i'}
 {\{g(X)-g(X_{-i},X_i')\}^2
  \over g(X)g(X_{-i},X_i')}.                        \tag{NR.5}
```

## 2. Exact reverse de Bruijn identity

### Theorem NR.1 (reverse density is paid by lifetime)

For every finite `q` and full-support row product `r`,

```math
\boxed{
D(r\Vert q)=\int_0^\infty\mathfrak K_r(P_t^rh)\,dt.} \tag{NR.6}
```

Consequently, if

```math
\mathcal I^{\leftarrow}(q)
=\inf_{p=\otimes_i p_i}D(p\Vert q),                 \tag{NR.7}
```

then the optimizer-free inequality

```math
\boxed{
\mathcal I^{\leftarrow}(q)
\le D(r\Vert q)
=\int_0^\infty\mathfrak K_r(P_t^rh)\,dt}            \tag{NR.8}
```

holds for every declared row product `r`.

*Proof.*  Put

```math
R(t)=-E_r\log h_t=D(r\Vert h_t r).
```

Since `partial_t h_t=-mathcal L_rh_t`,

```math
\begin{aligned}
R'(t)
 &=-E_r{\partial_t h_t\over h_t}\\
 &=\sum_i\left\{1-E_r{E_ih_t\over h_t}\right\}\\
 &=-\sum_iE_{r_{-i}}
   \left[(E_ih_t)(E_i h_t^{-1})-1\right]\\
 &=-\mathfrak K_r(h_t).                              \tag{NR.9}
\end{aligned}
```

The row-refresh flow is ergodic, so `h_t->E_rh=1` and `R(t)->0`.
Integration of (NR.9) proves (NR.6).  Since `r` itself is admissible in
(NR.7), (NR.8) follows. `square`

In the actual-child notation, take `r` to be the existing canonical inverse
row product and `q` the full negative-disorder escort.  Then

```math
D(r\Vert q)=\mathcal J,                              \tag{NR.10}
```

so (NR.6) is an exact multiscale representation of the canonical
negative-tilt product certificate.  It does not assume or reconstruct the
optimal reverse product `p^*`.

## 3. Walsh/Hoeffding lifetime cancellation

The sense in which (NR.6) removes raw multiplicity is already visible at
quadratic order.  Let

```math
h=1+\varepsilon g,
\qquad
g=\sum_{\varnothing\ne S\subseteq[m]}g_S           \tag{NR.11}
```

be the orthogonal Hoeffding decomposition under `r`, where `g_S` genuinely
depends on every row in `S`.  Then

```math
P_t^rg_S=e^{-|S|t}g_S,                               \tag{NR.12}
```

and finite Taylor expansion gives

```math
\mathfrak K_r(P_t^rh)
=\varepsilon^2\sum_{S\ne\varnothing}
 |S|e^{-2|S|t}\|g_S\|_{L^2(r)}^2
+O_{r,g}(\varepsilon^3e^{-3t}).                    \tag{NR.13}
```

Here and below the expansion is for a fixed finite product space and fixed
bounded `g`, with `epsilon` small enough that `1+epsilon*g>0`.  The remainder
is integrable in `t`: every nonconstant Hoeffding component of `P_t^r g`
decays at least as `e^(-t)`.

Therefore

```math
\int_0^\infty\mathfrak K_r(P_t^rh)dt
={\varepsilon^2\over2}\sum_{S\ne\varnothing}
 \|g_S\|_2^2+O(\varepsilon^3).                      \tag{NR.14}
```

The instantaneous factor `|S|` is cancelled by lifetime `1/(2|S|)`.
This agrees with
`D(r||(1+epsilon g)r)=epsilon^2||g||_2^2/2+O(epsilon^3)`.
The identity (NR.6), unlike this expansion, is exact at arbitrary amplitude.

## 4. Benchmark A: fixed-projective row factor

For the sharpness channel

```math
Q_{ij}=X_i z_j,
```

the forward likelihood is exactly

```math
p(B)=\prod_i p_i(B_i).
```

Its full negative escort is therefore

```math
q(B)\propto p(B)^{-\lambda}
=\prod_i p_i(B_i)^{-\lambda}.                       \tag{NR.15}
```

The canonical inverse row product is exactly this same law: `r=q`.  Hence

```math
\boxed{
h=1,
\qquad
\mathfrak K_r(P_t^rh)=0\ \ (t\ge0),
\qquad
\mathcal J=\mathcal I^{\leftarrow}=0.}              \tag{NR.16}
```

This holds although the raw negative-tilt cavity overlap converges to a
strictly positive constant.  Thus the row-base lifetime flow removes the
row-local false positive exactly; the fair-base reverse entropy does not.

## 5. Benchmark B: pure parity

Let `U` be fair on `d` bits, let `chi` be their full parity, and put

```math
{dq\over dU}={e^{-a\chi}\over\cosh a}
=1-\vartheta\chi,
\qquad \vartheta=\tanh a.                           \tag{NR.17}
```

Here the natural product base is `r=U`.  Under ordinary bit-noise with
Walsh eigenvalue `e^(-t)` per bit,

```math
h_t=1-\vartheta e^{-dt}\chi.                         \tag{NR.18}
```

Every coordinate contributes the same amount, so

```math
\mathfrak K_U(h_t)
={d\vartheta^2e^{-2dt}
  \over1-\vartheta^2e^{-2dt}}.                       \tag{NR.19}
```

Consequently

```math
\boxed{
\int_0^\infty\mathfrak K_U(h_t)dt
=-{1\over2}\log(1-\vartheta^2)
=\log\cosh a
=D(U\Vert q).}                                      \tag{NR.20}
```

The row-refresh version is identical with `d` replaced by the number of
rows met by the parity, and gives the same integral.  In contrast, the raw
time-zero bit Fisher cost is

```math
\mathfrak K_U(h_0)=d\sinh^2a.                        \tag{NR.21}
```

For the physical parity ceiling `a=lambda beta/sqrt(N)` and
`d=Theta(N^2)`, (NR.20) is `Theta(N^(-1))`, whereas (NR.21) is
`Theta(N)`.  Lifetime renormalization removes exactly the spurious
quadratic coordinate multiplicity identified by the parity falsifier.

### Theorem NR.2 (the late row-noise tail is uniformly finite-context)

The refresh representation also gives an exact depth truncation.  With
probability `(1-e^(-t))^m`, every row has been refreshed, and that component
of `h_t` is exactly `E_rh=1`.  Positivity of every other component gives

```math
h_t\ge(1-e^{-t})^m.
```

Using (NR.6) from time `T` onward therefore yields

```math
\boxed{
\int_T^\infty\mathfrak K_r(h_t)dt
=D(r\Vert h_Tr)
\le-m\log(1-e^{-T}).}                              \tag{NR.21a}
```

In particular, after `T=log m+u` the omitted tail is at most
`e^(-u)/(1-e^(-u)/m)`.  Hence only a logarithmic row-noise horizon is ever
needed for additive constant accuracy.  This does **not** compress the early
part of the path: for `t=O(1)`, `h_t` can still retain marginals involving a
linear number of rows.

The latter limitation is exact rather than heuristic.  If `p=e^{-t}`, the
product refresh construction expands as

```math
h_t(x)=\sum_{S\subseteq[m]}p^{|S|}(1-p)^{m-|S|}
             E_r[h\mid X_S=x_S].                  \tag{NR.21b}
```

Thus at every fixed positive time the represented subset has typical size
`pm=Theta(m)`.  Row noise replaces one full density by an average of its
linear-order rooted marginals; it does not reduce the exact state to bounded
or sublinear row order.

### Corollary NR.3 (extensive lifetime has a bounded-time witness)

Suppose `m<=N`, `0<eta<=1`, and

```math
D(r\Vert q)\ge\eta N.
```

Set `T_eta=log(4/eta)`.  Since
`-log(1-eta/4)<=eta/3`, (NR.21a) gives

```math
\int_0^{T_\eta}\mathfrak K_r(h_t)dt
\ge {2\eta N\over3}.
```

Consequently some `t in [0,T_eta]` satisfies

```math
\boxed{
\mathfrak K_r(h_t)\ge {2\eta N\over3\log(4/\eta)}.} \tag{NR.21c}
```

Hence an extensive canonical product error cannot hide exclusively at a
noise time diverging with system size: it has an extensive witness in a
fixed window depending only on its density.  Equation (NR.21b) also shows
why this is not yet a low-information certificate.  At every time in that
fixed window the exact witness can still use `Theta(N)`-row marginals.

## 6. What this can replace, and what it cannot

The fair-base logarithmic-Sobolev step in Theorem 37.50 gives

```math
\mathcal I^{\leftarrow}
\le D(U\Vert q)
```

and then bounds the right side by raw coordinate cavity mass.  Equation
(NR.8), with the canonical row base, gives the different structural
certificate

```math
\boxed{
\mathcal I^{\leftarrow}
\le\mathcal J
=\int_0^\infty\mathfrak K_r(P_t^r(dq/dr))dt.}        \tag{NR.22}
```

It requires neither raw overlap decay nor the unknown product optimizer.
The sole operational conclusion supplied by this note is the first
inequality in (NR.22): any independently proved upper bound on `J` is also an
upper bound on `I^leftarrow`.

This is not yet a Level-6 estimate.  The smoothed density `P_t^r(dq/dr)`
still contains the exact actual-child interaction, and (NR.22) equals the
canonical error `J`; it may overestimate the best reverse projection by the
entire coherent-retuning gap `J-I^leftarrow`.  Thus the identity supplies a
possible rare-event/renormalization diagnostic, not its optimizer-specific
bound.  The unresolved actual-child question exposed by the identity is

```math
\int_0^\infty\mathfrak K_r(P_t^r(dq/dr))dt=o(N),     \tag{NR.23}
```

or whether a fixed part of this lifetime mass survives after every canonical
row factor.  Neither alternative is established here, and (NR.23) is not
claimed to be operationally simpler than controlling the original
actual-child interaction.

For a noncanonical declared base this limitation is unavoidable even on
products: if both `q` and `r` are row products but `q!=r`, then
`I^leftarrow(q)=0` while the integral equals
`D(r||q)=sum_iD(r_i||q_i)`, which may be linear.  The canonical base avoids
this particular false positive because it equals `q` whenever `q` is row
product, but no theorem here shows that it removes an approximately
row-product background at sublinear error.
