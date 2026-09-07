# Weighted endpoint harvesting and the exact grouped descent threshold

Scratch memo for Wave 11.  Nothing here uses the false stopping inequality
`(10.435)`.

## 1. A genuinely new weighted consequence of `K <= 4`

Fix an endpoint partition tree for a matrix `A`.  For a directed bucket
`b=(v,X,sigma)`, let

```math
c_b=(\mu_{vX}^{\sigma})_+,
\qquad 0\le a_b\le c_b
```

be the capacity and the allocation supplied by `(10.473)`.  Let `e(b)` be
the underlying unoriented tree edge.  The path-cover theorem supplies a
finite measure `nu` on oriented root-to-leaf chains such that

```math
|\nu|=K\le4t_0\le4,
\qquad
\nu\{\pi:b\in\pi\}\ge \frac{a_b}{c_b}
```

for every positive-capacity bucket (with the usual zero convention).
Consequently, for **every** family of nonnegative weights `omega_b`,

```math
\boxed{
\sum_b\omega_ba_b
\le
\int\sum_{b\in\pi}\omega_bc_b\,d\nu(\pi)
\le
K\sup_\pi\sum_{b\in\pi}\omega_bc_b
\le4\sup_\pi\sum_{b\in\pi}\omega_bc_b.
}
\tag{W1}
```

This is stronger than merely applying `(10.341)` with all weights one.
It follows directly from coordinatewise domination, so no compatibility or
purification assumption is hidden in it.

Suppose the weights depend only on the underlying tree edge.  Along a chain,
write the split as `U_t=D_t sqcup X_t`, with `X_t=U_{t+1}`, and put

```math
d_t=Q(U_t)-Q(X_t)\ge0.
```

The selected capacity is bounded by the all-successor layer, and the local
step in the proof of `(10.341)` gives

```math
c_b\le \mathcal L_t\le d_t+Q(D_t).
```

Thus `(W1)` gives the exact weighted harvest

```math
\boxed{
\sum_b\omega_{e(b)}a_b
\le
K\sup_\pi
\sum_{t\in\pi}\omega_t\bigl(d_t+Q(D_t)\bigr).
}
\tag{W2}
```

Equation `(W2)`, rather than the unweighted `12Q`, is the strongest formal
weighted conclusion of the present endpoint certificate.

### A scale-local corollary

Fix `r>=1` and `s>=1`, and retain only edges for which

```math
|U_t|\ge r,
\qquad |D_t|\le s.
```

Give such an edge weight `|U_t|^{-3/2}` and every other edge weight zero.
Along one chain the `D_t` are disjoint.  Hence

```math
\sum_{\rm selected}d_t\le Q(A),
\qquad
\sum_{\rm selected}Q(D_t)
\le\sum_t|D_t|(|D_t|-1)
\le(s-1)n.
```

Therefore

```math
\boxed{
\sum_{b:\ |U_b|\ge r,\ |D_b|\le s}
\frac{a_b}{|U_b|^{3/2}}
\le
\frac K{r^{3/2}}
\bigl(Q(A)+(s-1)n\bigr).
}
\tag{W3}
```

For the dyadic block class `s<=|D_t|<2s`, the last term is instead
`(2s-2)n`.  In particular, when `r=rho*n` and `s=o(sqrt(n))`, the internal
block part of `(W3)` is `o(1)` after normalization.  What `(W3)` does **not**
say is that a definite portion of the endpoint allocation lies in this
size/order window.

There is also a useful intrinsic weighting.  Retain the initial chain edges
whose child and parent norms are both at least `q_*>0`, and take
`omega_t=1/Q(U_t)`.  Then

```math
\sum_t\frac{d_t}{Q(U_t)}
=\sum_t\left(1-\frac{Q(U_{t+1})}{Q(U_t)}\right)
\le\log\frac{Q(A)}{q_*}.
```

Thus

```math
\boxed{
\sum_b\frac{a_b}{Q(U_b)}
\le K\sup_\pi\left[
\log\frac{Q(A)}{q_*}
+\sum_t\frac{Q(D_t)}{Q(U_t)}
\right]
}
\tag{W4}
```

for buckets above the stopping threshold.  This is a logarithmic weighted
Carleson estimate, but again it applies to the endpoint allocation, not to a
prescribed temporal replenishment process.

## 2. Exact fixed-cardinality grouped identity

Let `A` be an exact order-`n` minimizer in doubled normalization,

```math
Q(A)=q_n.
```

Orient and switch an absolute ground to `1`, put `r=A1`, and choose a uniform
`h`-set `H`.  Write `T=H^c`, `m=n-h`, and

```math
e_H=1_T^TA[T]1_T,
\qquad d_H=q_n-Q(A[T]),
\qquad g_H=Q(A[T])-e_H,
\qquad \varepsilon_H=Q(A[T])-q_m.
```

The one-block identity `(10.328)` is

```math
d_H+g_H=2R_H-h_H,
```

where `R_H=sum_{i in H}r_i` and `h_H=1_H^TA[H]1_H`.  Uniform averaging gives

```math
\mathbb ER_H=\frac hnq_n,
\qquad
\mathbb Eh_H=\frac{h(h-1)}{n(n-1)}q_n.
```

Define

```math
\alpha_{n,h}
=\frac{2h}{n}-\frac{h(h-1)}{n(n-1)}.
```

Then

```math
\boxed{
q_n-q_m
=\alpha_{n,h}q_n-\mathbb E(g_H-\varepsilon_H).
}
\tag{G1}
```

There is an important audit simplification.  Pointwise,

```math
g_H-\varepsilon_H=q_m-e_H.
```

If

```math
\beta_{n,h}=\frac{m(m-1)}{n(n-1)},
```

then `E e_H=beta_{n,h}q_n` and `alpha_{n,h}+beta_{n,h}=1`.  Hence

```math
\boxed{
\mathbb E(g_H-\varepsilon_H)
=q_m-\beta_{n,h}q_n.
}
\tag{G2}
```

Thus `(G1)` is an exact bookkeeping identity, not by itself a new
replenishment theorem.

Put `rho=m/n` and `p=h/n`.  The exact coefficient that would preserve the
`3/2` scale is

```math
\begin{aligned}
\Delta_{n,h}
&=\alpha_{n,h}-\bigl(1-\rho^{3/2}\bigr)\\
&=\rho^{3/2}-\rho^2+\frac{p\rho}{n-1}\\
&=\rho^{3/2}(1-\sqrt\rho)+\frac{p\rho}{n-1}.
\end{aligned}
\tag{G3}
```

Equations `(G1)`--`(G3)` show the equivalence

```math
\boxed{
\mathbb E(g_H-\varepsilon_H)
\le\Delta_{n,h}q_n+E_{n,h}
\quad\Longleftrightarrow\quad
\frac{q_m}{m^{3/2}}
\le\frac{q_n}{n^{3/2}}+\frac{E_{n,h}}{m^{3/2}}.
}
\tag{G4}
```

Therefore the excess-adjusted `Delta` inequality is the logically weakest
scalar condition, but it is exactly normalized monotonicity in disguise.
It must not be advertised as a consequence of the chain cover.

## 3. Exact group-size constants

Ignoring only the displayed finite-population correction in `(G3)`, put

```math
\Delta(p)=(1-p)^{3/2}-(1-p)^2.
```

Then:

```math
\Delta(p)
=\frac p2-\frac{5p^2}{8}+\frac{p^3}{16}+O(p^4).
\tag{S1}
```

Consequently:

* `h=1`: `Delta*q_n=[(1-1/n)^{3/2}-1+2/n]q_n
  =(1/(2n)+3/(8n^2)+O(n^{-3}))q_n`;
* `h=o(n)`: the leading permissible net replenishment is
  `(h/(2n))q_n`, with the exact corrections in `(S1)` and `(G3)`;
* `h=lambda*sqrt(n)`: for `q_n=Cn^{3/2}`, the leading budget is
  `(C lambda/2)n`, while the next term is
  `-(5C lambda^2/8)sqrt(n)`;
* `h=n/4`: `Delta=(6sqrt(3)-9)/16`, plus `3/[16(n-1)]`;
* `h=n/2`: `Delta=(sqrt(2)-1)/4`, plus `1/[4(n-1)]`;
* `h=3n/4`: `Delta=1/16`, plus `3/[16(n-1)]`.

The macroscopic coefficient is maximized when `sqrt(rho)=3/4`, namely

```math
\boxed{
rho=\frac9{16},\qquad p=\frac7{16},\qquad
\max\Delta(p)=\frac{27}{256}.
}
\tag{S2}
```

Thus even the best macroscopic group permits only about `0.10547 q_n` of
net replenishment.

## 4. Bernoulli and field-weighted groups

For independent deletion indicators with arbitrary probabilities `p_i`, put

```math
a=p^Tr,
\qquad c=p^TAp,
\qquad K=|H|,
\qquad m=n-K.
```

Equation `(10.329)` and terminal excess give the exact identity

```math
\boxed{
q_n-\mathbb E q_m
=2a-c-\mathbb E(g_H-\varepsilon_H).
}
\tag{B1}
```

The exact scale-preserving budget is

```math
\boxed{
\mathfrak B(A,p)
=2a-c-q_n\,\mathbb E\left[1-(m/n)^{3/2}\right].
}
\tag{B2}
```

As in the fixed-size case,

```math
\mathbb E(g_H-\varepsilon_H)\le\mathfrak B(A,p)+E
```

is equivalent to

```math
\mathbb E q_m
\le q_n\mathbb E(m/n)^{3/2}+E.
```

For the field-proportional choice

```math
p_i=\frac{\lambda r_i}{n-1}
```

(or its restriction to a heavy set),

```math
a=\frac\lambda{n-1}\sum_i r_i^2,
\qquad
c=\frac{\lambda^2}{(n-1)^2}r^TAr,
\qquad
\mathbb EK=\frac{\lambda q_n}{n-1}.
\tag{B3}
```

If the row-square mass is at its Cauchy minimum and `c` is lower order,
then `(B2)` has leading value

```math
\frac{\lambda q_n^2}{2n(n-1)},
\tag{B4}
```

which is order `n` when `q_n=Theta(n^{3/2})`.  Hence the natural
field-proportional `Theta(sqrt(n))` group needs an order-`n` replenishment
bound.  The unlocalized `12q_n` theorem is larger by order `sqrt(n)`.

## 5. Why a naive use of congestion four cannot close the coefficient

Suppose, optimistically, that all replenishment in a fixed-size group could
be embedded as endpoint allocation and that scale-local harvesting yielded

```math
G\le K(D+E),
```

where `D` is the true norm decrement, `E` is an internal-block error, and
`K=4`.  Since the grouped identity has `D+G=alpha*q_n`, even with `E=0` this
would give only

```math
D\ge\frac{\alpha}{K+1}q_n=\frac{\alpha}{5}q_n.
```

The required decrement is `(1-rho^{3/2})q_n`.  In the continuum limit,

```math
\frac{1-\rho^2}{1-\rho^{3/2}}
=\frac{1+u+u^2+u^3}{1+u+u^2},
\qquad u=\sqrt\rho,
```

which lies strictly between `1` and `4/3`.  Therefore a naive congestion
argument would require

```math
K\le\frac{1-\rho^2}{1-\rho^{3/2}}-1<\frac13.
```

Even `K=1` would not have the needed coefficient; `K=4` is not close.
This no-go concerns the naive scalar charging only, not a future argument
that retains signed internal energy, terminal excess, or unused layer gain.

The global `12q_n` bound is still farther away.  By `(S2)`, its ratio to the
largest possible macroscopic budget is at least

```math
\frac{12}{27/256}=\frac{1024}{9}>113.
```

For `h=o(n)` the ratio is asymptotic to `24n/h`; at the natural
`h=Theta(sqrt(n))` scale it diverges like `Theta(sqrt(n))`.

## 6. The order-nine audit: raw replenishment is too strong

For the exact order-nine minimizer `(10.298)`, every order-eight child has

```math
Q(A[-i])=24,
\qquad q_8=20,
\qquad \varepsilon_i=4,
```

whereas `q_9=24`.  Averaging the singleton identity gives

```math
\mathbb Eg_i=\frac{2q_9}{9}=\frac{16}{3}.
```

The raw scale budget is only

```math
\Delta_{9,1}q_9
=\frac89(16\sqrt2-21)
\approx1.446.
```

Thus the stronger finite claim `E g_H<=Delta*q_n` is false.  The corrected
net value is

```math
\mathbb E(g_i-\varepsilon_i)=\frac43,
```

which satisfies the scale comparison.  Any viable theorem must retain the
terminal excess or allow an asymptotic/grouped error; it cannot simply bound
raw one-step replenishment at the sharp coefficient.

## 7. Exact missing lemma and the required uniformity

For a peeling segment starting from an exact minimizer and stopping at a
deterministic order `m=rho*n`, `(10.340)` says

```math
q_n-q_m
=2\mathbb E\sum_{t<L}a_t
-\mathbb E\sum_{t<L}c_t
-\mathbb E\sum_{t<L}(g_t-\varepsilon_L).
```

The exact scale-preserving condition is therefore

```math
\boxed{
\mathbb E\left[\sum_{t<L}g_t-\varepsilon_L\right]
\le
2\mathbb E\sum_{t<L}a_t
-\mathbb E\sum_{t<L}c_t
-\bigl(1-\rho^{3/2}\bigr)q_n
+E_{n,m}.
}
\tag{C1}
```

It gives

```math
\frac{q_m}{m^{3/2}}
\le\frac{q_n}{n^{3/2}}+\frac{E_{n,m}}{m^{3/2}}.
\tag{C2}
```

Condition `(C1)` is the weakest scalar condition but, by `(10.340)`, is
equivalent to `(C2)`.  The genuinely missing structural lemma must therefore
do more than repackage `(C1)`: it must couple the temporal obligations
`g_t` to the endpoint allocations in `(W1)`--`(W4)` **in the same order/size
window**, retain the terminal excess, and recover the centered
`(1-rho^{3/2})q_n` baseline.  A bound using only total chain mass four and
the uncentered inequality `L_t<=d_t+Q(D_t)` cannot do this, by Section 5.

Take `E_{n,m}>=0` (replace any signed error by its positive part).  A bare
error `E_{n,m}=o(n^{3/2})`, uniform for
`rho in [1/2,1]`, is not enough by itself to bridge a multiplicatively sparse
liminf subsequence: repeated dyadic comparisons may accumulate arbitrarily
many `o(1)` normalized errors.  One sufficient exact uniformity condition is

```math
\boxed{
\Omega(N)
=\sup_{M\ge N}
\sum_{j=0}^{J(M,N)-1}
\frac{E_{n_j,n_{j+1}}}{n_{j+1}^{3/2}}
\longrightarrow0,
}
\tag{C3}
```

where `n_0=M`, `n_J=N`, and every allowed comparison has
`n_{j+1}/n_j in [1/2,1)`.  Under `(C3)`, compare an arbitrary order `N`
directly down from the next member `M` of a liminf subsequence, using
factor-two steps and one final allowed step that lands exactly at `N`.
Then `L_N<=L_M+Omega(N)`, so the normalized sequence converges.

For singleton steps, the familiar sufficient specialization is

```math
\sum_n\frac{E_{n,1}}{n^{3/2}}<\infty.
```

For the **error accumulation along** a fixed-ratio chain `rho_0`, it is enough
to have a normalized envelope `eta(n)` with

```math
\sum_{j\ge0}\eta(\lfloor rho_0^{-j}N\rfloor)\longrightarrow0.
```

Polynomial decay, or `eta(n)=O((log n)^{-1-delta})`, makes this chain-tail
vanish; merely `eta(n)->0` does not.  A theorem available only at one fixed
ratio still permits log-periodic oscillation and does not by itself prove the
full limit.  One also needs a ratio-dense family of allowed comparisons (for
example every ratio in `[1/2,1)`, or ratios tending to one).
