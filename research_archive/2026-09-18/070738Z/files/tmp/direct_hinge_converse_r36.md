# Wave 36 main attack: a row-sensitive converse to rare cut coverage

## Status

The fixed-slice upper-tail inequality and its consequence for the normalized
restriction gap are **Verified**, conditional only on the already proved
quadratic decoupling/Hanson--Wright step used in (10.690).  The proof below
replaces the Bernstein estimate for the linear Hoeffding component by a
global Hoeffding bound, so no `R_infty` term is paid.  The result does not
prove convergence.  It shows that the arbitrary-cut/entropic-hinge target
cannot be realized as a generic rare fluctuation across a genuine constant
normalized gap: at its advertised row and entropy scales it already forces
the direct `q_m` recurrence.

## 1. A fixed-slice row-sensitive upper tail

Fix an exact order-`n` minimizer `A`, an oriented cut
`d=(sigma,x)`, and `m`, with

```math
p=m/n,
\qquad p_2=(m)_2/(n)_2,
\qquad
B=\sigma D_xAD_x.
```

Thus `B` is symmetric with zero diagonal and entries of modulus one,

```math
\|B\|_F^2=n(n-1),
\qquad
\|B\|_{\rm op}=\|A\|_{\rm op},
\qquad
\|B\mathbf1\|_2^2=R_2(d),
\qquad
|\mathbf1^TB\mathbf1|\le q_n.
\tag{C36.1}
```

For a selector indicator `xi`,

```math
X_d(S)=\xi^TB\xi-p_2\mathbf1^TB\mathbf1.
```

First take independent `eta_i~Ber(p)` and put `zeta=eta-p1`.  The exact
decomposition is

```math
X_d(\eta)
=(p^2-p_2)\mathbf1^TB\mathbf1
+2p(B\mathbf1)^T\zeta+\zeta^TB\zeta,
\qquad
p^2-p_2=\frac{p(1-p)}{n-1}.
\tag{C36.2}
```

Let

```math
b_0=(p^2-p_2)q_n,
\qquad
v=(u-b_0)_+.
```

Global Hoeffding for the independent linear summands gives

```math
\Pr\{2p(B\mathbf1)^T\zeta\ge v/2\}
\le\exp\{-c v^2/R_2(d)\},
\tag{C36.3}
```

with the usual zero-row interpretation.  The same decoupling and Gaussian
linearization used to prove (10.690), equivalently the standard bounded
Hanson--Wright tail, gives

```math
\Pr\{\zeta^TB\zeta\ge v/2\}
\le2\exp\left\{-c\min\left(
\frac{v^2}{n^2},\frac v{\|A\|_{\rm op}}
\right)\right\}.
\tag{C36.4}
```

Conditioning on `sum eta_i=m` gives the uniform slice.  Since `m` is a mode
of `Bin(n,m/n)`, the conditioning probability is at least `1/(n+1)`.
Combining (C36.2)--(C36.4) proves, with a universal `c>0`,

```math
\boxed{
U_m\{X_d\ge u\}
\le3(n+1)\exp\left\{-c\min\left(
\frac{v^2}{R_2(d)},
\frac{v^2}{n^2},
\frac v{\|A\|_{\rm op}}
\right)\right\}.
}
\tag{C36.5}
```

This differs importantly from the Bernstein presentation (10.691): the
linear term is bounded by its full squared row norm at all deviations and
does not introduce `R_infty`.

## 2. Coverage forces the direct restriction recurrence

Let

```math
h_d(S)=\widehat\ell(S,d)=Y_A(S)-X_d(S),
\qquad
G_{n,m}=q_m-p^{3/2}q_n.
```

Since every principal signing has `Q(A[S])>=q_m`, one has pointwise

```math
Y_A(S)\ge G_{n,m}.
\tag{C36.6}
```

Suppose a cut satisfies

```math
R_2(d)\le R,
\qquad
U_m\{h_d\le t\}\ge e^{-L}.
\tag{C36.7}
```

The event in (C36.7) is contained in
`{X_d>=G_(n,m)-t}`.  Comparing its lower probability with (C36.5), and
putting `H=L+O(log n)`, gives the general deterministic consequence

```math
\boxed{
G_{n,m}
\le t+b_0+C\max\left\{
\sqrt{RH},\ n\sqrt H,\ \|A\|_{\rm op}H
\right\}.
}
\tag{C36.8}
```

All constants are uniform for fixed ratio windows.  The exact-minimizer
bounds `q_n=O(n^(3/2))` and `||A||_op=O(n^(3/4))` make
`b_0=O(n^(1/2))`.

At the arbitrary-cut scales

```math
R=O(n^{9/4-c}),
\qquad L=O(n^{3/4-c}),
\qquad t=O(n^{3/2-c}),
\qquad 0<c<1/4,
```

the three inverse-tail terms are respectively

```math
O(n^{3/2-c}),
\qquad
O(n^{11/8-c/2})=O(n^{3/2-c}),
\qquad
O(n^{3/2-c}).
```

Therefore

```math
\boxed{
q_m\le(m/n)^{3/2}q_n+O(n^{3/2-c}).
}
\tag{C36.9}
```

This recurrence, with the same ratio-window uniformity and exact landing as
the restriction scheme, already proves convergence directly.  Thus (10.795)
is still a valid sufficient lemma, but its advertised tail cannot be supplied
by a generic large deviation while a fixed normalized `q_m/m^(3/2)` versus
`q_n/n^(3/2)` gap persists.

More quantitatively, if `G_(n,m)>=delta n^(3/2)`, `t=o(n^(3/2))`, and
`L=O(n^(3/4-c))`, then (C36.8), the Frobenius scale, and the exact-minimizer
operator bound force

```math
R_2(d)=\Omega(n^{9/4+c}).
\tag{C36.10}
```

Hence a cut covering `e^{-O(n^(3/4-c))}` selectors across a constant gap
must violate the (10.795) row budget by a factor `n^(2c)`.

## 3. Entropic-hinge version

The same proof has a transport form.  Cauchy--Schwarz, global Hoeffding for
the linear part, and the quadratic exponential bound behind (C36.4) give,
for `0<theta<c/||A||_op`,

```math
\log\mathbb E_{U_m}e^{\theta(X_d-b_0)}
\le O(\log n)
+\frac{C\theta^2(R_2(d)+n^2)}{1-C\theta\|A\|_{\rm op}}.
\tag{C36.11}
```

For any selector law `w`, entropy duality followed by optimizing `theta`
therefore yields, with `H_w=D(w||U_m)+O(log n)`,

```math
\mathbb E_wX_d
\le b_0+C\max\left\{
\sqrt{R_2(d)H_w},\ n\sqrt{H_w},\ \|A\|_{\rm op}H_w
\right\}.
\tag{C36.12}
```

Since `h<=t+[h-t]_+`, (C36.6) gives

```math
G_{n,m}
\le t+\mathbb E_w[h_d-t]_++\mathbb E_wX_d.
\tag{C36.13}
```

Consequently the Wave 35 entropic-hinge cost at
`lambda_cut=Theta(n^(-3/4))`, entropy `O(n^(3/4-c))`, row
`O(n^(9/4-c))`, and threshold `O(n^(3/2-c))` already implies (C36.9),
even before its strict soft-to-hard buffer is invoked.  The buffer remains
necessary only for extracting the literal hard-tail statement.

## 4. Scope and strategic use

This is a converse/diagnostic, not a proof of (10.795).  It does not make the
tail target false, and it does not replace minimizer-specific structure.
It says precisely what such structure must accomplish: it must rule out the
normalized restriction gap itself, rather than manufacture the needed mass
from ordinary quadratic fluctuations.  A Wave 36 continuation should either
find a signing-specific mechanism which directly forces (C36.9), or prefer a
route (such as conditional replacement or harmonic cancellation) whose
intermediate lemma contains additional structural information rather than a
generic reverse-tail estimate.
