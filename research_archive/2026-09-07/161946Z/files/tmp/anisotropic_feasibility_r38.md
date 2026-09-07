# Wave 38 follow-up: anisotropic feasibility is ruled out at project entropy

## Verdict

The proposed anisotropic product-noise certificate does **not** evade the
uniform spectral wall.  There are two independent no-go statements.

1. As written, (R38.A5) is empty: every structured complement flip satisfies

   ```math
   |\rho^TC\rho|\le3q_n,
   ```

   so its requirement `rho^T C rho>=4q_n` is impossible.

2. More robustly, coordinatewise Pinsker and the structured spectral bound
   give

   ```math
   \boxed{|\rho^TC\rho|\le6\sqrt{2q_n}\,D(\rho).}
   \tag{R38.F1}
   ```

   Thus even mean energy `kappa q_n` for any fixed `kappa>0` requires
   `D(rho)=Omega(sqrt(q_n))=Omega(n^(3/4))`.  The desired budget
   `D=O(R/q_n)=O(n^(3/4-c))` is smaller by a power.  At that budget the
   energy mean is `o(q_n)`, its variance is `o(q_n^2)`, and the probability
   of reaching energy `q_n` is `o(1)`, not a fixed constant.

This kills the constant-probability independent-product-noise compression
criterion, not merely the particular constants in (R38.A5).  A viable
successor must leave independent product tilting, abandon the project
entropy budget, or exploit a non-product conditioning whose entropy is not
controlled by `sum_i d(rho_i)`.

## 1. Exact box obstruction to the displayed `4q_n` threshold

Write `P=P_S`,

```math
B_S=-A+2PAP,
```

and put `u=D_x rho`.  Gauge and orientation do not change absolute quadratic
values, so

```math
|\rho^TC\rho|=|u^TB_Su|
\le |u^TAu|+2|(Pu)^TA(Pu)|. \tag{R38.F2}
```

For every zero-diagonal symmetric matrix, a quadratic form is affine in each
coordinate separately.  Therefore its maximum and minimum on the box
`[-1,1]^n` occur at cube vertices.  Since `A` is an exact order-`n`
minimizer,

```math
|u^TAu|\le q_n.
```

The same box argument and the principal-cap inequality `Q(A[S])<=q_n` give

```math
|(Pu)^TA(Pu)|\le Q(A[S])\le q_n.
```

Substitution in (R38.F2) proves

```math
\boxed{|\rho^TC\rho|\le3q_n.}
\tag{R38.F3}
```

Consequently no selector is anisotropically regularizable under the literal
definition (R38.A5).  The implication (R38.A8) remains formally true but is
vacuous.

Changing `4` to a smaller constant does not rescue the project-scale method,
as the next section shows.

## 2. Spectral--entropy no-go for every fixed positive mean

The exact-minimizer spectral estimate and the structured-flip identity give

```math
\lVert C\rVert_{op}=\lVert B_S\rVert_{op}
\le3\lVert A\rVert_{op}
\le3\sqrt{2q_n}. \tag{R38.F4}
```

For the Bernoulli relative entropy in (R38.A4),

```math
d(r)\ge\frac{r^2}{2},\qquad -1<r<1. \tag{R38.F5}
```

Indeed both sides are even, and for `r>=0` the derivative difference is
`atanh(r)-r>=0`.  Hence

```math
\lVert\rho\rVert_2^2\le2D(\rho). \tag{R38.F6}
```

Rayleigh, (R38.F4), and (R38.F6) now prove (R38.F1).  In particular,

```math
\rho^TC\rho\ge\kappa q_n
\quad\Longrightarrow\quad
D(\rho)\ge\frac{\kappa}{6\sqrt2}\sqrt{q_n}. \tag{R38.F7}
```

If `D(rho)<=K R/q_n` for a fixed constant `K`, then (R38.F7) forces

```math
\boxed{R\ge\frac{\kappa}{6\sqrt2K}q_n^{3/2}.}
\tag{R38.F8}
```

Since `q_n=Theta(n^(3/2))`, no fixed positive mean fraction of `q_n` is
possible at `R=O(n^(9/4-c))`.

This argument retains the coordinate box: it uses the exact product entropy,
not a Gaussian relaxation, and does not assume `max_i |rho_i|<=1/2`.

## 3. Exact variance and the event probability

The variance identity (R38.A2) also yields a bound without any small-bias
assumption.  Since `0<=v_i<=1`,

```math
\begin{aligned}
4\sum_i v_i(C\rho)_i^2
&\le4\lVert C\rho\rVert_2^2
\le72q_n\lVert\rho\rVert_2^2
\le144q_nD(\rho),\\
2\left\{\left(\sum_iv_i\right)^2-\sum_iv_i^2\right\}
&\le2n(n-1).
\end{aligned}
```

Therefore

```math
\boxed{
\operatorname{Var}(X)
\le144q_nD(\rho)+2n(n-1),
\qquad X=\sigma(x')^TB_Sx'.}
\tag{R38.F9}
```

At the advertised budget

```math
D(\rho)=O(R/q_n)=O(n^{3/4-c}),
```

(R38.F1) gives `|E X|=O(n^(3/2-c))=o(q_n)`, while (R38.F9) gives
`Var(X)=O(n^(9/4-c)+n^2)=o(q_n^2)`.  For all sufficiently large `n`,
Chebyshev thus gives

```math
\Pr\{X\ge q_n\}
\le\frac{\operatorname{Var}(X)}{(q_n-|\mathbb EX|)^2}
=O\left(\frac{D(\rho)}{q_n}+\frac{n^2}{q_n^2}\right)
=o(1). \tag{R38.F10}
```

So one cannot repair (R38.A5) merely by lowering its mean constant and using
a sharper one-sided concentration inequality.  At project entropy the
product law does not put constant mass on the certificate event at all.

## 4. Parent-row lower bound with the baseline retained

There is a second generalized-Rayleigh obstruction when the desired mean
constant `kappa` exceeds `2`.  From (R38.F2),

```math
|u^TB_Su|\ge\kappa q_n
\quad\Longrightarrow\quad
|u^TAu|\ge(\kappa-2)q_n.
```

Put `s=||u||_2^2=||rho||_2^2`.  Cauchy--Schwarz gives

```math
\lVert Au\rVert_2^2
\ge\frac{(\kappa-2)^2q_n^2}{s}. \tag{R38.F11}
```

Combining this with the **exact** baseline term in (R38.A3) gives the
necessary scalar feasibility condition

```math
\boxed{
\mathbb E R_2(x')
\ge\frac{(\kappa-2)^2q_n^2}{s}
 +(n-1)(n-s),
\qquad
0<s\le\min\{n,2D(\rho)\}.}
\tag{R38.F12}
```

If `D<=K R/q_n` and the proposed row condition is
`E R_2<=R/8`, then even after discarding the nonnegative baseline,

```math
\frac R8
\ge\frac{(\kappa-2)^2q_n^3}{2KR},
```

and hence

```math
\boxed{
R\ge\frac{2(\kappa-2)}{\sqrt K}\,q_n^{3/2}.}
\tag{R38.F13}
```

Keeping the baseline only strengthens this to

```math
\frac R8
\ge\frac{(\kappa-2)^2q_n^3}{2KR}
 +(n-1)\left(n-\frac{2KR}{q_n}\right)
\tag{R38.F14}
```

whenever `2KR/q_n<=n`.  The baseline is `Theta(n^2)` at project scale and
is not the leading obstruction; the directional parent-row term recreates
the `q_n^(3/2)` wall.

## 5. Scope

The exact rational checker for (R38.A1)--(R38.A3) was rerun and passes.  The
failure found here is not in those identities.  It is in the feasibility
claim following them:

- the literal `4q_n` energy clause contradicts the exact box cap (R38.F3);
- every fixed positive fraction of `q_n` contradicts the project entropy
  budget by (R38.F7)--(R38.F8); and
- at that budget, even the unconditioned chance of crossing `q_n` vanishes
  by (R38.F10).

This is a proved asymptotic no-go for the independent anisotropic
product-noise implementation.  It does not falsify complement-flip
fractional compression by a different, correlated law.
