# Would genuinely linear weighted BH rescue nonlinear filters?

2026-09-16. Status: **proved; independently reconstructed by root and the paper verifier**.
This is a scoped coefficient-certificate result, not a theorem about the
moment law of exact minimizers or convergence of the original minimum.

## 1. Normalization

For `s>=0`, define

```math
\mathcal W_s(f)^2=\sum_{r\ge1}r^{-2s}
 \left(\sum_{|S|=r}|\widehat f(S)|^{2r/(r+1)}\right)^{(r+1)/r}.
```

The paper uses `s=5`; the independently refined square weights `r^-5`
use `s=5/2`. Let

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\quad
Q(A)=\|H_A\|_\infty,\quad
Y_A=H_A/n^{3/2},\quad
N=\binom n2,\quad w_s=2^{-s-3/4}.
```

Every signing has

```math
\mathcal W_s(Y_A)
=2^{-s}N^{3/4}/n^{3/2}
=w_s(1-1/n)^{3/4}.                              \tag{1.1}
```

Consider the hypothetical **all-degree, all-dimension** inequality

```math
\mathcal W_s(f)\le C m\|f\|_\infty,
\qquad \deg f\le m.                            \tag{1.2}
```

The conclusions below also hold if (1.2) is restricted to globally even
Walsh polynomials. An asymptotic-only coefficient of `m`, with a finite
additive intercept, is a different statement; see Section 5.

## 2. Uniform nonlinear remainder: no spectral hypothesis

Let `P(t)=sum_(j=0)^d a_j t^j` satisfy `||P||_[-c,c]<=1` for a fixed
`c>0`. Bernstein--Walsh and Cauchy give

```math
|a_j|\le\left(\frac{ed}{cj}\right)^j\quad(j\ge1),
\qquad |a_1|=|P'(0)|\le d/c.                    \tag{2.1}
```

These estimates hold for complex coefficients as well. For the second
one, the usual Bernstein derivative bound at the center of the interval
can also be applied to the appropriately phase-rotated real part.

Degree-two hypercontractivity, uniformly over **all** full signings,
gives

```math
\|Y_A\|_{2j}\le(2j-1)\sqrt N/n^{3/2}
\le\sqrt2 j/\sqrt n.                            \tag{2.2}
```

No bound on the operator norm or on the cap is used here. Layerwise
norm comparison and Parseval give, for every Walsh polynomial `g`,

```math
\mathcal W_s(g)^2
\le en\sum_{r\ge1}r^{-2s-1}\|\widehat g^{=r}\|_2^2
\le en\|g\|_2^2.                               \tag{2.3}
```

Since `W_s` is a seminorm, equations (2.1)--(2.3) imply, with
`z=sqrt(2)e d/(c sqrt(n))<1`,

```math
\mathcal W_s\left(\sum_{j=2}^d a_jY_A^j\right)
\le\sqrt{en}\sum_{j=2}^d z^j
\le\sqrt{en}\frac{z^2}{1-z}.                    \tag{2.4}
```

Therefore, whenever `d=o(sqrt(n))`,

```math
\mathcal W_s(P(Y_A))
=w_s(1-1/n)^{3/4}|P'(0)|+O_c(d^2/\sqrt n),       \tag{2.5}
```

uniformly in the signing and in every admissible real or complex `P`.
The assertion is a two-sided estimate: apply the seminorm triangle
inequality in both directions to the linear term plus the remainder.

For a fixed `P`, let `n` tend to infinity first. Equation (2.5) says that
the entire weighted functional tends to `w_s|P'(0)|`, independently of
the signing. The stronger diagonal assertion allows `P` to vary
arbitrarily with `n`, provided its degree is `o(sqrt(n))`.

## 3. Exact normalized efficiency and cap threshold

If the univariate degree of `P` is exactly `delta` and `2delta<=n`, its
reduced Walsh degree is exactly `2delta`. Indeed the top coefficient on
each `2delta`-set is the nonzero leading coefficient of `P` times
`delta! haf(A[S])/n^(3delta/2)`, and every such hafnian is an odd integer.
Thus there is no hidden opportunity to divide by a smaller actual degree
in the range under consideration.

Define

```math
\mathcal E_{n,d}(A;c)=
\sup_{\substack{1\le\deg P\le d\\\|P\|_{[-c,c]}\le1}}
 \frac{\mathcal W_s(P(Y_A))}{2\deg P}.
```

Using (2.5) with the actual univariate degree `delta`, Bernstein's
bound, and `P(t)=t/c` for the lower bound, we get

```math
\mathcal E_{n,d}(A;c)
=\frac{w_s}{2c}+O_c\left(\frac d{\sqrt n}+\frac1n\right),
\qquad d=o(\sqrt n),                            \tag{3.1}
```

uniformly over **all** full signings. For a fixed odd degree, the same
leading efficiency is achieved by `P(t)=T_d(t/c)` since
`|T_d'(0)|=d`. For even degree, the odd part has degree at most `d-1`,
so its derivative cannot improve the efficiency.

Under a proposed cap `Q(A)<=c n^(3/2)`, inequality (1.2) supplies the
necessary tests `E_(n,d)(A;c)<=C`. Equation (3.1) says their asymptotic
threshold is exactly

```math
c_{\rm test}(C,s)=\frac{w_s}{2C}.                \tag{3.2}
```

This is **the same threshold obtained by applying (1.2) directly to
the quadratic form**. The Chebyshev derivative grows by `d`, but the
degree penalty grows by the same factor. This conclusion is stronger
than merely saying that a particular low-cap construction passes the
tests: the normalized functional has the same limit for every signing.

## 4. Concrete restrictions on an all-degree constant C

The order-six full signing stored in
`computations/results/exact_m6.json` has `N=15` and cap `5`; the cap and
all spin values were replayed by the companion power-check script.
Applying (1.2) at degree two forces

```math
C\ge\frac{2^{-s}15^{3/4}}{10}.
```

Consequently (3.2) can never exceed

```math
c_{\rm test}(C,s)
\le\frac5{30^{3/4}}
=0.3900578865534527\ldots.                      \tag{4.1}
```

This is below the currently verified lower constant
`0.4333221116640807`. It applies even if the putative linear inequality
is restricted to globally even functions: the order-six witness is
already quadratic and globally even. Thus an all-degree linear weighted
BH inequality of this form cannot improve the current cap lower bound
through these nonlinear filters.

For comparison, the quadratic CHSH polynomial
`x_1x_2+x_1x_3+x_4x_2-x_4x_3` has cap `2` and coefficient
`ell_(4/3)` norm `2sqrt(2)`. It forces `C>=2^(-s-1/2)` and the weaker
ceiling `2^(-5/4)=0.4204482076268573...`. The order-six example is
strictly more restrictive.

## 5. Do not confuse an all-degree constant with an asymptotic slope

Suppose instead that one could prove

```math
\mathcal W_s(f)\le b_m\|f\|_\infty,
\qquad b_m/m\longrightarrow C_\infty,
```

or a bound with a finite additive intercept. A fixed degree-two witness
does **not** by itself force `C_infinity` to satisfy the order-six floor.
One must not use (4.1) to exclude this distinct refinement.

For unrestricted polynomials, there is nevertheless the elementary
lower slope `C_infinity>=1`. For each fixed odd integer `m`, let

```math
f_n(x)=T_m\left(\frac{x_1+\cdots+x_n}{n}\right).
```

Its cube norm is at most one and its Walsh degree at most `m`.
As `n->infinity` **first for this fixed m**, its linear Fourier
`ell_1` norm tends to `|T_m'(0)|=m`: the linear Taylor term contributes
`m/n` per coordinate and the higher odd terms contribute `O_m(n^-2)`
per coordinate. Hence `W_s(f_n)>=m-o_n(1)` because the level-one weight
is one. Choosing `n` sufficiently large for each odd `m->infinity`
proves `C_infinity>=1`.

Accordingly, an unrestricted asymptotic-slope theorem still has filter
threshold at most `w_s/2<=2^(-7/4)=0.2973017787506803...` for `s>=0`.
The relevant possible escape would have to exploit a genuinely
globally-even or otherwise tailored class, not merely optimize the
unrestricted weighted growth exponent.

For the globally even class, a separate full-sign lower slope is valid.
Take the regular tensor controls

```math
R_4=J_4-2I_4,\quad B_k=R_4^{\otimes k},\quad n=4^k,\quad
A_k=B_k-(-1)^kI.
```

They satisfy `||A_k||op=sqrt(n)+1` and exactly
`Q(A_k)=n(sqrt(n)+1)/2`. The upper bound is spectral; it is attained
by the all-ones tensor or by a tensor with one balanced `R_4` eigenvector,
depending on the diagonal sign. Thus `Q(A_k)/n^(3/2)->1/2`.

For each fixed odd `d`, normalize by the **actual** cap and take
`f_k=T_d(H_(A_k)/Q(A_k))`. Its norm is at most one, degree `2d`, and
(2.5), with the normalization tending to `c=1/2`, gives

```math
\mathcal W_s(f_k)\longrightarrow 2w_s d.
```

Taking `k->infinity` first and then odd `d->infinity`, or choosing a
sufficiently slow diagonal `d=o(sqrt(n))`, proves the even-class lower
slope `C_infinity>=w_s`. Thus its filter threshold is at most `1/2`.

More intrinsically, let `c_infinity=liminf M_n/n^(3/2)`. Along a
liminf-realizing family, normalize by the actual cap and use the same
uniform all-sign estimate (no spectral hypothesis is needed). It yields

```math
C_\infty^{\rm even}\ge\frac{w_s}{2c_\infty}.
```

This is a necessary constant relation, not a new lower bound or
convergence mechanism. An even-class asymptotic-slope theorem with a
sharp constant could in principle be useful; its useful constant would
have to encode a sharper uniform cap bound. The present polynomial
growth theorem supplies neither that constant nor such a slope theorem.

## 6. Research conclusion

For a genuinely linear **all-degree** weighted inequality, fixed-degree
and `o(sqrt(n))`-degree nonlinear filters have exactly the same leading
efficiency as the original quadratic test. The explicit order-six cap
forces their best possible threshold below the current lower frontier.
For an **asymptotic-only** even-class linear slope, that finite-degree
ceiling is not justified; the remaining obstruction is the full-sign
constant relation above. These quantifiers are materially different.
