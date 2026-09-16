# Boolean Bohr radii and flat quadratic extrema: a degree-localization separation

2026-09-16. Status: **primary-source reconstruction and derived corollaries;
independent proof-verifier audit PASS**. This bounded connection investigates one
specific variational parameter, not ordinary degree-two BH or power tests.

## 1. Primary sources and exact domains

The relevant primary source is Defant--Mastylo--Perez,
[Bohr's phenomenon for functions on the Boolean cube](https://arxiv.org/html/1707.09186),
Definitions 1.1--1.2, Theorem 2.1, Theorem 3.1, Lemma 3.3, and the proof
of Corollary 3.2. Its scalar field is **real** and its supremum is over
the real Boolean cube. The proof of the homogeneous upper bound expressly
uses coefficients `+1` or `-1` on **every** set of the chosen degree.

For contrast, Bayart--Pellegrino--Seoane-Sepulveda,
[The Bohr radius of the n-dimensional polydisk](https://arxiv.org/html/1310.2834),
defines the complex-polydisc radius using arbitrary analytic monomials
and the supremum on the complex polydisc. That is not the same norm.
The real Boolean Sidon invariant in Defant et al.,
[Asymptotic insights for projection, Gordon--Lewis and Sidon constants](https://arxiv.org/html/2302.00233),
Section 2 and Theorem 5.1, ranges over arbitrary real coefficient
magnitudes; its comparisons do not preserve the nonlinear flat slice.

For real cube functions, the radius of a class is the largest `r` such
that `sum_S |fhat(S)|r^|S|<=||f||infinity` for every member. The all-real
class has exact radius `2^(1/n)-1`, whereas the union of homogeneous
classes has radius asymptotic to `sqrt(log(n)/n)`. Thus one must not
import the complex all-polynomial radius into the real all-function class.
Nor is changing only the scalar field harmless: for complex-valued
functions on the same real cube, the all-function radius is zero.
Indeed, `f(x)=1+i epsilon x_1` has norm `sqrt(1+epsilon^2)`, whereas
the weighted absolute coefficient sum is `1+r epsilon`; letting
`epsilon` decrease to zero excludes every positive `r`.

## 2. Flat homogeneous radius: the exact variational parameter

Define the complete-support degree-`m` Littlewood minimum and its radius by

```math
M_{n,m}=\min_{a_S\in\{\pm1\}}
 \left\|\sum_{|S|=m}a_Sx^S\right\|_\infty,
\qquad
r^{\rm flat}_{n,m}=
 \left(\frac{M_{n,m}}{\binom nm}\right)^{1/m}.
```

Let `r_(n,m)` be the radius for all real `m`-homogeneous Walsh
polynomials, and set

```math
R_n=\min_{1\le m\le n}r_{n,m},\qquad
R_n^{\rm flat}=\min_{1\le m\le n}r^{\rm flat}_{n,m}.
```

The constants and the zero polynomial do not affect these minima.
The definitions give `r_(n,m)<=r_flat_(n,m)`.

DMP's homogeneous theorem, including its **flat-sign upper proof**, gives

```math
c_m\left(\frac n{\binom nm}\right)^{1/(2m)}
\le r_{n,m}\le r^{\rm flat}_{n,m}
\le C_0^{1/m}\left(\frac n{\binom nm}\right)^{1/(2m)},          \tag{2.1}
```

uniformly for `1<=m<=n`, with `c_m->1` and absolute `C_0`.
Specifically, its Lemma 3.3 permits `C_0=6sqrt(log 2)` by choosing
`c_S=1` on the whole degree-`m` level. The source's original lower factor
is `c_m=m^(-1/(2m)) exp[-O(sqrt(log(m)/m))]`.

It follows by squeezing that

```math
R_n^{\rm flat}\sim R_n\sim\sqrt{\frac{\log n}{n}}.           \tag{2.2}
```

This is a legitimate positive implication **for complete-support real
sign polynomials**, but it minimizes over their degrees. It is not a
statement about fixed quadratic degree. Also, the radius-level relation
`r_flat/r->1` for `m->infinity` is an `m`th-root relation; it does not
assert that the corresponding flat and unrestricted Sidon constants
have ratio tending to one.

## 3. The discriminating theorem: all near-minimizing degrees are logarithmic

Put `L=log n`, `t=m/L`, and

```math
\Phi(t)=\frac12\left(\log t-1+\frac1t\right).
```

From (2.1) and `binom(n,m)<=(en/m)^m`,

```math
\frac{r_{n,m}}{\sqrt{L/n}}
\ge c_m\exp(\Phi(m/L)).                                      \tag{3.1}
```

The function `Phi` is nonnegative, tends to infinity at both endpoints
of `(0,infinity)`, and vanishes only at `t=1`. Taking
`m=floor(log n)` in the flat upper bound and using Stirling gives the
matching upper value `sqrt(L/n)(1+o(1))`.

**Degree localization.** If either `r_(n,m_n)` or `r_flat_(n,m_n)`
is asymptotic to `sqrt(log(n)/n)`, then

```math
\frac{m_n}{\log n}\longrightarrow1.                          \tag{3.2}
```

Proof: bounded degrees have radius larger by a diverging factor, so
`m_n->infinity` and `c_(m_n)->1`. Equation (3.1) then forces
`Phi(m_n/log n)->0`, which is equivalent to (3.2).

In particular, for **each fixed D**, removing every degree `1,...,D`
from either variational minimum leaves the minimum **exactly unchanged
for all sufficiently large n**. Every such fixed-degree radius divided
by `sqrt(log(n)/n)` diverges, while degree `floor(log n)` provides a
bounded ratio. This is a precise separation: the aggregate Bohr-radius
asymptotic eventually does not even consult the quadratic entry.

## 4. What polynomial BH really improves here

Suppose the Boolean BH constants obey `B_m<=K m^a`, with `K>=1`
and fixed `a>=0`. The new paper gives `a=27`; the independently
refined exponent changes constants but not the rate below. Repeating
the source's homogeneous lower-bound proof replaces `c_m` in (2.1) by

```math
\widetilde c_m=
\exp\left[-\frac{(a+1/2)\log m+\log K}{m}\right].            \tag{4.1}
```

Consequently both radius minima obey the quantitative estimate

```math
\frac{R_n}{\sqrt{\log(n)/n}},\quad
\frac{R_n^{\rm flat}}{\sqrt{\log(n)/n}}
=1+O\left(\frac{\log\log n}{\log n}\right).                 \tag{4.2}
```

Here is a uniform proof that includes all integer degrees. Let
`A=a+1/2`, `B=log K`, and write

```math
\log\frac{r_{n,m}}{\sqrt{L/n}}
\ge\Phi(m/L)-\frac{A\log m+B}{m}.                            \tag{4.3}
```

If `m<=L^2`, put `L'=L-4A log L-2B>0`. Since `log m<=2log L`,
the right side is at least

```math
\Phi(m/L')+\frac12\log(L'/L)
\ge\frac12\log(L'/L)
=-O(\log L/L).
```

If `m>L^2`, then `Phi(m/L)>=(log L-1)/2`, while the subtracted
term is at most `(2A log L+B)/L^2` for large `L`; hence it is positive.
This proves the uniform lower half of (4.2). For the upper half use
`m=floor L` and

```math
\log\left[
\frac{C_0^{1/m}(n/\binom nm)^{1/(2m)}}{\sqrt{L/n}}
\right]
=\Phi(m/L)+\frac{\log(2\pi m)}{4m}
 +\frac{\log C_0}{m}+O(m^{-2}+m/n)
=O(\log L/L).
```

Combining (4.3) with the upper bound shows that every **exact** minimizing
degree, for either radius, belongs to the window

```math
m=\log n+O\big(\sqrt{\log n\,\log\log n}\big).              \tag{4.4}
```

First use (3.2) to put `m/L` in a fixed neighborhood of one; then the
penalty in (4.3) is `O(log L/L)` and `Phi(t)` is comparable to
`(t-1)^2` there. General `1+o(1)`-near minimizers need only satisfy
(3.2), not necessarily the sharper window (4.4).

Equation (4.2) is a derived quantitative consequence, not a claim of
literature novelty. It improves the rate supplied by the older
subexponential-constant proof, not the leading Bohr-radius constant and
not the original quadratic minimax bound.

## 5. Why this does not settle the quadratic question

For the original problem, `M_n=M_(n,2)` exactly, and

```math
r^{\rm flat}_{n,2}=\sqrt{M_n/\binom n2},\qquad
n^{1/4}r^{\rm flat}_{n,2}
=\sqrt{\frac{2(M_n/n^{3/2})}{1-1/n}}.                         \tag{5.1}
```

Thus convergence of this **fixed-degree flat** radius is equivalent to
the original convergence question. The established scale is
`r_flat_(n,2)=Theta(n^-1/4)`, much larger than the degree-minimized
`sqrt(log(n)/n)`. By Section 3, even exact knowledge of the latter
minimum for all large dimensions does not recover the quadratic entry
through this variational identity: that entry is eventually never selected.
This is not an impossibility theorem for a future additional relation
between the degree-two and logarithmic-degree flat classes.

For arbitrary coefficients the homogeneous identity is
`r_(n,2)=Sidon(B_=2^n)^(-1/2)`; the flat restriction need not attain that
Sidon constant. The 2024 source's comparison with a projection constant
has a fixed degree-dependent factor. It is neither an equality on the
flat slice nor an asymptotically lossless replacement in fixed degree.

This coefficient-class gap is already strict in dimension four.
Parseval gives `M_(4,2)>=sqrt(6)`; every value of a six-term signing is
an even integer, so `M_(4,2)>=4`. The full signing in Section 6 attains
four. Thus the flat Sidon ratio is exactly `6/4=3/2`. In contrast,
the sparse CHSH polynomial

```math
x_1x_3+x_1x_4+x_2x_3-x_2x_4
```

has coefficient absolute sum four and cube norm two, proving that the
unrestricted real quadratic Sidon constant is at least `2>3/2`.
Equivalently, `r_(4,2)<=1/sqrt(2)<sqrt(2/3)=r_flat_(4,2)`.

There is also no remedy by degree complementation. Multiplication by
`x_1...x_n` sends a full degree-two signing to a full degree-`n-2`
signing with the same cap, so

```math
r^{\rm flat}_{n,n-2}
=\big(r^{\rm flat}_{n,2}\big)^{2/(n-2)}\longrightarrow1.
```

That high-degree encoding remains far from the logarithmic-degree
variational minimum.

## 6. An exact full-sign real/complex norm separation

The complete four-variable signing

```math
P(z)=z_1z_2-z_3z_4+(z_1+z_2)(z_3+z_4)
```

has six coefficients of absolute value one. Its real Boolean norm is
`4`, but its complex-polytorus norm is exactly `2sqrt(5)>4`.
For the latter, write `z_1=e^(i(alpha+u))`, `z_2=e^(i(alpha-u))`
and similarly `z_3,z_4` with `beta,v`. After removing a unit phase,

```math
P=2i\sin(\alpha-\beta)+4\cos u\cos v,
\qquad |P|^2\le4+16=20,
```

with equality when `u=v=0` and `alpha-beta=pi/2`. On real signs,
checking whether each pair agrees gives only magnitudes `0,2,4`.
An exact integer check over all sixteen cube points gives the value set
`{-4,-2,0,2,4}`; the torus point `(i,i,1,1)` gives `-2+4i`.
Thus even complete-support real sign coefficients do not identify the
complex and Boolean supremum norms.

## 7. Conclusion

The Boolean Bohr parameter does have a genuine flat-sign implication and
a better convergence rate under polynomial BH. Its minimizing degree is
logarithmic, however, and fixed quadratic degree is eventually absent
from the defining minimum. The precise missing input remains a theorem
for the fixed-degree flat radius (5.1), not the known aggregate Bohr limit.
No original-value bound or convergence mechanism follows from this
connection.
