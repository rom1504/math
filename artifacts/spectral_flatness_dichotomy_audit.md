# Spectral-flatness dichotomy audit

Status: exact inequalities proved; the proposed dichotomy is **not reduced to
a useful cap theorem**. The natural Gaussian route stops at normalized edge
cap `1/pi`, and the trace-four route is weaker than the existing universal
lower bound.

Write

```math
\kappa(A)=\max_{x\in\{\pm1\}^n}|x^{\mathsf T}Ax|=2\operatorname{cap}(A)
```

for a symmetric zero-diagonal sign matrix `A`.

## 1. Exact flatness identity

Every such matrix has `(A^2)_{ii}=n-1`. Therefore, with

```math
D=A^2-(n-1)I,
```

one has the exact identity

```math
\|D\|_F^2
=\operatorname{tr}(A^4)-n(n-1)^2
=\sum_{i\ne j}(A^2_{ij})^2.                          \tag{1}
```

The right side is nonnegative, so the minimum possible fourth spectral
moment is `n(n-1)^2`; equality is exactly the conference identity. Thus
`tr(A^4)=n(n-1)^2+o(n^3)` is precisely Frobenius-near-conference, but says
only `||D||_F=o(n^(3/2))`. It does not by itself give a small edit distance
to an exact conference matrix.

## 2. Exact Gaussian polynomial rounding and its barrier

For every real `t`, the matrix

```math
R_t=\frac{(A+tI)^2}{n-1+t^2}
```

is a correlation matrix: it is positive semidefinite and has diagonal one.
Gaussian hyperplane rounding gives a Boolean `x` satisfying

```math
\kappa(A)\ge
\left|\frac2\pi\sum_{i\ne j}
 \arcsin\!\left(\frac{a_{ij}D_{ij}+2t}{n-1+t^2}\right)\right|.      \tag{2}
```

This follows from
`E[sign(g_i)sign(g_j)]=(2/pi)arcsin((R_t)_{ij})` and oddness of
`arcsin`; no independence approximation is used.

For an exact conference matrix, `D=0`. Optimizing the common correlation in
(2) gives `t=sqrt(n-1)` and only

```math
\operatorname{cap}(A)
\ge \frac1\pi n(n-1)\arcsin\!\frac1{\sqrt{n-1}}
=\left(\frac1\pi+o(1)\right)n^{3/2}.                \tag{3}
```

The same barrier is robust in the flat branch. At `t=sqrt(n-1)`, every
argument in (2) stays uniformly away from `+1` and `-1` for large `n`, since
`|D_ij|<=n-2`. The Lipschitz error from replacing `D` by zero is at most a
universal constant times

```math
\frac{\sum_{i\ne j}|D_{ij}|}{n}
\le \|D\|_F.                                         \tag{4}
```

Consequently (1)--(2) prove only

```math
\operatorname{tr}(A^4)=n(n-1)^2+o(n^3)
\quad\Longrightarrow\quad
\operatorname{cap}(A)\ge(1/\pi-o(1))n^{3/2}.        \tag{5}
```

This is below the already known universal constant
`0.3364933644...`; it does not make the near-conference branch easier.

## 3. A trace-four cap inequality and the localization loss

There is also a direct exact inequality, but its constants expose the other
barrier. Let

```math
L^2=\max_i(A^4)_{ii}.
```

In the bipartite vector relaxation, take the left vector at `i` to be row
`i` of `A^2`, divided by `L` and padded orthogonally to unit length, and the
right vector at `j` to be row `j` of `A`, divided by `sqrt(n-1)`. Their
objective is exactly

```math
\frac{\operatorname{tr}(A^4)}{L\sqrt{n-1}}.          \tag{6}
```

Real Grothendieck rounding, followed by the elementary same-sign reduction,
therefore gives

```math
\operatorname{cap}(A)
\ge \frac{\operatorname{tr}(A^4)}
 {4K_G\sqrt{(n-1)\max_i(A^4)_{ii}}}.                 \tag{7}
```

For completeness, the same-sign reduction costs the factor two: if
`beta(A)=max_{x,y}|x^T A y|`, write `x=u+v`, `y=u-v` with disjoint ternary
`u,v`. Then `x^TAy=u^TAu-v^TAv`; each ternary quadratic value is at most
`kappa(A)` by randomly filling its zero coordinates. Hence
`beta(A)<=2kappa(A)`.

Inequality (7) depends not only on spectral spread but on how the fourth
moment is distributed among coordinates. Even under the favorable condition

```math
\max_i(A^4)_{ii}\le C\,\operatorname{tr}(A^4)/n,
```

it yields only

```math
\operatorname{cap}(A)
\ge \frac1{4K_G\sqrt C}
\sqrt{\frac n{n-1}\operatorname{tr}(A^4)},           \tag{8}
```

whose coefficient near the flat threshold is about `1/(4K_G)`, far below
the current universal bound. Without the coordinate-spread hypothesis, the
maximum diagonal term in (7) can absorb the entire fourth-moment gain.

## 4. Research judgment

The exact spectral identity (1) gives a clean algebraic dichotomy, but not
the needed mathematical one:

- the near-flat branch is only Frobenius-near an algebraic conference
  identity; no rounding/edit theorem maps it to an exact conference family,
  and the canonical Gaussian certificate stops at `1/pi`;
- the spread branch needs a new delocalization-sensitive Boolean rounding
  theorem. Trace `(A^4)` alone, through the available exact relaxation,
  neither reaches normalized cap `1/2` nor improves the known lower bound.

Accordingly spectral flatness is not presently a reduction of the landing or
convergence obligation. Progress would require a new theorem converting
fourth-moment excess into same-switch Boolean energy with essentially sharp
constant, or a rigidity theorem turning `||A^2-(n-1)I||_F=o(n^(3/2))` into a
cap-controlled exact structured approximation. The calculations above give
neither, so this route should stop rather than be counted as a new leading
mechanism.
