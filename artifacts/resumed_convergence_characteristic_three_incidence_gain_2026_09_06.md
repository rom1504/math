# A nonflat characteristic-three lower bound from projective incidence

Date: 2026-09-06. Joint derivation with the director. This is an actual
signing theorem; no spectral-flatness, row-sum, or Gaussian independence
hypothesis is imposed. The original unrestricted convergence problem
remains open.

## 1. Main theorem

For every sequence of hollow real symmetric additive-Cayley signings on
`F_(3^r)`, with `q=3^r -> infinity`,

\[
\boxed{\quad
\liminf\frac{Q(A)}{q^{3/2}}
\ge\frac{2080}{9\sqrt{269441}}
\approx0.4452346798594428 .\quad}                   \tag{1}
\]

Here `Q(A)=max_(x Boolean)|x^T A x|/2`. The theorem covers genuinely
nonflat spectra. It improves the unconditional `4/9` theorem in
`resumed_convergence_cayley_character_mobius_bound_2026_09_06.md`.

An explicit finite bound holds for `q>81`. Put

\[
\rho_q=\frac{q-29}{26(q-3)},\qquad
\gamma_q=\frac{q-81}{81(q-1)}.
\]

Then

\[
\frac{Q(A)}{q^{3/2}}
\ge\frac{4/9}
{\sqrt{1-\rho_q/10+\rho_q^2/(400\gamma_q)}}.
\tag{2}
\]

This statement concerns the Cayley TARGET class in characteristic
three. It neither asserts (1) for arbitrary sign matrices nor reaches
the conjectural one-half threshold for all Cayley signings.

## 2. Exact spectral accounting supplied by the sign kernel

Write `lambda(a)` for the real even Fourier eigenvalues and set

\[
K=\frac{2Q(A)}{q^{3/2}},\quad
t=\frac{\lambda(0)}{\sqrt q},\quad
h(a)=\frac{8\lambda(a)+\lambda(0)}{9\sqrt q}
\quad(a\ne0).
\tag{3}
\]

The function `h` is constant on the two nonzero vectors of a projective
point. Each `h(a)` is the normalized energy of the Boolean three-point
character square wave, so `|h(a)|<=K`. The constant Boolean vector also
gives `|t|<=K`.

Because the kernel has zero diagonal and every other entry is a sign,

    sum_a lambda(a)=0,
    sum_a lambda(a)^2=q(q-1).

Consequently, averaging uniformly over the `L=(q-1)/2` projective
points gives the EXACT identity

\[
\mathbb E h^2=\frac{64}{81}+\gamma_q t^2.
\tag{4}
\]

For `q>=81` this immediately yields `K>=8/9`, hence the earlier
four-ninths lower bound. The remaining argument excludes quantitative
near-equality by using dependencies among two-dimensional subspaces.

## 3. A uniform monochromatic-line density from three-dimensional incidence

Color every projective point by
`s(a)=sign(h(a))`, with either choice at zero. A projective line is the
set of four one-dimensional subspaces in a two-dimensional vector
subspace. Let `p_0,p_1,p_2` be the proportions of projective lines whose
four colors split as `4/0`, `3/1`, and `2/2`.

Inside every three-dimensional vector subspace, its projective plane
has 13 points and 13 lines, each point lying on exactly four lines.
The product of the four signs along a line is negative exactly for a
`3/1` split. Multiplying these products over all 13 lines gives `+1`,
because each point occurs four times. Thus the number of `3/1` lines
is even and at most 12. At least one line has another split. Averaging
over all three-dimensional vector subspaces yields

\[
p_0+p_2\ge\frac1{13}.                            \tag{5}
\]

There is a stronger elementary inequality in each projective plane:
if `N_0,N_2` count its `4/0` and `2/2` lines, then

\[
3N_0+N_2\ge3.                                    \tag{5a}
\]

Only the case `N_0=0` needs proof. Both color classes then meet every
line (they are blocking sets), and neither contains a complete line.
Such a blocking set must have at least six points. A set of at most
three points meets at most 12 lines. For a blocking set of four points
with no full line, let `n_i` count lines meeting it in exactly `i`
points, `i=1,2,3`. Counting lines, incidences, and pairs would give
`n_1+n_2+n_3=13`, `n_1+2n_2+3n_3=16`,
`n_2+3n_3=6`, forcing `n_3=3,n_2=-3`, impossible.
For a five-point blocking set the same equations instead force
`(n_1,n_2,n_3)=(9,1,3)`. But three distinct three-point subsets of a
five-point set cannot have pairwise intersections at most one: the
first two must meet in one point and have union all five; a third can
then contain at most two of those points without meeting one of the
first two in two points. Distinct projective lines do have pairwise
intersection one, a contradiction.

Thus a coloring with `N_0=0` splits its 13 points as six and seven.
Its total sign sum has square one. The exact line-pair count below,
applied inside that plane, gives `N_2=3`, proving (5a). Averaging it
over three-dimensional subspaces gives

\[
3p_0+p_2\ge\frac3{13}.                           \tag{5b}
\]

No independence of line colors has been assumed. To use the exact
pairwise incidence as well, let `r_pair` be the average product of the
colors of two distinct projective points. If their total color sum is
`M`,

\[
r_{\rm pair}=\frac{M^2-L}{L(L-1)}\ge-\frac1{L-1}.
\]

Every distinct pair determines exactly one projective line, so uniform
line-and-pair sampling gives exactly

\[
16p_0+4p_1
=\mathbb E_{\ell}\left(\sum_{a\in\ell}s(a)\right)^2
=4+12r_{\rm pair}.
\]

Using `p_0+p_1+p_2=1` gives `p_2=3p_0-3r_pair`. Combining with (5b),

\[
p_0\ge\frac1{26}-\frac{1}{2(L-1)}
=\rho_q.                                         \tag{6}
\]

This lower bound is positive for `q>29`. We use `q>81` so that the
square-identity coefficient `gamma_q` is also positive. It is an elementary incidence
bound, not an invocation of a large Ramsey threshold.

## 4. The Boolean nine-point test on a monochromatic line

Fix a projective line `ell`. Its four frequencies are the four lines
in a copy of `F_3^2` in the dual group. Choose any Boolean function on
`F_3^2` having five plus signs and four minus signs. Its squared mean
is `1/81`; its other Fourier weights sum to `80/81`.

Average its pullback energy over all invertible linear changes of its
two coordinates. The group `GL_2(F_3)` is transitive on the four
projective frequency lines, so each receives weight `20/81`. These
pullbacks are genuine Boolean vectors on the original group: two
linearly independent characters give a uniform quotient onto `F_3^2`.
Their averaged normalized energy is therefore exactly

\[
\frac{t+20\sum_{a\in\ell}\lambda(a)/\sqrt q}{81}
=\frac5{18}\sum_{a\in\ell}h(a)-\frac t9.
\tag{7}
\]

Its absolute value is at most `K`, because every energy being averaged
has that bound. There is no averaging over hypothetical independent
Fourier values.

Set `r(a)=K-|h(a)|>=0`. On a monochromatic line of color `sigma`, (7)
and its orientation by `sigma` imply

\[
\frac{10K}{9}-\frac5{18}\sum_{a\in\ell}r(a)
-\frac{\sigma t}{9}\le K.
\]

Hence

\[
\sum_{a\in\ell}r(a)
\ge\frac25(K-\sigma t)
\ge\frac25(K-|t|).                               \tag{8}
\]

Every projective point belongs to the same number of lines, so the
average over all lines of the left side is `4 E r`. Since `r>=0`,
(6)--(8) give

\[
\mathbb E r\ge\frac{\rho_q}{10}(K-|t|).
\tag{9}
\]

## 5. The retained spectral-defect inequality and the constant

The square deficit factors pointwise:

\[
D:=K^2-\mathbb E h^2
=\mathbb E[(K-|h|)(K+|h|)]
\ge K\mathbb E r.
\]

Substitute (4) and (9):

\[
\frac{64}{81}
\le K^2\left(1-\frac{\rho_q}{10}\right)
+\frac{\rho_q K|t|}{10}-\gamma_q t^2.
\tag{10}
\]

Completing the square in `|t|` bounds the last two terms by
`rho_q^2 K^2/(400 gamma_q)`, proving (2). As `q` tends to infinity,
the denominator factor tends to

\[
1-\frac1{260}+\frac{81}{4\cdot260^2}
=\frac{269441}{270400}.
\]

Taking the square root and dividing `K` by two gives exactly (1).

## 6. Why the second-moment-only two-dimensional route stalls

The independent exact verifier
`computations/resumed_convergence_characteristic_three_incidence_verify_2026_09_06.py`
enumerates all Boolean profiles on `F_3^2`. Their 24 distinct spectral
weight vectors have denominator 81. For a centered four-line symbol
whose entries are all `+1` or `-1`, the profile norm is `80/81` if all
four agree, and `8/9` otherwise.

In particular, choose a uniformly random `3/1` pattern and an independent
overall sign. Its four coordinates have mean zero and covariance `I`,
yet every sample has Boolean profile norm exactly `8/9`. Thus a method
using only the four spectral second moments cannot improve four-ninths.
This is a LOCAL LAW falsifier, not an actual global signing. The
13-line product identity above is exactly the global compatibility
constraint it violates if imposed at every projective plane.

The verifier also checks the exact `13 points / 13 lines / degree 4`
incidence, the parity obstruction, and `3N_0+N_2>=3` for all `2^13`
point colorings. The last inequality is sharp among individual planes:
the pairs `(N_0,N_2)=(0,3)` and `(1,0)` both occur.
These finite checks corroborate the elementary proof; no floating
integration or unproved asymptotic independence enters the theorem.
