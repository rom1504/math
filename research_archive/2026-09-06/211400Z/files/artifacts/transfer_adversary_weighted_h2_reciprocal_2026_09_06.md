# A concrete reciprocal-spectrum test for the weighted H2 seed

Date: 2026-09-06. New exact saturation identities and finite falsifiers.
No asymptotic strict `R<T` gap is proved. This isolates a specific binary
array-pair constraint for the smallest unresolved weighted seed, beyond
the generic request to realize an orthogonal polar action.

## 1. The seed, majorant, and exact defect

Use

```math
B=\begin{pmatrix}-1&2\\2&4\end{pmatrix},\quad
P=\sqrt2\begin{pmatrix}1&0\\0&4\end{pmatrix},\quad
M=P^{-1/2}BP^{-1/2}
=\frac1{\sqrt2}\begin{pmatrix}-1&1\\1&1\end{pmatrix}.
```

Then `M^2=I`, so `P` majorizes both signs of `B`. This proves
`T(B)<=5/sqrt2`. Equality can be checked exactly with the cut covariance
`C=[[1,3/4],[3/4,1]]`: direct multiplication gives `BCB=PCP`.
Consequently `|C^(1/2) B C^(1/2)|=C^(1/2) P C^(1/2)` and its nuclear
norm is `tr PC=5sqrt2`. The cut-covariance lower argument in the seed-loss
artifact gives `T(B)=5/sqrt2`. The original same-spin unamplified value is
`q(B)=7/2`, since its two possible quadratic values are `7/2` and `-1/2`.

Let `U=H_s/sqrt(s)` be the standard Walsh transform, with the normalized
counting inner product. For Boolean input functions `f_1,f_2`, put
`a=Uf_1`, `b=Uf_2`. For Boolean output functions `g_1,g_2`, define

```math
D=\mathbb E\left(a-\frac{-g_1+2g_2}{\sqrt2}\right)^2
 +\mathbb E\left(2b-\frac{g_1+2g_2}{\sqrt2}\right)^2.
```

Both vector fields in this squared distance have exact squared norm five.
Expanding the square therefore gives the identity

```math
\frac{\beta(H_s\otimes B)}{2s^{3/2}}
=\frac5{\sqrt2}-\frac1{2\sqrt2}\min_{f_1,f_2,g_1,g_2}D.       (1)
```

Changing all output signs accommodates the absolute value in `beta`.
Independent row/column signing equivalence of the standard even-dimensional
Walsh matrix and the regular order-four generator makes (1) applicable to
the prescribed `R` at `s=4^r`.

## 2. Saturation forces reciprocal spectra at correlation three quarters

If `D=0`, then pointwise

```math
a=\frac{-g_1+2g_2}{\sqrt2},\qquad
b=\frac{g_1+2g_2}{2\sqrt2},\qquad ab=\frac34.                (2)
```

Parseval gives `E a^2=E b^2=1`. Equation (2) then forces
`E g_1g_2=3/4`, and `E f_1f_2=E ab=3/4`. More explicitly, the large
absolute value of `a` occurs on `1/8` of the coordinates:

```math
|a|\in\{1/\sqrt2,3/\sqrt2\},\qquad
|b|\in\{1/(2\sqrt2),3/(2\sqrt2)\},
```

with reciprocal matching of the two levels and the same spectral sign.

Equivalently the periodic cross-correlations of the Boolean inputs on the
elementary abelian group satisfy

```math
\frac1s\sum_x f_1(x)f_2(x+t)=\frac34\,1_{\{t=0\}}.
```

Thus saturation requires a nonzero perfect binary array pair with prescribed
in-phase correlation AND these two amplitude levels, on an elementary
abelian 2-group. A scalar bent-function theorem does not supply this pair.

At an allowed finite order `s=4^r`, equality is already arithmetically
impossible: `U` has rational entries, whereas every nonzero target coordinate
in (2) is irrational. This is only finite nonattainment; the dyadic lattice
becomes dense and gives no uniform asymptotic gap.

## 3. The stable consequence uses L1, not a fourth-moment shortcut

Put `a_0=(-g_1+2g_2)/sqrt2` and
`b_0=(g_1+2g_2)/(2sqrt2)`, so `a_0b_0=3/4` pointwise. Since
`||a||_2=1`, `||b_0||_infinity<=3/(2sqrt2)`,
`||a-a_0||_2<=sqrt D`, and `||b-b_0||_2<=sqrt D/2`,

```math
\mathbb E|ab-3/4|
\le\left(\frac12+\frac3{2\sqrt2}\right)\sqrt D.              (3)
```

In particular the input correlation tends to `3/4` whenever `D->0`.
The output correlation obeys
`|E g_1g_2-3/4|<=sqrt D+D/2`, by comparing `||a_0||_2^2` with one.

Equation (3) is an actual necessary stability constraint for `R(B)=T(B)`.
It does NOT imply convergence of `E a^2b^2` to `9/16`; small L2 errors can
have rare large fourth-moment spikes. Any attempt to use a lower bound on
that fourth-order quantity alone must separately control this issue.

## 4. Exhaustive finite probe

`computations/transfer_adversary_reciprocal_pair_probe_2026_09_06.py` enumerates
every pair of Boolean functions with correlation exactly `3/4`, up to their
common reversal. All Fourier operations and comparisons are integer exact.

| Walsh order | Pairs examined | Minimum `E a^2b^2` | Minimum `E|ab-3/4|` |
| --- | ---: | ---: | ---: |
| 8 | 1,024 | `3/2` | `3/4` |
| 16 | 3,932,160 | `15/16` | `3/8` |

Order eight is an odd-dimensional Walsh diagnostic, not itself an allowed
`H4` catalyst. Order sixteen is allowed. Neither finite floor is asserted
uniform in order or in correlations approaching, but not exactly, `3/4`.

## 5. Exact primary-literature mapping and limitation

Arasu, Goyal and Puri,
[Binary sequence/array pairs via difference set pairs: A recursive approach,
Transactions on Combinatorics 6 (2017), 19--36](https://toc.ui.ac.ir/article_21466_7872b534dc27cd9c3fa2ab7a0cb15ba8.pdf),
Definitions (1.2), (1.5)--(1.6), and Proposition 2.2 identify the cross-perfect
condition with a difference-set pair in the underlying abelian group.
The theorem's group must be matched, not only its cardinality.

Their Theorem 3.2 gives balanced pairs with zero off-shift correlations in
every even-order abelian group. But its in-phase correlation is ALSO zero:
its parameters have `e=lambda`, and the in-phase minus off-phase difference
is `4(e-lambda)`. Thus that apparently broad existence theorem does not
supply (2). The later displayed recursive constructions multiply by odd
groups of orders `4m-1` or `4m+1`; they do not land in growing elementary
abelian 2-groups. No nonzero reciprocal pair satisfying (2), and no uniform
stable nonexistence theorem, has been imported from this source.

## 6. A larger exact witness defeats extrapolation of the order-sixteen floor

The bounded search
`computations/transfer_adversary_reciprocal_search_2026_09_06.py`
at order 64, batch 256, 10,000 iterations, seed 260906 found a pair
with exact input correlation `3/4` and

```math
\mathbb E |(Uf)(Ug)-3/4|=1312/4096=41/128<3/8.
```

Its final integer recomputation verifies both equalities. Thus the exact
order-sixteen minimum `3/8` is NOT a dimension-free lower bound. The
search is reproducible and prints the full witness and its integer Walsh
spectra. No optimality claim at order 64 is made. This witness does not
match the two-amplitude palette and does not imply a new lower bound on
the weighted seed's `R`.

A separate direct weighted-seed search at order 64, batch 256, 100 cycles,
seed 260906, found no score above the known `7/2`; this is only a bounded
search observation, not a certificate.

The exact elementary-2-group terminology occurs in Mao, Jiang, Zhao and
Zhou, *Research on perfect dyadic binary sequence pair*, Journal of
Electronics (China) 23 (2006), 361--364,
[public author manuscript](https://www.researchgate.net/publication/245462979_Research_on_perfect_dyadic_binary_sequence_pair),
DOI `10.1007/s11767-004-0163-6`. Their Theorem 3 gives the same Fourier
product identity, with unnormalized Walsh product `N-2d`, where `d` is
Hamming distance. The source supplies divisibility/weight conditions and
elementary sign/shift transformations. Its displayed order-sixteen pair
has `g=-f` and every squared Walsh coefficient of `f` equals sixteen, so
it only supplies correlation `-1`, not a nontrivial intermediate value.
No quantitative stable nonexistence result is asserted from this source.
