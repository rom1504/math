# Fractional reservoirs: finite-phase trimming and diffuse localized flips

**Status:** proved draft; no canonical edits.

Literal common-correct intersections can vanish after four positive shell
words.  This note replaces the intersection by a fractional measure which
has positive correlation with every fixed anchor.  There are two results.

1. A finite-dimensional polyhedral trimming theorem proves that such a
   diffuse fractional reservoir always exists for fixed anchor count `K`.
2. Independent Bernoulli inclusion from that reservoir gives a genuinely
   no-replacement version of the exact localized-flip argument.

Together they imply more than finite extension: every exact minimizer has a
slowly growing positive-shell packing whose pairwise separation is
`(1/4-o(1))M_n`.  The fractional constant deteriorates with the anchor
count, but it affects the number of sampled edges rather than the final
separation.  Thus the argument crosses the finite negative-holonomy ceiling
and proves a growing energy-scale packing.  It still does not change the
fixed-edge-scale lemma `L_projective`.

## 1. The finite-phase fractional-reservoir lemma

Let `r_e=(r_{1,e},...,r_{K,e})` be a finite multiset of vectors in
`{+1,-1}^K`.  Repetitions are allowed.

### Theorem FR.1 (bounded fractional trimming)

For every integer `K>=1` there is a finite constant `C_K` with the following
property.  If

```math
\sum_e r_{i,e}\ge m>0\qquad(1\le i\le K),         \tag{FR.1}
```

then there are weights `0<=w_e<=1` such that

```math
m\le W:=\sum_e w_e\le C_Km,
\qquad
\sum_e w_er_{i,e}\ge m\quad(1\le i\le K).        \tag{FR.2}
```

One may take the explicit, nonoptimal bound

```math
\boxed{
C_K=\max\left\{1,
\max_{2\le p\le K}{p^{1+p/2}\over2^{p-1}}
\right\}.}                                       \tag{FR.3}
```

For `K=1`, the inner maximum is omitted.

Equivalently, `p_e=w_e/W` is a probability distribution satisfying

```math
\max_e p_e\le {1\over m},
\qquad
\mathbb E_p r_{i,e}\ge {1\over C_K}\quad(1\le i\le K).       \tag{FR.4}
```

#### Proof

Group equal vectors.  Index the `2^K` sign patterns by `sigma`, let `n_sigma`
be their multiplicities, and let `B` be the `K`-by-`2^K` matrix whose
columns are the patterns.  After dividing by `m`, condition (FR.1) says that

```math
q=(n_\sigma/m)_\sigma\in
\mathcal Q_K:=\{x\in\mathbb R_+^{2^K}:Bx\ge\mathbf1\}.          \tag{FR.5}
```

The pointed polyhedron `\mathcal Q_K` has finitely many vertices and
recession cone

```math
\operatorname{rec}(\mathcal Q_K)
=\{d\ge0:Bd\ge0\}.
```

By the finite-dimensional Minkowski--Weyl decomposition,

```math
q=v+d,\qquad
v\in\operatorname{conv}(\operatorname{vert}\mathcal Q_K),
\quad d\in\operatorname{rec}(\mathcal Q_K).       \tag{FR.6}
```

In particular `0<=v<=q` coordinatewise and `Bv>=\mathbf1`.

It remains to bound the mass of every vertex.  If a vertex `v` has `p`
positive coordinates, then `p<=K`: otherwise the active row inequalities
cannot remove all directions supported on those coordinates.  More
precisely, one can choose `p` linearly independent active rows so that the
corresponding `p`-by-`p` sign matrix `D` is nonsingular and

```math
Dv_{\rm supp}=\mathbf1.
```

Cramer's rule expresses each coordinate as a quotient of determinants of
`p`-by-`p` sign matrices.  For `p>=2`, subtracting the first row from every
other row shows that every nonzero sign-matrix determinant is divisible by
`2^{p-1}`, hence has absolute value at least `2^{p-1}`.  Every numerator is
at most `p^{p/2}` by Hadamard's inequality.  Therefore

```math
\|v\|_1\le {p^{1+p/2}\over2^{p-1}}\le C_K.       \tag{FR.7}
```

for `p>=2`; the `p=1` bound is one.  The same bound holds for every convex
combination of vertices, in particular for the `v` in (FR.6).

Choose this `v` from (FR.6), put `x_\sigma=mv_\sigma`, and distribute total
weight `x_\sigma<=n_\sigma` arbitrarily among the copies of pattern `sigma`,
using unit weights plus at most one fractional remainder.  This gives
`0<=w_e<=1`, the upper bound in (FR.2), and all weighted row sums at least
`m`.  The lower bound `W>=m` follows because every weighted row sum is at
most `W`.  Dividing by `W` proves (FR.4). `square`

### Proposition FR.2 (the dependence on `K` cannot be uniform)

Let `K` be odd and take one copy of every sign vector having exactly
`(K+1)/2` positive coordinates.  Every such vector has coordinate sum one.
By symmetry every row sum equals

```math
m={1\over K}{K\choose(K+1)/2}.
```

If weights satisfy the `K` inequalities in (FR.2), summing those inequalities
over the rows gives

```math
W=\sum_e w_e\sum_{i=1}^K r_{i,e}\ge Km.
```

But the total available mass is exactly `Km`, so every admissible trimming
has `W=Km`.  Hence the optimal universal constant obeys `C_K^*>=K` for odd
`K`.  For even `K`, embed the construction for the first `K-1` rows and
duplicate one of those rows as the last row; then `C_K^*>=K-1`.

Thus no `K`-independent margin in (FR.4) is possible from the row-sum
hypothesis alone.  The determinant bound (FR.3) is not claimed sharp.

## 2. A no-replacement weighted localized-flip theorem

Let `A` be an exact order-`n` minimizer, put

```math
E={n\choose2},\qquad M=Q(A)=M_n,
```

and let `z_1,...,z_K` be positive augmented-cut words satisfying

```math
\langle a,z_i\rangle\ge M-2s=:m>0.               \tag{FR.8}
```

Write `r_{i,e}=a_e(z_i)_e`.  These are precisely the phase vectors to which
FR.1 applies.

### Theorem FR.3 (diffuse exact flips add a jointly separated word)

Let `C>=C_K` be a valid constant in FR.1 and let `r` be a positive integer
such that

```math
r<m,\qquad 4r<M,                                  \tag{FR.9}
```

and

```math
2^{n+1}\exp\left\{-{3r\over28C^2}\right\}
+\exp\left\{-{r\over3}\right\}<1.              \tag{FR.10}
```

Then there is a positive augmented-cut word `z_(K+1)` of deficit at most
`4r` such that, simultaneously for every `i<=K`,

```math
\boxed{
d_{\rm P}(z_{K+1},z_i)
\ge\min\left\{{m\over4},\ M-s-2r\right\}.}     \tag{FR.11}
```

The random construction below selects a set of physical edges.  No edge is
ever sampled twice, and no duplicate-removal error is hidden in (FR.11).

#### Proof

Apply FR.1 to (FR.8), obtaining weights with `m<=W<=Cm`.  Put

```math
q_e={r w_e\over W}.
```

Because `r<m<=W` and `w_e<=1`, all `q_e` lie in `[0,1]`, and
`\sum_e q_e=r`.  Include each edge independently with probability `q_e` and
call the resulting set `F`.  This is Bernoulli inclusion, not sampling with
replacement.

For every augmented cut `z`, define its weighted and sampled responses

```math
G_w(z)=\sum_e w_ea_ez_e,
\qquad
S_F(z)=\sum_{e\in F}a_ez_e.
```

Then

```math
\mathbb E S_F(z)={r\over W}G_w(z).
```

The variance sum is at most `r`, independently of `z`.  Bernstein's
inequality, with `t=r/(2C)`, gives

```math
\Pr\left\{
\left|S_F(z)-{r\over W}G_w(z)\right|>{r\over2C}
\right\}
\le2\exp\left\{-{3r\over28C^2}\right\}.         \tag{FR.12}
```

There are at most `2^n` augmented cuts.  Also,

```math
\Pr\{|F|>2r\}\le\exp\{-r/3\}.                   \tag{FR.13}
```

The union bound (FR.10) therefore supplies one set `F` for which:

1. every augmented cut obeys

   ```math
   \left|S_F(z)-{r\over W}G_w(z)\right|\le {r\over2C};          \tag{FR.14}
   ```

2. `|F|<=2r`.

Flip the edges of this fixed `F`.  Since `A` is an exact minimizer,
`Q(A^F)>=M`.  Orient an augmented cut `z` so that

```math
\langle a^F,z\rangle\ge M.
```

Its response under the original signing obeys

```math
H:=\langle a,z\rangle
=\langle a^F,z\rangle+2\sum_{e\in F}a_ez_e
\ge M-2|F|\ge M-4r>0.                            \tag{FR.15}
```

Thus `z` is positive for `A`, and its deficit `d=M-H` is at most `4r`.
Moreover

```math
S_F(z)=\sum_{e\in F}a_ez_e\le-{d\over2}\le0.    \tag{FR.16}
```

Apply the uniform approximation (FR.14) to this response.  Since
`S_F(z)<=0`,

```math
{r\over W}G_w(z)\le {r\over2C},
\qquad\text{hence}\qquad
G_w(z)\le {W\over2C}\le {m\over2}.              \tag{FR.17}
```

On the other hand, FR.1 gives `G_w(z_i)>=m` for every old anchor.  If
`h_i=d_H(z,z_i)`, then

```math
{m\over2}
\le G_w(z_i)-G_w(z)
=\sum_e w_ea_e\big((z_i)_e-z_e\big)
\le2\sum_{e:z_e\ne(z_i)_e}w_e
\le2h_i.
```

Thus `h_i>=m/4`.  Notice that `C` has disappeared from the separation; it
controls only how large `r` must be for uniform approximation.

Both words are positive under `A`, so the elementary two-response bound
also gives

```math
E-h_i\ge{H+\langle a,z_i\rangle\over2}
\ge M-s-2r.                                      \tag{FR.18}
```

Taking the minimum of `h_i` and `E-h_i` proves (FR.11). `square`

### Constants in the sampling step

For completeness, (FR.12) uses the two-sided Bernstein bound

```math
\Pr\{|X-\mathbb EX|>t\}
\le2\exp\left\{-{t^2\over2(r+t/3)}\right\}
```

with `t=r/(2C)` and `C>=1`; its exponent is at least `3r/(28C^2)`.
Independent inclusion makes `F` a genuine set.  In contrast, `r`
independent draws from `p_e=w_e/W` would have nonnegligible collisions at
the relevant `r=Theta_K(n)` scale and is not used.

## 3. Fixed and slowly growing packings

### Corollary FR.4 (arbitrarily large fixed energy-scale packings)

For every fixed integer `L`, there is a deterministic
`D_{n,L}=O_L(n)=o(M_n)` such that the positive `D_{n,L}`-deficit shell of
every exact order-`n` minimizer contains `L` words at pairwise projective
distance at least

```math
{M_n-D_{n,L}\over4}=(1/4-o(1))M_n.
```

With

```math
\widehat C_L=\max\left\{1,
\max_{2\le p\le L}{p^{1+p/2}\over2^{p-1}}
\right\},
\qquad
r_{n,L}=\left\lceil
100\widehat C_L^2\big((n+1)\log2+1\big)
\right\rceil,                                    \tag{FR.19}
```

one may take

```math
D_{n,L}=4r_{n,L}.                                 \tag{FR.20}
```

for all sufficiently large `n`.

#### Proof

Start with one positive ground word.  Always apply FR.3 to the same base
signing `A`, never to a previously flipped signing.  At every stage use the
common shell parameter `s=2r_{n,L}` and the harmless common constant
`C=\widehat C_L`.  Every old word has deficit at most `4r_{n,L}`, and every
new word does as well.  Since `L` is fixed and `M_n=Theta(n^(3/2))`, the
conditions (FR.9)--(FR.10) hold eventually.  Equation (FR.11) is then at
least `(M_n-4r_{n,L})/4` for every new-old pair.  Previously established
distances are unchanged.  After `L-1` applications the result follows.
`square`

This theorem has no common-intersection premise.  In particular, the
negative-holonomy quartet from TL.4 does not stop the next application:
the fractional measure is allowed to mix several phase cells while paying
all anchors jointly.

### Corollary FR.5 (a legitimate but nonuniform diagonal)

Fix any `0<alpha<1/2` and put

```math
L_n=\left\lfloor{\alpha\log n\over\log\log n}\right\rfloor.
```

Using (FR.19) with `L=L_n` gives

```math
\widehat C_{L_n}\le n^{\alpha/2+o(1)},
\qquad
D_{n,L_n}\le n^{1+\alpha+o(1)}=o(M_n).            \tag{FR.21}
```

Therefore every exact minimizer has a positive `o(M_n)` shell containing
`L_n\to\infty` words with pairwise distance at least

```math
{M_n-D_{n,L_n}\over4}
=(1/4-o(1))M_n.                                   \tag{FR.22}
```

This diagonal is mathematically valid because all constants and union
bounds were kept uniform up to `L_n`.  It does **not** justify replacing a
fixed `L` by an arbitrarily faster-growing sequence.  Proposition FR.2
shows that the fractional margin must deteriorate at least as `1/L` under
the bare row-sum hypothesis, and uniform response approximation costs order
`C_L^2n` flipped edges.  Crucially, once that approximation cost remains
`o(M_n)`, the final separation in (FR.22) does not deteriorate with `L`.

## 4. Frontier effect

The finite negative-holonomy ceiling is crossed: a slowly growing number of
anchors can be constructed, without a literal common-correct coordinate and
without accumulating shell deficit, at uniform `(1/4-o(1))M_n` separation.
This strictly strengthens the four-witness theorem TL.2 and completes the
previously missing **growing energy-scale packing** step.

The selected `L_projective` still does not follow.  Even a constant multiple
of `M_n=Theta(n^(3/2))` is a vanishing fraction of `E=Theta(n^2)`.  The
result therefore supplies a fractional-reservoir theorem and a genuine
growing mesoscopic packing, not a fixed-edge-scale compiler.

The next question is whether exact-minimizer cut structure improves the
fractional constant beyond the arbitrary phase-pattern lower bound in FR.2,
allowing a faster-growing packing, or whether a different amplifier can
convert the growing but mesoscopic packing (FR.22) into fixed-edge-scale
contextual information.
