# Independent audit: fractional reservoirs and localized flips

**Verdict:** PASS.

**Audited source:** `extremal_information/drafts/fractional_reservoir_localized_flip.md`,
SHA-256
`45628880b288e7f94aec91beba438f6f4f4fa8349499e0bc2e5cb76d92a49b0b`.

This audit checks FR.1--FR.5 independently.  In particular it checks the
uniform Bernstein event before the optimizing response is selected, the
exact-minimality sign conversion, both lifts in projective distance, the
nonaccumulating iteration, and the growing diagonal.  No canonical file is
edited here.

## 1. FR.1 polyhedral trimming: PASS

After grouping equal sign types, write `A` for the phase matrix and `n` for
the capacity vector.  The normalized capacity belongs to

```math
\mathcal Q_K=\{x\ge0:Ax\ge\mathbf1\}.
```

This polyhedron is pointed and has recession cone

```math
\{d\ge0:Ad\ge0\}.
```

Minkowski--Weyl therefore writes `n/m=v+d` with `v` a convex combination
of vertices and `d` in the recession cone.  Hence

```math
0\le mv\le n,\qquad A(mv)\ge m\mathbf1.
```

If a vertex has `p` positive coordinates, its zero-coordinate facets leave
`p` support directions, so `p` independent active row constraints are
required.  Thus `p\le K` and its positive coordinates solve

```math
Dv_{\rm supp}=\mathbf1
```

for an invertible `p`-by-`p` sign matrix `D`.  For `p\ge2`,
`|\det D|\ge2^{p-1}`: subtract the first row from the other rows and factor
`2` from each.  Every Cramer numerator is another sign determinant and is
at most `p^{p/2}` by Hadamard.  Therefore

```math
\|v\|_1\le {p^{1+p/2}\over2^{p-1}}.
```

The `p=1` value is one.  Convex combinations preserve the mass bound.
Distributing each aggregate mass among its occurrences gives weights in
`[0,1]`.  Finally every weighted row sum is at most the total mass `W`, so
`W\ge m`.  All constants and inequality directions in FR.1 are correct.

The exact-arithmetic checker
`experiments/finite_phase_fractional_reservoir_constants.py` independently
confirms the basic-support calculation through `K=5`.  This finite check is
not needed for the theorem.

## 2. Bernoulli inclusion and FR.12: PASS

FR.1 gives

```math
m\le W\le Cm,\qquad
G_w(z_i)=\sum_e w_ea_e(z_i)_e\ge m.
```

With `q_e=rw_e/W`, the hypotheses `r<m\le W` and `w_e\le1` imply
`0\le q_e<1`, while

```math
\sum_eq_e=r.
```

Thus independent Bernoulli inclusion really produces a set; no edge can be
sampled twice.

For each fixed augmented cut `z`,

```math
\mathbb ES_F(z)={r\over W}G_w(z).
```

Every centered summand has absolute value at most one, and the variance sum
is at most `\sum_eq_e=r`.  Two-sided Bernstein with
`t=r/(2C)` gives exponent

```math
{t^2\over2(r+t/3)}
={3r\over4C(6C+1)}
\ge {3r\over28C^2},
```

where the final inequality is exactly `C\ge1`.  Hence FR.12 has the right
constant and its factor two is correct.

There are at most `2^n` augmented-cut words, so the uniform-response failure
probability is at most

```math
2^{n+1}\exp\{-3r/(28C^2)\}.
```

Also `\mathbb E|F|=r`, and the standard multiplicative Chernoff bound at
factor two gives

```math
\Pr\{|F|>2r\}\le(e/4)^r\le e^{-r/3}.
```

Thus FR.10 supplies one fixed set satisfying FR.14 for **every** augmented
cut and `|F|\le2r`.  This uniformity is essential because the later
maximizing response depends on `F`; the proof has it.

## 3. Exact-minimality conversion: PASS

For the signing obtained by flipping `F`,

```math
\langle a^F,z\rangle
=\langle a,z\rangle-2S_F(z).
```

Exact order-`n` minimality gives `Q(A^F)\ge M`.  Orient a maximizing
augmented cut so that `\langle a^F,z\rangle\ge M` and put
`H=\langle a,z\rangle`.  Since `|F|\le2r`,

```math
H=\langle a^F,z\rangle+2S_F(z)
\ge M-2|F|\ge M-4r>0.
```

The response is therefore positive for the original signing.  Since
`Q(A)=M`, its deficit `d=M-H` is nonnegative and at most `4r`.  Rearranging
the flipped-response inequality gives

```math
S_F(z)\le-{d\over2}\le0,
```

which is FR.16 with the correct sign.

Applying the already-uniform FR.14 to this adaptively selected `z` gives

```math
{r\over W}G_w(z)\le {r\over2C},
\qquad
G_w(z)\le{W\over2C}\le{m\over2}.
```

No independence between `z` and `F` is assumed.

## 4. Signed and projective separation: PASS

For each old anchor,

```math
G_w(z_i)-G_w(z)\ge {m\over2}.
```

Only mismatching coordinates contribute to this difference, and their
individual contributions are at most `2w_e`.  Since `w_e\le1`,

```math
{m\over2}
\le2\sum_{e:z_e\ne(z_i)_e}w_e
\le2d_H(z,z_i).
```

Thus the actual Hamming distance obeys the `K`- and `C`-independent bound
`h_i\ge m/4`.  The constant `C` affects the sampling cost, not this final
distance.

For the complementary lift, the two response vectors agree on exactly
`E-h_i` coordinates.  Hence

```math
H+\langle a,z_i\rangle\le2(E-h_i).
```

Using `H\ge M-4r` and
`\langle a,z_i\rangle\ge M-2s` yields

```math
E-h_i\ge M-s-2r.
```

Taking `\min\{h_i,E-h_i\}` proves FR.11 exactly.  Positivity, the factor
`1/4`, and the complementary-distance constant are all correct.

## 5. Fixed-`L` iteration: PASS

At every stage FR.4 uses the same base exact minimizer and the same sample
size `r=r_{n,L}`.  Setting the common FR.3 shell parameter to

```math
s=2r
```

is correct because all previously generated words have deficit at most
`4r=2s`.  Thus

```math
m=M-2s=M-4r.
```

The first term of FR.11 is

```math
{m\over4}={M-4r\over4},
```

while the complementary term is `M-s-2r=M-4r` and is larger.  Every new
word again has deficit at most `4r`.  Since each perturbation starts from
the original `A`, shell width does not accumulate.

For

```math
r_{n,L}=
\left\lceil100\widehat C_L^2((n+1)\log2+1)\right\rceil,
```

the first exponent in FR.10 is at least
`(300/28)((n+1)\log2+1)`.  This dominates the `2^{n+1}` union factor.
For fixed `L`, `r=O_L(n)=o(M_n)`, so `r<m` and `4r<M_n` eventually.
Using `C=\widehat C_L\ge C_K` at every intermediate anchor count
`K\le L` is valid.  FR.4 follows after `L-1` finite applications.

## 6. Growing diagonal: PASS

Let

```math
L_n=\left\lfloor{\alpha\log n\over\log\log n}\right\rfloor,
\qquad0<\alpha<1/2.
```

The explicit determinant constant obeys

```math
\log\widehat C_{L_n}
\le\left(1+{L_n\over2}\right)\log L_n
=\left({\alpha\over2}+o(1)\right)\log n.
```

Conversely, taking `p=L_n` in the defining maximum gives

```math
\log\widehat C_{L_n}
\ge\left(1+{L_n\over2}\right)\log L_n-(L_n-1)\log2
=\left({\alpha\over2}+o(1)\right)\log n.
```

Therefore the exponent notation in the source is exact:

```math
\widehat C_{L_n}\le n^{\alpha/2+o(1)},\qquad
r_{n,L_n}=n^{1+\alpha+o(1)}.
```

Because `\alpha<1/2` and `M_n=\Theta(n^{3/2})`,
`r_{n,L_n}=o(M_n)`.  All FR.3 positivity conditions are consequently
uniform for every stage `K\le L_n`.  FR.10 is also uniform because its
choice of `r` already absorbs `\widehat C_{L_n}^2`; there is no additional
union bound over the sequential existence applications.

With `D_{n,L_n}=4r_{n,L_n}=o(M_n)`, the common pairwise lower bound is

```math
{M_n-D_{n,L_n}\over4}=(1/4-o(1))M_n,
```

and `L_n\to\infty`.  The strict threshold `\alpha<1/2` is exactly what
makes the flip budget smaller than the `n^{3/2}` energy scale.

## 7. Scope and source hygiene

| Item | Verdict |
|---|---|
| FR.1 Minkowski--Weyl domination | PASS |
| determinant divisibility and equation FR.3 | PASS |
| Bernoulli no-replacement construction | PASS |
| uniform FR.12 Bernstein bound | PASS |
| FR.10 union bound | PASS |
| exact-minimality response conversion | PASS |
| `m/4` signed distance | PASS |
| complementary/projective distance | PASS |
| fixed-`L` iteration with `s=2r` | PASS |
| `L_n` diagonal for `\alpha<1/2` | PASS |
| growing energy-scale packing | PASS |
| implication to fixed-scale `L_projective` | correctly denied |

The theorem supplies a growing packing with guaranteed separation
`\Omega(M_n)=\Omega(n^{3/2})`, while the guaranteed scale remains a
vanishing fraction of
`E=\Theta(n^2)`.  It therefore advances the mesoscopic frontier but does
not prove `L_projective`.

No theorem-level typo or missing hypothesis remains in the frozen source.
The only cosmetic issue is a line-wrap break in the prose following FR.7;
it does not affect rendering or meaning.
