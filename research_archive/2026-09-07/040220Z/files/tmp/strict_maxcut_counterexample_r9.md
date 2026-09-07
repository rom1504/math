# Strict weighted failure of the maximum-cut averaging target

This is an exact audit of (10.442), using undoubled edge weights throughout.
The verification script is `tmp/audit_strict_maxcut_r9.py`.

## Witness

Let `S={0,1,2,3}` and `T={4,5,6,7,8}`.  The six edges inside `S` all
have weight `9/112`.  Every row of the `S`-to-`T` block is

```text
(1/24, 1/24, 1/24, 55/672, 29/672).
```

Inside `T`, the triangle on `{4,5,6}` has edge weight `83/672`, each
of the six edges from `{4,5,6}` to `{7,8}` has weight `-55/672`, and
edge `{7,8}` has weight `113/672`.

Exact enumeration of all 256 projective parent cuts gives `c=1`, with
the all-one cut uniquely of weight zero and `S|T` uniquely of weight one.
Every other projective cut lies in `[1/112,111/112]`; hence both parent
endpoints are strict.

The child data are

```math
(a_S,m_S,M_S)=\left(\frac{27}{56},0,\frac9{28}\right),
\qquad
(a_T,m_T,M_T)=\left(\frac1{21},-\frac{55}{112},\frac{19}{112}\right).
```

Both children are strictly positive-dominant, since

```math
a_S-m_S-M_S=\frac9{56},
\qquad
a_T-m_T-M_T=\frac{31}{84}.
```

## Every maximum-cut pair

After fixing the first coordinate on each shore to `+1`, all maximum-cut
spins on `S` are

```text
(1,-1,-1, 1)
(1,-1, 1,-1)
(1, 1,-1,-1)
```

and all maximum-cut spins on `T` are

```text
(1,-1,-1,-1, 1)
(1,-1,-1, 1,-1)
(1,-1, 1,-1, 1)
(1,-1, 1, 1,-1)
(1, 1,-1,-1, 1)
(1, 1,-1, 1,-1).
```

There are therefore 18 projective maximizing pairs.  For every one of
them,

```math
z=x^{\mathsf T}By=0.
```

This also follows without checking the 18 pairs: all four rows of `B` are
identical, while each listed `S` maximum-cut spin is balanced, so
`x^{\mathsf T}B=0`.

Consequently the best possible left side of (10.442) is

```math
2(M_S+M_T)+|z|=\frac{55}{56},
```

whereas its right side is

```math
a_S+a_T-(m_S+m_T)=\frac{49}{48}.
```

Thus every maximum-cut pair fails (10.442), by the exact undoubled gap

```math
\frac{49}{48}-\frac{55}{56}=\frac{13}{336}.
```

As a factor-of-two check, the favorable relative sign of a maximum-cut
pair has doubled parent energy

```math
2(a_S+a_T)-4(M_S+M_T)-2|z|=-\frac{19}{21},
```

while the desired comparison level is

```math
2(m_S+m_T)=-\frac{55}{56}.
```

The former is larger by `13/168`, exactly twice the undoubled failure.

## Finite `\{\pm1\}` scope

The separate script `tmp/audit_maxcut_averaging_r9.py` exhausts switching
gauges and every actual parent endpoint pair, then enters precisely the
same-positive hard case and enumerates all child maximum-cut pairs.  It
finds:

- order six: all 1,024 switching classes pass, covering 322 hard endpoint
  pairs; the minimum undoubled slack is 3;
- order seven: all 32,768 switching classes pass, covering 73,297 hard
  endpoint pairs; the minimum undoubled slack is 1.

Random audits of 2,000 switching-gauge draws at each order 8 through 12 also
found no failure.  These random counts are numerical evidence only.

The weighted witness shows that no argument using only the parent cut-cone
constraints, child extremality, and positive dominance can prove (10.442).
Its missing finite-size condition is that an actual signing has every atomic
edge weight in `\{\pm1\}`.  That condition does not rescue the universal
stopping inequality asymptotically: the independently audited strict weighted
witness can be lifted by a biased random clone construction to complete sign
matrices while preserving a macroscopic failure of (10.435).

In fact the same lift preserves failure of (10.442).  Clone each of the nine
vertices `L` times and sample each clone edge independently with mean equal to
the displayed strict weight (using mean zero within a clone class).  Uniform
cut concentration gives `O(L^{3/2})=o(L^2)` error in all child cut extrema and
in every cross bilinear form.  Any sequence of actual child maximum cuts has,
after taking clone-density subsequences, a limit in the maximum face of the
corresponding multilinear child cut function.  Such a limiting density is the
mean of a distribution supported on Boolean maximum cuts.  Every Boolean
maximum spin on `S` is balanced, so every point of that maximum face has
coordinate sum zero.  Since the four rows of the mean cross block coincide,
its bilinear form vanishes against the entire `S` maximum face.  Consequently
the largest `|z|` over actual maximum-cut pairs is `o(L^2)`, while

```math
\bigl(a_S+a_T-m_S-m_T-2M_S-2M_T\bigr)L^2
=\frac{13}{336}L^2.
```

Thus sufficiently large complete sign-matrix blow-ups also violate (10.442)
for every actual pair of child maximum cuts, not merely its weighted
relaxation.
