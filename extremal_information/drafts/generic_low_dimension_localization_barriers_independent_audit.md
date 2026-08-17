# Independent audit: generic low-dimensional localization barriers

**Verdict: theorem-level PASS after minor proof and scope repairs.**

This audit freezes `generic_low_dimension_localization_barriers.md` at
SHA-256
`21effa64a5978e1bc09776692c89d0639d39b61c336256b94c1f7f0ae57a9112`.
GL.1--GL.7 are mathematically correct.  GL.6 needs one missing line to
exclude the antipodal mate when proving `D>0`; GL.1's “common positive
fraction” consequence should either state the obvious constant-generalized
version or retain the explicit `1/3` hypothesis; and several construction-
barrier phrases should be limited to the exact architectures covered.

The requested repairs were applied after the audit.  The repaired theorem
draft has SHA-256
`5c97efbe43c46d90ab060c3bc66391f868f9f469757d04433f6af303cefa35bd`:
it excludes both `c_0` and `c_0+mathbf1` in GL.6, retains the explicit
aggregate one-third hypothesis in GL.1, and limits the hierarchy statement
to genuine direct-sum nodes.

The main new result is GL.6--GL.7: for `k=o(t)=o(N)`, every deep hole has an
`o(t)` shell with projective diameter at least `(2-o(1))t`.  This is a genuine
low-dimensional theorem beyond PP.1, although it still reaches only the
covering-radius-deficit scale rather than fixed `Theta(N)` scale.

## 1. Dictionary and endpoint conventions: PASS

For `a=(-1)^y` and `z_c=(-1)^c`,

```math
\max_{c\in C}a\mathbin\cdot z_c
=N-2d(y,C).
```

Taking the minimum over `y` gives

```math
\Delta(C)=N-2\rho(C),\qquad t(C)=\Delta(C)/2.
```

The all-one word makes the response antipodal and forces
`rho<=N/2`; for odd `N`, `t` is naturally a half-integer.  Every later
strict condition `r<t` is compatible with this convention.  No parity
normalization is missing.

## 2. GL.1 direct-product activation: PASS

Covering radii and distances add under direct sums, so (GL.3) is exact.  The
mixed product word indexed by a block set `J` differs from the all-zero-choice
product in exactly

```math
w_J=\sum_{j\in J}d_j
```

coordinates.  Greedy addition reaches `N/3` because `sum d_j>=N/3`, and its
overshoot is at most `max d_j<=max N_j=o(N)`.  Hence

```math
N/3\le w_J\le N/3+o(N)<2N/3,
```

which yields the stated projective distance.  Both endpoints lie in the
radius-`rho+b` product shell.  This remains true if some selected local pair
is antipodal (`d_j=N_j`), because the global subset sum is deliberately kept
away from `0` and `N`.

Every positive-length factor containing its local all-one word has dimension
at least one, so `b<=k`.  Thus `b=o(t)` in the target scaling.

### Scope repair

The sentence saying that alternatives changing a “common positive fraction”
of every factor are covered is true but is not literally a consequence of
the displayed `N/3` hypothesis when that fraction is below `1/3`.  Either
say “at least one third in aggregate,” or state the immediate generalized
form: if `sum d_j>=alpha N`, greedy to `alpha N/2` gives projective separation
at least `min(alpha/2,1-alpha/2)N-o(N)`.  Any fixed `alpha>0` then suffices.

Likewise, the hierarchy conclusion applies at nodes whose code operation is
an actual direct sum and whose activation mass meets the hypothesis.  It does
not obstruct a coupled/non-product hierarchical extension merely because it
is described informally as nested.

## 3. GL.2 split discrepancy: PASS, including endpoints

For a length-`m`, dimension-`q` restricted code, the two-sided Rademacher
tail and a union bound over at most `2^q` words give failure probability

```math
2^q\,2\exp[-(q+2)\log2]=1/2<1.
```

Thus the constant in (GL.11) is valid.  When `m=0`, the response is
identically zero and the displayed square root is zero, so the endpoint is
covered separately without invoking a tail inequality.

If `U` is supported on `S`, it lies in the kernel of puncturing `C` to
`S^c`; rank--nullity gives complement-image dimension at most `k-ell`.
Choosing the two sign restrictions independently and using the triangle
inequality proves (GL.8).  Dividing by two via (GL.1) proves (GL.9).  No
independence between the two code images is needed.

## 4. GL.3 projective linear caps: PASS

Each projective pair `{w,w+1}` has exactly one member of weight at most `D`:
existence follows from (GL.12), while uniqueness follows from `2D<N`.
For low representatives `u,v`,

```math
\operatorname{wt}(u+v)\le2D<N-D
```

because `3D<N`; the high-lift option in (GL.12) is therefore impossible.
The low representatives form a subspace `U`, and
`W=U direct-sum <1>`.

For every coordinate in `supp U`, the coordinate functional is a nonzero
linear functional on `U`, hence equals one on exactly half its words.  Thus

```math
|\operatorname{supp}U|/2
=\mathbb E_{u\in U}\operatorname{wt}(u)\le D.
```

The strict `D<N/3` hypothesis is used exactly where needed.

## 5. GL.4 affine projective lifting: PASS

This was the most delicate quotient step, and it is correct.  Let
`bar c+bar W` lie in the shell image.  For every `bar w in bar W`, compare a
shell representative of `bar c` with one of `bar c+bar w`.  Their difference
lies in the projective class `bar w`, and its projective weight is exactly
their projective distance, at most `D`.  Therefore every projective class in
the inverse-image subspace

```math
W=\pi^{-1}(\bar W)
```

has a lift of weight at most `D`.  Here `dim W=r+1` and `1 in W`.  Since
`D=o(N)`, GL.3 applies eventually and produces an `r`-dimensional subcode
`U<=C` supported on at most `2D=o(N)` coordinates.

Applying GL.2 with `S=supp U` gives a first term

```math
O(\sqrt{Dk})=o(\sqrt{Nk})=o(N^{3/4}).
```

Because `Delta=2t=Theta(N^(3/4))`, the complement term forces
`k-r+2=Omega(sqrt N)`, hence `k-r=Omega(k)`.  Squaring the quantitative
version gives exactly

```math
k-r+2\ge(2c^2/\log2-o(1))\sqrt N;
```

the factor `2` is correct.

The diameter symbol `D` is used both for the projective cap and, in GL.2,
for support size.  Since the latter is only bounded by `2D`, this does not
affect any `o(N)` conclusion or (GL.17), but renaming it would make the proof
less error-prone.

## 6. GL.5 opposite-lift gap: PASS

For error supports of size at most `rho+s`, a pair using the large lift
satisfies

```math
N-d_{pr}=|E\triangle E'|
\le2(\rho+s)=N-2(t-s),
```

so `d_pr>=2(t-s)`.  If the entire diameter is strictly below this value,
no pair relative to `c_0` can use the large lift.  At the exact Hamming tie
`N/2`, the small and large lifts coincide, so (GL.19) remains valid.  Strict
inequality is correctly retained.

## 7. GL.6 local cap-cover entropy: PASS after one proof repair

Let `E_0=supp(y+c_0)`, so `|E_0|=rho` and `|E_0^c|=m`.  For every
`r`-set `F subseteq E_0^c`, covering-radius maximality supplies `c_F` with
`d(y+F,c_F)<=rho`.  The triangle inequality places `c_F` in `L_r(y)`.
Deepness gives `|E_F|>=rho`, and

```math
rho\ge|E_F\triangle F|
=|E_F|+r-2|E_F\cap F|
```

gives `|E_F cap F|>=ceil(r/2)=h` with the correct rounding.

Under (GL.20), GL.5 ensures that every shell word uses the ordinary lift from
`c_0`, so

```math
|G_c|\le d(c,c_0)=d_{pr}(c,c_0)\le D.
```

Because `F subseteq E_0^c`, one has the needed identity

```math
F\cap E_F=F\cap G_{c_F}.
```

Thus the `G_c` family half-covers every `r`-subset.

For an ordered sample without replacement, a prescribed set of `h` sample
positions lies in `G_c` with probability

```math
(|G_c|)_h/(m)_h\le(D/m)^h.
```

A union bound over the `binom(r,h)` position sets and then over shell words
proves (GL.21).  Taking base-two logarithms and `|L_r|<=2^k` proves (GL.22).
In fact `r<t` implies the shell contains at most one word per antipodal pair,
so `|L_r|<=2^(k-1)` yields the slightly sharper `k-1` version; the stated
weaker inequality is valid.

### Missing line in the frozen proof

The sentence “`c_F` cannot equal `c_0`, thus `D>0`” overlooks the possibility
`c_F=c_0+1`, which is different but projectively identical.  That possibility
is also impossible, using the already assumed `r<t`:

```math
d(y+F,c_0+1)
=N-(\rho+r)
=\rho+2t-r>\rho.
```

Hence `c_F` is not in the projective class of `c_0`, and `D>0` follows.
This is a one-line proof repair, not a counterexample.

All endpoint quantities are safe: `r<t<=m`, `h<=r`, and
`D<=N/2<m`, so the logarithm in (GL.22) is defined once the preceding repair
is added.

## 8. GL.7 quantifiers and constants: PASS

The hypotheses `k=o(t)` and `t->infinity` permit an integer sequence with
`k=o(r)` and `r=o(t)`; `r=ceil(sqrt(kt))` is one example.  Such `r` tends to
infinity and satisfies `1<=r<t` eventually.

If `D<2(t-r)`, GL.6 and the elementary bounds
`|L_r|<=2^k`, `binom(r,h)<=2^r` give

```math
D/m\ge2^{-(k+r)/h}.
```

Since `h=(1/2+o(1))r` and `k=o(r)`, the right side is
`1/4-o(1)`.  Also

```math
m=N/2+t=(1/2+o(1))N,
```

so this would force `D>=Theta(N)`, contradicting
`D<2(t-r)=o(N)`.  Therefore

```math
D\ge2(t-r)=(2-o(1))t.
```

The conclusion holds for **every** integer sequence with the two scale
relations, not only for one selected sequence.  In the target scaling,
`sqrt(kt)=Theta(N^(5/8))` is correct.

The constant `2` is the natural ceiling of this common-lift proof: PP.2 has
two opposite caps at projective separation `2t` and dimension comparable to
`t`, so it demonstrates the entry mechanism without contradicting GL.7.
This does not prove that every code saturating the bound actually has two
caps.

## 9. Independent finite checks

As a diagnostic independent of the proof, I enumerated every antipodal
binary linear code through length six and every deep hole.  GL.6 had 32
applicable `(code,deep-hole,r)` cases at length five and 64 at length six;
(GL.21) and `D>0` held in all 96.  I also enumerated every projectively
`D<N/3` subspace through length six; GL.3's closure, dimension, and
`|supp U|<=2D` conclusions held in every case.  These finite checks are not
used in the proofs.

## 10. Archive novelty and scope

* GL.1 is an elementary direct-sum amplification of PP.1's radius-one
  stability.  Its value is architectural: many small independently active
  blocks destroy localization before excess `t`.
* GL.2 is the standard random-sign discrepancy union bound with a new split
  application.  GL.3 is an elementary linear-anticode support lemma.  Their
  combination in GL.4 is new to the repository and rigorously eliminates a
  near-full-dimensional affine cap carrier.
* GL.5 is the two-word positivity/opposite-lift threshold already visible in
  AO.20 and PP.4, translated exactly into code radii.
* GL.6 is the substantive new generic theorem: it upgrades coordinate cover
  to a half-incidence covering-design inequality using all local `r`-flip
  obligations.
* GL.7 is a genuine unconditional low-dimensional shell-diffusion theorem at
  scale `Theta(t)`.  For the augmented-cut scaling it yields only normalized
  projective scale `Theta(N^(-1/4))` in generic `N` notation (equivalently
  `Theta(n^(-1/2))` when `N=Theta(n^2)`), still far below fixed scale.

The results do not prove that all hierarchical codes fail, that every
nonlinear shell contains a large affine subspace, or that `Theta(t)` spread
amplifies to `Theta(N)`.  The proposed nonlinear-cluster question is a
strictly narrower generic obligation, but an affirmative answer is not
suggested by GL.1--GL.7 alone.

## 11. Final classification

| Component | Verdict |
|---|---|
| GL.1 | PASS; generalize or narrow the “positive fraction” corollary |
| GL.2 | PASS with correct constants and endpoint handling |
| GL.3 | PASS; strict `D<N/3` is essential |
| GL.4 | PASS; projective affine lifting and dimension loss are valid |
| GL.5 | PASS |
| GL.6 | PASS after excluding `c_0+1` when proving `D>0` |
| GL.7 | PASS with correct quantifiers, `1/4`, and final factor `2` |
| generic fixed-scale diffusion | open exactly as stated |

No pathological or endpoint counterexample was found to the repaired
statements.  The draft contains theorem-level progress, not merely a
reformulation of PP.3, but its new unconditional scale remains `Theta(t)`.
