# Radial shell sufficiency and half-covering at a near-minimizer

Date: 2026-08-17.

Status: proof draft.  The theorem is an exact consequence of global
near-minimality.  It identifies which part of the augmented-cut landscape
answers every bounded Hamming edit, but it does not bound the size of that
part and therefore is not yet a composable state.

## 1. Setup

Let `E=binom([n],2)` and write a hollow signing as
`a in {+-1}^E`.  Let

```math
\mathcal Z_n=
\{\sigma(x_ix_j)_(i<j):x\in\{+-1\}^n,\ \sigma\in\{+-1\}\}
```

be the augmented cut code in sign coordinates.  Thus

```math
Q(a)=\max_(z\in\mathcal Z_n)\langle a,z\rangle.
```

For `F subseteq E`, let `a^F` be obtained by reversing the signs on `F`,
put `r=|F|`, and define

```math
d_a(z)=Q(a)-\langle a,z\rangle,
\qquad
s_F(z)=\sum_(e\in F)a_ez_e.
```

The exact radial response identity is

```math
Q(a^F)-Q(a)
=\max_(z\in\mathcal Z_n)\{-d_a(z)-2s_F(z)\}.       \tag{RS.1}
```

## 2. Thin-shell sufficiency

### Theorem RS.1 (near-minimal radial shell theorem)

Suppose

```math
Q(a)\le M_n+\eta.
```

For every `F subseteq E` of size `r`, an optimizer `z_F` for `Q(a^F)` can
be chosen with

```math
d_a(z_F)\le \eta+2r                                      \tag{RS.2}
```

and

```math
s_F(z_F)\le {\eta-d_a(z_F)\over2}\le {\eta\over2}.       \tag{RS.3}
```

Consequently, on the entire Hamming ball `|F|<=r`, formula (RS.1) is
unchanged if its maximum is restricted to the one common shell

```math
\mathcal S_a(\eta+2r)
=\{z\in\mathcal Z_n:d_a(z)\le\eta+2r\}.                  \tag{RS.4}
```

Equivalently, if

```math
D_z=\{e:a_ez_e=-1\},
```

then the shell disagreement sets satisfy the threshold-covering law

```math
\boxed{
\forall F\in\binom E r\quad
\exists z\in\mathcal S_a(\eta+2r):
|F\cap D_z|
\ge {r\over2}-{\eta\over4}+{d_a(z)\over4}.}              \tag{RS.5}
```

For an exact minimizer (`eta=0`), every `r`-edge direction therefore has a
`2r`-shell witness disagreeing with at least half of its edited edges.

#### Proof

Since `Q(a^F)>=M_n>=Q(a)-eta`, choose an optimizer `z_F` and use

```math
Q(a^F)=\langle a,z_F\rangle-2s_F(z_F)
=Q(a)-d_a(z_F)-2s_F(z_F).
```

It follows that

```math
d_a(z_F)+2s_F(z_F)\le\eta.                                \tag{RS.6}
```

As `s_F(z_F)>=-r`, this proves (RS.2), while (RS.3) is (RS.6)
rearranged.  A term with `d_a(z)>eta+2r` is strictly smaller than `-eta`
in (RS.1), whereas an optimizer has value at least `-eta`; hence the same
shell suffices simultaneously for every `|F|<=r` (using the larger shell
at radius `r`).  Finally

```math
s_F(z)=r-2|F\cap D_z|,
```

so (RS.3) gives (RS.5). `square`

## 3. What the theorem does and does not reduce

RS.1 is stronger than a scalar shell-cardinality statement: it gives an
exact response formula for all edit contexts in a Hamming ball and a
uniform covering obligation on the same shell.  It is, however, an exact
repackaging of the previously proved global flip inequality in
`nearmin_deterministic_inequalities.md`, not a new frontier theorem.  A thin
energy window need not have a small description or small response entropy.

However, the theorem supplies neither a small presentation of (RS.4) nor a
congruence under block composition.  In fact, retaining all affine
functionals in (RS.1) is precisely the local upper-response roof.  Calling
that roof a state would merely relocate the optimization unless its metric
entropy or algebraic structure can be bounded independently.

Thus RS.1 is currently:

* a proved Level-5 necessary law;
* a clean bridge from near-minimality to a declared family of future
  contexts;
* a target for lower/upper bounds on shell response entropy;
* **not** by itself the missing transfer lemma.

There is also a sharp scale warning.  For every signing and every
`|F|<=r`, entrywise Lipschitzness gives

```math
|Q(a^F)-Q(a)|\le2r.                                      \tag{RS.7}
```

Hence the zero-information predictor `Q(a)` already answers the whole
radius-`r` context class to `o(n^(3/2))` error whenever
`r=o(n^(3/2))`.  To make radial edit responses differ on a fixed positive
`n^(3/2)` scale one needs `r=Omega(n^(3/2))`; at exactly that scale the shell
width `eta+2r` in (RS.4) is no longer `o(n^(3/2))`.  Thus this particular
local-context route has a scale incompatibility: its thin-shell regime is
automatically trivial at the target response accuracy, while its
target-response regime does not force a vanishing shell.  A useful transfer
would need amplification by a low-cap physical continuation or a different
query class, not merely sharper metric entropy for sparse edge edits.

The next discriminating question is whether the shell in (RS.4), for some
`r` large enough to matter under a physical compiler but still
`eta+2r=o(n^(3/2))`, has a sub-landscape response quotient forced by
near-minimality.  Any proposed quotient must be tested against the universal
low-local-field affine cubes and the sparse-flip contextual packing.
