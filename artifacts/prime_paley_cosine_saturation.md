# Prime-Paley conference matrices have a saturating subsequence

Date: 2026-07-31. This is an agent-authored non-tensor design-family audit.
It proves an asymptotic theorem; it is not an extrapolation from the exact
orders 6, 14, and 18.

## 1. Exact Paley/Fourier mapping

Let `p=1 mod 4` be prime, let `chi` be the quadratic character of `F_p`, and
let `C_p` be the symmetric Paley conference matrix of order `p+1`. Its finite
core is

```math
(L_p)_{xy}=\chi(x-y),qquad x,y\in F_p,                \tag{PP1}
```

with `chi(0)=0`; the infinity row and column are all `+1`. For a Boolean
finite spin `f` and infinity spin `s`,

```math
H_{C_p}(s,f)=s\sum_x f(x)+{1\over2}f^{\mathsf T}L_pf. \tag{PP2}
```

Use the Fourier convention

```math
\widehat f(a)=\sum_{x\in F_p}f(x)e^{-2\pi iax/p}.
```

The quadratic Gauss-sum identity diagonalizes (PP1): for `a != 0`, its
additive character has eigenvalue `sqrt(p) chi(a)`. Therefore

```math
{f^{\mathsf T}L_pf\over p\sqrt p}
={1\over p^2}\sum_{a\ne0}\chi(a)|\widehat f(a)|^2.    \tag{PP3}
```

This checks the exact normalization against the project's Hamiltonian. The
right side is the difference between the fractions of Fourier energy on
quadratic-residue and nonresidue frequencies.

## 2. Cosine-threshold spin

Define

```math
f_p(x)=\operatorname{sign}\cos(2\pi x/p).             \tag{PP4}
```

No cosine vanishes for odd `p`. For `p=1 mod 4`, exactly `(p+1)/2` values are
positive, so

```math
\sum_xf_p(x)=1.                                       \tag{PP5}
```

For every fixed integer `a`, the Riemann-sum limit is

```math
{1\over p}\widehat f_p(a)
\longrightarrow
\int_0^1\operatorname{sign}(\cos2\pi t)e^{-2\pi iat}\,dt. \tag{PP6}
```

The continuous square wave has Fourier support on the nonzero odd
frequencies. For positive odd `l`, the combined energy of `+l` and `-l` is

```math
{8\over\pi^2l^2},qquad
\sum_{l\ge1,\ l\text{ odd}}{8\over\pi^2l^2}=1.         \tag{PP7}
```

Consequently, for every `epsilon>0`, there is an odd `L` such that, for all
sufficiently large `p`, at least `1-epsilon` of the normalized Fourier energy
of `f_p` lies on

```math
\{a:0<|a|\le L,\ a\text{ odd}\}.                     \tag{PP8}
```

This conclusion uses only (PP6), the finite partial sum in (PP7), and
Parseval; it makes no unproved pseudorandomness assumption about Paley
graphs.

## 3. Arithmetic subsequence

Let

```math
M_L=8\prod_{\substack{r\le L\\r\text{ odd prime}}}r.  \tag{PP9}
```

Dirichlet's theorem supplies infinitely many primes

```math
p=1\pmod {M_L}.                                       \tag{PP10}
```

For every odd prime `r<=L`, quadratic reciprocity and (PP10) give

```math
\left({r\over p}\right)
=\left({p\over r}\right)=1.                          \tag{PP11}
```

Also `(-1/p)=1`. Hence every positive and negative odd frequency in (PP8) is
a quadratic residue modulo `p`.

Choose `L_j -> infinity`, and for each `j` choose a sufficiently large prime
`p_j=1 mod M_(L_j)`. Equations (PP3), (PP8), and (PP11) imply

```math
{f_{p_j}^{\mathsf T}L_{p_j}f_{p_j}\over p_j\sqrt{p_j}}
\longrightarrow1.                                    \tag{PP12}
```

The infinity contribution in (PP2) is only one by (PP5). Thus

```math
\operatorname{cap}(C_{p_j})
\ge\left({1\over2}-o(1)\right)p_j^{3/2}.              \tag{PP13}
```

On the other hand, the exact conference identity `C_p^2=pI` gives

```math
\operatorname{cap}(C_p)\le{p+1\over2}\sqrt p.        \tag{PP14}
```

Combining (PP13)--(PP14) proves the theorem

```math
\boxed{
{\operatorname{cap}(C_{p_j})\over(p_j+1)^{3/2}}
\longrightarrow\frac12.}                             \tag{PP15}
```

## 4. Consequence for a structured landing family

Prime-field Paley conferences are materially different from the square-field
family: their spectral eigenvalue `sqrt(p)` is irrational, so there is no
exact Boolean eigenvector. At the known exact orders their normalized caps
are well below `1/2`:

```text
order 6:  cap 5;    order 14: cap 21;    order 18: cap 33.
```

Equation (PP15) proves that this finite advantage cannot remain uniformly
bounded away from `1/2` over the prime-Paley design family. The mechanism is
not tensor continuation: it is concentration of a concrete Boolean spin in
one additive-character eigenspace along an arithmetic subsequence.

The usual order-filling definition—take the next admissible prime conference
and a principal restriction—contains every native order `p+1`, including
the subsequence (PP15). Therefore its landing clause still requires

```math
M_{p_j+1}=\left(\frac12+o(1)\right)(p_j+1)^{3/2}       \tag{PP16}
```

along that subsequence. If the true limiting constant is strictly below
`1/2`, this family has a linear `b`-scale gap there. If the true constant is
`1/2`, the family may land, but proving that is exactly the unknown lower
constant problem.

Thus prime-field Paley designs do not provide a simpler alternative to the
square-field landing obligation. A surviving non-tensor design family must
prevent both exact Boolean eigenvectors and the approximate Fourier
concentration mechanism (PP4)--(PP12), while still admitting a quantitative
all-order filler.
