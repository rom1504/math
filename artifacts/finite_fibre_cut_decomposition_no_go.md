# Recurrence-first finite-fibre cut decomposition and its exact no-go

## Status

This note derives a theorem-level certificate before searching for a finite
gadget.  The certificate is exact and genuinely bounded-state: a finite
signed decomposition of one microkernel would imply a uniform lifting
inequality for **every** child signing, without optimizing any parent.

The outcome is a scalable no-go for that certificate interface.  Every
full-sign microkernel has a compulsory coefficient strictly larger than the
scale-preserving multiplier `k^(3/2)`.  Thus an AlphaEvolve-like search using
this verifier can never prove uniform amplification.  This does **not** rule
out nonuniform, seed-dependent, or nonlocal lifts.

## 1. Exact lift and parent energy decomposition

Use the one-copy normalization

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
\operatorname{cap}(A)=\max_x|H_A(x)|.                 \tag{1}
```

Let `A` be an arbitrary zero-diagonal signing of order `n`.  Let `R` be a
symmetric full-sign matrix of order `k`, including diagonal signs, and let
`D` be a zero-diagonal signing of order `k`.  Define the uniform fibre lift

```math
B=A\otimes R+I_n\otimes D.                            \tag{2}
```

This is an order-`nk` signing: `R` fills every pair between distinct macro
fibres, while `D` fills each diagonal fibre.

Choose one representative from each antipodal pair of microspins,

```math
U_k=\{u\in\{\pm1\}^k:u_1=1\},\qquad q=2^{k-1},       \tag{3}
```

and put

```math
K_R(u,v)=u^{\mathsf T}Rv.                             \tag{4}
```

Every parent spin has a unique form `z_i u_(alpha_i)` with `z_i in {+-1}`
and `u_(alpha_i) in U_k`.  Direct expansion gives the exact identity

```math
H_B(z,\alpha)
=\sum_{i<j}a_{ij}z_iz_jK_R(\alpha_i,\alpha_j)
 +\sum_i H_D(u_{\alpha_i}).                           \tag{5}
```

In particular, the second term is bounded by `n cap(D)`.  The first term is
the entire inter-fibre obligation.

## 2. A finite signed cut-cone certificate

For `h in [-1,1]^q`, multilinearity gives

```math
\left|\sum_{i<j}a_{ij}z_iz_jh_{\alpha_i}h_{\alpha_j}\right|
\le \operatorname{cap}(A).                           \tag{6}
```

Indeed, for fixed `z,alpha` the left side is the absolute value of the
multilinear extension of `H_A` at the point
`(z_i h_(alpha_i))_(i=1)^n in [-1,1]^n`; its maximum occurs at a Boolean
corner.

Consequently, if the finite microkernel has a signed rank-one decomposition

```math
K_R=\sum_{h\in\{-1,0,1\}^q}c_hhh^{\mathsf T},        \tag{7}
```

then (5)--(7) prove the uniform all-child theorem

```math
\boxed{
\operatorname{cap}(A\otimes R+I\otimes D)
\le L\operatorname{cap}(A)+n\operatorname{cap}(D),
\qquad L=\sum_h|c_h|.}                               \tag{8}
```

The diagonal-fibre term is lower order for fixed `k` and growing `n`.  Thus
`L<=k^(3/2)` would imply the ledger's uniform amplification lemma along the
factor `k`, and a family with a power-saving or summable excess would be a
recurrence-first gadget mechanism.

Because a finite symmetric matrix has a decomposition of the form (7), define

```math
L(R)=\min\left\{
\sum_h|c_h|:K_R=\sum_hc_hhh^{\mathsf T},
\ h\in\{-1,0,1\}^q
\right\}.                                             \tag{9}
```

This is a finite rational linear program.  Its dual is

```math
L(R)=\max_Y\langle Y,K_R\rangle_\triangle
\quad\text{subject to}\quad
|\langle Y,hh^{\mathsf T}\rangle_\triangle|\le1
\quad\text{for every ternary }h,                     \tag{10}
```

where `triangle` means that upper-triangular coordinates are each counted
once.  This convention is used consistently by the verifier; it avoids an
off-diagonal factor of two.

The state in (7) is strictly simpler than full parent minimization: it has
size depending only on `k`, and once certified, (8) holds for arbitrary
order `n`, arbitrary signing `A`, and every parent spin.

## 3. Exact universal obstruction

The diagonal entries of (7) immediately force a much stronger obstruction
than testing isolated lifted caps.  For any two microtypes `u,v`, subtracting
the corresponding diagonal identities gives

```math
u^{\mathsf T}Ru-v^{\mathsf T}Rv
=\sum_hc_h(h_u^2-h_v^2).
```

Since `|h_u^2-h_v^2|<=1`,

```math
L(R)\ge
\max_{u\in U_k}u^{\mathsf T}Ru
-\min_{u\in U_k}u^{\mathsf T}Ru.                    \tag{11}
```

Antipodal representatives lose no values because a quadratic form is even.
Write `S=R-diag(R)`.  The diagonal part is the constant `tr(R)` on Boolean
vectors, so the right side of (11) is the range of `u^T S u`.

For a zero-diagonal symmetric matrix `S`, put

```math
W(S)=\frac12\left(
\max_x x^{\mathsf T}Sx-\min_x x^{\mathsf T}Sx
\right).                                              \tag{12}
```

This `W` uses the doubled quadratic form; it is twice the centered half-range
of the one-copy Hamiltonian in (1).

If `x,y` differ exactly on a set `T`, direct subtraction gives

```math
x^{\mathsf T}Sx-y^{\mathsf T}Sy
=4\sum_{i\in T,j\notin T}s_{ij}x_ix_j.              \tag{13}
```

Conversely every choice of signs on the two shores occurs from such a pair
`x,y`.  Therefore

```math
W(S)=2\max_{T\subset[k]}\|S_{T,T^c}\|_{\infty\to1},
\qquad
\operatorname{range}_{x}(x^{\mathsf T}Sx)
=4\max_T\|S_{T,T^c}\|_{\infty\to1}.                \tag{14}
```

Equations (11)--(14) prove

```math
L(R)\ge4\max_T\|S_{T,T^c}\|_{\infty\to1}.           \tag{15}
```

Now fix a split with `|T|=r`, `|T^c|=s`.  For uniformly random signs `y` on
`T^c`, choose each sign on `T` to match its row sum.  Since every row contains
`s` signs,

```math
\max_{x,y}x^{\mathsf T}S_{T,T^c}y
\ge\mathbb E_y\sum_{i\in T}
\left|\sum_{j\in T^c}s_{ij}y_j\right|
=r\mu_s,                                              \tag{16}
```

where

```math
\mu_s=\mathbb E|\epsilon_1+\cdots+\epsilon_s|.
```

The sharp `p=1` Khintchine inequality gives `mu_s>=sqrt(s/2)`.  Taking
`r=floor(2k/3)` and `s=k-r` yields

```math
L(R)\ge4r\sqrt{s/2}>k^{3/2}                          \tag{17}
```

for every `k>=3` except `k=4`.  The cases `k=2,4` follow from a one-column
cut in (15):

```math
L(R)\ge4(k-1)>k^{3/2}.                               \tag{18}
```

For completeness, square (17) and write `k=3a+b`, `b in {0,1,2}`.
For `b=0,1,2`, respectively, the positive differences
`8r^2s-k^3` are

```math
5a^3,
\qquad 5a^3+5a^2-9a-1,
\qquad 5a^3+10a^2+4a,                                \tag{19}
```

with `r=floor(2k/3)` and `s=k-r`.  They are positive for the applicable
`a>=1`, except that the middle expression vanishes at `a=1`, exactly the
special case `k=4`.  The reproducible verifier also evaluates the stronger
exact binomial mean in (16) through `k=256`; its smallest ratio to
`k^(3/2)` is greater than one.

The gap is uniform, not merely strict.  Direct differentiation of the three
rational functions `8r^2s/k^3` (after clearing their positive denominators)
shows that the first is constant, the second is increasing for `a>=2`, and
the third is increasing for `a>=1`.  The least squared ratios in the three
residue classes are therefore

```math
{32\over27},\qquad {384\over343},\qquad {144\over125},
```

with the middle class starting at `a=2` because `k=4` is handled by (18).
All three exceed `1.05^2`; the special `k=2,4` ratios are larger.  Hence the
safe all-order form of the obstruction is

```math
\boxed{L(R)>1.05\,k^{3/2}\quad(k\ge2).}              \tag{20}
```

Thus varying the fibre size cannot make this certificate's normalized excess
vanish, let alone make it geometrically summable.

Combining (8), (11), and (17)--(18):

> **Verified scalable no-go.** For every fibre size `k>=2` and every
> symmetric full-sign microkernel `R`, the signed cut-cone certificate (7)
> necessarily has `L(R)>k^(3/2)`.  It therefore incurs a fixed leading
> multiplier loss at every iteration and cannot establish uniform
> amplification or a summable recurrence defect.

This conclusion concerns the proof interface (7), not the actual Boolean cap
of `A tensor R`, which may exploit cancellations not visible after the
triangle inequality over atoms.

### Asymmetric/Krivine-style atoms do not evade the range obstruction

One might enlarge (7) to symmetric left/right atoms

```math
K_R=\sum_t c_t\,{p_tq_t^{\mathsf T}+q_tp_t^{\mathsf T}\over2},
\qquad p_t,q_t\in[-1,1]^q.                            \tag{21}
```

For a parent assignment, one atom contributes

```math
{1\over2}P^{\mathsf T}AQ,
```

where `P_i=z_i p_(alpha_i)` and `Q_i=z_i q_(alpha_i)`.  By multilinearity
it suffices to take Boolean `P,Q`.  Polarization for symmetric zero-diagonal
`A` gives

```math
|P^{\mathsf T}AQ|\le4\operatorname{cap}(A),          \tag{22}
```

so (21) produces the all-child coefficient `2 sum_t|c_t|`.

On the other hand, the difference of two diagonal entries of one atom in
(21) has absolute value at most two.  Therefore

```math
\sum_t|c_t|\ge{1\over2}
\operatorname{range}_{u\in U_k}u^{\mathsf T}Ru,
```

and after multiplying by the polarization cost in (22),

```math
2\sum_t|c_t|
\ge\operatorname{range}_{u\in U_k}u^{\mathsf T}Ru
>1.05\,k^{3/2}.                                      \tag{23}
```

Thus the same all-`k` theorem rules out the entire bounded symmetric
left/right box-decomposition interface as well.  This includes a direct
finite-dimensional Krivine-style decomposition when each left/right channel
is bounded separately by polarization.  A successful asymmetric construction
must retain cancellations between channels beyond the triangle inequality;
merely replacing `hh^T` by `sym(pq^T)` does not help.

## 4. Recurrence-first finite search

The verifier enumerates every symmetric full-sign `R` at `k=2,3,4`, solves
(9), rationalizes the best candidate's primal and dual solutions, and
independently checks that candidate's equalities and dual inequalities in
exact `Fraction` arithmetic.  It finds

| `k` | matrices `R` | ternary atoms | `min_R L(R)` |
|---:|---:|---:|---:|
| 2 | 8 | 4 | 4 |
| 3 | 64 | 40 | 9 |
| 4 | 1,024 | 3,280 | 16 |

The displayed `min_R` values are **exhaustive numerical LP minima with exact
certificates for the attaining candidates**.  The script does not rationalize
a dual lower certificate for every nonattaining `R`, so the word "minimum" in
the table is not a standalone exact global certificate.  No all-`k` formula
is claimed.  The all-`k` no-go is the independent analytic theorem
(15)--(18), which does not use this table.

The search therefore did not optimize isolated `cap(B)`.  Every objective it
evaluated was exactly the coefficient in the universal theorem (8), and the
dual certificates diagnose why that theorem cannot scale.

## 5. Scope and next falsifiable interface

This audit rules out a broad bounded-state proof that expands one fixed
microtype interaction kernel into scalar cut directions and applies the child
cap to each direction separately.  It includes arbitrary signed ternary
atoms, not just spectral eigenvectors, tensor-product spins, or one particular
Hadamard alphabet.

It does not rule out:

1. edge-dependent kernels `R^(ij)` coordinated nonlocally with the seed;
2. seed-dependent decompositions whose atom costs cancel before absolute
   values are taken;
3. constructions using several child signings or a growing bridge state; or
4. a theorem that controls the joint vector-valued child action without
   decomposing it into separately bounded scalar cuts.

The next exact test should therefore demand one of those features explicitly.
For edge-dependent kernels, a certificate that applies `cap(A)` separately to
each atom is no longer valid unless its coefficient across macro edges is a
scalar multiple of `A`; allowing arbitrary coefficients asks for the full
weighted cut response of the child.  Thus a proposed generalization must
exhibit a compressed joint cancellation law, not merely replace `R` by a
table of `R^(ij)`.

## Reproduction

```text
.venv/bin/python computations/finite_fibre_cut_decomposition_certificate.py \
  --output computations/results/finite_fibre_cut_decomposition_certificate.json
```

The result file contains the exact rational primal/dual supports for each
saved minimum and a canonical payload hash.  Floating-point LP is used only
to locate candidates; every reported optimum has a matching exact rational
primal and dual certificate.
