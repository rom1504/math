# Partial matchings as a high-rank cut-code state: exact audit

Date: 2026-08-15. Independent scratch report; not a tracked project file.

## 1. Verdict

The partial-matching orbit has exactly the hidden size suggested by the
growing-state heuristic:

```math
|\mathcal M_l|
={n!\over 2^l l!(n-2l)!}
=\exp\{\alpha n\log n+O(n)\},\qquad l=\alpha n.
```

Its normalized add/delete coefficient is

```math
c_l={\sqrt{(l+1){n-2l\choose2}}\over E},
\qquad E={n\choose2},
```

and the radial matching path has asymptotic spectral scale

```math
2\max_l c_l
={4\over3\sqrt3}{1+o(1)\over\sqrt n}
=0.769800\ldots\,n^{-1/2}.
```

Thus the orbit evades the `exp(O(n))` support obstruction and genuinely
reaches the required order.  But the two most natural realizations are now
closed by exact theorems.

1. A scalar amplitude-square on the matching orbit loses the entire
   `n^{-1/2}` gain under the cut-code twirl: its rooted certificate is at
   most `2/sqrt(E)=O(1/n)`.
2. A kernel whose positive operator Fourier coefficients are supported
   directly on matchings has `lambda<=1/(n-1)`, at every matrix rank.  This
   follows by marginalizing pointwise positivity to one vertex star.

The main surviving `S_n`-equivariant partial-matching architecture has
Fourier support on **differences of matchings** (alternating paths and even
cycles) and uses a non-rank-one positive Gram block.  Perfect-matching fibers themselves are
multiplicity-free, so their invariant algebra is commutative.  Noncommuting
blocks first arise from multiplicities across partial-matching levels.  This
is a real algebraic opening, but the arbitrary root signing destroys the
symmetry that block-diagonalizes the unrooted operator.  No cited
association-scheme theorem supplies the required uniform signed root bound.

## 2. Cut-code fibers are fixed-boundary perfect matchings

Use the augmented cut code

```math
C_n^+=\{(\sigma x_ix_j)_{i<j}:\sigma,x_i\in\{\pm1\}\}.
```

For an edge set `S`, the restriction of its character to `C_n^+` is
determined by

```math
\pi(S)=(|S|\bmod2,\partial S),
```

where `partial S` is the set of odd-degree vertices.  If `M` is an
`l`-matching, then

```math
\pi(M)=(l\bmod2,U(M)),
```

where `U(M)` is its `2l`-vertex covered set.  Therefore two matchings have
the same cut-code syndrome exactly when they are perfect matchings of the
same covered set `U`.  The fiber size is

```math
p_l=(2l-1)!!,
\qquad |\mathcal M_l|={n\choose2l}p_l,
\qquad {p_{l+1}\over p_l}=2l+1.                    \tag{2.1}
```

Two matchings in one fiber differ by a disjoint union of alternating even
cycles.  This is exactly the perfect-matching association scheme inside
that fiber.

## 3. Exact no-go for a scalar matching amplitude square

Let `v_M=w_l>=0` on all `l`-matchings and zero elsewhere.  Put

```math
g(a)=\sum_Mv_M\chi_M(a),\qquad K(a)=g(a)^2.
```

The Fourier coefficients of `K` are the nonnegative autocorrelation
`q=v*v`.  If normalized cube adjacency `A` satisfies

```math
Av\ge\lambda v,
```

then `Aq=(Av)*v>=lambda q`; hence both `K` and
`F=(tau-lambda)K` are positive definite, while `K>=0` pointwise.  This is
an admissible scalar rooted kernel.

Set

```math
\alpha_l=\sqrt{{n\choose2l}p_l}\,w_l,
\qquad
c_l={\sqrt{(l+1){n-2l\choose2}}\over E}.
```

The super-eigenvector condition becomes

```math
c_{l-1}\alpha_{l-1}+c_l\alpha_{l+1}
\ge\lambda\alpha_l.                               \tag{3.1}
```

For an edge signing `a`, cut-code orthogonality and the fiber description
give the exact rooted mass

```math
{T_a\over|C_n^+|}
=\sum_lw_l^2\sum_{|U|=2l}
 \operatorname{haf}(a[U])^2,                       \tag{3.2}
```

where

```math
\operatorname{haf}(a[U])
=\sum_{M\text{ perfect on }U}\prod_{e\in M}a_e.
```

In particular

```math
{T_1\over|C_n^+|}
=D:=\sum_lp_l\alpha_l^2.                           \tag{3.3}
```

The cut-code fiber-sum vector has normalized level amplitude

```math
\beta_l=\sqrt{p_l}\,\alpha_l.
```

The quotient syndrome graph has coefficient

```math
\bar c_l
={\sqrt{{n-2l\choose2}{2l+2\choose2}}\over E}
=\sqrt{2l+1}\,c_l.                                 \tag{3.4}
```

Consequently, if `rho_bar` is the Rayleigh quotient of this fiber-sum
vector,

```math
\rho_{\rm bar}D
=2\sum_lp_{l+1}c_l\alpha_l\alpha_{l+1},
\qquad {J\over T_1}=\rho_{\rm bar}-\lambda.         \tag{3.5}
```

Multiplying (3.1) by `p_l alpha_l` and summing yields

```math
\lambda D
\le\sum_l(p_l+p_{l+1})c_l\alpha_l\alpha_{l+1}.
```

Subtracting (3.5) gives

```math
\lambda-{J\over T_1}
=2\lambda-\rho_{\rm bar}
\le {2\sum_lp_lc_l\alpha_l\alpha_{l+1}\over D}.
```

Using `beta_l=sqrt(p_l)alpha_l` and (2.1), the right side is the
Rayleigh quotient of a tridiagonal path with off-diagonal coefficients

```math
d_l={c_l\over\sqrt{2l+1}}.
```

Therefore

```math
\lambda-{J\over T_1}
\le 2\max_ld_l
\le {2\over\sqrt E}.                               \tag{3.6}
```

Finally every Fourier coefficient of `K` is nonnegative, so
`T_a<=T_1`; also `J>=0`.  Hence, whenever `T_a>0`, the actual rooted
lower bound obeys

```math
\boxed{
\lambda-{J\over T_a}\le {2\over\sqrt E}=O(n^{-1})
\quad\text{for every }a\text{ with }T_a>0.}        \tag{3.7}
```

If `T_a=0`, this kernel gives no rooted certificate at that root, so the
same no-go conclusion holds.  The proof of (3.7) assumes that the scalar
amplitude is `S_n`-invariant (`v_M` depends only on `|M|`) and nonnegative;
it does not cover an arbitrary nonradial matching amplitude.

This is not merely a bad finite constant.  The factorial fiber collision
`p_{l+1}/p_l=2l+1` cancels the full leading matching spectral scale.

## 4. Exact rank-independent no-go for direct matching Fourier support

There is a distinct possible construction: put the Fourier coefficients of
an operator kernel directly on matchings,

```math
K(a)=\sum_{M\text{ matching}}\chi_M(a)A_M,
\qquad A_M\succeq0,                                 \tag{4.1}
```

and require `K(a)>=0` pointwise and
`F=(tau-lambda)K` operator-positive-type.

This too cannot reach `n^{-1/2}`.  Marginalize (4.1) over every edge except
the `n-1` edges incident with a fixed vertex `i`.  Since a matching
contained in a star has size at most one, pointwise positivity gives

```math
A_\varnothing+\sum_{j\ne i}x_jA_{ij}\succeq0
\quad(x_j\in\{\pm1\}).
```

Taking every `x_j=-1` gives

```math
\sum_{j\ne i}A_{ij}\preceq A_\varnothing.
```

Summing over `i`,

```math
2\sum_eA_e\preceq nA_\varnothing.                 \tag{4.2}
```

At the empty Fourier index, positivity of `F` requires

```math
{1\over E}\sum_eA_e\succeq\lambda A_\varnothing.  \tag{4.3}
```

Combining (4.2)--(4.3), on the support of `A_empty`, proves

```math
\boxed{\lambda\le {n\over2E}={1\over n-1}.}        \tag{4.4}
```

This theorem allows arbitrary matrix rank and noncommuting `A_M`.  Thus
neither rank nor the association scheme rescues Fourier support that is
literally confined to partial matchings.

## 5. What representation theory does and does not compress

For a fixed `2l`-set `U`, the perfect matching permutation module is
multiplicity-free:

```math
\mathbb C[\operatorname{PM}(U)]
\cong\bigoplus_{\lambda\vdash l}S^{2\lambda}.
```

Equivalently `(S_{2l},S_2\wr S_l)` is a Gelfand pair, and its
Bose--Mesner algebra is commutative.  Its orbitals are indexed by the even
cycle type of the union of two perfect matchings.  Thus an invariant kernel
inside one boundary fiber decomposes into scalar spherical channels; there
are no noncommuting blocks there.

The module of all `l`-matchings on `[n]` is

```math
\mathbb C[\mathcal M_l]
=\operatorname{Ind}^{S_n}_{S_{n-2l}\times(S_2\wr S_l)}\mathbf1,
```

with Frobenius characteristic

```math
h_{n-2l}h_l[h_2]
=h_{n-2l}\sum_{\lambda\vdash l}s_{2\lambda}.        \tag{5.1}
```

Pieri's rule shows that this module is not generally multiplicity-free.
Already for `n=6,l=2`, `S^(4,2)` occurs twice: it is obtained from both
`2(2)=(4)` and `2(1,1)=(2,2)` by adding a horizontal 2-strip.  Hence the
direct sum over partial-matching levels does have genuine matrix-valued
multiplicity spaces, and add/delete incidence acts nontrivially on them.

This is the surviving algebraic opening.  It is, however, only an
**unrooted** compression.  A general edge signing multiplies each matching
basis vector by

```math
a_M=\prod_{e\in M}a_e.
```

That diagonal twist does not commute with `S_n` and mixes the isotypic
blocks in (5.1).  The perfect-matching association scheme diagonalizes
cycle-type convolution, not the arbitrary signed vectors
`(a_M)_{M\in PM(U)}`.  The trivial fiber component is the hafnian in (3.2);
nontrivial components are signed matching immanants/harmonics.

There is also a useful positivity warning.  A projector onto a nontrivial
perfect-matching eigenspace annihilates the all-ones matching vector at
`a=1`.  If its rooted energy came by itself from an edge-cube positive
definite kernel, Fourier positivity would imply `T_a<=T_1=0`, so it would
vanish identically.  Nontrivial matching harmonics can therefore enter only
jointly with enough trivial/baseline mass to restore edge-cube positivity;
they are not free positive channels.

## 6. Exact surviving target and falsification test

A partial-matching route now has to construct a positive Gram kernel whose
edge-cube Fourier support includes symmetric differences of matchings.
Those differences are unions of alternating paths and even cycles.  It must
simultaneously satisfy

```math
K(a)\succeq0\quad\text{pointwise},
\qquad K\succeq_{\rm PD}0,
\qquad (\tau-\lambda)K\succeq_{\rm PD}0,
\qquad \lambda=\Theta(n^{-1/2}),                    \tag{6.1}
```

and prove a root inequality after the arbitrary diagonal matching twist.
An `S_n`-equivariant decomposition may compress the unrooted Gram budget
to multiplicity matrices, but it counts as a real escape only if the rooted
quantity also closes in those matrices.  Merely carrying every signed
principal hafnian/immanant or every labeled alternating-cycle coefficient
is the full growing root profile in another basis.

Reject a proposed matching construction if any of the following occurs.

1. Its Fourier support is directly on matchings: (4.4) applies.
2. It is a single nonnegative matching amplitude square: (3.7) applies.
3. It diagonalizes each perfect-matching fiber and pays the scalar
   spherical channels separately.
4. It uses association-scheme eigenvalues only for the unrooted numerator
   while leaving the arbitrary signed root vector uncontrolled.
5. Its compressed state contains all labeled principal matching
   harmonics and therefore reconstructs an exponential root profile.

The smallest genuinely open statement is a cross-level multiplicity-space
Gram identity for which the complete alternating-cycle remainder is
positive only jointly, together with a uniform generalized-eigenvalue or
determinant bound after every edge-sign diagonal twist.

Outside the `S_n`-equivariant route, (3.7) also leaves open a nonradial
rank-one amplitude on the matching set.  Such a choice forfeits the proposed
association-scheme compression and would need a separate uniform
support-versus-fiber-collision theorem.  The arbitrary-rank theorem (4.4),
by contrast, has no radiality, equivariance, or commutativity assumption.

## 7. Literature boundary

Primary sources checked:

- M. K. Srinivasan, *The perfect matching association scheme*, Algebraic
  Combinatorics 3 (2020), 559--591,
  https://doi.org/10.5802/alco.104 .  It gives the multiplicity-free
  decomposition and the orbital/eigenvalue algorithms for fixed-boundary
  perfect matchings.
- M. K. Srinivasan, *Incidence matrices for matchings*, Special Matrices 6
  (2018), 297--300, https://doi.org/10.1515/spma-2018-0024 .  It proves the
  matching up maps have full row rank in the semi-Peck range, but supplies
  no signed-root inequality.
- D. Gijswijt and S. Polak, *Semidefinite lower bounds for covering codes*,
  arXiv:2504.01932v2 (2026), https://arxiv.org/abs/2504.01932 .  Its
  polynomial block diagonalization uses the full Hamming automorphism group
  and optimizes over an arbitrary covering code.  It does not determine the
  covering radius of the fixed cut code, whose root has only `S_n`
  symmetry.

These sources verify that the association-scheme machinery is real, but
none supplies the missing uniform control of the signed cut-code root.
