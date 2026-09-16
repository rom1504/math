# Critical-scale regular-block lift of the strict nine-type obstruction

This is a scratch proof memo.  It is not a tracked ledger entry.

## Claim

Let `w=(w_ij)` be the strict nine-type rational certificate (10.448),
with projective cut minimum/maximum gaps

```math
\delta=1/112,
```

and shores `S={0,1,2,3}`, `T={4,5,6,7,8}`.  There is a fixed
constant `kappa>0` and, for every sufficiently large `L` in an infinite
arithmetic progression, a complete zero-diagonal sign matrix `A_L` of
order `n=9L` such that

1. the block-constant positive and negative macro spins are the exact
   projective endpoints of `A_L`;
2. both endpoint children are positive-dominant;
3. (10.435), (10.437), and (10.440) fail at the root; and
4. `Q(A_L)=O(n^{3/2})`.

In particular, the dense `Theta(n^2)` blow-up from Section 10.61 can be
replaced by a genuinely critical-scale `O(n^{3/2})` obstruction.  The
average sign in block `ij` is `kappa w_ij/sqrt(L)+O(1/L)`.

## Regular signed blocks

Restrict to `L=1 mod 4`.  Fix `kappa`, to be chosen after the spectral
constant below.  For every `i<j`, let `r_ij` be a nearest odd integer to

```math
kappa w_ij sqrt(L).
```

Thus `|r_ij-kappa w_ij sqrt(L)|<=1`, and for large `L`,
`0<|r_ij|<L`.  Put

```math
d_ij=(L-|r_ij|)/2.
```

Choose an `L` by `L` zero-one matrix `H_ij` with every row and column
sum equal to `d_ij`, and define

```math
B_ij=sgn(r_ij)(J-2H_ij).
```

Then `B_ij` is a sign matrix with every row and column sum `r_ij`.  If

```math
E_ij=B_ij-(r_ij/L)J,
```

then

```math
E_ij=-2 sgn(r_ij)(H_ij-(d_ij/L)J).
```

Fix `alpha=1/2` and `m=2` in Theorems A and B below.  Theorem B of
Tikhomirov--Youssef, *The spectral gap of dense random
regular graphs*, applies because `d_ij=L/2-O(sqrt(L))`.  It implies that
for all sufficiently large `L` there is such a regular matrix with

```math
||E_ij||_op <= Gamma sqrt(L),
```

where `Gamma` is an absolute constant independent of `i,j,L` (enlarge
one theorem constant to cover the finitely many blocks).  Their model
is exactly all zero-one square matrices with fixed row and column sum;
loops are allowed, as needed here.  The arXiv source is
<https://arxiv.org/abs/1610.01765>.

Inside every clone class, choose a simple `(L-1)/2`-regular graph `G_i`
and put

```math
B_ii=J-I-2G_i.
```

This is symmetric, zero-diagonal, has sign off-diagonal entries, and
has every row sum zero.  Theorem A of the same paper supplies a choice
whose nontrivial adjacency eigenvalues are `O(sqrt(L))`.  Hence, after
enlarging `Gamma`,

```math
||B_ii||_op <= Gamma sqrt(L),
qquad B_ii 1=0.
```

The congruence `L=1 mod 4` makes `(L-1)/2` even, so these regular simple
graphs exist.  Since there are only finitely many blocks, the
high-probability conclusions of Theorems A and B can be intersected, or
the blocks can simply be selected separately.

Assemble the symmetric matrix `A_L` from the `B_ij`, with
`B_ji=B_ij^T`.  Let `P` be orthogonal projection onto vectors constant
on each clone class.  Write

```math
A_L=D_L+E_L,
```

where `(D_L)_ij=(r_ij/L)J` off the block diagonal and `(D_L)_ii=0`.
Every block of `E_L` has zero row and column sums.  Therefore

```math
E_L P=P E_L=0,
qquad
||E_L||_op <= 9 Gamma sqrt(L).
```

The last estimate follows by comparing the nine by nine matrix of
block operator norms with `Gamma sqrt(L) J_9`.

## Exact endpoint stability

Put

```math
beta_ij=r_ij/sqrt(L).
```

Then `beta_ij=kappa w_ij+O(L^{-1/2})`.  For all large `L`, every
nontrivial macro cut has `beta`-weight at least `kappa delta/2`; the
rounding error in a macro cut is at most `20/sqrt(L)`.

Let `R` be a micro cut, replace it by its complement if needed, and put

```math
d_i=|R cap V_i|,
theta_i=d_i/L,
d=sum_i d_i <= 9L/2,
s=d/L.
```

If `Z_i` are independent Bernoulli variables of parameters `theta_i`,
the cut contributed by `D_L` is

```math
C_D(R)=L^{3/2} E g_beta({i:Z_i=1}).
```

The expected number of disagreeing macro pairs is

```math
D(theta)=8s-s^2+sum_i theta_i^2
        >=8s-(8/9)s^2
        >=4s.
```

A nonconstant nine-bit vector has at most twenty disagreeing pairs.
Consequently,

```math
Pr(Z is nonconstant)>=D(theta)/20>=s/5,
```

and hence

```math
C_D(R)>=(kappa delta/10) sqrt(L) d.
```

For the spin `x` associated with `R`, write `x=Px+u`.  Since
`E_LP=0`, the exact cut identity is

```math
C_A(R)=C_D(R)-(1/4)u^T E_L u.
```

Moreover,

```math
||u||_2^2
=sum_i L(1-(1-2d_i/L)^2)
<=4d.
```

It follows that

```math
|(1/4)u^T E_L u|<=9 Gamma sqrt(L)d.
```

Thus the block-constant all-one spin is the unique projective positive
endpoint whenever

```math
kappa delta/10>9 Gamma.
```

For example it is enough to take

```math
kappa>90 Gamma/delta=10080 Gamma.
```

Let `nu` be the block-constant `S|T` endpoint.  The signing

```math
A_L^-=-diag(nu) A_L diag(nu)
```

has the same regularity and spectral estimates.  Its macro cut weights
are `beta^-_ij=-nu_i nu_j beta_ij`, and

```math
g_{beta^-}(Y)=c_beta-g_beta(T triangle Y).
```

The strict upper gap in (10.449), again stable under rounding, says
that every nontrivial `beta^-` cut has weight at least
`kappa delta/2`.  Applying the preceding argument to `A_L^-` proves

```math
C_{A_L}(T)-C_{A_L}(T triangle R)=C_{A_L^-}(R)>0
```

for every nontrivial projective `R`.  Thus `nu` is the exact negative
endpoint.

This is the point at which regularity is essential.  Independent
biased edges have `Theta(sqrt(L))` noise in every single-vertex row and
cannot survive a union bound over the `Theta(L)` cuts at distance one.
Exact row and column sums make all first-order centered terms vanish;
the remaining error is bilinear in deviations and is bounded by
`O(sqrt(L)d)` with no logarithmic loss.

## Child stability and the stopping gap

For a shore `X`, let `a_X(beta),m_X(beta),M_X(beta)` be the finite
macro total, minimum cut, and maximum cut.  Let the superscript `(L)`
denote the corresponding micro quantities.  Group-constant cuts are
represented exactly, and every fractional macro cut is a convex
combination of macro vertex cuts.  If `X` contains `t<=5` macro types,
then

```math
||E_L[X]||_op<=t Gamma sqrt(L),
qquad ||u_X||_2^2<=tL.
```

Therefore, with `eta=25 Gamma/4`,

```math
L^{3/2}(m_X(beta)-eta) <= m_X^(L) <= L^{3/2}m_X(beta),
L^{3/2}M_X(beta) <= M_X^(L) <= L^{3/2}(M_X(beta)+eta),
a_X^(L)=L^{3/2}a_X(beta).
```

The exact macro positive-dominance gaps are `9/56` on `S` and `31/84`
on `T`.  Since `beta -> kappa w`, choosing the already fixed `kappa`
large enough compared with `Gamma` makes

```math
a_X^(L)-m_X^(L)-M_X^(L)>0
```

on both children for all sufficiently large `L`.  Indeed, the
one-sided bounds give

```math
a_X^(L)-m_X^(L)-M_X^(L)
>=L^{3/2}(a_X(beta)-m_X(beta)-M_X(beta)-eta).
```

More importantly, the hard stopping gap cannot be decreased by the
centered blocks: since a group-constant macro minimizer is an available
micro cut,

```math
m_X^(L)<=L^{3/2}m_X(beta).
```

Hence

```math
a_S^(L)+a_T^(L)-m_S^(L)-m_T^(L)-c_L
>=L^{3/2}
  [a_S(beta)+a_T(beta)-m_S(beta)-m_T(beta)-c_beta].
```

The bracket tends to `kappa/48`, so it is positive for all large `L`.
This proves failure of (10.440), and positive dominance gives failure
of (10.437).

It also proves the actual stopping failure (10.435), not only its
stronger sufficient form.  Here `a_X^(L)>0` for large `L`.  In doubled
quadratic normalization,

```math
c_beta-[a_S(beta)+a_T(beta)] -> kappa (79/168)>0,
c_beta+[a_S(beta)+a_T(beta)] -> kappa (257/168)>0.
```

Thus the intended parent maximum and minimum energies have the required
opposite signs for all large `L`.  Consequently, in doubled quadratic
normalization,

```math
I_L=4(a_S^(L)+a_T^(L)),
Q_X=2a_X^(L)-4m_X^(L),
b_X=2c_L-Q_X+2a_X^(L)=2c_L+4m_X^(L).
```

Consequently the exact identity

```math
I_L-(b_S+b_T)
=4[a_S^(L)+a_T^(L)-m_S^(L)-m_T^(L)-c_L]
>0
```

holds.

## Critical quadratic norm

On the block-constant subspace, `D_L` is represented by the nine by
nine matrix `(r_ij)`, so

```math
||D_L||_op<=9 max_{i<j}|r_ij|=O(kappa sqrt(L)).
```

Together with `||E_L||_op<=9 Gamma sqrt(L)`, this yields

```math
||A_L||_op=O(sqrt(L)).
```

Since `n=9L`,

```math
Q(A_L)=max_x |x^T A_L x|
      <=n||A_L||_op
      =O(n^{3/2}).
```

Thus `kappa` is large but fixed, and the construction remains exactly
at the critical scale.

## Hypotheses imported from the literature

Only the following existence facts are imported.

1. A dense uniform `d`-regular directed graph (equivalently, a square
   zero-one matrix with every row and column sum `d`) has second
   singular value `O(sqrt(d))` with probability tending to one,
   uniformly for `L^alpha<=d<=L/2`.
2. A dense uniform simple undirected `d`-regular graph has every
   nontrivial adjacency eigenvalue `O(sqrt(d))` with probability
   tending to one in the same range.

These are exactly Theorems B and A, respectively, in Tikhomirov and
Youssef, arXiv:1610.01765 (Annals of Probability 47 (2019), 362--419).
The paper explicitly allows loops in the directed model and identifies
its adjacency matrices with all zero-one matrices having the prescribed
constant row and column sums.  Our degrees are `L/2-O(sqrt(L))`, so all
range and nonemptiness hypotheses hold for sufficiently large
`L=1 mod 4`.
