# Trivial-channel obstruction for nonabelian fusion schemes

Date: 2026-07-31. This is an agent-authored theorem-level audit of a
genuinely scale-dependent nonabelian/association-scheme proposal. It does not
assume a fixed tensor generator.

## 1. Exact signing mapping

Let `X` have size `n`, and let

```math
I=A_0,A_1,\ldots,A_d
```

be the adjacency matrices of a homogeneous association scheme. Thus the
off-diagonal relations partition `X x X`, and every `A_j` has a constant row
sum `k_j`. Choose signs `epsilon_j` and define

```math
S=\sum_{j=1}^d\epsilon_jA_j.                           \tag{NF1}
```

When the scheme is symmetric, (NF1) is a symmetric zero-diagonal signing of
the complete graph. Its definition may change completely with `n`; neither a
fixed generator nor tensor closure is assumed.

The nonabelian group case is exact. For a finite group `G`, fuse inverse
conjugacy classes `C_j` and put

```math
S_{u,v}=\epsilon_j
\quad\text{when }u^{-1}v\in C_j,quad u\ne v.          \tag{NF2}
```

Then (NF2) is (NF1) for the conjugacy-class scheme, with `k_j=|C_j|`.
Arbitrary scale-dependent choices of groups, class fusions, and signs are
allowed.

## 2. Exact regular-channel lower bound

Every matrix in a homogeneous Bose--Mesner algebra preserves the trivial
idempotent. In elementary terms,

```math
S\mathbf1=r\mathbf1,qquad
r=\sum_{j=1}^d\epsilon_jk_j\in\mathbb Z.              \tag{NF3}
```

The all-ones vector is Boolean, so the project's two-sided absolute
Hamiltonian obeys

```math
\boxed{
\operatorname{cap}(S)\ge {n\over2}|r|.}               \tag{NF4}
```

This is an exact signing inequality, not a one-sided cut or frustration
quantity.

## 3. Exact and approximate orthogonal designs saturate `1/2`

Suppose first that the fusion produces a conference signing,

```math
S^2=(n-1)I.                                           \tag{NF5}
```

Apply (NF5) to `1` and use (NF3). Then

```math
r^2=n-1.                                              \tag{NF6}
```

Consequently `n-1` must be a square and

```math
\boxed{
\operatorname{cap}(S)={n\over2}\sqrt{n-1}.}           \tag{NF7}
```

Equality follows because the ordinary spectral bound supplies the matching
upper estimate. Thus every conference signing lying in a homogeneous
association scheme has a Boolean extremizer in the trivial representation
and exactly saturates the conference constant.

The conclusion is stable. Suppose a scale-dependent family satisfies

```math
\left\|S_n^2-(n-1)I\right\|_{op}\le\eta_n n,
\qquad \eta_n\longrightarrow0.                        \tag{NF8}
```

Since `1` is an eigenvector, (NF8) gives

```math
|r_n^2-(n-1)|\le\eta_n n.                             \tag{NF9}
```

Equations (NF4) and (NF9) prove

```math
\boxed{
\operatorname{cap}(S_n)
\ge\left(\frac12-o(1)\right)n^{3/2}.}                 \tag{NF10}
```

Hence any homogeneous fusion family whose tractable certificate is
asymptotic orthogonality cannot have a normalized cap bounded below `1/2`.
This includes scale-dependent nonabelian groups with irreducible dimensions
tending to infinity: the high-dimensional representations do not remove the
trivial one-dimensional representation.

## 4. Difference-set and group-developed Hadamard formulation

The same obstruction can be written without association-scheme language.
Let `phi:G->{+1,-1}` be inverse-symmetric and set

```math
Q_{u,v}=\phi(u^{-1}v).                                 \tag{NF11}
```

If a nonabelian difference-set construction makes `Q` Hadamard,

```math
QQ^{\mathsf T}=nI,                                    \tag{NF12}
```

then its common row sum

```math
R=\sum_{g\in G}\phi(g)                                \tag{NF13}
```

satisfies `R^2=n`. Remove the diagonal to form the complete-graph signing
`A=Q-phi(e)I`. The all-ones spin gives

```math
\operatorname{cap}(A)
\ge {n\over2}|R-\phi(e)|
\ge {n\over2}(\sqrt n-1).                             \tag{NF14}
```

Thus nonabelian difference sets, cocyclic Hadamards that remain
group-developed, and their scale-dependent fusions inherit the same
`1/2-o(1)` obstruction.

## 5. Composition and landing consequence

Wreath products, conjugacy-class fusions, and scale-dependent class signings
can give a compact eigenmatrix state and closed multiplication rules. But if
their proposed composition theorem preserves (NF8), its exact landing bound
is already constrained by (NF10). Combined with the project's all-order
Paley upper construction,

```math
M_n\le\left(\frac12+o(1)\right)n^{3/2},               \tag{NF15}
```

the scheme can at best reproduce the existing `1/2` family. Its landing
condition is again the unknown assertion that the true optimum has leading
constant `1/2`. If the true limiting constant is strictly smaller, (NF10)
implies a linear gap after taking the `2/3` power.

This is a scalable obstruction rather than a claim that nonabelian groups
are irrelevant. It identifies the exact necessary escape:

> A useful nonabelian fusion must leave the homogeneous Bose--Mesner algebra
> or abandon near-orthogonality as its cap certificate.

The first option requires a multi-fiber coherent configuration with
nonconstant row sums and a Boolean theorem controlling its quotient and
nontrivial modules simultaneously. The second requires a cap-specific
inequality proving Boolean slack despite a larger operator norm. Merely
choosing groups with growing nontrivial representation dimensions does not
remove the trivial Boolean channel and therefore does not produce a new
landing mechanism.
