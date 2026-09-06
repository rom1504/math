# Wave 32 root route: row-square noise and the local-proximity wall

## Status

The noise identities and the complete-signing construction below are
verified.  They show that trace, operator norm, and degree-two noise alone do
not put an arbitrary completion within sublinear Hamming distance of an
`O(n^2)` row center.  The construction is not an exact signing minimizer and
its distinguished cut is not asserted to be a favorable principal
completion, so this is a mechanism wall rather than a falsifier of the Wave
31 forest route.

## 1. Exact independent-noise drift

Let `B=A^2`, `R(x)=x^T Bx`, and let `epsilon_i` be independent signs with
`E epsilon_i=theta`.  Put `y=x odot epsilon` and
`v=1-theta^2`.  Since `B_ii=n-1`, direct Walsh expansion gives

```math
E R(y)=theta^2 R(x)+(1-theta^2)n(n-1).                 (N1)
```

Thus degree-two noise contracts only the *excess above the uniform mean*.
It does not create an additive target-scale improvement independent of the
starting row cost.

There is also an exact variance.  Write

```math
q_ij=2B_ij x_i x_j  (i<j),
ell_i=sum_(j ne i) q_(min(i,j),max(i,j)).
```

Expanding `epsilon_i=theta+eta_i`, the centered linear and quadratic Walsh
parts are orthogonal, as are distinct quadratic monomials.  Hence

```math
Var R(y)
=v theta^2 sum_i ell_i^2+v^2 sum_(i<j)q_ij^2.          (N2)
```

This identifies a possible minimizer-specific input—control of the local
`B=A^2` fields `ell_i`—but supplies no reverse drift beyond (N1).

## 2. Exact fixed-sphere drift

If `F` is uniform among the `d`-subsets of `[n]` and `x^F` flips `x` on
`F`, then for distinct `i,j`

```math
E[(-1)^(1_F(i)+1_F(j))]
=gamma_(n,d)
=((n-2d)^2-n)/(n(n-1))
=1-4d(n-d)/(n(n-1)).
```

Therefore

```math
E_(|F|=d) R(x^F)
=gamma_(n,d)R(x)+(1-gamma_(n,d))n(n-1).                (N3)
```

In particular some radius-`d` point is no worse than the right side.  But
for `d=o(n)`, `1-gamma=O(d/n)`: a row cost of order `n^(5/2)` remains of
that order throughout the mean.  At the forest radius
`d=O(n^(3/4-c)/log n)`, the contraction fraction is only
`O(n^(-1/4-c)/log n)`.

## 3. An explicit correct-operator-scale complete-signing wall

Let `C` be a symmetric conference matrix of order `g`, so

```math
C=C^T,
diag(C)=0,
C_ij in {+-1},
C^2=(g-1)I_g.
```

Paley's construction supplies infinitely many such orders `g=q+1` with
`q` a prime power congruent to `1 mod 4`.  Put `h=g`, `N=g^2`, and define
the `N` by `N` complete signing

```math
A=C tensor J_g+I_g tensor (J_g-I_g).                   (N4)
```

Off-diagonal blocks are constant signs, diagonal blocks are positive
cliques, and the diagonal is zero.  On block-constant vectors `A` is
`gC+(g-1)I`; on the within-block zero-sum space it is `-I`.  Consequently

```math
||A||_op=g sqrt(g-1)+g-1=Theta(N^(3/4)).                (N5)
```

For the all-plus spin, using `||C 1||^2=(g-1)g`,

```math
R(1)=g ||gC1+(g-1)1||_2^2=Theta(g^5)=Theta(N^(5/2)).   (N6)
```

The triangle inequality gives, for every spin `y` at Hamming distance at
most `d` from either projective representative of `1`,

```math
sqrt(R(y))
>=sqrt(R(1))-2||A||_op sqrt(d).                         (N7)
```

Hence for every `d=o(N)`, uniformly on that projective ball,

```math
R(y)=Theta(N^(5/2)).                                    (N8)
```

In particular the full forest radius `k_0=O(N^(3/4-c_0)/log N)` contains no
`O(N^2)` center and no target-row cut.  This is stronger than an abstract
PSD wall: it is an actual zero-diagonal complete signing with exactly the
operator scale available for minimizers.

The construction is deliberately scoped.  It is not claimed to minimize
`Q(A)`, and the high-row all-plus cut is not tied to any principal child
ground.  It proves that a Wave 32 proximity theorem must use the joint
structure of exact minimization and favorable principal completions; the
generic identities (N1)--(N3), trace, and the minimizer-scale operator norm
cannot prove it.

## 4. Sharp remaining local input

A necessary singleton input is that a favorable-completion fiber cannot be
trapped in a high-row basin such as (N8).  The full sufficient local theorem
is stronger: for every relevant batch group, choose one favorable completion
`x^a` per selector and an `R_2=O(n^2)` center `z`, and join the whole family
by a rooted tree whose **total over all edges** is `O(k_0)`.  Along that tree,
the exact class-switch drift (10.914) must remove every high-row excess using
only `O(k_0)` collective coordinate changes.  Per-selector proximity alone
does not provide this group clustering, and neither (N1) nor uniform low-row
abundance provides even the required conditional-fiber input.
