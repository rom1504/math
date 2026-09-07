# Wave 31 Route 1: exchange identities and low-row forest closure

## Scope and conclusion

This memo audits whether exact principal-ground exchanges control the two
batch statistics in (10.893)--(10.898).  The useful positive result is a
new deterministic closure theorem: a family joined by a short projective
Hamming forest to an auxiliary low-row center automatically has both few
coordinate signatures and the required `C_sig` bound.  The center need not
be favorable for any selector.

The principal-ground exchange identities themselves are only first-order
signed boundary identities.  They do not control the quadratic boundary
susceptibility that appears in `C_sig`.  A five-vertex weighted example
below makes this separation exact.  It is deliberately not a complete
signing or an asymptotic counterexample; it shows that a proof must add a
genuinely minimizer-specific quadratic input rather than re-sum the same
cycle identities.

## 1. Exact cross-regret cocycle

For an oriented full cut `d`, write

```math
c_S(d)=\langle A[S],d[S]\rangle.
```

Let `d^a` be a favorable completion for selector `S_a`, with

```math
\delta_a=Q(A[S_a])-c_{S_a}(d^a)\ge0.
```

For any two batch indices define the cross-regret

```math
\rho_b(a)=Q(A[S_b])-c_{S_b}(d^a)\ge0.
```

For a directed cycle `a_1,...,a_t`, with `a_(t+1)=a_1`, direct
telescoping gives

```math
\boxed{
\sum_{\ell=1}^t\rho_{a_{\ell+1}}(a_\ell)
=\sum_{\ell=1}^t\delta_{a_\ell}
 +\sum_{\ell=1}^t
  \{c_{S_{a_\ell}}(d^{a_\ell})
    -c_{S_{a_{\ell+1}}}(d^{a_\ell})\}.
}
```

For exact principal grounds all `delta_a` vanish.  If
`S=U union {i}` and `T=U union {j}`, its two-cycle specialization is

```math
\boxed{
\rho_T(S)+\rho_S(T)
=2\sum_{u\in U}\left[
A_{iu}(d^S_{iu}-d^T_{iu})
+A_{ju}(d^T_{ju}-d^S_{ju})
\right].
}
```

Thus cycles see signed boundary sums on the exchanged stars.  Principal
optimality is used only to assert nonnegativity of the regrets.  Global
minimality of the complete signing does not add a term to this identity.

## 2. Exact class-switch susceptibility

Fix a base spin `x`, gauge `w_ij=A_ij x_i x_j`, and put

```math
r_i=\sum_{j\ne i}w_{ij}=x_i(Ax)_i.
```

For a shore `F`, define the full signed boundary field at each vertex by

```math
a_i(F)=\sum_{j:\,1_F(i)\ne1_F(j)}w_{ij}.
```

If `x^F` is obtained by flipping `x` on `F`, then its gauged row field is
exactly `r_i-2a_i(F)`, and hence

```math
\boxed{
R_2(x^F)
=\sum_i(r_i-2a_i(F))^2
=R_2(x)+4\sum_i\{a_i(F)^2-r_i a_i(F)\}.
}
```

When `F` ranges over unions of coordinate-signature classes, the maximum
of the left side is exactly `C_sig`.  By contrast, comparison with an exact
child ground only says that the relevant signed shore total is
nonnegative.  That shore total is a linear projection of `a(F)` (for a
full selector it is `(1/2) sum_i a_i(F)`); the displayed row cost needs the
whole squared vector and its correlation with `r`.  The missing term is
therefore the quadratic boundary susceptibility

```math
\sum_i a_i(F)^2-\sum_i r_i a_i(F),
```

uniformly over signature-class unions.  It is not another cycle sign.

## 3. Low-row rooted-forest closure theorem

Use projective Hamming distance

```math
d_{\rm pr}(x,y)=\min\{d_H(x,y),d_H(x,-y)\}.
```

Consider one family `F={x^1,...,x^q}` of favorable completion spins and
an auxiliary spin `x^0`, not required to be favorable.  Join the vertices
`0,1,...,q` by a rooted tree and put

```math
D=\sum_{uv\text{ in the tree}}d_{\rm pr}(x^u,x^v).
```

Orient representatives successively from the root so every tree edge
realizes its projective distance.  Let `U` be the union of the edge
disagreement sets.  Then

```math
|U|\le D,
```

and every witness agrees with the root off `U`.  Consequently the
coordinate-signature partition of the enlarged family
`{x^0} union F` has

```math
\boxed{J\le1+|U|\le1+D.}
```

Every hybrid spin in its coarsest signature coset agrees, up to one global
sign, with `x^0` off `U`.  It can therefore be written

```math
h=\epsilon x^0+v,
\qquad \operatorname{supp}v\subseteq U,
\qquad \|v\|_2\le2\sqrt D.
```

The triangle inequality and `(a+b)^2<=2a^2+2b^2` give the exact uniform
bound

```math
\boxed{
C_{\rm sig}(F\cup\{x^0\})
\le 2R_2(x^0)+8\|A\|_{\rm op}^2D.
}
```

Adding the auxiliary center to the signature family does not hurt: the
containing coset includes the center and every favorable witness, while
only witnesses are used to hit selectors.

For every complete signing, uniform averaging gives

```math
\mathbb E_{x\sim U_{\{\pm1\}^n}}R_2(x)
=\operatorname{tr}A^2=n(n-1),
```

so at least one admissible auxiliary center has `R_2(x^0)<=n(n-1)`.
Moreover Markov gives

```math
\Pr\{R_2(x)\le2n(n-1)\}\ge\frac12,
```

so low-row centers are abundant up to a harmless factor two.  The open
issue is finding one of these centers close to the chosen favorable
witnesses.  Abundance by itself supplies no forest and no collision.

For an exact minimizer, `||A||_op^2<=2q_n=O(n^(3/2))`.  At the structural
scale

```math
L_0=n^{3/4-c_0},
\qquad k_0=Theta(L_0/log n),
```

the conditions

```math
R_2(x^0)\le 2n(n-1),
\qquad D\le C k_0
```

therefore imply

```math
J=O(L_0/log n),
\qquad C_{\rm sig}=O(n^{9/4-c_0}).
```

Wave 30's constrained refinement theorem then constructs one eligible
row-good common coset for the favorable witnesses in the component.

## 4. Exact adversarial-law sufficient event, including many groups

Let `w` be an arbitrary law on the selector slice and let

```math
r_0=ceil(n log(2k_0)/L_0).
```

A sufficient structural event for an iid `w`-batch is:

1. the `r_0` indices can be partitioned into `T<=n^eta` nonempty groups;
2. every selector in every group has a completion with
   `delta_S<=B_(n,m)+O(n^(3/2-c_0))`;
3. for each group, its completion spins and an auxiliary center of row
   cost at most `2n(n-1)` admit a rooted tree of total projective Hamming
   length `D_t<=C k_0`.

If this event has probability at least `exp{-O(r_0 L_0)}` uniformly for
every `w`, the forest theorem gives `T` common row-good cosets covering the
realized batch.  The partial/group collision reduction then costs
`exp{O(TL_0)}`.  For `eta<c_0`, this is

```math
\exp\{O(n^{3/4-(c_0-eta)})\},
```

while the row and tolerance bounds at `c_0` are stronger than those needed
at `c'=c_0-eta`.  Hence this event is sufficient for convergence with the
degraded positive saving `c'`.  Every probability and constant here is
uniform against the original arbitrary selector law; no uniform-slice
replacement is made.

This criterion is stronger than bare common-coset incidence but weaker
than requiring one coherent completion family for the whole batch.  Its
remaining minimizer-specific input is now concrete: a sub-`n^c0` number of
low-row-centered Hamming-tree clusters for favorable principal
completions.

## 5. Why exchange signs alone do not imply the forest input

Take the symmetric zero-diagonal real weighted matrix on five vertices
with nonzero unordered weights

```text
A01=10, A02=10, A12=20, A34=20,
A03=A04=1, A13=A14=-1.
```

Let `x=(1,1,1,1,1)` and flip `F={0,1,2}`.  Both spins are exact positive
grounds: the two dense positive pieces force constancy inside `F` and
`{3,4}`, and the total cross weight is zero, so the relative piece sign is
free.  Their energies are both `120`; every cross-regret and shore-energy
exchange between them is zero.  Nevertheless

```math
a(F)=(2,-2,0,0,0),
\qquad R_2(x)=2968,
\qquad R_2(x^F)=3048.
```

Thus a zero signed exchange can carry a strictly positive quadratic row
change.  Orthogonal sums of high-gap components can also realize all eight
coordinate patterns of three ground witnesses while all energy exchanges
between component sign choices vanish.  Therefore the pair/cycle algebra
does not forbid a shattered triple; `VCdim<=2` would be a new higher-order
selection theorem.  At most, one shattered triple forces a constant amount
of Hamming-tree length, and many triples can reuse the same eight
coordinates, so there is no useful aggregate charge.

This example is a scoped algebraic separation only: it has real/zero
weights and is not a complete-signing minimizer.  It does not rule out a
new theorem that uses the discrete global minimality of `A` beyond the
known edge-flip witnesses.  The current edge-flip certificate merely
asserts existence of a near-ground witness for each flipped edge; it gives
no inequality on the boundary vectors `a(F)` of the selected principal
completions.  That is the precise missing signed/quadratic term.
