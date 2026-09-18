# Independent cross-audit of the Wave 17 mixed-product memo

## Verdict

The order-eighteen certificate and growing Sylvester/Hadamard no-go are
correct.  The exact checker passes.  The factor `2` in the ground-pair
energy, the order `18s`, and the final normalization all agree.  The
fibrewise replacement estimate `36m` is also correct, with the scope caveat
below.

## Independently checked identities

For `H=R+E`, `P=A+D`,

```math
H\otimes P-\operatorname{diag}(H\otimes P)
=R\otimes(A+D)+E\otimes A.
```

The cancelled diagonal is exactly `E tensor D`.  For fibres `x_a`, the
quadratic energy is

```math
\sum_{a\ne b}r_{ab}x_a^TPx_b
+\sum_a e_ax_a^TAx_a.
```

For product states, the tensor in the memo's ordering is `v tensor u`
(factor state first, seed state second), and

```math
(v\otimes u)^T\mathcal T(v\otimes u)
=(u^TAu+\delta)(v^TRv+\varepsilon)-\delta\varepsilon.
```

The scalar formula is correct; only the displayed tensor order `u tensor v`
in the memo should be swapped or declared permutation-equivalent.

For `B_D=H_2 tensor P-diag`, `tr(H_2)=0`.  If `p,n` are opposite `A_9`
grounds, the Boolean fibre states `(p,n)` and `(p,-n)` have energies

```math
48\pm2p^TPn.
```

Thus the exact cross-Gram lower bound `15` gives `Q(B_D)>=48+30=78`.
The checker confirms:

- positive/negative projective ground counts `10,15`;
- `512*10*15=76,800` completed cross-Gram entries;
- maximum histogram `15:7, 17:82, 19:204, 21:165, 23:48, 25:6`;
- all `2^17` projective order-eighteen states for every `D`, and in fact all
  four order-two diagonal choices `E` (`2048` completions total);
- exact minimum norm `78`, attained four times over all `(D,E)` and attained
  for the Hadamard diagonals `E=(1,-1),(-1,1)`.

Since `78>18 sqrt(17)`, an order-18 conference signing already proves
nonminimality.

For `s=4^r`, `H_s=H_4^{tensor r}`, and `y=w^{tensor r}`,
`H_sy=sqrt(s)y` and `y^TH_sy=s sqrt(s)`.  Associativity gives the unpunctured
matrix `H_s tensor K_D`.  Because `tr K_D=0`, subtracting its diagonal changes
no Boolean quadratic value.  A norm-78 order-eighteen witness therefore
gives

```math
Q(\mathcal T_m)\ge78s\sqrt s,
\qquad N=9m=18s.
```

The normalized constant is exactly

```math
78/(18 sqrt(18))=1.0213764617...>1.
```

The standard Paley-principal-submatrix bound
`q_N<=(1+o(1))N^{3/2}` applies along `N=18*4^r`, so the products are
nonminimal for all sufficiently large `r` for every seed diagonal `D`.

For fibrewise replacement, the perturbation is a direct sum of `m` copies of
`(G_3-A_3) direct-sum (G_6-A_6)`, up to harmless fibre signs.  Hence

```math
Q(perturbation)
\le m[Q(G_3-A_3)+Q(G_6-A_6)]
\le m[(6+6)+(10+14)]=36m.
```

The replicated *within-fibre* captured excess is exactly `4m`.  This is not
the full captured excess of the enlarged `3m+6m` partition; the
`R tensor (A_r+D_r)` terms can create leading enlarged-shore excess, as the
memo correctly discusses immediately afterward.

## Corrections and scope

1. Equations (2) and (3) are missing visible plus signs between their two
   sums in the Markdown source.
2. In (4), use `v tensor u` for the factor-first Kronecker ordering
   `H tensor P`, or state that `u tensor v` is after the perfect-shuffle
   permutation.
3. The converse classification in Section 1 is valid for genuinely
   Kronecker-separable rules, i.e. rules obtained by completing `R` and `A`
   on their diagonals and multiplying the two entries.  Merely prescribing
   `r_ab a_ij` when both coordinates differ leaves many nonseparable signs on
   `a=b` or `i=j`; those are not all encoded by two diagonal vectors `E,D`.
4. Describe `4m` as the replicated within-fibre excess, not the full
   enlarged-shore captured excess.
5. The exact checker does not audit the final Paley coordinate-ascent table.
   Those claims remain numerical, as labeled.  The current
   `mixed_product_r17.py` entry point reruns only outer order `198`, not the
   full list printed in the memo.

No correction is needed to the factor counts, the `78` certificate, the
Sylvester scaling, the constant `78/(18sqrt(18))`, the asymptotic comparison,
or the `36m` perturbation bound.
