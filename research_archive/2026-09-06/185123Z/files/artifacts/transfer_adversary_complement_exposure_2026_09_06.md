# Complement exposure can reveal the entire compressed weave

Date: 2026-09-06. Exact finite obstruction to one proposed conditional-net
strategy. No general exclusion of conditional exposures, adaptive clusters,
or convergence is asserted.

## 1. Construction and the relevant randomness

Let `H_i`, `1<=i<=m`, be arbitrary Hadamard matrices of order `m`, and
let `S` be any symmetric full sign matrix of order `m`. Form the weave

```math
W_{(i,a),(j,b)}=S_{ij}H_i(a,j)H_j(b,i),\qquad W^2=m^2I.
```

In each fibre retain `k` rows and delete `l=m-k`, where `0<k,l<m`.
After this partition write `W=[[K,L],[L^T,D]]`; its full order is
`M=m^2`, and its retained order is `N=mk`. The conclusions below hold
pointwise, so apply to every law of the seed signs, row selectors,
recursive row gates, and Hadamard bases. No conditional independence
assumption is made.

The middle spectral space of `K` is `E=range L`, since
`K^2=m^2I-LL^T`. Thus revealing `L` is a natural proposal for choosing
the small-dimensional candidate net before applying an energy tail bound.
The following identities show exactly how much randomness that proposal
can consume.

## 2. Revealing both cross-block and deleted block leaves no randomness

Choose one deleted anchor `alpha_i` in each fibre. Then every retained
entry is recovered by the three-sign identity

```math
K_{(i,a),(j,b)}=
 L_{(i,a),(j,\alpha_j)}
 L_{(j,b),(i,\alpha_i)}
 D_{(i,\alpha_i),(j,\alpha_j)}.                 (1)
```

Substituting the weave formula cancels each deleted Hadamard entry
twice and leaves the desired retained entry. The identity also holds
when `i=j`, including the diagonal. Consequently conditioning on `(L,D)`
determines `K` exactly. All conditional cap events have probability zero
or one; the original independent edge signs and fresh row gates cannot
be reused as if they remained unconditioned.

In fact only the `m`-by-`m` deleted-anchor matrix
`epsilon_ij=D_((i,alpha_i),(j,alpha_j))` is needed beyond `L`.

## 3. What L alone reveals

Choose also one retained anchor `b_j` in each fibre. The observed `L`
determines the following sign matrices, of sizes `k`-by-`m` and
`l`-by-`m`, respectively:

```math
R_i(a,j)=L_{(i,a),(j,\alpha_j)},\qquad
C_i(\alpha,j)=
 L_{(j,b_j),(i,\alpha)}L_{(j,b_j),(i,\alpha_i)}.            (2)
```

Writing `z_ij=H_i(alpha_i,j)`, direct substitution gives

```math
R_i(:,j)=S_{ij}z_{ji}H_i(T_i,j),\qquad
C_i(:,j)=z_{ij}H_i(T_i^c,j),\qquad
\epsilon_{ij}=S_{ij}z_{ij}z_{ji}.
```

Each anchor row of `C_i` is all positive. Column orthogonality is exactly

```math
R_i^TR_i+\operatorname{diag}(\epsilon_i)
 C_i^TC_i\operatorname{diag}(\epsilon_i)=mI_m.             (3)
```

For two completions having the same `L`, put
`delta_ij=epsilon'_ij epsilon_ij`. Their ratio signs therefore obey

```math
\delta_{ia}\delta_{ib}=1
\quad\hbox{whenever }(C_i^TC_i)_{ab}\ne0.                 (4)
```

Let `c_i` be the number of connected components of the graph on the
`m` columns with an edge wherever the off-diagonal deleted Gram entry
is nonzero. The row `delta_i` is constant on each component. Distinct
components span mutually orthogonal nonzero subspaces of `R^l`, so
`c_i<=l`. It follows that the total number of possible retained matrices
given `L` is at most

```math
2^{\sum_i c_i}\le2^{ml}=2^{M-N}.                          (5)
```

Symmetry and any prescribed seed diagonal can only reduce this count.
In particular, under any ensemble law, the conditional Shannon entropy
of `K` given `L` is at most `(M-N) log 2`. This is a support bound, not
a claim that the surviving signs are independent or unbiased.

## 4. Odd deletion gives an exact cap-measurability obstruction

If `l` is odd, every inner product of two sign columns of `C_i` is an
odd integer, hence nonzero. All the graphs in Section3 are complete.
Thus `delta_ij=u_i` for every `j`, and symmetry gives `u_i=u_j` for
every pair of fibres. There are at most two retained completions: `K`
and `-K`. The same conclusion holds whenever all the deleted Gram graphs
are connected, even if `l` is even. Therefore

```math
Q(\operatorname{hollow}K)
\quad\hbox{is a measurable function of }L.                (6)
```

An explicit recovery in the odd case avoids searching any completions.
Use column 1 as reference, put `q_i1=1`, and for `j!=1` set

```math
q_{ij}=-\frac{(R_i^TR_i)_{1j}}{(C_i^TC_i)_{1j}}.
```

Equation (3) proves these ratios are exactly signs and equal
`epsilon_i1 epsilon_ij`. Setting
`epsilon_recovered_ij=q_1i q_ij` recovers the deleted-anchor matrix up
to its single diagonal sign `epsilon_11`. Inserting it in
`K_ij=epsilon_ij R_i(:,j)R_j(:,i)^T` recovers `K` up to that same
global sign. If the seed diagonal is fixed, even that ambiguity is fixed.

This parity case is compatible with any limiting retention `p in (0,1)`:
choose the nearest odd `l` to `(1-p)m`, changing `k/m` by only `O(1/m)`.
This observation does not replace the specified ensemble by fiat. It
shows that a proposed exposure argument relying only on the limiting
retention cannot assume fresh cap randomness survives `L`.

## 5. Exact cap and partition consequences, and remaining options

For odd `l`, conditional on `L`, every threshold `u` satisfies

```math
\Pr\{Q(\operatorname{hollow}K)>u\mid L\}
=\boldsymbol1_{\{Q(\operatorname{hollow}K)>u\}}.
```

Moreover the two-sided hollow partition

```math
Z_+(b)+Z_-(b)
=\sum_x\left(e^{bH_K(x)}+e^{-bH_K(x)}\right)
```

where `H_K(x)=x^T hollow(K)x/2`,
is invariant under the one remaining global sign and is therefore also
`L`-measurable. Consequently revealing `L`, choosing a small spectral
net from `E=range L`, and then applying the unconditional row/edge
annealed bound to those candidates is circular in this case. The
conditioning has already determined the very cap and partition being
estimated. Revealing `(L,D)` has this defect for every positive `l`.

These statements do NOT rule out a coarser statistic of the complement,
a partial exposure leaving useful randomness, or a joint unconditional
counting estimate for low-dispersion, high-energy spins. They identify
the exact failure of the full-complement and full-cross-block proposals.
For general even `l`, (5) does not itself give a conditional tail bound;
the law of the surviving completions still needs to be proved.

Reproduction:
`computations/transfer_adversary_complement_exposure_2026_09_06.py`
checks (1)--(5) and explicit odd-deletion recovery on finite actual weaves.
