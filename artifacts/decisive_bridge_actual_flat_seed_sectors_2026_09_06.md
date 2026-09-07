# Actual all-order flat-spectrum sectors retain the Boolean seed cap

Date: 2026-09-06. Status: exact restricted-sector computation, not an upper
bound for all spins. It contrasts with the Gaussian diagnostic, which sees
seed singular values in the fully correlated-column family.

Take k=4 and the regular Hadamard H_4=J_4-2I_4. Its tensor powers have order
m=4^r and row and column sums sqrt(m). Write m=k d, with d=4^(r-1).
Use this same regular basis in every group of the literal vector seed weave,
and retain all rows (ell=m). Fix c in {+-1}^k. In EVERY group alpha choose
physical spins `x_(i,alpha)=c_i 1_m`. Their spectra are exactly
`h_(i,alpha)=c_i sqrt(m) 1_m`. Thus every normalized outgoing block is the
same matrix `X=c 1^T`, with ||X||²=k². This is an actual joint spectrum,
not an arbitrary prescribed vector source. Independent block permutations
do not change it.

For a symmetric full sign seed A,

`<X,R_A X>=sum_ij A_ij c_i c_j=c^T A c`.

After averaging the independent shared sign on each inter-group edge, the
positive kernel weight is therefore exactly

`K_(R_A)(X,X)=exp(-2t k²)cosh(2t c^T A c)`.

The inter-group defect partition of this one fixed spin sector is exactly
this number to the power binom(d,2). Its normalized logarithm tends to

`-t k²+(1/2)logcosh(2t c^T A c)`.

Maximizing over the finitely many common c produces precisely the full
Boolean cap `max_c |c^T A c|`. One may use the centrally symmetric source
(delta_X+delta_-X)/2 instead: the folded kernel is constant on this support,
so its entropic selftransport has the same value. No transport optimization
or numerical Sinkhorn error enters this example.

For k=4 with diagonal +1, compare A_all=J_4 and A_one obtained by negating
only entries12 and21. Their full Boolean caps are respectively16 and12.
For c=1 and t=1 the sector pressures are

`-16+(1/2)logcosh(32)` and `-16+(1/2)logcosh(24)`.

The former exceeds the latter by strictly less than4 and strictly more
than `4-(1/2)exp(-48)`. Thus the better Boolean seed suppresses this
particular actual sector by essentially 4 d² in log partition. This is the
correct seed-cap direction, unlike the correlated Gaussian spectral probe.

LIMITATIONS: There are only 2^k such common-row sign choices, so their
entropy at scale d² is zero. They need not dominate the spin union bound.
Moreover ell=m means full retention; the universal weave matching obstruction
still prevents this example from implying an improved normalized cap. The
calculation concerns inter-group factors; including previously discarded
within-group defects can only decrease the sector weight. The result is a
precise seed-sensitive test sector and a lower floor on the corresponding
inter-group spin-summed positive partition, not a whole-ensemble equality.

## Positive-entropy thickening, with a uniform all-spectra estimate

Fix epsilon in (0,1/2). Instead of the single flat tuple, allow each of the
k d physical rows to differ from c_i 1_m in at most floor(epsilon m) spin
positions. The number of tuples is

`[sum_(j<=floor(epsilon m)) binom(m,j)]^(k d)
 =exp[k² h(epsilon)d²+o(d²)]`.

Let X collect all normalized outgoing blocks on directed inter-group slots,
and Z the corresponding flat blocks. Parseval and full retention give

`||X||² <= k² d²`, `||Z||² <= k² d²`,
`||X-Z||² <=4 epsilon k² d²`.

These are total squared Euclidean norms over all directed slots; no bound
on an individual block is asserted or needed. The estimates survive every
block permutation, because the flat block does not depend on its slot.

For arbitrary two arrays X,Z define F(X)=sum_edges logK_R(X_ij,X_ji).
The exponential-norm part obeys

`|sum_directed ||X_ij||²-sum_directed ||Z_ij||²|
 <= ||X-Z|| (||X||+||Z||)`.

Since logcosh is 1-Lipschitz, symmetrize the difference of each bilinear
term as one half of `(X-Z)^T R(X+Z)` in the two orientations. Cauchy then
gives the SAME bound, times t, for the entire logcosh contribution. Hence

`|F(X)-F(Z)| <=2t ||X-Z|| (||X||+||Z||)
                  <=8t k² sqrt(epsilon) d²`.

This is a deterministic global bound even for concentrated spectra. It is
not the invalid claim that deleting few energetic edges has small log loss.
It controls a genuine total spectral L2 perturbation of a flat reference.

Consequently the log of the spin-summed inter-group partition restricted to
this Hamming neighborhood, divided by d², lies between

`-t k²+(1/2)logcosh(2t c^T A c)+k²h(epsilon)
                       +/- 8t k²sqrt(epsilon) + o(1)`.

For two seeds with a strict flat-sector pressure gap, sufficiently small
fixed epsilon preserves a strict separation while adding strictly positive
entropy density. At k=4,t=1 the allplus/one-negative-edge gap is almost4;
epsilon=10^-6 gives error at most .128 for each seed, far smaller than the
gap. This remains an upper AND lower estimate only for the specified
positive-entropy sector. It neither establishes its dominance over the
remaining spins nor removes the full-retention matching obstruction.
