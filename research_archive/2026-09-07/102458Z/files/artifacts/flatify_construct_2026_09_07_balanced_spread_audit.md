# Balanced spread-port construction: exact audit and seed-transfer loss

This audits the root's new actual-sign construction. Let r be a power of2,
k=r^d and m=1+(k-1)/(r-1). A spread partitions the nonzero Walsh
characters of F2^(d log2 r) into m-1 subspaces of dimension log2 r.
Each port map g to F2^(log2 r) is uniform, and distinct ports are pairwise
independent. Give each opposite fibre one port. Let S be any r-by-r
full-sign seed whose row AND column sums vanish. Set

 C_ij[u,v]=S[g_ij(u),g_ji(v)]

with transpose on reversing the ordered fibre pair. This is symmetric
and full-sign across fibres; the seed itself need not be symmetric.

## 1. Operator audit

Put V_g[u,a]=sqrt(r/k)1_{g(u)=a}. The image under V_g of1-perp has
dimension r-1, and these images are orthogonal for distinct ports. Balance
annihilates the shared constant mode, so

 C_ij=(k/r)V_g S V_h^T

acts independently on the two port spaces belonging to that edge. Thus
the cross operator norm is exactly(k/r)||S||op, and

 Q(C_cross)/N^(3/2)<=sqrt(r-1)||S||op/(2r)+o(1),
 N=mk.

Good internal fibre fillings cost O_r(N^(5/4)). However this operator
bound can NEVER give a coefficient below1/2: rank(S)<=r-1 and
||S||F=r force ||S||op>=r/sqrt(r-1). This is a limitation of this spectral
bound, not an impossibility theorem for the actual construction.

## 2. Matching witness audit

For a matched fibre pair i,j, choose Boolean h,l attaining
beta(S)=max_{h,l}|h^T S l|, and put x_i(u)=h(g_ij(u)),
x_j(v)=l(g_ji(v)). Every other port sees the constant conditional mean of
this spin function, by pairwise independence. Seed balance kills its
interaction there. Hence nonmatching cross edges vanish EXACTLY, not only
after a random sign average. Each matching edge contributes magnitude
(k/r)^2 beta(S).

The common cross polarity can be selected without changing any internal
quadratic energy. Therefore arbitrary internal sign fillings still obey

 Q(C)/N^(3/2)>=sqrt(r-1) beta(S)/(2r²)-o(1).

For the one unmatched fibre, if present, averaging its whole-fibre spin
sign removes remaining interactions; this does not weaken the leading
bound. Internal terms are handled by choosing the common matching
polarity with the sign of their sum.

## 3. Order4 seeds are already rank-two CHSH tiles

Every balanced Boolean vector of length4 is plus/minus one of the three
nonconstant Hadamard characters. Since these characters are linearly
independent, column balance of S forces equal positive/negative counts
of each row type. Four rows therefore use either one type four times
(rank1), or two types once with each sign(rank2).

In the rank2 case, the singular values are sqrt8,sqrt8 and beta(S)=8.
Writing the two row selectors as two Boolean characters gives exactly
the CHSH expression(1/2)[f0,f1] H2 [h0,h1]^T. Thus this is not a new
higher-rank local primitive. The spread port has three nonconstant modes
but the seed uses only two of them. Its matching coefficient is sqrt3/4
and its operator coefficient sqrt6/4. The extra unused port mode means
it is not identical in dimension accounting to the earlier full-rank
rank-two weave, which used k approximately2m rather than3m.

## 4. The usual balanced augmentation loses at leading order

Let A be any hollow n-by-n sign child and fill its diagonal arbitrarily
to a full-sign F. Use the obvious balanced seed

 S=[[F,-F],[-F,F]], r=2n.

Then beta(S)=4 beta(F): differences of two Boolean vectors lie in
[-2,2]^n, and opposite copies attain the bilinear maximum. Also
beta(F)>=max_x |x^T F x|>=2Q(A)-n. The matching lower bound becomes

 sqrt(2n-1) beta(F)/(2n²)
     >=sqrt2 Q(A)/n^(3/2)-O(n^(-1/2)).

Thus retaining an actual child through this balanced augmentation incurs
a fixed factor at least sqrt2 in its normalized cap. For a child near
the active upper coefficient .4936, the lower witness is already about
.698. This is an actual witness, not just an inefficient norm upper
bound. The augmentation cannot deliver favorable child-cap transfer.

## 5. Independent binary features have no higher-rank alternative

We independently checked the root's representation classification. If
u^T R v lies in{+/-1} for EVERY pair of independent Boolean vectors u,v,
then expanding its square forces orthogonal distinct rows, orthogonal
distinct columns, every2-by-2 permanent zero, and ||R||F=1. Disjoint
nonzero entries force the other two entries of their rectangle nonzero.
The support is therefore one complete rectangle. A complete2-by-3
rectangle is impossible: its three row ratios would have to be pairwise
negatives. Orthogonality says the numbers of nonzero rows and columns
both equal rank. Hence the only possibilities are one signed unit entry,
or a2-by-2 block with all magnitudes1/2 and negative sign product.

The spread construction legitimately escapes this classification through
its constrained port codebooks, not through more independent Boolean
features. A seed-Q-based global bound exploiting those constraints has
not been proved here.

## 6. Coset ports remove the augmentation factor, but retain a half floor

Root proposed the following more economical encoding. Let n=2^s,
m=2^d and k=nm. Fibre i has coordinates(a,b) with a in F2^s and b in F2^d.
Use the n-character coset indexed by opposite fibre j and set

 C_ij[(a,b),(c,e)]=F_ac (-1)^(b dot j+e dot i).

The full matrix including i=j blocks is, up to an index permutation,
F tensor W_m, where

 W_m[(i,b),(j,e)]=(-1)^(b dot j+e dot i).

The cross-only construction omits i=j and inserts good internal signs.
This indeed improves the dimension accounting: its spectral coefficient
is ||F||op/(2sqrt n), and matching floor beta(F)/(2n^(3/2)).

Nevertheless an archived half-floor theorem applies exactly on the
power-four seed subsequence. Let H4=J4-2I4. For m=2,
diag(-1,1,1,1) W2 diag(-1,1,1,1)=-H4. Tensoring bit coordinates shows
W_m is switching-equivalent to(-H4)^(tensor d). Consequently, when
n=4^s and m>=2sqrt n is dyadic, the explicit rectangular witness and
same-spin conversion from
`decisive_independent_prescribed_hadamard_half_floor_2026_09_07.md`
give

 q(F tensor W_m)>= (1/2)n^(3/2)m³,                    (3)

where q includes the diagonal. This works for every full symmetric sign
seed F, not just a favorable one. The archive's prescribed Walsh-only
theorem does NOT by itself cover seed orders2 times a power of4.

### Removing the old diagonal fibre blocks is a genuine issue

It is not valid simply to treat the omitted i=j blocks as a small
perturbation. Their individual ranks are low and their total cap can
have leading order. The following witness averaging resolves this issue.

Write w=(i,b) in F2^(2d), and equip this space with the nondegenerate
alternating form B(w,w')=b dot j+e dot i. Every symplectic linear map T
preserves W_m. Thus applying T to the physical coordinate of a witness
leaves its FULL energy in F tensor W_m unchanged. The symplectic group
is transitive on nonzero vectors: extend any nonzero vector to a
symplectic basis and map such bases to one another. Therefore, for any
w!=w', under a uniform symplectic map, their difference becomes uniform
nonzero. Its first d coordinates vanish with probability

 (m-1)/(m²-1)=1/(m+1).

This is exactly the probability that the corresponding edge lies inside
one fibre. Edges with w=w' always lie inside. Their full quadratic energy
has absolute value at most m² q(F), since for each w the n-spin seed
contribution is bounded by q(F).

Take a Boolean witness for(3), choose sigma=+/-1 so that sigma times its
full quadratic energy is positive, and average these automorphic
witnesses. This is the polarity in the absolute objective, not a claim
that reversing the whole spin vector changes a quadratic energy's sign.
There is one whose signed cross-only energy
is at least

 Q(C_cross)>= [m/(m+1)] [(1/2)n^(3/2)m³-m²q(F)].      (4)

This argument changes the witness, not the matrix or the authorized
construction. No arbitrary orthogonal basis is substituted for W_m.

If each fibre is filled by a sign matrix of cap O(k^(3/2)), all internal
fillings together cost O(m k^(3/2)). With N=nm², (4) gives

 Q(C_full)/N^(3/2)>=
   [m/(m+1)] [1/2-q(F)/(n^(3/2)m)]-O(m^(-1/2)).       (5)

For F obtained by completing an actual minimizer's diagonal,
q(F)=O(n^(3/2)), so the coefficient tends to1/2 uniformly as m grows,
including growing seeds n=4^s with m>=2sqrt n. Thus good replacement of
the internal fibre blocks does not rescue this coset-port seed transfer
on that subsequence. The proof does not assert a half floor for all
arbitrary port codebooks or for every odd-logarithm Walsh seed order.
