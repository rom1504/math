# Haar audit of the marked-tree and elementary unmarked identities

Date: 2026-09-06. Scope: the exact statements below, at an averaged random
root. No claim of deterministic uniform-root estimates or operator-Gram
universality. The Haar cap estimate is archived mathematics; the purpose
here is to identify which successful response identities it also obeys.

## 1. Model and the traffic theorem actually used

Let `n` be even and

    U=O diag(I_(n/2),-I_(n/2)) O^T,       O Haar orthogonal.

Then `U^2=I`, `||U||op=1`, and with probability tending to one
`epsilon_n=max_ij|U_ij|=n^(-1/2+o(1))`. Put `C=U-diag(U)` when a hollow
matrix is needed. The operator distance `||C-U||op` tends to zero. The
injective tree formulas below are unchanged by hollowing because their
individual tree embeddings have no loops.

For every fixed connected multigraph `alpha`, the primary
[Gorini--Jones--Kunisky--Pesenti theorem](https://arxiv.org/pdf/2604.11729),
Theorem 4.2 (PDF page 27, printed page 24), gives

    lim n^-1 E z_alpha(U) = product_(cycles sigma) kappa_|sigma|

if `alpha` is a bridgeless cactus, and zero otherwise. Here `z` sums
injective vertex labelings, and `kappa` are the free cumulants of
`(delta_1+delta_-1)/2`. Thus `kappa_1=kappa_3=0`, `kappa_2=1`, and
`kappa_4=-1`; all odd cumulants vanish. Its hypothesis, convergence of
tracial moments in L2, is exact in this model. This is the unpunctured
orthogonally invariant theorem, not the deterministic projected theorem.

A loopless cactus with `v` vertices and `e` edges satisfies

    e <= 2(v-1),

with equality precisely when it is a doubled tree. Indeed each cycle of
length ell contributes ell edges and ell-1 vertices, and ell>=2.

## 2. The old marked-tree local Gaussian law survives

Use precisely the old tree class and normalizations in
`fresh_limit_hierarchical_tree_energy_2026_09_05.md`: root degree one,
every other vertex of odd degree, one iid spin at every nonroot vertex,
injective embeddings, division by the square root of the rooted
automorphism count. Write the resulting field `X_(T,i)`. The edge tree is
`G_i=sum_(j!=i)U_ij S_j`.

For every fixed finite collection of these fields and every fixed joint
moment, averaging over `U`, the input, and a uniform root `I` gives the
same limit as independent standard Gaussians indexed by the rooted tree
types. This holds for Gaussian inputs and for Rademacher inputs.

Proof. Expand the moment as a finite sum over identifications among all
vertex positions. Each individual tree remains injective. Let `D` be the
total number of nonroot positions, equivalently the total number of edge
occurrences. Nonzero input moments require every marked label to occur
at least twice. Since the common root is not a marked label in any copy,
the quotient graph has `v<=D/2+1` and `e=D`. There are no loops inside
any tree copy. Noncactuses vanish by the traffic theorem, while a
loopless cactus requires `D<=2(v-1)`. Thus equality holds throughout:
every marked label has exactly two occurrences and the quotient is a
doubled tree.

At a marked vertex only two tree copies meet. Doubling every incident
edge forces those two copies to agree on each neighboring vertex and
then, by propagation through their connected nonroot parts, on the
entire tree. The copies therefore pair in whole isomorphic rooted
trees. Automorphism normalization makes each whole-tree pair have
weight one. This is exactly Gaussian Wick pairing in the tree-type
coordinates. Input blocks of size four or more lose a free vertex and
vanish, explaining why Gaussian and Rademacher input give the same
answer. All finite moments are controlled by the same finite expansion.

The root spin is exactly absent from all old fields, so the Rademacher
root spin is independent of them even before taking a limit.

The same calculation with two independent, colored Gaussian input
replicas proves concentration of conditional covariance and variance at
a random root. In particular,

    E_(U,I) (E_N X_(T,I)^2-1)^2 ->0,
    E_(U,I) E_N X_(T,I)^4 ->3.

The nonnegative conditional fourth cumulant in each fixed Wiener chaos
therefore tends to zero in mean. This is useful when adjoining other
Gaussian-chaos fields; annealed fourth moments alone would not imply
the needed conditional statement.

## 3. The old paired-energy identity survives

For fixed polynomials `F,H` of the old coordinates,

    n^-1 E F(X)^T U [S H(X)]
        -> sum_T E[partial_T F(Z)] E[h_T(Z) H(Z)].               (6)

Here `h_T` is the same even normalized child-Hermite product as in the
original marked-tree creation identity, and `Z` is the independent old
Gaussian family.

The combinatorial reason is the same as in Section 2. The energy adds
one bridge edge and one explicit marked spin at its second endpoint.
If the old factors contribute `D` edges/spin positions, there are
`D+1` edges and `D+1` marked positions in total. There are at most
`(D+1)/2+1` quotient vertices. Thus the only surviving quotient graphs
are again doubled trees with every marked block paired. The bridge
must double the distinguished top edge of one tree occurrence on the
F-side. Removing that occurrence exposes its child Hermite product on
the H-side; all other whole-tree pairs yield the derivative and
expectations in (6). No fourth free cumulant survives this count.

The diagonal term `i=j` is negligible separately by
`||diag U||op=o(1)` and the polynomial moment bounds. Formula (6)
extends to bounded Gaussian-a.e.-continuous feasible responses by
polynomial L2 approximation and `||U||op=1`.

Consequently the entire old finite marked-tree variational functional
has the relaxation ceiling

    sup_(finite old feasible F,H)
        sum_T E[partial_T F] E[h_T H] <= sqrt(15)/8.             (7)

To see this, use the two feasible means `m_+=F+S H` and `m_-=-F+S H`,
with `|F|+|H|<=1`. Half their normalized expected energy difference is
the left side of (6). The maximum absolute feasible energy bounds that
half-difference. Apply the Haar cap upper bound from
`resumed_convergence_involution_local_ceiling_2026_09_06.md`, Section 4.
In particular, taking a supremum over arbitrary but fixed finite tree
families cannot make the old functional tend to 1/2.

## 4. Elementary unmarked channels also survive locally

For a fixed odd integer `r>=3`, first use Gaussian input `N` and set

    G=UN,       Z_r=U h_r(G),       h_r=He_r/sqrt(r!).

Conditional on `U`, the coordinates of `G` are EXACTLY independent
standard Gaussians. Hence

    Cov_N(Z_r,Z_s)=1_(r=s) I                                   (8)

exactly. A fixed root of `Z_r` is a weighted sum of independent
`h_r(G_j)`; its fourth cumulant is a constant depending on r times
`sum_j U_ij^4 <=epsilon_n^2`. It is therefore asymptotically Gaussian.

Its covariance with an old tree of different input degree is exactly
zero. For an old tree with `d(T)=r`, expand the equal-degree covariance:
the graph is a cone over the old tree, with an additional hub k joined
to the output root and to every marked old vertex. If the hub is an old
vertex or the root, there is a loop and its limiting traffic weight is
zero. If all are distinct, every old nonroot edge and its two hub
spokes form a genuine triangle. A noncactus vanishes; a cactus
containing such a triangle has factor `kappa_3=0` and also vanishes.
The same argument for the squared conditional covariance uses two
independent colored input replicas, glued at the output root. Any
quotient still contains a loop or one of the original genuine
triangles, so the conditional covariance tends to zero in averaged L2.

Combine this covariance conclusion with the conditional old fourth-
cumulant conclusion of Section 2 and the new fourth-cumulant estimate.
The fixed-chaos multivariate fourth-moment theorem gives joint
independent Gaussian limits for any fixed old family and finitely many
`Z_r`. A primary proof and quantitative formulation is Noreddine--
Nourdin, [arXiv:1009.1310](https://arxiv.org/pdf/1009.1310), Theorems
1.1 and 1.5; these were checked directly. Conditioning on `(U,I)` and
subsequence compactness avoids assuming unconditional variables with
random kernels belong to a single fixed Gaussian chaos.

Here is a useful exact input-transfer calculation. Write `u_j` for row
j of U and suppress the fixed Hermite normalization. The r-chaos
kernel of the new channel is

    K_(r,i)=sum_j U_ij u_j^(tensor r).

It has squared norm one. Its squared norm on the slice where the first
two input positions coincide is, for r>=3,

    sum_(j,l) U_ij U_il (sum_a U_ja^2 U_la^2)
                          <u_j,u_l>^(r-2)
      =sum_j U_ij^2 sum_a U_ja^4 <=epsilon_n^2.                  (9)

Its influence at an input coordinate a is

    sum_j U_ij^2 U_ja^2 <=epsilon_n^2.                          (10)

Equation (9) controls Gaussian repeated-position deletion. For the
literal Rademacher expression `U h_r(US)` one must additionally compare
`h_r(u_j.S)` with the evaluation of its squarefree Gaussian kernel.
This is done BEFORE multiplication by U. Expanding a fixed one-variable
Hermite polynomial and grouping repeated input indices gives L2 error
at most `C_r epsilon_n` uniformly in j: pair blocks cancel the Hermite
subtractions; every remaining repeated block has size at least three,
and its squared-coefficient or even-block sum contains an extra factor
`max_a |U_ja|^2`. The other singleton/pair sums are bounded by powers of
`sum_a U_ja^2=1`. Multiplication by U preserves normalized L2 error.
This step is needed; a small Gaussian repeated-position kernel norm
alone does not justify evaluation of that remainder on Boolean inputs.

For the resulting squarefree kernels, fixed-degree input replacement
costs o(1). Old injective degree-d kernels have maximum influence bounded
by `C_d n^(d-1) epsilon_n^(2d)=n^(-1+o(1))`. Standard finite-degree
Taylor replacement (or the corresponding moment expansion) now transfers
the joint law to Rademacher inputs. The own output input is removed at
the same influence cost. There is no Rademacher replacement of a
macroscopic coherent first-chaos atom in this argument.

## 5. A scope boundary needed for the resumed nonlinear theorem

Sections 2--4 establish old local laws, old paired energies, exact
Gaussian-input covariance (8), and old/elementary-new local Gaussian
independence for the Haar model. They already show that old creation
identities plus elementary new local covariance cannot characterize
the flat Boolean optimum at 1/2.

They do NOT by themselves establish the full nonlinear `(U),(S)` tests
of `resumed_bound_audit_restricted_channel_center_update_2026_09_06.md`.
In particular one must still identify the mixed law or tested
projection of `B R_F`, not simply observe that the new channel itself
is Gaussian. The fixed-first-order theorem in the companion artifact
does cover any implementation satisfying its explicit equivariance,
Lipschitz-approximation, and finite-depth hypotheses. Passing from every
injective forest expression to that class has not been asserted here.

This distinction matters: a statement that all flat matrices have cap
1/2 is not falsified by Haar matrices, which are not sign-flat. The
original-class statement actually proved in the companion artifact is
an algorithmic ceiling on exact incoherent involutions, not a ceiling
on their optimum.
