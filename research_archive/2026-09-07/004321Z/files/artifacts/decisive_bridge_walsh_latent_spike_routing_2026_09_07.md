# Finite latent means can be routed on diagonal spikes in a Walsh ensemble

Date: 2026-09-07. Status: detailed lower-bound theorem and proof for the
specified WALSH + FULL COLUMN PERMUTATION ensemble. It is not a theorem
about arbitrary randomized Hadamard bases or actual signing minima.

## 1. Ensemble and proposed obstruction to a seed-transfer upper

Fix a power-of-two seed order k, a symmetric full sign matrix A with A_ii=1,
a density p in (0,1], and t>0. Let m=k n tend through powers of two, with n
even. Each of n groups uses a Walsh Hadamard of order m with an independent
uniform permutation of ALL m columns. Use a common uniform ell=floor(pm)
row selector within each group. All k retained physical spin rows are summed
over. Between group pairs use one independent shared sign, giving the exact
positive kernel K_(R_A) on the normalized outgoing k-by-k blocks.

Let Z_spin denote this spin-summed, selector/column-permutation/sign-averaged
FULL positive defect partition. The exact self-block routing at the end of
section 3 restores within-group factors; without that extra routing the same
argument proves the inter-group version.
Write nu_p=(1-p)delta_0+(p/2)(delta_(-1/sqrt(p))+delta_(1/sqrt(p))). Then

`liminf n^-2 log Z_spin >= k² [p log2+E_t(nu_p)]`.

The mechanism is a lower bound for the ENTIRE spin-summed partition, not
an upper bound for a chosen sector. It prevents this precise Walsh ensemble
from obtaining a first-moment certificate smaller than the scalar E value
by inserting a fixed seed A. It does not identify its pressure exactly:
joint latent channels could make it larger. It also does not rule out the
more elaborate recursively randomized bases used in the campaign's best
upper construction. Walsh character algebra is essential below.

## 2. Joint channel, common mask, and exact entropy normalization

At one physical row let Y=(M sigma_1,...,M sigma_k), with M Bernoulli(p)
and independent fair sigma_i. Thus H(Y)=h(p)+pk log2. Start with any finite
symmetric scalar channel L_i|Y_i, and apply it independently across i
conditional on Y. Then

`I(Y;L_1,...,L_k) <= sum_i I(Y_i;L_i)=k I_scalar`.

This follows from H(L_vector)<=sum H(L_i) and conditional channel
independence. Joint conditioning also gives

`E Var(Y_i|L_vector) <= E Var(Y_i|L_i)=p D_scalar`,

where D_scalar refers to the normalized scalar source nu_p.

Symmetry lets us write the joint labels as (T,epsilon_1,...,epsilon_k),
where T carries the possibly correlated absolute labels and epsilon_i are
independent fair row orientations, independent of T. A zero signed label
can be given a redundant independent orientation. Conditional means have

`f_i(T,epsilon)=epsilon_i a_i(T)`.

In particular their cross second moments vanish. The average residual
covariance is diagonal; row symmetry makes all diagonal entries equal to
some D_joint<=p D_scalar. Correlation of the labels through the COMMON
MASK does not destroy this property.

First approximate the finite amplitude-label distribution and conditional
source law by rational/dyadic ones, retaining the sign symmetries, with
arbitrarily small changes in information, distortion, and p. Limiting p
can be restored at the end by O(delta m) physical-row repairs. Encode T by
r independent fair bits, mapping their 2^r states to the desired dyadic
amplitude distribution. Redundant encoding adds no information about Y.

Place these r bits and the k orientation bits on independent Walsh
characters of the physical row index g. Their joint law over all m rows is
EXACTLY uniform once m is sufficiently large. For each label cell impose a
conditional empirical type of the physical vector Y(g) approximating its
target conditional law. This keeps all coordinates simultaneously zero or
simultaneously nonzero, hence defines a COMMON selector. The weighted number
of (selector, spin) configurations with these types has exponent

`m[H(Y|L_vector)-h(p)]
 =m[pk log2-I(Y;L_vector)]`.

Indeed the conditional multinomial count gives m H(Y|L); dividing by the
uniform selector count binom(m,ell) subtracts m h(p). This subtraction is
essential and is not an additional arbitrary entropy penalty.

## 3. Disjoint Walsh supports and diagonal routing

Expand each mean function f_i in the label Walsh basis. Its support consists
only of frequencies containing orientation bit epsilon_i and NONE of the
other orientation bits, together with an arbitrary subset of the r shared
amplitude bits. Thus these k Fourier supports are PAIRWISE DISJOINT, even
though the mask and amplitudes are shared. Let S=2^r, padding coefficients
with zeros. Label the support of row i by (i,s), s in [S].

Choose S disjoint perfect matchings of the n group vertices, one per s;
these exist for n>S. In each group, condition the full column permutation
to place Walsh frequency (i,s) in column i of the block belonging to its
s-matching partner. This prescribes only kS columns. Its probability is
at least m^(-kS), so over n groups its log cost is O_(k,S)(n log m)=o(n²).
All unprescribed columns remain uniformly permuted.
Decompose that residual permutation into uniform unordered block formation,
independent within-block ordering, and finally a uniform permutation of the
formed blocks among the unreserved destinations. The good empirical-profile
event in section 4 depends only on formation and not on the last destination
permutation. Conditional on the profile, those destination permutations are
still uniform and independent across groups, as required by positive transport.

Exact conditional physical types fix ALL Fourier coefficients in the label
span. On the prescribed block s, the mean spectral matrix is DIAGONAL:
entry ii is (m/sqrt(ell)) times the s-th Fourier coefficient of f_i; the
offdiagonal entries are exactly zero. The same coefficients occur at both
ends of the matching edge. Therefore R_A acts on this block only through
A_ii=1, and the two spike matrices match exactly. Each such folded kernel
is at least 1/2. All S n/2 matching edges thus cost only O_(S)(n), independent
of the offdiagonal seed signs. This is the actual rare-spike compatibility
argument; no deletion of energetic unmatched spikes is being performed.

There is an exact way to include all within-group defects as well. The full
label Walsh span has 2^(r+k) columns, while the union of nonzero row-mean
cosets has at most k 2^r columns. Every remaining label-span column has
coefficient zero in EVERY physical row, by the exact conditional types.
For k>=2 there are at least k such columns; for k=1 at least one exists.
Prescribe k of these common zero columns as the group's self block. This
adds only k prescribed columns to the permutation event. Consequently
`h_(i,alpha)((j,alpha))=0` for all i,j, so every within-group squared
defect vanishes identically, regardless of the within-group signs. The
routing cost remains O_(k,S)(n log m), and no energy tail estimate is needed.
For parity/rounding impose exact sign-orbit conditional counts on sufficiently
divisible Walsh orders; general densities follow by finite-type approximation.

## 4. Residual spectrum is an isotropic Gaussian empirical profile

Remove the S prescribed blocks from each group. Its remaining normalized
outgoing block profile converges in W2 to N(0,v I_(k²)), where
v=D_joint/p. Here are the needed elementary details.

Before exact-type repair, generate bounded independent physical vectors Y(g)
with their prescribed conditional distribution in each label cell. After
subtracting f(L(g)), every Fourier coefficient outside the finite label
span is a centered sum of m bounded independent vectors divided by
sqrt(ell). Conditional covariance functions depend only on the finitely
many label bits. Thus cross-covariances of two spectral coordinates vanish
unless the difference of their frequencies lies in that finite label span.
For generic distinct column frequencies their covariance is v times the
identity; mean offdiagonal row covariance vanishes by independent orientations.

Uniform random column grouping makes the fraction of exceptional one-block
and two-block choices O_(k,r)(1/m). For all other choices, the joint
characteristic-function Taylor expansion has error O_(k,r,u)(m^-1/2),
because each summand is bounded and scaled by m^-1/2. Hence the one-block
law tends to the stated Gaussian and the two-block law to two independent
copies. The empirical integrals of every bounded continuous test function
therefore converge in probability, by the second-moment argument.

Repair physical conditional counts within their label cells to exact target
counts. The initial deviations are O(sqrt(m log m)) per cell with probability
tending to one. Editing that many WHOLE physical vectors preserves the
common-mask support. Parseval makes the average squared spectral perturbation
o(1). The repair has exp(o(m)) preimages, so the good exact-type configurations
still have the entropy exponent in section 2. Finally Parseval and the fixed
spike coefficients show the residual empirical second moment tends exactly
to k² v. Weak convergence plus this moment convergence gives W2 convergence.
Finite exceptional zero coefficients or blocks do not change the limit.

## 5. Residual transport pressure, and conclusion

For two directed block arrays X,Z and any orthogonal R, the deterministic
global estimate proved in the flat-sector artifact is

`|sum_edges logK_R(X_ij,X_ji)-sum_edges logK_R(Z_ij,Z_ji)|
 <=2t ||X-Z|| (||X||+||Z||)`.

Couple two row permutations using an optimal matching of their empirical
profiles. This gives a two-sided partition comparison in W2; the total L2
matching cost and energies are deterministic. Hence all profiles sufficiently
W2-close to the same Gaussian have residual log partition per n² uniformly
close to that homogeneous Gaussian's transport value. To justify the latter
for unbounded Gaussian profiles, approximate them in W2 by a fixed bounded
finite alphabet, use the exact positive homogeneous theorem, then use this
same global L2 inequality. Deleting S prescribed matching edges causes only
O(n) error for the bounded approximation, and the W2 comparison transfers
the conclusion back. More explicitly, for a bounded alphabet compare a row
with n-1-S slots to a full row by adding S dummy colors. In a uniform full-row
arrangement delete the S missing-edge positions and minimally recolor at most
S surviving slots to restore the target residual counts. A permutation-
equivariant minimal-repair coupling makes the repaired row uniformly arranged
with its target counts. At most O(Sn) endpoint colors and edges change over
the graph, so bounded logK gives O(Sn) log partition error. Thus this is a
row-count coupling argument, not a claim that dropping factors alone compares
the two partition laws. It is not an uncontrolled tail-deletion argument.

The folded Gaussian transport value is k² g_t(v), independent of A, by the
exact folding identity and isotropic Gaussian formula. Combining source
counts, subleading routing cost, matched-spike kernels, and residual pressure
gives

`liminf n^-2 log Z_spin
 >= k[pk log2-I(Y;L_vector)] + k² g_t(D_joint/p)
 >= k²[p log2-I_scalar+g_t(D_scalar)]`.

The last inequality uses decreasing g_t. Approximate arbitrary finite scalar
channels and take their supremum to obtain section 1. Limits are in the
order: fix channel and finite coding approximation; send n to infinity;
then improve the coding approximation and channel. S never grows with n.

## Audit targets and nonclaims

The substantive interfaces to recheck independently are: exact selector
normalization in section 2; disjoint orientation cosets in section 3; and
uniform W2 empirical control after conditional-type repair in section 4.
The result concerns annealed full positive defect pressure for a specified
Walsh ensemble. It does not prove that every optimized recursive Hadamard
ensemble has this lower bound, nor does it yield an actual minimax lower
bound. Its finite-dimensional mean routing uses genuine spectrum/character
structure rather than assuming arbitrary sign gauges preserve the cube.
