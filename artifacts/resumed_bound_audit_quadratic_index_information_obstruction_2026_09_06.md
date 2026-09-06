# Quadratic index maps cannot match independent Fourier clouds

Date: 2026-09-06. Status: exact counting and concentration theorem,
independently reconstructed by the director. This is a scalable
obstruction to ONE specific realization strategy, not R<T and not
nonconvergence of the original signing minima. It leaves correlated
profiles, nonquadratic index maps, and sufficiently long realizable
compositions open.

## 1. Ambient-independent count of quadratic dual-index maps

Fix an index dimension m. Consider any quadratic vectorial dual-bent
system with index permutation sigma on F_2^m, fixing zero. The ambient
bent dimension d may be arbitrarily large. The independently audited
polar identity is

    M_alpha N_(sigma(alpha))=I_d,   alpha!=0,
    M_alpha=sum_i alpha_i M_i,  N_beta=sum_j beta_j N_j.

Every nonzero M_alpha and N_beta is nonsingular. In particular the map
beta -> N_beta is injective.

For each matrix entry (s,t), form the polynomial equation

    sum_(i,j) alpha_i beta_j (M_i N_j)_(st) + delta_st = 0.

Its coefficient vector lies in a fixed binary vector space of dimension

    D=m^2+1,

with basis the m^2 bilinear monomials alpha_i beta_j and the constant
monomial. Let L be the span of all d^2 entry equations. Their common
zero set is exactly M_alpha N_beta=I_d. For each alpha!=0 this has
one and only one solution beta: M_alpha is invertible, N is injective,
and the required sigma(alpha) exists. Thus L uniquely determines sigma,
regardless of the ambient dimension d.

There are at most 2^(D^2) subspaces of F_2^D: encode each by its unique
reduced-row-echelon basis, padded with zero rows to a D-by-D matrix.
Consequently the number of possible quadratic index permutations is

    #Sigma_m <= 2^((m^2+1)^2).                     (1)

This bound already includes every ambient dimension and every choice
of polar matrices. Linear reparametrizations of index coordinates do
not evade it, because their transformed entry equations belong to the
same coefficient space. The bound does not claim that every such
subspace is realizable.

For affine matrix pencils M_alpha=M_0+sum alpha_i M_i and similarly N,
the identical proof applies whenever the target polar pencil is
injective and every source matrix is invertible. The coefficient space
then has dimension (m+1)^2. Without target polar injectivity, this
particular proof of uniqueness must not be asserted.

There is a broader affine count when the polar pencil is not injective.
Write the affine quadratic families as q_x and r_y, and suppose the
target functions y->r_y are injective modulo additive constants. This
holds in particular when r_y=r_0+y dot Psi and the logical index map
Psi is uniform. Besides N_y M_x=I, impose, for every ambient basis
vector e_l,

    r_y(M_x e_l)+r_y(0)=q_x(e_l)+q_x(0).            (1a)

The inverse-polar identity fixes the quadratic part of r_y; the
translation identity fixes its differences on the basis {M_x e_l}.
Thus these equations say exactly that r_y=q_x* modulo an additive
constant, and injectivity modulo constants makes y unique.

The left side of (1a) has degree at most two in x and at most one in
y, because M_x is affine and r_y is quadratic in its ambient argument
and affine in y. In Boolean polynomial normal form all equations
therefore lie in a coefficient space of dimension at most

    D'=(m+1)(1+m+binom(m,2)).

Their span again uniquely determines the index map, giving at most
2^((D')^2)=2^(O(m^6)) maps, independently of ambient dimension. This
slightly weaker count is still subexponential in the number 2^m of
logical labels. It requires the stated function-level injectivity;
mere affine nonsingularity is not silently substituted for it.

## 2. An independent-profile concentration lemma

Let F be an N-by-k random matrix with independent rows X_i in
{+1,-1}^k. Assume E X_i=0; the coordinates within a row may be arbitrarily
correlated. Let Y be any deterministic N-by-k real matrix and U any
deterministic orthogonal N-by-N matrix. For 0<epsilon<k,

    P{N^(-1)||UF-Y||_F^2 <= epsilon}
       <= exp[-N(k-epsilon)/(2k)].                (2)

The same bound holds conditionally if Y and a finite candidate family
of U's depend on arbitrary randomness independent of F.

Proof: put u=N^(-1)||Y||_F^2. The Frobenius norm of F is exactly sqrt(Nk).
If u=0, the event is impossible. Otherwise it implies

    Z=<F,U^T Y>_F >= N(k+u-epsilon)/2.

Writing h_i for row i of U^T Y, the independent centered summands
X_i dot h_i have ranges of length at most 2sqrt(k)||h_i||. Their total
squared range is at most 4kNu. Hoeffding's elementary exponential bound
therefore gives

    P{Z>=N(k+u-epsilon)/2}
       <=exp[-N(k+u-epsilon)^2/(8ku)].

Since (k-epsilon+u)^2/u>=4(k-epsilon), this proves (2). No Gaussian
approximation or independence of transformed Fourier rows is used.

For a fixed family U_N of size exp(o(N)), the union bound yields

    P{min_(U in U_N) N^(-1)||UF-Y||_F^2<=epsilon}->0. (3)

In particular, with probability tending to one the defect remains at
least k/2 simultaneously for the entire family. This conclusion is
uniform in the deterministic target Y, even if its norm is unbounded.

## 3. Application to the independent Walsh-cloud construction

The Fourier-matching majorant theorem starts from two INDEPENDENT
bounded Boolean profile samples F and G, whose Walsh clouds converge
to covariance-matched Gaussian laws. Its polar target is Y=GA, for a
fixed covariance-preserving feature map A. Conditional on G, Section 2
applies to every fixed operator H P_sigma H.

At N=2^m, (1) gives

    log #Sigma_m = O(m^4)=o(N).

Hence no data-dependent selection of a single quadratic dual-bent
index map from this entire class makes the independent-cloud polar
defect vanish with high probability. The count is uniform over the
ambient physical bent dimension, so allowing enormous realizations
for each fixed logical index dimension does not avoid this conclusion.

The exact-balanced version has profiles
F(b,z)=(-1)^b X_z and G(b,z)=(-1)^b Y_z. There are N/2 independent
sample rows. For an operator preserving the active frequency hyperplane,
the preceding argument applies to its induced orthogonal action on those
rows with N replaced by N/2. Even without preservation, compressing
the target to this balanced subspace is a contraction, and the same
Hoeffding proof gives the corresponding positive exponential bound.
Thus the balanced convention does not remove the obstruction.

More generally, a candidate family formed by words of length L_m in
quadratic index maps remains subexponential if L_m m^4=o(2^m), and is
still excluded for independent-profile matching. This is only a word-
count consequence; it assumes no such composition is physically
available, and does not forbid longer realizable words.

## 4. Why this differs from the bilinear matching theorem

The director's deterministic compatibility graph u dot v=0 has dense
edges, no large empty rectangles, and enough length-five augmenting
paths to match any two W2-convergent clouds. That theorem independently
passes. It provides a vast class of compatible permutations.

Quadratic inverse-polar realizability is much more restrictive: the
whole index graph must be cut out by a subspace of bilinear equations
and be single-valued. Its subexponential count is information-theoretic
data absent from the necessary bilinear edge rule. Thus dense local
compatibility does not imply that the matching belongs to the quadratic
realization family.

This is not a bound for specially correlated or jointly designed
Boolean profiles. It does not count arbitrary nonquadratic dual-bent
maps, arbitrary logical row/column sign modifications, or unrestricted
operator circuits. In particular it does not prove R<T, nor rule out
profile-specific finite gates such as affine quadratic Toffoli maps.
The original minimizing sequence M_n/n^(3/2) is not changed by this
obstruction.

## 5. Arbitrary Fourier signs do not rescue independent clouds

The director proposed the following stronger argument, independently
checked here. It removes the need to count possible index-dependent
Arf signs. Let H be normalized Walsh and, for a fixed permutation sigma,
put

    D_sigma=N^(-1/2) min_(D diagonal signs)
                    ||P_sigma H F-D H G A||_F.

Assume F,G are independent row samples of the same symmetric Boolean
law of covariance C, and A^T C A=C. Optimizing the signs row by row
gives the exact squared cost as the sum of squared row norms minus
twice the absolute row inner products. The one-row Walsh CLT and
uniform fourth-moment bounds therefore give, uniformly in sigma,

    E D_sigma^2 -> d(C)
       :=2tr(C)-2E|Z dot Y|,

where Z,Y are independent Gaussian vectors of covariance C. This
expectation is the same for every sigma even at finite N, because
the samples F,G are independent and symmetry makes every Walsh row
have the same marginal law. Diagonalizing C with eigenvalues lambda_i
and applying the triangle inequality gives

    E|Z dot Y|=E|sum_i lambda_i g_i h_i|
       <=(2/pi)sum_i lambda_i.

Hence d(C)>=2(1-2/pi)tr(C)>0; for cut covariance tr(C)=k. No
independence among the Walsh rows within one cloud is asserted.

Here is an elementary concentration proof requiring only product-row
independence before Walsh transformation. If two inputs differ in h
of their 2N bounded rows, orthogonality, the triangle inequality, and
the contraction property of taking the minimum over signs give

    |D_sigma(z)-D_sigma(z')| <= C_A sqrt(h/N),

with C_A independent of N and sigma. Write R=2N. For any measurable
set S in this product space with probability at least 1/2, its Hamming
distance d_H(z,S) is 1-Lipschitz in each coordinate. McDiarmid's lower
tail at zero yields

    E d_H(z,S) <= sqrt(R log(2)/2).

The upper tail then gives

    P{d_H(z,S)>=ell}
       <=exp[-(2/R)(ell-sqrt(R log(2)/2))_+^2].

Apply this first to S={D_sigma<=median} and then to
S={D_sigma>=median}. The preceding square-root Hamming modulus shows

    P{|D_sigma-median|>=t}
       <=2exp[-(2/R)
          (Nt^2/C_A^2-sqrt(R log(2)/2))_+^2].

For each fixed t>0 this is exp(-c_t N). Since D_sigma is uniformly
bounded, its squared expectation differs from its squared median by
o(1), uniformly in sigma. The medians are therefore uniformly separated
from zero. A union bound over exp(o(N)) index permutations still works,
even though every one of the 2^N Fourier signs has already been optimized.

For the exactly balanced doubled profiles, only N/2 input rows per
cloud are independent. The same modulus and proof retain the same
exponential order. Frequencies moved from the active hyperplane into
the inactive one contribute zero cross correlation, so their presence
can only increase the limiting expected squared cost relative to a
permutation preserving the active set.

Thus the joint quadratic-index count and this phase-eliminated
concentration theorem exclude independent-cloud matching even when
arbitrary index-dependent Fourier signs are granted for free. They
still do not exclude deliberately correlated profiles or much larger
families of realizable composite transformations.
