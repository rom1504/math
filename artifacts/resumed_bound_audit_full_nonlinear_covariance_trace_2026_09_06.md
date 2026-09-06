# Full nonlinear old-response covariance in normalized trace norm

Date: 2026-09-06. Status: proved, independently reconstructed in
`resumed_full_nonlinear_channel_director_audit_2026_09_06.md` and
`resumed_response_full_channel_independent_audit_2026_09_06.md`.
This supplies the covariance/mean-standard-deviation module requested
after the restricted center-update proof. It deliberately does NOT assert
raw covariance convergence in operator norm.

## 1. Statement

Let B=A/sqrt(n-1) be a symmetric hollow signing with fixed ||B||op<=L,
and Q=B^2. Let X_i be a fixed finite family of normalized injective
marked odd-tree fields. Let R be a fixed odd polynomial of the limiting
independent old Gaussian coordinates, with zero first local Gaussian
chaos projection. Decompose it by LOCAL Hermite degree:

    R=sum_{odd k>=3} R_k,       w_k=E R_k(G)^2,
    tau^2=E R(G)^2=sum_k w_k.

For actual Rademacher inputs define

    C_n=E[R(X) R(X)^T],
    S_n=sum_{odd k>=3} w_k Q^(circ k).

Then the proposed conclusion is

    ||C_n-S_n||_*/n -> 0,                             (1)

where ||.||_* is nuclear/trace norm. In particular, with

    u_i=(B C_n B)_ii=E[(B R(X))_i^2],

one has

    n^-1 sum_i sqrt(u_i) >= tau-o(1).                (2)

Both statements are for every fixed finite response and fixed L. There
is no claim that (1) is o_op(1), and no claim that total nonlinear noise
is independent of every old variable before its cross-covariances are
checked separately.

## 2. Kernel inputs actually needed

For each old tree T let X_T,i=I_dT(K_T,i) under Gaussian input, with
K_T,i the exact normalized symmetric injective kernel. The following
three bounds are already available under a fixed operator cap:

* Row Hilbert norms are O(1), and every fixed-root proper marked-slot
  flattening has operator norm O(n^-1/2).
* The global root map i -> K_T,i has operator norm O(1).
* The actual old covariance matrices satisfy
  Gamma_TU=E[X_T X_U^T]=1_{T=U}Q+o_op(1).

The first two are the direct global/fixed-root tree-cut bounds. The
third is proved, including its output-collision correction, in
`fresh_limit_injective_input_gram_2026_09_05.md`, Sections 1--6:
the exact own-spin-marked inputs have covariance delta_TU I+o_op,
and X_T=B V_T plus an error whose covariance operator is O(1/n).
Gaussian and Rademacher covariance of these injective old kernels is
identical, because each surviving input label occurs exactly twice.

All normalization constants involving the fixed degrees and tree
automorphisms are fixed and harmless in estimates below. They are
reinserted when identifying the whole-branch main term.

## 3. First replace local Hermite monomials by Wick forests

For a normalized local Hermite monomial indexed by multiplicities m_T,
write k=sum_T m_T>=3 and D=sum_T m_T d_T. Its Gaussian version differs
in uniformly vanishing L2 from the pure D-chaos kernel

    P_m,i = I_D(sym tensor_T K_T,i^(tensor m_T))
               /sqrt(product_T m_T!).

This is the local Wick/forest replacement already used in the nonlinear
transport lemma. Whole equal-branch contractions are canceled by the
Hermite polynomials; whole distinct-branch contractions and proper
contractions vanish. Let P be the same finite linear combination of
these forest variables as R is of the local Hermite monomials.

The covariance operator of P is uniformly bounded. Indeed the global
root map of an unsymmetrized product of branch kernels has Gram matrix
equal to a Schur product of their positive semidefinite global Gram
matrices. The Schur multiplier inequality bounds its operator norm by
a constant. Symmetrization and the finite sum preserve this bound.

For the moment work with the covariance C_P of this pure Gaussian
forest vector. Its product-formula expansion is a finite sum of slot
pairings BETWEEN the two forests. No marked slots pair internally
within one forest covariance factor.

## 4. The branch-matching graph of a slot pairing

Fix output roots i,j and a pairing between two forests of equal total
marked degree D. Make a bipartite simple graph whose vertices are the
old branches on the i side and on the j side. Join two branch vertices
when at least one marked slot from them is paired. Record the number
of such paired slots as the positive edge weight. Every branch vertex
has at least one incident edge. A vertex is split if it has at least
two distinct neighbors.

For fixed i,j, this is an ordinary tensor network: each branch carries
its fixed-root kernel and each graph edge contracts its indicated slot
block. The common output labels i and j are fixed parameters in this
entrywise calculation. They are not additional contracted slots.

An entire closed tensor-network contraction has absolute value at most
the product of its tensor Hilbert norms: repeatedly merge two connected
tensors, contracting all their shared indices, and use Hilbert
Cauchy--Schwarz. No uncontracted self-loop is created because all shared
indices are contracted at the merge.

More sharply, contracting a split branch A with any one neighbor B
along their shared slot block gives a merged Hilbert norm at most

    ||flatten_shared(A)||op ||B||_HS=O(n^-1/2).

The cut of A is proper because another neighbor remains. This is the
basic gain used below.

## 5. Every genuinely partial pairing has normalized Frobenius decay

There are three component types.

1. A component consisting of one edge joins two whole branches. Summing
   its internal slot bijections gives an old covariance factor Gamma_TU.

2. A nontrivial star component has exactly one split branch, the center.
   Its contraction is O(n^-1/2): contract one leaf into the center using
   the proper-flattening gain and bound everything else in Hilbert norm.
   Its number p of leaves is ODD and at least three. Indeed the center
   marked degree is the sum of the leaves' marked degrees, all of which
   are odd; the center degree is odd as well.

3. Every connected component which is neither one edge nor a star
   contains a simple P4 (a path on four distinct vertices). This follows
   for connected bipartite graphs: a graph with no P4 is a star. Choose
   the two disjoint outer edges of this P4. Each inner endpoint is split
   because it also meets the middle edge. Contract along the two outer
   edges FIRST. The contractions use disjoint pairs of tensors and each
   gains O(n^-1/2). Bound the rest of the network by the product of
   Hilbert norms. This component's entry is therefore O(1/n).

These estimates are uniform in i,j. They also cover repeated branch
types; distinctness of types is not being assumed.

If a pairing has a component of type 3, its entire entry is O(1/n).
If it has at least two components of type 2, the same is true. In
either case the resulting n-by-n matrix has Frobenius norm O(1), hence
nuclear norm at most O(sqrt(n)).

The only remaining partial case has exactly one nontrivial star and
all other components are whole-branch edges. There MUST be at least
one such whole edge. In fact if the star center is on the left and
has p>=3 leaves on the right, the left response has 1+t branches and
the right response p+t branches, where t counts the other whole edges.
Both local Hermite degrees are at least three, so t>=2.

The star gives an entrywise O(n^-1/2) matrix, multiplied entrywise by
the whole-branch Gamma factors. Keep one Gamma factor, whose operator
norm is O(1) and thus whose Frobenius norm is O(sqrt(n)); all other
factors are entrywise bounded. Their total product consequently has
Frobenius norm O(1) again.

Thus EVERY pairing which is not a collection of whole-branch matches
has Frobenius norm O(1), with constants depending on the fixed forest
degrees. There are finitely many pairings. This proves that their
sum has nuclear norm o(n), but does not prove its operator norm
vanishes. A dense O(1/n) entry bound alone would not justify the latter.

## 6. Whole-branch pairings give the Schur main term

For whole-branch pairings, sum all internal slot bijections of each
paired branch first. This produces its exact Gamma_TU covariance,
including its d_T! factor. The remaining count is exactly the Wick
pairing count of the LOCAL Gaussian Hermite monomials.

Replacing Gamma_TU by delta_TU Q costs o_op(1) in each finite Schur
product. To justify this step even for cross-type Gram factors, note
that a matrix Gamma_TU has a Gram factorization
Gamma_TU(i,j)=<a_i,b_j> with uniformly bounded row-vector norms.
Its Schur multiplier norm is therefore uniformly bounded, by tensor
compression. Thus Schur multiplication by other such factors preserves
an o_op(1) error; entrywise smallness alone is not used.

For two distinct local normalized Hermite monomials the whole-branch
main term is zero: at least one branch type must mismatch. For identical
monomials with total local degree k, the number of matching branches
of each repeated type is m_T!, exactly canceled by the two factors
sqrt(m_T!). The main term is Q^(circ k). Consequently

    ||C_P - sum_k w_k Q^(circ k)||_*/n -> 0.

The main-term error is o_op(1), hence o(n) in nuclear norm; the
genuinely partial pairings have nuclear norm O(sqrt(n)) by Section 5.
If two forest terms have unequal original input-chaos degrees, their
covariance is exactly zero and no slot pairing is present.

## 7. Raw responses, repeated slots, and Rademachers

For any two random vectors U,V, covariance Cauchy--Schwarz in nuclear
norm gives

    ||E[U V^T]||_* <= sqrt(E||U||_2^2 E||V||_2^2).

For example, regard the random vectors as Hilbert--Schmidt operators
from the probability L2 space and use ||AB*||_*<=||A||_HS||B||_HS.
Thus an averaged L2 error o(1), together with O(1) averaged second
moments, changes covariance by o(n) in nuclear norm. This applies to
the local Gaussian Wick replacement of Section 3. No covariance-
operator estimate for its raw error is required.

Next delete coincident marked slots from each forest kernel. At a
fixed local root, each coefficient has size O(n^-D/2) and only
O(n^(D-1)) tuples are deleted. The Hilbert error is O(n^-1/2).
Its effect on covariance is therefore o(n) in nuclear norm, by the
same argument. The resulting homogeneous kernels are squarefree, so
their Gaussian and Rademacher covariance matrices agree EXACTLY.

Finally, the direct Rademacher local Hermite-to-disjoint-forest moment
argument gives uniformly vanishing L2 error: the allowed leading
pairings are precisely those encoded by the Hermite subtractions,
while every smaller-label or nonempty-parity pattern vanishes.
Applying nuclear covariance Cauchy--Schwarz one last time proves (1)
for the actual raw R(X_i) on signs.

This sequence explicitly handles the potential operator-sized
accumulation of tiny local collision errors. Such accumulation cannot
be ignored in o_op, but it is negligible in normalized nuclear norm.

## 8. Mean standard deviation without a matrix square-root continuity theorem

Let E_n=C_n-S_n, a symmetric matrix. For the unit-norm row b_i of B,

    |b_i^T E_n b_i| <= b_i^T |E_n| b_i.

For nonnegative scalars a,b, |sqrt(a)-sqrt(b)|<=sqrt(|a-b|).
Consequently

    n^-1 sum_i |sqrt((B C_n B)_ii)-sqrt((B S_n B)_ii)|
      <= sqrt(n^-1 Tr[Q |E_n|])
      <= L sqrt(||E_n||_*/n) -> 0.

This uses only scalar square-root continuity and the exact trace-norm
estimate; no operator-Lipschitz assertion about square root is needed.

If tau>0, S_n/tau^2 is a convex mixture of odd Schur powers of Q.
The already proved odd-Schur mean-standard-deviation inequality gives

    n^-1 sum_i sqrt((B S_n B)_ii) >= tau.

Combine the two inequalities to prove (2). If tau=0, the conclusion
is immediate.

## 9. Remaining distinction before using all tau as a center innovation

The theorem controls the TOTAL variance of the Gaussianized nonlinear
transport B R. To use tau in the center-update bound in place of the
smaller edge-only nonlinear coefficient, one still needs the relevant
cross-covariances with the finite old family to vanish (or another
argument showing that any old-correlated portion remains usable under
the chosen mask). Joint Gaussianity alone only gives a possible old
linear drift plus an independent residual, whose variance may be
smaller than total variance.

This note therefore banks (1)--(2) as a covariance module submitted for
audit. It does not equate total Gaussian variance with fresh independent
innovation variance, and it does not assert the stronger numerical
center update before that last distinction is addressed.

## 10. Subsequent endpoint-graph lemma kills the old drift

The distinction in Section 9 admits a separate actual-tree argument.
Let P_j be one fully injective odd forest with k>=3 child branches and
q marked vertices. Let V_T,a be an exact old input tree: all its
vertices, INCLUDING its root a, carry spins; its root degree is even
and all other degrees are odd. Its output old tree satisfies
X_T=B V_T+E_T with Cov(E_T) operator norm O(1/n).

Claim:

    ||E[V_T P^T]||_F=O(1).                            (3)

Different marked degrees have exactly zero covariance. For equal
marked degree q, a surviving covariance identifies the two marked
sets bijectively. The label a must therefore be a marked vertex of P,
whereas the unmarked forest root j is excluded from BOTH marked sets.
In particular the diagonal a=j is exactly zero. There are q-1 free
labels after fixing a,j. The input V has q-1 edges and the forest P
has q edges, so each pattern has the prefactor n^-1/2 times a bounded
signed graph density.

In its parity-edge graph the only odd vertices are a,j. Every free
vertex combines two odd degrees; root a combines an even V degree
with an odd P degree. Crucially, the parity degree at j is EXACTLY
k>=3: V contains no label j, and P's injective root neighbors are
distinct, so none of its k edges at j can cancel.

If a free-free parity edge survives, the usual bilinear discrepancy
bound gives density O(beta(A)/n^2)=O(n^-1/2). With the prefactor,
the entire matrix has entries O(1/n) and Frobenius norm O(1).

If no free-free parity edge survives, the graph consists of an
optional edge a-j and p length-two paths a-u-j; all other vertices
are isolated. Write e=0 or 1 for the optional edge. The exact degree
at j gives p=k-e>=2. Dropping finite injectivity exclusions changes
the density by O(1/n). Thus its matrix is, up to fixed constants and
Frobenius-o(1) error, either

    n^-1/2 (Q^(circ p)-I),          e=0,
    B circ Q^(circ p),              e=1.

The subtraction of I enforces the exact zero diagonal a=j; the graph
formula was derived at distinct roots. The first has operator norm
O(n^-1/2) by the correlation Schur bound and this diagonal subtraction.
The second does as well: its absolute row sum is at most
n^-1/2 sum_j |Q_aj|^p<=L^2 n^-1/2. Their Frobenius norms are O(1).
This proves (3). The forbidden edge-only graph would have p=0; it
is excluded by the actual forest root degree k>=3, not by a generic
Gaussian orthogonality assertion.

The fully injective forest's covariance operator is O(1), since marked
collision deletion is an orthogonal projection on the common marked
tensor space and cannot increase its global root-map norm. Therefore
Cov(E_T,P) has operator norm O(n^-1/2), hence Frobenius norm O(1).
It follows from (3) that

    ||Cov(BP,X_T)||_F=O(1).                           (4)

Consequently its same-root cross-covariances vanish in averaged absolute
value, by Cauchy--Schwarz on the diagonal. Raw local Wick errors do not
alter that conclusion after B, by averaged L2 Cauchy--Schwarz.

Thus the Gaussianized full nonlinear transport has no surviving
same-root old drift in the averaged limit. More strongly, (4) gives
the squared cross-covariance control over ALL pairs of output roots
needed in the next section.

## 11. General full-noise contraction into exceptional first-chaos transport

Let Z_i=(B P)_i for one fully injective odd forest P of marked degree
q and local branch count k>=3. We must not infer independence from
B X_T merely from (4), because B X_T need not be Gaussian. Consider
every mixed kernel contraction. Proper-Z contractions are negligible
by Z's small proper flattenings (in averaged roots after collision
deletion). Only contraction of ALL q slots of Z into B X_T needs
another argument.

Use the exact-input representation B X_T=Q V_T+small error. Write
V_T,a=N_a H_a under Gaussian input. The even H_a is an exact injective
child forest, excludes N_a, and has marked degree e=d(T)-1. The
variable H_a is therefore ALREADY a pure degree-e Gaussian multiple
integral. Apply the positive covariance-of-squares contraction identity
to this exact pure chaos and the exact degree-q channel component,
before using any approximation by the raw local polynomial h_T(X_a).
That approximation is used only to identify limiting joint moments,
not to assert contraction positivity for a sum of chaoses. The
root-not-hit part of a full Z contraction is

    sum_a Q_ia N_a U_a,i,

where U_a,i contracts the full Z_i kernel into H_a. The creation
inequality bounds its squared L2 norm by a fixed constant times
sum_a Q_ia^2 ||U_a,i||_2^2.

The old-at-a / Z-at-i joint Gaussian approximation, together with the
positive covariance-of-squares contraction identity, bounds the latter
contraction by

    ||U_a,i||_2^2 <= C sum_{old coordinates T}
                            |Cov(X_T,a,Z_i)|^2 + epsilon_i.

Here n^-1 sum_i epsilon_i->0, uniformly in a. To justify that averaged
form, the collision-free transported forest has uniformly small proper
flattenings. Its marked-slot collision deletion has vanishing averaged
Hilbert error and uniformly bounded total Hilbert norms. Thus the
Gaussian-chaos approximation errors may depend on i but have vanishing
average, uniformly over the other old output root a. Replace H_a by
its local Hermite polynomial only in this fixed-moment argument; its
uniform local Wick error is harmless.

By (4), the sum of the squared cross-covariances over all a,i is O(1).
Since |Q_ia|<=1 and n^-1 sum_ia Q_ia^2=O(1), averaging the preceding
bound proves the root-not-hit class is o_L2 in averaged roots. No
unproved uniform independence of every old/new root pair is required.

The ROOT-HIT class has the form

    diag[Q J B],

where J(a,j) is the random partial contraction of the P_j kernel with
one marked slot fixed to a, the remaining q-1 slots contracted into
H_a. It has e-q+1 uncontracted stochastic slots and is zero if this
number is negative. We need only

    E||J||_F^2=o(n),                                (5)

not a small operator norm or a special star derivative formula.

Here is a direct doubled-diagram proof of (5). It suffices to fix each
of the finitely many choices of the marked slot forced to a, and each
finite slot contraction pattern; Cauchy--Schwarz handles their sum.
In the square, take two P copies rooted at j and two H copies rooted
at a. In each pair, q-1 marked P slots are contracted to H slots; the
e-q+1 residual H slots are paired between the two H copies by the
Gaussian inner product. Each P copy has one distinguished marked
vertex fixed to a. All labels a,j are summed in the Frobenius norm.

The denominator has exponent q+e: there are 2(q+e) edge occurrences,
each carrying a normalized B entry. Leading diagrams have

    2 + 2(q-1) + (e-q+1) = q+e+1

distinct labels, including a,j. Thus their nominal absolute size is
O(n). Any additional collision loses at least one free label and
contributes O(1); this includes a residual H label equal to j.

In a leading diagram, the contracted-label sets from the two P copies
are disjoint. Every neighbor of their common root j is either their
fixed marked vertex a or one of those contracted labels; no residual
H slot is a P neighbor. At most the common edge j-a can be doubled.
Each P root has k>=3 distinct neighbors, so at least two, indeed at
least four for a fixed identical slot choice, parity incidences at j
survive. The parity graph is necessarily NONEMPTY.

All its vertex labels are summed here, including a,j. Choose any
surviving edge and apply the bilinear discrepancy estimate after
fixing the other labels. The leading signed density is bounded by
beta(A)/n^2=O(n^-1/2), so its O(n) nominal contribution is O(sqrt(n)).
Lower-label diagrams total O(1), and there are finitely many patterns.
Therefore E||J||_F^2=O(sqrt(n))+O(1)=o(n), proving (5).

Finally,

    n^-1 E||diag[Q J B]||_2^2
       <= L^6 n^-1 E||J||_F^2 -> 0.

This closes the full-Z root-hit case for arbitrary higher odd local
forests. It uses neither generic conditional Gaussianity nor an
unproved generalization of the star's explicit derivative matrix.

All mixed contractions between the full nonlinear forest vector and
the exceptional first-chaos transports now vanish in averaged roots.
The Gaussian product formula therefore gives their needed asymptotic
moment independence jointly with the local old family. The exceptional
first-chaos transports themselves may remain non-Gaussian.

## 12. Sign-input transfer with averaged, rather than uniform, influences

For a full forest channel, marked-slot deletion need only yield a
vanishing AVERAGE of maximum influences; insisting on a uniform
O(n^-1/2) rate is unnecessary. Let

    c_i,a=O(n^-1/2)+||D_a Z_i||_2,
    rho_i=max_a c_i,a.

For the squarefree forest main kernel, sum_a c_i,a^2 is uniformly
bounded and n^-1 sum_i rho_i^2->0 by the proper-flattening argument
and its averaged Hilbert deletion error.

In the unmarked test Y_i f_i, Y_i=(B X_T)_i enters linearly and is
squarefree, while f_i is a fixed polynomial of old fields and the
forest channel. Under every hybrid law, repeated a-coordinate
derivatives of f_i have fixed Lp bounds O(c_i,a^r). The fourth-order
replacement error is bounded by

    C[sum_a c_i,a^4
        + (sum_a ||D_aY_i||_2^2)^(1/2)
             (sum_a c_i,a^6)^(1/2)]
       <= C rho_i^2.

Uniform fixed-degree moment bounds control all constants. Averaging
therefore makes this input transfer vanish. High influences of Y_i
are still harmless because it enters only linearly.

For the direct Boolean root-spin test, remove the own input from the
channel at vanishing averaged cost. Its off-diagonal double-derivative
bound needs only the channel's total influence, not its maximum:

    (n sqrt(n))^-1 sum_ia |Q_ia| ||D_a Z_i||_2
       <= C/sqrt(n).

Use Cauchy--Schwarz in a, bounded row norms of Q, and bounded total
influences. Thus the marked identity extends as well.

Together, Sections 10--12 fill the old-drift and exceptional
first-chaos gaps for the full nonlinear finite-polynomial channel.
They are submitted for the director's independent audit before any
new original-problem numerical constant is claimed.

## 13. Packaged full-noise center consequence, including bounded F

For a fixed odd polynomial F, write F=P1F+R, set

    K_F=U^{-1}P1F=sum_T E[F G_T] h_T,
    tau_F^2=||R||_2^2=E F^2-sum_T(E[F G_T])^2,
    Z=B R(X).

This channel is ALIGNED EXACTLY by the identity B F=B(P1F)+Z; no
separate weighted projection theorem or covariance slope is needed.
The preceding sections give its local Gaussian variance u_i, asymptotic
independence from the old local family, and the restricted test identities

    avg E[M_even(X_i) psi_odd(Z_i) (BF)_i]
      = E M avg E[sqrt(u_i)N psi_odd(sqrt(u_i)N)]+o(1),

    avg E[S_i M_even(X_i) psi_even(Z_i) (BF)_i]
      = E[M K_F] avg E[psi_even(sqrt(u_i)N)]+o(1).

Deterministic bounded weights are allowed. Main forest variances are
uniformly bounded for each fixed polynomial and L; their raw local
Wick errors can be removed in averaged L2. Section 8 gives
avg sqrt(u_i)>=tau_F-o(1).

Here and in the common-density approximation below, u_i is the variance
of the squarefree FOREST MAIN, not necessarily the raw variance of
B R(X). Normalized nuclear covariance control alone does not bound every
raw variance uniformly. First run the fixed-delta, fixed-softsign
argument on the forest main. Its averaged L2 distance from B R(X)
then transfers the test, since softsign at fixed width is Lipschitz.
The mean-standard-deviation bound transfers by
avg |sd(Z_raw,i)-sd(Z_main,i)| <= (avg E|Z_raw,i-Z_main,i|^2)^(1/2).

For any bounded even H in [0,1], the auxiliary-Gaussian smoothing and
softsign argument in Sections 7--9 of
`resumed_bound_audit_restricted_channel_center_update_2026_09_06.md`
now applies with slope one and channel variance u_i. This is an analytic
inequality even when F itself is an unbounded polynomial:

    liminf avg E[H(X_i)|(BF(X))_i|]
        >= E[H(X) E_N |K_F(X)+tau_F N|].             (6)

No Boolean rounding of a polynomial F is involved in (6). The channel
smoothing first adds delta times an independent Gaussian, so its test
variance is bounded below by delta^2; all polynomial L2 approximations
are made with fixed delta and fixed L. The softsign regret inequality
then removes delta and the softsign width uniformly, even at zero u_i.
Convexity in sqrt(u_i), followed by the mean-standard-deviation lower
bound, supplies tau_F on the right of (6).

Now let the actual F be bounded, odd and Gaussian-a.e.-continuous on
finitely many old coordinates, and let H be bounded, even and similarly
continuous, with |F(x)|+H(x)<=1 POINTWISE. Approximate F in Gaussian
L2 by odd polynomials P; these need not be feasible. At fixed L,

    |avg E H|BF| - avg E H|BP||
       <= L (avg E|F(X_i)-P(X_i)|^2)^(1/2),

whose matrix limit is bounded by L||F-P||_2. Orthogonal projections
give ||K_F-K_P||_2<=||F-P||_2 and
|tau_F-tau_P|<=||F-P||_2. Thus the right side of (6) is continuous in
this approximation. Take the matrix limit for each FIXED P first and
then improve P. There is no need to establish a uniform Gaussian law
or a uniform variance bound for the actual nonpolynomial residual
channel B(F-P1F).

For the actual bounded feasible F,H, the exact conditional means

    mu_plus = F+H sign(BF),
    mu_minus=-F+H sign(BF)

give Q(A)/(n sqrt(n-1))>=avg E H|BF| by symmetry, hollow independent
rounding, and the exact half energy difference. Therefore the submitted
full-response center theorem is

    liminf_n M_n/n^(3/2)
      >= E[H(X) E_N |K_F(X)+tau_F N|],              (7)

with tau_F^2=E F^2-sum_T(E[F G_T])^2 and every finite pointwise-feasible
pair just described. The same value holds for every fixed operator cap;
principal Grothendieck deletion then extends it to all signings exactly
as in the frozen minimal numerical proof.

Every construction, polynomial approximation, and smoothing width is
fixed before its matrix-size limit. The theorem does not run an infinite
matrix iteration, does not assert Gaussianity of B(P1F), and does not
turn covariance alone into a conditional expectation. Those were the
three distinct issues addressed by Sections 10--12.

## 14. Infinite Gaussian closure and the banked finite-anchor application

The full-response functional extends continuously to every globally odd
bounded F and globally even 0<=H<=1 on the countable old Gaussian
space, with |F|+H<=1. Put C=U^{-1}P1, K_F=C F, and
tau_F=||(I-P1)F||_2. The map C is a contraction, and

    ||K_F||_2^2+tau_F^2=||F||_2^2.

Writing J_c(F,H)=E H E_N|K_F+tau_F N|, one has for bounded feasible
pairs

    |J_c(F,H)-J_c(F',H')|
      <= (1+E|N|)||F-F'||_2+||H-H'||_2.

Indeed use the absolute-value Lipschitz inequality for the first term;
for the second use Cauchy--Schwarz and
||E_N|K_F'+tau_F' N|||_2<=||F'||_2<=1.

Finite-coordinate conditional expectations preserve odd/even parity
and the JOINT feasibility constraint:

    |E[F|J]|+E[H|J] <= E[|F|+H|J] <= 1.

They converge in L2 as the ancestor-closed finite coordinate sets grow.
Applying the same finite-dimensional Gaussian Markov smoothing to both
functions preserves this joint inequality and produces bounded smooth
finite responses. Thus every infinite Gaussian feasible pair is
approximated by legitimate finite constructions before taking matrix
limits. This does not create a growing-depth matrix iteration.

For the banked anchor mask, V=U g is standard Gaussian, H=1{|V|<=alpha},
W=U H has variance p and covariance E VW=w, and the old response is
F=sign(W)(1-H). Let r=E[F V] and J=E[F W]. Because F is a function of
the Gaussian pair (V,W), its first Gaussian-chaos projection is

    P1F=a V+b W,
    a=(p r-w J)/(p-w^2),
    b=(J-w r)/(p-w^2).

Consequently

    K_F=a g+b H,
    tau_F^2=(1-p)-(a^2+2ab w+b^2 p).

On the center H=1, the full-response certificate therefore reduces to

    E[1{|V|<=alpha} E_N |a g+b+tau_F N|].

The anchor construction supplies g as an explicit finite polynomial in
its finite anchor coordinates and its fixed-point innovation, whose
joint distribution is the standard product Gaussian law. Thus this is
a well-defined finite Gaussian integral even though V itself is a
first-chaos variable supported on infinitely many tree coordinates.
No decimal evaluation or strict improvement is asserted in this note.
