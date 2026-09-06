# Resumed adversarial reconstruction of the 0.4306581794055286 bound

Date: 2026-09-06. Audit started shortly after 01:14 UTC. This note reconstructs
the non-strict numerical lower bound from its actual arguments, without
counting or relying on earlier positive audit verdicts. The exact arithmetic
replay completed at 01:22 UTC. No source report or global document was changed.

## Scope and present verdict

For a real symmetric hollow signing A of order n, set

    H_A(x) = x^T A x / 2,
    Q(A) = max_{x in {+1,-1}^n} |H_A(x)|,
    M_n = min_A Q(A).

The following implication reconstructs without a structural gap:

    liminf_n M_n / n^(3/2) >= 0.4306581794055286.

This audit does not claim convergence, optimality of the Gaussian method,
the stronger strict escape above the entire scalar hierarchy, or an
unbounded-operator nonlinear transport theorem. None is needed below.

The only substantial analytic/combinatorial ingredients are proved in
Sections 1--5 below. The actual finite certificate is checked in Sections
6--8. In particular, an earlier audit's verdict is not a lemma in this
dependency chain.

## 1. Principal spectral deletion, with no lower-bound dependency

Write beta(A)=max_{x,y signs}|x^T A y|. The same maximum over two cubes
[-1,1]^n is equal to beta(A), by separate multilinearity. Since A is
hollow, H_A on the cube is an average of its Boolean values under
independent coordinate rounding. Thus |H_A(u)|<=Q(A) throughout the cube.
For x,y Boolean, let u=(x+y)/2 and v=(x-y)/2. Symmetry gives

    x^T A y = 2[H_A(u)-H_A(v)],

so beta(A)<=4Q(A). No optimizer structure enters this inequality.

Here is the needed diagonal majorant, including its proof. Let
c=asinh(1) and K=pi/(2c). For arbitrary unit vectors u_i,v_j, use odd
tensor powers with absolute coefficient weights from sin(ct), assigning
the alternating coefficient signs to the second vector family. The new
vectors have squared norms sinh(c)=1 and cross inner products
sin(c<u_i,v_j>). Common Gaussian hyperplane rounding then has correlation
(2/pi)arcsin(sin(c<u_i,v_j>))=(2c/pi)<u_i,v_j>, since c<pi/2. Hence

    |sum_ij A_ij <u_i,v_j>| <= K beta(A).

The semidefinite program minimizing sum_i d_i subject to
diag(d)>=A and diag(d)>=-A has dual

    max Tr A(X-Y),  X,Y>=0,  diag(X+Y)=1.

If X and Y are Gram matrices of a_i and b_i, respectively, then
u_i=(a_i,b_i) and v_i=(a_i,-b_i) are unit and have cross Gram matrix X-Y.
The preceding inequality bounds the dual objective by K beta(A).
A sufficiently large scalar diagonal is strictly primal feasible;
X=Y=I/2 is strictly dual feasible. Standard finite-dimensional SDP duality
therefore supplies D=diag(d)>=+/-A with Tr D<=K beta(A). Each d_i>=0
because A_ii=0. Alternatively, attainment of the primal follows directly
from nonnegativity and compactness of a bounded trace sublevel set.

Delete coordinates with d_i>K beta(A)/(epsilon n). Fewer than epsilon n
coordinates are deleted. The retained principal signing A_R, of order r,
satisfies

    ||A_R||op <= K beta(A)/(epsilon n)
              <= 4K Q(A)/(epsilon n).

If Q(A)<=C n^(3/2), then r>=(1-epsilon)n and, for each fixed epsilon>0,
||A_R||op/sqrt(r-1) has a fixed bound depending only on C and epsilon.
No coherence, independence of rows, or involution property is asserted.

Principal monotonicity is exact: fix the retained spins and average all
omitted independent unbiased spins. Every omitted edge has zero mean,
so the average full energy equals the retained energy. Therefore
Q(A)>=Q(A_R), without a bridge penalty.

## 2. The injective odd-tree Gaussian limit

Let T be a finite tree with an external root o of degree one and every
other vertex of odd degree. Set d(T)=|V(T)|-1 and a(T)=|Aut_o(T)|.
The handshake lemma makes d(T) odd. With m=n-1, independent Rademachers
S, and a fixed output i, define

    X_T,i = [a(T)^(-1/2) m^(-d(T)/2)]
            sum_{f injective, f(o)=i}
              product_{uv edge} A_f(u),f(v)
              product_{v != o} S_f(v).

Injectivity includes the output label. Thus the entire finite local
field family is exactly independent of S_i. There is no approximation
or subsequent own-spin removal in this construction.

Assume beta(A)=o(n^2). Expand any fixed mixed moment with D free-vertex
occurrences. A nonzero expectation requires every spin label to occur
an even number of times. Labels equal to the root are prohibited exactly.
Fewer than D/2 free labels cost O(1/n) by absolute counting. Odd D gives
zero. At the leading order, the occurrences form pairs; injectivity
forbids a pair from lying within one copy.

Reduce the quotient graph's edge multiplicities modulo two. Every free
quotient vertex comes from two odd-degree vertices, so its parity degree
is even. The root also has even parity degree: the number of copies is
even when D is even. A nonempty simple Eulerian graph has a free-free
edge. Fix every label except this edge's two endpoints u,v. The remaining
summand is A_uv f(u)g(v) with |f|,|g|<=1. Forbidden fixed labels are
enforced by unary zeros, and u!=v is enforced by A_uu=0. The normalized
contribution is bounded by beta(A)/n^2=o(1), uniformly in output i.

If the parity graph is empty, let V=D/2 be the number of free quotient
vertices and E the number of distinct quotient edges. Connectivity gives
E>=V; even edge multiplicity and D=2V edge occurrences give E<=V.
Consequently the quotient is a tree with every edge repeated twice.
The root edge of one original copy must pair with that of another copy.
Their first free label has exactly these two occurrences. Every incident
edge must therefore occur in both copies. This propagates through the
tree, and cannot switch to a third copy because the two occurrences of
each encountered label are exhausted. Entire rooted-isomorphic copies
pair, and no partial-copy matching survives.

For each such pairing, rooted isomorphism counts cancel the displayed
a(T)^(-1/2) factors. Its free labeling count (n-1)_V/m^V tends to one.
Thus all mixed moments converge to those of independent standard normals
G_T indexed by rooted-isomorphism classes. Fixed higher-moment bounds give
tightness and uniform integrability; Gaussian moment determinacy gives
the joint CLT. The estimates are uniform in the output coordinate.

For bounded normalized operator norm ||A/sqrt(m)||op<=L, the hypothesis
holds automatically: beta(A)<=n||A||op<=Ln sqrt(m). Thus this CLT does
not need the desired numerical lower bound as an input.

## 3. Exact leading classification for the paired energy

Put B=A/sqrt(m). For fixed polynomials F,H of a finite tree-field family,
expand

    n^(-1) E sum_ij F(X_i) B_ij S_j H(X_j).

It is enough to consider monomials. Let D be the number of tree free
vertex occurrences. There are D+1 spin occurrences after including the
explicit S_j, and the bridge adds one edge to the D tree edges.
The denominator is n m^((D+1)/2), asymptotic to n^((D+3)/2).

Nonzero expectations have at most (D+1)/2 spin labels. The output i can
contribute one additional label. A leading diagram therefore has exactly
V=(D+3)/2 total labels: every spin occurrence belongs to a pair, and i
belongs to no spin pair. In particular j is paired with a tree free
vertex, not with i. All lower-label patterns vanish by absolute counting.

If any odd-multiplicity edge remains, BOTH its endpoints are summed:
unlike a root moment, there is no fixed output root in this energy sum.
Fixing other labels gives the same bilinear discrepancy bound as above,
and the normalized contribution is O(beta(A)/n^2). No Eulerian condition
is needed here. This argument controls cyclic and noncyclic parity
graphs equally; it does not silently discard cycles without cancellation.

For an empty parity graph, connectivity and exactly 2(V-1) edge
occurrences force a tree with every edge repeated twice. The bridge's
partner must be the top edge of an F copy: i occurs nowhere internally,
so only F top edges can meet i. The selected F copy's top vertex is j,
and its spin occurrence is paired with the explicit S_j. This exhausts
the two spin occurrences at j.

Now consider all other F copies. Their top edges pair among themselves
at i. From each paired top vertex, the same two-occurrence propagation
as in Section 2 forces whole-copy rooted pairings. They cannot meet j,
whose two spin occurrences are already exhausted. Remove these paired
copies. What remains at root j is the selected F copy's child branches
and the H copies. All their free labels again occur in pairs, so the
same propagation forces whole-branch pairings. Two selected child
branches cannot pair with each other: they belong to one original
injectively labeled F copy. They can only pair to H copies; remaining
H copies pair among themselves.

This argument is the classification, not an assumption built into an
enumeration. In particular it excludes transfers of a pairing to a
third tree, child-child contractions, and a selected branch crossing
back to the i side.

Let m_U(T) count child branches of type U below the top of T. Every
sum_U m_U(T) is even. Define

    h_T(G) = product_U He_{m_U(T)}(G_U)/sqrt(m_U(T)!).

The automorphism identity
a(T)=product_U m_U(T)! a(U)^(m_U(T)) shows that a selected T copy is
a product of normalized child branches with the additional denominator
product_U sqrt(m_U(T)!). Child-child-excluded Wick pairings give exactly
h_T. Selecting the F copy supplies the formal derivative. The limit is

    sum_T E[partial_T F(G)] E[h_T(G) H(G)].             (1)

All labels of a surviving graph are globally distinct at quotient level,
so its count is (n)_V/[n m^(V-1)] -> 1. No additional n, two, or
automorphism factor is missing from (1).

## 4. Gaussian creation and bounded-operator L2 realization

On the countable product Gaussian space indexed by these trees, the h_T
are exactly the normalized finite Hermite monomials of even TOTAL degree.
Indeed, attaching any finite even multiset of smaller rooted trees below
a new top vertex gives precisely one such tree T. Conversely deleting
the root edge gives that multiset. This is a bijection; it is not a
recursive assertion requiring a fixed point of tree types.

Thus h_T form an orthonormal basis of the globally even subspace, and

    U H = sum_T <h_T,H> G_T

is an isometry from that subspace onto the first Gaussian chaos. The
infinite sum converges in L2. It is a Gaussian random variable: its
coefficients are deterministic scalar inner products, not random
nonlinear coefficients. Its variance is ||H||_2^2.

Gaussian integration by parts rewrites (1), for even H, as E[F U H].
Terms outside the finite coordinates used by F have zero covariance
with F, even when U H uses additional coordinates.

Fix L<infinity and a bounded-operator matrix sequence. For bounded
Gaussian-a.e.-continuous finite responses F,H, approximate them by
polynomials P,Q in Gaussian L2, preserving odd/even parity as needed.
The uniform root CLT and fixed higher polynomial moments imply

    max_i E|F(X_i)-P(X_i)|^2 -> E|F(G)-P(G)|^2,

and the analogous statement for H,Q. Boundedness plus higher polynomial
moments gives uniform integrability of every term in this assertion.
For ||K||_(2,n)^2=n^(-1) sum_i E K_i^2, matrix Cauchy--Schwarz gives

    |n^(-1) E[F^T B(S H)-P^T B(S Q)]|
      <= L ( ||F-P||_(2,n)||H||_(2,n)
              + ||P||_(2,n)||H-Q||_(2,n) ).

Multiplication by S is an exact pointwise Euclidean isometry, regardless
of dependencies between coordinates. First take n->infinity for each
fixed P,Q; only then improve the polynomial approximation. On the
Gaussian side E[F U H] is continuous by the isometry and ordinary
Cauchy--Schwarz. No Sobolev norm, derivative estimate, cross-root CLT,
or unbounded-operator extension is required.

For a finite-coordinate even mask 0<=H<=1, let J be a finite
ancestor-closed coordinate set containing its coordinates, and put
W_J=projection_J U H. Use F=sign(W_J)(1-H), or its smooth bounded
sign approximation. The means

    mu_plus = F + S H,       mu_minus = -F + S H

lie in [-1,1]^n pointwise. Independent coordinate rounding preserves
each conditional quadratic energy exactly because B is hollow. The
algebra is

    E H_B(mu_plus)-E H_B(mu_minus) = 2 E F^T B(S H).

Each rounded energy is in [-Q(B),Q(B)], so Q(B)>=E F^T B(S H).
By (1) and its L2 extension the normalized limiting right side is
E[F U H]=E[|W_J|(1-H)]. The omitted Gaussian tail is independent of
the J coordinates and centered. Let J increase after the matrix limit;
W_J->U H in L2. This gives J(H)=E|U H|(1-H).

For an arbitrary infinite-coordinate even mask, finite-coordinate
conditional expectations preserve both parity and [0,1], and converge
in L2. Gaussian smoothing gives bounded continuous finite masks if
desired. The functional satisfies

    |J(H)-J(K)| <= 2||H-K||_2

by isometry, |H|,|K|<=1, and Cauchy--Schwarz. Consequently every such
mask gives, for EVERY fixed L, the SAME lower bound J(H) for every
growing matrix sequence with normalized operator norm at most L.

All masks, finite coordinate sets, and polynomial approximants are fixed
before their associated matrix-size limit. The approximants can depend
on L and target accuracy, which causes no problem. Only their limiting
value J(H) must be independent of L.

## 5. Returning from the operator cap to all signings

The preceding normalized conclusion is initially
Q(A)/(n sqrt(n-1))>=J(H)-o(1); replacing n sqrt(n-1) by n^(3/2)
does not change the limit.

To prove universality, suppose there were an unbounded-order sequence
with Q(A)/n^(3/2)<=J(H)-eta for some eta>0. It is a low-cap sequence
Q(A)<=C n^(3/2) with fixed C. Fix epsilon>0, apply Section 1, and use
Section 4 on the retained growing orders r. These orders need not
contain every positive integer. The theorem applies to every growing
sequence with the fixed operator cap obtained from this epsilon.
Principal monotonicity gives

    liminf Q(A)/n^(3/2) >= (1-epsilon)^(3/2) J(H).

Choosing a sufficiently small FIXED epsilon contradicts the assumed
gap eta. Equivalently, take the matrix limit first and then let
epsilon decrease to zero. No minimizing-sequence assumption, all-order
upper bound, or claimed lower bound is used to obtain regularization.

## 6. Anchored contraction really constructs a same-space Gaussian

Choose a finite ancestor-closed anchor family A. Then every inverse
feature h_T for T in A is a function only of the independent anchor
coordinates G=(G_T:T in A). Let rho_T be fixed, R^2=sum rho_T^2<1,
and s=sqrt(1-R^2). Suppose the finite polynomial h(G,z) is jointly even,
has standard product-Gaussian norm one, is orthogonal to every anchor
feature h_T(G), and has D_z=E|partial_z h|^2<1.

Let S be the unit sphere in first chaos orthogonal to the anchor
coordinates. It is nonempty and complete in the L2 metric, though not
convex. Define R(Z)=U h(G,Z). Every Z in S is standard Gaussian and
independent of G. Therefore the defining norm and orthogonality
conditions imply R(Z) is again in S.

For two innovations Z,Z' in S with correlation q, their JOINT pair is
independent of G, because all variables lie in the same Gaussian first
chaos and both are orthogonal to the finite anchor span. Expanding in
the z-Hermite basis gives

    E h(G,Z)h(G,Z') = K(q) = sum_{ell>=0} w_ell q^ell,
    w_ell>=0, sum w_ell=1, sum ell w_ell=D_z.

For -1<=q<=1, 1-q^ell <= ell(1-q), by the finite geometric sum.
Hence

    ||R(Z)-R(Z')||_2^2 = 2[1-K(q)]
                       <= D_z ||Z-Z'||_2^2.

Banach contraction applies to a complete metric space and does not
require convexity. Its fixed point Z_* remains in the closed Gaussian
first chaos, and remains independent of G. Define

    g = sum_{T in A} rho_T h_T + s h(G,Z_*),
    V = sum_{T in A} rho_T G_T + s Z_* = U g.

Then ||g||_2=1 and V is a standard Gaussian on the original countable
space. This is a deterministic-coefficient Gaussian Hilbert-space
fixed point, not a nonlinear Gaussian transform declared Gaussian and
not a matrix iteration whose depth is taken to infinity with n.

## 7. Finite polynomial and its scalar evaluation

For alpha>0 let chi=1{|sum rho_T G_T+sZ|<=alpha}. Write
p=2Phi(alpha)-1 and v=(rho,s). Its coefficient in the normalized
multivariate Hermite basis at total even degree d and multiindex k is

    c_k = beta_d sqrt(d!/product k_j!) product v_j^(k_j),
    beta_0=p,
    beta_d=-2 phi(alpha) He_(d-1)(alpha)/sqrt(d!)  (d>=2).

This follows by the generating function for a Hermite polynomial of a
unit linear combination and one-dimensional integration by parts for
the central-indicator coefficients.

Retain degrees <=D, delete exactly the anchor feature monomials, divide
each remaining coefficient by a+ell where ell is its innovation degree,
then normalize. The resulting h is a genuine finite polynomial; all
orthogonalities hold exactly. If c_T=<chi,h_T>, the grouped squared
coefficients are

    B_ell = sum_{even d<=D, d>=ell}
              beta_d^2 binom(d,ell) R^(2(d-ell)) s^(2ell)
             - 1_{ell=0} sum_{T in A} c_T^2.

The multinomial theorem gives this identity. Every deleted feature has
innovation degree zero because the anchor family is ancestor closed.
In particular B_ell are nonnegative coefficient-square sums. Set

    N=sum B_ell/(a+ell)^2,
    D_z=[sum ell B_ell/(a+ell)^2]/N,
    w=sum rho_T c_T+s[sum B_ell/(a+ell)]/sqrt(N).

For H=1{|V|<=alpha}, the same-space pair (V,W=U H) is jointly Gaussian
with variances 1,p and covariance <g,H>=w. The standard product law of
(G,Z_*) makes the finite formula above apply exactly. If 0<w<sqrt(p),
put sigma=sqrt(p-w^2). Then

    J(H) = E|W|1{|V|>alpha}
         = 2w phi(alpha)[2Phi(w alpha/sigma)-1]
           + 4sqrt(p) phi(0) barPhi(alpha sqrt(p)/sigma).       (2)

An independent derivation of (2): condition on V=v, giving folded-normal
mean 2sigma phi(wv/sigma)+wv[2Phi(wv/sigma)-1], and integrate over
|v|>alpha. Integration by parts of the linear term contributes
2w phi(alpha)[2Phi(w alpha/sigma)-1] plus a Gaussian tail with coefficient
4w^2/sqrt(p). The first term contributes the same tail with coefficient
4sigma^2/sqrt(p). Their sum is 4sqrt(p), verifying the factor in (2).

## 8. Concrete arithmetic verification and replay

Read completely:

* `computations/fresh_finite_anchor_fixed_point_certificate.py`;
* `computations/fresh_limit_rooted_lower_certificate.py`;
* `computations/fresh_limit_hierarchical_fixed_point_certificate.py`.

The 21 child tuples are distinct, have even lengths, and only reference
earlier indices. This supplies actual distinct ancestor-closed trees,
not arbitrary Gaussian anchors with innovation-dependent inverse images.
The parameters are alpha=361/500, a=17/5, D=200 and

    rho=(8108,-3110,1441,1662,-691,-1088,-758,-874,-444,
         331,638,496,572,502,357,563,392,452,230,286,330)/10000.

The script groups c_k^2 after factoring out 4phi(alpha)^2. Starting at
total degree two already deletes the edge constant; only the 20 other
anchor squares are subsequently subtracted from innovation degree zero.
The removed square for child multiplicities m_j is exactly
He_(sum m_j-1)(alpha)^2 product rho_j^(2m_j)/product m_j!.
The common factor cancels in D_z and leaves the positive factor
2phi(alpha) in the normalized covariance. No infinite Hermite tail is
asserted absent: the construction itself is the finite degree-200
polynomial, so no truncation error in its definition is needed.

Interval primitives round outward on a 60-decimal rational grid.
Addition, multiplication, division away from zero, and integer-square-root
enclosures are valid. Pi is enclosed by Machin's identity and alternating
arctangent tails. The exp and Gaussian-integral series have explicit
tails. For the extended Gaussian integral, initial Taylor terms need not
decrease; x^2/2<degree+2 ensures that the omitted alternating tail does,
which is the only monotonicity the first-omitted-term bound requires.

A fresh execution used `.venv/bin/python -B`, redirected stdout to an
in-memory stream, and intercepted the script's report-file open with a
second in-memory stream. Thus no saved report was overwritten. Parsing
the replay gave exact equality with the existing JSON report. An
additional exact Fraction assertion checked that its lower endpoint
exceeds 4306581794055286/10^16, a stronger target than the script's own
internal assertion value 0.4306.

The replay obtained

    D_z in
    [.990872910662428464868565727570052064530720110786980942487571,
     .990872910662428464868565727570052064530720110786980942487572],

    w in
    [.700415430882898878903274666519237791299666284716117727237748,
     .700415430882898878903274666519237791299666284716117727237840],

    J(H) in
    [.430658179405528602724053804634711026327173238336325190455581,
     .430658179405528602724053804634711026327173238336325190455840].

The script also verifies positive innovation variance, nonnegative grouped
squares, positive covariance, and positive p-w^2. Therefore the contraction
and nondegenerate Gaussian formula hypotheses hold with strict margins.

## Adversarial conclusion

The delicate bridge and child contractions were reconstructed explicitly,
including the reason they cannot change partners or form hidden cycles.
Own-spin independence is exact for the chosen injective fields. The
Gaussian fixed point is a valid construction in closed first chaos, with
deterministic coefficients and anchor/innovation independence maintained.
Bounded-operator L2 approximation is sufficient for the final bounded
mask, so none of the stronger Sobolev or unmarked-response extensions is
a numerical-bound dependency. Finally the deletion fraction is fixed
before the matrix limit and only subsequently tends to zero.

No structural counterexample or gap was found in this minimal proof.
The finite certificate and original absolute same-spin normalization
therefore support the stated lower endpoint. This is not a conclusion
about the unresolved convergence problem or about stronger claims in
the surrounding campaign.
