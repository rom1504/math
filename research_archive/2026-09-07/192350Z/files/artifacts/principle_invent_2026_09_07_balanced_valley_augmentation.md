# Actual balanced augmentation and a necessary complexity of near-minimizers

2026-09-07. **PASS: independently reconstructed in
`principle_construct_2026_09_07_balanced_augmentation_audit.md`.**
This is an actual full-sign extension preserving every old edge. It gives
a structural exclusion for asymptotically optimal signings, not convergence.
The elementary version uses signature pairings. A primary vector-balancing
theorem strengthens its range to sublinear factorization complexity.

## 1. Statement

Write H_A(x)=sum_(i<j) A_ij x_i x_j and Q(A)=max_x |H_A(x)|. Suppose
A_n is a sequence of symmetric hollow full signings and there are sets
G_n of Boolean reference patterns and fixed constants kappa>0, 1<alpha<2
such that, for EVERY Boolean x,

    Q(A_n)-|H_(A_n)(x)|
       >=kappa n^(3/2) [d_H(x,G_n)/n]^alpha.               (1)

Assume the sign matrix whose rows are G_n admits a factorization

    G_n=L_n V_n,
    max_i ||(V_n)_:i||_2<=1,
    max_(g in G_n) ||(L_n)_g:||_2<=gamma_n,
    gamma_n^2=o(n).                                     (2)

Then for every sufficiently small fixed epsilon>0 and q=floor(epsilon n)
there exists a symmetric hollow full signing A_n^+ of order n+q, whose
old n by n principal block is EXACTLY A_n, and

    limsup_n [Q(A_n^+)-Q(A_n)]/n^(3/2) <=D(epsilon),
    D(epsilon)=o(epsilon) as epsilon decreases to zero.   (3)

One valid rate, with constants depending only on kappa and alpha, is

    D(epsilon)=O((epsilon log(1/epsilon))^(alpha/(2alpha-2))
          +epsilon^(2alpha/(2alpha-1))+epsilon^(3/2)).     (4)

In particular (2) holds whenever |G_n|=o(n): take V=G/sqrt(|G|)
and L=sqrt(|G|)I. Thus sublinearly many uniformly robust valleys are
enough for the construction. Factorization permits more centers when
their coordinate patterns have low-dimensional structure.

The absolute energy in (1) includes both energy polarities. It automatically
requires G_n to contain every exact maximizer of |H_A|, including its
global reversal. There is no silent assumption about the cross-block
sign after reversing a child.

## 2. Elementary signature-pairing precursor

If G has R reference patterns, partition old vertices by their R-bit
signatures across G. Pair vertices within every signature cell. At most
ell<=2^R vertices remain unpaired. For each new vertex independently,
give the two old vertices in every pair opposite independent random
cross signs; give every unpaired vertex its own independent fair sign.
This is a genuine full +/-1 cross column v.

For any g in G each paired contribution to v dot g vanishes. If x has
Hamming distance d from g, at most d pairs contribute nonzero coefficients
to v dot x, each of magnitude two; the unpaired coefficients have magnitude
one. Thus

    E exp(lambda v dot x)
       <=exp[lambda^2(4d+ell)/2].                         (5)

The leftovers are RANDOM, not paid by a deterministic q ell bound.
Consequently this elementary construction already proves (3) whenever
2^R=o(n), and more generally when the actual number of odd signature
cells is o(n) and log R=o(n). At the centers, its cross contribution
is of scale q sqrt(ell), hence normalized scale epsilon sqrt(ell/n).

## 3. Stronger cross-column distribution from Gram--Schmidt balancing

The imported result is Theorem 1.4 of
[Bansal, Dadush, Garg and Lovett, The Gram--Schmidt Walk (2019)](https://theoryofcomputing.org/articles/v015a021/v015a021.pdf).
For vectors of norm at most one and initial coloring zero, it supplies
a random full sign coloring whose vector discrepancy has MGF bounded
by exp(40||theta||_2^2/2). This is an all-directions MGF statement, not
only coordinate tail bounds. The theorem and its exact convention were
checked in the primary publication.

Apply it to the vectors

    w_i=((V_n)_:i,e_i)/sqrt(2),                          (6)

whose norms are at most one. Let v be the resulting full sign vector.
For x and a nearest g in G, use the test vector

    theta=(sqrt(2)(L_n)_g:, sqrt(2)(x-g)).

Its inner product with sum_i v_i w_i is exactly v dot x, while
||theta||_2^2<=2gamma_n^2+8d_H(x,G). Therefore

    E exp(lambda v dot x)
       <=exp[lambda^2(80gamma_n^2+320d_H(x,G))/2].         (7)

The identity coordinates in (6) are essential: they preserve the small
Hamming-increment control while the other coordinates balance the reference
family. Independent samples of v form the q actual cross columns.
No real-valued cross matrix, scalar covariance rounding, or matrix-sign
replacement is being used.

## 4. The center count is automatically subexponential

For completeness, (2) implies log|G_n|=o(n); this need not be an extra
assumption. The sign family G_n has VC dimension at most gamma_n^2.
Indeed if d coordinate positions are shattered, for every sign vector
epsilon on them there is a row vector l_epsilon of norm at most gamma_n
with <l_epsilon,V_i>=epsilon_i on those positions. Hence

    d <=gamma_n ||sum_i epsilon_i V_i||_2.

Average over independent fair epsilon_i and use Cauchy--Schwarz and
||V_i||<=1 to obtain d<=gamma_n sqrt(d), so d<=gamma_n^2.

The elementary Sauer recursion bounds a binary family of VC dimension
D by sum_(j<=D) binom(n,j): splitting on the last coordinate gives the
recursion f(n,D)<=f(n-1,D)+f(n-1,D-1), with the usual boundary cases.
Thus for gamma_n^2=o(n),

    log|G_n|<=n h2(floor(gamma_n^2)/n)+o(n)=o(n).          (8)

When gamma_n^2 is below one, nonempty sign families already force gamma_n
at least one by Cauchy--Schwarz, so no exceptional denominator is hidden.

## 5. One simultaneous bridge bound for all old and new spins

Let C be the n by q cross matrix with independent columns from (7).
For fixed x,y, with x at distance d from G and y in {+-1}^q,
x^T C y is subgaussian with variance proxy

    q(80gamma_n^2+320d).

At distance d there are at most |G| binom(n,d) old configurations. Union
over every shell d=0,...,n, all old configurations, both tails and all
2^q new spin words gives, with probability at least 1/2, simultaneously,

    |x^T C y| <=T_d,
    T_d=sqrt(2q(80gamma_n^2+320d)
          log[4|G|(n+1)2^q binom(n,d)]).                 (9)

There is no unnecessary confidence term of order n in the logarithm:
constant positive existence probability suffices. By (2), (8), and
log binom(n,d)<=n h2(d/n), the normalized envelope converges uniformly
in r=d/n to at most

    F_epsilon(r)=sqrt(640epsilon r[h2(r)+epsilon log2]).  (10)

Choose any new q by q hollow signing B_q with Q(B_q)<=q^(3/2). Such a
signing exists by the elementary independent-edge tail and union over
all Boolean words; no existing strict upper theorem is required here.
Set A_n^+ to have blocks A_n,C,C^T,B_q. For every x,y,

    |H_(A_n^+)(x,y)|
       <=Q(A_n)-kappa n^(3/2)r^alpha+T_d+Q(B_q).          (11)

This inequality pays the genuine absolute cross contribution. In particular
it respects the child-reversal identity; no cancellation between separately
bounded channels is assumed.

## 6. Why the added cost is smaller than the dilution gain

Since h2(r)<=r log(e/r), (10) is at most a universal constant times

    sqrt(epsilon) r sqrt(log(e/r))+epsilon sqrt(r).       (12)

Split the negative term kappa r^alpha equally between the two summands.
For the first term, take a cutoff

    r_0=C_(alpha,kappa)
          (epsilon log(1/epsilon))^(1/(2alpha-2)).

For r>=r_0, the first summand is at most (kappa/2)r^alpha
once the constant is large and epsilon is small. For r<=r_0, the
function r sqrt(log(e/r)) is increasing, so its value at r_0 bounds
the whole positive excess by
O((epsilon log(1/epsilon))^(alpha/(2alpha-2))).

The ordinary one-variable maximum of C epsilon sqrt(r)-(kappa/2)r^alpha
is O(epsilon^(2alpha/(2alpha-1))). Together with Q(B_q), this proves
(4). Both powers before the child term exceed one precisely as needed:
alpha/(2alpha-2)>1 when alpha<2, and 2alpha/(2alpha-1)>1 always.
All constants are independent of n and of the center family.

First fix epsilon, then let n grow so gamma_n^2/n and log|G_n|/n vanish;
only afterward decrease epsilon. There is no assumption of a growing
center family being uniform in an unproved finite-dimensional limit.

## 7. A necessary geometry of asymptotically optimal signings

Let c_inf=liminf_n M_n/n^(3/2)>0. Suppose a subsequence of signings has
Q(A_n)/n^(3/2)->c_inf and satisfies (1)--(2) with fixed kappa,alpha.
Choose a sufficiently small fixed epsilon so

    D(epsilon)<c_inf[(1+epsilon)^(3/2)-1].

Then (3) yields actual signings along the larger orders with

    limsup Q(A_n^+)/(n+floor(epsilon n))^(3/2)
       <=[c_inf+D(epsilon)]/(1+epsilon)^(3/2)<c_inf,

a contradiction to the definition of c_inf. Thus no liminf-realizing
sequence can have a uniformly subquadratic robust-valley landscape with
sublinear factorization complexity as in (1)--(2).

In particular asymptotic minimizers cannot have only o(n) reference
valleys with a common kappa>0 and exponent alpha<2. They must violate
at least one of: small center complexity, uniform coercivity, or the
subquadratic exponent. This is a structural theorem about the ORIGINAL
full-sign problem obtained from an explicit successful augmentation.

It is not a proof that arbitrary minimizers admit the favorable geometry;
in fact it proves they cannot. The theorem explains a genuine obstruction
to simple seed extension: near-optimal signs must retain enough competing
near-extremal structure to prevent balanced cross columns from making
dilution profitable. No convergence conclusion is asserted.

## 8. A multiscale covering version, stronger than a thick-level statement

The same contradiction holds under the following more flexible hypothesis.
There are fixed C>0, 1<alpha<2 and delta_0>0 such that for EVERY fixed
delta in (0,delta_0], the near-extreme level

    E_(n,delta)={x: Q(A_n)-|H_(A_n)(x)|<=delta n^(3/2)}

is covered by Hamming balls of radius C delta^(1/alpha)n, whose centers
form a family G_(n,delta) with factorization square o_delta(n).
The centers need not be the same at different levels.

To prove this, fix the augmentation fraction epsilon first. Take finitely
many dyadic levels from delta_0 down to delta_min<=epsilon^2. Form the
union G of their center families. A finite union still has factorization
square o(n): stack its finitely many feature matrices divided by the
square root of the number of levels, and multiply the corresponding row
factor by that square root. The VC argument still gives log|G|=o(n).

If a spin has normalized energy deficit e in (delta_min,delta_0], use
the least dyadic level delta>=e. Its distance to G is at most
C(2e)^(1/alpha)n. If e<=delta_min, its distance is at most
C delta_min^(1/alpha)n. Thus on these levels there is the relaxed
coercivity

    e>=c_(C,alpha) [d_H(x,G)/n]^alpha-delta_min.           (13)

Spins with e>delta_0 are harmless for a sufficiently small fixed epsilon:
the identity coordinates in the GS embedding give the universal cross
bound O(sqrt(epsilon))n^(3/2), strictly less than delta_0 n^(3/2).
On the remaining spins, (13) adds at most delta_min<=epsilon^2 to the
previous o(epsilon) augmentation cost. The same dilution contradiction
follows. Only finitely many levels are used at each fixed epsilon; no
uniformity in an infinite depth of covers is assumed.

This profile exclusion is not the generic observation that a whole
fixed-width superlevel has linear gamma2 complexity. Such a superlevel
already contains a Boolean face of dimension Omega(delta n), as proved
in `principle_invent_2026_09_07_thick_level_complexity_caution.md`.
That face fits inside ONE Hamming ball of radius C delta^(1/alpha)n
when delta is small and alpha>1. The multiscale conclusion instead
requires competing centers or collective excursions beyond these local
faces, at some scales relative to the square-root energy radius.

## 9. Explicit full-sign non-vacuity example

The favorable robust-valley hypothesis is not restricted to weighted
matrices. Let q be an odd prime power, n=q^2, and take t=(q+3)/2
one-dimensional directions in F_q^2. For distinct vertices u,v put
A_uv=+1 if u-v lies on a selected direction, and A_uv=-1 otherwise.
Set A_uu=0. This is an actual symmetric hollow full signing.

Each selected line contributes q-1 neighbors, so the constant-vector
eigenvalue is

    d=2t(q-1)-(q^2-1)=2q-2.

For a nontrivial additive character, its sum over a line minus zero
is q-1 on its unique annihilated direction and -1 on every other
direction. Consequently all remaining eigenvalues of A are exactly
q-2 or -q-2. This proves the spectrum directly by finite Fourier
diagonalization; no random spectral theorem is being imported.

For q>=9, the unique absolute maximizing spins are +1 and -1, with

    Q(A)=(q-1)q^2,       Q(A)/n^(3/2)=1-1/q.

Indeed if r=d_H(x,{+-1})/n<=1/2, decomposing x into its constant and
orthogonal components yields

    Q(A)-H_A(x)>=2q n r(1-r)>=q n r,
    Q(A)+H_A(x)>=(q-4)n/2>=q n/4.

Therefore

    Q(A)-|H_A(x)|>=(1/2)n^(3/2) r.

The two centers have a factorization with gamma^2=1. Since r^alpha<=r,
the theorem applies with uniform kappa=1/2 for every fixed 1<alpha<2.
It supplies an actual old-edge-preserving dilution improvement of this
explicit infinite family. The family is not competitive with the best
known cap, and this example does NOT prove sharpness of the alpha=2
threshold. It establishes that the successful augmentation hypothesis
has genuine full-sign realizations.

Regression check: direct integer construction for q=11,13,17 (orders
121,169,289) had the stated constant row sums and the three predicted
eigenvalues to numerical residual below 1e-9. For each matrix, 1000
independent Boolean test vectors satisfied the displayed coercivity
inequality. These checks are diagnostic only; the finite Fourier and
spectral inequalities above are the proof.

## 10. Endpoint multiscale cover obstruction from the critical barrier

The tuned GS construction in
`principle_construct_2026_09_07_critical_barrier_augmentation.md` improves
the limiting bridge envelope to

    sqrt(320epsilon r[h2(r)+epsilon log2]).              (14)

I independently checked its feature tuning, all-directions MGF, uniform
shell limit, the exact coefficient 80/kappa at a kappa r^2 log(e/r)
barrier, and its mesoscopic-excitation consequence. The following is a
covering consequence of the same calculation.

Let c_inf be the original liminf. Fix a>0 with

    40a^2<(3/2)c_inf.                                  (15)

A liminf-realizing sequence CANNOT have, for every sufficiently small
fixed eta>0, a cover of E_(n,eta) by Hamming balls of radius

    a sqrt(eta/log(e/eta)) n,                           (16)

whose center family has factorization square o_eta(n).

Proof. Choose a small geometric mesh ratio 1+theta. For any small fixed
epsilon use only finitely many levels down to eta_min<=epsilon^2, and
form their center union as in Section 8. The tuned feature construction
still applies to this finite union. If a spin has deficit eta between
adjacent mesh levels, it has a center at relative distance no larger
than r_delta=a sqrt(delta/log(e/delta)), where delta<=(1+theta)eta.

As delta decreases to zero,

    r_delta^2 log(e/r_delta)=(a^2/2)(1+o(1))delta.        (17)

The error is uniform on all levels below a sufficiently small fixed top
level. Since r sqrt(log(e/r)) is increasing on [0,1], the first part
of (14) is consequently bounded by

    a sqrt(160(1+theta)(1+zeta)epsilon eta),

for arbitrarily small fixed zeta>0 after choosing that top level small.
Subtract eta and optimize over eta>=0. The result is at most

    40a^2(1+theta)(1+zeta)epsilon.                       (18)

The second part of the split square root in (14) is at most a constant
times epsilon times the square root of the top covering radius; this
is an arbitrarily small additional multiple of epsilon. The bottom
level contributes O(sqrt(epsilon eta_min))=o(epsilon), and the new
child costs epsilon^(3/2). Spins below the top energy level cannot win
once the universal O(sqrt(epsilon)) bridge is less than their fixed
energy deficit. These statements use n-limit first for each fixed
epsilon and finite family of centers.

Under (15), choose theta,zeta and the top level so the total linear
coefficient stays below (3/2)c_inf. The actual extension then contradicts
liminf dilution exactly as before. This proves the endpoint profile
obstruction without asserting a single common center set at every scale.

The audited positive lower bound c_inf>.433 suffices, for example, to
take a=1/10 in (16). This is stronger than the earlier alpha<2 covering
criterion. It remains distinct from the generic thick-level Boolean
face, whose radius is only O(eta)n and therefore fits inside (16) for
small eta. Some additional competing-center or collective-excursion
structure is required at these larger radii.

Sharpness at the level of the envelope only: for f(r)=kappa r^2 log(e/r),
choose r=r_epsilon with r sqrt(log(e/r))=sqrt(80epsilon)/kappa.
Then h2(r)/(r log(e/r))->1 and the epsilon log2 term is negligible,
so the supremum of (14)-f has asymptotic slope exactly80/kappa. Thus
that scalar coefficient cannot be improved without improving the bridge
estimate or using more information. This is NOT an actual full-sign
example proving optimality of the geometric threshold; that question
remains open.
