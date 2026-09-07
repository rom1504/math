# Actual critical landscapes at every cap coefficient above one half

2026-09-07. **Proof complete; independent audit requested.** This strengthens
the actual paired construction in `principle_invent_2026_09_07_actual_critical_landscape.md`.
It remains a non-minimizing family: the audited original-sign upper is
strictly below one half. The point is a concrete full-sign construction
with an exact, controlled energy landscape, not a claimed convergence
theorem or an optimizer-specific sharpness statement.

For EVERY fixed c>1/2 there is a sequence of symmetric hollow full-sign
matrices A_n such that

    Q(A_n)/n^(3/2) -> c,

whose only absolute ground states are the two constant words, and whose
absolute deficit has matching critical order around those grounds:

    Q(A_n)-|H_A(x)| >= kappa_c n^(3/2) psi(r)

for every Boolean x, where r=d_H(x,{+1,-1})/n and
psi(r)=r^2 log(e/r). There are scales r_j tending to zero and actual
words x_j for which

    0<Q(A_n)-|H_A(x_j)|<=K_c n^(3/2) psi(r_j).

All constants are positive and independent of n. Thus no fixed power
coercivity kappa r^alpha with alpha<2 holds for this family, despite its
two-center gamma2 square being one and its uniform macroscopic isolation.

The new operation uses a large finite-field block, not a two-spin pair.
Its nonconstant modes carry an almost-isometric full-sign filling, while
a regular graph is inserted in the constant mode. A small residual
constant-mode interaction is controlled, rather than incorrectly set
to zero.

## 1. Coarse hierarchical graphs with a fixed negative spectral gap

Let D=110 and let m=4^ell tend to infinity. Use the exact dyadic-community
construction of Sections1--3 of the paired artifact, with the following
modification of each initial community piece.

Start with its100-regular bipartite double cover, whose second largest
eigenvalue is at most21 and whose bipartition classes have equal size.
Within EACH bipartition class add the simple10-regular cyclic graph
joining each vertex to its five successors and five predecessors. These
internal edges do not overlap the bipartite edges. On the span of the
two part-constant vectors, the resulting adjacency has eigenvalues110
and -90. On the orthogonal complement its norm is at most21+10=31.
Consequently the piece has degree110, second eigenvalue at most31,
smallest eigenvalue at least -90, and cut expansion at least39. It
retains the original bipartite perfect matching for the switches.

Use communities of masses r_j=2^(-j-1), 1<=j<=J=floor(log_2(m)/4),
and a core containing the remaining at least half of the vertices.
With a=1/1000, remove ell_j=ceil(a m psi(r_j)) disjoint matching edges
from community j and the same number of previously unused core matching
edges, then reconnect them crosswise. Exactly as in the paired artifact,
this gives a simple EXACTLY D-regular connected graph G with external
degree at most one and

    cut_G(S)>=c_H m psi(|S|/m),  |S|<=m/2,
    c_H=9a/32,                                      (1)

and

    cut_G(C_j)<=4a m psi(r_j)                        (2)

uniformly j<=J for large m. Every row of the difference from the
initial block-diagonal adjacency has absolute sum at most two, so

    lambda_min(G)>=-92=-(D-18).                     (3)

There is also a useful Cheeger bound that retains the smallest scale.
Write r_*=r_J and

    h_m=(3a/4)r_* log(e/r_*).

Then, for all large m,

    cut_G(S)>=h_m min(|S|,m-|S|).                    (4)

To verify this rather than infer it from an untruncated asymptotic
profile, use the same majority rounding as in the paired artifact.
If the within-community minority count t is at least |S|/4, the cut
is at least |S|/4. Otherwise the selected rounded satellite union has
a largest community of mass r_max>=3|S|/(8m), and r_max>=r_*.
Its external cut divided by |S| is at least
(3a/4)r_max log(e/r_max), which is at least h_m. This proves (4).
In particular h_m has order m^(-1/4) log m.

The sole imported graph input remains Friedman's expander existence,
read and normalized in the paired artifact from
[Bordenave, Theorem1](https://arxiv.org/pdf/1502.04482). All matching,
degree preservation, negative-gap, and cut-profile steps here are exact.

## 2. A masked regular Hadamard signing

Let H_4=J_4-2I_4 and H_m=H_4 tensor ... tensor H_4 (ell factors).
It is symmetric, has entries +/-1, satisfies H_m^2=mI, and every row
sum is sqrt(m). Its diagonal is the constant (-1)^ell. Set

    R_ij=(H_m)_ij if i!=j and ij is NOT an edge of G,
    R_ij=0 otherwise.

This is a symmetric signing of the complement of G, with

    L_R=||R||op<=sqrt(m)+D+1,
    M_R=||R1||infinity<=sqrt(m)+D+1.                 (5)

Indeed deleting the diagonal costs operator norm one, and deleting the
G-supported entries costs at most their maximum absolute row sum D.
The same deletion estimate applies to every row sum. Set K=L_R+M_R.
The R correction is not assumed random or row-regular.

## 3. Finite-field block and the order subsequence

Let q=9^k. On the finite field F_q define C_{uv}=chi(u-v), where chi
is the quadratic character and chi(0)=0. Since q is1 modulo4, C is
symmetric, hollow, and has all off-diagonal entries +/-1. Elementary
character identities give

    C1=0,  C^2=qI-J.                                (6)

For completeness the off-diagonal identity reduces by translation and
scaling to sum_t chi(t(t-1))=-1. Completing the square reduces that to
sum_t chi(t^2-1)=-1: the number of pairs (t,y) with y^2=t^2-1 equals
q-1 because (t-y)(t+y)=1 and2 is invertible. The diagonal identity is
the sum of q-1 squares of nonzero characters. Thus C has eigenvalues
+/-sqrt(q) on the orthogonal complement of the constant vector.
The matrix C+I is a FULL sign matrix, including its diagonal, and has
constant row sum one.

Choose a subsequence of m=4^ell and q=9^k for which

    q/m -> rho=(2c/D)^2.                            (7)

Such a subsequence exists because log9/log4 is irrational. Density of
the fractional parts of k log9/log4 lets k log9-ell log4 approach any
prescribed real number, with k,ell tending to infinity. Irrationality
follows from unique factorization: a nonzero integer power of3 cannot
equal an integer power of2. No prime-gap or degree-repair result is
being imported here. We use only these subsequence orders.

## 4. The actual full-sign matrix and its exact invariant modes

There are m fibres of q physical vertices each, so n=mq. Inside every
fibre use the hollow signing C. Between fibres i and j use

    J_q                         if ij is an edge of G,
    R_ij (C+I_q)                otherwise.

This defines a symmetric hollow FULL sign matrix A: the potentially
problematic matching-coordinate entries in a cross-fibre block are the
diagonal entries of C+I, hence are +/-1, not zeros.

For each Boolean fibre word x_i write

    x_i=a_i1+y_i,  a_i=(1/q)sum_v (x_i)_v,
    y_i perpendicular to1.

Then a_i belongs to[-1,1] and ||y_i||^2=q(1-a_i^2). The constant-fibre
and nonconstant-fibre subspaces are EXACTLY invariant. On normalized
constant vectors, A acts as qG+R. On the nonconstant subspace it acts as

    R tensor (C+I)+I_m tensor C,

whose norm is at most

    L=(sqrt(q)+1)L_R+sqrt(q).                       (8)

Therefore, putting u=sum_i(1-a_i^2),

    H_A(x)=(q^2/2)a^T G a+(q/2)a^T R a+H_bulk,
    |H_bulk|<=(Lq/2)u.                              (9)

By (5)--(7), L/(qD) tends to1/(2c)<1. Fix any mu>0 smaller than
1-1/(2c); for all large orders, qD-L>=mu qD.

The energy at the positive constant word is

    Q_0=(q^2 Dm/2)+(q/2)1^T R1,
    Q_0/n^(3/2)->(D/2)sqrt(rho)=c.                 (10)

We next prove this is the absolute cap, without treating the small
R correction as zero or assuming all maximizing fibres are aligned.

## 5. A soft-spin cut lemma

For a in[-1,1]^m, choose its global sign so that

    r=sum_i(1-a_i)/(2m)<=1/2.

Let Y_i be independent signs with mean a_i, and let T count their
negative entries. Its mean is mr and its variance is at most mr. Hence

    E[T(m-T)]=mr(m-mr)-Var(T)
              >=mr(m/2-1).

Because T(m-T)<=m min(T,m-T), for m>=4,

    E[min(T,m-T)]>=mr/4.                            (11)

The function psi is increasing and convex on[0,1/2]:
psi''(r)=-1-2log r>0 there. Jensen and (11) give

    E psi(min(T,m-T)/m)>=psi(r/4)>=psi(r)/16.        (12)

Applying (1) and (4) to each realized sign word yields BOTH

    E cut_G(Y)>=c_H m psi(r)/16,
    E cut_G(Y)>=h_m mr/4.                           (13)

The link to the soft-spin energy is the exact independent-rounding
identity

    a^T L_G a+D sum_i(1-a_i^2)
       =Dm-a^T G a=4 E cut_G(Y).                    (14)

This is a quantitative rounding lemma with a gain, not an assumption
that every averaged fibre can be replaced by a sign without cost.

## 6. Positive deficit, including the residual Hadamard mode

Let z=(1-a)/2. Then0<=z_i<=1 and sum_i z_i=mr. The exact identity

    1^T R1-a^T R a=4z^T R1-4z^T Rz

and (5) imply

    1^T R1-a^T R a>=-4Kmr.                          (15)

Indeed ||z||^2<=sum_i z_i, so the operator bound on R is sufficient
for the second term. This estimate uses the actual bounded row sums
of the masked regular Hadamard in the first term.

Subtracting (9) from (10), using qD-L>=mu qD, and then (14)--(15),

    Q_0-H_A(x)
       >=(mu q^2/2)[a^T L_G a+D u]-2qKmr
        =2mu q^2 E cut_G(Y)-2qKmr.                 (16)

Here q h_m/K tends to infinity, since q has order m, h_m has order
m^(-1/4)log m, and K has order sqrt(m). Thus one half of the positive
term in (16), using the second estimate in (13), absorbs the entire
negative term for all large m, uniformly r. The other half, using the
first estimate in (13), gives

    Q_0-H_A(x)>= (mu c_H/16)q^2 m psi(r).           (17)

The physical Hamming fraction of x from its nearer constant word is
EXACTLY r: averaging the fibre means does not change this distance.
Also q^2m/n^(3/2)=sqrt(q/m) tends to the positive constant sqrt(rho).
Thus (17) is the required uniform positive-polarity critical barrier.

## 7. Absolute polarity and matching actual soft states

The negative spectral gap (3) and the bulk estimate (9) give

    Q_0+H_A(x)
       >=(q^2/2)[18||a||^2+(D-L/q)u]-(qK/2)m
       >=(q^2m/2)min(18,mu D)-(qK/2)m.              (18)

Since K/q tends to zero, this is a fixed positive multiple of q^2m,
hence of n^(3/2). Combining (17)--(18) proves absolute coercivity and
shows that Q(A)=Q_0 and the only absolute grounds are the two constant
words.

For the word x_j constant within fibres and negative exactly on C_j,
the bulk vanishes and

    Q_0-H_A(x_j)=2q^2 cut_G(C_j)+2q cut_R(C_j).      (19)

The signed R cut satisfies |cut_R(C_j)|<=Kmr_j by the same indicator
calculation as (15), now applied in both directions. Consequently

    Q_0-H_A(x_j)
       <=8a q^2m psi(r_j)+2qKmr_j
        <=(8a+o(1))q^2m psi(r_j),                  (20)

uniformly j<=J, because K/[q r_*log(e/r_*)] tends to zero. The deficit
is strictly positive by (17) and is far below Q_0, so H_A(x_j)>0 and
the same estimate holds for the absolute deficit. Since r_J tends to
zero, these are genuine critical-scale soft collective excitations.

## 8. Scope and relation to the campaign

This theorem supplies actual critical landscapes at every fixed cap
coefficient above one half, arbitrarily close to that threshold. It
does not supply such a family below one half, where the optimum is
already known to lie. The parameter mu, and thus the coercivity
constant, degenerates as c decreases to one half.

The mechanism is constructive and distinguishable from a relaxation:
large quadratic-character blocks fill all cross-fibre entries with
signs, a bounded-row-sum Hadamard remainder controls the residual
constant-mode interaction, and the soft-spin rounding lemma controls
ALL physical words. The dimension and finite-field subsequence are
part of the construction; no all-orders statement is asserted.

### Enlarged centers remove this particular critical obstruction

There is an important distinction between sharpness around the TWO
actual grounds and sharpness of a criterion that may choose a larger
low-complexity center family. The present family proves only the former.
Indeed let F contain ALL fibre-constant words. Its cardinality is2^m,
but its gamma2 square is at most m=o(n): use one coordinate feature for
each fibre, with each physical column equal to that fibre's unit vector.
Each center row is then represented by a sign vector of norm sqrt(m).

The proof even gives linear coercivity transverse to F. After absorbing
the R term in (16), the remaining term is mu q^2 E cut_G(Y), which by
(14) is at least mu q^2 D u/4. Writing d=d_H(x,F), one has

    d=(q/2)sum_i(1-|a_i|),  u>=2d/q.

Thus the positive deficit is at least (mu qD/2)d; the leading negative
gap supplies the same form for the absolute deficit after reducing its
constant. Since qD has order sqrt(n), this is a uniform linear barrier
in the normalized distance from F. The previously proved actual
balanced-center augmentation therefore applies and dilutes this family
while preserving all its old edges.

Consequently these examples are NOT sharpness examples for the
optimizer-specific multiscale covering obstruction, nor a lower bound
on the best possible extension cost. They demonstrate an actual critical
two-ground landscape and, simultaneously, why adaptive enlargement of
the balanced center code matters. Claiming the stronger sharpness from
the two-ground calculation would miss this explicit low-complexity code.
