# An actual full-sign family with a critical quadratic-log landscape

2026-09-07. **PASS; independently reconstructed in
`principle_construct_2026_09_07_actual_critical_landscape_audit.md`.** This is an
explicit existence construction in the original hollow full-sign class,
not a weighted surrogate or merely a sharp example for a scalar upper
envelope. Its normalized cap is a fixed, noncompetitive constant. It
does not prove that minimizers have this landscape, or settle convergence.

Write psi(r)=r^2 log(e/r), with psi(0)=0. There are constants
0<c<C<infinity and a sequence of hollow symmetric full-sign matrices
A_n, n tending to infinity, with the following properties:

1. Q(A_n)/n^(3/2) tends to 101/sqrt(2).
2. The only absolute ground states are the two constant words.
3. For EVERY Boolean x, writing r=d_H(x,{+1,-1})/n,

       Q(A_n)-|H_A(x)| >= c n^(3/2) psi(r).

4. There are scales r_j=2^(-j-1), for j=1,...,J_n with J_n tending to
   infinity, and words x_j at distance r_j n from the constant words,
   such that

       0<Q(A_n)-|H_A(x_j)| <= C n^(3/2) psi(r_j).

In particular this actual family cannot satisfy a fixed coercive power
barrier kappa r^alpha around its two grounds for ANY alpha<2. On the
other hand it is uniformly isolated away from its grounds at every
fixed positive Hamming fraction. The ground-center gamma2 square is one.

The matching order psi does not assert sharpness of the constant80/kappa
in the augmentation theorem, nor optimizer-specific sharpness. It shows
that the critical landscape scale can genuinely occur for actual full
signings of bounded normalized cap.

## 1. Fixed-degree pieces with exact perfect matchings

Fix d=100. For every sufficiently large integer L there exists a simple
d-regular graph on L vertices whose nonconstant eigenvalues have absolute
value at most21. This is the existence consequence of Friedman's theorem
in [Bordenave, Theorem1](https://arxiv.org/pdf/1502.04482), which I read
directly: it applies to uniformly chosen simple d-regular graphs along
every admissible sequence of orders, and 2sqrt(99)+1<21.

Take its bipartite double cover. This is a simple d-regular graph on
2L vertices. Its second largest eigenvalue is at most21: the extra
eigenvalue -d does not affect this one-sided assertion. Thus every set
U of size at most L has at least

    (d-21)|U|(1-|U|/(2L)) >= (79/2)|U|

edges leaving it. The inequality follows directly by decomposing the
indicator of U into its constant and orthogonal parts. The double cover
also has a perfect matching: for a set in one bipartition, counting its
d|U| incident edges proves Hall's condition.

Consequently there are simple d-regular pieces of EVERY sufficiently
large even order with a perfect matching and edge expansion at least39.
Deleting an arbitrary subset of that perfect matching reduces every
cut by at most its smaller side. The retained internal expansion is
therefore at least38; below we only use the weaker constant2.

## 2. Exact-regular dyadic communities

Let M=2^L tend to infinity and put J=floor(L/4). The coarse graph H has
M vertices, partitioned into communities C_0,C_1,...,C_J, with

    |C_j|=M r_j,  r_j=2^(-j-1)  (1<=j<=J),
    |C_0|=M-sum_(j=1)^J |C_j| >=M/2.

All these sizes are even and tend to infinity. On each community choose
one of the pieces from Section1, together with a perfect matching. Set
a=1/1000 and

    ell_j=ceil(a M psi(r_j)).

For all large M, ell_j<=|C_j|/2 and sum_j ell_j<=|C_0|/2. To check these
claims uniformly, psi(r_j)/r_j=r_j log(e/r_j)<1, while the smallest
community size tends to infinity. Also sum_j psi(r_j)<=sum_j r_j<=1/2,
so sum_j ell_j<=aM/2+J<M/4.

For each j remove ell_j distinct matching edges from C_j and ell_j
previously unused matching edges from C_0. Pair the removed edges
arbitrarily. Replacing a pair {u,v} in C_j and {s,t} in C_0 by {u,s}
and {v,t} is a degree-preserving switch. No loops or repeated edges are
created: communities were initially disjoint and every switched vertex
is used only once. Let H be the resulting simple graph.

H is EXACTLY d-regular. Every vertex has at most one external edge.
Each induced community retains internal expansion at least38. There
are exactly

    e_j=2ell_j

edges joining C_j to C_0, no other intercommunity edges, and every e_j
is positive. Thus H is connected. The full cut of C_j is exactly e_j.

## 3. The coarse cut profile, including arbitrary nonaligned sets

Let S be any set of size s<=M/2, and put r=s/M. Round S within each
community to its majority choice, obtaining a union U of entire
communities. Let t=|S symmetric_difference U|, the sum of the minority
sizes in the communities. Internal expansion gives at least38t cut
edges. Since the external graph has maximum degree one, changing t
vertices changes its cut by at most t. Hence

    cut_H(S) >=37t+cut_external(U).                  (1)

If t>=s/4, this is at least s/4, which is at least M psi(r)/4 because
r log(e/r)<=1 for r<=1/2.

Otherwise min(|U|,M-|U|)>=3s/4. One of U and its complement excludes
the core C_0, and is therefore a union of selected satellite communities.
Write v for its normalized size and r_* for the largest selected
community mass. The dyadic sizes imply v<=2r_*, so r_*>=3r/8. Since
psi is increasing on [0,1/2],

    cut_external(U) >=2a M psi(r_*)
                    >=(9a/32) M psi(r).             (2)

The last inequality follows directly if r_*<=r from the logarithmic
factor; if r_*>r it follows from monotonicity. Thus, with
c_H=9a/32,

    cut_H(S)>=c_H M psi(|S|/M)                      (3)

for every S of size at most M/2.

The entire communities give the matching upper bound

    cut_H(C_j)=2ceil(a M psi(r_j))
              <=4a M psi(r_j)                       (4)

for all large M, uniformly j<=J. Indeed M psi(r_J) tends to infinity,
at least at order M^(1/2) log M. This also checks that rounding at the
smallest retained scale does not swamp the asserted profile.

## 4. A clique blowup supplies degree sqrt(order) and a negative gap

Replace each vertex of H by a clique of b=M vertices, and each edge of
H by a complete bipartite graph between its two cliques. Denote this
lexicographic blowup by G. Its order and degree are

    m=M^2,  q=b(d+1)-1=101M-1.

Its adjacency matrix is

    G=(H+I_M) tensor J_b - I_m.

It is connected and q-regular. Its eigenvalues are b(lambda+1)-1 for
eigenvalues lambda of H, together with -1 on the within-clique
zero-sum directions. In particular

    lambda_min(G)>=-b(d-1)-1.                       (5)

This elementary bound is enough; H need not have a uniform global
spectral gap and could even be bipartite.

The cut profile transfers to ALL sets in G. If S occupies a_v vertices
of the clique over v, put f_v=a_v/b and round f_v at1/2. Let F be the
rounded set of coarse vertices and t=sum_v min(f_v,1-f_v). The exact
cut formula is

    cut_G(S)/b^2
       =sum_v f_v(1-f_v)
        +sum_{uv in E(H)}(f_u+f_v-2f_uf_v).          (6)

Write r=|S|/m<=1/2. If t>=rM/4, the first term is at least t/2, hence
at least M psi(r)/8. Otherwise min(|F|,M-|F|)>=3rM/4. Each edge crossing
F contributes at least1/2 to the second sum in (6). Applying (3),
and using psi(3r/4)>=(9/16)psi(r), gives

    cut_G(S)>=c_G b^2 M psi(r),
    c_G=min(1/8,9c_H/32)>0.                         (7)

For the full blowup of a community C_j, equality of the coarse cut
scaling and (4) give

    cut_G(C_j times [b])<=4a b^2 M psi(r_j).         (8)

## 5. An exact paired full-sign compiler

Choose a symmetric hollow matrix R supported on the NONEDGES of G,
with independent fair signs on its nonzero unordered entries. There
exists such a choice with

    ||R||op<=8sqrt(m).                              (9)

Here is a self-contained sufficient estimate. For a fixed unit vector v,
v^T Rv is a sign sum with coefficient-square sum
4 sum_(i<j,nonedge) v_i^2v_j^2<=2. Therefore
Pr(|v^T Rv|>=t)<=2exp(-t^2/4). A1/4-net of the unit sphere has size
at most9^m; its maximum quadratic value is at least half the operator
norm. Setting t=4sqrt(m) makes the union probability at most
2 exp[(log9-4)m]<1 for large m. This proves (9) uniformly for G.

Create n=2m physical vertices, in pairs. The edge within each pair is
+1. Between two pairs use the2x2 block J_2 when their coarse vertices
are adjacent in G. On a nonedge use

    R_ij [[1,-1],[-1,1]].

Every off-diagonal entry is now exactly +/-1, the matrix A is symmetric,
and its physical diagonal is zero. For a physical Boolean word whose
pair entries are (a_i,b_i), put

    z_i=(a_i+b_i)/2,  w_i=(a_i-b_i)/2.

Then z_i,w_i belong to {-1,0,1}, z_i w_i=0, and z_i^2+w_i^2=1. The
energy identity, including the within-pair edge, is EXACT:

    H_A(x)= z^T(2G+I)z + w^T(2R-I)w.                (10)

Since q>2||R|| for all large M, (10) is at most Dm, where D=2q+1,
with equality only when w=0 and z is a constant word. On the other
hand (5) and (9) give

    Dm+H_A(x)
      >=(4b-2)||z||^2+(2q-2||R||)||w||^2
      >=2b m.                                      (11)

Thus no negative-polarity word is an absolute ground, and

    Q(A)=Dm=(202M-1)M^2,
    Q(A)/(2M^2)^(3/2) ->101/sqrt(2).                (12)

The only absolute grounds are the two physical constant words.

## 6. Coercivity for every physical Boolean word

The positive deficit in (10) can be written

    Q(A)-H_A(x)
      =2 z^T L_G z+w^T((2q+2)I-2R)w
      >=2 sum_{ij in E(G)}(z_i-z_j)^2+q t,          (13)

where t=||w||^2 is the number of misaligned pairs. Switch x globally
so that its distance from the positive constant word is at most n/2.
Let s be the number of pairs with z_i=-1. Its normalized distance is

    r=(2s+t)/(2m)<=1/2.

If t>=rm, (13) is at least qmr, hence at least qm psi(r).
Otherwise s>=rm/2 and s<=rm. Each edge leaving the set {z=-1}
contributes at least one to (z_i-z_j)^2. Applying (7) to that set,

    Q(A)-H_A(x)>=2c_G b^2 M psi(s/m)
                >=(c_G/2)b^2 M psi(r).              (14)

Because b=M and n^(3/2)=2sqrt(2) M^3, (13)--(14) prove a fixed
positive multiple of n^(3/2) psi(r). The negative deficit (11) is a
fixed positive multiple of n^(3/2), so it supplies the same bound after
decreasing the constant. This proves the absolute coercivity claim.

Finally, take the physical word that reverses BOTH entries of every
pair in the full blowup of C_j, and no other entries. Then w=0,
r=r_j, and its exact positive deficit is

    Q(A)-H_A(x_j)=8cut_G(C_j times [b])
                 <=32a M^3 psi(r_j).                (15)

It is positive because G is connected. It is much smaller than Q(A)
at these scales, so H_A(x_j)>0 and the same identity holds for the
absolute deficit. Equations (14)--(15) establish matching critical
order, with r_J tending to zero.

## 7. Scope of the result

The paired compiler is a reusable actual-sign operation: a regular
nonnegative graph landscape can be inserted exactly in one pair sector,
while a bounded-operator signing of its complement controls the other.
This differs from replacing zeros by independent signs and hoping their
cut norm is negligible. Their contribution is not small; it is confined
to a disjoint spin sector and dominated there.

For the present construction the cap coefficient is intentionally large.
There is no claimed counterexample to convergence, no near-minimizer,
and no claim that the constant in the critical augmentation bound is
optimal. What is now realized in the original class is the geometric
phenomenon: two low-complexity grounds, a uniform quadratic-log lower
barrier, and actual soft collective excitations of matching order at
arbitrarily small macroscopic scales.
