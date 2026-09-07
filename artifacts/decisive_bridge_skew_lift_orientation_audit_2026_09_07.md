# Skew-orientation doubling: exact interface and finite obstruction

This audits the director's proposed width-to-balanced-cap doubling. It
does not prove or disprove the desired asymptotic orientation theorem.

Let A be a hollow symmetric full-sign matrix of order n and let C be
skew-symmetric with C_ii=0 and C_ij in {+-1}. Define

    L=[[A,C],[-C,-A]].

L is symmetric and has only n missing undirected edges, the matching
between corresponding coordinates in the two copies. Swapping the two
blocks sends L to -L exactly, so Q(L)=W(L). Filling the missing matching
by arbitrary signs changes each oriented endpoint by at most n, hence
costs at most n in cap and at most n in midpoint.

## 1. Exact directed-cut representation

For Boolean x,y put z=(x+y)/2 and w=(x-y)/2. At each coordinate exactly
one of z,w is nonzero, and that value is a sign. Direct expansion gives

    H_L(x,y)=2 z^T(A-C)w.

Orient edge i--j from i to j exactly when C_ij=-A_ij. Then A-C has
entry 2A_ij in the forward direction and zero in the reverse direction.
Writing

    D(A,orientation)=max_(S disjoint T, S union T=[n])
                    max_(z in {+-1}^S,w in {+-1}^T)
                    |sum_(i in S,j in T,i->j) A_ij z_i w_j|,

we therefore have the exact identity

    Q(L)=4D(A,orientation).                               (1)

The desired scale-preserving doubling from the width is precisely

    min_orientation D(A,orientation)
       <= W(A)/sqrt(2)+o(n^(3/2)).                         (2)

It is not enough to cancel the signed sum on one selected extremal cut.
The same orientation must satisfy all cut partitions and all sign profiles.

There is a stronger exact warning against a cancellation argument. Swapping
z and w preserves z^TAw and negates z^TCw. Since both cut directions are
admissible,

    Q(L)=2 max_disjoint(z,w) (|z^TAw|+|z^TCw|),
    D=(1/2) max_disjoint(z,w) (|z^TAw|+|z^TCw|).           (1a)

The proposed improvement therefore requires a JOINT magnitude bound:
profiles with large A response must have small C response and conversely.
Even a cleverly chosen signed cancellation is defeated by cut reversal.

For every orientation, D>=W(A)/2. To see this, width equals the maximum
full bilinear norm across a vertex cut: the lower inequality follows by
independent block sign reversal, and the top/bottom extremizers supply
equality by splitting their agreements and disagreements. On any fixed
cut the two directed parts sum to its full bilinear form, each bounded
by D. This proves the asserted lower bound and shows that the target
factor in (2) is not contradicted by the inherited seed witness alone.

## 2. Exact finite failures at actual width minimizers

The checker first exhausts every switching class of A (through n=7) to
select an actual global width minimizer. For this fixed A it then imposes
ALL directed-cut inequalities on the orientation variables and minimizes
the integer D with CP-SAT. It also evaluates the returned lift on every
Boolean spin in dimension 2n. Results for the selected actual minimizers:

    n       W(A)      min D       W(A)/sqrt(2)
    5         4          3           2.828427...
    6         5          5           3.535534...
    7         8          6           5.656854...

All three integer optimizations closed with status OPTIMAL. For n=6 the
separate direct enumeration of all 32768 orientations independently gives
the same minimum D=5. The n=6 seed has upper-triangle signs, in lexicographic
edge order,

    (1,1,1,1,1,-1,-1,1,1,1,-1,1,1,-1,-1).

It is the symmetric-conference width/cap minimizer. Every lift in this
orientation family has cap at least 20 BEFORE filling its six missing edges.

These disprove an exact finite factor-1/sqrt(2) theorem even on actual
width minimizers. They do not disprove (2): an o(n^(3/2)) error can absorb
all fixed-order failures. No tensor repetition of these finite witnesses
is assumed to preserve their actual orientation optimum.

## 3. The existing spectral obstruction applies only to certificates

For every choice of C,

    ||L||_F^2=4n(n-1),
    ||L||op>=sqrt(2(n-1)).

Thus the operator-norm cap certificate n||L||op, divided by (2n)^(3/2),
is at least (1/2)sqrt(1-1/n). Since the actual width minimizer has a
strict-subhalf all-order upper available, such a spectral certificate
cannot prove its scale-preserving (2). This is a certificate floor,
NOT an actual Boolean lower bound of 1/2 for the skew-lift family.

The archive's related skew-Clifford calculation in
joint_finite_fibre_action_audit.md section 5 likewise proves an operator
certificate floor; it does not establish a Boolean floor for this precise
A/−A skew-cross lift. No stronger scalable obstruction was found in this
bounded audit.

Reproduction:
computations/decisive_bridge_skew_lift_orientation_2026_09_07.py
--n 6 --seconds 30 --brute

The executable stores no unverified random-orientation assumption. It
checks all partitions and sign profiles exactly, including both directions
of each cut, and fixes two redundant shore signs only because absolute
value removes their common sign action.

## 4. Cross-order inequalities which do follow without the missing theorem

For skew C, its full rectangular Boolean norm satisfies

    beta(C)=2 max_disjoint(z,w)|z^TCw|.

Indeed write any full Boolean x,y as x=z+w,y=z-w, so x^TCy=-2z^TCw;
the converse parametrization is exact. Thus (1a) gives, for EVERY A,C,

    Q(L)<=2W(A)+beta(C).

After filling the matching, this proves the genuine cross-order inequality

    M_(2n)<=2W_n+min_(C skew full sign) beta(C)+n.         (3)

At an order admitting a skew conference matrix, beta(C)<=n sqrt(n-1),
so (3) becomes M_(2n)<=2W_n+n sqrt(n-1)+n. This statement is conditional
on that order; no unverified all-order skew-conference existence is used.

Without any existence input, independent random skew signs give an
explicit weaker all-order comparison. For any fixed x,y, x^TCy has
independent Rademacher coefficients x_i y_j-x_j y_i in {0,+2,-2}.
Writing tau_i=x_i y_i, precisely |{tau=+1}| |{tau=-1}| coefficients
are nonzero. Its variance proxy is at most n^2. There are at most
2^(2n+1) signed tests including the absolute-value sign. The elementary
Rademacher log-mgf bound and optimized soft maximum therefore give

    E beta(C)<=n sqrt(2(2n+1)log2).

Consequently

    M_(2n)<=2W_n+n sqrt(2(2n+1)log2)+n.                  (4)

Both inequalities are valid but noncompetitive with the known direct
all-order cap upper bound. Neither supplies (2) or original convergence.

The newly audited bipartite lower bound gives beta(C)>=2c_* n^(3/2)
asymptotically even when C is skew. Through (1a) this only yields
D>=c_* n^(3/2)/2. Applying the original lower bound directly to the
completed lift is stronger: D>=c_* n^(3/2)/sqrt(2)-o(n^(3/2)). Neither
inequality disproves the requested factor relative to W(A).
