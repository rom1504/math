# Full-forest root-hit contraction: independent director proof

Date: 2026-09-06. This is a proof of one deterministic/random-kernel module,
not by itself a full nonlinear-response theorem. All forests and degrees
below are fixed before the order tends to infinity.

Let A be a symmetric hollow signing, B=A/sqrt(n-1), Q=B^2, with fixed
||B||op<=L. Then beta(A)<=L n sqrt(n-1). Coefficients suppressed below
depend only on the finite forests and their Gaussian-chaos normalizations.

## 1. The contraction matrix

An unmarked response forest R has an external root j, k>=3 odd branches
attached there, and D odd marked nonroot vertices. Every marked vertex
has odd graph degree. Its fully injective coefficient kernel has scale
(n-1)^(-D/2). Let P be an even child forest rooted at a, with d even marked
nonroot vertices and scale (n-1)^(-d/2). Its marked vertices have odd
degree, and its root has even degree. Different branches may be of
different old-tree types.

In the root-hit term, choose one marked slot of R and set its label to a.
Contract its other D-1 slots into distinct marked slots of P. If d<D-1,
the term is zero. Otherwise E=d-D+1 even marked slots of P remain. The
result J_(a,j) is a fixed-degree homogeneous multilinear Gaussian chaos
in these E slots, with deterministic internal sums over the D-1 contracted
labels. Its normalization is (n-1)^(-(D+d)/2). There are only finitely many
choices of slot matches; prove the estimate for one and sum them.

We claim

    E ||J||_F^2 = O_L(sqrt(n)).                         (1)

The same bound holds with Rademacher inputs because all remaining slots
are distinct within each summand. The homogeneous Gaussian-chaos convention
only changes fixed factorial constants.

## 2. Leading free-label count in the squared Frobenius norm

First exclude j from the remaining marked slots of P. It is already excluded
from the D-1 contracted slots because those belong to the injective R.
Within one coefficient graph, the two roots a,j, the contracted labels, and
the remaining marked labels are now all distinct.

When two copies of J_(a,j) are multiplied, expectation pairs their E
remaining marked labels by a bijection. At the leading label count there
are no further identifications: the D-1 internal contracted labels in the
first copy and those in the second copy are distinct from each other and
from the paired remaining labels. Including the summed roots a,j, the
number of distinct labels is

    V=2+2(D-1)+E=D+d+1.

The squared normalization has exponent D+d. Therefore absolute counting
alone gives O(n), not a sufficient estimate. Any additional identification
loses at least one label and contributes O(1). There are only finitely many
such equality patterns, with fixed degrees. This includes the previously
excluded case in which j is one of the remaining marked labels: its marked
partner is forced to be j in the other copy, so at least one free label is
lost. Terms with a=j vanish under injectivity of R.

## 3. A parity edge necessarily survives at the leading count

Delete paired edges modulo two in the doubled coefficient graph. At root j,
the first copy of R has k distinct incident edges. Their other endpoints
are all contracted labels, except possibly a. They cannot be remaining
marked slots of P: every R slot other than the designated a was contracted.

The second copy has its own distinct set of D-1 internal contracted labels.
Thus the only possible common edge incident to j between the two copies
is ja. At least 2(k-1)>=4 incident edge occurrences survive modulo two.
No edge of either P can cancel these edges, because P has no vertex with
label j in the leading pattern. Hence the reduced simple parity graph is
nonempty.

All labels, including a and j, are summed in the Frobenius norm. Choose any
surviving parity edge uv, fix the other labels, and factor the remaining
sign product into A_uv f(u)g(v), with |f|,|g|<=1. The distinctness exclusions
from already fixed labels are unary zeros. The remaining exclusion u!=v
is enforced by A_uu=0. Consequently the sum over u,v is at most beta(A)
in absolute value. Equivalently, after normalized averaging, the factor is
at most beta(A)/n^2=O_L(n^(-1/2)).

Multiplying by the O(n) nominal leading scale proves O_L(sqrt(n)) for
every leading pattern. The lower label patterns cost O(1), proving (1).
There is no assumption of generic cancellation or Gaussian independence.

## 4. The root-hit term vanishes at the required averaged scale

Let the transported forest be Z_i=sum_j B_ij R_j. In a full contraction
of its kernel into Y_i=sum_a Q_ia N_a P_a, the class hitting N_a once is,
up to a fixed factor,

    sum_(a,j) Q_ia J_(a,j) B_ji = [Q J B]_ii.

It follows from (1) that

    (1/n) E sum_i |[Q J B]_ii|^2
      <= (||Q||op^2 ||B||op^2/n) E||J||_F^2
      =O_L(n^(-1/2))=o(1).                         (2)

This is sufficient for averaged mixed-contraction vanishing. It is neither
a uniform-root estimate nor a claim of small operator norm for J. In
particular, no new matrix-derivative theorem is needed for this root-hit
class. The root-not-hit contractions, local Gaussianization, nonlinear
covariance identification, and bounded feasible response limits remain
separate modules in the full theorem.
