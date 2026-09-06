# Gap-two orientation dichotomy and response-code geometry

Date: 2026-09-06. Exact consequences for a hollow symmetric signing that
is locally minimizing under every single coefficient flip. Exact global
minimizers satisfy this hypothesis. The cap is `Q=max_x |H_A(x)|`, with
`H_A(x)=sum_(i<j) a_ij x_i x_j`.

The first section independently audits the director's strengthened
gap-two separation argument. The later sections give an orientation
dichotomy and its general antipodal-code counterpart. These statements
do not supply the common near-state law or bridge discrepancy needed for
original convergence.

## 1. The square-root separation threshold

Switch a chosen oriented ground state to the positive all-one state of
a signing `B`. Thus `H_B(1)=Q`, and its local fields

\[
\ell_i=\sum_{j\ne i}B_{ij}\ge0,
\qquad \sum_i\ell_i=2Q.
\]

For an integer `r<=n/2`, assume

\[
n(2r-1)<Q,\qquad Q>r(n-r)+1.                       \tag{1}
\]

There is an oriented state of gap at most two whose projective Hamming
distance from the chosen ground exceeds `r`.

Suppose not. Represent every such state by a flipped set `S`, `|S|<=r`.
A negative-oriented near-state would satisfy
`Q<=cut_B(S)+1<=r(n-r)+1`, contrary to (1). All near-states therefore have
positive orientation.

At any positive near-state of gap `g<=2`, the local field measured in its
own spin coordinates is at least `-g/2>=-1`: otherwise a single spin flip
would exceed the cap `Q`. For `i in S` this local field equals

\[
\ell_i^S=-\ell_i+2\sum_{j\in S\setminus\{i\}}B_{ij}.
\]

Consequently every flipped vertex satisfies

\[
\ell_i\le 2(|S|-1)+1\le2r-1.                      \tag{2}
\]

The set `T={i:ell_i>2r-1}` is fixed at `+1` in every gap-two state. It
must contain a positive edge: indeed,

\[
\sum_{i\in T}\ell_i\ge2Q-n(2r-1)>Q,
\]

whereas if all its internal edges were negative then
`sum_T ell_i <= cut_B(T) <= Q`. The final inequality follows by comparing
`H_B(1^T)=Q-2cut_B(T)` with the lower cap `-Q`.

Flipping a positive edge inside `T` decreases every gap-two oriented
energy by two. Every other oriented state has gap at least four, by
parity, so even its possible increase by two leaves it below `Q`. This
contradicts coefficient-flip local minimality.

Thus for any fixed universal lower coefficient `c>0`, along exact
minimizers with `Q>=(c-o(1))n^(3/2)`, every oriented ground has a gap-two
state at projective distance at least

\[
\left(\frac c2-o(1)\right)\sqrt n.                 \tag{3}
\]

This improves the earlier quadratic internal-energy threshold, which
only gave `n^(1/4)` separation. The gain uses local fields at the
near-state itself, not just at the ground.

## 2. Exact orientation dichotomy

Let `mathcal N` be the entire oriented gap-two window, identifying spins
under global negation. Assume `Q>2`, so a spin cannot occur there in both
orientations. Retain the switched positive all-one ground.

Either:

1. There is a **negative-oriented** state in `mathcal N`, at projective
   distance `s` satisfying

   \[
   s(n-s)\ge Q-1,
   \qquad s\ge\frac{Q-1}{n};                       \tag{4}
   \]

or:

2. Every state in `mathcal N` has positive orientation. If its number is
   `k`, then

   \[
   2^{k-1}\ge\frac{n^2}{2Q+n},
   \qquad
   k\ge1+\log_2\frac{n^2}{2Q+n}.                  \tag{5}
   \]

   Moreover, after selecting any representatives of the near-states
   relative to the all-one ground, the union of the coordinates they
   change has size at least `(n-1)/2`.

For (4), write the negative near-state as `1^S`, `|S|=s<=n/2`. Its
energy is `H_B(1^S)=-Q+g`, `0<=g<=2`, so
`cut_B(S)=Q-g/2>=Q-1`. The cut has at most `s(n-s)` edges.

For the second alternative, choose representatives
`x^(1)=1,x^(2),...,x^(k)` and assign vertex `i` its literal response code

\[
c_i=(x_i^{(1)},\ldots,x_i^{(k)})\in\{1\}\times\{\pm1\}^{k-1}.
\]

Vertices with the same code form a class `C`. Every edge inside such a
class is negative. Otherwise a positive internal edge is aligned in
every gap-two state, and flipping it produces the contradiction used in
Section 1.

Independently choose one random sign for each class and assign it to all
vertices of that class. Every interclass term averages to zero and all
intraclass terms remain negative. Therefore

\[
Q\ge\sum_C\binom{|C|}{2},
\qquad \sum_C|C|^2\le2Q+n.                        \tag{6}
\]

Cauchy--Schwarz gives at least `n^2/(2Q+n)` nonempty classes. There are
at most `2^(k-1)` literal codes, proving (5). No antipodal identification
of vertex codes is needed: each starts with the anchored coordinate
`+1`, and independent class switches are legitimate Boolean spins.

The unchanged coordinates are one such negative-clique class `C_0`.
For any of its vertices, ground local-field positivity gives

\[
0\le\ell_i\le-(|C_0|-1)+(n-|C_0|)
=n-2|C_0|+1.
\]

Thus `|C_0|<=(n+1)/2`, proving the linear collective-support bound.

Along exact minimizers, the all-order spectral construction gives
`Q<=(1/2+o(1))n^(3/2)`. The single-orientation alternative consequently
forces

\[
k\ge\left(\frac12-o(1)\right)\log_2n.              \tag{7}
\]

This is in the fixed **gap-two** window. It is a conditional improvement
over the general logarithm-over-logarithm count in an `O(sqrt(n))` window;
it does not replace that unconditional theorem when both orientations
occur.

## 3. General orientation vector and antipodal code pairs

The common-code argument has an exact version without a common
orientation. Enumerate all gap-two states as `(sigma_a,x^(a))`, with the
first state the switched ground, and put

\[
\sigma=(\sigma_1,\ldots,\sigma_k),\qquad
c_i=(x_i^{(1)},\ldots,x_i^{(k)}),\qquad
\tau(c)=\sigma\odot c.
\]

The first coordinate of every code and of `sigma` is `+1`. If two
distinct vertices have codes `c` and `tau(c)`, then

\[
\sigma_a x_i^{(a)}x_j^{(a)}=1
\quad\text{for every }a.
\]

Their edge must therefore be negative by single-edge local minimality.
When orientations are mixed, `tau` is a fixed-point-free involution on
the anchored code space. It pairs classes, and every edge between paired
classes is negative. Empty classes are allowed.

Let `m_c=|C_c|`. Then

\[
\sum_{\{c,\tau(c)\}}m_cm_{\tau(c)}\le Q.           \tag{8}
\]

To prove it, choose an independent random sign per unordered `tau` pair,
use that sign on one class, and use either the same or the opposite sign
on its partner. Averaging kills every interaction between different
pairs. The sum `D` of all within-class edge energies is fixed. Choosing
all partners the same gives expectation `D-L`; choosing all opposite
gives `D+L`, where `L` is the left side of (8). Both expectations lie in
`[-Q,Q]`, so `|D|+L<=Q`.

Equivalently the empirical code law `p_c=m_c/n` obeys

\[
\sum_c p_c p_{\tau(c)}\le\frac{2Q}{n^2}.           \tag{9}
\]

At the `n^(3/2)` cap scale this overlap is `O(n^(-1/2))`. With all
orientations positive, `tau` is the identity and the diagonal correction
recovers (6), not (8).

The distinction matters: for a nontrivial involution, a probability law
can be supported on a transversal of its pairs and have zero overlap
even with finitely many atoms. Therefore (9) alone does not force many
response codes in the mixed-orientation case. Replacing it by the
single-orientation collision bound would be incorrect.

## 4. Scope of the macroscopic-response conclusion

The collective support in Alternative 2 is linear, but it need not be
contained in any one distant state. It can be distributed over many
near-states. The opposite-orientation alternative only guarantees
square-root separation, and (3) has the same scale without orientation
conditions.

None of these facts yields a linear Hamming separation, an extensive
packing, a covariance lower bound, or unavoidable discrepancy against
every new coefficient row. In particular they do not yet produce the
critical-window response theorem needed for an upper-preserving
cross-order insertion. They are exact constraints on actual minimizers,
not on arbitrary asymptotic near-minimizers that may tolerate a
two-unit cap improvement.
