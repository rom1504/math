# Affine-MUB block conjugation: actual signs and an exact spike-sector obstruction

Date: 2026-09-07. Status: **Pending independent audit**. The construction and
normalizations below are exact. The asymptotic lower uses the stated standard
bounded-rank hypergraph edge-coloring theorem. This is not a favorable
flatification theorem, a bound on the original liminf, or nonconvergence.

## 1. Concrete proposed construction

Let A be a hollow symmetric signing of order k. Let q admit a real Hadamard
matrix and a net consisting of k parallel classes on q² points: each class
partitions the points into q lines of size q, and lines in distinct classes
meet in exactly one point. A finite affine plane of order q supplies this
when q is a power of two and k<=q+1.

For each class i and each of its lines L, choose ANY q-by-q Hadamard H_(i,L),
with its columns indexed by points of L. Define the q²-by-q² orthogonal
matrix U_i whose row (L,a) has entry H_(i,L)(a,p)/sqrt(q) at p in L and zero
elsewhere. Orthogonality follows within each line from the Hadamard identity,
and between lines from disjoint support. Moreover

```math
 |(U_iU_j^T)_{ab}|=1/q\quad(i\ne j).
```

Thus the order-N=kq² matrix with cross-fibre blocks

```math
 C_{ij}=q a_{ij}U_iU_j^T\quad(i\ne j),\qquad C_{ii}=0
```

has actual full sign cross blocks. Arbitrary hollow full sign matrices D_i
may be inserted in its diagonal fibre blocks to make a genuine full signing.
Before those insertions, C is orthogonally similar to q(A tensor Id_(q²)).
This preserves seed spectral information without merely permuting the cube.
The diagonal-block triangle-inequality payment is O(k q³) when low-cap D_i
are selected, hence O(k^(-1/2)) after normalization by N^(3/2). This observation
alone does not control the cross-fibre cap.

For Boolean fibre vectors x_i and y_i=U_i^T x_i, the exact cross energy is

```math
 H_C(x)=q\sum_{p}H_A(y_1(p),\ldots,y_k(p)).                 (1)
```

The tempting lossless bound q³Q(A) is FALSE in general. Its failure cannot
be repaired just by selecting arbitrary diagonal-fibre full signs.

## 2. The finite partial-support linear program

For sigma in {+1,-1}, define

```math
 \kappa_\sigma(A)=
 \max\left\{\sum_z\lambda_z\sigma H_A(z):
       \lambda_z\ge0,\quad
       \sum_z\lambda_z z_i^2\le1\ (1\le i\le k),\quad
       z\in\{0,1,-1\}^k\right\}.                         (2)
```

Zero and singleton vectors have zero objective and can be omitted. Equivalently,
use one variable per nonempty support S and objective
P_sigma(A_S)=max_(z in signs on S) sigma H_(A_S)(z). Only patterns with positive
objective are useful. The LP dual is

```math
 \kappa_\sigma(A)=\min\left\{\sum_i d_i:
      d_i\ge0,\quad
      \sigma H_A(z)\le\sum_i d_i z_i^2
           \text{ for all ternary }z\right\}.             (3)
```

This is a finite polyhedral diagonal-majorant relaxation, NOT a PSD
majorant: requiring the same inequality for every real vector is stronger.
It lies between the one-polarity Boolean cap and the corresponding
correlation-matrix SDP value. Indeed sum_z lambda_z zz^T is PSD with diagonal
at most one; add a nonnegative diagonal to obtain a correlation matrix,
which does not change its pairing with hollow A.

## 3. Exact upper bound on the pure-column spike sector

On every physical line block choose x_(i,L)=epsilon_(i,L) H_(i,L)(:,p_(i,L))
for one point p_(i,L) in L. Then y_i equals epsilon sqrt(q) at the selected
point of each line and vanishes at every other point. Write y_i(p)=sqrt(q)z_i(p).
The vectors z(p) are ternary and satisfy

```math
 \sum_p z_i(p)^2=q,\qquad
 H_C(x)=q^2\sum_p H_A(z(p)).                              (4)
```

Taking lambda_z equal to its point frequency divided by q proves, EXACTLY,

```math
 \max_{x\text{ in the spike sector}}\sigma H_C(x)
       \le q^3\kappa_\sigma(A).                           (5)
```

The following lower proves that (5) is asymptotically sharp in this sector
when the diagonal fibre blocks are zero. It also proves the SAME lower
against every arbitrary full diagonal-fibre completion.

## 4. General asymptotic lower, robust to arbitrary fibre completion

**Theorem.** Fix k and A. Along any sequence of the above q-nets with q->infinity,
and uniformly over the choices of line Hadamards and arbitrary hollow full
diagonal-fibre completions D_i, the assembled full signing F satisfies

```math
 \max_x\sigma H_F(x)\ge
       [\kappa_\sigma(A)-o(1)]q^3.                       (6)
```

For zero diagonal-fibre blocks, combining (5) and (6) identifies the exact
asymptotic spike-sector cap. There is NO matching upper asserted for the
whole Boolean cube: non-spike states can do better.

### 4.1 Hypergraph and its checked hypotheses

Fix an optimal LP solution with one signed pattern z^S per support S of size
at least two, and write its weights lambda_S. Choose fixed L>0 with
L>=max_S lambda_S. For every point p and useful S, independently retain the
hyperedge

```math
 e(p,S)=\{(i,L_i(p)):i\in S\}
```

with probability lambda_S/L. Hypergraph vertices are the kq physical lines.
Distinct points give distinct edges for a fixed S because two lines from
different classes determine their intersection. Edge sizes are between 2 and
k. Every vertex (i,L_i) has expected degree (q/L)sum_(S contains i)lambda_S,
at most q/L. Standard Bernoulli concentration and a union bound over kq
vertices give maximum degree at most q/L+o(q) with probability tending to one.
At least one constraint is tight for a nonzero optimum, so maximum degree is
also q/L-o(q). Each pair of vertices has codegree at most 2^(k-2), because
its two lines determine at most one point and there are only that many
support choices. These constants depend on fixed k, never on q or the D_i.

The bounded-rank small-codegree edge-coloring theorem therefore partitions
the retained hyperedges into at most (1+o(1))q/L matchings. A directly checked
primary statement is Jeff Kahn, *Asymptotically good list-colorings*, JCTA 73
(1996), 1–59, DOI 10.1006/jcta.1996.0001: bounded edge size, maximum degree D,
and maximum pair degree o(D) imply list chromatic index (1+o(1))D. Ordinary
edge coloring is the weaker conclusion used here. The author's university
publication record spells out all hypotheses and the limit:
[Rutgers primary record](https://www.researchwithrutgers.org/en/publications/asymptotically-good-list-colorings/).

### 4.2 Boolean law associated to one matching

For every matched hyperedge e(p,S), choose its line blocks to be the Hadamard
columns at p, with signs epsilon_e z_i^S for i in S. The epsilon_e are
independent fair signs. Every unmatched physical line independently chooses
a uniform Hadamard column and a fair overall sign.

All cross-energy interactions between different components have mean zero.
The contribution inside matched component e(p,S) is exactly
q² H_A(z^S). No component contains two physical lines from the same fibre,
so all D_i interactions between different physical lines also average to zero.
Only the internal energy of each selected physical line remains.

Write b_(i,L)(p)=H_(D_i restricted to physical line L)(H_(i,L)(:,p)). Since
the principal line block is hollow,

```math
 \sum_{p\in L}b_{i,L}(p)=0,\qquad
 |b_{i,L}(p)|\le q(q-1)/2.                                (7)
```

Unmatched lines have zero average internal energy. Thus the expected signed
full energy for a matching M is the sum over e(p,S) in M of

```math
 q^2\sigma H_A(z^S)+
       \sigma\sum_{i\in S}b_{i,L_i(p)}(p).                (8)
```

### 4.3 Average over the color classes

Sum (8) over all retained edges. Its expectation over the initial independent
edge sampling is

```math
 (q^4/L)\kappa_\sigma(A).                                 (9)
```

The b terms cancel exactly in expectation: each i,L has coefficient
(sum_(S contains i)lambda_S)/L multiplying the zero sum in (7). Each summand
in (8) is O_k(q²), and there are O_k(q²) independent inclusion variables.
Consequently the variance of the total is O_k(q^6), uniformly over all line
Hadamards and D_i. Chebyshev implies that with probability tending to one,
the total in (9) has error o(q^4), simultaneously with the degree bounds.
Fix one realization having both properties.

At least one of its at most (1+o(1))q/L color classes has expected full energy
at least [kappa_sigma(A)-o(1)]q³. Some Boolean assignment in that finite law
attains its expectation or more. This proves (6).

## 5. The actual optimal pentagon seed: kappa_+=kappa_-=5>Q=4

Take the pentagon Seidel signing A_5, with negative cycle edges and positive
chords (switching-equivalent to the archived exact order-five witness).
Five triangles have edge-sign product +1 and five have product -1. In each
family every vertex belongs to three triangles. A positive-product triangle
has positive cap 3; a negative-product triangle has negative cap 3.
Assign weight 1/3 to each of the five triangles of the chosen polarity.
This is feasible in (2) and has objective 5.

Conversely every ternary z obeys |H_(A5)(z)|<=sum_i z_i². Check by support size:
0,1 give zero; 2 gives at most1; 3 gives at most3; every four-vertex restriction
has cap4; and the full seed has cap4. Thus d_i=1 is feasible in BOTH duals,
and kappa_+(A5)=kappa_-(A5)=5 exactly. The value Q(A5)=4 is elementary: for
spins x, the cycle/chord signing's quadratic energies lie between -4 and4,
as direct enumeration of the sixteen projective spins verifies.

For this example random thinning is unnecessary. The hypergraph using the
five favorable triangle patterns at every point is 3-uniform, exactly
3q-regular, and has codegree at most3. Pippenger–Spencer directly gives an
edge coloring with (3+o(1))q colors. There are 5q² edges. Summing local line
energies over all edges cancels EXACTLY by (7), while the total signed cross
energy is 15q^4. Hence some color class and actual Boolean assignment attain
(5-o(1))q³ against ANY full diagonal-fibre completion.

The exact primary hypothesis statement is in Pippenger–Spencer,
*Asymptotic behavior of the chromatic index for hypergraphs*, JCTA51 (1989),
24–42, DOI 10.1016/0097-3165(89)90074-5:
[author's institutional record](https://scholarship.claremont.edu/hmc_fac_pub/1041/).

Since N=5q², this family is forced to have normalized cap at least 1/sqrt(5).
It cannot retain the finite seed's row-normalized coefficient
Q(A5)/(5sqrt(4))=2/5. The gap is specific to this infinite affine-MUB
construction class, with a fixed finite seed. It does not refute favorable
flatification of comparable selected large minimizing children.

## 6. Finite diagnostics and research scope

`computations/flatify_construct_2026_09_07_mub_diagnostic.py` builds the q=4
net using F4, checks every exact cross-block sign and basis identity, and
performs a reproducible Boolean local search. For k=5 it finds cross energy
324, above q³Q(A5)=256 and even above the pure-spike ceiling q³kappa(A5)=320.
This last fact illustrates why the spike functional is not a full-cube upper.
For k=4 the finite search reaches256 but finds no larger value; that failure
is not an upper certificate.

`computations/flatify_construct_2026_09_07_partial_packing.py` solves the
finite support LPs for archived exact children of orders3 through10. The
solver outputs are numerical unless separately rationally reconstructed.
The analytic pentagon equality above has no floating-point dependency.

The first attempted elementary H2 amplification was separately checked in
`computations/flatify_construct_2026_09_07_tensor_diagnostic.py`: its ratios
already exceed2sqrt(2) on several finite optimal witnesses, and can equal4.
Those finite failures do not themselves give an asymptotic obstruction.

No no-inflation condition on actual minimizing seeds is promoted as a new
proof target here. The intended positive construction has produced a concrete,
scalable obstruction and an exact diagnostic functional, not a recurrence.
