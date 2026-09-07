# Exact actual-sign map and the remaining positive obligation

2026-09-07. Status: exact identities plus exploratory alternatives. No
original convergence theorem is claimed.

## 1. The constrained compiler must preserve whole orthogonal row codes

Let each H_i be k selected rows of a Hadamard matrix of order m. For every
row spin x_i in {+-1}^k, write

```math
u_i=H_i^T x_i/\sqrt k,\qquad \|u_i\|_2^2=m,
\qquad \mathcal C_i=\{H_i^Tx_i/\sqrt k:x_i\in\{\pm1\}^k\}.
```

For a symmetric signing S including freely chosen diagonal entries, define
the actual full symmetric sign matrix of order N=mk by

```math
W_{(i,a),(j,b)}=S_{ij}H_i(a,j)H_j(b,i).
```

The diagonal blocks in this formula are rank-one full-sign blocks. For
both sigma=+-1, the exact all-incidence defect is

```math
D_\sigma(u)=\sum_{i,j}(u_i(j)-\sigma S_{ij}u_j(i))^2
=2\left(m^2-\frac{\sigma x^TWx}{k}\right).
\tag{1}
```

Thus a favorable actual cap is a lower bound on the distance between the
PRODUCT OF THE WHOLE ROW CODES C_i and the two signed incidence-matching
subspaces. Deleting W's diagonal costs at most N/2 in the quadratic energy.
The global child-reversal issue is retained by requiring both polarities.

Replacing C_i by all words of its empirical type enlarges the spin set.
The hard compiler proves that this enlarged set contains many nearly
matched arrays. It cannot be turned into a favorable parent upper bound:
the inequality has the wrong direction. Actual recovery needs a compiler
constrained to the precise orthogonal-cube code, or a comparison proving
that its controlled enlargement preserves the relevant distance.

There is also an exact incidence-gauge obstruction. For arbitrary
`g_i(j) in {+-1}`, replacing

```math
H_i\mapsto H_i\operatorname{diag}(g_i),\qquad
S_{ij}\mapsto S_{ij}g_i(j)g_j(i)
```

leaves W unchanged. In particular a factorization `S_ij=z_i(j)z_j(i)`
followed by matching column rephasings can cancel S completely. A balanced
factorization does not by itself preserve an optimizing child's energy.

## 2. Controlled pair routing has an exact spectral calculus

Pair the m incidence coordinates within each row, and retain a Hadamard
row subspace antisymmetric on each pair; for m even it has dimension m/2
and a flat Hadamard basis obtained from `H_(m/2) tensor [1,-1]`, followed
by the required pair permutations and phases. The signed edge flip and
the within-row pair involution form alternating circuits.

On the antisymmetric subspace, a circuit of length ell compresses the
edge flip to one half of a signed cycle adjacency. Its eigenvalues are
`cos((2pi r+theta)/ell)`, where theta is 0 or pi according to the cycle
holonomy. For an even four-cycle with negative holonomy, the spectral
radius is 1/sqrt(2). With retention p=1/2, the normalized spectral cap is

```math
\frac{1}{2\sqrt p}\frac1{\sqrt2}=\frac12.
```

This is an actual routing operation, but it only recovers the spectral
half scale. Odd cycles always have an absolute eigenvalue one; longer
frustrated even cycles have larger spectral radius. A strict-subhalf
extension must retain Boolean non-saturation of these eigenspaces, not
only their eigenvalues. No minimizer-preserving routing law was found in
this bounded exploration.

## 3. An independent analytic alternative for signed Euler counting

For binary hard matching, choose one independent edge sign z_e; exact
balance is a signed-incidence lattice probability. Its Fourier product is
`product_(i<j) cos(theta_i+S_ij theta_j)`. Put
`c_i=cos(2theta_i)`, `s_i=sin(2theta_i)`, and `cbar=m^-1 sum c_i`. Then
the following identity is exact:

```math
\sum_{i<j}\sin^2(\theta_i+S_{ij}\theta_j)
=\frac14\left[(m-1)\sum_i s_i^2
+m\sum_i(c_i-\bar c)^2+s^TSs\right].
```

Therefore `||S||op=o(m)` forces global localization near the two saddle
classes in which every theta is near 0 or every theta is near pi/2,
modulo pi. The two contributions have relative sign equal to the parity
of the number of positive edges, matching the exact incidence parity
obstruction. Bounded cap implies the required operator estimate through
the audited fourth-moment theorem.

A full signed local central limit theorem was NOT completed. The root's
transition-system count is stronger and avoids this analytic tail problem.
Adjacent primary sources, consulted as possible methods rather than
imported signed theorems, were:

- Barvinok--Hartigan, *Maximum entropy Edgeworth estimates of the number of
  integer points in polytopes*, https://arxiv.org/abs/0910.2497.
- Barvinok--Hartigan, *The number of graphs and a random graph with a given
  degree sequence*, https://arxiv.org/abs/1003.0356.
- Isaev, *Asymptotic behavior of the number of Eulerian orientations of
  graphs*, https://arxiv.org/abs/1110.2598.

## 4. Comparison with the archive

The raw independent-edge type-repair method was already proved in
`decisive_bridge_positive_homogeneous_pressure_2026_09_06.md` and
`transfer_reconstruction_homogeneous_edge_count_2026_09_06.md`; no novelty
is claimed for that step. New candidate contributions of this branch are:

- an arbitrary-heterogeneous canonical dual and explicit uniform edge-atom
  lower bound;
- fixed finite-type seed-pressure Lipschitz comparison in beta(S-T), via a
  smooth separated-kernel lemma;
- independent audit of the root's entropy-rich signed Euler compiler; and
- a multicolor hard compiler using load-balanced alternating-path degree
  repair, with exact finite regression.

The first independent freeze's broad all-order generator and lower-tail
variational targets remain open architectures. They have not acquired a
recovery theorem, and are not presented as strict reductions of the original
problem. The positive frontier is constrained row-code realization in (1),
not another enlargement to local-type pressure.
