# Leading discrepancy from the entire switched-permuted convex orbit hull

Status: proved elementary obstruction, independently reconstructed. It excludes
small **difference norm**, not favorable cap comparison or selected-child
flatification. The row-l1 strengthening below was suggested by the constructive
agent after the adversarial agent noted strict-convexity failure of exact
orbit averaging.

For a hollow symmetric real matrix E define

```math
Q(E)=\max_{x\in\{\pm1\}^N}|x^TE x|/2,
\qquad \beta(E)=\max_{x,y\in\{\pm1\}^N}|x^TEy|.
```

## 1. General row-deficit theorem

Suppose D is hollow symmetric and every row satisfies

```math
\sum_{j\ne i}|D_{ij}|\le(1-\varepsilon)(N-1).
```

Then for every hollow full signing C,

```math
Q(C-D)\ge {\varepsilon\over4\sqrt2}N\sqrt{N-1}.       (1)
```

Proof: E=C-D has row l1 at least epsilon(N-1), hence by Cauchy--Schwarz
its row l2 is at least epsilon sqrt(N-1). For uniform Rademacher x,

```math
\beta(E)\ge\mathbb E\|Ex\|_1
\ge {1\over\sqrt2}\sum_i\|E_{i,*}\|_2
\ge {\varepsilon\over\sqrt2}N\sqrt{N-1},
```

using the sharp real Khintchine inequality. To check the polarization factor,
write z=(x+y)/2,w=(x-y)/2. They belong to the continuous cube and symmetry
gives x^TEy=z^TEz-w^TEw. Hollowness and independent Boolean rounding imply
|v^TEv|<=2Q(E) for every cube vector v. Therefore beta(E)<=4Q(E).

A completely elementary weaker constant, independent of sharp Khintchine, is
epsilon/(4sqrt(3)). Indeed for Z=sum a_j x_j, E Z^4<=3(E Z^2)^2; Hölder gives
E Z^2<= (E|Z|)^(2/3)(E Z^4)^(1/3), so E|Z|>=||a||_2/sqrt(3).
Thus even that self-contained bound proves a positive leading-scale defect.

## 2. Application to every convex orbit average of two children

Let N=m+l and

```math
B=\operatorname{diag}\left(
\sqrt{(N-1)/(m-1)}A_m,\sqrt{(N-1)/(l-1)}A_l\right),
```

where the children are arbitrary full signings, including exact minimizers.
The row-l1 norms are exactly sqrt((N-1)(m-1)) or sqrt((N-1)(l-1)).
Let D be any convex combination of matrices UBU^T, where U is a signed
permutation matrix; optionally allow global edge negation in the orbit too.
The triangle inequality preserves the common row bound

```math
\sum_j|D_{ij}|\le(N-1)\sqrt{(\max(m,l)-1)/(N-1)}.
```

Consequently (1) applies with

```math
\varepsilon_N=1-\sqrt{(\max(m,l)-1)/(N-1)}.
```

For comparable splits m,l in [N/4,3N/4], epsilon_N>=1-sqrt(3)/2>0.
For equal children epsilon_N tends to 1-1/sqrt(2). In particular the equal
split asymptotic discrepancy constant is at least
`(sqrt(2)-1)/8`, approximately .0517767.

The optimization is over the ENTIRE convex switched-permuted orbit hull, and
the conclusion is uniform over the full signing C. Thus even mixing many
different partitions does not make a subsequent small-norm correction possible.
No special property of the child sign patterns enters the proof.

## 3. Exact scope

Every orbit member has cap Q(B), so convexity gives Q(D)<=Q(B). The theorem
proves that completing a construction by the triangle inequality
Q(C)<=Q(D)+Q(C-D) necessarily budgets a leading error if its only available
estimate for the remainder is its norm. It does **not** show
Q(C)>=Q(D)+Omega(N^(3/2)): norm cancellation can be substantial.

The selected-child target explicitly permits joint reorganization of internal
and bridge responses, so it survives. This is a stronger version of the
already recorded fixed-bridge difference-norm obstruction, extending its
scope to arbitrary convex orbit mixing; it is not an original-value lower
bound or a nonconvergence mechanism.

## 4. Independent finite numerical diagnostic

`tmp/flatify_adversary_2026_09_07_orbit_hull_lp.py` enumerates all eight
first-row-gauged order-four parents. The two-order-two-child orbit consists
of the twelve signed perfect matchings, each scaled by sqrt(3). For each
parent a linear program minimizes Q(C-D) over their full convex hull.
All parent spins are included as constraints, and the returned mixture is
independently reevaluated. Its JSON preserves every mixture and discrepancy.

For the six optimal parent signings, the numerical minimum is
1.26794919243112 (consistent with 3-sqrt(3)); for the two cap-six parents
it is twice this value. These LP optima are numerical checks, not additional
exact theorems. All satisfy the independently proved bound (1).
