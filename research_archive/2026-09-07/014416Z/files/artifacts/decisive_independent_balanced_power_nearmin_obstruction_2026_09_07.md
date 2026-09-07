# Balanced pointwise powered aggregation fails on original near-minimizers

Status: complete proof, independently audited PASS by `transfer_seeds`
(reverse Minkowski, separate majorant cores, drift subsequence, all error
scales, and exact/near-minimizer quantifiers). Date: 2026-09-07.

This concerns a PRESCRIBED balanced principal partition of an asymptotically
minimizing parent. It does not concern an existentially selected partition,
an exact minimizing parent after modification, or minimum values alone.

## 1. Statement

There are even orders N_j tending to infinity, actual hollow signings A'_j,
and prescribed equal shores S_j,T_j, such that

```math
Q(A'_j)=M_(N_j)+o(N_j^(3/2)),
Q(A'_j)/N_j^(3/2) < 1/2-epsilon
```

for a fixed epsilon>0, but

```math
Q(A'_(j,S_j))^(2/3)+Q(A'_(j,T_j))^(2/3)-Q(A'_j)^(2/3)
 >= c N_j/log N_j.
```

Here c>0 is absolute. Thus any error O(N^(1-delta)), delta>0, fails for
pointwise powered aggregation on arbitrary original normalized near-minimizers,
even at exactly equal shores and with a strict-subhalf parent. The children
are not asserted strict-subhalf, nor exactly minimizing. The partition can
be chosen before the modification; no claim is made about all partitions of
the modified signing.

## 2. Two elementary ingredients

For r=2/3 the positive-cone reverse Minkowski inequality says

```math
[(a+d)^r+(b+d)^r]^(1/r)
 >= (a^r+b^r)^(1/r)+2^(1/r)d,       a,b,d>=0.          (1)
```

One direct proof puts u=a^r,v=b^r and observes that
g(u)=(u^(1/r)+d)^r is convex: for u>0,
g'(u)=(1+d/u^(1/r))^(r-1), which is increasing. Jensen gives
g(u)+g(v)>=2g((u+v)/2), equivalent to (1).

The second ingredient is the banked Grothendieck diagonal-majorant bound:
for any hollow symmetric signing B there is a nonnegative diagonal D with

```math
-D <= B <= D,                 tr D <= 4 K_G Q(B).     (2)
```

In particular, when |B|=N/2 and Q(B)<=C N^(3/2), at least N/4 indices
have D_ii<=16 K_G C sqrt N. Every k-subset I of those indices satisfies

```math
Q(B_I) <= 8 K_G C k sqrt N.                           (3)
```

The two shores may use DIFFERENT majorants. No unjustified intersection of
one global half-core with both shores is required.

## 3. Dyadic exact parents with an almost nonpositive value drift

Let m_N=M_N/N^(3/2), u_j=m_(2^j)^(2/3). This sequence is bounded above
and below by positive constants. Hence there is an infinite subsequence of
j for which

```math
e_j=max(0,u_j-u_(j-1)) -> 0.                          (4)
```

Otherwise the increments would eventually be bounded below by a positive
constant, contradicting boundedness. Set N=2^j and choose an EXACT minimizing
parent A_N. Fix ANY partition into two shores of size N/2. In N^(3/2) energy
units write

```math
P=Q(A_N)/N^(3/2),  a=Q(A_S)/N^(3/2),  b=Q(A_T)/N^(3/2).
```

Since both children dominate their order minimum,

```math
a^r+b^r >= m_(N/2)^r >= P^r-e_j.                    (5)
```

All of P,a,b have uniform positive lower bounds and finite upper bounds:
the lower bounds use the original universal lower cap at their respective
orders, and a,b<=P follows from principal-cap monotonicity.

## 4. Gauged clique additions

Define

```math
delta_j=max(1/log N, sqrt(e_j)),
k=floor(sqrt(2 delta_j) N^(3/4)),
d=k(k-1)/(2N^(3/2)).                                  (6)
```

Then delta_j->0, e_j<=delta_j^2, and d=delta_j+o(delta_j).
For large N, k<N/4. On EACH shore choose k indices in its core from (3).
Let x be a ground state of that shore, and let s in {+1,-1} be its energy
orientation, so s q_B(x)=Q(B). Replace the edges of the selected k-block by

```math
a'_il=s x_i x_l,             i<l in that block.        (7)
```

Every entry remains a sign. At the OLD child ground state, the new block
has signed energy exactly k(k-1)/2. By (3), removal of the old block costs
at most L k sqrt N, where L=8 K_G C can be fixed uniformly. Hence each new
child satisfies

```math
Q(A'_S)/N^(3/2) >= a+d-eta,
Q(A'_T)/N^(3/2) >= b+d-eta,
eta=L k/N = O(sqrt(delta_j) N^(-1/4)).                 (8)
```

The parent triangle inequality, with NO presumed cancellation between the
two gadgets and with arbitrary child orientations, gives

```math
Q(A'_N)/N^(3/2) <= P+2d+2eta.                        (9)
```

Because delta_j>=1/log N, eta=o(delta_j). Thus the new actual parent is an
original normalized near-minimizer and remains strictly subhalf by the banked
all-order bound limsup m_N<1/2.

## 5. Powered defect

All arguments of the power function stay in a fixed positive compact
interval. Errors eta and e_j therefore change the following powered
expressions by O(eta+e_j), uniformly. Apply (1) and (5):

```math
(a+d)^r+(b+d)^r
 >= [(a^r+b^r)^(1/r)+2 sqrt(2)d]^r
 >= [P+2 sqrt(2)d]^r-O(e_j).                         (10)
```

Equations (8)--(10) imply that the powered defect, divided by N, is at least

```math
[P+2 sqrt(2)d]^r-[P+2d]^r-O(e_j+eta)
 >= c_0 d-O(e_j+eta),                               (11)
```

where c_0>0 is uniform. For example, once P+2 sqrt(2)d<=1, the mean-value
theorem permits c_0=(2/3)(2 sqrt(2)-2). Since e_j<=delta_j^2,
eta=o(delta_j), and d/delta_j->1, (11) is at least c_0 delta_j/2 for all
sufficiently large j. This proves the claimed N/log N lower defect.

## 6. Exact boundary

The starting parent is exact; the MODIFIED parent need not be. Its excess
over M_N is O(delta_j N^(3/2)), and no assertion allows that excess to be
discarded on the much smaller O(N) energy scale corresponding to an
O(sqrt N) powered error. This theorem therefore does not refute the actual
minimum-value reverse-Fekete conjecture, a selectable exact-parent theorem,
or the existence of a better partition of A'_N. It is a balanced, original
near-minimizer pointwise obstruction, not an original convergence result.
