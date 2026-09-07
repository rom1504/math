# H2 instability persists with a bounded normalized operator norm

Status: complete actual-sign construction, pending independent audit.
This strengthens the subhalf part of the earlier H2 planting theorem.
It does NOT establish fixed-operator-norm near-minimizers for the ORIGINAL
problem; the restricted-class optimality scope below is essential.

Let `L(A)=[[A,A+I],[A+I,-A]]`. There are a fixed finite K, eta>0 and
actual hollow signings A_j of unbounded orders N_j such that

```math
Q(A_j)<\tfrac12N_j^{3/2},\qquad
\|A_j\|_{op}\le K\sqrt{N_j},\qquad
Q(L(A_j))\ge2\sqrt2 Q(A_j)+\eta N_j^{3/2}.              (1)
```

The proof uses the established all-order restricted-weave upper bound.
It does not assume the original sequence converges.

## 1. A fixed operator class closed under the canonical H2 lift

Fix any L greater than sqrt(32/31), for example L=2. At large orders let

```math
\mathcal C_N^L=\{A:\ A\text{ is a hollow signing},\
 \|A+I\|_{op}\le L\sqrt N-5\},\qquad
b_N=\min_{A\in\mathcal C_N^L}Q(A)/N^{3/2}.              (2)
```

These sets are nonempty for every sufficiently large N and satisfy
`b_N<=U<1/2` eventually. Here is the operator check on the banked weave
construction, so this assertion does not introduce an extra near-minimizer
regularization premise. The full restricted weave K has representation
`K=T^T S T`, where the block fibre transform obeys `T^T T=mI`, and S is
a signed coordinate-swap involution. Thus `||K||op<=m`. Its order is
`N=mk` with k/m tending to p=31/32. The hollow signing is A=K-diag(K),
so `||A+I||op<=m+2`. The all-order principal-restriction step changes
the ratio of orders by 1+o(1), retaining the bound
`||A+I||op/sqrt(N)<=1/sqrt(p)+o(1)`. A fixed L above this limit absorbs
the additive margin5. The cap upper bound is the independently audited
strict all-order certificate, so U can be chosen strictly below1/2.

The class (2) is EXACTLY closed under L(A):

```math
L(A)+I=(A+I)\otimes H_2+\operatorname{diag}(0,2I),
```

and hence

```math
\|L(A)+I\|_{op}
\le L\sqrt{2N}-5\sqrt2+2\le L\sqrt{2N}-5.
```

Consequently, for an exact restricted-class minimizer A,
`Q(L(A))>=b_(2N)(2N)^(3/2)`. Also `||A||op<=L sqrt(N)`.
The bounded dyadic sequence b_N has infinitely many increments
`b_(2N)-b_N` arbitrarily close to nonnegative.

## 2. A lift maximizer has two macroscopic supports

Fix a universal original lower constant ell>0 with
`2 sqrt(2) ell>2U`, for instance any certified ell between .43 and the
banked .4333221116640807. At all sufficiently large orders,
`M_(2N)>=(ell-o(1))(2N)^(3/2)`.

The lift has symmetric positive and negative extrema under
`(x,y)->(y,-x)`, so choose a positive extremizer (x,y). Set
`S={i:x_i=y_i}`, `T={i:x_i=-y_i}`. If h=min(|S|,|T|), the exact
cross-support identity and the operator norm give

```math
Q(L(A))\le2Q(A)+2Lh\sqrt N+N.                         (3)
```

For a small fixed kappa>0 depending only on L,U,ell, (3) and the universal
lower bound therefore force `|S|,|T|>=2 kappa N` eventually. For example
any kappa strictly less than
`(2 sqrt(2) ell-2U)/(8(L+1))` suffices. Select k=floor(kappa N) indices
from each part, without changing the maximizing spin.

## 3. Sparse edge flips add the coherent direction without a large eigenvalue

Put rho=sqrt2-1 and g0=(1+rho²)/2. On the selected 2k vertices, choose
the real hollow TARGET F which, after gauging by x, has entry+1 within
the selected S part, entry-1 within the selected T part, and entry rho
across them. It is zero outside this block. The target is only an
intermediate mean direction: the final output remains a flat signing.

As in the earlier coherent gadget calculation,

```math
Q(F)\le g_0k^2,
\quad q_{L_0(F)}(x,y)=2(1+\rho)k^2-2k
=R g_0k^2-2k,\qquad R=2(1+\sqrt2),                    (4)
```

where `L_0` omits the cross-diagonal completion and is linear in F.

Fix theta>0 and set epsilon=theta/sqrt(N), eventually less than1.
Independently flip each selected edge A_ij with probability
`epsilon(1-A_ij F_ij)/2`. Outside edges remain untouched. This produces
an actual signing A' with

```math
\mathbb E A'=A+\epsilon(F-A_{\rm selected}).            (5)
```

Let W=A'-E A'. Each selected upper-triangle entry is independent,
centered, bounded by2, and has variance at most4 epsilon. Scalar Bernstein
plus a union bound over the `2^(2k)` block spins gives, with probability
1-o(1),

```math
Q(W)\le20(\sqrt\epsilon\,k^{3/2}+k).                  (6)
```

For example use Bernstein's deviation parameter4k: its bound is at most
`8 sqrt(epsilon) k^(3/2)+(16/3)k`, and the union bound decays exponentially
in k. Separately each row changes on at most `4 epsilon k` edges with
probability1-o(1), by scalar Chernoff and union bound; its mean is at most
`2 epsilon k`. Therefore, on a common successful realization,

```math
\|A'-A\|_{op}\le\|A'-A\|_{\infty\to\infty}
\le8\epsilon k\le8\theta\kappa\sqrt N.                (7)
```

No matrix concentration theorem is required for this bounded-op conclusion.
Both estimates also hold for theta_N tending to zero provided
`theta_N>=c/log N`, as used below.

The original cap triangle inequality, (4)--(6), and
`Q(A_selected)<=L k sqrt(N)` imply

```math
Q(A')\le Q(A)+g_0\theta\kappa^2N^{3/2}
 +O_{L,\kappa,\theta}(N^{5/4}+N).                       (8)
```

The OLD lift maximizer is retained without changing a single spin. Its
mean gain is (4) minus an old-block contribution bounded by
`6 Q(A_selected)`. Moreover `Q(L_0(W))<=6Q(W)` by polarization. Thus

```math
Q(L(A'))\ge Q(L(A))+R g_0\theta\kappa^2N^{3/2}
 -O_{L,\kappa,\theta}(N^{5/4}+N).                       (9)
```

More uniformly the errors are
`O_(L,kappa)(sqrt(theta) N^(5/4)+(1+theta)N)`.

## 4. Dyadic extraction and precise near-minimality scope

Put delta=g0 theta kappa² and choose it small enough that U+delta<1/2.
For restricted-class minimizing A, subtracting `2 sqrt2` times (8) from
(9) gives

```math
Q(L(A'))-2\sqrt2 Q(A')
\ge2\sqrt2(b_{2N}-b_N)N^{3/2}+2\delta N^{3/2}
-o(N^{3/2}).                                          (10)
```

On infinitely many dyadic orders the first term is no more negative than
`-(delta/2)N^(3/2)`. Equations (7)--(10) prove (1), with a fixed finite
operator bound L+8 theta kappa.

One can instead choose dyadic orders with negative increment tending to
zero, set `delta_j=max(1/log N_j,sqrt(max(0,b_(N_j)-b_(2N_j))))`, and
take theta_j=delta_j/(g0 kappa²). Then (10) has excess at least
`N_j^(3/2)/log N_j`, while

```math
Q(A'_j)\le b_{N_j}N_j^{3/2}+o(N_j^{3/2}),\qquad
\|A'_j+I\|_{op}\le(L+o(1))\sqrt{N_j}.
```

This is asymptotic optimality RELATIVE TO THE CLASS (2), whose optimal cap
may exceed the original M_N. No equality `b_N-m_N->0` is assumed or proved.
Consequently this second statement does not establish bounded-op ORIGINAL
near-minimizers with the quantitative lift failure. It establishes that
bounded operator norm does not repair the subhalf-cap implication, while
leaving the original selected-exact-minimizer question open.
