# Subhalf caps do not imply H2 lift stability

Status: complete scalable actual-sign argument, pending independent audit.
This does NOT falsify a selectable EXACT-minimizer H2 theorem, or a merely
vanishing-normalized-loss H2 theorem for asymptotic minimizers.

Write `q_A(x)=x^T A x/2`, `Q(A)=max_x|q_A(x)|`,
`m_n=M_n/n^(3/2)`, and define the actual hollow signing

```math
L(A)=\begin{pmatrix}A&A+I\\A+I&-A\end{pmatrix}.
```

Assume only a proved all-order bound `limsup m_n<U<1/2`, available in this
campaign. Then:

1. There are actual signings A_j of unbounded orders N_j and a fixed eta>0
   for which `Q(A_j)<N_j^(3/2)/2`, but
   `Q(L(A_j))>=2 sqrt(2) Q(A_j)+eta N_j^(3/2)`.
2. There are actual asymptotic minimizers
   `Q(A_j)=M_(N_j)+o(N_j^(3/2))`, still strictly subhalf, for which
   `Q(L(A_j))-2 sqrt(2) Q(A_j)>=N_j^(3/2)/log N_j` eventually.

Thus subhalf cap alone cannot give even vanishing normalized H2 error.
Asymptotic minimality alone cannot give O(N) error, or any prescribed
power-saving error, for EVERY near-minimizing sequence.

## 1. A two-block flat gadget with an excessive H2 ratio

Put `rho=sqrt(2)-1`, and for each k choose a k-by-k signing C satisfying

```math
\|C-\rho J\|_{\infty\to1}\le4k^{3/2}.                 (1)
```

Such actual signings exist for all sufficiently large k: choose independent
entries of mean rho. For each Boolean pair, Hoeffding bounds the centered
bilinear sum tail at `4 k^(3/2)` by `2 exp(-8k)`; union bound over `4^k`
pairs suffices. Maxima over cubes and Boolean pairs agree.

Define the hollow signing on 2k vertices

```math
B_k=\begin{pmatrix}J-I&C\\C^T&-J+I\end{pmatrix},
\qquad g_0={1+\rho^2\over2}.
```

For two block spins with normalized sums a,b, its energy is
`(k^2/2)(a^2-b^2+2 rho ab)` plus an error bounded by `4 k^(3/2)`.
The maximum absolute value of the displayed quadratic on `[-1,1]^2`
is `1+rho^2`, attained at `(a,b)=(1,rho)` or its sign counterparts.
Rounding the magnetization b costs O(k). Consequently

```math
Q(B_k)=g_0k^2+O(k^{3/2}).                              (2)
```

To lower-bound its lift, take p equal to 1 on the first block and 0 on
the second, and r equal to 0 on the first and 1 on the second. Use lift
spins `(x,y)=(p+r,p-r)`. The diagonal-completion term x dot y is zero.
The other energy equals

```math
2[q_{B_k}(p)-q_{B_k}(r)+p^TB_kr]
=2(1+\rho)k^2+O(k^{3/2}).
```

Therefore, with `R=2(1+sqrt(2))`,

```math
Q(L(B_k))\ge R g_0 k^2-O(k^{3/2}),\qquad
R-2\sqrt2=2.                                          (3)
```

## 2. Planting preserves a lift extremizer up to O(k sqrt N)

Start with an EXACT minimizing signing A of order N. For all sufficiently
large N, `Q(A)<=U N^(3/2)`. The elementary polarization estimate gives

```math
Q(L(A))\le6Q(A)+N\le6U N^{3/2}+N.                     (4)
```

Also the lift energy changes sign under `(x,y)->(y,-x)`. Choose a positive
lift extremizer z. Every signed local field
`ell_i=z_i(L(A)z)_i` is nonnegative, because flipping coordinate i cannot
increase its energy, and `sum ell_i=2Q(L(A))`.

For each original index i, pair the two fields i and N+i. At least three
quarters of original indices satisfy

```math
\ell_i+\ell_{N+i}\le(48U+8)\sqrt N.                   (5)
```

Independently, the simultaneous Grothendieck diagonal-majorant theorem
supplies `D>=+/-A` and `Tr D<=4K_G Q(A)`. Hence at least three quarters of
indices have `D_ii<=16K_G U sqrt(N)`. The intersection with (5) has size
at least N/2. Select any set S of size 2k inside it, where `2k<=N/2`.

The induced old block obeys

```math
\|A_S\|_{op}\le16K_G U\sqrt N,\qquad
Q(A_S)=O_U(k\sqrt N),\qquad
Q(L(A_S))=O_U(k\sqrt N).                               (6)
```

Replace A_S by B_k, leaving every other edge unchanged, and call the
result A'. The original triangle inequality and (2),(6) give

```math
Q(A')\le M_N+g_0k^2+O_U(k^{3/2}+k\sqrt N).             (7)
```

For clarity, the lift lower estimate does NOT assume that old and new
extremizers happen to agree. Change the 4k selected coordinates of z to
an extremizer of L(B_k). If T is the subset of coordinates actually flipped,
the exact spin-flip identity is

```math
q_{L(A)}(z^{T})-q_{L(A)}(z)
=-2\sum_{i\in T}\ell_i+4q_{L(A)_T}(z_T).
```

By (5) the first loss is O_U(k sqrt(N)), and by principal monotonicity
and (6) the second loss is too. Replacing the old selected principal block
then loses at most another `Q(L(A_S))` and gains `Q(L(B_k))`. The identity
of the selected principal block with L(A_S) includes its cross-diagonal
completion; those unchanged diagonal-completion edges cancel exactly when
the replacement is made. Thus

```math
Q(L(A'))\ge Q(L(A))+R g_0k^2
-O_U(k^{3/2}+k\sqrt N).                                (8)
```

In particular, since `Q(L(A))>=M_(2N)`, equations (7)--(8) yield

```math
Q(L(A'))-2\sqrt2 Q(A')
\ge 2\sqrt2\,(m_{2N}-m_N)N^{3/2}+2g_0k^2
-O_U(k^{3/2}+k\sqrt N).                                (9)
```

This is an original-signing statement: no zero coefficients, fractional
final entries, matrix relaxation, or seed realization assumption remains.

## 3. A bounded dyadic sequence supplies the necessary orders

For any bounded sequence `(m_(2^j))`,

```math
\limsup_{j\to\infty}(m_{2^{j+1}}-m_{2^j})\ge0.
```

Otherwise its increments are eventually uniformly negative, contradicting
boundedness below. Hence there are arbitrarily large N=2^j with the
increment in (9) as close to nonnegative as desired. No convergence of m
is being assumed here.

For the first assertion, fix `0<delta<1/2-U` and set
`k=floor(sqrt(delta/g_0) N^(3/4))`. Along infinitely many dyadic orders
with `m_(2N)-m_N>=-delta/(4 sqrt(2))`, the error in (9) is o(N^(3/2)).
Equation (7) gives `Q(A')<N^(3/2)/2` eventually, while (9) gives a lift
excess at least `delta N^(3/2)` eventually.

For the second assertion choose a dyadic subsequence N_j along which
`b_j=max(0,m_(N_j)-m_(2N_j))->0`, and put

```math
\delta_j=\max\{1/\log N_j,\sqrt{b_j}\},\qquad
k_j=\lfloor\sqrt{\delta_j/g_0}\,N_j^{3/4}\rfloor.
```

Then delta_j tends to zero, the dyadic increment is at least `-delta_j^2`,
and the error in (9) is `o(delta_j N_j^(3/2))`, since
`k_j/sqrt(N_j)->infinity`. The excess is eventually at least
`delta_j N_j^(3/2)>=N_j^(3/2)/log N_j`. Equation (7), together with
minimality as a lower bound, proves asymptotic minimality of A'_j.

## Scope and the surviving exact-minimizer question

The modified matrices need not be exact minimizers. This theorem therefore
does not falsify the possibility that EACH order admits an EXACT minimizer
A with `Q(L(A))<=2 sqrt(2) M_N+O(N)`, or another selected hereditary
construction with a summable error. The second family has a vanishing
normalized lift excess, so it also does not falsify a qualitative o(N^1.5)
claim for all near-minimizers. Its point is that neither a strict subhalf
cap nor unquantified near-minimality supplies the desired stronger premise.

The good dyadic orders need not be algorithmically known; existence of the
signings is unconditional given the stated all-order subhalf bound. No
nonconvergence of M_n is claimed.
