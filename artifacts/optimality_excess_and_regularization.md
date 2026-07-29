# Optimality excess versus spectral and row tails

## 1. Near-minimality does not imply unpeeled uniform integrability

Let
\[
q_n=\min_A Q(A)
\]
in the doubled quadratic normalization, and choose a subsequence
\(m\to\infty\) such that
\[
\frac{q_m}{m^{3/2}}\longrightarrow
c_{\inf}:=\liminf_n\frac{q_n}{n^{3/2}}.
\]
Take an exact minimizer \(B_m\), switch it and replace it by its negative
if necessary so that
\[
\mathbf1^\top B_m\mathbf1=q_m.
\]

Let \(k_m\to\infty\) with \(k_m=o(\sqrt m)\), and adjoin \(k_m\)
universally positive vertices:
\[
\widetilde B_m=
\begin{pmatrix}
J_k-I_k&J_{k,m}\\
J_{m,k}&B_m
\end{pmatrix}.
\]
The exact calculation gives
\[
Q(\widetilde B_m)
=q_m+2km+k(k-1).
\]
For \(n=m+k\),
\[
\frac{Q(\widetilde B_m)}{n^{3/2}}\to c_{\inf}.
\]
Thus this is an asymptotically minimizing sequence.

Nevertheless, applying \(\widetilde B_m\) to the normalized indicator
of the new vertices gives
\[
\|\widetilde B_m\|_{\rm op}\ge(1-o(1))\sqrt{km},
\]
and hence
\[
\frac{\|\widetilde B_m\|_{\rm op}}{\sqrt n}
\ge(1-o(1))\sqrt k\longrightarrow\infty.
\]
At its all-one maximizer, every new row has field \(n-1\). Therefore the
normalized row-square tail
\[
\frac1n\sum_i
\left(\frac{r_i}{\sqrt n}\right)^2
\mathbf1_{\{|r_i|>K\sqrt n\}}
\]
is at least \((1-o(1))k\) for every fixed \(K\), and diverges.

Consequently neither unpeeled spectral regularity nor row-square uniform
integrability can be bounded by the normalized optimality excess. The
exceptional set has only \(k=o(n)\) vertices, so a post-peeling theorem
remains possible. Notice also that
\[
\frac{(\sqrt{kn})^2}{n^{3/2}}
=\frac{k}{\sqrt n},
\]
the same scale as the normalized excess added by the construction. This
suggests that a Schatten-square spectral-tail charge, after localization,
is more plausible than a row-square charge.

## 2. A constructive excess--regularity tradeoff

There is a useful positive statement in the other direction. Let \(A_n\)
be an exact minimizer and fix \(\varepsilon>0\).
Grothendieck--Pietsch factorization supplies a principal set
\[
U\subset[n],\qquad |U|\ge(1-\varepsilon)n,
\]
such that
\[
\|A_n[U]\|_{\rm op}
\le \frac{4K_Gq_n}{\varepsilon n}
=O\!\left(\frac{\sqrt n}{\varepsilon}\right).
\]
Let \(S=[n]\setminus U\), \(h=|S|\).

Fill \(S\) with a conference-type signing of operator norm
\(O(\sqrt h)\), and fill the \(U\times S\) rectangle with independent
signs. A standard union bound gives a realization for which
\[
\max_{x\in\{\pm1\}^{U},\,y\in\{\pm1\}^{S}}
|x^\top By|
=O\!\left(\sqrt{nh(n+h)}\right)
=O(\sqrt\varepsilon\,n^{3/2}),
\]
while random-matrix norm bounds give
\[
\|B\|_{\rm op}=O(\sqrt n+\sqrt h)=O(\sqrt n).
\]
For the resulting full signing \(A'_n\),
\[
Q(A'_n)
\le Q(A_n[U])+O(h^{3/2})
+O(\sqrt\varepsilon\,n^{3/2})
\le q_n+O(\sqrt\varepsilon\,n^{3/2}),
\]
and
\[
\|A'_n\|_{\rm op}
=O\!\left(\frac{\sqrt n}{\varepsilon}\right).
\]

Equivalently, writing \(K=1/\varepsilon\), there exist full-order
signings with
\[
\boxed{
Q(A'_n)\le q_n+O(K^{-1/2}n^{3/2}),
\qquad
\|A'_n\|_{\rm op}=O(K\sqrt n).
}
\]
Taking \(K\to\infty\) arbitrarily slowly produces spectrally controlled
asymptotic near-minimizers. The tradeoff is not strong enough for the
current coefficient-one joint-selection argument, but it is a rigorous
two-limit regularization theorem.

## 3. Exact block lower bound relevant to replenishment

For
\[
A=\begin{pmatrix}D&B\\B^\top&E\end{pmatrix},
\]
fix \(y\) in the \(E\)-block. Averaging over a uniform Boolean
\(x\) in the \(D\)-block and then using the relative global sign of
\(x\) gives
\[
\boxed{
Q(A)\ge
|y^\top Ey|
+2\,\mathbb E_x|x^\top By|.
}
\]
By the sharp Khintchine inequality,
\[
\mathbb E_x|x^\top By|
\ge\frac1{\sqrt2}\|By\|_2.
\]
Thus
\[
Q(A)\ge
\max_y\left(
|y^\top Ey|+\sqrt2\,\|By\|_2
\right).
\]
This charges cross structure visible to a core energy witness. The
remaining replenishment difficulty is precisely that the high singular
directions of \(B\) may avoid every near-ground-state layer of \(E\).

There is a second, sharper form on the exact ground-state layer.  Let
\[
\operatorname{GS}(E)
=
\{y\in\{\pm1\}^{|E|}:|y^\top Ey|=Q(E)\}.
\]
For any \(y\in\operatorname{GS}(E)\), orient \(y\) so that
\(y^\top Ey=Q(E)\), and choose
\[
x=\operatorname{sign}(By).
\]
Changing the global sign of \(x\) changes only the cross term.  Hence
\[
Q(A)\ge Q(E)-Q(D)+2\|By\|_1,
\]
and therefore
\[
\boxed{
Q(A)-Q(E)
\ge
2\max_{y\in\operatorname{GS}(E)}\|By\|_1-Q(D).
}
\tag{3.1}
\]

Now apply (3.1) to a nested peeling tower.  Write
\[
A_t=
\begin{pmatrix}
D_t&B_t\\
B_t^\top&A_{t+1}
\end{pmatrix},
\qquad
d_t=Q(A_t)-Q(A_{t+1}),
\]
where the peeled vertex blocks are disjoint, and put
\[
V_t=
\max_{y\in\operatorname{GS}(A_{t+1})}\|B_ty\|_1.
\]
Then
\[
2V_t\le d_t+Q(D_t).
\]
The decrements telescope:
\[
\sum_t d_t\le Q(A_0).
\]
Also
\[
\sum_tQ(D_t)\le2Q(A_0).
\]
Indeed, choose an absolute ground state on each diagonal block; retain
the common energy-sign class carrying at least half of
\(\sum_tQ(D_t)\), and independently randomize the global signs of those
block states.  Some choice makes all cross terms nonnegative in the
chosen orientation.  Consequently
\[
\boxed{
\sum_tV_t\le\frac32Q(A_0).
}
\tag{3.2}
\]
Thus successive cross blocks have only \(O(Q(A_0))\) total
\(\ell_1\)-visibility on the successive exact ground-state frames.

## 4. The pointwise inverse is false, even for a regular core

It is tempting to hope that spectral regularity converts a large
replenishment gap \(g_t\) into comparable visibility \(V_t\).  This is
false even for a singleton deletion and
\(\|A_t\|_{\rm op}<1.5\sqrt{|A_t|}\).

Take the following order-nine signing:

\[
T=
\begin{pmatrix}
0&-1&1&1&-1&1&1&1&-1\\
-1&0&-1&-1&-1&-1&1&1&1\\
1&-1&0&-1&-1&1&-1&-1&1\\
1&-1&-1&0&-1&1&1&1&-1\\
-1&-1&-1&-1&0&-1&-1&-1&-1\\
1&-1&1&1&-1&0&-1&-1&-1\\
1&1&-1&1&-1&-1&0&1&1\\
1&1&-1&1&-1&-1&1&0&-1\\
-1&1&1&-1&-1&-1&1&-1&0
\end{pmatrix}
=
\begin{pmatrix}0&b\\b^\top&E\end{pmatrix}.
\]
Here
\[
b=(-1,1,1,-1,1,1,1,-1).
\]
Exact Boolean enumeration gives
\[
Q(T)=Q(E)=28.
\]
The full matrix has a positive ground state
\[
x=(+,-,+,+,-,+,+,+,-)
\]
whose restriction has \(E\)-energy \(12\).  Thus its replenishment gap
is
\[
g=28-12=16.
\]
On the other hand, the complete absolute ground-state set of \(E\) is
the pair
\[
\pm(-,+,-,+,+,-,-,-),
\]
and \(b^\top y=0\) on this pair.  Therefore
\[
\boxed{g=16,\qquad V=0.}
\tag{4.1}
\]

The regularity assertion has an exact rational certificate.  All
leading principal minors of \(81I-4T^2\) are positive:
\[
\begin{aligned}
&49,\ 2385,\ 115425,\ 4688865,\ 203295825,\\
&6488617905,\ 251693598465,\ 6095015474625,\
197299523840625.
\end{aligned}
\]
Sylvester's criterion implies
\[
81I-4T^2\succ0,
\qquad
\|T\|_{\rm op}<\frac92=\frac32\sqrt9.
\]

So no inequality of the form
\[
g_t\le C(K)V_t
\]
can hold pointwise under \(\|A_t\|_{\rm op}\le K\sqrt{|A_t|}\), even
for \(K=3/2\).  The cumulative estimate (3.2) remains useful, but any
inverse theorem must group several scales, enlarge the visible layer
from exact ground states to near-ground states, or charge a different
quantity.

## 5. An existential regularized adaptivity-gap theorem

Although the pointwise visibility inverse fails, spectral regularity
does control cumulative replenishment for a suitable deletion order.

**Theorem.**  Suppose \(A\) has order \(n\) and
\[
\|A\|_{\rm op}\le K\sqrt n.
\]
There is a singleton principal-deletion order, with an arbitrary
choice of positive ground state at each suffix, for which
\[
\boxed{
\sum_{t=0}^{n-2}g_t
\le
2Kn^{3/2}-Q(A).
}
\tag{5.1}
\]
Equivalently, choosing the next deleted vertex uniformly at random
from the current core gives (5.1) in expectation.

To prove this, let \(S\) be the current core, \(m=|S|\), and let
\(x\) be a positive ground state:
\[
x^\top A[S]x=Q(A[S]).
\]
For \(i\in S\), put
\[
\ell_i=x_i(A[S]x)_i.
\]
Boolean one-flip optimality gives \(\ell_i\ge0\).  If \(i\) is deleted,
the old state's oriented energy on the new core is
\[
e_i=Q(A[S])-2\ell_i.
\]
Writing
\[
d_i=Q(A[S])-Q(A[S\setminus\{i\}]),
\qquad
g_i=Q(A[S\setminus\{i\}])-e_i,
\]
gives the exact singleton identity
\[
2\ell_i=d_i+g_i.
\tag{5.2}
\]

For uniform \(i\in S\),
\[
\mathbb E_i\ell_i^2
=\frac1m\|A[S]x\|_2^2
\le\|A[S]\|_{\rm op}^2
\le\|A\|_{\rm op}^2
\le K^2n.
\]
Hence
\[
\mathbb E_i\ell_i\le K\sqrt n.
\]
Summing the conditional expectations along the random deletion chain,
using (5.2) and
\[
\sum_td_t=Q(A),
\]
proves (5.1).  The deterministic version follows by choosing at every
step a coordinate with
\(\ell_i\le(\frac1m\sum_j\ell_j^2)^{1/2}\le K\sqrt n\).

If \(Q(A)\ge c_0n^{3/2}\) is any universal lower bound, then
\[
\sum_tg_t
\le
\left(\frac{2K}{c_0}-1\right)Q(A).
\tag{5.3}
\]
Thus spectrally regular signings have an existential
\(C(K)\)-bounded suffix-FTL adaptivity gap.

The quantifier on the order matters.  The construction preferentially
deletes a coordinate with a small current local field.  It does not
show that an order constrained to delete current heavy coordinates has
bounded replenishment, and therefore does not by itself finish the
heavy-core peeling argument.

## 6. Two attempted upgrades and their exact losses

### 6.1 Capped field-biased deletion

On a current order-\(m\) core with local fields \(\ell_i\ge0\), consider
\[
w_i=\min\left(1,\frac{\ell_i}{H\sqrt m}\right),
\qquad
p_i=\frac{w_i}{\sum_jw_j}.
\]
Even under the scale-local assumption
\(\|A[S]\|_{\rm op}\le K\sqrt m\), the direct estimates are only
\[
\sum_iw_i
\le\frac{\sum_i\ell_i}{H\sqrt m}
\le\frac{Km}{H}
\tag{6.1}
\]
and, while at least one \(H\)-heavy coordinate exists,
\[
\mathbb E_p\ell_i
\le
\frac1{H\sqrt m}\sum_i\ell_i^2
\le\frac{K^2m^{3/2}}H.
\tag{6.2}
\]
Here the denominator in \(p_i\) is bounded below only by \(1\).
Summing (6.2) through \(\Theta(n)\) deletions costs
\(O(K^2n^{5/2}/H)\), a factor \(n/H\) above the required
\(n^{3/2}\) scale.

If there is just one heavy coordinate, (6.1) gives selection
probability only
\[
\frac{H}{Km}.
\]
This supplies constant hazard over \(m/H\) successive states only when
the same heavy coordinate persists.  It does not control a coordinate
that becomes heavy late, or a tower in which replenishment continually
moves the heavy field.  Thus capped linear bias alone does not yield a
useful probability that a large stopped core is heavy-free.

### 6.2 Average one-vertex deletion

For the first uniformly random deletion from a positive ground state,
\[
\mathbb E_i\ell_i
=\frac1n\sum_i x_i(Ax)_i
=\frac{Q(A)}n.
\]
Together with \(2\ell_i=d_i+g_i\), this gives the exact identity
\[
\frac1n\sum_i\bigl(Q(A)-Q(A_{-i})\bigr)
+\frac1n\sum_i g_i
=\frac{2Q(A)}n.
\tag{6.3}
\]
In particular, the automatic conclusion is the lower bound
\[
\frac1n\sum_iQ(A_{-i})
\ge\left(1-\frac2n\right)Q(A).
\]
The convergence-useful upper bound
\[
\frac1n\sum_iQ(A_{-i})
\le\left(1-\frac{3}{2n}+o(n^{-1})\right)Q(A)
\]
is exactly equivalent to the new estimate
\[
\frac1n\sum_i g_i
\le\left(\frac12+o(1)\right)\frac{Q(A)}n.
\tag{6.4}
\]
The cumulative theorem (5.1) does not localize sharply enough to imply
(6.4).

Regularity alone cannot prove it.  The matrix \(T\) in Section 4
satisfies
\[
Q(T_{-i})=Q(T)=28
\quad\text{for every }i,
\]
although \(\|T\|_{\rm op}<1.5\sqrt9\).  There is also an order-seven
signing with \(Q=18\), the exact order-seven minimum in the doubled
normalization, for which all seven one-vertex principal norms equal
\(18\).  Therefore any surviving version of (6.4) must use
large-order asymptotic minimality, not merely spectral regularity or
finite-order exact minimality.
