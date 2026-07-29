# Spectral peeling: a dense-sign counterexample

## Claim that fails

The following plausible statement is false:

> If \(A_n\) is a symmetric zero-diagonal sign matrix and
> \(Q(A_n)=O(n^{3/2})\), then there are sets \(S_n=o(n)\) such that
> \(\|A_n[S_n^c]\|_{\rm op}=O(\sqrt n)\), with one uniform implied
> constant.

The obstruction is an integrable tower of positive cliques whose
normalized sizes have unbounded support.

## Construction

For convenience take \(n=N^2\). Fix
\[
p_j=\frac{j^{-3}}{\sum_{\ell\ge1}\ell^{-3}},
\qquad
K_j=aj,
\]
where \(a>0\) is fixed. Let \(J_N\to\infty\) sufficiently slowly, for
example \(J_N=\lfloor N^{1/5}\rfloor\). For \(1\le j\le J_N\), form
\[
b_j=\left\lfloor\frac{p_jn}{k_j}\right\rfloor
\quad\text{blocks of size}\quad
k_j=\lfloor K_jN\rfloor.
\]
Put any unassigned \(o(n)\) vertices into additional blocks of size at
most \(O(N)\).

Set every edge inside every block equal to \(+1\). Set the edges between
different blocks to independent signs, symmetrically.

The rounding and truncation errors in the partition are \(o(n)\), since
\[
\sum_{j\le J_N} k_j=O(NJ_N^2)=o(n),
\qquad
\sum_{j>J_N}p_j=o(1).
\]

## The quadratic norm is \(O(n^{3/2})\)

For a block \(B\) of size \(k\),
\[
x_B^\top(J_k-I_k)x_B
=\left(\sum_{i\in B}x_i\right)^2-k,
\]
so its absolute contribution is at most \(k^2\). Moreover,
\[
\sum_B |B|^2
\le
\left(a\sum_{j\ge1}p_jj+o(1)\right)n^{3/2}.
\]

Let \(C\) be the matrix containing only the random between-block edges.
For each fixed Boolean \(x\), \(x^\top Cx\) is a sum of at most
\(\binom n2\) independent variables of magnitude \(2\). Hoeffding's
inequality followed by a union bound over the \(2^n\) Boolean vectors
shows that, with positive probability,
\[
\max_x|x^\top Cx|\le C_0n^{3/2}
\]
for an absolute constant \(C_0\). Hence a deterministic realization
exists for which
\[
Q(A)
\le Q(C)+\sum_B|B|^2
=O(n^{3/2}).
\]

## No \(o(n)\) deletion leaves a uniformly regular core

Suppose \(|S_n|=o(n)\), and fix any proposed constant \(C>0\). Choose a
fixed \(j\) with \(K_j>2(C+1)\). Type \(j\) occupies
\[
b_jk_j=(p_j+o(1))n
\]
vertices.

To make every remaining type-\(j\) block have at most
\((C+1)\sqrt n\) vertices, one must delete at least
\[
b_j\bigl(k_j-(C+1)\sqrt n\bigr)
\ge \left(\frac{p_j}{2}+o(1)\right)n
\]
vertices. This contradicts \(|S_n|=o(n)\). Thus some type-\(j\) block
retains more than \((C+1)\sqrt n\) vertices.

The vector supported uniformly on those retained vertices has Rayleigh
quotient \(r-1>C\sqrt n\), since the corresponding principal block is
\(J_r-I_r\). Therefore
\[
\|A[S_n^c]\|_{\rm op}>C\sqrt n.
\]
As \(C\) was arbitrary, no \(o(n)\) deletion gives a uniform
\(O(\sqrt n)\) operator norm.

## Arbitrarily small fixed leading perturbation

The same obstruction can be implanted in any signing \(B_n\) with
\(Q(B_n)=O(n^{3/2})\). Keep every between-block edge of \(B_n\), and
overwrite only each block interior with \(+1\). If \(F\) edges are
changed, then
\[
Q(A_n)\le Q(B_n)+4F
\le Q(B_n)+2\sum_B|B|^2.
\]
Since
\[
\sum_B|B|^2
=\left(a\sum_jp_jj+o(1)\right)n^{3/2},
\]
the added normalized cost can be made arbitrarily small, but fixed, by
choosing \(a\) small. Meanwhile \(K_j=aj\) remains unbounded, so the
failure of uniform regularization persists.

It cannot be made an \(o(n^{3/2})\) perturbation by this construction
while preserving the obstruction. If
\(\sum_jp_jK_j=o(1)\), Markov's inequality makes the vertex mass in
blocks with \(K_j>C\) equal to \(o(1)\) for every fixed \(C\), and those
vertices can be deleted.

## Consequence for the research route

An absolute peeling theorem based only on
\(Q(A)=O(n^{3/2})\) is impossible. A viable replacement must charge the
spectral-tail mass to the excess
\[
Q(A)-(\text{asymptotically minimal value})\,n^{3/2},
\]
or use an optimizer-specific replenishment statistic.

The tower inequality
\[
2\sum_tR_t\le3Q(A_0)+\sum_tg_t
\]
only bounds an analogue of the first moment
\(\sum_jp_jK_j\). The construction has finite first moment but unbounded
support. Thus such a bound alone permits a positive (though arbitrarily
small) linear fraction of vertices above every fixed spectral
threshold. Controlling the cumulative replenishment by the *optimality
excess*, or proving uniform integrability of the tower tail, is the
remaining necessary ingredient.

## The optimal qualitative replacement

Grothendieck--Pietsch factorization does give an
\(\varepsilon\)-regular core. Put
\[
B=\|A\|_{\infty\to1}
=\max_{x,y\in\{\pm1\}^n}|x^\top Ay|
\le2Q(A).
\]
There are probability weights \((\mu_i)\), \((\nu_j)\) such that
\[
|u^\top Av|
\le K_GB
\left(\sum_i\mu_i u_i^2\right)^{1/2}
\left(\sum_j\nu_j v_j^2\right)^{1/2}.
\]
Delete every index for which either
\[
\mu_i>\frac2{\varepsilon n}
\quad\text{or}\quad
\nu_i>\frac2{\varepsilon n}.
\]
At most \(\varepsilon n\) indices are deleted. On the remaining
principal set \(U\),
\[
\boxed{
\|A[U]\|_{\rm op}
\le\frac{2K_GB}{\varepsilon n}
\le\frac{4K_GQ(A)}{\varepsilon n}.
}
\]
Thus \(Q(A)=O(n^{3/2})\) always gives
\[
|U|\ge(1-\varepsilon)n,
\qquad
\|A[U]\|_{\rm op}=O_\varepsilon(\sqrt n).
\]
The clique-tower construction shows why the constant must be allowed to
diverge as \(\varepsilon\downarrow0\); its \(1/\varepsilon\)-type
dependence is qualitatively sharp under only a first-moment budget.
