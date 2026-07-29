# Concentration--compactness for competitive signings

Let
\[
M(A)=\max_{x\in\{\pm1\}^n}
\left|\sum_{i<j}a_{ij}x_ix_j\right|
\]
and suppose \(M(A)=O(n^{3/2})\).

This note separates two different concentration phenomena:

1. endpoint-energy bubbles, visible through large switched local fields;
2. operator bubbles, removed uniformly by Grothendieck--Pietsch
   factorization.

It then proves that both can be purified at arbitrarily small leading cost.
The remaining obstruction is not concentration but all-order realization of
the compact core profile.

## 1. Exact endpoint bubble profile

Suppose the positive endpoint is dominant and switch a maximizing spin to
\(\mathbf1\).  Write the switched signing as \(D\) and its row fields as
\[
r_i=(D\mathbf1)_i.
\]
Flipping coordinate \(i\) changes the quadratic energy by \(-2r_i\) in the
\(H=\sum_{i<j}\) normalization.  Global maximality therefore gives
\[
\boxed{r_i\ge0\quad\text{for every }i.}
\tag{1}
\]
Moreover
\[
2P(A)=\mathbf1^\top D\mathbf1=\sum_i r_i.
\tag{2}
\]
Define the empirical field measure
\[
\mu_n=\frac1n\sum_{i=1}^n\delta_{r_i/\sqrt n}.
\]
Then
\[
\boxed{
\frac{P(A)}{n^{3/2}}
=\frac12\int_0^\infty t\,d\mu_n(t).}
\tag{3}
\]
The first moments are uniformly bounded.  Along a subsequence the measures
converge weakly, but their first moments may lose a mass
\[
\beta
=\lim_n\int t\,d\mu_n-\int t\,d\mu
\ge0.
\tag{4}
\]
This is the exact endpoint bubble energy.

For a threshold \(K\), put
\[
S_K=\{i:r_i>K\sqrt n\}.
\]
Since all fields are nonnegative,
\[
|S_K|\le\frac{2P(A)}{K\sqrt n}=O(n/K).
\tag{5}
\]
Thus all escaped first moment lives on vanishing vertex sets as
\(K\to\infty\).

If \(R=S_K^c\), and
\[
D=\begin{pmatrix}D[S_K]&C\\C^\top&D[R]\end{pmatrix},
\]
then the endpoint has the exact decomposition
\[
2P(A)
=\mathbf1^\top D[S_K]\mathbf1
+\mathbf1^\top D[R]\mathbf1
+2\mathbf1^\top C\mathbf1.
\tag{6}
\]
Flipping all spins in \(S_K\) reverses the last term and preserves the two
internal terms.  Since the original spin is a global positive maximizer,
\[
\mathbf1^\top C\mathbf1\ge0.
\tag{7}
\]
Equations (3)--(7) are an exact concentration--compactness decomposition of
one endpoint.

## 2. Exact block energy aggregation, and its limit

For every vertex partition \(V_1,\ldots,V_k\),
\[
\sum_jP(A[V_j])\le P(A),\qquad
\sum_jQ(A[V_j])\le Q(A).
\]
Hence
\[
\boxed{
\sum_jW(A[V_j])\le W(A)\le M(A),}
\tag{8}
\]
where \(W=(P+Q)/2\).

For every union \(S\) of partition blocks,
\[
\boxed{
\|A[S,S^c]\|_{\infty\to1}\le W(A)\le M(A).}
\tag{9}
\]

These are the strongest scalar aggregation statements available.  There is
no quadrature or additive combination of (8) and (9): the order-five
counterexample in `pythagorean_centered_width_block.md` has a nonzero cross
block whose entire discrepancy fits inside the internal endpoint slack and
adds exactly zero width.

## 3. Bubble scales

If a bubble has size \(s=n^a=o(n)\), then:

- its trivial internal energy ceiling is \(O(s^2)\), which reaches the
  \(n^{3/2}\) scale at \(a=3/4\);
- a coherent \(s\times(n-s)\) cross block can have norm \(O(sn)\), which
  reaches the \(n^{3/2}\) scale already at \(a=1/2\);
- a random cross signing has
  \[
  \|C\|_{\infty\to1}
  =O\!\left(\sqrt{sn(s+n)}\right)
  =O(n\sqrt s)=o(n^{3/2}).
  \tag{10}
  \]

Thus vanishing bubbles can carry leading energy at several exponents.  A
single scalar atom at infinity does not encode their internal/cross
geometry.

## 4. Grothendieck--Pietsch core/bubble decomposition

Cube polarization gives
\[
\|A\|_{\infty\to1}^{\rm matrix}\le4M(A)=O(n^{3/2}).
\]
Grothendieck--Pietsch factorization implies that, for every fixed
\(0<\varepsilon<1\), there is a set \(S\) with
\[
|S|\le\varepsilon n
\tag{11}
\]
such that
\[
\boxed{
\|A[R]\|_{\rm op}\le C\frac{\sqrt n}{\varepsilon},
\qquad R=[n]\setminus S.}
\tag{12}
\]
Unlike the endpoint tail set, this single hub removes all positive and
negative operator spikes simultaneously.

## 5. Purification theorem

### Theorem

For every fixed \(\varepsilon>0\) and every competitive signing \(A\), there
is a signing \(\widetilde A\) of the same order such that
\[
M(\widetilde A)
\le M(A)+C\sqrt\varepsilon\,n^{3/2}+o(n^{3/2})
\tag{13}
\]
and
\[
\boxed{
\|\widetilde A\|_{\rm op}\le C_\varepsilon\sqrt n.}
\tag{14}
\]

### Proof

Take the split \(R\cup S\) from (11)--(12).

1. Keep \(A[R]\).
2. Fill \(S\) with a principal block of a conference matrix of order
   \((1+o(1))|S|\).  The resulting signing \(E\) obeys
   \[
   M(E)=O(|S|^{3/2}),\qquad
   \|E\|_{\rm op}=O(\sqrt{|S|}).
   \]
3. Fill \(S\times R\) by an independent random signing \(C\).  With positive
   probability,
   \[
   \|C\|_{\infty\to1}
   =O\!\left(\sqrt{|S||R|n}\right)
   =O(\sqrt\varepsilon\,n^{3/2})
   \]
   and
   \[
   \|C\|_{\rm op}=O(\sqrt n).
   \]

For
\[
\widetilde A=
\begin{pmatrix}E&C\\C^\top&A[R]\end{pmatrix},
\]
the block triangle inequality gives
\[
\begin{aligned}
M(\widetilde A)
&\le M(E)+M(A[R])+\|C\|_{\infty\to1}\\
&\le O(\varepsilon^{3/2}n^{3/2})
+M(A)
+O(\sqrt\varepsilon\,n^{3/2}),
\end{aligned}
\]
which is (13).  The operator-norm triangle inequality and (12) give (14).

### Consequence

For every \(\eta>0\), uniformly over all orders,
\[
\frac{F_n}{n^{3/2}}
\]
can be approximated within \(\eta\) by signings whose normalized
\(2\to2\) norm is bounded by a constant depending only on \(\eta\).
Vanishing bubbles are therefore variationally unnecessary for the minimum.

## 6. Action continuity on the purified class

Let \(T_n=\widetilde A_n/\sqrt n\), with
\(\|T_n\|_{2\to2}\le C_\varepsilon\).  If \(T_n\) action-converges, then the
Boolean objective
\[
\Phi(T_n)=\sup_{f\in\{\pm1\}}|\langle f,T_nf\rangle|
\]
is continuous along the sequence.

Indeed, the uniform \(L^2\) output bound makes the products
\(f\,T_nf\) uniformly integrable.  A maximizing Boolean input has a
subsequential one-profile limit, and its energy integral converges.
Conversely, a limiting profile whose input marginal is supported on
\(\{\pm1\}\) can be approximated by finite inputs and rounded back to signs;
the \(2\to2\) bound makes the rounding error in the quadratic energy tend to
zero.

Thus purification repairs exactly the uniform-integrability failure in
`action_convergence_boolean_spikes.md`.

## 7. Why convergence still does not follow

Compactness and continuity provide a limit object along a liminf
subsequence.  They do not show that this signed limit object is realizable
at every sufficiently large order with \(o(1)\) objective loss.

The logical gap admits a compact abstract model.  Let the compact profile
space be an interval, let the objective be the identity, and let the
admissible set at order \(n\) consist of
\[
c_n=1+\delta\sin\log\log(n+e^e).
\]
Then:

- every sequence has convergent subsequences;
- the objective is continuous;
- \(d(c_{n+1},c_n)\to0\);
- after multiplying by \(n^{3/2}\), the costs can be chosen convex and to
  satisfy the same superadditive, padding, and rectangular gluing
  inequalities recorded in `centered_width_rectangular_system.md`;
- nevertheless \(c_n\) oscillates.

This is a near-minimizer-compatible slow-drift countermodel with no bubbles
at all.

## 8. Verdict

The bubble problem can be regularized away at arbitrarily small leading
cost.  This is genuine progress: minimization reduces to uniformly
integrable, \(2\to2\)-bounded action profiles.  But the remaining problem is
an amplification/realization theorem:

> Every subsequential extremal signed action profile must be approximable,
> with the same normalization and objective, by sign matrices at every
> sufficiently large order.

Neither action compactness nor concentration--compactness supplies this
theorem, and without it slow drift remains possible.

