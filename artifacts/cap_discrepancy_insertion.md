# Cap discrepancy and one-vertex insertion

## Status

This note audits the direct one-vertex route to convergence.  It first
records the exact scale which an insertion theorem would have to meet.
In particular, a total \(o(\sqrt n)\) insertion increment is impossible
uniformly along minimizers, and an \(o(\sqrt n)\) *error beyond the
scale-correct derivative* is still not, by itself, enough to force a
limit.

For two exact affine cap faces, the insertion problem becomes an
explicit signed bipartite-incidence discrepancy problem.  The
frustration index of that signed graph gives a deterministic bound on
the fields of every exact top and bottom state.  The connected
components simultaneously give exact zero-cut principal closures.
This isolates the direct route's remaining analytic obstruction:
controlling the thick cap after the exact endpoint fields have been
balanced.

Throughout,
\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
M(A)=\max_x|H_A(x)|,
\]
and \(M_n=\min_A M(A)\).

## 1. Exact insertion identity

Adjoin a vertex whose incident signs are
\(b=(b_1,\ldots,b_n)\in\{\pm1\}^n\).  For the old spin vector \(x\)
and the new spin \(y\),
\[
H_{A,b}(x,y)=H_A(x)+y\,b\cdot x.
\]
For real \(h,t\),
\[
\max_{y=\pm1}|h+yt|=|h|+|t|.
\]
It follows that
\[
\boxed{
E(A):=
\min_b\max_{x,y}|H_{A,b}(x,y)|
=
M(A)+\min_b\Delta_A(b),
}
\tag{1.1}
\]
where
\[
g_A(x)=M(A)-|H_A(x)|,\qquad
\boxed{
\Delta_A(b)=\max_x\bigl(|b\cdot x|-g_A(x)\bigr).
}
\tag{1.2}
\]
In particular,
\[
M_{n+1}\le E(A)
\]
for every order-\(n\) signing \(A\), with equality after minimizing
over \(A\).

Thus an insertion bound \(\Delta_A(b)\le d\) is exactly the
simultaneous weighted-discrepancy system
\[
\boxed{|b\cdot x|\le d+g_A(x)\quad\text{for every }x.}
\tag{1.3}
\]
Equivalently, with
\[
d_\pm(b,x)=\min\{d_H(b,x),d_H(b,-x)\},
\]
it is the intersection-of-Hamming-caps condition
\[
\boxed{
d_\pm(b,x)\ge \frac{n-g_A(x)-d}{2}
\quad\text{for every }x.
}
\tag{1.4}
\]
The exact endpoint constraint is ordinary discrepancy of the union of
the positive and negative ground-state families.  The nonendpoint
states enter with precisely their energy slack.

## 2. The scale-correct differential target

Put
\[
f_n=\frac{M_n}{n^{3/2}},\qquad
D_n=M_n\left[\left(1+\frac1n\right)^{3/2}-1\right].
\tag{2.1}
\]
If an optimal order-\(n\) signing admits a row with
\[
\Delta_A(b)\le D_n,
\tag{2.2}
\]
then
\[
\frac{M_{n+1}}{(n+1)^{3/2}}
\le
\frac{M_n}{n^{3/2}}.
\tag{2.3}
\]
Thus (2.2) eventually at every order would prove convergence by
eventual monotonicity of \(f_n\).

There is a robust summable-error form.  If
\[
\Delta_A(b)\le D_n+r_n,\qquad r_n\ge0,
\tag{2.4}
\]
for an optimal \(A\) at every sufficiently large \(n\), and
\[
\boxed{\sum_n\frac{r_n}{n^{3/2}}<\infty,}
\tag{2.5}
\]
then \(f_n\) converges.  Indeed,
\[
(f_{n+1}-f_n)_+
\le\frac{r_n}{(n+1)^{3/2}},
\]
so \(f_n\) has finite total positive variation; since it is bounded
below, it has a limit.

By contrast, merely replacing (2.5) by
\[
r_n=o(\sqrt n)
\tag{2.6}
\]
does not imply convergence.  After normalization, (2.6) only says
\[
f_{n+1}\le f_n+o(1/n).
\]
Bounded nonconvergent sequences such as
\[
f_n=c+\varepsilon\sin(\log\log n)
\]
satisfy this one-sided local condition.  A one-step theorem therefore
needs an exact sign, a summable error, or additional nonlocal control
on multiplicative intervals.

An equivalent nonlocal sufficient statement is: if
\[
\boxed{
M_N\le\left(\frac Nn\right)^{3/2}M_n+\varepsilon_nN^{3/2}
\quad(N\ge n),\qquad \varepsilon_n\to0,
}
\tag{2.7}
\]
then \(f_n\) converges.  Take \(n\) along a subsequence realizing the
liminf and then let \(N\ge n\) be arbitrary.

## 3. A total \(o(\sqrt n)\) increment is impossible

Suppose that for every sufficiently large \(n\), some optimal
order-\(n\) signing had a row with
\[
\Delta_A(b)=o(\sqrt n)
\]
uniformly in the tail.  Then (1.1) would give
\[
M_{n+1}-M_n=o(\sqrt n).
\]
Summation would imply \(M_n=o(n^{3/2})\), contradicting the universal
positive lower bound.

The present rigorous cage makes this quantitative.  Write
\[
c_*=0.336493364431\ldots
\]
for the proved lower constant and use
\(\limsup M_n/n^{3/2}\le1/2\).  For dyadic intervals,
\[
\begin{aligned}
M_{2N}-M_N
&\ge
\left(2^{3/2}c_*-\frac12-o(1)\right)N^{3/2}\\
&=
(0.451746\ldots-o(1))N^{3/2}.
\end{aligned}
\tag{3.1}
\]
Consequently, in every sufficiently large dyadic interval there is
some \(n\in[N,2N)\) for which
\[
\boxed{
M_{n+1}-M_n
\ge
\left(
\frac{2^{3/2}c_*-1/2}{\sqrt2}-o(1)
\right)\sqrt n
=
(0.3194\ldots-o(1))\sqrt n.
}
\tag{3.2}
\]
For every optimal \(A\) at such an order,
\[
\min_b\Delta_A(b)\ge M_{n+1}-M_n.
\]
Thus the direct route cannot aim to make the *total* insertion
increment \(o(\sqrt n)\).  It must hit the derivative-scale target in
Section 2.

## 4. Two affine cap faces: exact signed-incidence reduction

Assume the positive and negative exact ground families contain affine
families.  Applying exact affine closure separately to the two
orientations enlarges them to the full type-constant families
\[
\mathcal G_+
=
\{x_i=\alpha_i u_{p(i)}:u\in\{\pm1\}^{P}\},
\qquad
\mathcal G_-
=
\{y_i=\beta_i v_{q(i)}:v\in\{\pm1\}^{Q}\}.
\tag{4.1}
\]
Here \(p(i)\) and \(q(i)\) are the occupied positive- and
negative-face types.

Construct a bipartite multigraph \(G\) with left vertices \(P\), right
vertices \(Q\), and one edge \(i\) joining \(p(i)\) to \(q(i)\).  Give
edge \(i\) the sign
\[
\gamma_i=\alpha_i\beta_i.
\tag{4.2}
\]
For a proposed insertion row put
\[
c_i=b_i\alpha_i.
\tag{4.3}
\]
Its type imbalances are
\[
L_p(c)=\sum_{i:p(i)=p}c_i,\qquad
R_q(c)=\sum_{i:q(i)=q}\gamma_ic_i.
\tag{4.4}
\]
Because all type signs in (4.1) are independently available,
\[
\boxed{
\max_{x\in\mathcal G_+}|b\cdot x|
=\sum_{p\in P}|L_p(c)|,\qquad
\max_{y\in\mathcal G_-}|b\cdot y|
=\sum_{q\in Q}|R_q(c)|.
}
\tag{4.5}
\]
Thus the exact two-face insertion discrepancy is
\[
\boxed{
\kappa(G,\gamma)
=
\min_{c\in\{\pm1\}^{E(G)}}
\max\left\{
\|L(c)\|_1,\|R(c)\|_1
\right\}.
}
\tag{4.6}
\]
The relevant norm is \(\ell_1\), not merely maximum vertex
discrepancy.  This is forced by the Cartesian type-constant ground
closure and prevents residuals in disconnected components from being
hidden.

## 5. Frustration controls exact-cap insertion

Call a signed bipartite graph balanced if there are vertex signs
\(s_p,t_q\) such that
\[
\gamma_i=s_{p(i)}t_{q(i)}
\tag{5.1}
\]
on every edge.  Equivalently, the product of the edge signs on every
cycle is \(+1\).  Let
\[
\tau(G,\gamma)
=
\min_{s,t}
|\{i:\gamma_i\ne s_{p(i)}t_{q(i)}\}|
\tag{5.2}
\]
be its frustration index.

The affine parity theorem says that each of the two type partitions
has at most one odd class.  Therefore, if \(n\) is even, every vertex
of \(G\) has even degree.  If \(n\) is odd, there is one odd vertex on
each side; the two lie in the same connected component.

### Proposition 5.1

Under this parity hypothesis,
\[
\boxed{\kappa(G,\gamma)\le 2\tau(G,\gamma)+1.}
\tag{5.3}
\]
For even \(n\), the \(+1\) can be omitted.

#### Proof

Choose vertex signs \(s,t\) realizing (5.2), and call the exceptional
edge set \(F\).  Gauge the edge colors by
\[
d_i=s_{p(i)}c_i.
\]
Then
\[
L_p(c)=s_p\sum_{i\ni p}d_i,
\]
while
\[
R_q(c)
=t_q\left(
\sum_{i\ni q}d_i
-2\sum_{\substack{i\in F\\i\ni q}}d_i
\right).
\tag{5.4}
\]

If all degrees are even, alternate \(d_i=\pm1\) around an Euler tour
in each connected component.  The ordinary incidence sum
\(\sum_{i\ni v}d_i\) then vanishes at every vertex.  Hence
\[
\|L(c)\|_1=0,\qquad
\|R(c)\|_1
\le2|F|=2\tau(G,\gamma).
\]

If there are two odd vertices, use an alternating Euler trail in
their common component and alternating Euler tours elsewhere.  The
ordinary incidence imbalance has magnitude one at the two endpoints
and vanishes elsewhere.  Equation (5.4) gives
\[
\|L(c)\|_1\le1,\qquad
\|R(c)\|_1\le1+2|F|.
\]
This proves (5.3). \(\square\)

In particular,
\[
\boxed{
\tau(G,\gamma)=o(\sqrt n)
\Longrightarrow
\max_{\mathcal G_+\cup\mathcal G_-}|b\cdot x|
=o(\sqrt n)
}
\tag{5.5}
\]
for an explicit deterministic insertion row.

This controls only the exact endpoints.  To reach (1.3), one must
show that the same row has
\[
|b\cdot x|\le D_n+g_A(x)
\]
throughout the thick cap.  Ordinary entropy union bounds are known to
fail at the needed constants.

## 6. Connected components are simultaneous principal closures

Every connected component \(C\) of the signed-incidence graph is a
union of positive types and also a union of negative types.  In the
positive endpoint gauge, every inter-positive-type block has total
zero.  Consequently
\[
\sum_{ij\in\delta(C)}a_{ij}x_i^+x_j^+=0.
\tag{6.1}
\]
The zero-cut principal-closure theorem therefore applies to every
union of incidence components.  In particular,
\[
\boxed{
W(A)\ge W(A[C])+W(A[C^c]),
}
\tag{6.2}
\]
with Cartesian positive-ground closure and cross-block annihilation
on the principal ground spans.

Thus the affine insertion branch has an exact trichotomy:

1. \(\tau=o(\sqrt n)\), in which case both exact cap fields admit an
   \(o(\sqrt n)\) simultaneous balancing row;
2. one incidence component has size \(n-o(n)\), in which case its
   principal restriction is a scale-preserving descent candidate;
3. frustration \(\Omega(\sqrt n)\) remains distributed over a
   genuinely multipart component decomposition.

The third case is the same double-centered residual obstruction that
appears in the affine refill analysis.  Resolving it requires a
same-signing theorem turning distributed frustration into either an
additional \(n^{3/2}\)-scale width payment or a nearly macroscopic
principal block.

## 7. Exact endpoint balance does not control the thick cap

There is a concrete signing showing that no deterministic implication
of the form
\[
\max_{\{|H_A|=M(A)\}}|b\cdot x|\ \hbox{small}
\quad\Longrightarrow\quad
\Delta_A(b)=O(\sqrt n)
\tag{7.1}
\]
can hold for the same row \(b\), even when the endpoint field is
exactly zero.

Let \(k\ge8\) be even.  Let \(B\) be the \(k\times k\) circulant sign
matrix whose first row consists of \(k/2\) consecutive \(+1\)'s
followed by \(k/2\) consecutive \(-1\)'s.  Then
\[
B\mathbf1=B^\mathsf T\mathbf1=0,\qquad
\|B\|_{\rm op}=\frac2{\sin(\pi/k)}<k-2.
\tag{7.2}
\]
Define the order-\(2k\) signing
\[
A_0=
\begin{pmatrix}
J_k-I_k&B\\
B^\mathsf T&J_k-I_k
\end{pmatrix}.
\tag{7.3}
\]
For spins \(x,y\in\{\pm1\}^k\), put
\(s=\mathbf1^\mathsf Tx\), \(t=\mathbf1^\mathsf Ty\).  Since the
cross block annihilates constants,
\[
x^\mathsf TBy
=
(x-s\mathbf1/k)^\mathsf TB(y-t\mathbf1/k).
\]
Writing \(X=k^2-s^2\), \(Y=k^2-t^2\), one obtains
\[
H_{A_0}(x,y)
\le
k(k-1)-\frac{X+Y}{2}
+\frac{\|B\|_{\rm op}}k\sqrt{XY}.
\tag{7.4}
\]
Because \(\|B\|_{\rm op}<k\), equality at the upper endpoint is
possible only when \(X=Y=0\).  Thus the exact positive grounds of
\(A_0\) are the two-block-constant states and their common energy is
\[
P_0=k(k-1).
\tag{7.5}
\]
Also
\[
H_{A_0}(x,y)\ge-k-k\|B\|_{\rm op},
\tag{7.6}
\]
whose absolute value is strictly below \(P_0\) by (7.2).

Now flip one cross coefficient \(B_{i_0j_0}=a\).  Call the resulting
signing \(A\).  Every energy changes by at most two.  Hence
\[
M(A)=P_0+2,
\]
and its only absolute-ground pair is, up to global sign, the
two-block-constant state
\[
z=(\mathbf1,-a\mathbf1).
\tag{7.7}
\]
Indeed this state gains two, while the other two-block orientation
\[
z'=(\mathbf1,a\mathbf1)
\]
loses two and has slack exactly four:
\[
M(A)-|H_A(z')|=4.
\tag{7.8}
\]

Take the insertion row
\[
b=(\mathbf1,a\mathbf1).
\tag{7.9}
\]
It is exactly orthogonal to every absolute ground state,
\[
\max_{\{|H_A|=M(A)\}}|b\cdot x|=0,
\]
but
\[
|b\cdot z'|=2k=n.
\]
Therefore
\[
\boxed{\Delta_A(b)\ge n-4.}
\tag{7.10}
\]

This does not say that no *other* row works for this noncompetitive
example.  It proves the precise point needed in the insertion audit:
endpoint discrepancy of a chosen Euler/balancing row has no
deterministic thick-cap continuation.  A successful construction must
choose the row using the near-cap hierarchy itself, not only the exact
top and bottom faces or their ranks.

## 8. What remains open in the direct route

The direct route now has two separate missing estimates.

1. **Exact-face inverse.**  Outside the affine branch, prove a
   rank-versus-discrepancy theorem for the two opposite exact cap
   faces.  Low joint cap rank should yield a row with small endpoint
   fields; high rank must be converted, using the same-\(A\)
   four-point laws, into a larger current width.
2. **Thick-cap transfer.**  Even a row which balances every exact
   ground state must satisfy the slack-weighted inequalities (1.3).
   A successful theorem must use overlap geometry or a chained
   cap-intersection argument; counting the states in each energy layer
   separately is quantitatively insufficient.

Finally, any claimed one-step theorem must be checked against Section
2: \(o(\sqrt n)\) endpoint fields are useful input, but an
\(o(\sqrt n)\) total increment is impossible, and an
\(o(\sqrt n)\) derivative error does not alone imply convergence.
