# Switching-minimal graphs and the missing two-sided condition

## 1. Exact graph translation

Let

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad a_{ij}\in\{\pm1\}.
\]

Choose a sign of \(A\) and switch vertices so that

\[
H_A(\mathbf1)=M(A)=\max_x|H_A(x)|=:M.
\]

For a vertex set \(S\), let \(x^S\) be \(-1\) on \(S\) and \(+1\)
outside, and define the signed cut

\[
C_A(S)=\sum_{i\in S,\ j\notin S}a_{ij}.
\]

Then

\[
H_A(\mathbf1)-H_A(x^S)=2C_A(S).
\]

The absolute-ground-state choice gives the **two-sided** cut condition

\[
\boxed{
0\le C_A(S)\le M
\qquad(S\subseteq[n]).
}
\tag{1.1}
\]

Let \(G\) be the graph of negative edges.  Since

\[
C_A(S)=|S|(n-|S|)-2e_G(S,S^c),
\]

(1.1) is equivalent to

\[
\boxed{
\frac{|S|(n-|S|)-M}{2}
\le e_G(S,S^c)
\le\frac{|S|(n-|S|)}2.
}
\tag{1.2}
\]

Also

\[
\boxed{
M=\binom n2-2|E(G)|.
}
\tag{1.3}
\]

The upper inequality in (1.2) says exactly that \(G\) is
switching-minimal.  The lower inequality is equally important: it
records that the chosen positive ground state is an **absolute**
ground state.

## 2. One-sided switching-minimality cannot give an \(n^{3/2}\) bound

The graph literature usually calls \(G\) switching-minimal when only

\[
e_G(S,S^c)\le\frac{|S|(n-|S|)}2
\tag{2.1}
\]

is assumed.  This condition by itself permits only a linear total
deficit.

### Theorem 2.1

The maximum number of edges in a switching-minimal graph is

\[
\boxed{
\max |E(G)|=
\begin{cases}
\dfrac{n(n-2)}4,&n\ \text{even},\\[2mm]
\dfrac{(n-1)^2}4,&n\ \text{odd}.
\end{cases}
}
\tag{2.2}
\]

Equivalently, the minimum possible one-sided signed total is

\[
\boxed{
\binom n2-2|E(G)|=\left\lfloor\frac n2\right\rfloor.
}
\tag{2.3}
\]

#### Proof

First suppose \(n\) is even.  Apply (2.1) to a singleton.  Every
degree is at most \(\lfloor(n-1)/2\rfloor=(n-2)/2\), so

\[
|E(G)|\le\frac12n\frac{n-2}{2}.
\]

Equality is attained by

\[
G=K_{n/2}\ \dot\cup\ K_{n/2}.
\]

To verify switching-minimality directly, choose a balanced sign vector
\(u\) and put

\[
a_{ij}=-u_iu_j.
\]

Then

\[
H_A(x)=\frac{n-(u\cdot x)^2}{2}
\le\frac n2=H_A(\mathbf1),
\]

so every signed cut is nonnegative.  The negative edges are exactly
the two equal cliques.  \(\square\)

Now suppose \(n\) is odd.  Put

\[
r_i=n-1-2d_G(i).
\]

Singleton cuts give \(r_i\ge0\), and every \(r_i\) is even.  If
\(r_i=r_j=0\), the two-vertex cut inequality

\[
C_A(\{i,j\})=r_i+r_j-2a_{ij}\ge0
\]

forces \(a_{ij}=-1\), so all zero-row vertices form a clique.  Each
has degree \((n-1)/2\), hence there are at most \((n+1)/2\) of them.
Every nonzero \(r_i\) is at least two, and consequently

\[
\sum_i r_i\ge2\left(n-\frac{n+1}{2}\right)=n-1.
\]

Since \(\sum_i r_i=2(\binom n2-2|E(G)|)\), this proves the odd-order
upper bound in (2.2).

The rank-one construction with \(|\sum_i u_i|=1\) gives equality: its
negative graph is the union of cliques of sizes \((n-1)/2\) and
\((n+1)/2\), with

\[
H_A(\mathbf1)=\frac{n-1}{2}.
\]

Thus, in every order, one-sided switching-minimality permits

\[
\binom n2-2|E(G)|=\left\lfloor\frac n2\right\rfloor.
\tag{2.4}
\]

The construction fails the original problem maximally:

\[
H_A(u)=\frac{n-n^2}{2}=-\binom n2.
\]

Equivalently, its switching class contains a graph with nearly all
edges.  This is precisely what the missing lower inequality in (1.2)
forbids.

This distinction is consistent with the graph-theory literature.
Jelínek--Jelínková--Kratochvíl define switching-minimality by minimum
edge count in a switching class and record the singleton maximum-degree
bound; their work concerns decision complexity, not a two-sided
switching-class-width estimate:
https://arxiv.org/abs/1603.00254

## 3. Switching-class width is the exact remaining graph problem

Switching \(G\) by \(S\) changes its number of edges by

\[
|E(G^S)|-|E(G)|=C_A(S).
\tag{3.1}
\]

Hence (1.1) says that the entire switching class has edge counts in

\[
\boxed{
|E(G)|
\le |E(G^S)|
\le |E(G)|+M.
}
\tag{3.2}
\]

Since a uniformly random switch makes every edge present with
probability \(1/2\), the mean switching-class edge count is
\(\binom n2/2\).  Equation (1.3) centers (3.2) as

\[
\frac12\binom n2-\frac M2
\le |E(G^S)|
\le
\frac12\binom n2+\frac M2.
\tag{3.3}
\]

Thus the original min--max problem is exactly:

> How narrow can the edge-count distribution of a Seidel switching
> class be?

Ordinary switching-minimality controls only the left endpoint.  The
rank-one example proves that no theorem using only that endpoint can
produce even \(\omega(n)\), much less the conjectural
\(\frac12n^{3/2}\).

## 4. Exact vertex-deletion repair lemma

Although it does not close the asymptotics, switching-minimality has a
useful exact recursion.

Let \(G\) be switching-minimal and fix a vertex \(v\).  Put

\[
H=G-v,\qquad
r_v=n-1-2d_G(v)\ge0.
\]

Choose \(T\subseteq V(H)\) which minimizes the number of edges in the
switching class of \(H\), and write

\[
\Delta=|E(H)|-|E(H^T)|\ge0.
\]

Switching \(T\) in the full graph changes the incident edges at \(v\)
by

\[
\alpha=|T|-2e_G(v,T).
\]

Switching \(V(H)\setminus T\) performs the same internal switch and
changes the incident edges by \(r_v-\alpha\).  Full
switching-minimality gives

\[
\alpha-\Delta\ge0,
\qquad
r_v-\alpha-\Delta\ge0.
\]

Therefore

\[
\boxed{
0\le\Delta\le\frac{r_v}{2}.
}
\tag{4.1}
\]

In signed-total form, the switching-minimal representative of the
induced class has

\[
\boxed{
M-r_v
\le
\left[\binom{n-1}{2}-2|E(H^T)|\right]
=M-r_v+2\Delta
\le M.
}
\tag{4.2}
\]

This is a clean deletion/repair theorem: deleting \(v\) loses the row
sum \(r_v\), and re-minimizing the induced switching class can replenish
at most that same amount.

The obstruction is also exact.  Equation (4.1) supplies no positive
fractional loss; \(\Delta=r_v/2\) is allowed.  Consequently scalar
iteration gives monotonicity but not an \(O(\sqrt n)\) increment or a
scale-contraction theorem.

## 5. Multi-cut and quotient identities

Let \(\mathcal P=(V_1,\ldots,V_m)\) be a vertex partition and put

\[
w_{ab}=\sum_{i\in V_a,\ j\in V_b}a_{ij}
\qquad(a<b).
\]

For every \(I\subseteq[m]\), the union cut obeys

\[
\boxed{
0\le
\sum_{a\in I,\ b\notin I}w_{ab}
\le M.
}
\tag{5.1}
\]

Thus every contracted weighted quotient inherits the same two-sided
cut cap.  This is an exact structural recursion, but its edge weights
are arbitrary integers rather than flat signs; the flatness rigidity
is lost under contraction.

For two cuts \(S,T\), let \(D(S,T)\) be the signed sum over edges which
cross both cuts.  The cut parity identity

\[
\delta(S)\triangle\delta(T)=\delta(S\triangle T)
\]

gives

\[
\boxed{
D(S,T)
=
\frac{
C_A(S)+C_A(T)-C_A(S\triangle T)
}{2}.
}
\tag{5.2}
\]

Hence

\[
\boxed{
-\frac M2\le D(S,T)\le M.
}
\tag{5.3}
\]

In the four-cell partition generated by \(S,T\), (5.2) controls the
sum of the two opposite inter-cell blocks.  Cyclically permuting the
three cuts controls the other pairings.  These are the exact
triangle/multi-cut constraints beyond one-cut positivity.

They do not by themselves improve the scale: a fixed number of cuts
produces a fixed-dimensional weighted quotient, while the rank-one
one-sided obstruction is eliminated only by the full upper cap over
all \(2^{n-1}\) cuts.

## 6. Cut-cone/localizing hierarchy

Use the doubled quadratic normalization

\[
q_A(x)=x^\top Ax,\qquad Q=q_A(\mathbf1)=2M.
\]

The full condition is

\[
0\le Q-q_A(x)\le2Q
\qquad(x\in\{\pm1\}^n).
\tag{6.1}
\]

For every Boolean polynomial \(p\),

\[
0\le
\mathbb E[(Q-q_A(x))p(x)^2]
\le
2Q\,\mathbb E[p(x)^2].
\tag{6.2}
\]

At degree one, \(p(x)=u\cdot x\), this gives

\[
\boxed{
-\frac Q2I\preceq A\preceq\frac Q2I.
}
\tag{6.3}
\]

This recovers only the spectral \(M\ge\sqrt{n-1}\) scale.  Higher
degree produces the localizing/SOS hierarchy for the nonnegative cut
polynomial.  Degree two introduces the signed triangle operator on
edge characters; fixed degree still sees only finitely many
switching-invariant cycle statistics.  Reaching the \(n^{3/2}\)
extreme-value scale requires degree growing with \(n\), where the
hierarchy becomes the original covering-radius problem.

This explains why ordinary cut-cone duality cannot solve the problem:
after relaxing \(a_{ij}\in\{\pm1\}\) to an interval, the flat
integrality effect disappears.  The \(n^{3/2}\) deficit is not a
consequence of one-cut linear inequalities.

## 7. Verdict

The switching-minimal reformulation is useful only after retaining its
two-sided cap.

New rigorous conclusions from this audit are:

1. one-sided switching-minimal graphs have only a linear forced
   deficit, with the even-order extremum (2.2);
2. the original problem is the minimum possible width of the
   edge-count distribution of a switching class;
3. vertex deletion has the exact replenishment bound (4.1)--(4.2);
4. two cuts give the exact four-cell constraint (5.2)--(5.3).

No improved \(n^{3/2}\) constant follows.  Any successful graph
argument must use the upper and lower cut inequalities simultaneously
at a number of correlated cuts growing with \(n\).  A theorem stated
only for switching-minimal graphs cannot address the original limit.
