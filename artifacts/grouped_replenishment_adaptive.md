# Adaptive ground-state closure for replenishment gaps

Checkpoint date: 2026-07-25.

## 1. Verdict

The false pointwise inverse between a replenishment gap and successor
ground-layer visibility is not needed if the deletion blocks may be
chosen adaptively.

After an arbitrary initial restriction, repeatedly delete the smaller
disagreement set between the inherited state and a one-sided ground
state of the current core.  Along this **ground-state closure chain**,
the next replenishment gap is at most one half of the current gap:

\[
\boxed{
g_{j+1}\le \frac12g_j.
}
\tag{1.1}
\]

Here \(P(C)=\max_x x^\top Cx\) is the one-sided positive maximum in one
fixed orientation.  Consequently

\[
\boxed{
\sum_jg_j\le2g_1
\le2\bigl(P(A)+N(A)\bigr)
\le4Q(A),
}
\tag{1.2}
\]

where \(N(A)=-\min_xx^\top Ax\).  If the inherited core energy is
nonnegative in this orientation, the right side improves to \(2Q(A)\).
This is an exact \(O(Q)\) grouped replenishment theorem.  It does not
bound the sum of gaps along a prescribed singleton order; the order and
block sizes are changed adaptively.

## 2. One-sided notation and monotonicity

For a symmetric zero-diagonal matrix \(C\), put

\[
E_C(x)=x^\top Cx,\qquad
P(C)=\max_{x\in\{\pm1\}^{V(C)}}E_C(x).
\]

The mean of \(E_C(x)\) over the cube is zero, so \(P(C)\ge0\).
If \(C'\) is a principal submatrix of \(C\), then

\[
P(C')\le P(C). \tag{2.1}
\]

Indeed, take a positive maximizer on \(C'\), independently complete
the deleted coordinates with mean-zero signs, and average.  The
expected full energy is the energy on \(C'\), so one completion has at
least that energy.

For disjoint vertex blocks \(D_1,\ldots,D_r\) inside a common ambient
matrix \(A\),

\[
\boxed{
\sum_{j=1}^rP(D_j)\le P(A).
}
\tag{2.2}
\]

Choose on each block a vector attaining \(P(D_j)\).  Independently
multiply the entire vector on every block by a random sign, and assign
independent random signs to any unassigned vertices.  All inter-block
and unassigned-edge terms have mean zero, while the internal
contribution has mean \(\sum_jP(D_j)\).  Some realization therefore
has full energy at least this sum.

More strongly, \(P\) is superadditive across every two-block
partition:

\[
\boxed{
P\!\begin{pmatrix}D&B\\B^\top&C\end{pmatrix}
\ge P(D)+P(C).
}
\tag{2.3}
\]

Choose positive ground states \(u,v\) of \(D,C\).  The larger of the
energies of \((u,v)\) and \((-u,v)\) is

\[
P(D)+P(C)+2|u^\top Bv|.
\]

The same argument after reversing the matrix orientation gives

\[
N\!\begin{pmatrix}D&B\\B^\top&C\end{pmatrix}
\ge N(D)+N(C). \tag{2.4}
\]

## 3. Closure construction

Start with a core \(C_1\) and an inherited Boolean state \(y_1\).
Fix one orientation of the matrix for the whole construction.  Let

\[
p_j=P(C_j),\qquad
g_j=p_j-E_{C_j}(y_j)\ge0.
\]

Choose \(z_j\) with \(E_{C_j}(z_j)=p_j\).  Since replacing \(z_j\) by
\(-z_j\) does not change its energy, choose its global sign so that
the disagreement set

\[
S_j=\{v:y_j(v)\ne z_j(v)\}
\]

has size at most \(|V(C_j)|/2\).  If \(g_j>0\), then \(S_j\) is
nonempty.  Delete \(S_j\), put

\[
C_{j+1}=C_j[V(C_j)\setminus S_j],
\qquad
y_{j+1}=y_j|_{V(C_{j+1})}.
\]

By construction,

\[
y_{j+1}=z_j|_{V(C_{j+1})}. \tag{3.1}
\]

The procedure terminates after finitely many steps, either with
\(g_j=0\) or with a core of order at most one.

This is a legitimate adaptive peeling process: at every step the
state inherited by the new core is exactly the restriction of a
positive ground state of the preceding core.

## 4. Exact half-contraction identity

Write

\[
C_j=
\begin{pmatrix}
D_j&B_j\\
B_j^\top&C_{j+1}
\end{pmatrix}
\]

according to \(S_j\) and its complement.  Since \(y_j=-z_j\) on
\(S_j\) and \(y_j=z_j\) on its complement, their energies inside
\(D_j\) agree.  Set

\[
a_j=z_{j,S_j}^\top D_jz_{j,S_j}
    =y_{j,S_j}^\top D_jy_{j,S_j}.
\]

If

\[
c_j=z_{j,S_j}^\top B_jz_{j,V(C_{j+1})},
\]

then the cross term changes sign between \(z_j\) and \(y_j\).  Hence

\[
g_j
=E_{C_j}(z_j)-E_{C_j}(y_j)
=4c_j. \tag{4.1}
\]

Using (3.1) and

\[
E_{C_{j+1}}(y_{j+1})=p_{j+1}-g_{j+1},
\]

the ground-state energy has the exact decomposition

\[
p_j
=a_j+p_{j+1}-g_{j+1}+\frac12g_j.
\]

Therefore

\[
\boxed{
g_{j+1}
=\frac12g_j+a_j-(p_j-p_{j+1}).
}
\tag{4.2}
\]

By (2.3),

\[
p_j-p_{j+1}\ge P(D_j),
\]

while \(a_j\le P(D_j)\).  Substitution in (4.2) proves the pure
contraction (1.1).

The identity is sharper than the earlier energy-layer inequality

\[
\frac12g_j\le
\bigl(p_j-p_{j+1}\bigr)+Q(D_j)+g_{j+1},
\]

which points in the wrong direction for telescoping.  Retaining the
signed internal energy and the full one-sided decrement reveals that
the positive induced-block term cancels completely.

## 5. Telescoping

Let the nonzero gaps be \(g_1,\ldots,g_L\), and set \(g_{L+1}=0\).
The pure contraction gives immediately

\[
\boxed{
\sum_{j=1}^{L}g_j
\le g_1\sum_{r=0}^{L-1}2^{-r}
<2g_1.
}
\tag{5.1}
\]

For reference, summing the exact identity (4.2) also gives

\[
\sum_{j=1}^{L}g_j
=2g_1+2\sum_{j=1}^{L}a_j
-2(p_1-p_{L+1}). \tag{5.2}
\]

Since \(g_1=p_1-E_{C_1}(y_1)\), this is

\[
\sum_{j=1}^{L}g_j
=2\left(
\sum_{j=1}^{L}a_j+p_{L+1}-E_{C_1}(y_1)
\right). \tag{5.3}
\]

The blocks \(D_1,\ldots,D_L\) and the terminal core \(C_{L+1}\)
are disjoint.  Applying (2.2) to all of them gives

\[
\sum_ja_j+p_{L+1}
\le
\sum_jP(D_j)+P(C_{L+1})
\le P(A).
\]

For the universal bound, one-sided monotonicity gives
\[
g_1
=P(C_1)-E_{C_1}(y_1)
\le P(A)+N(A),
\]
which proves (1.2).  If \(E_{C_1}(y_1)\ge0\), then
\(g_1\le P(C_1)\le P(A)\), and the cumulative bound improves to
\(2Q(A)\).

## 6. Orientation-free initial deficit

The block inequalities retain

\[
|E_{C_1}(y_1)|
\]

and lose the absolute deficit

\[
\bar g_1
=Q(C_1)-|E_{C_1}(y_1)|.
\]

Choose a sign \(s\in\{\pm1\}\) such that
\(P(sC_1)=Q(C_1)\), and orient the entire ambient matrix by the same
sign.  Then the one-sided starting gap is

\[
g_1
=Q(C_1)-sE_{C_1}(y_1)
\ge Q(C_1)-|E_{C_1}(y_1)|
=\bar g_1. \tag{6.1}
\]

All later steps use the fixed orientation \(sA\) and its positive
maxima.  Thus (1.2) controls a one-sided closure whose initial gap
dominates the absolute deficit.  It does **not** assert that every
later absolute deficit is dominated by the corresponding one-sided
gap; the correct way to plug the theorem into the block recursion is
to use the one-sided block inequalities in the next section.

## 7. One-sided block inequalities

The absolute-value block lemmas have exact one-sided versions.  Write

\[
A=\begin{pmatrix}D&B\\B^\top&C\end{pmatrix}
\]

and fix \(y\in\{\pm1\}^{V(C)}\).  Put

\[
e=E_C(y),\qquad
L=\|By\|_2,\qquad
H=\|By\|_1,\qquad
N(D)=-\min_zE_D(z).
\]

Then

\[
\boxed{
P(A)\ge e+\sqrt2\,L
}
\tag{7.1}
\]

and

\[
\boxed{
P(A)\ge e+\Phi(H,N(D)).
}
\tag{7.2}
\]

For (7.1), average a uniform \(z\) on the deleted block and compare
\((z,y)\) with \((-z,y)\).  Their larger (not absolute) energy is
\[
E_D(z)+e+2|z^\top By|.
\]
The internal term has mean zero, and sharp Khintchine gives (7.1).

For (7.2), independently round \(Z\) with
\[
\mathbb EZ=t\,\operatorname{sign}(By).
\]
Then
\[
\mathbb EE_A(Z,y)
\ge e+2tH-t^2N(D),
\]
and optimization in \(t\in[0,1]\) gives (7.2).

These formulas avoid every orientation problem: orient the original
matrix so that \(P(A)=Q(A)\), and keep that orientation throughout the
whole deletion tree.

There is also a sharper one-sided cumulative visibility theorem.  If
\[
V_j^+=\max_{y:E_{C_{j+1}}(y)=P(C_{j+1})}\|B_jy\|_1
\]
and \(d_j=P(C_j)-P(C_{j+1})\), then
\[
\boxed{
2V_j^+\le d_j+N(D_j).
}
\tag{7.3}
\]
For disjoint deleted blocks,
\[
\boxed{
\sum_jV_j^+
\le\frac{P(A)+N(A)}2
\le Q(A).
}
\tag{7.4}
\]
Thus the earlier \(3Q/2\) absolute-ground visibility bound improves to
\(Q\) when the recursion is kept in one fixed orientation.

## 8. Exact plug into the block recursion

At a closure step, (4.1) implies
\[
\|B_jy_{j+1}\|_1\ge \frac14g_j.
\]
Define the harvested one-sided gain
\[
\Gamma_j
=\Phi\left(\frac14g_j,N(D_j)\right).
\]
Equations (7.2) and (4.2) give
\[
p_j\ge p_{j+1}+\Gamma_j-g_{j+1}. \tag{8.1}
\]

There is an exact local payment inequality.  For every \(b,d\ge0\),
\[
\Phi(b/4,d)+d\ge b/2. \tag{8.2}
\]
When \(b/4\le d\), the difference between the two sides is
\[
\frac{(d-b/4)^2}{d};
\]
when \(b/4\ge d\), equality holds.  Combining (8.2) with
\[
g_{j+1}\le\frac12g_j
\]
gives
\[
\boxed{
g_{j+1}
\le
\Gamma_j+N(D_j).
}
\tag{8.3}
\]

Consequently, when inequalities (8.1) are telescoped, every
replenishment gap after the first is charged exactly once:
\[
\boxed{
\Gamma_j-g_{j+1}
\ge-N(D_j).
}
\tag{8.4}
\]
The total remaining charge is bounded by
\[
\sum_jN(D_j)
\le N(A)
\le Q(A). \tag{8.5}
\]

More explicitly, suppose an initial selected block gives
\[
P(A)\ge p_1+\Gamma_0-g_1.
\]
Following it by closure steps and a terminal core \(R\) gives
\[
\boxed{
P(A)
\ge
P(R)+\Gamma_0-g_1
-\sum_jN(D_j).
}
\tag{8.6}
\]
Positive unused \(\Gamma_j\) terms can be retained, so (8.6) is the
worst-case form.

## 9. The precise \(3/2\)-scale barrier

Equation (8.6) removes the unbounded adaptivity gap, but its present
coefficient does not prove scale-preserving amplification.  The negative-norm
charge in (8.5) can itself be of leading order.

If \(s_j=|D_j|\), then the elementary bound
\[
N(D_j)\le s_j(s_j-1)
\]
gives the explicit scale-sensitive form
\[
\sum_jN(D_j)
\le
\sum_js_j^2. \tag{9.1}
\]
Therefore the closure is asymptotically lossless at the \(n^{3/2}\)
scale under the concrete condition
\[
\boxed{
\sum_js_j^2=o(n^{3/2}).
}
\tag{9.2}
\]
In particular, (9.2) holds if every closure block has
\(o(\sqrt n)\) vertices.  A dyadic block-size class
\(s\le |D_j|<2s\) contributes at most \(O(ns)\), because the blocks
are disjoint.

The missing theorem is now sharply geometric:

> A mesoscopic initial deletion must admit a ground-state closure
> satisfying (9.2), or a macroscopic disagreement block must itself
> yield an additional quadratic witness that pays its negative norm.

The first alternative is not automatic.  In the certified order-15
replenishment example, deleting one vertex can produce a first
one-sided gap \(20\) whose closest successor ground state differs on
five of the remaining fourteen vertices.  Thus a microscopic deletion
can trigger a macroscopic closure block even when the favorable
orientation is chosen.

This is the surviving coefficient barrier.  The replenishment losses
contract geometrically, but controlling the **sizes or negative norms**
of the adaptive disagreement blocks is still required when the
generic gains \(\Gamma_j\) are themselves telescoped in a
\(3/2\)-scale deletion inequality.

## 10. What this proves, and what it does not

### Proved

* Every initial restriction admits an adaptive sequence of
  disagreement-block deletions with exact geometric contraction
  \(g_{j+1}\le g_j/2\).  Its cumulative one-sided replenishment is at
  most \(2g_1\), hence at most \(4Q(A)\), improving to \(2Q(A)\) when
  its inherited energy is nonnegative.
* The contraction is exact and geometric; no entropy hypothesis,
  spectral regularity, or ground-layer visibility inverse is used.
* One-sided block extraction and cumulative visibility avoid all
  orientation changes, and (8.4) charges each later gap exactly once
  using only the negative norm of its disagreement block.
* A scalable universal-row tower cannot refute this adaptive grouped
  statement.  Such a tower only concerns a prescribed singleton
  deletion history, whereas closure changes the subsequent blocks.

### Still needed for the main limit problem

The exact generic-gain recursion (8.6) still pays the negative norms
of the closure blocks.  They are globally \(O(Q(A))\), but this is not
\(o(n^{3/2})\) for a near-minimizing sequence.  Closing the main limit
problem by this route therefore requires either (9.2), a strict
coefficient recovery from the unused \(\Gamma_j\), or a new witness
forced by a macroscopic two-replica disagreement block.
