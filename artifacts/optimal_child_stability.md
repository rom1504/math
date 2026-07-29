# Optimal-child stability under one-vertex puncturing

## Status

This note attacks the defect left open by the near-cap insertion
identity.  If \(A\) is a global order-\(n\) minimizer, \(B_i=A[-i]\),
and
\[
d_i=M_n-M(B_i),\qquad e_i=M(B_i)-M_{n-1},
\]
then
\[
e_i+d_i=M_n-M_{n-1}.
\]
Controlling \(d_i\) alone is therefore insufficient.

The first result below is a certified obstruction at the exact
arithmetic orders \(n\equiv1\pmod4\): even there a global minimizer
need not have any optimal child.  At order nine the defect is
simultaneously positive at every vertex, while every deletion cost
vanishes.  This kills optimal-puncture and average-puncture arguments
in their strongest forms.  It leaves a sharper possible target:
whether some child always has \(e_i\le2\), equivalently whether a
puncture of a deepest augmented-cut-code hole is within one unit of
the child covering radius.

## 1. Exact order-nine obstruction

Let
\[
A=\begin{pmatrix}
0&1&-1&1&1&1&1&-1&-1\\
1&0&1&-1&-1&-1&1&-1&-1\\
-1&1&0&1&1&1&1&1&-1\\
1&-1&1&0&1&1&1&-1&1\\
1&-1&1&1&0&1&-1&1&-1\\
1&-1&1&1&1&0&-1&-1&-1\\
1&1&1&1&-1&-1&0&1&1\\
-1&-1&1&-1&1&-1&1&0&-1\\
-1&-1&-1&1&-1&-1&1&-1&0
\end{pmatrix}.
\tag{1.1}
\]
Direct enumeration of the \(2^8\) spins modulo global negation gives
\[
M(A)=12.
\tag{1.2}
\]
For each \(i\in[9]\), an independent enumeration of the \(2^7\)
spins of the principal child gives
\[
\boxed{M(A[-i])=12\quad\hbox{for all }i.}
\tag{1.3}
\]
The certified exact values \(M_9=12\) and \(M_8=10\) therefore imply
that \(A\) is a global minimizer and
\[
\boxed{d_i=0,\qquad e_i=2\quad(i=1,\ldots,9).}
\tag{1.4}
\]

The calculation was performed by two independent evaluators:
`optimal_child_stability_verify.cpp`, which evaluates the displayed
matrix from scratch, and a direct Python evaluator reconstructed from
the same matrix.  Both return parent norm \(12\) and the child-norm vector
\[
(12,12,12,12,12,12,12,12,12).
\]

### Consequences

1. The statement “every order \(4k+1\) minimizer has an optimal
   order-\(4k\) principal child” is false.
2. Even the weaker statement “some deletion realizing the best
   near-cap insertion cost has an optimal child” is false: here every
   deletion has the best possible cost \(d_i=0\), but no child is
   optimal.
3. Averaging does not repair the issue:
   \[
   \frac1n\sum_i d_i=0,\qquad
   \frac1n\sum_i e_i=2.
   \]
4. Since \(d_i=0\), a child absolute ground has zero field at the
   deleted vertex.  Its two extensions are parent absolute grounds
   differing only at coordinate \(i\).  Thus every coordinate
   direction occurs in the graph induced by the parent absolute
   ground set.  Dense exact-face coverage is therefore compatible
   with a uniformly nonoptimal puncture profile.

The obstruction is stronger than a parity artifact: it occurs at
\(n=9\equiv1\pmod4\), exactly where edge-flip replacement witnesses
upgrade to exact face witnesses.

## 2. Covering-radius translation

Let \(N_n=\binom n2\), let \(\mathcal C_n\) be the augmented cut code,
and put
\[
\rho_n=\rho(\mathcal C_n)=\frac{N_n-M_n}{2}.
\]
For an order-\(n\) signing \(A\), let
\[
r_i=\operatorname{dist}(A[-i],\mathcal C_{n-1})
    =\frac{N_{n-1}-M(A[-i])}{2}.
\]
Then the child defect is exactly
\[
\boxed{e_i=2(\rho_{n-1}-r_i).}
\tag{2.1}
\]
Thus (1.4) says that a deepest order-nine hole can have every
puncture exactly one unit below the order-eight covering radius.

The next sharp finite target is consequently
\[
\boxed{
\max_{A:\,\operatorname{dist}(A,\mathcal C_n)=\rho_n}
\ \min_i\bigl(\rho_{n-1}-r_i\bigr)\le1
}
\tag{2.2}
\]
for \(n\equiv1\pmod4\), or at least an asymptotic/summable-error
version after choosing a convenient deepest hole.  In energy
language this is \(\min_i e_i\le2\).

## 3. A basic extension inequality and why it does not prove (2.2)

Take a closest augmented cut codeword to the puncture \(A[-i]\).
After restoring vertex \(i\), choose its sign to minimize the number
of star disagreements.  Since the restored star has \(n-1\) edges,
\[
\rho_n\le r_i+\left\lfloor\frac{n-1}{2}\right\rfloor.
\tag{3.1}
\]
For odd \(n\), if every puncture were at least \(k\) below the child
covering radius, this would only imply
\[
M_n-M_{n-1}\ge2k.
\tag{3.2}
\]
Hence the general extension inequality rules out radius defect
\(\ge k\) only when the scalar increment is \(<2k\).  At the natural
\(\Theta(\sqrt n)\) increment scale it gives no constant puncture
stability.  Any proof of (2.2) must use compatibility among the
different child approximants, global edge-flip optimality, or special
triangle relations of the augmented cut code.

## 4. Current frontier

The exact order-nine obstruction changes the useful question from
“is there an optimal child?” to one of the following.

1. **Radius-one puncture:** prove (2.2), perhaps only for a
   lexicographically chosen deepest hole.
2. **Exchange component:** prove that every optimal order-\(n\)
   switching class is connected, through norm-preserving batch
   exchanges, to a minimizer with an \(e_i\le2\) child.
3. **Asymptotic stability:** prove
   \[
   \min_i e_i=O(1)
   \quad\text{or}\quad
   \sum_n\frac{\min_i e_i}{n^{3/2}}<\infty
   \]
   for a compatible sequence of global minimizers.

The first statement is the cleanest and is not falsified by any exact
value through order ten.  The order-nine example shows that its
constant \(2\) would be sharp.

## 5. Exact local-to-global decoding identity

The natural agreement approach can be made completely exact.  It
also exposes why ordinary agreement estimates do not yet see a
constant child defect.

Let \(A\) be a deepest order-\(n\) hole.  For each \(i\), choose a
closest augmented cut codeword
\[
c_i\in\mathcal C_{n-1}
\quad\hbox{to}\quad A[-i],
\qquad
d(A[-i],c_i)=r_i.
\]
Define the global repair cost of this family by
\[
\mathcal R(c_1,\ldots,c_n)
=
\min_{c\in\mathcal C_n}
\sum_{i=1}^n d(c[-i],c_i).
\tag{5.1}
\]
Then
\[
\boxed{
(n-2)\rho_n
\le
\sum_i r_i+\mathcal R(c_1,\ldots,c_n).
}
\tag{5.2}
\]
Indeed, for every global \(c\in\mathcal C_n\),
\[
d(A[-i],c[-i])
\le r_i+d(c_i,c[-i]).
\]
After summing over \(i\), every global edge disagreement is counted
exactly \(n-2\) times.  Minimize the repair term and use that every
global codeword is at distance at least \(\rho_n\) from the deep
hole \(A\).

Substituting
\[
r_i=\rho_{n-1}-\frac{e_i}{2}
\]
and cancelling the edge-count terms gives the energy form
\[
\boxed{
(n-2)M_n
\ge
nM_{n-1}+\sum_i e_i-2\mathcal R.
}
\tag{5.3}
\]
In particular, if the closest child codewords are exactly compatible
and hence are the punctures of one global codeword, then
\(\mathcal R=0\) and
\[
(n-2)(M_n-M_{n-1})
\ge
2M_{n-1}+\sum_i e_i.
\tag{5.4}
\]
Thus exact agreement forces a scalar increment of coefficient at
least \(2\), not the scale-correct coefficient \(3/2\).  The child
approximants of an asymptotically smooth minimizing sequence must
therefore be genuinely incompatible.  For the order-nine matrix in
Section 1, exact compatibility would give
\[
7(12-10)\ge2(10)+9(2),
\]
which is false; its nine closest child codewords necessarily fail to
glue.

### 5.1 Robust agreement repair

There is also a quantitative version.  Write
\[
\Gamma
=
\sum_{i<j}
d\!\left(
c_i\big|_{E([n]\setminus\{i,j\})},
c_j\big|_{E([n]\setminus\{i,j\})}
\right).
\tag{5.5}
\]
For \(n\ge12\), a direct majority decoder and the triangle test for the
augmented cut code give an absolute constant \(C\) such that
\[
\boxed{\mathcal R\le C\,\Gamma/n.}
\tag{5.6}
\]
One may take the nonoptimized constant \(C=100\).

Here is a self-contained proof.  Represent an augmented cut codeword
as
\[
c_{uv}=t+z_u+z_v\pmod2.
\]
First retain the majority value of \(t\) among the \(c_i\)'s.  If
two local words have opposite \(t\), their difference on
\(r=n-2\) common vertices is a complemented cut and has weight at
least
\[
D_r=\binom r2-\left\lfloor r^2/4\right\rfloor.
\tag{5.7}
\]
Thus the total cost of discarding the minority-\(t\) local words,
even at the trivial cost \(\binom{n-1}2\) per word, is
\(O(\Gamma/n)\).

For the remaining same-\(t\) words, take the edgewise majority \(w\)
of their local predictions.  If \(L\) is the total number of local
edge disagreements with \(w\), pairwise majority counting gives
\[
L\le 6\Gamma/n.
\tag{5.8}
\]
Every triangle on which \(w\) has parity different from \(t\) must
contain a local edge disagreement in every local view containing that
triangle.  Double counting triangle--view incidences shows that the
number of bad triangles is at most \(4L\).  Choose a root contained
in at most \(12L/n\) bad triangles and set
\[
z_u=w_{ru}+t.
\]
The global word \(c_{uv}=t+z_u+z_v\) then differs from \(w\) on at
most \(12L/n\) edges.  Summing its discrepancies over all local views
costs at most
\[
L+(n-2)\frac{12L}{n}\le13L.
\]
Together with the opposite-\(t\) estimate this proves (5.6), with
room in the displayed constant.

Combining (5.3) and (5.6) yields the exact dichotomy
\[
\boxed{
(n-2)(M_n-M_{n-1})
\ge
2M_{n-1}+\sum_i e_i-\frac{2C}{n}\Gamma.
}
\tag{5.9}
\]
Consequently, if the scalar increment is near the \(3/2\)-homogeneous
derivative and all child defects are at least \(4\), then necessarily
\[
\Gamma=\Omega(n^{5/2}).
\tag{5.10}
\]
Because two same-\(t\) local augmented cuts either agree or differ on
at least \(n-3\) overlap edges, (5.10) forces at least
\(\Omega(n^{3/2})\) incompatible pairs, unless a substantial part of
the mass comes from opposite-\(t\) pairs.

This is a rigorous local-to-global reduction, but also a scale wall.
The difference between \(e_i=2\) and \(e_i=4\) contributes only
\(2n\) to (5.9), whereas the unavoidable leading repair budget is of
order \(n^{3/2}\).  A generic agreement theorem, without a
leading-term cancellation specific to global minimizers, cannot
resolve the radius-one question.

## 6. Exact equivalences and the weakest useful asymptotic form

Let
\[
\delta_n=M_n-M_{n-1}.
\]
For a deepest parent,
\[
e_i+d_i=\delta_n,
\]
so radius-one puncturing is exactly the max-deletion statement
\[
\boxed{
\min_i e_i\le2
\iff
\max_i d_i\ge\delta_n-2.
}
\tag{6.1}
\]
It is **not** equivalent to the stronger residue-class increment bound
\(\delta_n\le2\).  The latter implies (6.1) trivially, but (6.1)
allows \(d_i\) and \(\delta_n\) to be of order \(\sqrt n\).

The minimum in (6.1), by itself, does **not** combine with a
small-\(d_i\) deletion chosen by the near-cap argument.  In fact
\[
e_i=\delta_n-d_i
\tag{6.2}
\]
is perfect anticorrelation: the child with smallest defect is exactly
the child with largest deletion cost.  Separate estimates on
\(\min_i e_i\) and \(\min_i d_i\) need not concern the same vertex.

For the convergence campaign, suppose first that the near-cap argument
supplies vertices \(i(n)\) with
\[
d_{i(n)}
\le\frac{3M_n}{2n}+O(1).
\tag{6.3}
\]
The exact joint child-stability input needed by that one-step proof is
\[
\boxed{
\sum_n\frac{e_{i(n)}}{n^{3/2}}<\infty.
}
\tag{6.4}
\]
Indeed
\[
M_n-M_{n-1}
\le\frac{3M_n}{2n}+e_{i(n)}+O(1),
\]
and the normalized positive variation is summable.  Uniform
\(e_i=O(1)\) suffices, as does
\[
e_i=O\!\left(\frac{\sqrt n}{(\log n)^{1+\varepsilon}}\right).
\]
In contrast, an unquantified \(e_i=o(\sqrt n)\) is not enough: after
division by \(n^{3/2}\) it can still produce a divergent harmonic-type
error series.

There are three stronger but cleaner sufficient formulations.

1. If the near-cap theorem gives only the average estimate
   \[
   \overline d_n=\frac1n\sum_i d_i
   \le\frac{3M_n}{2n}+O(1),
   \]
   then it is enough to prove
   \[
   \sum_n\frac{\overline e_n}{n^{3/2}}<\infty,
   \qquad
   \overline e_n=\frac1n\sum_i e_i,
   \tag{6.5}
   \]
   because \(\delta_n=\overline d_n+\overline e_n\).
   In particular, the mean-radius-one theorem
   \[
   \boxed{\frac1n\sum_i e_i\le2}
   \tag{6.6}
   \]
   would close the recurrence.  In covering language, (6.6) is
   \[
   \frac1n\sum_i r_i\ge\rho_{n-1}-1.
   \]

2. If only one vertex with (6.3) is obtained, the uniform theorem
   \[
   \max_i e_i=O(1)
   \tag{6.7}
   \]
   suffices.

3. Radius-one puncturing (6.1) becomes sufficient when supplemented
   by bounded profile spread:
   \[
   \max_i d_i-\min_i d_i=O(1).
   \tag{6.8}
   \]
   Indeed the \(e\)-profile has exactly the same range, in reverse
   order, so (6.1) and (6.8) imply (6.7).

Thus a compatible-class replacement of radius-one puncturing may
target (6.4) or (6.5): for each order, choose a global minimizer and a
high-traffic deletion jointly so that the resulting child defects
form a summable sequence.  This is strictly weaker than asking every
deep hole, or even every order, to have a radius-one puncture.

## 7. A failed residue-class shortcut

A tempting strengthening was:

> Every order \(m\equiv0\pmod4\) signing \(B\) has an insertion row
> whose affine norm is at most \(M(B)+2\).

This would imply \(M_{4k+1}\le M_{4k}+2\), and hence radius-one
puncturing trivially.  It is false.  The signing
\[
B=\begin{pmatrix}
0&1&1&1&-1&-1&1&-1&1&-1&1&-1\\
1&0&-1&1&-1&-1&1&1&1&1&-1&1\\
1&-1&0&-1&1&1&1&-1&1&-1&-1&1\\
1&1&-1&0&1&1&-1&1&1&-1&-1&-1\\
-1&-1&1&1&0&1&1&1&-1&-1&-1&-1\\
-1&-1&1&1&1&0&-1&-1&1&1&-1&-1\\
1&1&1&-1&1&-1&0&1&-1&1&-1&-1\\
-1&1&-1&1&1&-1&1&0&-1&-1&-1&1\\
1&1&1&1&-1&1&-1&-1&0&1&-1&1\\
-1&1&-1&-1&-1&1&1&-1&1&0&-1&-1\\
1&-1&-1&-1&-1&-1&-1&-1&-1&-1&0&-1\\
-1&1&1&-1&-1&-1&-1&1&1&-1&-1&0
\end{pmatrix}
\tag{7.1}
\]
has
\[
M(B)=18,
\]
whereas exhaustive enumeration of all \(2^{11}\) insertion rows and
all \(2^{11}\) core spins gives
\[
\boxed{
\min_{b\in\{\pm1\}^{12}}
\max_x\bigl(|H_B(x)|+|b\cdot x|\bigr)=24.
}
\tag{7.2}
\]
Thus the fixed-core insertion gap is \(6\), even at a multiple-of-four
order.  This calculation is implemented deterministically in
`optimal_child_stability_verify.cpp`; the displayed matrix was
independently checked by a vectorized exact evaluator.

It is not currently known whether (7.1) is a global order-twelve
minimizer.  Therefore (7.2) does not disprove the residue-class bound
for a *specially chosen compatible minimizer*.  It does prove that
any such argument must exploit global optimality and the choice of
switching class; parity alone, or a universal completion lemma, cannot
work.
