# Bent/parity nonconvergence audit

## 0. Question and normalization

Let

\[
F(n)=
\min_{a_{ij}=\pm1}
\max_{x\in\{\pm1\}^n}
\left|
\sum_{i<j}a_{ij}x_ix_j
\right|.
\]

This note tests whether \(F(n)/n^{3/2}\) could fail to converge for an
arithmetic reason analogous to bent versus semi-bent Walsh spectra.
The conclusion is negative for every finite congruence, divisibility,
Seidel-integrality, and conference-existence obstruction checked here.
A sparse multiplicative-scale obstruction is not ruled out.

---

## 1. Padding theorem

### 1.1 Monotonicity and one-vertex padding

For every \(n\),

\[
\boxed{
F(n)\le F(n+1)\le F(n)+n.
}                                                            \tag{1.1}
\]

For the lower inequality, restrict any signing of \(K_{n+1}\) to
\(K_n\).  Given a Boolean vector \(x\) on the restriction, the two
extensions have energies

\[
H_A(x)\pm h(x).
\]

At least one has absolute value at least \(|H_A(x)|\).  Minimizing over
the parent signing gives \(F(n+1)\ge F(n)\).

For the upper inequality, append one arbitrarily signed vertex to an
optimal order-\(n\) signing.  The new linear field has absolute value
at most \(n\), giving \(F(n+1)\le F(n)+n\).

It follows immediately that

\[
\frac{F(n+1)}{(n+1)^{3/2}}-\frac{F(n)}{n^{3/2}}\longrightarrow0. \tag{1.2}
\]

In particular, **even and odd orders cannot have different
subsequential constants**.

### 1.2 Random rectangular padding

For arbitrary \(n,h\), choose optimal signings \(A,D\) of orders
\(n,h\), and fill the \(n\times h\) cross block \(B\) with independent
signs.  For fixed Boolean \(x,y\), \(x^\top By\) is a sum of \(nh\)
independent signs.  Hoeffding and a union bound over at most
\(2^{n+h}\) pairs give a deterministic \(B\) with

\[
\|B\|_{\infty\to1}
\le
\sqrt{2nh(n+h+2)\log2}.                                     \tag{1.3}
\]

Therefore

\[
\boxed{
F(n+h)
\le
F(n)+F(h)+
\sqrt{2nh(n+h+2)\log2}.
}                                                            \tag{1.4}
\]

A separate random-sign union bound gives the explicit estimate

\[
F(h)\le
\sqrt{2\binom h2(h+2)\log2}
=O(h^{3/2}).                                                 \tag{1.4a}
\]

Combining this with monotonicity yields the scale-local continuity
theorem

\[
\boxed{
h=o(n)
\quad\Longrightarrow\quad
F(n+h)-F(n)=o(n^{3/2}).
}                                                            \tag{1.5}
\]

Consequently,

\[
\boxed{
h=o(n)
\quad\Longrightarrow\quad
\frac{F(n+h)}{(n+h)^{3/2}}
-
\frac{F(n)}{n^{3/2}}
\longrightarrow0.
}                                                            \tag{1.6}
\]

For completeness, the normalized comparison is two-sided.  Put
\(m=n+h\) and \(a_k=F(k)/k^{3/2}\).  Since \(F(k)=O(k^{3/2})\),
monotonicity and (1.5) give

\[
\begin{aligned}
|a_m-a_n|
&\le
\frac{F(m)-F(n)}{m^{3/2}}
+
\frac{F(n)}{n^{3/2}}
\left|1-\left(\frac nm\right)^{3/2}\right|\\
&=o(1)+O(h/n)=o(1).
\end{aligned}                                                \tag{1.6a}
\]

### 1.3 Dense-subsequence transfer

Let \(\mathcal S\subseteq\mathbb N\) have relative gaps tending to
zero: for every large \(n\), some \(s\in\mathcal S\) satisfies

\[
n\le s=n+o(n)
\]

uniformly.  Equations (1.5)--(1.6) imply that the full normalized
sequence and its restriction to \(\mathcal S\) have the same liminf
and limsup.  Hence

\[
\boxed{
\lim_n\frac{F(n)}{n^{3/2}}\ \text{exists}
\quad\Longleftrightarrow\quad
\lim_{\substack{n\to\infty\\n\in\mathcal S}}
\frac{F(n)}{n^{3/2}}\ \text{exists},
}
                                                               \tag{1.7}
\]

and the limits agree when they exist.

Indeed, map any sequence \(n_k\to\infty\) to such points
\(s_k\in\mathcal S\).  Equation (1.6a) says \(a_{s_k}-a_{n_k}\to0\),
so every cluster point of the full sequence is a cluster point of the
restricted sequence.  The reverse inclusion is immediate because
\(\mathcal S\) is a subsequence.  This proves equality of the complete
cluster sets, not only equality conditional on convergence.

This rules out distinct constants on:

* even versus odd orders;
* any fixed residue classes modulo \(m\);
* any finite collection of fixed \(2\)-adic classes;
* any other syndetic or \(o(n)\)-gap arithmetic family.

Using the prime number theorem in arithmetic progressions, it also
recovers the equivalence with orders \(q+1\), where
\(q\equiv1\pmod4\) is prime: the next such order is \(n+o(n)\).

Thus a genuine arithmetic nonconvergence mechanism would have to live
on a sparse sequence with multiplicative-size gaps.  Finite parity or
congruence cannot do it.

---

## 2. Exact cut-code divisibility

Let \(N=\binom n2\).  In binary notation define the augmented cut code

\[
\mathcal C_n=
\left\{
(t+z_i+z_j)_{i<j}:
t,z_i\in\mathbb F_2
\right\}.
\]

If a signing is represented by \(a\in\mathbb F_2^N\), then

\[
\max_x|H_a(x)|
=N-2d(a,\mathcal C_n).
\]

Therefore

\[
\boxed{
F(n)=N-2\rho(\mathcal C_n),
}                                                            \tag{2.1}
\]

where \(\rho\) is covering radius.

The weights of \(\mathcal C_n\) are exactly

\[
s(n-s)
\quad\text{and}\quad
N-s(n-s),
\qquad 0\le s\le n.                                         \tag{2.2}
\]

Let \(\Delta_n\) be the greatest common divisor of all nonzero codeword
weights.  Then

\[
\boxed{
\Delta_n=
\begin{cases}
2,&n\equiv1\pmod4,\\
1,&\text{otherwise}.
\end{cases}
}                                                            \tag{2.3}
\]

Proof:

* If \(n\) is even, the cut weights \(n-1\) and \(2(n-2)\) have gcd
  \(1\).
* If \(n\) is odd, all cut weights are even and those same two weights
  have gcd \(2\).
* When \(n\equiv3\pmod4\), the all-one word has odd weight
  \(N\), reducing the augmented-code gcd to \(1\).
* When \(n\equiv1\pmod4\), \(N\) is even, so the gcd remains \(2\).

This sharply separates the present code from the first-order
Reed--Muller code.  An affine truth table of length \(L=2^m\) has
weight \(0,L/2\), or \(L\), so \(\operatorname{RM}(1,m)\) has growing
divisibility \(L/2\).  By contrast, the augmented cut code has
divisibility at most \(2\), uniformly in \(n\).  Therefore the direct
bent/semi-bent divisibility mechanism has no analogue here.

The dual description is equally local:

\[
\boxed{
\mathcal C_n^\perp
=
\{\text{Eulerian subgraphs of \(K_n\) having even edge count}\}.
}                                                            \tag{2.4}
\]

For \(n\ge4\), its minimum distance is \(4\), attained by a four-cycle.
Equivalently, codewords of \(\mathcal C_n\) are exactly edge labelings
whose parity sum on every triangle is the same bit \(t\).  These are
fixed-size local constraints, not constraints whose divisibility grows
with \(n\).

Finally, (2.1) itself gives only the unavoidable parity lattice

\[
\boxed{
F(n)\equiv N\pmod2.
}                                                            \tag{2.5}
\]

Its spacing is \(2/n^{3/2}\) after normalization.

---

## 3. Energy congruences inside a switching class

For a fixed signing \(A\), let

\[
\sigma(A)=\sum_{i<j}a_{ij}.
\]

If \(x\) corresponds to a vertex set \(S\), \(|S|=s\), then

\[
H_A(x)
=
\sigma(A)-2\sum_{ij\in\delta(S)}a_{ij}.
\]

Since a sum of \(s(n-s)\) signs has that same parity,

\[
\boxed{
H_A(x)
\equiv
\sigma(A)-2s(n-s)
\pmod4.
}                                                            \tag{3.1}
\]

Hence:

* if \(n\) is odd, every energy in a switching class has one residue
  modulo \(4\);
* if \(n\) is even, there are at most two residues, according to the
  parity of \(s\).

This is the complete universal mod-\(4\) restriction obtained from
switching.  It changes any extremal energy by at most \(O(1)\), hence
by \(o(n^{3/2})\).

Triangle products are switching invariants.  If

\[
\tau(A)=
\sum_{i<j<k}a_{ij}a_{jk}a_{ki},
\]

then

\[
\boxed{
\operatorname{tr}A^3=6\tau(A),
\qquad
\tau(A)\equiv\binom n3\pmod2.
}                                                            \tag{3.2}
\]

Thus the triangle imbalance has lattice spacing \(2\), and
\(\operatorname{tr}A^3\) has lattice spacing \(12\).  Even if exact
conference orthogonality asks for \(\tau=0\), a parity obstruction
forces only \(|\tau|\ge1\), which is negligible at every leading
spectral or \(n^{3/2}\) scale.

---

## 4. Seidel characteristic-polynomial parity

Every zero-diagonal sign matrix satisfies

\[
A\equiv J-I\equiv J+I\pmod2.
\]

The matrix determinant lemma, applied as a polynomial identity, gives

\[
\det(\lambda I-A)
\equiv
\begin{cases}
(\lambda+1)^n,&n\ \text{even},\\
\lambda(\lambda+1)^{n-1},&n\ \text{odd}
\end{cases}
\pmod2.                                                      \tag{4.1}
\]

Thus all characteristic-polynomial parity data depend only on the
ordinary parity of \(n\).  Equation (1.2) already proves that this
cannot create two normalized asymptotic constants.

Exact symmetric conference matrices obey

\[
A^2=(n-1)I.
\]

Their two eigenvalues have equal multiplicity, so exact equality
requires even \(n\), and further arithmetic conditions restrict the
known orders.  But exact conference nonexistence does **not** imply a
leading spectral gap.  If \(C\) is a conference matrix of order
\(N\ge n\) and \(A\) is an \(n\times n\) principal submatrix, then

\[
\|A\|_{\rm op}\le\sqrt{N-1}
\]

and hence

\[
M(A)\le\frac n2\sqrt{N-1}.                                  \tag{4.2}
\]

Paley conference orders \(N=q+1\) with prime
\(q\equiv1\pmod4\) satisfy \(N=n+o(n)\) for a suitable next prime.
Therefore, in **every** residue class of \(n\),

\[
\frac{M(A)}{n^{3/2}}\le\frac12+o(1).                         \tag{4.3}
\]

This rules out any determinant, eigenvalue-multiplicity, or exact
conference-existence obstruction that would force a positive
constant spectral gap on one congruence class.

The same fact can be phrased in fourth-moment stability terms.  For a
principal \(n\times n\) block of a conference matrix of order
\(N=n+r\),

\[
A^2-(n-1)I=rI-BB^\top,
\]

where \(B\) is the deleted \(n\times r\) cross block.  If \(r=o(n)\),
the normalized conference defect is \(o(1)\) at the scale needed for
any leading \(n^{3/2}\) stability correction.  Congruence can forbid
zero defect, but cannot force macroscopic defect.

---

## 5. What this proves, and what remains possible

### Ruled out

The following cannot cause nonconvergence of \(F(n)/n^{3/2}\):

1. even/odd order;
2. any fixed residue class or fixed \(2\)-adic class;
3. the weight divisibility of the augmented cut code;
4. the mod-\(4\) energy lattice of a switching class;
5. triangle-product parity;
6. characteristic-polynomial parity;
7. exact nonexistence of conference matrices at particular orders.

The first two are ruled out by the padding theorem, not merely by
scale heuristics.  The code calculation shows why the Reed--Muller
bent analogy does not transfer: the relevant divisibility is bounded
rather than growing.

### Not ruled out

A nonconvergent sequence could still oscillate on multiplicatively
separated scales.  For example, (1.5) does not compare order \(n\) to
order \(2n\) with \(o(n^{3/2})\) loss.  An obstruction confined to a
sparse family such as \(2^{2k+1}\) would propagate only to
\(o(n)\)-width neighborhoods of those orders, leaving macroscopic gaps
between them.

However, being a power of two supplies no universal structure to an
arbitrary edge signing.  Bent, semi-bent, Paley, and group-circulant
phenomena constrain selected algebraic construction families; they
do not give lower bounds for the minimum over all
\(\binom n2\) independent edge signs.

Thus the audit finds no credible arithmetic route to nonconvergence.
The only surviving possibility is a genuinely scale-dependent
geometric phenomenon, not a finite congruence or divisibility
phenomenon.  Ruling that out still requires the missing
proportional-scale amplification/deletion theorem.
