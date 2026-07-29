# Cut-triangle midpoint balancing for width-minimizing signings

## Status

Let

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
P(A)=\max_xH_A(x),\qquad
L(A)=\min_xH_A(x),
\]

and let

\[
W(A)=\frac{P(A)-L(A)}2.
\]

This note records two exact structural lemmas for a signing minimizing
\(W(A)\):

1. every edge is certified by one of the two one-step endpoint layers;
2. after switching a top state to \(\mathbf1\), a maximum signed cut
   induces a weighted domination inequality on both of its sides.

These are stronger than the earlier one-objective edge-flip
certificate, but they do not yet prove that the energy midpoint is
\(O(n)\) or \(o(n^{3/2})\). The remaining obstruction is explicit:
one near-endpoint cut can certify quadratically many edges by
cancellation between its positive and negative members.

## 1. Exact two-sided edge-flip certificate

### Theorem 1.1

Suppose \(A\) is a global minimizer of \(W(A)\) among all order-\(n\)
signings. For every edge \(e=ij\), at least one of the following holds:

\[
\exists x\in\{\pm1\}^n:
\quad
H_A(x)\ge P(A)-2,\qquad
a_{ij}x_ix_j=-1,
\tag{1.1}
\]

or

\[
\exists y\in\{\pm1\}^n:
\quad
H_A(y)\le L(A)+2,\qquad
a_{ij}y_iy_j=+1.
\tag{1.2}
\]

#### Proof

Flip only the coefficient \(a_{ij}\), obtaining \(A'\). For a spin
state \(z\),

\[
H_{A'}(z)=H_A(z)-2a_{ij}z_iz_j.
\tag{1.3}
\]

Assume neither (1.1) nor (1.2) holds. If
\(a_{ij}z_iz_j=+1\), then (1.3) lowers the energy by two. If
\(a_{ij}z_iz_j=-1\), failure of (1.1) and the common energy lattice
give

\[
H_A(z)\le P(A)-4,
\]

so (1.3) is still at most \(P(A)-2\). Therefore

\[
P(A')\le P(A)-2.
\tag{1.4}
\]

Similarly, states with \(a_{ij}z_iz_j=-1\) are raised by two, while
failure of (1.2) implies that every state with
\(a_{ij}z_iz_j=+1\) has energy at least \(L(A)+4\). Hence

\[
L(A')\ge L(A)+2.
\tag{1.5}
\]

Equations (1.4)--(1.5) imply

\[
W(A')\le W(A)-2,
\]

contradicting minimality. \(\square\)

### Coding formulation

Let \(N=\binom n2\), encode \(A\) by an edge word \(a\), and let
\(C_n\) be the cut code. Put

\[
f=d(a,C_n),\qquad
g=d(a,\mathbf1+C_n).
\]

Minimizing \(W\) is equivalent to maximizing \(f+g\). A top state with
gap at most two corresponds to a word of weight at most \(f+1\) in
the coset \(a+C_n\). A bottom state with gap at most two corresponds,
after taking its agreement set, to a word of weight at most \(g+1\)
in \(a+\mathbf1+C_n\).

Thus Theorem 1.1 is equivalent to:

\[
\boxed{
\bigcup_{\substack{u\in a+C_n\\|u|\le f+1}}\operatorname{supp}u
\ \cup\
\bigcup_{\substack{v\in a+\mathbf1+C_n\\|v|\le g+1}}
\operatorname{supp}v
=E(K_n).
}
\tag{1.6}
\]

Every coordinate is covered by a leader or near-leader from one of the
two antipodal cosets.

This is the exact two-sided local certificate absent from a generic
deep-hole argument.

## 2. Cut formulation

Switch a top state to \(\mathbf1\), so

\[
H_A(\mathbf1)=P(A).
\]

For \(S\subseteq[n]\), define

\[
C(S)=\sum_{i\in S,\ j\notin S}a_{ij}.
\]

Then

\[
H_A(x^S)=P(A)-2C(S).
\]

Consequently,

\[
\boxed{
0\le C(S)\le W(A)\qquad(S\subseteq[n]),
}
\tag{2.1}
\]

and some cut \(T\) has \(C(T)=W(A)\).

For a positive edge \(e\), Theorem 1.1 says more concretely:

\[
\boxed{
\begin{array}{l}
\text{either some cut crossing }e\text{ has }C(S)\le1,\\
\text{or some cut avoiding }e\text{ has }C(S)\ge W(A)-1.
\end{array}
}
\tag{2.2}
\]

The analogous statement after orienting at a bottom state supplies the
complementary constraints on negative edges.

## 3. Exact symmetric-difference identity

For two vertex sets \(S,T\), let

\[
D(S,T)
=\sum_{\substack{ij:\,ij\text{ crosses }S\\
ij\text{ crosses }T}}a_{ij}.
\]

The cut identity

\[
\delta(S)\mathbin\triangle\delta(T)
=\delta(S\mathbin\triangle T)
\]

gives

\[
\boxed{
C(S\mathbin\triangle T)
=C(S)+C(T)-2D(S,T).
}
\tag{3.1}
\]

If \(S,T\) are both low cuts, \(C(S),C(T)\le1\), then

\[
D(S,T)\le1
\tag{3.2}
\]

because the left side of (3.1) is nonnegative. If \(S,T\) are both
high cuts, applying (3.1) to complements and using
\(C(S\triangle T)\le W\) gives the corresponding four-cell
restrictions. These identities control signed sums of opposite blocks
in the Venn partition, but not their unsigned sizes.

## 4. Domination inside a maximum cut

Let \(T=U\subset[n]\) be a maximum signed cut:

\[
C(U)=W.
\]

For \(R\subseteq U\), write

\[
A_U(R)
=\sum_{i\in R,\ j\in U\setminus R}a_{ij},
\]

and

\[
B_U(R)
=\sum_{i\in R,\ j\in U^c}a_{ij}.
\]

The cuts \(R\) and \(U\triangle R=U\setminus R\) have values

\[
C(R)=A_U(R)+B_U(R)\ge0
\tag{4.1}
\]

and

\[
C(U\triangle R)
=W-B_U(R)+A_U(R)\le W.
\tag{4.2}
\]

Therefore:

\[
\boxed{
|A_U(R)|\le B_U(R)
\qquad(R\subseteq U).
}
\tag{4.3}
\]

In particular \(B_U(R)\ge0\) for every \(R\subseteq U\). Since
\(B_U(R)\) is additive in the indicator of \(R\), this is equivalent
to nonnegative cross row sums

\[
b_i=\sum_{j\in U^c}a_{ij}\ge0
\qquad(i\in U).
\tag{4.4}
\]

Moreover,

\[
\boxed{
\left|
\sum_{i\in R,\ j\in U\setminus R}a_{ij}
\right|
\le
\sum_{i\in R}b_i.
}
\tag{4.5}
\]

The same theorem holds on \(U^c\), with its cross column sums.

Thus each induced signing is dominated, on every cut, by a
nonnegative vertex measure whose total mass is

\[
\sum_{i\in U}b_i=W.
\tag{4.6}
\]

### Midpoint identity

Let

\[
I_U=\sum_{\{i,j\}\subset U}a_{ij},\qquad
I_{U^c}=\sum_{\{i,j\}\subset U^c}a_{ij}.
\]

Since the total edge sum is \(P(A)\) and the maximum cross sum is
\(W(A)\),

\[
\boxed{
d(A):=\frac{P(A)+L(A)}2=P(A)-W(A)
=I_U+I_{U^c}.
}
\tag{4.7}
\]

Singleton instances of (4.5) give only

\[
|I_U|\le\frac W2,\qquad
|I_{U^c}|\le\frac W2,
\]

and hence the trivial \(|d|\le W\). Improving this to \(O(n)\) requires
using the entire cut-domination family (4.5), not merely singleton
rows.

## 5. Current obstruction

The two-family certificate (1.6) is a support-cover statement. A
single low cut can cross \(\Theta(n^2)\) positive edges while keeping
signed value \(0\) or \(1\), because nearly as many negative edges
cross it. Likewise, a single high cut can avoid quadratically many
edges.

The triangle identity (3.1) controls signed sums in intersections of
certificate cuts, but it does not directly control their unsigned
coverage. Therefore a naive union bound over low/high certificate cuts
cannot bound the midpoint.

The sharpened remaining target is:

> Use (1.6), the full symmetric-difference closure (3.1), and the
> weighted cut domination (4.5) to show
> \[
> |I_U+I_{U^c}|=O(n)
> \quad\text{or at least }o(n^{3/2})
> \]
> for a width-minimizing signing.

The next two results delimit this target.  First, even the exact
two-sided edge certificate together with all triangle identities is
insufficient.  Second, even the correct \(n^{3/2}\) width scale
together with the full cut-domination family is insufficient.  A
successful argument must use *global* width minimality in a genuinely
nonlocal way.

## 6. An all-orders local trap

### Proposition 6.1

The all-positive signing of \(K_n\) is an edgewise flat local minimum
of \(W\): flipping any one edge leaves \(W\) unchanged.  Its midpoint
and width are

\[
d_+(n)=\frac{n^2-2n+\epsilon_n}{4},
\qquad
W_+(n)=\frac{n^2-\epsilon_n}{4},
\qquad
\epsilon_n=n\bmod 2.
\tag{6.1}
\]

In particular, an edgewise local minimum may have midpoint
\(\Theta(n^2)\).

#### Proof

For the all-positive signing,

\[
H_+(x)=\frac{\left(\sum_i x_i\right)^2-n}{2}.
\tag{6.2}
\]

Thus

\[
P_+=\binom n2,\qquad
L_+=\frac{\epsilon_n-n}{2},
\tag{6.3}
\]

which gives (6.1).

Now flip the edge \(12\).  The new energy is

\[
H'(x)=H_+(x)-2x_1x_2.
\tag{6.4}
\]

The all-equal states show \(P'\ge P_+-2\).  If \(x_1x_2=+1\), then
\(H'(x)\le P_+-2\).  If \(x_1x_2=-1\), the largest possible absolute
magnetization is \(n-2\), so

\[
H'(x)\le \frac{(n-2)^2-n}{2}+2\le P_+-2
\qquad(n\ge3).
\tag{6.5}
\]

Hence \(P'=P_+-2\).

There is a minimum-energy balanced or nearly balanced state having
\(x_1x_2=+1\).  It gives \(L'\le L_+-2\).  Conversely, (6.4) is at
least \(L_+-2\) when \(x_1x_2=+1\), and at least \(L_++2\) when
\(x_1x_2=-1\).  Therefore \(L'=L_+-2\), and

\[
\frac{P'-L'}2=\frac{P_+-L_+}2.
\]

Every edge is equivalent by symmetry. \(\square\)

The all-positive signing satisfies the local disjunction of Theorem
1.1: every edge is certified on the bottom side by a balanced or
nearly balanced state.  It is also a genuine complete-graph signing,
so all cut-triangle identities hold exactly.  Consequently no theorem
using only single-edge local optimality, the support cover (1.6), and
formal triangle closure can prove midpoint balancing.

### Exact small-order audit

The program `width_local_minima_audit.cpp` exhaustively enumerates the
\(2^{\binom{n-1}{2}}\) switching classes, evaluates every spin state,
and tests all \(\binom n2\) original edge flips.  It gives:

| \(n\) | global \(W_n\) | edgewise local minima | strict local minima | largest \(|d|\) among local minima |
|---:|---:|---:|---:|---:|
| 3 | 2 | 2 | 0 | 1 |
| 4 | 4 | 8 | 0 | 2 |
| 5 | 4 | 34 | 12 | 4 |
| 6 | 5 | 94 | 12 | 6 |
| 7 | 8 | 10,284 | 720 | 9 |
| 8 | 10 | 267,366 | 4,200 | 12 |

For \(n\ge5\), the largest-midpoint examples in this audit are
nonglobal local traps; the all-positive switching class is one of
them.

## 7. Correct-scale obstruction to cut geometry alone

The preceding example has width \(\Theta(n^2)\), so one might hope
that the known upper bound \(W_n=O(n^{3/2})\), combined with the full
cut identities, forces centering.  It does not.

### Proposition 7.1

There is a sequence of complete-graph signings \(A_n\) such that

\[
W(A_n)=O(n^{3/2})
\tag{7.1}
\]

but

\[
|d(A_n)|=\Theta(n^{3/2}).
\tag{7.2}
\]

After switching a top state to \(\mathbf1\), these signings satisfy
the full cut interval \(0\le C(S)\le W(A_n)\), every
symmetric-difference and triangle identity, and the maximum-cut
domination theorem (4.3)--(4.5).

#### Construction

Choose

\[
k=\lfloor n^{3/4}\rfloor+O(1),\qquad m=n-k,
\]

with \(m\) even, and write \(m=2r\).  On the first \(k\) vertices put
the all-positive signing.  On the other \(2r\) vertices use a chiral
signing

\[
D=
\begin{pmatrix}
B&C\\
C&-B
\end{pmatrix},
\tag{7.3}
\]

where \(B\) is a symmetric zero-diagonal sign matrix and \(C\) is a
symmetric full sign matrix.  Finally join the two blocks by a
\(k\times m\) sign matrix \(R\):

\[
A_n=
\begin{pmatrix}
J_k-I_k&R\\
R^\mathsf T&D
\end{pmatrix}.
\tag{7.4}
\]

Let

\[
S=
\begin{pmatrix}
0&I_r\\
-I_r&0
\end{pmatrix}.
\]

Then \(S\) is a signed permutation of the Boolean cube and

\[
S^\mathsf TDS=-D.
\tag{7.5}
\]

Therefore the energy multiset of \(D\) is centrally symmetric and

\[
d(D)=0.
\tag{7.6}
\]

Moreover \(B,C\) can be chosen so that

\[
\max_{z\in\{\pm1\}^{m}}|z^\mathsf TDz|=O(m^{3/2}).
\tag{7.7}
\]

Indeed, choose their independent upper-triangular entries randomly.
For each fixed \(z=(x,y)\), \(z^\mathsf TDz\) is a Rademacher sum
whose squared coefficient sum is \(O(r^2)\).  Hoeffding's inequality
followed by a union bound over \(2^{2r}\) states proves (7.7) with an
absolute constant.

Independently choose \(R\) randomly.  The same argument, now for
\(u^\mathsf TRv\), gives a choice satisfying

\[
\|R\|_{\infty\to1}
:=\max_{u,v}|u^\mathsf TRv|
\le 2\sqrt{km(k+m)}
=O(n^{11/8})
=o(n^{3/2}).
\tag{7.8}
\]

For \(u\in\{\pm1\}^k,v\in\{\pm1\}^m\),

\[
H_{A_n}(u,v)
=H_+(u)+H_D(v)+u^\mathsf TRv.
\tag{7.9}
\]

Adding a uniformly bounded perturbation changes each endpoint by at
most its uniform norm.  Hence

\[
W(A_n)
\le W_+(k)+W(D)+\|R\|_{\infty\to1}
=O(k^2+m^{3/2}+n^{11/8})
=O(n^{3/2}),
\tag{7.10}
\]

while, using (6.1), (7.6), and (7.8),

\[
\begin{aligned}
d(A_n)
&=d_+(k)+d(D)+O(\|R\|_{\infty\to1})\\
&=\frac{k^2}{4}+o(n^{3/2})
=\left(\frac14+o(1)\right)n^{3/2}.
\end{aligned}
\tag{7.11}
\]

This proves the proposition. \(\square\)

Proposition 7.1 is not a counterexample to midpoint balancing for
*global width minimizers*: the constructed \(A_n\) need not minimize
\(W\).  Its role is sharper.  It proves that the following information
still does not imply centering:

1. the correct \(O(n^{3/2})\) width scale;
2. every exact cut and triangle identity;
3. the complete maximum-cut domination family.

Thus the only remaining source of a midpoint theorem is comparison
with other signings, not the intrinsic cut geometry of one signing.

## 8. Revised frontier

The hoped-for bound

\[
F_n\le W_n+o(n^{3/2})
\tag{8.1}
\]

remains open.  The route has isolated the necessary extra ingredient:
a global replacement or purification lemma.

A useful theorem would have the form:

> If a signing with \(W=O(n^{3/2})\) has
> \(|d|\ge\varepsilon n^{3/2}\), then one can replace a structured
> vertex block or a positive fraction of its edges and reduce its
> half-range by \(c(\varepsilon)n^{3/2}\).

The planted construction in Proposition 7.1 shows that such a lemma
cannot merely detect a large induced signed sum.  It must also use
that *every* global replacement fails for a width minimizer.

At present the exact status is:

\[
\boxed{
\begin{array}{l}
\text{edge-local optimality + triangle identities: insufficient;}\\
\text{correct scale + full cut domination: insufficient;}\\
\text{global width minimality: the only unexploited midpoint input.}
\end{array}}
\]
