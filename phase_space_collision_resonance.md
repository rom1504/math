# Phase-space collision versus Boolean resonance

## Setup and normalization

Let \(V=\mathbb F_2^d\), \(n=|V|=2^d\), and let
\[
U_{u,s}=n^{-1/2}(-1)^{u\cdot s}
\]
be the normalized Walsh transform.  Suppose
\[
\pi:V\setminus\{0\}\longrightarrow V\setminus\{0\}
\]
is a permutation satisfying
\[
r\cdot\pi(r)=0 \qquad(r\ne0).
\]
For signs \(\epsilon_r\in\{\pm1\}\), define the symmetric zero-diagonal
sign matrix
\[
A_{uv}=\epsilon_{u+v}(-1)^{v\cdot\pi(u+v)}
\qquad(u\ne v).
\]
We use
\[
Q(A)=\max_{x\in\{\pm1\}^n}|x^\top Ax|,
\]
which is twice the one-copy objective.

## 1. Weyl representation and exact Fourier dual

For \(r\ne0\), define
\[
(P_rf)(u)=(-1)^{u\cdot\pi(r)}f(u+r).
\]
The orthogonality condition \(r\cdot\pi(r)=0\) implies
\[
P_r^\top=P_r,\qquad P_r^2=I,
\]
and
\[
A=\sum_{r\ne0}\epsilon_rP_r.
\]
Writing
\[
w_r=(r,\pi(r))\in V\oplus V
\]
and equipping \(V\oplus V\) with the symplectic form
\[
\langle(r,p),(s,q)\rangle=r\cdot q+s\cdot p,
\]
one has
\[
P_rP_s=(-1)^{\langle w_r,w_s\rangle}P_sP_r.
\]
Thus \(P_r,P_s\) commute exactly when
\[
r\cdot\pi(s)+s\cdot\pi(r)=0.
\]

A direct Walsh calculation gives, for \(s\ne t\),
\[
(UAU)_{st}
=
\epsilon_{\pi^{-1}(s+t)}
(-1)^{s\cdot\pi^{-1}(s+t)},
\]
and its diagonal is zero.  Hence
\[
D=-UAU
\]
is again a zero-diagonal sign matrix, with parameters
\[
\pi_D=\pi^{-1},\qquad
\epsilon^D_q=-\epsilon_{\pi^{-1}(q)}.
\]
So this phase-space family is exactly closed under one Walsh
anti-conjugation.

## 2. The doubled parent: spectrum and fourth-moment regularization

Set
\[
G=\begin{pmatrix}A&\sqrt n\,U\\
\sqrt n\,U&D\end{pmatrix}
=
\begin{pmatrix}A&W\\W&D\end{pmatrix}.
\]
Conjugating by \(\operatorname{diag}(I,U)\) gives
\[
\operatorname{diag}(I,U)\,G\,\operatorname{diag}(I,U)
=
\begin{pmatrix}A&\sqrt n\,I\\
\sqrt n\,I&-A\end{pmatrix}.
\]
Consequently
\[
G^2=
\begin{pmatrix}
A^2+nI&0\\
0&U(A^2+nI)U
\end{pmatrix},
\]
and if \(\lambda_1,\dots,\lambda_n\) are the eigenvalues of \(A\), the
eigenvalues of \(G\) are
\[
\pm\sqrt{\lambda_i^2+n}\qquad(1\le i\le n).
\]
In particular,
\[
\|G\|_{\rm op}=\sqrt{\|A\|_{\rm op}^2+n}.
\]

The conference defect obeys the exact identity
\[
\boxed{\;
\|G^2-(2n-1)I\|_F^2
=2\|A^2-(n-1)I\|_F^2.
\;}
\]
Thus the normalized fourth-moment defect is divided asymptotically by
four in one doubling.  This would be a useful regularizer if an exact
dualizing Hadamard could be found again at the new order.

## 3. Exact nonlinear collision formula

Let
\[
\mathcal P=\bigl\{\{r,s\}:r<s,\ 
\langle w_r,w_s\rangle=0\bigr\}
\]
be the set of commuting unordered pairs and put
\[
C(\pi)=|\mathcal P|.
\]
For a commuting pair, \(P_rP_s\) is, up to a sign, the Weyl operator
with phase-space label \(w_r+w_s\).  Weyl operators with different
labels are orthogonal in Frobenius inner product.  Therefore
\[
A^2-(n-1)I
=
2\sum_{\{r,s\}\in\mathcal P}
\epsilon_r\epsilon_sP_rP_s.
\]
For every nonzero phase-space label \(t\), let
\[
\mathcal P_t=\{\{r,s\}\in\mathcal P:w_r+w_s=t\}.
\]
There are fixed signs \(\sigma_{rs}\in\{\pm1\}\), determined only by
the Weyl convention, such that
\[
\boxed{\;
\|A^2-(n-1)I\|_F^2
=4n\sum_t
\left(
\sum_{\{r,s\}\in\mathcal P_t}
\sigma_{rs}\epsilon_r\epsilon_s
\right)^2.
\;} \tag{3.1}
\]

Three consequences are immediate.

First, averaging over independent \(\epsilon_r\) gives
\[
\boxed{\;
\mathbb E_\epsilon\|A^2-(n-1)I\|_F^2
=4nC(\pi).
\;} \tag{3.2}
\]

Second, if \(U(\pi)\) is the number of labels \(t\) for which
\(|\mathcal P_t|=1\), then for every choice of signs
\[
\boxed{\;
\|A^2-(n-1)I\|_F^2\ge4nU(\pi).
\;} \tag{3.3}
\]

Third, put
\[
T(\pi)=\sum_t\binom{|\mathcal P_t|}{2}.
\]
Since
\[
\sum_{k_t\ge2}k_t\le
2\sum_t\binom{k_t}{2},
\]
one has
\[
U(\pi)\ge C(\pi)-2T(\pi).
\]
Consequently
\[
\boxed{\;
\|A^2-(n-1)I\|_F^2
\ge4n\bigl(C(\pi)-2T(\pi)\bigr)_+.
\;} \tag{3.4}
\]
The quantity \(T(\pi)\) counts the nonlinear additive parallelograms
\[
r+s=p+q,\qquad
\pi(r)+\pi(s)=\pi(p)+\pi(q)
\]
formed by two distinct commuting pairs.

## 4. A universal polar-graph lower bound

Every \(w_r\) is singular for the hyperbolic quadratic form
\[
q(r,p)=r\cdot p.
\]
The nonzero singular vectors form the hyperbolic quadric
\(Q^+(2d-1,2)\).  Its collinearity graph joins two distinct singular
vectors when their symplectic product is zero.  Here is a direct
derivation of the spectral data used below.

There are \(n-1\) singular vectors with first coordinate zero, and for
each nonzero first coordinate there are \(n/2\) choices of the second
coordinate.  Thus the number of vertices is
\[
v=\left(\frac n2+1\right)(n-1). \tag{4.0a}
\]
The orthogonal group is transitive on these vertices, so fix
\(w=(e_1,0)\).  A perpendicular singular vector \((x,y)\) has
\(y_1=0\).  There are \(n\) choices when \(y=0\), and \(n/2\) choices
of \(x\) for each of the \(n/2-1\) nonzero possibilities for \(y\).
After removing \(0\) and \(w\), the degree is
\[
k=2\left(\frac n4+1\right)\left(\frac n2-1\right). \tag{4.0b}
\]

For two nonperpendicular singular vectors, use the representative
\((e_1,0),(0,e_1)\).  A common perpendicular vector has
\(x_1=y_1=0\), so counting singular vectors in dimension
\(2(d-1)\), and then removing zero, gives
\[
\mu=\frac{n^2}{8}+\frac n4-1.
\]
For two distinct perpendicular singular vectors, use
\((e_1,0),(e_2,0)\).  Now \(y_1=y_2=0\); the same elementary count,
followed by removal of \(0\) and the two fixed vertices, gives
\[
\lambda=\frac{n^2}{8}+\frac n2-3.
\]
The strongly regular graph identity
\[
\mathcal A^2=(k-\mu)I+(\lambda-\mu)\mathcal A+\mu J
\]
shows that on the orthogonal complement of the constant vector the
two eigenvalues are the roots of
\[
t^2-(\lambda-\mu)t-(k-\mu)=0.
\]
They are
\[
\frac n2-1
\qquad\hbox{and}\qquad
s=-\left(\frac n4+1\right). \tag{4.0c}
\]
In particular, (4.0c) proves directly that \(s\) is the least
adjacency eigenvalue.

For a vertex subset of size \(m=n-1\), the least-eigenvalue bound gives
\[
2C(\pi)
\ge
\frac{k m^2}{v}
+s\left(m-\frac{m^2}{v}\right).
\]
Substitution and simplification yield
\[
\boxed{\;
C(\pi)\ge
\frac{(n-1)(n+4)(n-4)}{8(n+2)}
=\left(\frac18+o(1)\right)n^2.
\;} \tag{4.1}
\]

Combining (3.2) and (4.1),
\[
\boxed{\;
\mathbb E_\epsilon
\frac{\operatorname{tr}A^4}{n(n-1)^2}
\ge
1+
\frac{(n+4)(n-4)}{2(n+2)(n-1)}
=\frac32-o(1).
\;} \tag{4.2}
\]
The Haar-involution/conference value of the normalized fourth moment
is \(1\).  Thus random coefficient signs in this family cannot be a
fourth-moment-unbiased conference construction.  More strongly, any
Sidon-like phase-space graph with \(T(\pi)=o(n^2)\) satisfies, for
every coefficient signing,
\[
\|A^2-(n-1)I\|_F^2
\ge\left(\frac12-o(1)\right)n^3.
\]

## 5. An explicit admissible Sidon family

The preceding Sidon obstruction is attained by a natural nonlinear
orthogonal permutation.

Identify \(V\) with \(\mathbb F_{2^d}\) and use the trace pairing
\[
x\cdot y=\operatorname{Tr}_{\mathbb F_{2^d}/\mathbb F_2}(xy).
\]
Choose \(\alpha\ne0\) with \(\operatorname{Tr}(\alpha)=0\), and set
\[
\pi(r)=\frac{\alpha}{r}\qquad(r\ne0).
\]
Then \(\pi\) is a permutation and
\[
r\cdot\pi(r)=\operatorname{Tr}(\alpha)=0.
\]

The phase-space graph
\[
S_\alpha=\{(r,\alpha/r):r\ne0\}
\]
is Sidon for sums of two distinct points.  Indeed, if
\[
r+s=p+q=t\ne0
\]
and
\[
\frac1r+\frac1s=\frac1p+\frac1q,
\]
then
\[
\frac{t}{rs}=\frac{t}{pq},
\]
so \(rs=pq\), and the unordered pairs \(\{r,s\}\) and \(\{p,q\}\)
are the same roots of the same quadratic.

Hence every commuting pair occupies a singleton collision class:
\[
U(\pi)=C(\pi),\qquad T(\pi)=0,
\]
and (3.1) becomes independent of the coefficient signs:
\[
\boxed{\;
\|A^2-(n-1)I\|_F^2=4nC(\pi)
\quad\hbox{for every }\epsilon.
\;} \tag{5.1}
\]

The commuting-pair count is explicit up to a Kloosterman sum.  Put
\[
K(\alpha)=
\sum_{t\in\mathbb F_{2^d}^{*}}
(-1)^{\operatorname{Tr}(\alpha t+\alpha/t)}.
\]
For \(r=ts\), the commutation condition is
\[
\operatorname{Tr}\!\left(\alpha(t+t^{-1})\right)=0.
\]
Therefore
\[
C(\pi)=
\frac{(n-1)(n-3+K(\alpha))}{4}.
\]
By the Weil bound \(|K(\alpha)|\le2\sqrt n\),
\[
C(\pi)=\left(\frac14+O(n^{-1/2})\right)n^2.
\]
Consequently, for every coefficient signing,
\[
\boxed{\;
\frac{\operatorname{tr}A^4}{n(n-1)^2}=2+O(n^{-1/2}).
\;} \tag{5.2}
\]
This explicit nonlinear family is Wigner-like rather than
conference/Haar-like already at the fourth moment.

## 6. The opposite wall: linear maps force Boolean eigenvectors

Suppose
\[
\pi(r)=Lr
\]
with \(L\) invertible.  The condition \(r\cdot Lr=0\) for every \(r\)
means that \(L\) is alternating (symmetric with zero diagonal over
\(\mathbb F_2\)); in particular, this is possible only for even \(d\).

Choose a quadratic form \(q\) with polar form
\[
q(u+r)+q(u)+q(r)=u\cdot Lr.
\]
Switching vertices by \(s_u=(-1)^{q(u)}\) turns \(A\) into the Cayley
matrix
\[
A'_{uv}=a_{u+v},\qquad
a_0=0,\quad
a_r=\epsilon_r(-1)^{q(r)}.
\]
The Walsh characters are Boolean eigenvectors of \(A'\), with
eigenvalues
\[
\widehat a(k)=\sum_{r\in V}a_r(-1)^{k\cdot r}.
\]
It follows exactly that
\[
\boxed{\;
Q(A)=n\max_k|\widehat a(k)|.
\;} \tag{6.1}
\]
Parseval gives
\[
\max_k|\widehat a(k)|\ge\sqrt{n-1},
\]
and hence
\[
\boxed{\;
Q(A)\ge n\sqrt{n-1}.
\;} \tag{6.2}
\]
Conversely, when \(d\) is even, take a bent Boolean function \(b\)
and set \(a_r=b_r\) for \(r\ne0\), \(a_0=0\).  Since every Walsh
coefficient of \(b\) has magnitude \(\sqrt n\),
\[
\max_k|\widehat a(k)|\le\sqrt n+1.
\]
Thus the best doubled normalized Boolean maximum in the linear
subfamily tends exactly to
\[
\frac{Q(A)}{n^{3/2}}\longrightarrow1,
\]
or \(1/2\) in the original one-copy normalization.

So the most additive end of the family permits spectral flattening,
but simultaneously installs a complete Boolean eigenbasis and forces
the spectral-ceiling constant.

## 7. Why BSG/Freiman does not close the intermediate regime

Suppose the conference defect is \(o(n^3)\).  From (3.1), the number
of odd collision classes is \(o(n^2)\).  Together with (4.1), this
forces
\[
T(\pi)\ge \frac12C(\pi)-o(n^2)
\ge\left(\frac1{16}-o(1)\right)n^2.
\]
This is only a speed-\(n^2\) additive-energy statement.

For a set \(S\) of size \(m\asymp n\) in the ambient group
\(\mathbb F_2^{2d}\), whose size is \(n^2\), random-scale additive
energy is already \(\Theta(n^2)\).  In the usual BSG notation
\[
E(S)\ge \frac{|S|^3}{K},
\]
an \(E(S)=\Theta(n^2)\) hypothesis corresponds to
\[
K=\Theta(n).
\]
Balog--Szemerédi--Gowers then extracts at best a subset of size
\[
\frac{n}{K^{O(1)}},
\]
which is not macroscopic, and Freiman/PFR information on that subset
does not produce a leading-order Boolean witness.  The loss is a full
power of \(n\), not a removable logarithm.

Therefore collision count alone cannot prove the desired inverse
statement.  A successful inverse theorem must use the signed
cancellation equations in (3.1), simultaneously for almost every
phase-space sum, not merely the unsigned additive energy.  The exact
remaining structural question is:

> If signs \(\epsilon_r\) make the signed collision sums in (3.1)
> have total square \(o(n^2)\), must the signed graph
> \(\{(r,\pi(r),\epsilon_r)\}\) contain a macroscopic quadratic/affine
> component that yields a Boolean resonance witness?

Ordinary BSG and Freiman theory do not see the signs and are
quantitatively too weak at this threshold.

## 8. Boolean facts for the doubled parent

For all Boolean \(x,y\),
\[
(x,y)^\top G(x,y)
=x^\top Ax+y^\top Dy+2x^\top Wy.
\]
Replacing \(x\) by \(-x\) leaves the internal terms fixed and reverses
the cross term, so
\[
\boxed{\;
Q(G)=
\max_{x,y}
\left(
|x^\top Ax+y^\top Dy|+2|x^\top Wy|
\right).
\;} \tag{8.1}
\]
In particular,
\[
Q(G)\ge2\|W\|_{\infty\to1}.
\]
If \(d\) is even, a bent sign vector \(y\) has \(Uy\in\{\pm1\}^n\).
Taking \(x=Uy\) cancels the two child energies and gives
\[
\boxed{\;
Q(G)\ge2n^{3/2}.
\;} \tag{8.2}
\]

The exact Walsh-dual property need not recurse.  A necessary condition
for an \(N\times N\) signing \(B\) to have a Hadamard \(H\) with
\(H^\top BH/N\) zero-diagonal is that the zero-energy shell
\[
\{x\in\{\pm1\}^N:x^\top Bx=0\}
\]
contain \(N\) mutually orthogonal sign vectors (the columns of \(H\)).
For the order-eight parent obtained from the unique \(d=2\)
orthogonal permutation, exact enumeration gives:

- \(Q(G)=20\);
- there are \(12\) zero-energy sign lines modulo \(x\sim-x\);
- their orthogonality graph has clique number \(6\).

Thus this first parent has no dualizing Hadamard of any equivalence
class.  Exact one-step anti-conjugacy is real, but exact recursive
closure fails already at order eight.

For reproducibility, fix the first coordinate of each sign vector to
\(+1\), encode a minus sign in coordinate \(i\) by bit \(i\), and use
\(\pi(1)=2,\pi(2)=1,\pi(3)=3\) with all three coefficient signs
negative.  The twelve zero-energy masks are
\[
34,39,40,45,68,71,72,75,114,116,123,125.
\]
Direct dot products give an orthogonality graph with a clique of size
six, while checking its \(\binom{12}{7}=792\) seven-vertex subsets
shows that none is a clique.  This is a complete finite certificate
of the claimed clique number.
