# Structured flat cross blocks and exact energy-layer bounds

## 0. Normalization

For a symmetric zero-diagonal sign matrix \(A\) of order \(n\), write

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j=\frac12x^\top Ax,
\qquad
M(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|.
\]

Thus \(Q(A)=2M(A)\) in the doubled normalization used elsewhere in the
campaign.  If

\[
G=\begin{pmatrix}A&B\\B^\top&D\end{pmatrix},
\]

then

\[
H_G(x,y)=H_A(x)+H_D(y)+x^\top By
\]

and, because replacing \(x\) by \(-x\) changes only the cross term,

\[
M(G)=\max_{x,y}\bigl(
|H_A(x)+H_D(y)|+|x^\top By|
\bigr).                                                     \tag{0.1}
\]

The aim here is to understand the restricted cross norm

\[
S_{p,q}(B)=
\max_{\substack{H_A(x)=p\\H_D(y)=q}}|x^\top By|
\]

when \(B\) is a Hadamard block.

---

## 1. The exact anti-conjugate layer ellipse

Let \(B=\sqrt n\,U\), where \(U\) is orthogonal (so a sign \(B\) is a
Hadamard matrix).  Suppose first that

\[
D=-U^\top A U.                                               \tag{1.1}
\]

For Boolean \(x,y\), put \(z=Uy\).  Then

\[
H_A(x)+H_D(y)
=\frac12\bigl(x^\top Ax-z^\top Az\bigr),\qquad
x^\top By=\sqrt n\,x^\top z.
\]

After replacing \(x\) by \(-x\), assume \(x^\top z\ge0\), and set

\[
t=\frac{x^\top z}{n}\in[0,1].
\]

Since \(\|x\|_2=\|z\|_2=\sqrt n\),

\[
\|x-z\|_2\|x+z\|_2=2n\sqrt{1-t^2}.
\]

Symmetry of \(A\) gives the exact polarization

\[
x^\top Ax-z^\top Az=(x-z)^\top A(x+z).
\]

Consequently,

\[
\boxed{
|H_A(x)+H_D(y)|
\le n\|A\|_{\rm op}\sqrt{1-t^2},
\qquad
|x^\top By|=n^{3/2}t.
}                                                            \tag{1.2}
\]

In particular, on exact energy layers,

\[
\boxed{
S_{p,q}(B)
\le n^{3/2}
\sqrt{
1-\left(\frac{|p+q|}{n\|A\|_{\rm op}}\right)^2
}
}                                                            \tag{1.3}
\]

whenever the displayed ratio is at most one.  Thus the joint
internal-energy/cross-energy profile lies inside a literal ellipse.
Maximizing the support function of this ellipse gives

\[
\boxed{
M(G)\le n\sqrt{\|A\|_{\rm op}^2+n}.
}                                                            \tag{1.4}
\]

This is also immediate after the orthogonal change of variables
\(\operatorname{diag}(I,U)\):

\[
\operatorname{diag}(I,U)\,G\,\operatorname{diag}(I,U^\top)
=
\begin{pmatrix}A&\sqrt n I\\ \sqrt n I&-A\end{pmatrix},
\]

whose square is

\[
\begin{pmatrix}A^2+nI&0\\0&A^2+nI\end{pmatrix}.
\]

The layer proof is stronger than this spectral observation because it
identifies exactly where cross energy must disappear as internal energy
grows.

### Approximate anti-conjugacy

More generally, write

\[
D=-U^\top A U+E.
\]

Then

\[
\boxed{
S_{p,q}(B)\le n^{3/2}
\sqrt{
1-
\left(
\frac{(|p+q|-\frac n2\|E\|_{\rm op})_+}
{n\|A\|_{\rm op}}
\right)^2
}
}                                                            \tag{1.5}
\]

and

\[
\boxed{
M(G)
\le n\sqrt{\|A\|_{\rm op}^2+n}
+\frac n2\|E\|_{\rm op}.
}                                                            \tag{1.6}
\]

In total-order normalization \(N=2n\), with
\(r=\|A\|_{\rm op}/\sqrt n\),

\[
\frac{M(G)}{N^{3/2}}
\le
\frac{\sqrt{1+r^2}}{2\sqrt2}
+\frac{\|E\|_{\rm op}}{4\sqrt{2n}}.                          \tag{1.7}
\]

Because \(\|A\|_{\rm op}\ge\|A\|_F/\sqrt n=\sqrt{n-1}\), this
particular deterministic certificate has a hard asymptotic value
\(1/2\), attained at the conference spectral scale.  Anti-conjugacy is
therefore the right way to prevent the *sum* of the internal and cross
spectral costs, but the operator-norm ellipse by itself cannot certify a
constant below \(1/2\).

---

## 2. A sharp Hadamard exponential-moment inequality

Let \(B\) be an \(n\times n\) Hadamard matrix and \(U=B/\sqrt n\).
For every fixed Boolean \(y\) and every real \(t\),

\[
\boxed{
\mathbb E_x
\exp\!\left(\frac{t}{\sqrt n}x^\top By\right)
\le(\cosh t)^n.
}                                                            \tag{2.1}
\]

Indeed, put \(z=Uy\), so \(\sum_i z_i^2=n\).  Conditioning on \(y\),

\[
\mathbb E_x e^{t x^\top z}=\prod_i\cosh(tz_i).
\]

The function

\[
f(s)=\log\cosh(t\sqrt s),\qquad s\ge0,
\]

is concave: its derivative is

\[
f'(s)=\frac{t^2}{2}\frac{\tanh u}{u},
\qquad u=t\sqrt s,
\]

and \(\tanh u/u\) decreases on \(u\ge0\).  Jensen therefore gives

\[
\sum_i\log\cosh(tz_i)
\le n\log\cosh t.
\]

Equality holds precisely when \(|(Uy)_i|=1\) for every \(i\), i.e. at
a Boolean-to-Boolean (bent) Hadamard pair.

There is also a rectangular version.  If \(B\) is \(m\times n\) and
\(BB^\top=nI_m\), then

\[
\mathbb E_x
\exp\!\left(\frac{t}{\sqrt n}x^\top By\right)
\le
\left[
\cosh\!\left(t\sqrt{\frac nm}\right)
\right]^m.                                                   \tag{2.2}
\]

This follows from \(\|By/\sqrt n\|_2^2\le n\) and the same concavity
argument.

---

## 3. Exact switched-block free-energy recursion, and its limitation

For a Hamiltonian \(A\), define

\[
Z_A(\lambda)=\sum_xe^{\lambda H_A(x)},\qquad
\mathcal Z_A(\lambda)=Z_A(\lambda)+Z_A(-\lambda).
\]

Randomize a fixed Hadamard cross block by independent row and column
switches:

\[
B_{s,t}=\operatorname{diag}(s)B\operatorname{diag}(t).
\]

For every fixed \(x,y\), the switched vectors
\(\operatorname{diag}(s)x,\operatorname{diag}(t)y\) are independent
uniform cube points.  If

\[
\varphi_B(\lambda)=
\mathbb E_{u,v}e^{\lambda u^\top Bv},
\]

then

\[
\mathbb E_{s,t}Z_G(\lambda)
=\varphi_B(\lambda)Z_A(\lambda)Z_D(\lambda)
\]

and the same identity holds at \(-\lambda\).  Hence some switching
satisfies

\[
\boxed{
\mathcal Z_G(\lambda)
\le
\varphi_B(\lambda)\,
\mathcal Z_A(\lambda)\mathcal Z_D(\lambda).
}                                                            \tag{3.1}
\]

Define the two-sided scaled free energy

\[
f_A(\beta)=\frac1n
\log\mathcal Z_A\!\left(\frac{\beta}{\sqrt n}\right).
\]

For equal child orders, (2.1) and (3.1) give a switching with

\[
\boxed{
f_G(\beta)
\le
\frac12f_A(\beta/\sqrt2)
+\frac12f_D(\beta/\sqrt2)
+\frac12\log\cosh(\beta/\sqrt2).
}                                                            \tag{3.2}
\]

If this construction is iterated on a dyadic hierarchy and only
(2.1) is used to control every cross block, then

\[
\limsup_k f_{A_k}(\beta)
\le
\log2+
\frac12\sum_{j\ge1}
\log\cosh\!\left(\frac{\beta}{2^{j/2}}\right).
                                                               \tag{3.3}
\]

Here the seed order is taken to infinity before the number of
amplification levels; a fixed seed contributes the harmless additional
term \((\log2)/n_{\rm seed}\) coming from the two-sided partition
function.

Since

\[
\frac{M(A_k)}{|A_k|^{3/2}}
\le\frac{f_{A_k}(\beta)}{\beta},
\]

the best constant certified by this recursion is

\[
\inf_{\beta>0}
\frac1\beta
\left[
\log2+\frac12\sum_{j\ge1}
\log\cosh\!\left(\frac{\beta}{2^{j/2}}\right)
\right]
=0.7666897549\ldots                                         \tag{3.4}
\]

at \(\beta=2.194\ldots\).  This is far above \(1/2\).  Thus
independent row/column randomization plus a one-parameter annealed
free-energy profile is not a viable flat amplification proof.  A
successful flat construction must correlate the children with the
cross basis (as in Section 1), not merely randomize their relative
switching.

---

## 4. An infinite exact Walsh anti-conjugate sign family

The exact condition (1.1) is not empty.  It has a useful phase-space
description.

Let \(V=\mathbb F_2^d\), \(n=2^d\), and let

\[
W_{u,v}=(-1)^{u\cdot v}
\]

be the Walsh Hadamard matrix.  Let

\[
\pi:V\setminus\{0\}\longrightarrow V\setminus\{0\}
\]

be a permutation satisfying

\[
r\cdot\pi(r)=0\qquad(r\ne0),                                 \tag{4.1}
\]

and let \(\varepsilon_r\in\{\pm1\}\).  Define

\[
A_{u,u}=0,\qquad
A_{u,v}=
\varepsilon_{u+v}
(-1)^{v\cdot\pi(u+v)}
\quad(u\ne v).                                               \tag{4.2}
\]

Condition (4.1) makes \(A\) symmetric.  A direct Walsh transform gives

\[
\frac1n(W^\top A W)_{s,t}
=
\begin{cases}
0,&s=t,\\
\varepsilon_r(-1)^{s\cdot r},
&s\ne t,\quad r=\pi^{-1}(s+t).
\end{cases}                                                  \tag{4.3}
\]

Therefore \(W^\top A W/n\) is again a symmetric zero-diagonal sign
matrix.  Taking

\[
D=-\frac1nW^\top A W
\]

produces an exact sign parent satisfying (1.1).

There are infinite such examples.  In every even dimension \(d\), take
\(\pi(r)=Jr\), where \(J\) is an invertible alternating matrix over
\(\mathbb F_2\); then \(r\cdot Jr=0\).

Equivalently,

\[
A=\sum_{r\ne0}\varepsilon_r\,M_{\pi(r)}T_r,
\]

a signed sum of symmetric Weyl operators.  This makes clear both the
large algebraic family and its rigidity: it has only \(n-1\) freely
chosen signs once \(\pi\) is fixed, rather than \(\binom n2\).

### Linear \(\pi\) is completely resonant

The simplest infinite subfamily cannot beat the \(1/2\) barrier.  Let
\(\pi(r)=Jr\) for an invertible alternating matrix \(J\).  The Weyl
operators

\[
W_r=M_{Jr}T_r
\]

all commute, because

\[
W_rW_s=(-1)^{r\cdot Js}W_{r+s}
=W_sW_r.
\]

Choose a quadratic form \(q\) with polar form \(J\).  For every
\(s\in V\), the Boolean vector

\[
\xi_s(v)=(-1)^{q(v)+s\cdot v}
\]

is a simultaneous eigenvector, and

\[
A\xi_s
=
\left[
\sum_{r\ne0}
\varepsilon_r(-1)^{q(r)+s\cdot r}
\right]\xi_s.                                                \tag{4.4}
\]

Thus \(A\) has a complete Boolean eigenbasis.  If

\[
\eta_0=0,\qquad
\eta_r=\varepsilon_r(-1)^{q(r)}\quad(r\ne0),
\]

then, exactly,

\[
\boxed{
M(A)=\frac n2\max_s|\widehat\eta(s)|.
}                                                            \tag{4.5}
\]

Parseval gives

\[
\max_s|\widehat\eta(s)|\ge\sqrt{n-1},
\]

and hence

\[
\boxed{
\frac{M(A)}{n^{3/2}}
\ge\frac12\sqrt{1-\frac1n}.
}                                                            \tag{4.6}
\]

So every *linear* exact anti-conjugate family is forced back to the
spectral \(1/2\) scale by explicit Boolean stabilizer eigenvectors.
Any genuinely Haar-like anti-conjugate candidate must use a nonlinear
orthogonal permutation \(\pi\), equivalently a noncommuting Weyl graph.

For later use, the multiplication rule also gives a concrete measure
of this nonlinearity.  For two distinct displacements \(r,s\), the two
ordered terms in \(A^2\) cancel exactly when

\[
r\cdot\pi(s)+s\cdot\pi(r)=1,                                 \tag{4.7}
\]

and reinforce when this quantity is zero.  Thus near-conference
flatness in the nonlinear family is a signed collision/cancellation
problem among the phase-space sums

\[
(r+s,\pi(r)+\pi(s)).
\]

Absent structured collisions, this expansion has the Wigner-scale
fourth-moment count; linear graphs have maximal collisions but also the
Boolean eigenbasis (4.4).  A useful future theorem would quantify this
apparent tradeoff.

### Exact small-order audit

All statements below use the one-copy normalization \(M\).

* At \(d=2\), there is one admissible \(\pi\).  Exhausting its eight
  sign choices gives \(M(A)=6\), while every resulting
  anti-conjugate parent of order \(8\) has

  \[
  M(G)=10,\qquad \frac{M(G)}{8^{3/2}}=0.4419417\ldots.
  \]

  This happens to equal the known optimum \(F(8)=10\).

* At \(d=3\), there are exactly 24 admissible permutations and 3072
  pairs \((\pi,\varepsilon)\).  Exhaustion shows

  \[
  \min M(A)=14,
  \]

  and **every** resulting phase-space parent of order \(16\) has

  \[
  M(G)=32,\qquad \frac{M(G)}{16^{3/2}}=\frac12.
  \]

* The order-\(8\) parent from the \(d=2\) optimum cannot itself be
  anti-conjugated by any Hadamard basis.  Indeed it has only 12
  zero-energy Boolean rays, and exhaustive clique search finds no
  set of eight mutually orthogonal zero-energy rays whose pairwise
  \(A\)-bilinear values all have magnitude \(8\).  Those conditions
  are necessary for \(H^\top A H/8\) to be zero-diagonal sign.

Thus exact anti-conjugacy can produce an excellent isolated lift, but
the property is not automatically recursive.  The first nontrivial
next order returns to the \(1/2\) ceiling.

---

## 5. Relation to the bent/Walsh obstruction

The familiar Walsh-diagonal family

\[
A=\frac1{\sqrt n}W\operatorname{diag}(\lambda)W-\delta I
\]

with a bent sign mask \(\lambda\) is different from the phase-space
family in Section 4.  It always contains Boolean Walsh columns as
eigenvectors.  Consequently its one-copy Boolean maximum is already
at the spectral scale \(1/2+o(1)\), before any block lift is made.

Direct computations supplied by the parent campaign also show that
the obvious coherent parents

\[
\begin{pmatrix}A&W\\W&A\end{pmatrix},
\qquad
\begin{pmatrix}A&W\\W&-A\end{pmatrix}
\]

do not generate a low-norm recursion: at total orders \(32\) and
\(128\), greedy witnesses are already at or above roughly
\(0.51\)--\(0.62\) depending on the sign choice.  This agrees with the
structural point: exact cap pairs have low child energy, but
intermediate-overlap pairs create the obstruction.

---

## 6. Current verdict and the surviving lemma

The flat-block route now has a clean division.

1. **Independent/random relative orientation is ruled out as a proof
   mechanism.**  Bulk exact-energy layers retain essentially the full
   rectangular norm (previous note), and the sharp annealed recursion
   (3.2) is quantitatively much too expensive.

2. **Exact anti-conjugacy gives the correct layer geometry.**  The
   ellipse (1.3) is a sharp deterministic statement and completely
   cancels the additive spectral cost of the cross block.

3. **The remaining issue is genuinely arithmetic/discrete.**  To beat
   \(1/2\), one must improve the ellipse on Boolean points.  Equivalently,
   for an exact or approximate sign anti-conjugate pair, one needs a
   forbidden-band theorem saying that Boolean pairs cannot populate the
   support point

   \[
   t=\frac1{\sqrt{1+r^2}}
   \]

   of the spectral ellipse with simultaneously extremal polarization.
   The phase-space family (4.2) is a concrete test bed for precisely
   this question.

No scale-preserving recursion below \(1/2\) has been proved.  The
strongest reusable output is the exact layer ellipse (including its
robust error form) and the explicit infinite anti-conjugate sign
family on which any proposed Boolean improvement can be tested.
