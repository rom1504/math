# Global block replacement and midpoint balancing

## Status

Let

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
W(A)=\frac{\max H_A-\min H_A}{2},
\qquad
d(A)=\frac{\max H_A+\min H_A}{2},
\]

and let

\[
W_n=\min_A W(A).
\]

This note records the global optimality information that is absent
from single-edge flip arguments.

The main conclusions are:

1. an exact replacement identity reduces global optimality to the
   low-slack top/bottom-pair profiles;
2. every edge block has an approximate minimax dual certificate;
3. every induced subgraph of a global width minimizer satisfies an
   exact exchange inequality against \(W_k\);
4. this exchange inequality rules out the planted biased-block
   obstruction from `cut_triangle_midpoint_balancing.md`;
5. a macroscopic midpoint can therefore not be localized behind a
   small bilinear boundary, but a diffuse midpoint remains possible.

The last item is the surviving obstruction to

\[
|d(A)|=o(n^{3/2})
\]

for some width-minimizing signing.

## 1. Pair-profile representation

For an edge \(e=ij\), put

\[
\chi_e(x)=x_ix_j
\]

and, for a pair \(p=(x,y)\), define

\[
s_e(p)=\frac{\chi_e(x)-\chi_e(y)}2\in\{-1,0,1\}.
\tag{1.1}
\]

Then

\[
\boxed{
W(A)=\max_{p=(x,y)}\sum_e a_es_e(p).
}
\tag{1.2}
\]

For later use define the slack

\[
\delta_A(p)=W(A)-\sum_ea_es_e(p)\ge0.
\tag{1.3}
\]

If \(x\) and \(y\) are respectively a top and a bottom state, then
\(\delta_A(x,y)=0\).

## 2. Exact global replacement identity

Let \(T\subseteq E(K_n)\), and replace \(a_e\) by
\(\beta_e\in\{\pm1\}\) for \(e\in T\), leaving all other coefficients
fixed.  Denote the resulting signing by \(A^{T\to\beta}\).

### Theorem 2.1

\[
\boxed{
W(A^{T\to\beta})-W(A)
=
\max_p\left\{
\sum_{e\in T}(\beta_e-a_e)s_e(p)-\delta_A(p)
\right\}.
}
\tag{2.1}
\]

#### Proof

For every pair \(p\), its score after replacement is

\[
\sum_ea_es_e(p)
+
\sum_{e\in T}(\beta_e-a_e)s_e(p).
\]

Subtract \(W(A)\), use (1.3), and maximize over \(p\). \(\square\)

Thus \(A\) is globally width-minimizing only if, for every block
replacement \(\beta\), some pair profile satisfies

\[
\sum_{e\in T}(\beta_e-a_e)s_e(p)\ge\delta_A(p).
\tag{2.2}
\]

This is the nonlocal analogue of the one-edge endpoint certificate.
It also gives an exact entropy criterion.  If \(\boldsymbol\beta\) is
any random block replacement and

\[
\sum_p
\Pr\left[
\sum_{e\in T}(\boldsymbol\beta_e-a_e)s_e(p)
\ge\delta_A(p)
\right]<1,
\tag{2.3}
\]

then some deterministic realization strictly lowers \(W\).  Grouping
the sum by slack shells makes clear why the entropy of low-slack
profiles, rather than the number of exact endpoint pairs alone, is the
relevant statistic.

## 3. Approximate dual certificate for every block

Let \(v(T)\) be the number of vertices incident with \(T\), and put
\(m=|T|\).  For \(z\in[-1,1]^T\), define the continuous block
relaxation

\[
\Phi_T(z)=
\max_p\left\{
\sum_{e\notin T}a_es_e(p)+\sum_{e\in T}z_es_e(p)
\right\}.
\tag{3.1}
\]

Let

\[
V_T=\min_{z\in[-1,1]^T}\Phi_T(z).
\tag{3.2}
\]

### Lemma 3.1: rounding the continuous block

There is an absolute bound

\[
\min_{\beta\in\{\pm1\}^T}\Phi_T(\beta)
\le V_T+\eta_T,
\tag{3.3}
\]

where

\[
\boxed{
\eta_T=
\sqrt{2m\left(2v(T)\log2+1\right)}.
}
\tag{3.4}
\]

#### Proof

Round an optimizer \(z\) independently to signs
\(\boldsymbol\beta_e\) with
\(\mathbb E\boldsymbol\beta_e=z_e\).  For a fixed restriction of
\((x,y)\) to the vertices incident with \(T\),

\[
\sum_{e\in T}(\boldsymbol\beta_e-z_e)s_e(x,y)
\]

is a centered sum of independent variables with total squared range
at most \(4m\).  Hoeffding's inequality gives upper tail
\(\exp(-t^2/(2m))\).  There are at most \(4^{v(T)}\) restricted pair
profiles.  The union bound with \(t=\eta_T\) has probability less than
one, so some rounding increases no profile by more than \(\eta_T\).
\(\square\)

If \(A\) globally minimizes \(W\), then changing only \(T\) cannot
improve it.  Since \(a_T\) itself is feasible,

\[
\min_{\beta\in\{\pm1\}^T}\Phi_T(\beta)=W(A).
\]

Consequently,

\[
\boxed{
W(A)-\eta_T\le V_T\le W(A).
}
\tag{3.5}
\]

### Theorem 3.2: blockwise dual measure

For every edge block \(T\) of a global width minimizer \(A\), there is
a probability measure \(\mu_T\) on spin pairs \(p=(x,y)\) such that,
with

\[
m_e=\mathbb E_{\mu_T}s_e(p),
\]

one has

\[
\boxed{
\mathbb E_{\mu_T}\sum_ea_es_e(p)\ge W(A)-\eta_T
}
\tag{3.6}
\]

and

\[
\boxed{
2\sum_{e\in T}(a_em_e)_+\le\eta_T.
}
\tag{3.7}
\]

#### Proof

Write

\[
c(p)=\sum_{e\notin T}a_es_e(p).
\]

Finite-dimensional minimax gives

\[
\begin{aligned}
V_T
&=\min_{z\in[-1,1]^T}
  \max_{\mu}\mathbb E_\mu
  \left[c(p)+\sum_{e\in T}z_es_e(p)\right]\\
&=\max_{\mu}
\left\{
\mathbb E_\mu c(p)
-\sum_{e\in T}\left|\mathbb E_\mu s_e(p)\right|
\right\}.
\end{aligned}
\tag{3.8}
\]

Take a maximizing measure and write \(m_e=\mathbb E_\mu s_e\).
Then

\[
\begin{aligned}
\mathbb E_\mu\sum_ea_es_e-V_T
&=\sum_{e\in T}\left(|m_e|+a_em_e\right)\\
&=2\sum_{e\in T}(a_em_e)_+.
\end{aligned}
\tag{3.9}
\]

The left side is nonnegative, while every pair score is at most
\(W(A)\).  Combining (3.5) and (3.9) proves both claims. \(\square\)

For an induced \(k\)-vertex block,

\[
\eta_T=O(k^{3/2}).
\tag{3.10}
\]

For a \(k\)-by-\((n-k)\) cross block,

\[
\eta_T=O(n\sqrt{k})
\qquad(k\le n/2).
\tag{3.11}
\]

Theorem 3.2 is a genuine global certificate: \(\mu_T\) is concentrated
on average within \(O(\eta_T)\) of the top/bottom-pair layer, while its
mean width-gradient on \(T\) is almost never aligned with the existing
coefficients.

## 4. Exact induced-block exchange inequality

Partition the vertices into \(S\) and \(S^c\), with \(|S|=k\), and
write

\[
A=
\begin{pmatrix}
B&C\\
C^\mathsf T&D
\end{pmatrix}.
\tag{4.1}
\]

Put

\[
R(C)=\|C\|_{\infty\to1}
=\max_{u\in\{\pm1\}^k,\ v\in\{\pm1\}^{n-k}}
|u^\mathsf TCv|.
\tag{4.2}
\]

### Theorem 4.1

If \(A\) globally minimizes \(W\), then

\[
\boxed{
W(B)\le W_k+2R(C)
}
\tag{4.3}
\]

and, symmetrically,

\[
\boxed{
W(D)\le W_{n-k}+2R(C).
}
\tag{4.4}
\]

#### Proof

With the cross block deleted, the two vertex blocks are independent,
so their energy ranges add:

\[
W(B\oplus D)=W(B)+W(D).
\tag{4.5}
\]

The cross energy \(u^\mathsf TCv\) has uniform absolute value at most
\(R(C)\).  Therefore

\[
W(A)\ge W(B)+W(D)-R(C).
\tag{4.6}
\]

Replace \(B\) by a signing \(B_*\) attaining \(W_k\), while leaving
\(C,D\) fixed.  For the new full signing \(A_*\),

\[
W(A_*)\le W_k+W(D)+R(C).
\tag{4.7}
\]

Global minimality gives \(W(A)\le W(A_*)\).  Combining
(4.6)--(4.7) proves (4.3); the other side is identical. \(\square\)

This theorem is deterministic and has no rounding loss.

## 5. The planted biased block is not a minimizer

The obstruction in `cut_triangle_midpoint_balancing.md` planted an
all-positive block of order

\[
k\asymp n^{3/4}
\]

inside a centered \(O(n^{3/2})\)-width bulk, with a random cross block
satisfying

\[
R(C)=O(n^{11/8})=o(k^2).
\tag{5.1}
\]

For the all-positive block,

\[
W(B)=\frac{k^2-(k\bmod2)}4.
\tag{5.2}
\]

On the other hand, standard random-sign or conference constructions
give

\[
W_k=O(k^{3/2})=o(k^2).
\tag{5.3}
\]

Equations (5.1)--(5.3) violate (4.3) for all sufficiently large \(n\).
Thus:

\[
\boxed{
\text{the planted correct-scale midpoint obstruction cannot globally
minimize }W.
}
\tag{5.4}
\]

More generally, a block \(B\) with

\[
W(B)>W_k+2R(C)
\tag{5.5}
\]

is a deterministic certificate that the full signing is not a global
width minimizer.

## 6. What a large midpoint now forces

Switch a top state of \(A\) to \(\mathbf1\), and let
\(U|U^c\) be a maximum signed cut.  As in
`cut_triangle_midpoint_balancing.md`,

\[
d(A)=I_U+I_{U^c},
\tag{6.1}
\]

where \(I_U,I_{U^c}\) are the total signed internal edge sums.

For any signing \(B\) and any state \(z\),

\[
W(B)\ge\frac{|H_B(z)|}{2},
\tag{6.2}
\]

because the mean of \(H_B\) over the Boolean cube is zero and hence
lies between its two endpoints.

Suppose

\[
|d(A)|\ge\varepsilon n^{3/2}.
\tag{6.3}
\]

Then one side \(S\in\{U,U^c\}\) satisfies

\[
|I_S|\ge\frac{\varepsilon}{2}n^{3/2},
\]

and therefore

\[
W(A[S])\ge\frac{\varepsilon}{4}n^{3/2}.
\tag{6.4}
\]

Applying Theorem 4.1 gives the exact dichotomy

\[
\boxed{
\frac{\varepsilon}{4}n^{3/2}
\le W_{|S|}+2\|A_{S,S^c}\|_{\infty\to1}.
}
\tag{6.5}
\]

In particular, if \(|S|=o(n)\), then \(W_{|S|}=o(n^{3/2})\), so

\[
\|A_{S,S^c}\|_{\infty\to1}
\ge\left(\frac{\varepsilon}{8}-o(1)\right)n^{3/2}.
\tag{6.6}
\]

Thus a macroscopic midpoint of a global minimizer cannot sit in a
sublinear vertex block with a weak boundary.  It must be either

1. spread over a linear fraction of all vertices, or
2. protected by a boundary block having macroscopic bilinear norm.

For every partition,

\[
\|A_{S,S^c}\|_{\infty\to1}\le W(A),
\tag{6.7}
\]

because pair states differing by a flip of \(S\) realize every
bilinear orientation of that cross block.  Hence (6.6) is compatible
with the known \(W(A)=\Theta(n^{3/2})\) scale; it is structural but not
yet contradictory.

## 7. Why diffuse midpoint bias survives block replacement

Let \(x^+,x^-\) be exact top and bottom states and put

\[
r_e=\frac{\chi_e(x^+)+\chi_e(x^-)}2\in\{-1,0,1\}.
\tag{7.1}
\]

Then

\[
d(A)=\sum_ea_er_e.
\tag{7.2}
\]

For a uniformly random \(k\)-vertex set \(S\),

\[
\mathbb E_S\sum_{e\subset S}a_er_e
=
\frac{k(k-1)}{n(n-1)}\,d(A).
\tag{7.3}
\]

If \(d(A)=\varepsilon n^{3/2}\), the expected midpoint signal inside
the block is only

\[
\asymp \varepsilon\frac{k^2}{\sqrt n}.
\tag{7.4}
\]

By contrast, the unavoidable generic replacement/rounding scale for
an induced \(k\)-block is

\[
\eta_T=\Theta(k^{3/2}).
\tag{7.5}
\]

Their ratio is

\[
\Theta\!\left(\varepsilon\sqrt{\frac{k}{n}}\right).
\tag{7.6}
\]

Thus no argument that merely samples a sublinear induced block and
compares its average midpoint mass with generic block discrepancy can
close the theorem.  The planted \(k^2\)-sized spike is removable, but
a density-\(n^{-1/2}\) bias spread throughout \(\Theta(n^2)\) edges
lies below the natural \(k^{3/2}\) fluctuation scale on every
sublinear block.

This explains exactly what global replacement has achieved and what
it has not:

\[
\boxed{
\begin{array}{l}
\text{localized midpoint spikes are excluded;}\\
\text{diffuse }n^{-1/2}\text{-density midpoint bias is not.}
\end{array}}
\tag{7.7}
\]

## 8. Revised target

To prove midpoint balancing from here, one needs a replacement whose
gain adds coherently across many blocks while its discrepancy cost does
not add at the generic \(n\sqrt{k}\) scale.

Two precise sufficient routes are:

1. **entropy-sensitive replacement:** apply (2.3) using a replacement
   distribution whose upper-tail exponent beats the entropy of every
   low-slack pair shell;
2. **correlated multiblock replacement:** couple the replacements of
   many induced blocks so that their diffuse midpoint gains add
   linearly while their width fluctuations cancel.

The dual obstruction (3.6)--(3.7) gives the object either route must
defeat.  A proof cannot rely only on single-edge certificates or on
the total midpoint mass.

