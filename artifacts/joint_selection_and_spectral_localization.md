# Joint selection and spectral localization

Let \(A\) be a symmetric zero-diagonal sign matrix and put
\[
Q=Q(A)=\max_{x\in\{\pm1\}^n}|x^\top A x|.
\]
For the opposite-orientation field-plus-spin law, let
\[
R=\sigma X^\top A X,\qquad S=X^\top A^2X=\|AX\|_2^2,
\]
where \(\sigma\) is uniform on \(\{\pm1\}\).

## 1. What follows without any new concentration theorem

Suppose along a competing sequence
\[
Q=(c_*+o(1))n^{3/2},\qquad
\mathbb E R\ge(c_*-o(1))n^{3/2}.
\]
Since \(R\le Q\) pointwise,
\[
\mathbb E(Q-R)=o(n^{3/2}).
\]
Choose a sample with \(Q-R\le\mathbb E(Q-R)\). Then
\[
R=Q-o(n^{3/2})
\]
and Cauchy--Schwarz gives
\[
S\ge \frac{R^2}{n}
=(c_*^2-o(1))n^2.
\]
Applying the exact negative-field correction to the oriented matrix
\(\sigma A\) gives
\[
L_-^2\le Q(Q-R),
\]
for all sufficiently large \(n\), and consequently
\[
L_-=o(n^{3/2}).
\]
Thus the literal simultaneous target \(S=\Omega(n^2)\) is automatic.
The real unresolved target is the sharp coefficient
\(S\ge(1-o(1))n^2\).

## 2. Uniform-integrability selection lemma

The polarization identity implies
\[
\max_{x,y\in\{\pm1\}^n}|x^\top Ay|\le 2Q.
\]
Indeed, with \(u=(x+y)/2\) and \(v=(x-y)/2\),
\[
x^\top Ay=u^\top Au-v^\top Av,
\]
and randomized Boolean completion gives
\(|u^\top Au|,|v^\top Av|\le Q\). Hence, for every Boolean \(x\),
\[
S=\|Ax\|_2^2
\le (n-1)\|Ax\|_1
\le 2(n-1)Q.
\]

More importantly, if the rounded variables \(S/n^2\) are uniformly
integrable, then the exact average theorem
\[
\mathbb E S\ge n(n-1)
\]
and \(\mathbb E(Q-R)=o(n^{3/2})\) force a single sample with
\[
R=Q-o(n^{3/2}),\qquad
S\ge(1-o(1))n^2,\qquad
L_-=o(n^{3/2}).
\]
For a direct proof, truncate \(S\) at \(Kn^2\), discard the event
\(\{Q-R>\varepsilon n^{3/2}\}\), and then let
\(n\to\infty\), \(\varepsilon\downarrow0\), and \(K\to\infty\).

A particularly clean sufficient condition is
\[
\|A\|_{\mathrm{op}}=O(\sqrt n),
\]
because then \(S\le n\|A\|_{\mathrm{op}}^2=O(n^2)\) pointwise.

Quantitatively, put
\[
\eta=\frac{\mathbb E(Q-R)}{n^{3/2}},\qquad
U=\frac{\|A\|_{\mathrm{op}}^2}{n}.
\]
If \(U\eta=o(1)\), choosing the energy-deficit cutoff
\(a=\sqrt{U\eta}\) gives a sample satisfying
\[
Q-R\le a n^{3/2},\qquad
S\ge\bigl(1-o(1)-\sqrt{U\eta}\bigr)n^2.
\]

## 3. Why the two averaged inequalities alone are insufficient

Fix \(c>0\), put \(Q=cn^{3/2}\), and set
\[
p=\frac{1-c^2}{2c\sqrt n-c^2}.
\]
Consider the abstract two-point law
\[
(R,S,L_-)=
\begin{cases}
(Q,c^2n^2,0),&\text{with probability }1-p,\\
(0,2cn^{5/2},Q),&\text{with probability }p.
\end{cases}
\]
It has
\[
\mathbb ER=Q-\Theta(n)=Q-o(n^{3/2}),
\qquad
\mathbb ES=n^2,
\]
but every near-\(Q\) sample has only \(S=c^2n^2\). It also saturates
the known scalar constraints
\[
S\ge R^2/n,\qquad S\le2nQ,\qquad
L_-^2\le Q(Q-R).
\]
The corresponding local-field profiles can be modeled by
\(r_i=c\sqrt n\) in the good state, and by \(c\sqrt n\) entries \(+n\),
\(c\sqrt n\) entries \(-n\), and zeros in the rare state. This is not
claimed to come from one actual signing; it proves that moment
bookkeeping alone cannot establish coefficient-one selection.

## 4. Exact spectral bound and anomalous-mode localization

Let \(\lambda=\|A\|_{\mathrm{op}}\), and replace \(A\) by \(-A\) if
needed so that \(Av=\lambda v\) for a unit vector \(v\).
The bilinear polarization bound above also gives
\[
\lambda^2\le
\max_{x,y\in\{\pm1\}^n}|x^\top Ay|
\le2Q.
\]
For the first inequality, take \(x=\operatorname{sgn}v\), and take a
random sign vector \(y\) with
\(\mathbb Ey=v/\|v\|_\infty\). Then
\[
\mathbb E\,x^\top Ay
=\lambda\frac{\|v\|_1}{\|v\|_\infty}
\ge\lambda^2,
\]
because the eigenvector equation at a largest coordinate gives
\(\lambda\|v\|_\infty\le\|v\|_1\).
Thus every competing sequence satisfies
\[
\|A\|_{\mathrm{op}}=O(n^{3/4}).
\]

There is also a strong localization theorem. For any \(\theta>1\), set
\[
s^2=\frac{\lambda}{\theta Q},\qquad
T=\{i:|v_i|>s\},\qquad
\alpha=\sum_{i\in T}v_i^2.
\]
Then
\[
\boxed{
|T|\le\frac{\theta Q}{\lambda},
\qquad
\alpha\ge\frac{1-\theta^{-1}}3.
}
\]

Proof: write \(u=v_T\), \(w=v_{T^c}\). The eigenvector identity gives
\[
w^\top Aw
=\lambda(1-2\alpha)+u^\top Au
\ge\lambda(1-3\alpha).
\]
On the other hand, \(\|w\|_\infty\le s\), so \(w/s\in[-1,1]^n\).
Randomized Boolean completion therefore gives
\[
|w^\top Aw|\le Qs^2=\frac{\lambda}{\theta}.
\]
Combining the two inequalities proves the mass bound, while the
cardinality bound follows from \(\sum_i v_i^2=1\).

Consequently, if
\[
\lambda=L_n\sqrt n,\qquad L_n\to\infty,
\qquad Q=O(n^{3/2}),
\]
then a set of \(O(n/L_n)=o(n)\) vertices carries a fixed positive
fraction of the top eigenvector's mass. With \(\theta=2\), the explicit
figures are
\[
|T|\le\frac{2Q}{\lambda},
\qquad
\|v_T\|_2^2\ge\frac16.
\]
Thus failure of the regular-branch joint-selection theorem forces a
genuinely localized spectral anomaly, supplying a concrete target for
vertex peeling or conditioning.

## 5. The strongest capped-field consequence

Here is the exact conclusion available after a hypothetical peeling step.
Switch by a near-maximizing witness and write
\[
D=\operatorname{diag}(x)A\operatorname{diag}(x),\qquad
r=D\mathbf1,\qquad
q=\mathbf1^\top D\mathbf1.
\]
Assume, for fixed \(K,H\), that
\[
\|D\|_{\rm op}\le(K+o(1))\sqrt n,\qquad
q=Q-o(n^{3/2}),
\]
\[
\sum_i r_i^2=(s+o(1))n^2,\quad s\ge1,\qquad
\max_i|r_i|\le(H+o(1))\sqrt n,
\]
and that the negative mass is \(o(n^{3/2})\). The cap makes the squared
negative mass \(o(n^2)\), so replacing \(r\) by \(r_+\) changes its
normalized second and third spectral moments by \(o(1)\).

Put
\[
c=\frac{q}{n^{3/2}},\qquad
t=\frac{r^\top Dr}{n^{5/2}}.
\]
The spectral measure of \(\mathbf1\) has first three moments \(c,s,t\)
on \([-K,K]\). The polynomial majorant
\[
z^3\le (K+2a)z^2-(2Ka+a^2)z+Ka^2
\]
follows from \((K-z)(z-a)^2\ge0\). Optimizing in \(a\) gives
\[
t\le M_K(c,s)
:=Ks-\frac{(s-Kc)^2}{K-c}.
\]

For \(0\le\alpha\le1/H\), choose independent cut probabilities
\[
p_i=\frac{\alpha(r_i)_+}{\sqrt n}.
\]
Every cut has signed weight at most
\((Q+q)/4=(c/2+o(1))n^{3/2}\), whereas its expected weight is
\[
\left(\alpha s-\alpha^2t+o(1)\right)n^{3/2}.
\]
Consequently every limiting profile obeys
\[
\boxed{
\frac c2\ge
\max_{0\le\alpha\le1/H}
\left[\alpha s-\alpha^2M_K(c,s)\right],
\qquad
c\ge\frac{s}{H}.
}
\]
This is the strongest conclusion of the present capped-profile
calculation.

For the coefficient-one case \(s=1\), abbreviate
\[
M_K(c)=K-\frac{(1-Kc)^2}{K-c}.
\]
Then
\[
\max_{0\le\alpha\le1/H}
(\alpha-M_K(c)\alpha^2)
=
\begin{cases}
\dfrac1{4M_K(c)},&H\le2M_K(c),\\[6pt]
\dfrac1H-\dfrac{M_K(c)}{H^2},&H>2M_K(c).
\end{cases}
\]

At \(c=c_*=0.672986728863\ldots\), let
\[
M_*=M_K(c_*).
\]
For
\[
1\le K<K_{\rm crit}
=1.022079887507\ldots,
\qquad
2c_*M_*<1,
\]
the value \(c_*\) is impossible whenever
\[
H<h_+(K)
:=\frac{1+\sqrt{1-2c_*M_*}}{c_*}.
\]
Numerically,
\[
h_+(1)=1.941916296158\ldots,
\]
and \(h_+(K)\) decreases to
\[
\frac1{c_*}=1.485913402318\ldots
\]
at \(K=K_{\rm crit}\). For \(K\ge K_{\rm crit}\), this calculation
improves on \(c_*\) only under the elementary cap
\(H<1/c_*\).

Thus a spectrally flat core with \(K=1+o(1)\) cannot saturate the current
lower constant unless it contains a positive local field of at least
\((1.9419-o(1))\sqrt n\). The universal-positive-vertex construction
shows why such heavy coordinates must be peeled rather than merely
bounded by a scalar inequality.
