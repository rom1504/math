# Endpoint-shell large deviations under spectral regularization

## 1. Setup and normalization

Let \(A\) be a symmetric \(n\times n\) signing with zero diagonal, and write

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j=\frac12x^\top A x,
\qquad x\in\{\pm1\}^n.
\]

Set

\[
P(A)=\max_x H_A(x),\qquad m(A)=\min_x H_A(x).
\]

This note gives a rigorous exponential upper bound for the number of Boolean
vectors in a shell near either endpoint when \(A\) has controlled operator
norm.  It then compares the resulting exponent with the \(s=4\) traffic
threshold

\[
h_4(r,\rho)<\frac{5}{64}\frac{r^2}{1-\rho^2}.
\]

The conclusion is negative but useful: spectral regularization plus a
one-replica Rademacher-chaos large-deviation bound is quantitatively incapable
of supplying the required low endpoint entropy.  Moreover, the
delete-and-refill Pietsch regularization becomes asymptotically *less*
informative as its accuracy is sent to zero.

---

## 2. Shifted Hubbard--Stratonovich bound

Let the eigenvalues of \(A/\sqrt n\) lie in

\[
[-k_-,k_+],
\qquad k_-,k_+\ge0,
\]

and put

\[
L=k_-+k_+,\qquad v=k_-^2+\frac1n\operatorname{tr}(A^2)/n.
\]

For a signing,

\[
\operatorname{tr}(A^2)=n(n-1),
\]

so

\[
v=k_-^2+1-\frac1n.
\]

### Proposition 2.1

For every \(0<\theta<L^{-1}\),

\[
\frac1n\log
\mathbb E_x\exp\!\left(\frac{\theta}{\sqrt n}H_A(x)\right)
\le
\frac{v}{2L^2}
\bigl[-\log(1-L\theta)-L\theta\bigr].
\tag{2.1}
\]

Consequently, for every \(p\ge0\),

\[
\mathbb P_x\{H_A(x)\ge p n^{3/2}\}
\le
\exp\{-n I_{L,v}(p)\},
\tag{2.2}
\]

where

\[
\boxed{
I_{L,v}(p)
=
\frac pL-\frac{v}{2L^2}
\log\!\left(1+\frac{2Lp}{v}\right).
}
\tag{2.3}
\]

The lower tail follows by applying the same statement to \(-A\), with
\(k_-\) and \(k_+\) interchanged.

### Proof

Let \(c=-\lambda_{\min}(A)=k_-\sqrt n\) and

\[
M=A+cI\succeq0.
\]

Since \(x^\top x=n\),

\[
H_A(x)=\frac12x^\top Mx-\frac12cn.
\]

For \(\lambda=\theta/\sqrt n\), Gaussian linearization gives

\[
\begin{aligned}
\mathbb E_x e^{\lambda H_A(x)}
&=
e^{-\lambda cn/2}
\mathbb E_g\mathbb E_x
\exp\!\left(\sqrt\lambda\,g^\top M^{1/2}x\right)\\
&=
e^{-\lambda cn/2}
\mathbb E_g\prod_{i=1}^n
\cosh\!\left(\sqrt\lambda\,(M^{1/2}g)_i\right).
\end{aligned}
\]

Using \(\cosh u\le e^{u^2/2}\),

\[
\mathbb E_x e^{\lambda H_A(x)}
\le
e^{-\lambda cn/2}\det(I-\lambda M)^{-1/2},
\tag{2.4}
\]

provided \(\lambda\|M\|_{\mathrm{op}}<1\).

Let

\[
z_i=\frac{c+\lambda_i(A)}{\sqrt n}\in[0,L].
\]

Then

\[
\frac1n\sum_i z_i=k_-,
\qquad
\frac1n\sum_i z_i^2=v.
\]

Expanding the logarithm in (2.4), the linear term cancels the shift:

\[
\begin{aligned}
\frac1n\log\mathbb E_x e^{(\theta/\sqrt n)H_A(x)}
&\le
\frac12\sum_{\ell\ge2}
\frac{\theta^\ell}{\ell}
\left(\frac1n\sum_i z_i^\ell\right).
\end{aligned}
\]

Since \(0\le z_i\le L\),

\[
\frac1n\sum_i z_i^\ell\le L^{\ell-2}v,
\qquad \ell\ge2.
\]

Summing the series proves (2.1).  Chernoff optimization uses

\[
\theta_*=\frac{2p}{v+2Lp}
\]

and gives (2.3). \(\square\)

### Remark 2.2

This is a rigorous *upper* bound based on
\(\cosh u\le e^{u^2/2}\), not an assertion that the Ising quadratic form has
a Gaussian log determinant.  It therefore does not conflict with the
previous fourth-cumulant obstruction to Gaussian free-energy transfer.

---

## 3. Symmetric operator-norm corollary

Assume

\[
\|A\|_{\mathrm{op}}\le k\sqrt n.
\]

Then \(L\le2k\) and \(v\le k^2+1\).  Repeating the last series estimate
directly with \(2k\) and \(k^2+1\) gives the convenient uniform rate

\[
\boxed{
I_k(p)=
\frac p{2k}
-
\frac{k^2+1}{8k^2}
\log\!\left(1+\frac{4kp}{k^2+1}\right).
}
\tag{3.1}
\]

Thus

\[
\mathbb P_x\{|H_A(x)|\ge p n^{3/2}\}
\le 2e^{-nI_k(p)}.
\tag{3.2}
\]

For a signing, the Frobenius identity forces

\[
k\ge\sqrt{1-\frac1n}.
\tag{3.3}
\]

Also the spectral inequality gives

\[
\max_x|H_A(x)|
\le\frac n2\|A\|_{\mathrm{op}}
\le\frac k2 n^{3/2}.
\tag{3.4}
\]

---

## 4. Endpoint-shell entropy

Suppose

\[
P(A)=p_+n^{3/2},
\qquad
-m(A)=p_-n^{3/2}.
\]

For \(u\ge0\), define the upper and lower endpoint shells

\[
\mathcal E_+(u)
=
\{x:P(A)-H_A(x)\le u n^{3/2}\},
\]

\[
\mathcal E_-(u)
=
\{x:H_A(x)-m(A)\le u n^{3/2}\}.
\]

From (3.2),

\[
\boxed{
\frac1n\log|\mathcal E_+(u)|
\le
\log2-I_k((p_+-u)_+)+o(1),
}
\tag{4.1}
\]

\[
\boxed{
\frac1n\log|\mathcal E_-(u)|
\le
\log2-I_k((p_--u)_+)+o(1).
}
\tag{4.2}
\]

The same statements with \(I_{L,v}\) use the actual one-sided spectral
edges and can be sharper.

---

## 5. An intrinsic ceiling of this method

The rate \(I_k(p)\) is increasing in \(p\), while (3.4) restricts
\(p\le k/2\).  Therefore this argument can never certify more endpoint
rarity than

\[
I_k(k/2)
=
\frac14
-
\frac{k^2+1}{8k^2}
\log\!\left(1+\frac{2k^2}{k^2+1}\right).
\tag{5.1}
\]

At the conference scale \(k=1\),

\[
I_1(1/2)
=
\frac14-\frac14\log2
=0.07671320486\ldots.
\]

As \(k\to\infty\),

\[
I_k(k/2)
\longrightarrow
\frac14-\frac18\log3
=0.1126734639\ldots.
\]

In particular, throughout the entire spectrally allowed range this
one-replica method remains far below \(\log2\).  It cannot prove that an
exact endpoint shell is subexponential: even in its best limiting regime,
the certified entropy upper bound is no smaller than

\[
\log2-\left(\frac14-\frac18\log3\right)
=0.5804737\ldots.
\tag{5.2}
\]

This is a limitation of the certificate, not a lower bound on the true
endpoint entropy.

For a representative near-optimal scale \(p=0.4\), \(k=1\),

\[
I_1(0.4)
=0.2-\frac14\log1.8
=0.053053\ldots,
\]

so (4.1) only yields endpoint entropy at most
\(0.64009\ldots\) per vertex.

---

## 6. Comparison with the traffic threshold

For the exact-sum \(s\times s\) block at \(s=4\), the best nontrivial
quadratic coefficient is

\[
\frac1{\kappa_4}=\frac5{64}=0.078125.
\]

At zero overlap, a shell argument would need

\[
h_4(r,0)<\frac5{64}r^2.
\tag{6.1}
\]

At \(r=0\), the right side vanishes, so one needs a subexponential exact
endpoint family.  Equations (5.1)--(5.2) show that spectral LDP alone
cannot certify anything close to this.

For the relevant near-endpoint range \(0\le r\le p\le1/2\),

\[
\frac5{64}r^2\le\frac5{256}=0.01953125,
\]

whereas the best entropy bound obtainable from (4.1) remains above
\(0.58\).  Hence there is no positive-slack window in the natural
near-conference range where this estimate verifies (6.1).

The obstacle is exponential-scale dependence among near-ground states:
one-replica tail control only counts the rarity of a randomly sampled
spin, while the traffic argument needs structural control of an entire
correlated endpoint family.

---

## 7. Pietsch deletion and refill

Deleting \(\varepsilon n\) vertices can leave a principal block with

\[
\|A_{R,R}\|_{\mathrm{op}}
\le k_\varepsilon\sqrt n,
\qquad
k_\varepsilon=O(\varepsilon^{-1}).
\]

One must not apply (4.1) to this block after simply conditioning the
deleted spins: the cross edges induce a linear field

\[
h=A_{R,T}z
\]

whose size is not controlled by the principal-block operator norm.

The legitimate route is the previously established delete-and-refill
regularization: for fixed \(\varepsilon\), replace the deleted part and
obtain a full-order signing \(A'\) satisfying

\[
\|A'\|_{\mathrm{op}}\le k_\varepsilon\sqrt n,
\]

while changing the normalized width by at most

\[
O(\sqrt\varepsilon).
\]

The endpoint LDP applies to \(A'\), not automatically to the original
exact minimizer.

Unfortunately, for fixed \(p=O(1)\) and large \(k\),

\[
I_k(p)=\frac{p^2}{k^2}+O_p(k^{-3}).
\tag{7.1}
\]

Therefore

\[
I_{k_\varepsilon}(p)=O_p(\varepsilon^2).
\tag{7.2}
\]

Sending \(\varepsilon\to0\) removes the width perturbation but makes the
rarity exponent vanish.  Balancing an \(O(\sqrt\varepsilon)\) objective
error against an \(O(\varepsilon^2)\) rarity exponent cannot produce the
uniform low-shell entropy required by traffic compression.

---

## 8. Final verdict and next missing lemma

The shifted log-determinant argument provides a clean, explicit, and
fully rigorous endpoint-shell LDP for spectrally regular signings.  It
also proves a useful no-go statement:

> Operator-norm regularization plus any one-replica quadratic-form tail
> estimate of this type cannot establish the low-slack traffic criterion.

Progress now requires information absent from the one-spin law, such as:

1. a two-replica or cluster LDP for endpoint configurations;
2. a tangent-cone argument forcing covariance degeneracy of a large
   endpoint family;
3. a block-replacement theorem that converts such degeneracy into an
   improvement of the signing.

The next branch pursues precisely the second and third alternatives.
