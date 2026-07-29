# Projection-sign images and the cap-stable entropy problem

## 1. Exact projection-sign inverse lemma

Let \(U=C/\sqrt{n-1}\) be a symmetric conference involution,
\[
P_+=\frac{I+U}{2},\qquad P_-=\frac{I-U}{2},
\]
and define
\[
T(x)=\operatorname{sign}(P_+x).
\]
Ties may be resolved in favor of the current coordinate.

Put
\[
r(x)=\frac{x^\top Ux}{n}.
\]
If \(T(x)_i\ne x_i\), then
\[
x_i(P_+x)_i\le0.
\]
Since \(P_-x=x-P_+x\),
\[
x_i(P_-x)_i=1-x_i(P_+x)_i\ge1.
\]
Every changed coordinate therefore contributes at least one to
\(\|P_-x\|^2\), and
\[
\boxed{
d_H(x,T(x))
\le\|P_-x\|^2
=\frac{1-r(x)}2\,n.} \tag{1}
\]
Equivalently,
\[
\frac{x^\top Ux}{2n}
\le\frac12-\frac{d_H(x,T(x))}{n}. \tag{2}
\]

Moreover, \(T\) does not decrease projection energy.  Indeed,
\[
T(x)^\top P_+x=\|P_+x\|_1\ge x^\top P_+x,
\]
while Cauchy--Schwarz gives
\[
T(x)^\top P_+x
\le\sqrt{T(x)^\top P_+T(x)}\,
     \sqrt{x^\top P_+x}.
\]
Thus
\[
T(x)^\top P_+T(x)\ge x^\top P_+x. \tag{3}
\]

Equations (1)--(3) give the desired linear supporting-hyperplane profile
without guessing square waves, finite-cyclic gadgets, or multi-harmonic
templates.  The natural template library is the high-energy part of
\(\operatorname{Range}(T)\), or of the fixed-point set \(T(y)=y\).

## 2. Why the energy threshold is essential

The full projection-sign image and fixed-point sets are exponentially large.
Exhaustive enumeration on the antipodal cube gives:

| Paley order \(p\) | cube size \(2^{p-1}\) | stable \(T(y)=y\) | distinct \(T\)-images |
|---:|---:|---:|---:|
| 5 | 16 | 11 | 11 |
| 13 | 4096 | 586 | 703 |
| 17 | 65536 | 13397 | 15165 |

For \(p=17\), the stable-state energy histogram is:

| \(H\) | \(-8\) | \(-4\) | \(0\) | \(4\) | \(8\) | \(12\) | \(16\) | \(20\) | \(24\) | \(28\) | \(32\) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| count | 68 | 136 | 358 | 884 | 1768 | 2992 | 3281 | 2380 | 1224 | 272 | 34 |

The top \(34\) states form one affine orbit, while roughly twenty percent of
the entire cube is stable.  Any theorem counting all fixed points is
therefore unusable.  The relevant complexity is
\[
\Sigma_n(\varepsilon)
=\frac1n\log
\#\left\{y:T(y)=y,\ 1-r(y)\le\varepsilon\right\}. \tag{4}
\]
The precise live target is
\[
\boxed{
\lim_{\varepsilon\downarrow0}\limsup_{n\to\infty}
\Sigma_n(\varepsilon)=0.} \tag{5}
\]
An analogous statement for high-energy \(T\)-images is equally sufficient.

## 3. Local fields of a cap-near state

Write
\[
\ell_i=y_i(Uy)_i.
\]
Orthogonality of \(U\) gives the exact identities
\[
\frac1n\sum_i\ell_i=r(y),\qquad
\frac1n\sum_i\ell_i^2=1.
\]
Therefore
\[
\boxed{
\frac1n\sum_i(\ell_i-1)^2=2(1-r(y)).} \tag{6}
\]
If \(1-r(y)\le\varepsilon\), all but
\[
\frac{2\varepsilon}{\tau^2}n
\]
coordinates have \(|\ell_i-1|\le\tau\).  Stability additionally requires
\(\ell_i>-1\).  Thus a cap-near stable state has almost every switched row
sum within \(O(\tau\sqrt n)\) of \(+\sqrt n\).

This is a useful encoding input, but by itself it does not bound the number
of states: a quantitative Paley/arithmetic constraint is still needed to
encode the exceptional coordinates or the near-\(+\) eigenspace pattern.

## 4. Branching finite-cyclic gadgets do not refute (5)

Recursive \(\ell\)-adic gadgets can have exponentially many internal
branches as a function of their period \(q\).  This invalidates any claim
that a **fixed-defect** cap family must be finite or polynomial.

However, the known gadget-to-Paley transfer fixes \(q\) first and then sends
\(p\to\infty\).  If a period-\(q\) step word is sampled on \(p\) coordinates,
the ideal and sampled realizations differ only in \(O(q)\) boundary cells.
For the Paley core, \(\|C_p\|_{\rm op}=\sqrt p\).  If two Boolean vectors
differ in \(h\) coordinates, then
\[
\begin{aligned}
|x^\top C_px-y^\top C_py|
&\le
\|x-y\|_2\,\|C_p(x+y)\|_2\\
&\le4p\sqrt h,
\end{aligned}
\]
and hence
\[
\left|
\frac{x^\top C_px-y^\top C_py}{p^{3/2}}
\right|
\le4\sqrt{\frac hp}. \tag{7}
\]
With \(h=O(q)\), the uncontrolled transfer error is
\[
O\left(\sqrt{\frac qp}\right). \tag{8}
\]

Consequently a linearly growing period \(q=cp\) carries an \(O(\sqrt c)\)
error under the currently proved transfer.  To certify cap defect at most
\(\varepsilon\), this architecture requires
\[
\frac qp=O(\varepsilon^2). \tag{9}
\]
If the number of branching gadgets is \(\exp(\kappa q)\), their entropy per
Paley coordinate is then at most
\[
\kappa\frac qp=O(\kappa\varepsilon^2)\longrightarrow0. \tag{10}
\]

Thus recursively branching finite-period gadgets are consistent with
\(\Sigma(\varepsilon)\to0\).  They may force \(\Sigma(\varepsilon)>0\) at
every fixed \(\varepsilon\), but do not presently produce a positive limiting
entropy as \(\varepsilon\downarrow0\).

To refute (5), one would need a direct finite-\(p\) construction with
positive entropy rate and vanishing Paley defect, avoiding the
fixed-period step-sampling loss (8).

## 5. Current frontier

The projection-sign lemma reduces arithmetic template classification to a
single entropy statement:

> Count stable Paley sign patterns whose switched row sums are
> \(+\sqrt p+o(\sqrt p)\) in mean square, and prove their entropy rate
> vanishes with the allowed mean-square defect.

Square waves, multi-harmonic \(p=17\) fixed points, and finite-cyclic gadgets
all fit automatically into this formulation.  The exact lemma (1) supplies
the linear rigidity needed by correlated sparse perturbation once (5) is
available.

## 6. Paley-only predictor and cyclotomic reduction

There is a useful exact refinement for the cyclic Paley core.  Let
\(C_{st}=\chi(s-t)\) on \(\mathbb F_p\), \(U=C/\sqrt p\), and let
\(P_Q,P_N,P_0\) be the Fourier projections onto the nonzero QR, NQR, and
zero frequencies.  The multipliers of \(U\) are \(+1,-1,0\), respectively,
so every sign vector satisfies
\[
\boxed{\|y-Uy\|_2^2
=4\|P_Ny\|_2^2+\|P_0y\|_2^2.} \tag{11}
\]
Since \(U_{ii}=0\), \((Uy)_i\) is a leave-one-out linear predictor for
\(y_i\).  If \(\mathcal S\) is a translation-invariant family satisfying
the right side of (11) at most \(4\varepsilon p\), and \(Y\) is uniform on
\(\mathcal S\), translation invariance gives
\[
\mathbb E(Y_i-(UY)_i)^2\le4\varepsilon
\]
for every \(i\).  Therefore
\[
\Pr\{\operatorname{sign}((UY)_i)\ne Y_i\}\le4\varepsilon,\qquad
H(Y_i\mid Y_{-i})\le h_2(4\varepsilon). \tag{12}
\]
This is not by itself a joint entropy bound: the even-parity distribution
has all leave-one-out conditional entropies zero and joint entropy \(p-1\).
A Paley-specific causal reconstruction or cluster theorem is still needed.

The arithmetic form of that missing theorem is particularly clean.  Put
\[
\alpha_y=\sum_{j=0}^{p-1}y_j\zeta_p^j\in\mathbb Z[\zeta_p].
\]
All nonzero DFT coefficients are Galois conjugates:
\[
\widehat y(k)=p^{-1/2}\sigma_{-k}(\alpha_y),
\qquad \sigma_k(\zeta_p)=\zeta_p^k. \tag{13}
\]
Hence
\[
\|P_Ny\|_2^2\le\varepsilon p
\quad\Longleftrightarrow\quad
\sum_{k\in N}|\sigma_k(\alpha_y)|^2\le\varepsilon p^2. \tag{14}
\]
Thus (5) is equivalently a lattice-point estimate for Littlewood
cyclotomic integers whose NQR half of the Minkowski embedding is small.

Chang/large-spectrum bounds control only the coefficients of normalized
size at least a fixed \(\tau\): there are at most \(\tau^{-2}\) of them and
only \(p^{O_\tau(1)}=\exp(o(p))\) possible frequency supports.  They leave
uncontrolled a flat component of \(\Theta(p)\) coefficients of size
\(\Theta(p^{-1/2})\).  Controlling that component is precisely the missing
Paley inverse theorem.

An anisotropic theta-series relaxation and the full details of the
cyclotomic formulation are recorded in
`paley_stable_tail_cyclotomic_entropy.md`.
