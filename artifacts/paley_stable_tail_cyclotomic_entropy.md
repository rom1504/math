# Paley stable-tail entropy: predictor and cyclotomic reductions

## Status

This note records two exact reductions for the proposed Paley stable-tail
theorem.  Neither reduction by itself proves the desired entropy estimate.
They identify the additional theorem that is genuinely missing.

Throughout, \(p\equiv1\pmod 4\) is prime, \(\chi\) is the Legendre symbol, and
\[
  C_{st}=\chi(s-t),\qquad s,t\in\mathbb F_p,
\]
with \(\chi(0)=0\).  In the unitary DFT basis,
\[
  U:=\frac C{\sqrt p}
\]
has multiplier \(+1\) on the nonzero quadratic-residue frequencies, \(-1\)
on the nonresidue frequencies, and \(0\) at frequency zero.  Write
\(P_Q,P_N,P_0\) for the three corresponding orthogonal projections.

The order-\(p+1\) Paley conference matrix removes the zero mode and gives an
exact symmetric involution.  All statements below have an exact
order-\(p+1\) version; the cyclic order-\(p\) formulation is more useful for
the arithmetic reduction, and the one-dimensional zero mode is harmless at
the \(p\)-scale.

## 1. Exact leave-one-out predictor

For a sign vector \(y\in\{\pm1\}^p\),
\[
 \|y-Uy\|_2^2
 =4\|P_Ny\|_2^2+\|P_0y\|_2^2. \tag{1}
\]
This is immediate in the DFT basis from the multipliers \(1,-1,0\).
Moreover \(U_{ii}=0\), so \((Uy)_i\) is a linear function of all signs other
than \(y_i\).

Let \(\mathcal S\subseteq\{\pm1\}^p\) be invariant under cyclic translations,
let \(Y\) be uniform on \(\mathcal S\), and suppose every \(y\in\mathcal S\)
satisfies
\[
 4\|P_Ny\|_2^2+\|P_0y\|_2^2\le 4\varepsilon p. \tag{2}
\]
Translation invariance and (1) give, for every coordinate \(i\),
\[
 \mathbb E\bigl(Y_i-(UY)_i\bigr)^2\le4\varepsilon. \tag{3}
\]
If \(\widehat Y_i=\operatorname{sign}((UY)_i)\), then on the event
\(\widehat Y_i\ne Y_i\) one has
\(\lvert Y_i-(UY)_i\rvert\ge1\).  Consequently
\[
 \Pr(\widehat Y_i\ne Y_i)\le4\varepsilon,\qquad
 H(Y_i\mid Y_{-i})\le h_2(4\varepsilon), \tag{4}
\]
where the last inequality is in bits and \(h_2\) is binary entropy.

For the order-\(p+1\) conference matrix the same argument is exact without a
zero-mode term.  Averaging over coordinates gives (3)--(4), and the single
distinguished Paley vertex affects normalized quantities by \(o(1)\).

### Why this is not yet an entropy bound

Small leave-one-out conditional entropies do **not** imply small joint
entropy.  The even-parity distribution on \(\{\pm1\}^p\) has
\[
 H(Y_i\mid Y_{-i})=0\quad\hbox{for every }i,
 \qquad H(Y)=p-1.
\]
Thus (4) must be supplemented by a causal reconstruction theorem, a
cluster decomposition, or a Paley-specific assertion ruling out
high-dimensional parity-like dependencies.  Shearer's inequality has the
wrong direction for deriving an upper bound on \(H(Y)\) from (4).

## 2. Stationary covariance formulation

For the same translation-invariant uniform measure, the covariance
\[
 R=\mathbb E[YY^\top]
\]
is circulant.  Its DFT eigenvalues are
\[
 \lambda_k=\mathbb E|\widehat Y(k)|^2
\]
under the unitary DFT normalization.  Condition (2) implies
\[
 4\sum_{k\in N}\lambda_k+\lambda_0\le4\varepsilon p. \tag{5}
\]
Hence the desired counting theorem can equivalently be viewed as a
finite-cyclic, binary stationary-process theorem:

> A translation-invariant binary process whose covariance spectral mass is
> almost entirely on the Paley QR frequencies should have entropy rate
> tending to zero.

Ordinary stationary prediction theory does not immediately apply.  The QR
and NQR sets are equidistributed and interlaced, rather than converging to
two fixed arcs of the circle, so there is no fixed macroscopic spectral gap
from which to build a one-sided prediction filter.

## 3. Cyclotomic reformulation

Let \(\zeta=e^{2\pi i/p}\) and
\[
  A_y(z)=\sum_{j=0}^{p-1}y_j z^j,\qquad
  \alpha_y=A_y(\zeta)\in\mathbb Z[\zeta].
\]
For every \(k\ne0\),
\[
  A_y(\zeta^k)=\sigma_k(\alpha_y),
  \qquad \sigma_k(\zeta)=\zeta^k. \tag{6}
\]
Thus all nonzero DFT coefficients are Galois conjugates of one cyclotomic
integer.  In the unitary normalization,
\[
 \widehat y(k)=p^{-1/2}\sigma_{-k}(\alpha_y). \tag{7}
\]
The negative Paley energy condition is exactly
\[
 \sum_{k\in N}\bigl|\sigma_k(\alpha_y)\bigr|^2
 \le \varepsilon p^2, \tag{8}
\]
up to the convention-dependent factor in the definition of
\(\varepsilon\).

This makes the missing result an arithmetic lattice-point theorem:

> Count Littlewood cyclotomic integers
> \(\alpha_y=\sum y_j\zeta^j\), \(y_j=\pm1\), for which the NQR half of the
> Minkowski embedding has squared norm at most \(\varepsilon p^2\).

The map \(y\mapsto\alpha_y\) is injective except for the trivial collision
between the two constant sign vectors: if \(A_y(\zeta)=A_z(\zeta)\), then
\(A_y-A_z\) is an integer multiple of
\(\Phi_p(z)=1+\cdots+z^{p-1}\).

This formulation explains why ordinary coherence does not see the relevant
structure.  The exceptional square waves are precisely cyclotomic integers
whose small and large Galois conjugates are separated by the QR/NQR
partition.

## 4. Theta/Poisson route and its obstruction

The cube can be embedded into the odd integer lattice.  Since an odd
integer vector \(z\in(2\mathbb Z+1)^p\) obeys
\(\|z\|_2^2\ge p\), with equality exactly for sign vectors, for \(s,t>0\)
one obtains the rigorous relaxation
\[
\begin{split}
 \#\{y\in\{\pm1\}^p:\|P_Ny\|_2^2\le\varepsilon p\}
 \le{}& e^{sp+t\varepsilon p}\\
 &\times\sum_{z\in(2\mathbb Z+1)^p}
 e^{-s\|z\|_2^2-t\|P_Nz\|_2^2}. \tag{9}
\end{split}
\]
The right-hand side is an anisotropic theta series for
\(sI+tP_N\).  Poisson summation turns its continuous determinant term into
a strong volume bound, but also produces a dual theta series containing
integer vectors close to the opposite Paley eigenspace.  Those dual
near-eigenvectors are exactly the resonant templates that the desired
inverse theorem must classify.  Dropping the dual theta terms is invalid.

Thus Poisson summation does not remove the hard part; it cleanly splits it
into:

1. a bulk determinant contribution, which has the desired small-volume
   behavior as \(\varepsilon\downarrow0\); and
2. a resonant dual-lattice contribution, which requires a cluster/template
   theorem.

## 5. Large-spectrum methods: exact limit of the present argument

For \(a_k=p^{-1}\sum_j y_j\zeta^{-jk}\), Parseval gives
\(\sum_k|a_k|^2=1\).  A threshold
\[
 L_\tau(y)=\{k:|a_k|\ge\tau\}
\]
has \(|L_\tau(y)|\le\tau^{-2}\).  Consequently, the genuinely spiky part of
the spectrum has only \(O_\tau(1)\) frequencies, and its frequency support
has only \(p^{O_\tau(1)}=\exp(o(p))\) possible choices.

This does **not** control the remaining flat spectral mass.  Chang-type
large-spectrum lemmas only encode \(L_\tau(y)\); they do not prevent
\(\Theta(p)\) coefficients of size \(\Theta(p^{-1/2})\) from carrying most
of the QR energy.  Excluding that flat scenario is essentially the
Paley-to-random-orthogonal transfer theorem, which remains unproved.

## 6. Precise surviving target

A sufficient theorem is the following two-phase inverse statement.  For
every \(\delta>0\), there should be a finite-complexity template family
\(\mathcal T_{\delta,p}\) with
\[
 \log|\mathcal T_{\delta,p}|=o_\delta(p)
\]
such that every sign vector with
\(\|P_Ny\|_2^2\le\varepsilon p\) either

1. lies within \(o_\varepsilon(p)\) Hamming distance of a template in
   \(\mathcal T_{\delta,p}\), or
2. has flat spectrum and satisfies a uniform energy gap
   \(\|P_Ny\|_2^2\ge c(\delta)p\).

The square-wave and finite cyclic hierarchical gadgets must be included in
the first alternative.  The established gadget transfer estimates are
consistent with an \(o_\varepsilon(p)\) cluster entropy, but do not prove
that the template library is complete.

## Conclusion

The Paley stable-tail entropy theorem remains plausible, but it is not a
consequence of coherence, leave-one-out predictability, or a standard
large-spectrum lemma.  The exact missing input is a cyclotomic
near-conjugate inverse theorem (equivalently, control of the resonant terms
in the anisotropic theta series).
