# Banach foreign critique: action recovery, terminal drift, and pressure

## Bottom line

The foreign mechanism tested is **action-profile compactness followed by a \(\Gamma\)-limsup recovery sequence at every index**. Unlike the usual Banach passage from a quadratic polynomial to its polarized bilinear form, the one-profile law of \((v,Tv)\) keeps the same test function in both slots and preserves the Boolean diagonal functional without a polarization constant. This is useful conceptual leverage, but it creates no new route: compactness gives only a subsequence, while exact signed recovery at every sufficiently large order is the already-isolated AR obligation.

| Candidate | Polarization/projective loss | State secretly required | Judgment |
|---|---|---|---|
| fixed-\(C\) SR | none in the statement or action observable | only \(Q\) and \(\|A\|_{op}\); recovery is missing | genuine isolated B input, not a route alone |
| terminal \(L_{\rm drift}\) | none | top-two active-face incidence, not the full histogram | strict quotient, but Paley forces the sharp constant \(1/2\) |
| zero-temperature composition | none if the parent maximum is kept joint | all \(2^{m+n}\) bridge states | full parent optimization remains |
| fixed-temperature pressure | none if kept same-spin | full Gibbs/Laplace bridge response | softening does not compress bridge selection |
| tensor/factorization substitutes | fixed real polarization or projective-output loss | independent-input tensor norm | wrong norm or fixed leading coefficient |

## 1. Exact action normalization and the foreign experiment

For a symmetric hollow signing \(A\), let \(T_A=A/\sqrt n\) in the finite probability-space normalization. Backhausz--Szegedy action profiles give
\[
 \Phi(T)=\sup_{\|v\|_\infty\le1}|\mathbb E[v(Tv)]|,
 \qquad \Phi(T_A)=\frac1{n^{3/2}}\max_x|x^{\mathsf T}Ax|
 =\frac{2Q(A)}{n^{3/2}}.                           \tag{1}
\]
If \(T_{A_j}\to T\) in action distance and \(\sup_j\|T_{A_j}\|_{2\to2}<\infty\), convergence of one-profiles plus \(L^2\) uniform integrability gives \(\Phi(T_{A_j})\to\Phi(T)\). Equation (1) never introduces independent signs \(x,y\), real polarization, or \(K_G\).

Fixed-\(C\) SR would supply asymptotically minimizing \(A_n\) with \(\|A_n\|_{op}\le C\sqrt n\). Along a liminf subsequence, action compactness then produces \(T\) with \(\Phi(T)=2\liminf b_n\). This is progress only up to subsequential compactness.

Graph-limit sampling does not finish the argument. An arbitrary limiting \(P\)-operator need not be an ordinary kernel, and sampling it does not yield a symmetric hollow matrix with every coefficient exactly \(\pm1\), a fixed spectral bound, and the same \(\Phi\). Independent sign purification leaves order-one residual variance per edge and an order-\(n^{3/2}\) Boolean ground state. Lossless correlation of that residual is the old dependent signed-realization problem. Ordinary graphons are weaker still: competitive signed graphons tend to zero at the \(n^2\) scale.

## 2. Fixed-\(C\) SR: Banach loss audit

SR itself has no hidden polarization and recovers no coset histogram. It is genuinely stronger than the proved two-parameter regularization
\[
 Q(A'_n)\le M_n+O(K^{-1/2}n^{3/2}),\qquad
 \|A'_n\|_{op}=O(K\sqrt n),                         \tag{2}
\]
because fixing the operator bound fixes a leading error in (2).

The loss appears in attempted proofs. Spectral truncation destroys unimodularity. Independent rounding pays a new \(\Theta(n^{3/2})\) chaos; controlling it through \(\ell_\infty\to\ell_1\), Grothendieck factorization, or a projective tensor norm either changes to independent inputs or pays a fixed constant. Projective control is excessive: \(\ell_1^n\widehat\otimes_\pi\ell_1^n=\ell_1^{n^2}\) sees \(\Theta(n^2)\) coefficient mass, while the symmetric injective norm sought is \(\Theta(n^{3/2})\). Thus SR remains the cleanest local statement, but action compactness plus SR terminates at all-order signed recovery.

## 3. Terminal drift and the Paley collision

For a coset root \(U\), write
\[
 z_n(U)=\frac{N_n-2r_n(U)}{n^{3/2}},\qquad
 b_n(U)=|\{e:r_n(U+e)=r_n(U)+1\}|.                 \tag{3}
\]
There is no bilinearization here. The exact identity
\[
 b_n(U)=N_n-\left|\bigcup_{g_v\le1}N_v\right|      \tag{4}
\]
reads coverage only by the top two same-spin layers. Verified examples show \((r,b)\) does not recover the coset histogram; the danger is the proposed uniform law over every root, not a hidden tensor norm.

For infinitely many square-field Paley conference roots \(C_q\), \(n=q+1\), the supplied facts are
\[
 \boxed{b_n(C_q)=0,\qquad z_n(C_q)\to\frac12.}     \tag{5}
\]
If \(L_{\rm drift}\) held with \(b_n(U)/N_n=\beta(z_n(U))+o(1)\) and a unique zero \(c\), continuity and (5) force \(c=1/2\). Deepest roots also have \(b_n=0\), so
\[
 \frac{M_n}{n^{3/2}}\longrightarrow\frac12.        \tag{6}
\]
Thus terminal drift is information-theoretically compressed but contains the full sharp lower theorem \(\liminf b_n\ge1/2\). Its strongest falsifier is a scalable non-Paley terminal family with \(z\to z_*\ne1/2\), or two same-\(z\) families with separated \(b/N\). The exact order-ten dead end at \(z=15/(10\sqrt{10})\) already displays the finite metastability mechanism.

Standard action convergence cannot rescue (4). A one-edge change vanishes in action distance, whereas \(b/N\) is an empirical law of microscopic directional derivatives and constant-size energy gaps. A two-scale edge-marked tangent-action topology would be needed; if it retains every active cut and margin, it merely recovers the top-two incidence that (4) was meant to compress.

## 4. Pressure and composition

For a block completion,
\[
 Q\!\begin{pmatrix}A&D\\D^{\mathsf T}&B\end{pmatrix}
 =\max_{x,y}\bigl(|H_A(x)+H_B(y)|+|x^{\mathsf T}Dy|\bigr).       \tag{7}
\]
Keeping (7) intact pays no polarization loss but keeps the whole hard object: one bridge \(D\) must answer every \((x,y)\). Paying \(\|D\|_{\infty\to1}\) separately leaves a leading \(\Theta(\sqrt{mn(m+n)})\) channel. Grothendieck or vector-valued factorization adds a fixed constant and cannot correlate that channel with both child landscapes. The \(2/3\)-power, squared-cap, and compressed-lift scalarizations therefore leave the full parent response.

At fixed temperature, \(\log\overline Z_n(A,\beta)\) replaces the maximum by a complete same-spin Laplace sum. Exact bridge minimization still ranges over all sign bridges and every parent Gibbs state. One \(\beta\) is not literally the histogram, but the universally quantified completion has the same backward-response state, and all \(\beta\) determine the absolute-energy distribution. Random bridges contract child temperatures and leave a linear defect; Guerra--Toninelli does not survive the outer deterministic minimum.

Full tensor-power regularization is no cure. The \(\varepsilon\)-input to \(\pi\)-output theorem controls \(\|A\|_{\infty\to1}\), with
\[
 Q(A)\le\tfrac12\|A\|_{\infty\to1}\le2Q(A),       \tag{8}
\]
only to a root-exponential rate on one power subsequence. The fixed interval in (8), \(N^{o(1)}\) prefactor, and lack of all-order minimum control are fatal.

## 5. Strongest failure and strongest possible rescue

**Strongest failure reason.** At the \(n^{3/2}\) scale, topology gives compactness but forgets exact sign residuals and microscopic terminal derivatives. Recovering them requires the missing all-order signed realizer, exact joint bridge/Gibbs response, or top-two active-face incidence. Banach factorization exposes these only after fixed polarization/projective loss. Paley additionally forces the terminal quotient to prove the sharp constant \(1/2\).

**Strongest possible rescue theorem.** Given fixed-\(C\) SR, the minimal lossless rescue is:

> For every self-adjoint action limit \(T\) of \(C\)-spectrally bounded competitive signings, there exist at every sufficiently large order \(m\) symmetric hollow sign matrices \(B_m\) with \(\|B_m\|_{op}\le C'\sqrt m\) and \(\limsup_m\Phi(T_{B_m})\le\Phi(T)\).

Applied to a liminf action limit, (1) gives \(\limsup b_m\le\liminf b_n\). This signed \(\Gamma\)-limsup theorem needs neither polarization nor full action-profile recovery, only objective recovery. No existing action, graphon, Banach-factorization, or tensor-power theorem supplies exact sign realizers with this lossless every-order conclusion. It is the cleanest possible rescue, but presently another form of the missing cross-order obligation, not a new proved route.
