# Banach-space / Grothendieck / tensor-norm retrieval toolkit

Status: retrieval packet checked through 2026-08-15; not a proposed solution.

Retrieval date: 2026-08-15. This is a theorem-level source packet, not an attempted solution. The motivating quantity is

\[
M_n:=\min_{a_{ij}\in\{\pm1\}}\max_{x\in\{\pm1\}^n}
 \left|\sum_{1\le i<j\le n}a_{ij}x_ix_j\right|.
\]

The directly relevant literature below gives the scale (M_n\asymp n^{3/2}), but I found no primary theorem determining a limiting constant, proving that a normalized limit exists, or optimizing the complete-support same-spin problem. Most famous sharp constants in this neighborhood belong to a different norm.

## Norm dictionary: do not merge these problems

For a symmetric zero-diagonal matrix (A), put (q_A(x)=x^{\mathsf T}Ax), with the harmless factor-of-two convention fixed once and for all. Then

* **same-spin quadratic norm:** (\|q_A\|_{\Delta}=\max_{x\in\{\pm1\}^n}|x^{\mathsf T}Ax|);
* **independent-block bilinear norm:** (\|A\|_{\infty\to1}=\max_{x,y\in\{\pm1\}^n}|x^{\mathsf T}Ay|);
* **Hilbert/vector relaxation:** (\sup_{\|u_i\|=\|v_j\|=1}|\sum A_{ij}\langle u_i,v_j\rangle|);
* **coefficient inequality:** an estimate of an (\ell_p)-norm of all coefficients by one of the preceding function norms;
* **complete-support minimization:** the coefficients are constrained to be unimodular on *every* edge and one minimizes the function norm over the signing.

Always (\|q_A\|_{\Delta}\leq\|A\|_{\infty\to1}), but the reverse passage is polarization/decoupling and can lose a constant. Grothendieck controls the last two bilinear norms, not (\|q_A\|_{\Delta}). A generic BH inequality supplies a lower bound for every coefficient array; it does not construct, classify, or optimize a complete-support signing.

## A. Direct Boolean-cube and complete-support results

### 1. Defant–Galicer–Mansilla–Mastyło–Muro (2024)

**Andreas Defant, Daniel Galicer, Martín Mansilla, Mieczysław Mastyło, Santiago Muro, “Asymptotic Insights for Projection, Gordon–Lewis, and Sidon Constants in Boolean Cube Function Spaces,” IMRN 2024, 11239–11270.** [DOI](https://doi.org/10.1093/imrn/rnae083) · [arXiv](https://arxiv.org/abs/2302.00233)

* **Theorems/normalization.** For (\mathcal S\subset2^{[N]}), (\mathcal B^N_{\mathcal S}=\operatorname{span}\{\chi_S:S\in\mathcal S\}\subset C(\{\pm1\}^N)). Theorem 3.1 identifies its projection constant exactly as
  \[
  \lambda(\mathcal B^N_{\mathcal S})=2^{-N}\sum_x\left|\sum_{S\in\mathcal S}x^S\right|.
  \]
  For fixed (d), Theorem 4.1 gives (N^{-d/2}\lambda(\mathcal B^N_{=d})\) as the Gaussian (L_1)-norm of (h_d/d!). At (d=2), (\lambda(\mathcal B^N_{=2})/N\to\sqrt{2/(\pi e)}). Proposition 5.10 gives, for sufficiently large degree-(d) supports, Sidon size of order (\sqrt{|\mathcal S|/N}), with degree-dependent constants. Its KSZ step gives signs on an arbitrary prescribed support with
  \[
  \Big\|\sum_{S\in\mathcal S}\varepsilon_S\chi_S\Big\|_\infty
  \le C\sqrt{N|\mathcal S|}.
  \]
  Thus (\mathcal S=\binom{[N]}2) gives a **complete-support, same-spin** (O(N^{3/2})) signing.
* **Proof mechanism.** Rudin averaging gives the projection formula. For all (d)-sets, monomial decomposition or Beckner’s identity turns the Walsh kernel into a Hermite polynomial of (N^{-1/2}\sum x_i); CLT plus (L_2)-uniform integrability transfers expectations. Sidon comparisons factor through Gordon–Lewis and summing norms; the upper signing is probabilistic KSZ.
* **Boundary.** The explicit constant (\sqrt{2/(\pi e)}) is for the *projection constant/all-ones kernel in (L_1)*, not the minimum (L_\infty) signing. The KSZ constant is not claimed sharp and yields no normalized limit.

### 2. Defant–Mastyło–Pérez (2019)

**Andreas Defant, Mieczysław Mastyło, Antonio Pérez, “On the Fourier Spectrum of Functions on Boolean Cubes,” Math. Ann. 374 (2019), 653–680.** [DOI](https://doi.org/10.1007/s00208-018-1756-y) · [arXiv](https://arxiv.org/abs/1706.03670)

* **Theorem/normalization.** If (f:\{\pm1\}^N\to\mathbb R) has degree at most (d), then for an absolute (C),
  \[
  \left(\sum_{|S|\le d}|\widehat f(S)|^{2d/(d+1)}\right)^{(d+1)/(2d)}
  \le C^{\sqrt{d\log d}}\|f\|_\infty.
  \]
  The tetrahedral polynomial on ([-1,1]^N) is isometric to its restriction to the cube. At (d=2), (m) unimodular coefficients force (\|f\|_\infty\gtrsim m^{3/4}); for all edges (m=\binom N2), this is the direct (\Omega(N^{3/2})) scale.
* **Proof mechanism.** Blei’s mixed-norm inequality, Boolean hypercontractivity/Khinchine, symmetric multilinearization, and a multiaffine polarization estimate; the recursion is optimized at a block size (k\asymp\sqrt{d/\log d}).
* **Boundary.** This is a universal coefficient inequality. It neither selects a signing nor supplies the best (d=2) constant for complete support. Its real cube normalization is not the complex Steinhaus/polydisc normalization.

### 3. Defant–Frerick–Ortega-Cerdà–Ounaïes–Seip (2011)

**Andreas Defant, Leonhard Frerick, Joaquim Ortega-Cerdà, Myriam Ounaïes, Kristian Seip, “The Bohnenblust–Hille Inequality for Homogeneous Polynomials Is Hypercontractive,” Ann. Math. 174 (2011), 485–497.** [DOI](https://doi.org/10.4007/annals.2011.174.1.13) · [arXiv](https://arxiv.org/abs/0904.3540)

* **Theorem/normalization.** For a complex (m)-homogeneous (P(z)=\sum_{|\alpha|=m}a_\alpha z^\alpha) on the polydisc,
  \[
  \|(a_\alpha)\|_{2m/(m+1)}\le C_m\|P\|_{\mathbb D^n},
  \]
  with the explicit bound (C_m\le(1+1/(m-1))^{m-1}\sqrt m\,2^{(m-1)/2}). The exponent is dimension-free and sharp.
* **Proof mechanism.** Blei mixed norms, Bayart/Khinchine hypercontractivity, and Harris’s polarization estimates; the same package converts coefficient estimates into Sidon bounds.
* **Boundary.** This is complex Steinhaus/polydisc and generic-coefficient theory. It does not compare directly with real Boolean same-spin extrema, and it does not solve complete-support minimization.

### 4. Diniz–Muñoz-Fernández–Pellegrino–Seoane-Sepúlveda (2014)

**Diogo Diniz, Gustavo A. Muñoz-Fernández, Daniel Pellegrino, Juan B. Seoane-Sepúlveda, “Lower Bounds for the Constants in the Bohnenblust–Hille Inequality: the Case of Real Scalars,” Proc. AMS 142 (2014), 575–580.** [DOI](https://doi.org/10.1090/S0002-9939-2013-11791-0) · [arXiv](https://arxiv.org/abs/1111.3253)

* **Theorem/normalization.** For real (m)-linear (T:(\ell_\infty^N)^m\to\mathbb R), the BH coefficient exponent is (2m/(m+1)). Their explicit recursively constructed forms give lower bounds on the optimal real constants; at (m=2), together with Littlewood’s upper estimate, the sharp real constant is (\sqrt2):
  \[
  \left(\sum_{i,j}|T(e_i,e_j)|^{4/3}\right)^{3/4}\le\sqrt2\,\|T\|.
  \]
* **Proof mechanism.** Sparse recursively signed extremal forms and exact coefficient counting.
* **Boundary.** (T(x,y)) has independent blocks and includes arbitrary rectangular support. The sharp (\sqrt2) is not a same-spin quadratic constant.

### 5. Núñez-Alarcón–Pellegrino (2015/2016)

**Daniel Núñez-Alarcón, Daniel Pellegrino, “The Optimal Constants for the Real Hardy–Littlewood Inequality for Bilinear Forms on (c_0\times\ell_p),”** [arXiv:1508.02355](https://arxiv.org/abs/1508.02355)

* **Theorem/normalization.** For (A:\ell_p\times\ell_q\to\mathbb R), (p,q\ge2), the nested Hardy–Littlewood coefficient norm uses (\lambda=pq/(pq-p-q)). In the (c_0\times\ell_p) regime the sharp constant is (C_{p,\infty}=2^{1/2-1/p}) for (p\ge p_0/(p_0-1)\approx2.18), where (\Gamma((p_0+1)/2)=\sqrt\pi/2); the endpoint recovers the sharp mixed Littlewood constant (\sqrt2).
* **Proof mechanism.** Optimal Khinchine constants plus interpolation/mixed-sum inequalities, with explicit two-variable witnesses for sharpness.
* **Boundary.** These constants depend on the ordered mixed norm and on independent sequence-space factors. They are not constants for a symmetric zero-diagonal Boolean quadratic.

### 6. Jiménez-Rodríguez–Muñoz-Fernández–Murillo-Arcila–Seoane-Sepúlveda (2015/2016)

**Pedro Jiménez-Rodríguez, Gustavo A. Muñoz-Fernández, Marina Murillo-Arcila, Juan B. Seoane-Sepúlveda, “Sharp Values for the Constants in the Polynomial Bohnenblust–Hille Inequality,” Linear Multilinear Algebra 64 (2016).** [DOI](https://doi.org/10.1080/03081087.2015.1115810) · [arXiv](https://arxiv.org/abs/1502.02173)

* **Theorem/normalization.** The extreme points of the unit ball of real 2-homogeneous polynomials on (\ell_\infty^2) are classified, and the exact two-variable real polynomial BH constant is (D_{\mathbb R,2}(2)\approx1.837373); it is obtained by maximizing a one-variable expression over the extreme family (t x^2-t y^2\pm2\sqrt{t(1-t)}xy).
* **Proof mechanism.** Krein–Milman/extreme-point reduction followed by elementary one-dimensional optimization.
* **Boundary.** This illustrates how sharp polynomial constants are extracted, but the space has only two variables and diagonal monomials are allowed. It is neither the dimension-free constant nor the tetrahedral complete graph.

## B. Polarization and decoupling: the price of replacing one spin by two

### 7. Dimant–Galicer–Rodríguez (2022)

**Verónica Dimant, Daniel Galicer, Jorge Tomás Rodríguez, “The Polarization Constant of Finite Dimensional Complex Spaces Is One,” Math. Proc. Camb. Phil. Soc. 172 (2022), 105–123.** [DOI](https://doi.org/10.1017/S030500412100013X) · [arXiv](https://arxiv.org/abs/1908.08107)

* **Theorem/normalization.** For a (k)-homogeneous polynomial (P(x)=\check P(x,\ldots,x)), let (c(k,X)) be the best (\|\check P\|\le c(k,X)\|P\|). Universally (c(k,X)\le k^k/k!), and this is sharp in general; Hilbert spaces have equality of the two norms. Theorem 1.1 says (\limsup_k c(k,X)^{1/k}=1) for every finite-dimensional **complex** (X). The real statement fails; in particular (c(2,X)=1) forces a real (X) to be Hilbertian.
* **Proof mechanism.** Approximate a finite-dimensional space by a quotient of finite-dimensional (\ell_1), use exact complex (\ell_1^d) estimates, then Stirling; the real analysis is tied to Bochnak complexification and type/cotype.
* **Boundary / negative lesson.** At fixed degree two the generic real inequality only gives (\|\check P\|\le2\|P\|). The asymptotic-in-degree complex result does not erase a degree-two real polarization loss, and it says nothing special about zero-diagonal/tetrahedral forms.

### 8. O’Donnell–Zhao (2016)

**Ryan O’Donnell, Yu Zhao, “Polynomial Bounds for Decoupling, with Applications,” CCC 2016, 24:1–24:18.** [DOI](https://doi.org/10.4230/LIPIcs.CCC.2016.24) · [arXiv](https://arxiv.org/abs/1512.01603)

* **Theorem/normalization.** For a Banach-valued multilinear polynomial (f(x)=\sum_{|S|\le k}a_Sx^S), define one-block decoupling
  \[
  \breve f(y,z)=\sum_Sa_S\sum_{i\in S}y_i z^{S\setminus\{i\}}.
  \]
  If (f) is homogeneous, (\breve f(x,x)=k f(x)). Theorem 2.9 compares convex moments/tails with (C_k=O(k)) for Gaussians, (O(k^2)) for Rademachers, and (O(k^{3/2})) for homogeneous Rademachers; Corollary 2.12 gives (\|\breve f\|_\infty\le O(k^2)\|f\|_\infty). Full homogeneous decoupling has a (k^k/k!\) scale (and a stated ((2e)^k) general sup bound).
* **Proof mechanism.** Random restrictions and interpolation isolate a block; hypercontractivity and anti-concentration transfer tails.
* **Boundary.** This legitimizes comparison of same-spin and independent-block degree-two norms only with an explicit constant and normalization. It is not an equality and does not preserve an optimal complete-support constant.

## C. Real Grothendieck/Krivine toolkit (intrinsically bilinear)

### 9. Krivine (1979)

**Jean-Louis Krivine, “Constantes de Grothendieck et fonctions de type positif sur les sphères,” Adv. Math. 31 (1979), 16–30.** [DOI](https://doi.org/10.1016/0001-8708(79)90017-3)

* **Theorem/normalization.** For real matrices, the Hilbert/vector bilinear value is at most (K_G) times (\max_{\varepsilon_i,\delta_j=\pm1}|\sum a_{ij}\varepsilon_i\delta_j|), and
  \[
  K_G\le \frac{\pi}{2\log(1+\sqrt2)}=1.7822\ldots.
  \]
* **Proof mechanism.** The Gaussian sign identity (\mathbb E\operatorname{sgn}\langle g,u\rangle\operatorname{sgn}\langle g,v\rangle=(2/\pi)\arcsin\langle u,v\rangle); invert the correlation and realize its power series through tensor powers/nonlinear preprocessing, then random-hyperplane round.
* **Boundary.** Both sign families are independent. Symmetrizing (A) or setting (y=x) after the theorem is not licensed by Grothendieck’s inequality.

### 10. Braverman–Makarychev–Makarychev–Naor (2013)

**Mark Braverman, Konstantin Makarychev, Yury Makarychev, Assaf Naor, “The Grothendieck Constant Is Strictly Smaller than Krivine’s Bound,” Comm. Pure Appl. Math. 66 (2013).** [DOI](https://doi.org/10.1002/cpa.21398) · [arXiv](https://arxiv.org/abs/1103.6161)

* **Theorem.** (K_G<\pi/(2\log(1+\sqrt2))). The original proof establishes a positive but numerically tiny gap rather than a practically useful decimal improvement.
* **Proof mechanism.** A two-dimensional non-hyperplane Gaussian partition perturbs the hyperplane correlation; it is mixed with Krivine rounding so the inverse-series absolute-coefficient condition remains feasible. Hermite coefficients detect the improving direction.
* **Boundary.** This refutes Krivine’s conjectured value but still concerns independent-block bilinear rounding. The small perturbative gap is not a Boolean quadratic signing constant.

### 11. Naor–Regev (2014)

**Assaf Naor, Oded Regev, “Krivine Schemes Are Optimal,” Proc. AMS 142 (2014), 4315–4320.** [DOI](https://doi.org/10.1090/S0002-9939-2014-12169-1) · [arXiv](https://arxiv.org/abs/1205.6415)

* **Theorem/normalization.** For every projection dimension (k), there is a (k)-dimensional *oblivious mixed* Krivine scheme of quality ((1+O(1/k))K_G). Thus increasing-dimensional Gaussian partitions recover the true, still unknown bilinear constant.
* **Proof mechanism.** Duality first supplies a probability measure on pairs of sign partitions of (S^{k-1}). Its rotationally invariant Gaussian correlation (f_k(t)) is inverted. Rouché/Cauchy bounds control inverse coefficients, and nonlinear tensor-power maps (S,T:S^{\infty}\to S^{\infty}) realize the inverse correlation.
* **Boundary.** This is an existential dimension-transfer theorem for mixed schemes, not a finite explicit partition and not a same-spin theorem. It is the global input used by the 2026 lower-bound preprint below.

### 12. Heilman (2026 preprint)

**Steven Heilman, “An Upper Bound on Grothendieck’s Constant.”** [arXiv:2606.00247](https://arxiv.org/abs/2606.00247)

* **Theorems.** Explicit two-dimensional threshold schemes (f_\eta(x)=\operatorname{sgn}(x_2-\eta h_5(x_1))) prove the Braverman–et al. no-extra-randomization conjecture. The paper records analytic gaps (10^{-389}) (degree five) and (10^{-217}) (oppositely oriented degree three), and a **rigorous interval-arithmetic** certificate
  \[
  K_G<\frac\pi{2\log(1+\sqrt2)}-10^{-5}.
  \]
* **Proof mechanism.** Perturb the Gaussian correlation, analytically continue and invert it, and show its inverse absolute-coefficient sum crosses one beyond (\operatorname{arsinh}(1)); Taylor bounds give tiny analytic gaps, interval arithmetic gives the useful one.
* **Status/boundary.** Unrefereed 2026 preprint; distinguish the analytic theorem from the computer-assisted (10^{-5}) certificate. Numerically superseded by the August 2026 claim below if that claim survives review. Entirely bilinear.

### 13. Jones–Malavolta (2026 preprint)

**Chris Jones, Giulio Malavolta, “The Grothendieck Constant Is Strictly Larger than Davie–Reeds’ Bound.”** [arXiv:2603.30039](https://arxiv.org/abs/2603.30039) · [author PDF](https://chrisjones.space/assets/papers/grothendieck-strictly-larger.pdf)

* **Theorem/normalization.** With (K_{DR}) the classical Davie–Reeds Gaussian-operator lower bound, (K_G\ge K_{DR}+10^{-12}).
* **Proof mechanism.** Write the gap instance as a Hermite projection game (A=\sum c_k\Pi_k). A stability analysis classifies near optimizers of (A_{DR}=\Pi_1-\lambda_*I); all carry a uniform degree-three correlation (the paper obtains (\langle\Pi_3f,\Pi_3g\rangle\ge0.046) for the strip forms). Adding a small negative cubic projector improves the SDP/scalar ratio.
* **Status/boundary.** New unrefereed preprint; its numerical inequalities should be independently checked. Superseded as a numerical lower bound by the next source, but its perturb-near-extremizers mechanism remains useful. Bilinear (L_\infty(\gamma)\to L_1(\gamma)), not same-spin.

### 14. Saha–Li–Xue–Chaudhuri–Klivans–Kothari–Meka (August 2026 preprint)

**Rahul Saha, Alan Li, Anton Xue, Swarat Chaudhuri, Adam Klivans, Pravesh K. Kothari, Raghu Meka, “New Lower and Upper Bounds for the Grothendieck Constant,” arXiv v2, 12 Aug. 2026.** [arXiv:2608.11158v2](https://arxiv.org/abs/2608.11158v2) · [HTML](https://arxiv.org/html/2608.11158v2)

* **Main claim.**
  \[
  \frac{6\pi}{11}\le K_G\le\frac\pi{2\log(1+\sqrt2)}-3.47\cdot10^{-4}=1.7818\ldots.
  \]
  The abstract rounds the upper improvement to (10^{-4}).
* **Upper mechanism.** Normalize (H_{f,g}(t)=(\pi/2)\mathbb E[f(X)g(Y)]), write (H^{-1}(z)=\sum a_nz^n), and use (\sum|a_n|\gamma^n\le1\Rightarrow K_G\le\pi/(2\gamma)). “Limiting” schemes allow coordinate correlations (\rho(t)=\sum_{d\text{ odd}}c_dt^d), (\sum|c_d|\le1). Multivariate CLT plus Vitali and inverse-branch stability transfer them to honest high-dimensional schemes. The certified cubic–quintic choice has
  \[
  \rho(t)=\frac{t-s_3^2t^3+s_5^2t^5}{1+s_3^2+s_5^2},\quad
  (\eta,s_3,s_5)=(.136419125,.34101124,.05276111),
  \]
  and thresholds (\operatorname{sgn}(w\pm\vartheta\mathrm{He}_3(x))).
* **Lower mechanism.** For every odd sign pair, if (H=b_1t+b_3t^3+\cdots), prove the affine strip (b_3\ge2b_1-11/6). It survives mixing and coefficientwise limits. The first two inverse coefficients force every admissible (\gamma\le11/12); Naor–Regev optimality then gives (K_G\ge6\pi/11). The strip proof uses (h=(f+g)/2,k=(f-g)/2), cancellation of the (P_1k) term, one-dimensional rearrangement, and a certified ternary moment/fiber inequality.
* **Rigor flag.** This was posted three days before retrieval, is unrefereed, was AI-assisted, and both halves contain reproducible interval-arithmetic certificates (the authors state they independently verified every theorem). Treat it as a major **provisional claim**, not settled folklore or a machine-only fact. It remains a bilinear result.

## D. Tensor powers, completely bounded norms, and vector-valued factorization

### 15. Aubrun–Müller-Hermes (2025)

**Guillaume Aubrun, Alexander Müller-Hermes, “Limit Formulas for Norms of Tensor Power Operators,” J. Funct. Anal. 289 (2025), article 111113.** [DOI](https://doi.org/10.1016/j.jfa.2025.111113) · [arXiv](https://arxiv.org/abs/2410.23063)

* **Theorem/normalization.** For any bounded (\phi:X\to Y), with the convention (+\infty) for unbounded tensor maps,
  \[
  \lim_{k\to\infty}\|\phi^{\otimes k}:X^{\otimes_\varepsilon k}\to Y^{\otimes_\pi k}\|^{1/k}=\gamma_2^*(\phi),
  \]
  the 2-dominated norm. It is finite exactly for 2-dominated (\phi). For (\phi:X\to H), the (\varepsilon\to h) limit is the 2-summing norm (\pi_2(\phi)). Corollary 4.6 gives exact multiplicativity (\gamma_2^*(\phi^{\otimes k})=\gamma_2^*(\phi)^k).
* **Proof mechanism.** Pietsch/Hilbert factorization gives the limsup. Trace duality, Haar averaging, and random nearly entangled Hilbert tensors give the reverse inequality; finite-dimensional compressions extend it to general Banach spaces.
* **Boundary.** This is the correct way to prove existence of a regularized *full tensor-power* norm. Symmetric powers or diagonal evaluations are restrictions/quotients and may incur polarization loss, so this theorem does not imply existence or multiplicativity for (M_n/n^{3/2}).

### 16. Arunachalam–Dutt–Escudero Gutiérrez–Palazuelos (2025)

**Srinivasan Arunachalam, Arkopal Dutt, Francisco Escudero Gutiérrez, Carlos Palazuelos, “A cb-Bohnenblust–Hille Inequality with Constant One and Its Applications in Learning Theory,” Math. Ann. 392 (2025), 3367–3396.** [DOI](https://doi.org/10.1007/s00208-025-03142-5)

* **Theorem/normalization.** For complex (d)-homogeneous (P=\sum_{|\alpha|=d}a_\alpha z^\alpha),
  \[
  \|(a_\alpha)\|_{2d/(d+1)}\le\|P\|_{cb},
  \quad
  \|P\|_{cb}=\sup_{\mathcal H,\,\|Z_i\|\le1}
  \Big\|\sum a_\alpha Z_1^{\alpha_1}\cdots Z_n^{\alpha_n}\Big\|,
  \]
  and both constant one and exponent are optimal. Proposition 4.2 separately gives the exact (2^{(d-1)/d}) coefficient bound for **Boolean-valued outputs** (f:\{\pm1\}^n\to\{\pm1\}).
* **Proof mechanism.** Encode polynomial coefficients in a one-sided tensor, test it on explicit shift/contraction matrices, and combine mixed-norm Hölder with noncommutative Cauchy–Schwarz.
* **Boundary.** (\|P\|_{cb}\ge\|P\|_\infty), often strictly; constant one in the larger norm does not sharpen the scalar same-spin problem. “Boolean-valued output” is also different from a real-valued polynomial with Boolean input and unimodular coefficients.

### 17. Briët–Escudero Gutiérrez–Gribling (2024)

**Jop Briët, Francisco Escudero Gutiérrez, Sander Gribling, “Grothendieck Inequalities Characterize Converses to the Polynomial Method,” Quantum 8 (2024), 1526.** [DOI](https://doi.org/10.22331/q-2024-11-18-1526) · [arXiv](https://arxiv.org/abs/2212.08559)

* **Theorems/normalization.** For bounded bilinear forms, the worst one-query additive approximation error is exactly (1-1/K_G). In degree four, no uniform additive error (<1) exists for the analogous two-query converse; the obstruction is an unbounded separation between scalar and completely bounded norms.
* **Proof mechanism.** Tensor-norm duality identifies quantum-query forms with cb tensors; convex separation turns a universal approximation statement into the Grothendieck norm ratio.
* **Boundary.** The positive statement is bilinear with independent inputs; the negative quartic result warns against extrapolating a degree-two factorization. Neither theorem identifies the diagonal restriction (y=x).

### 18. Defant–Junge (1999)

**Andreas Defant, Marius Junge, “A Vector-Valued Grothendieck Inequality with an Application to ((p,q))-Completely Bounded Operators,” Indiana Univ. Math. J. 48 (1999), 295–310.** [DOI](https://doi.org/10.1512/iumj.1999.48.1692)

* **Theorem/normalization.** Their vector-valued Grothendieck inequality yields equivalence, up to a universal constant, between the operator-space cb norm and the a priori larger ((\infty,1))-completely bounded norm. The bridge is a vector-valued Maurey–Rosenthal factorization theorem.
* **Proof mechanism.** Turn a vector-valued family inequality into a weighted/factorized estimate through suitable (L_p)-type spaces, then apply scalar Grothendieck at the factor level.
* **Boundary.** This is the nearest rigorous source I found to the prompt’s “nonlinear factorization / jointly-canceling vector-valued” language; **“jointly-canceling” is not the paper’s terminology**, and no primary theorem with that exact descriptor was located. Operator-space factorization does not preserve complete support, symmetry, or a single shared sign vector.

## Method map for the project

1. **Getting the (n^{3/2}) scale.** Boolean BH gives the universal lower scale because (\binom n2) coefficients of modulus one have (\ell_{4/3})-norm (\binom n2^{3/4}). Boolean KSZ gives a full-edge signing with (O(\sqrt{n\binom n2})) norm. This is the clean direct sandwich; its constants are not matched.

2. **Dimension transfer.** Two reliable patterns recur. On the cube, rewrite a symmetric Walsh kernel as a polynomial of (n^{-1/2}\sum x_i), apply CLT, then prove uniform integrability before passing an (L_1) norm. For Gaussian rounding schemes, use Hermite blocks and multivariate CLT for pointwise correlation convergence, then Vitali for local uniform holomorphic convergence and Rouché/Cauchy estimates for inverse-series stability. Neither transfer automatically preserves a minimum over signings.

3. **Extracting sharp constants/extrema.** In fixed dimension, convexity reduces a polynomial norm problem to extreme points and then a finite/one-variable optimization (Jiménez-Rodríguez et al.). For dimension-free coefficient inequalities, mixed-norm interpolation and Khinchine constants identify the exponent and sometimes the constant. Complete-support minimization is a much smaller, nonconvex slice of the coefficient ball; generic extremizers may have sparse or diagonal support.

4. **Proving an asymptotic norm exists.** Full tensor powers are submultiplicative and, more strongly, their regularization has the single-letter value (\gamma_2^*). A graph-size or symmetric-diagonal sequence needs its own composition law with controlled interface error. Fekete cannot be invoked until such a law is proved; symmetrization and deletion of diagonal terms are not cost-free.

5. **Symmetry/polarization loss.** (q(x)=B(x,x)) forgets the off-diagonal values of (B(x,y)). Polarization reconstructs them by evaluating (q(x\pm y)), which expands the domain and costs up to (2) at degree two in a general real Banach space. Decoupling gives distributional/supremum comparisons with explicit constants but not equality. Therefore a bilinear Grothendieck or sharp Littlewood constant is, at best, a surrogate bound until a same-spin transfer theorem is supplied.

6. **Computer-assisted claims.** Heilman’s (10^{-5}) upper improvement is explicitly interval-certified. The August 2026 Saha et al. upper and part of the lower proof use reproducible interval certificates and the paper reports extensive AI assistance; the authors claim human verification, but the result is days old and unrefereed. Jones–Malavolta is also a new preprint. Keep all three labeled provisional when quoting the current numerical frontier.

## Bottom-line relevance boundaries

* **Directly on target:** the Boolean BH lower bound and Boolean KSZ/Sidon upper construction, both at order (n^{3/2}).
* **Useful but lossy:** polarization/decoupling from (q(x)) to (B(x,y)); any constant obtained after that passage includes the transfer loss.
* **Different optimization:** Grothendieck compares scalar and Hilbert **bilinear** values for a fixed matrix; it does not minimize over symmetric complete-support sign matrices.
* **Different field/domain:** Steinhaus/polydisc BH and cb-BH constants cannot be substituted for real cube constants.
* **No retrieved theorem:** existence or value of (\lim M_n/n^{3/2}), a sharp complete-support same-spin constant, or an exact “jointly-canceling vector-valued” inequality tailored to this class.
