# Spin-glass thermodynamic-limit toolkit for dense quadratic Ising energies

Retrieval date: 2026-08-15. This is a literature packet, not a proposed solution. It records theorem-level inputs and, especially, the point at which each input stops applying to deterministic optimization over the disorder.

## 0. Normalization dictionary

For a symmetric zero-diagonal coupling array \(J=(J_{ij})_{i<j}\), put

\[
Q_J(\sigma)=\sum_{i<j}J_{ij}\sigma_i\sigma_j,
\qquad \sigma\in\{-1,1\}^N.
\]

The usual SK normalization is \(H_N^J(\sigma)=N^{-1/2}Q_J(\sigma)\). Therefore

\[
\frac1N\max_\sigma H_N^J(\sigma)
=\frac{1}{N^{3/2}}\max_\sigma Q_J(\sigma).
\]

With i.i.d. variance-one disorder and the \(i<j\) convention,

\[
\mathbb E H_N(\sigma)H_N(\tau)
=\frac N2R(\sigma,\tau)^2+O(1),
\qquad R(\sigma,\tau)=N^{-1}\sum_i\sigma_i\tau_i.
\]

Thus the corresponding mixed-spin covariance function is \(\xi(r)=r^2/2\). Papers summing over all ordered pairs often use \(\xi(r)=r^2\); this changes inverse-temperature and ground-state constants by a factor of \(\sqrt2\). Papers integrating against the uniform probability measure \(2^{-N}\sum_\sigma\) subtract \(\log 2\) from pressures defined with counting measure.

For counting measure, if

\[
p_N^J(\beta)=\frac1N\log\sum_\sigma e^{\beta H_N^J(\sigma)},
\qquad G_N(J)=\frac1N\max_\sigma H_N^J(\sigma),
\]

then the exact log-sum-exp bounds are

\[
G_N(J)\le \frac{p_N^J(\beta)}\beta
\le G_N(J)+\frac{\log2}{\beta}.
\]

Hence a finite-temperature limit uniform enough in \(N\), followed by \(\beta\to\infty\), is the standard route to a ground-state limit.

The absolute-energy version is not a new one-sided Hamiltonian:

\[
A_N(J):=\max_\sigma|Q_J(\sigma)|
=\max\left\{\max_\sigma Q_J(\sigma),\max_\sigma Q_{-J}(\sigma)\right\}.
\]

For a sign-symmetric random law, the two one-sided maxima have the same marginal distribution, but they are dependent. More importantly, symmetry of a random law does not commute with an outer deterministic minimization over \(J\).

## 1. Executive synthesis

### What are the native objects?

The mature theory natively controls one of four objects:

1. quenched pressure \(N^{-1}\mathbb E_J\log Z_N(J)\);
2. a typical-sample pressure \(N^{-1}\log Z_N(J)\), after concentration around the quenched mean;
3. a typical ground state \(N^{-1}\max_\sigma H_N^J(\sigma)\), usually obtained from the pressure at zero temperature;
4. constrained pressures, where the constraint is on magnetization, self-overlap, or overlaps among replicas.

It does **not** natively control

\[
\inf_{J\in\{\pm1\}^{\binom N2}}\max_\sigma |Q_J(\sigma)|,
\]

because there the disorder is a decision variable rather than a random environment.

### Which limits are standard?

- For Gaussian SK and even mixed \(p\)-spin models, the quenched pressure, almost-sure pressure, quenched ground state, and almost-sure ground state all have thermodynamic limits. Guerra--Toninelli size interpolation is the clean existence proof.
- For independent mean-zero variance-one non-Gaussian couplings, Lindeberg replacement transfers the Gaussian pressure and one-sided ground-state limit under a finite-variance Lindeberg condition. This includes Rademacher disorder.
- For positive-semidefinite multi-species covariance, a Parisi formula and hence a limit are known. As of June 2026, convexity is no longer required for **centered Ising** multi-species models, but centeredness is essential to the new argument.
- Constrained free energies can still converge when exact constraints live on incompatible \(1/N\)-lattices. The standard repair is a near-superadditive inequality with a sublinear defect plus a restricted Fekete/de Bruijn--Erdős lemma.
- Heavy-tailed models require a different normalization and a sparse combinatorial interpolation; finite-variance universality is not a cosmetic technicality.

### Which mechanisms transfer values across system sizes?

1. **Guerra--Toninelli split interpolation:** compare one \(N=N_1+N_2\) Gaussian process with independent processes on the two blocks; covariance convexity signs the derivative and gives exact or approximate superadditivity.
2. **Cell decomposition:** first pin additive one-replica/self-overlap parameters to small cells, then interpolate within each cell. This prevents diagonal order parameters from drifting during the split.
3. **Near-superadditivity with a defect:** allow \(a_{N_1+N_2}\ge a_{N_1}+a_{N_2}-o(N)\), usually only for comparable block sizes; a generalized Fekete lemma still gives a limit.
4. **Cavity/ASS increments:** represent the limit through \(M\)-spin increments \(Q_{N+M}-Q_N\), with \(M\) fixed before \(N\to\infty\); the reservoir becomes a random overlap structure.
5. **Sparse edge interpolation:** add/delete or relocate random edges one at a time, used when heavy tails make the dense Gaussian covariance interpolation unavailable.

Every one of these uses a distributional consistency or covariance structure across sizes. An outer optimizer over deterministic \(J\) destroys that input unless an additional deterministic extension/restriction lemma is supplied.

### Quenched, annealed, constrained, and adversarial are different quantifiers

| Statement | What it gives | What it does not give |
|---|---|---|
| \(N^{-1}\mathbb E\log Z_N\to p\) | average random-sample value | a bound for every coupling matrix |
| \(N^{-1}\log Z_N\to p\) in probability/a.s. | typical-sample value; in particular, witnesses exist | the optimum over all samples or a limit for that optimum |
| universality under iid replacement | same typical value for two random ensembles | robustness under choosing an exceptional array after seeing all entries |
| constrained overlap/magnetization formula | control of selected spin/replica sectors | a constraint on the disorder array |
| \(\sup_m\inf_\pi\) or \(\sup_p\inf_q\) formula | saddle point over order parameters | an adversary selecting \(J\) |
| a disorder large-deviation principle | cost of rare random samples | a deterministic uniform theorem unless its speed and constant beat the cardinality of the disorder class |

For Bernoulli disorder there are \(2^M\) arrays, \(M=\binom N2\), and every array has mass \(2^{-M}\). Thus a random bound \(\mathbb P(E_N)\le e^{-cN^2}\) rules out **all** Bernoulli arrays in \(E_N\) only if it is sharp enough to be \(<2^{-M}\), asymptotically requiring \(c>(\log2)/2\) in this normalization. A merely positive quadratic rate is not enough. Gaussian large deviations do not transfer at this exponential resolution through ordinary Lindeberg universality.

### Closest currently available results to optimizing the disorder

- Huang--Sellke prove that suppressing the spherical Gaussian ground state below a threshold costs at least \(e^{-cN^2}\). This is the correct *direction* for an outer minimization and the correct quadratic speed, but it is spherical, Gaussian, and has no exact rate constant.
- Chen--Guionnet--Ko--Lacroix-A-Chez-Toine--Mourrat compute the \(e^{-N I(r)}\) **upper** tail of the Ising Gaussian ground state. It is an exact modern disorder LDP, but its tail direction is opposite to disorder minimization.
- Recent spherical universality theorems sharply identify when rare large couplings dominate. They show why typical universality cannot be presumed under either heavy tails or deliberate disorder selection.
- No primary theorem in this packet is uniform over all deterministic \(\{\pm1\}\) coupling arrays, and none proves a thermodynamic limit after an outer minimization over those arrays.

## 2. Primary theorem toolkit

### 1. Guerra--Toninelli: exact system-size interpolation

**Source.** Francesco Guerra and Fabio L. Toninelli, [*The Thermodynamic Limit in Mean Field Spin Glass Models*](https://arxiv.org/abs/cond-mat/0204280), CMP 230 (2002), 71--79.

**Hypotheses and normalization.** Ising spins; independent standard Gaussian \(J_{ij}\); SK Hamiltonian

\[
H_N(\sigma,h,J)=-N^{-1/2}\sum_{i<j}J_{ij}\sigma_i\sigma_j-h\sum_i\sigma_i,
\quad
\alpha_N=N^{-1}\mathbb E\log Z_N.
\]

**Conclusion.** For every split \(N=N_1+N_2\),

\[
N\alpha_N\ge N_1\alpha_{N_1}+N_2\alpha_{N_2}.
\]

Consequently \(\alpha_N\to\sup_N\alpha_N\). The analogous averaged ground-state energies are superadditive and converge. Gaussian concentration upgrades both pressure and ground-state convergence to almost sure convergence. Moreover the zero-temperature limits commute through the exact sandwich

\[
0\le \frac{1}{\beta N}\log Z_N-e_N(J,h)\le \frac{\log2}{\beta},
\]

with the paper's sign convention for \(e_N\).

**Mechanism.** Interpolate between the full \(N\)-spin field and independent fields on the two blocks. Gaussian integration by parts makes the derivative a multiple of

\[
-\Big\langle R_{12}^2-\frac{N_1}{N}(R_{12}^{(1)})^2-\frac{N_2}{N}(R_{12}^{(2)})^2\Big\rangle.
\]

Since \(R_{12}\) is the weighted average of its block overlaps, convexity of \(x^2\) gives the sign.

**Extensions in the paper.** Even \(p\)-spin models use convexity of \(x^p\). Symmetric non-Gaussian disorder with suitable moment bounds gives an approximate interpolation with a sublinear error; bounded disorder also has concentration.

**Boundary.** The theorem interpolates *quenched averages of compatible random ensembles*. It neither selects one deterministic optimizer at each size nor couples such optimizers coherently across sizes.

### 2. Guerra--Toninelli: generalized mean-field template and cell decomposition

**Source.** Francesco Guerra and Fabio L. Toninelli, [*The Infinite Volume Limit in Generalized Mean Field Disordered Models*](https://arxiv.org/abs/cond-mat/0208579), MPRF 9 (2003), 195--207.

**Hypotheses.** Product single-spin spaces with bounded additive order parameters. The Hamiltonian mean and covariance have, uniformly in configurations, the forms

\[
N^{-1}b_N(\sigma)=g(m_N(\sigma))+O(N^{-1}),
\qquad
N^{-1}c_N(\sigma,\tau)=f(Q_N(\sigma,\tau))+O(N^{-1}),
\]

where \(m_N,Q_N\) decompose as weighted averages under a block split, \(g\) is differentiable, and \(f\) is differentiable and convex on the relevant overlap domain. The diagonal/self-overlap coordinates are bounded.

**Conclusion.** The quenched free-energy density and ground-state energy density have infinite-volume limits for this class, covering generalized SK, even \(p\)-spin, non-Ising, and coupled-replica models.

**Mechanism.** Partition configuration space into small cells that pin the additive magnetization and diagonal overlap. Interpolate separately on each cell, use convexity for the off-diagonal covariance, control the cell count subexponentially, and then shrink the cells.

**Boundary.** Convexity of the covariance is structural here. The theorem allows deterministic mean interactions only through fixed additive order parameters; it does not allow the full \(O(N^2)\)-dimensional coupling array to be optimized.

### 3. Aizenman--Sims--Starr: cavity increments and random overlap structures

**Source.** Michael Aizenman, Robert Sims, and Shannon Starr, [*An Extended Variational Principle for the SK Spin-Glass Model*](https://arxiv.org/abs/cond-mat/0306386) (2003).

**Hypotheses and normalization.** Gaussian mean-field fields with covariance \(Nf(R)/2\), where \(f\) is convex (the paper treats SK and even mixed \(p\)-spin classes). A random overlap structure (ROSt) consists of random summable weights \((\xi_\alpha)\) and a positive-semidefinite overlap kernel \(q_{\alpha\alpha'}\) with diagonal one. Gaussian cavity fields have covariances determined by \(f'(q)\) and \(qf'(q)-f(q)\).

**Conclusion.** For the paper's \(M\)-spin ROSt trial functional \(G_M\),

\[
P_M\le \inf_{\mathrm{ROSt}}G_M\le P_U,
\qquad
P=\lim_{M\to\infty}\inf_{\mathrm{ROSt}}G_M,
\]

where \(P_U\) is the limiting upper pressure and \(P\) is the thermodynamic pressure once existence is invoked.

**Mechanism.** The upper comparison is Guerra interpolation. The reverse direction takes a large \(N\)-spin reservoir and adds \(M\) cavity spins. A superadditive-increment lemma converts \(Q_N/N\) into fixed-\(M\) increments \((Q_{N+M}-Q_N)/M\); the actual cavity ratio differs from the ideal ROSt functional by \(O(M/N)\).

**Boundary.** The order of limits is \(N\to\infty\) with \(M\) fixed, then \(M\to\infty\). The reservoir is itself a quenched random Gibbs object. This is not a deterministic recurrence for an optimizer over disorder.

### 4. Panchenko: exact overlap constraints with a sublinear size defect

**Source.** Dmitry Panchenko, [*A Note on the Free Energy of the Coupled System in the Sherrington--Kirkpatrick Model*](https://arxiv.org/abs/math/0405359) (2004).

**Hypotheses.** Two jointly Gaussian mixed \(p\)-spin Hamiltonians \(H_N^1,H_N^2\), with coefficient sequences \(a_p^\ell\), external fields \(h_1,h_2\), and covariance functions \(\xi_{\ell\ell'}\) convex on \([-1,1]\). Spins are constrained by the exact cross-overlap \(R(\sigma^1,\sigma^2)=u_N\), where \(u_N\in\{-1,-1+2/N,\ldots,1\}\) and \(u_N\to u\).

**Conclusion.** The constrained pressure

\[
F_N(u_N)=\frac1N\mathbb E\log
\sum_{R(\sigma^1,\sigma^2)=u_N}
e^{H_N^1(\sigma^1)+H_N^2(\sigma^2)+h_1\sum\sigma_i^1+h_2\sum\sigma_i^2}
\]

has a limit depending only on \(u\), not on the approximating lattice sequence. If \(U_{N,\varepsilon}\) is an \(\varepsilon\)-window around the target overlap, the paper proves an exact-to-window comparison of the form

\[
F_N(U_{N,\varepsilon})\le F_N(u_N)+L\sqrt\varepsilon.
\]

Its size interpolation yields superadditivity after subtracting \(A\sqrt N\), for comparable block sizes; a restricted de Bruijn--Erdős/Fekete lemma gives convergence. An ASS-type variational description for the coupled system also follows.

**Mechanism.** Approximate incompatible target overlaps in the two blocks, quantify the rounding error, and absorb it into a \(\sqrt N\) defect.

**Boundary.** This is an important model for how to handle lattice-incompatible constraints, but the constraint is on two replicas, not on the disorder. Gaussian independence and convex overlap covariance remain essential.

### 5. Carmona--Hu: direct universality for pressure and ground state

**Source.** Philippe Carmona and Yueyun Hu, [*Universality in Sherrington--Kirkpatrick's Spin Glass Model*](https://arxiv.org/abs/math/0403359), AIHP 42 (2006), 215--222.

**Hypotheses.** I.i.d. disorder \(\xi\) with \(\mathbb E\xi=0\), \(\mathbb E\xi^2=1\), and \(\mathbb E|\xi|^3<\infty\). The paper uses an all-\((i,j)\) summation convention, so constants should be converted before comparison with \(i<j\).

**Conclusion.** The normalized quenched pressure converges almost surely and in mean to the Gaussian SK limit. Quantitatively, in the paper's convention,

\[
|\alpha_N(\beta,\xi)-\alpha_N(\beta,g)|
\le 9\,\mathbb E|\xi|^3\,\frac{\beta^3}{\sqrt N}.
\]

For \(S_N(\xi)=\sup_\sigma\sum_{i,j}\xi_{ij}\sigma_i\sigma_j\), \(N^{-3/2}S_N(\xi)\) converges almost surely and in mean to the same Gaussian ground-state constant.

**Mechanism.** Approximate Gaussian integration by parts/one-coordinate replacement, followed by the finite-temperature-to-zero-temperature sandwich.

**Boundary.** This compares typical iid ensembles. It has no uniformity over all deterministic arrays and does not resolve a rare lower tail selected by an outer minimization.

### 6. Chatterjee: Lindeberg invariance at the native \(N^{3/2}\) maximum scale

**Source.** Sourav Chatterjee, [*A Simple Invariance Theorem*](https://arxiv.org/abs/math/0508213) (2005 preprint).

**Hypotheses.** Independent triangular arrays \(J_{ij}\) with mean zero, variance one, and the Lindeberg condition

\[
\forall\varepsilon>0,\qquad
N^{-2}\sum_{i<j}\mathbb E[J_{ij}^2;|J_{ij}|>\varepsilon\sqrt N]\longrightarrow0.
\]

I.i.d. mean-zero variance-one entries satisfy it with no third-moment assumption.

**Conclusion.** The SK pressure has the same asymptotics as the Gaussian model. At zero temperature,

\[
N^{-3/2}\max_{\sigma\in\{\pm1\}^N}\sum_{i<j}J_{ij}\sigma_i\sigma_j
\]

converges in probability and in mean to the Gaussian constant. With uniformly bounded third absolute moments, the pressure comparison error is \(O(N^{-1/2})\); for smooth test functions of the maximum, the smoothing argument gives \(O(N^{-1/6})\).

**Mechanism.** Replace entries one at a time and Taylor-expand through third order. For a log-sum-exp \(F=\alpha^{-1}\log\sum_f e^{\alpha f}\), the paper proves

\[
\lambda_2(F)\le3\alpha\lambda_2(\mathcal F),
\qquad
\lambda_3(F)\le13\alpha^2\lambda_3(\mathcal F).
\]

The hard maximum is approximated uniformly within \(\alpha^{-1}\log|\mathcal F|\).

**Boundary.** The result is distributional/averaged. Lindeberg replacement controls a smooth function under two product laws; it does not say that this function is close on every pair of deterministic input arrays. It is therefore unsafe at the outer-disorder-optimization quantifier.

### 7. Panchenko: positive-semidefinite multi-species Parisi formula

**Source.** Dmitry Panchenko, [*The Free Energy in a Multi-Species Sherrington--Kirkpatrick Model*](https://arxiv.org/abs/1310.6679), Ann. Probab. 43 (2015), 3494--3513.

**Hypotheses and normalization.** A finite species set \(\mathscr S\); partitions \(I_s\) with \(|I_s|/N\to\lambda_s\in(0,1)\); Ising spins; independent Gaussian \(g_{ij}\) with variance \(\Delta_{st}^2\) for \(i\in I_s,j\in I_t\); symmetric positive-semidefinite matrix \(\Delta^2\). The paper uses

\[
H_N(\sigma)=N^{-1/2}\sum_{i,j=1}^N g_{ij}\sigma_i\sigma_j,
\qquad
R_s(\sigma,\tau)=|I_s|^{-1}\sum_{i\in I_s}\sigma_i\tau_i.
\]

**Conclusion.** \(N^{-1}\mathbb E\log\sum_\sigma e^{H_N(\sigma)}\) converges to the infimum of the multi-species Parisi functional over a common cascade distribution \(\zeta\) and monotone species paths \(q^s\). The functional uses

\[
Q_\ell=\sum_{s,t}\Delta_{st}^2\lambda_s\lambda_tq_\ell^sq_\ell^t,
\qquad
Q_\ell^s=2\sum_t\Delta_{st}^2\lambda_tq_\ell^t.
\]

**Mechanism.** Guerra interpolation gives an error

\[
-\tfrac12\mathbb E\langle(\Delta^2(R-q),R-q)\rangle\le0.
\]

The ASS lower bound is closed using multi-species Ghirlanda--Guerra identities and synchronization: each species overlap becomes a deterministic monotone function of the total overlap.

**Boundary.** Positive semidefiniteness is exactly what signs the Guerra remainder. The paper is quenched Gaussian. The 2026 centered-Ising theorem below removes convexity only by a different mechanism and still does not optimize disorder.

### 8. Auffinger--Chen: the zero-temperature Parisi formula

**Source.** Antonio Auffinger and Wei-Kuo Chen, [*Parisi Formula for the Ground State Energy in the Mixed \(p\)-Spin Model*](https://arxiv.org/abs/1606.05335), Ann. Probab. 45 (2017), 4617--4631.

**Hypotheses.** Ising spins and a centered Gaussian mixed \(p\)-spin field

\[
H_N(\sigma)=\sum_{p\ge2}\frac{c_p}{N^{(p-1)/2}}
\sum_{i_1,\ldots,i_p}g_{i_1\cdots i_p}\sigma_{i_1}\cdots\sigma_{i_p}
+h\sum_i\sigma_i,
\]

with \(\sum_p2^pc_p^2<\infty\) and at least one nonzero coefficient; \(\xi(s)=\sum_pc_p^2s^p\). For the \(i<j\) SK normalization, use \(\xi(s)=s^2/2\).

**Conclusion.** If \(L_N=N^{-1}\max_\sigma H_N(\sigma)\), then \(L_N\to\mathrm{GSE}\) almost surely and

\[
\mathrm{GSE}=\inf_{\gamma\in\mathcal U}
\left\{\Psi_\gamma(0,h)-\frac12\int_0^1t\xi''(t)\gamma(t)\,dt\right\}.
\]

Here \(\mathcal U\) is the set of nonnegative, nondecreasing, right-continuous, integrable functions on \([0,1)\), and

\[
\partial_t\Psi_\gamma
=-\frac{\xi''(t)}2\left(\partial_{xx}\Psi_\gamma
+\gamma(t)(\partial_x\Psi_\gamma)^2\right),
\qquad \Psi_\gamma(1,x)=|x|.
\]

**Mechanism.** Take the zero-temperature limit of the finite-temperature Parisi formula. A stochastic-control representation gives compactness of the rescaled Parisi measures and controls the apparent singularity at \(t=1\); the PDE and linear terms cancel the endpoint mass correctly.

**Boundary.** This precisely identifies the typical Gaussian one-sided maximum. It is not a formula for an outer minimization over \(J\), and it does not by itself handle the joint pair \(J,-J\) in an absolute maximum.

### 9. Barbier--Macris: adaptive interpolation and its Nishimori boundary

**Source.** Jean Barbier and Nicolas Macris, [*The Adaptive Interpolation Method: A Simple Scheme to Prove Replica Formulas in Bayesian Inference*](https://arxiv.org/abs/1705.02780), PTRF 174 (2019), 1133--1185.

**Hypotheses.** Bayes-optimal rank-one symmetric matrix estimation: \(S_i\stackrel{\mathrm{iid}}\sim P_0\) with bounded support and observations

\[
W_{ij}=S_iS_j/\sqrt n+\sqrt\Delta\,Z_{ij},
\]

with independent Gaussian noise and the posterior using the true prior/channel. Let \(f_n=-n^{-1}\mathbb E\log Z\).

**Conclusion.** The free energy converges to the replica-symmetric scalar variational formula

\[
\lim_{n\to\infty}f_n=\min_{m\ge0}f_{\mathrm{RS}}(m;\Delta),
\qquad
f_{\mathrm{RS}}(m;\Delta)=\frac{m^2}{4\Delta}
+f_{\mathrm{den}}\!\left(\sqrt{\Delta/m}\right).
\]

**Mechanism.** A \(K\)-step interpolation replaces \(1/K\) of the dense matrix channel at a time by scalar Gaussian channels. A vanishing side channel enforces overlap concentration. Nishimori identities turn signal-overlap terms into replica overlaps. The trial \(m_k\) is chosen adaptively from the expected overlap so the sum-rule remainder has the desired sign/cancellation; then \(n\to\infty\), \(K\to\infty\), and the perturbation vanishes.

**Boundary.** This method's decisive identities come from the matched Bayesian/Nishimori setting. The paper explicitly contrasts this with low-temperature SK, where the needed scalar overlap concentration is unavailable in general. It is not an interpolation for generic or adversarial disorder.

### 10. Jagannath--Lopatto: thermodynamic limits beyond finite variance

**Source.** Aukosh Jagannath and Patrick Lopatto, [*Existence of the Free Energy for Heavy-Tailed Spin Glasses*](https://arxiv.org/abs/2211.09879), CMP 405 (2024), article 231.

**Hypotheses.** I.i.d. symmetric couplings with exact power tail

\[
\mathbb P(|J|\ge t)=C_0t^{-\alpha}\quad(t>1),
\qquad \mathbb E|J|<\infty,
\qquad 1<\alpha<2,
\]

and

\[
H_N(\sigma)=N^{-1/\alpha}\sum_{i<j}J_{ij}\sigma_i\sigma_j.
\]

**Conclusion.** For every \(\beta>0\), \(N^{-1}\mathbb E\log Z_N(\beta)\) has a finite limit. The free energy is self-averaging; specifically, for each \(\delta>0\),

\[
\mathbb P\left(N^{-1}|\log Z_N-\mathbb E\log Z_N|>t\right)
\le C N^{1-\alpha+\delta}/t^2.
\]

**Mechanism.** Truncate small couplings at \(N^{1/\alpha-\varepsilon}\), reduce the model to a sparse weighted multigraph, fix its edge count, and apply Bayati--Gamarnik--Tetali combinatorial edge interpolation. The resulting almost-superadditivity has an integrable sublinear defect. A martingale argument yields self-averaging.

**Boundary.** The normalization changes from \(N^{-1/2}\) to \(N^{-1/\alpha}\), and a few extreme edges govern the geometry. This is evidence that moment hypotheses encode a real localization transition, not a disposable proof detail.

### 11. Chen: a genuine max--min formula, but only over order parameters

**Source.** Hong-Bin Chen, [*Free Energy in Spin Glass Models with Conventional Order*](https://arxiv.org/abs/2401.10223) (2024; revised 2025).

**Hypotheses.** Vector spins \(\sigma_i\in\mathbb R^D\) drawn iid from a finite measure supported on the unit ball; centered Gaussian field with

\[
\mathbb E H_N(\sigma)H_N(\tau)
=N\xi(N^{-1}\sigma\tau^\top),
\]

where \(\xi\) satisfies the paper's regularity assumptions and is convex on the positive-semidefinite cone. A conventional order is \(m_N=N^{-1}\sum_i h(\sigma_i)\), with bounded measurable \(h\), and \(G\) is locally Lipschitz, not necessarily convex.

**Conclusion.** For the self-overlap-corrected model,

\[
F_N^{\mathrm{soc},G}=\frac1N\mathbb E\log\int
e^{H_N(\sigma)-\frac N2\xi(N^{-1}\sigma\sigma^\top)+NG(m_N)}\,dP_N,
\]

the limit is

\[
\sup_{m\in\mathbb R^d}\inf_{\pi\in\Pi}\inf_{x\in\mathbb R^d}
\left\{\mathscr P^h(\pi,x)-m\!\cdot x+G(m)\right\}.
\]

The paper also removes the self-overlap correction by adding the self-overlap matrix as a conventional order parameter.

**Mechanism.** Interpolate to the corrected, overlap-only model along a Hamilton--Jacobi PDE and apply a Hopf formula; an alternate proof conditions on small self-overlap cells.

**Boundary.** This is an important warning about terminology: its “max--min” is over a finite-dimensional conventional order and a Parisi overlap path. It is not a min--max game between spins and disorder, and convex covariance is still assumed.

### 12. Sawhney--Sellke: sharp light-tail universality for spherical models

**Source.** Mehtaab Sawhney and Mark Sellke, [*Free Energy Universality of Spherical Spin Glasses*](https://arxiv.org/abs/2408.13701) (2024 preprint).

**Hypotheses.** A finite mixed spherical \(p\)-spin model on \(\|\sigma\|_2=\sqrt N\), with independent symmetric-tensor entries normalized so each \(p\)-layer has the Gaussian covariance \(N\gamma_p^2R^p\). Under a uniform \((C,\varepsilon)\) condition, entries in the \(p\)-layer have bounded \((2p+\varepsilon)\)-moments. For the asymptotic iid theorem, mean zero, variance one, and finite \(2p\)-th moment are assumed in each active layer.

**Conclusion.** Uniformly over \(\beta\in[0,\infty]\), including the ground state,

\[
\mathbb E\big|F_{N,\beta}(J)-\mathbb EF_{N,\beta}(G)\big|
\le N^{1-c}
\]

for some \(c=c(P,\varepsilon)>0\) under the stronger uniform moments. For fixed iid laws with finite \(2p\)-moments, \(F_{N,\beta}/N\) converges almost surely to the Gaussian spherical Parisi value for every \(\beta\le\infty\). If an active pure \(p\)-layer has infinite \(2p\)-moment, the normalized ground state has infinite limsup.

**Mechanism.** Multiscale truncation separates huge entries; the sphere is covered by a subexponential family of delocalized subspheres; Lindeberg replacement works on each; rotational invariance controls induced external fields.

**Boundary.** This is spherical rather than Ising and remains a theorem for sampled independent disorder. Its sharp \(2p\) threshold is a concrete failure mode: localized coordinate spikes can dominate the ground state.

### 13. Huang--Sellke: exponentially many near-ground states and a quadratic lower tail

**Source.** Brice Huang and Mark Sellke, [*A Constructive Proof of the Spherical Parisi Formula*](https://arxiv.org/abs/2311.15495), revised 2024.

**Hypotheses.** Spherical Gaussian mixed \(p\)-spin model

\[
H_N(\sigma)=\sum_{p\ge1}\frac{\gamma_p}{N^{(p-1)/2}}
\sum g_{i_1\cdots i_p}\sigma_{i_1}\cdots\sigma_{i_p},
\qquad \sum_p2^p\gamma_p^2<\infty.
\]

**Conclusion, geometry.** For any finite chain of Parisi-support overlaps \(q_0<\cdots<q_D\), tolerances \(\delta,\varepsilon>0\), and branching \(k\le e^{cN}\), with probability \(1-e^{-cN}\) there is a \(k\)-ary ultrametric tree whose node energies meet the prescribed Parisi energy profile and whose leaves carry the full free energy in replicated bands. At zero temperature, this yields exponentially many ultrametrically arranged near-ground states.

**Conclusion, rare disorder.** Let \(GS_N=N^{-1}\max H_N\), let \(\xi^{\gamma_1\leftarrow0}\) remove the degree-one/external-field component, and let \(\mathcal Q\) be the spherical zero-temperature Parisi value. For every \(\varepsilon>0\),

\[
\liminf_{N\to\infty}-N^{-2}\log
\mathbb P\!\left(GS_N\le
\mathcal Q(\xi^{\gamma_1\leftarrow0})-\varepsilon\right)
\ge C_2(\xi,\varepsilon)>0.
\]

For 1RSB models without external field, the paper also gives the full upper-tail LDP at speed \(N\), with rate \(-\Theta_*(E)\).

**Mechanism.** Construct exponentially many nearly orthogonal near-maximizers, exactly orthogonalize a linear-size subfamily, and average their energies. The replicated maximum has Gaussian variance proxy \(O(N^{-2})\), so simultaneously suppressing all peaks costs \(e^{-cN^2}\).

**Boundary.** The quadratic lower-tail speed is the closest result here to deliberately finding unusually low-ground-state disorder. But it is spherical and Gaussian, the displayed theorem gives only a positive rate lower bound rather than an exact rate, and no ordinary universality theorem preserves \(N^2\)-speed probabilities.

### 14. Bates--Sohn: balanced nonconvex multi-species comparison

**Source.** Erik Bates and Youngtak Sohn, [*Balanced Multi-Species Spin Glasses*](https://arxiv.org/abs/2507.06522) (2025 preprint).

**Hypotheses.** Finite species set with \(\lambda_{s,N}\to\lambda_s>0\); Ising hypercube or a product of species spheres; Gaussian mixed \(p\)-spin coefficients \(\Delta_{s_1\ldots s_p}\ge0\); analytic decay

\[
\sum_{p\ge1}(1+\varepsilon)^p
\sum_{s_1,\ldots,s_p}\Delta_{s_1\ldots s_p}^2
\lambda_{s_1}\cdots\lambda_{s_p}<\infty.
\]

The balance condition is that \(\Delta_t^2\) is independent of \(t\) for \(p=1\), and, for every \(p\ge2\),

\[
\sum_{s_2,\ldots,s_p}\Delta_{t,s_2,\ldots,s_p}^2
\lambda_{s_2}\cdots\lambda_{s_p}
\]

is independent of the first species \(t\).

**Conclusion.** With

\[
\beta_p^2=\sum_{s_1,\ldots,s_p}\Delta_{s_1\ldots s_p}^2
\lambda_{s_1}\cdots\lambda_{s_p},
\]

the liminf of the multi-species free energy is at least the limiting free energy of the associated single-species model with parameters \((\beta_p)\); the analogous ground-state liminf inequality also holds. Equality follows at high temperature when the single-species model attains its annealed value, at all temperatures when \(\xi\) is convex on the positive orthant, and at zero temperature for specified balanced pure bipartite spherical models.

**Mechanism.** A Guerra-style one-sided interpolation plus a multi-species Talagrand positivity principle. Balance makes the interpolation error have a sign and identifies the optimal diagonal comparison path.

**Boundary.** In the general nonconvex case this paper proves a lower bound, not convergence. It explicitly recorded the bipartite Ising free energy as open in July 2025; the centered-Ising result below resolved that class in June 2026.

### 15. Chen--Issa--Mourrat: nonconvex centered-Ising multi-species limit

**Source.** Hong-Bin Chen, Victor Issa, and Jean-Christophe Mourrat, [*Free Energy of Non-Convex Multi-Species Spin Glasses with Centered Ising Spins*](https://arxiv.org/abs/2606.16636) (June 2026 preprint).

**Hypotheses and normalization.** Finite species set; proportions \(\lambda_N\to\lambda_\infty\in(0,1)^{\mathscr S}\); centered Ising spins; overlap

\[
R_{N,s}(\sigma,\tau)=N^{-1}\sum_{i\in I_{N,s}}\sigma_i\tau_i
\]

(note the denominator \(N\), not \(|I_{N,s}|\)); centered Gaussian field with covariance \(N\xi(R_N)\), where \(\xi:\mathbb R^{\mathscr S}\to\mathbb R\) admits an absolutely convergent power series. No convexity of \(\xi\) is assumed.

The corrected pressure uses probability measure and the paper's minus sign:

\[
\overline F_N(t,0)=-\frac1N\mathbb E\log
2^{-N}\sum_\sigma
e^{\sqrt{2t}H_N(\sigma)-Nt\xi(\lambda_N)}.
\]

Let \(\psi(q)=\sum_s\lambda_{\infty,s}\psi_\circ(q_s)\) be the centered-Ising cascade transform and

\[
\mathcal J_{t,q}(q',p)=
\psi(q')+\langle q-q',p\rangle_{L^2}
+t\int_0^1\xi(p(r))\,dr.
\]

**Conclusion.** For every \(t\ge0\) and \(q\in\mathcal Q_2^{\mathscr S}\),

\[
\lim_{N\to\infty}\overline F_N(t,q)
=\sup_{p\in\mathcal Q_\infty^{\mathscr S}}
\inf_{q'\in\mathcal Q_\infty^{\mathscr S}}
\mathcal J_{t,q}(q',p).
\]

The limit is the Lipschitz viscosity solution of

\[
\partial_t f-\int_0^1\xi(\partial_qf)=0,
\qquad f(0,\cdot)=\psi.
\]

Balanced models reduce to an associated one-species free energy.

**Mechanism.** Prove displacement convexity of the centered-Ising cascade transform; reduce rational species proportions to vector spins; use cavity/Ghirlanda--Guerra calculations to locate critical points; then Hamilton--Jacobi comparison and the Hopf formula force the liminf and limsup to coincide.

**Boundary.** Centered \(\pm1\) spins are essential. The variational statement is false for a one-species biased prior \(p\delta_1+(1-p)\delta_{-1}\) when \(\max(p,1-p)>(3+\sqrt3)/6\approx0.79\), so a deterministic external field is excluded. This remains a quenched Gaussian theorem, not an adversarial-disorder theorem. As of the retrieval date it is a new, unrefereed preprint.

### 16. Chen--Guionnet--Ko--Lacroix-A-Chez-Toine--Mourrat: exact Ising upper-tail disorder LDP

**Source.** Hong-Bin Chen, Alice Guionnet, Justin Ko, Bertrand Lacroix-A-Chez-Toine, and Jean-Christophe Mourrat, [*One-Sided Large Deviations for the Ground-State Energy of Spin Glasses*](https://arxiv.org/abs/2603.06368) (March 2026 preprint).

**Hypotheses.** Ising mixed Gaussian \(p\)-spin field with covariance \(N\xi(R)\), deterministic external field \(h\sum_i\sigma_i\), and a power series \(\xi(r)=\sum_{p\ge2}\beta_p^2r^p\) finite for every real \(r\). Put \(L_N=N^{-1}\max_\sigma H_N(\sigma)\) and \(gs=\lim\mathbb EL_N\).

**Conclusion.** For every \(r\ge gs\),

\[
\lim_{N\to\infty}-\frac1N\log\mathbb P(L_N\ge r)=\Lambda^*(r),
\]

where \(\Lambda^*\) is given explicitly by an infimum over bounded martingales \((\alpha_t)\) with \(|\alpha_1|\le1\) and the tail constraints

\[
\int_t^1\xi''(s)(\mathbb E\alpha_s^2-s)\,ds\ge0
\quad\text{for all }t.
\]

The same admissible martingales give an “un-inverted” variational formula for \(gs\). If \(h\ne0\), \(\Lambda^*(r)\) is globally comparable to \((r-gs)^2\); if \(h=0\), \(\Lambda^*(r)/(r-gs)^2\to\infty\) as \(r\downarrow gs\).

**Mechanism.** Start from a Parisi formula for fractional moments of \(Z_N\), take the zero-temperature Laplace limit, rewrite it through martingale convex duality, and apply Gärtner--Ellis/differentiability arguments.

**Boundary.** This theorem is explicitly one-sided **above** the typical maximum. Disorder minimization asks for the opposite, lower-deviation side, which is expected to have speed \(N^2\). It is Gaussian and, as of the retrieval date, unrefereed.

### 17. Kim: a 2026 universality/extremes phase diagram

**Source.** Taegyun Kim, [*A Sharp Universality Dichotomy for the Free Energy of Spherical Spin Glasses*](https://arxiv.org/abs/2601.08599) (January 2026 preprint).

**Hypotheses.** Pure spherical \(p\)-spin models with iid symmetric disorder whose tails are regularly varying with index \(\alpha\). The Hamiltonian is normalized by a tail-dependent scale built from the \(1/M_{N,p}\) upper quantile, \(M_{N,p}=\binom{N+p-1}{p}\), and by \(\sqrt N\) in the bulk regime.

**Conclusion.** The paper identifies three regimes:

- \(\alpha<2p\): finitely many extreme tensor entries dominate; free energy and ground state have Fréchet/Poisson-extreme limits rather than the Gaussian Parisi value.
- \(\alpha>2p\), with finite \(2p\)-moment: the Gaussian spherical Parisi limit holds.
- \(\alpha=2p\): depending on the ratio of the tail quantile to \(\sqrt N\), the limit is extreme-dominated or a random TAP-type variational mixture of a spike and Gaussian bulk.

**Mechanism.** Poisson point-process convergence for extreme monomials, nonintersection of the leading monomials, spherical slicing/TAP reduction, and the Sawhney--Sellke bulk universality theorem.

**Boundary.** This is spherical and concerns randomly sampled heavy-tailed arrays. Its relevance is diagnostic: deliberately choosing disorder can mimic an exceptional/extreme regime even when the nominal random ensemble is bounded or light-tailed, so product-law universality is not a uniform principle. New, unrefereed preprint.

### 18. Aizenman--Lebowitz--Ruelle: the annealed high-temperature region and its failure

**Source.** Michael Aizenman, Joel Lebowitz, and David Ruelle, [*Some Rigorous Results on the Sherrington--Kirkpatrick Spin Glass Model*](https://doi.org/10.1007/BF01217677), CMP 112 (1987), 3--20; [author-hosted PDF](https://www.ihes.fr/~ruelle/PUBLICATIONS/%5B90%5D.pdf).

**Hypotheses and normalization.** Zero external field, Gaussian SK with

\[
Z_N(\beta)=\sum_\sigma\exp\left(\frac\beta{\sqrt N}\sum_{i<j}g_{ij}\sigma_i\sigma_j\right).
\]

Then

\[
\log\mathbb EZ_N=N\log2+\frac{\beta^2}{4}(N-1).
\]

**Conclusion.** For \(\beta<1\), the total free-energy correction has the nondegenerate Gaussian limit

\[
\log\frac{Z_N}{\mathbb EZ_N}
\Rightarrow
\mathcal N\!\left(
\frac14[\log(1-\beta^2)+\beta^2],
-\frac12[\log(1-\beta^2)+\beta^2]
\right).
\]

In particular the pressure equals its annealed value

\[
\lim_N\frac1N\mathbb E\log Z_N=\log2+\frac{\beta^2}{4}
\]

for \(\beta<1\), and continuity gives the endpoint \(\beta=1\). For \(\beta>1\), the Parisi/replica-symmetry-breaking theory gives a strict gap below this annealed value.

**Mechanism.** Cluster expansion (and, for the pressure equality, a second-moment argument plus concentration). The second moment remains effective exactly to the zero-field SK critical point.

**Boundary and failure lesson.** Jensen always gives quenched \(\le\) annealed, but equality is a phase-specific fact. Above the critical temperature parameter, rare disorder samples dominate \(\mathbb EZ_N\), so annealed interpolation loses the typical free energy. With nonzero external field, the simple quenched-equals-annealed statement is already false. Consequently annealed estimates are especially unreliable for an outer optimization that intentionally searches atypical disorder.

## 3. Mechanism-to-use map

| Need | Strongest theorem/mechanism in this packet | Required input that must be checked |
|---|---|---|
| prove an ordinary quenched thermodynamic limit | Guerra--Toninelli split interpolation | compatible random ensembles; convex covariance or a signed replacement remainder |
| handle an exact discrete constraint across sizes | Panchenko coupled free energy | approximate block constraints; sublinear defect; restricted Fekete lemma |
| express the limit through size increments | Aizenman--Sims--Starr cavity principle | fixed cavity size; random reservoir; overlap-kernel positivity |
| move Gaussian to Rademacher iid disorder | Chatterjee Lindeberg invariance | independence, moment matching, Lindeberg tails; only typical/averaged output |
| pass pressure to a ground state | uniform log-sum-exp sandwich; Auffinger--Chen | control uniform enough to send \(\beta\to\infty\) |
| treat multiple species | Panchenko 2015; Chen--Issa--Mourrat 2026 | PSD covariance in the classical theorem; centered Ising in the nonconvex theorem |
| treat heavy tails | Jagannath--Lopatto sparse interpolation | \(\alpha\in(1,2)\), \(N^{-1/\alpha}\) scale; no Gaussian universality |
| understand a rare low-ground-state sample | Huang--Sellke quadratic lower-tail bound | spherical Gaussian model; only a positive rate bound |
| compute an exact disorder LDP | 2026 one-sided Ising theorem | upper tail only; Gaussian mixed \(p\)-spin |
| interpret a “min--max” spin-glass formula | Chen 2024/2025, Chen--Issa--Mourrat 2026 | optimization variables are order parameters, not couplings |

## 4. Impossibility and non-transfer lessons to keep explicit

1. **Expectation plus concentration is not uniformity.** Even exponentially good \(e^{-cN}\) concentration leaves exponentially many exceptional coupling arrays. The deterministic disorder class has \(e^{\Theta(N^2)}\) elements.
2. **Ordinary universality is too low-resolution for optimized disorder.** \(o(1)\) comparison of expected normalized maxima says nothing about probabilities of order \(e^{-cN^2}\), precisely the scale on which lower-ground-state disorder can live.
3. **A Parisi saddle point is not an adversarial game.** The variables are overlap paths, magnetizations, or self-overlaps induced by the Gibbs measure. They do not parameterize individual \(J_{ij}\)'s.
4. **Convexity failures are genuine.** Classical Guerra upper bounds can be wrong, not merely unproved, for nonconvex covariance. The 2026 centered-Ising repair uses displacement convexity of a special cascade transform and demonstrably fails for sufficiently biased spins.
5. **Annealed calculations can be governed by the wrong samples.** Above the SK critical point, \(\log\mathbb EZ\) is strictly larger than \(\mathbb E\log Z\). An outer minimization is also a rare-sample operation, but in the opposite direction; neither annealed nor typical quenched values identify it.
6. **Absolute energy couples \(J\) and \(-J\).** Sign symmetry identifies their separate laws under random sampling, not the joint maximum after a common deterministic choice.
7. **Size interpolation needs a coherent extension rule.** Quenched proofs freely resample independent block disorders. A deterministic optimizer at size \(N\) need not restrict to optimizers at \(N_1,N_2\), and optimizers at smaller sizes need not admit a low-cost completion. Without a deterministic gluing/restriction inequality with \(o(N^{3/2})\) error, Fekete-style conclusions do not follow.
8. **Spherical rare-disorder theorems are informative but not interchangeable with Ising/Rademacher.** Sphere geometry supplies exact orthogonalization and rotational invariance; both are absent on the hypercube.

## 5. Short priority reading order

1. Guerra--Toninelli 2002 for the exact system-size derivative and the zero-temperature sandwich.
2. Panchenko 2004 for near-superadditivity under incompatible discrete constraints.
3. Chatterjee 2005 for the strongest clean Gaussian-to-Rademacher transfer and its quantifier limitation.
4. Aizenman--Sims--Starr 2003 for the cavity increment principle.
5. Chen--Issa--Mourrat 2026 for the current nonconvex centered-Ising frontier.
6. Huang--Sellke 2024 and the 2026 one-sided LDP paper for the two different disorder-deviation speeds and directions.
