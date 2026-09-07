# Spin-glass domain report: deterministic antipodal pressure gluing

This report was frozen without consulting the project archive. It uses one
architecture only. No rare-disorder LDP is invoked.

## 1. Native translation and normalization

Let \(\mathcal A_n\) be the hollow symmetric \(\{\pm1\}\) arrays and put

\[
 h_A(\sigma)=\sum_{i<j}a_{ij}\sigma_i\sigma_j,
 \qquad
 \mathcal H_{n,A}(\sigma)=\frac{h_A(\sigma)}{\sqrt n}.
\]

The one-sided SK ground-state density of the deterministic sample \(A\) is

\[
 G_n^+(A)=\frac1n\max_\sigma \mathcal H_{n,A}(\sigma),
\]

whereas the problem asks for the antipodal, or two-channel, density

\[
 G_n^{\pm}(A)
 =\max\{G_n^+(A),G_n^+(-A)\}
 =\frac{Q(A)}{n^{3/2}}.
\]

Thus, with \(g_n=M_n/n^{3/2}\), the exact order of play is

\[
 g_n
 =\min_{A\in\mathcal A_n}\ \max_{s\in\{\pm1\},\,\sigma\in\{\pm1\}^n}
 \frac{s\,\mathcal H_{n,A}(\sigma)}{n}.
\]

The disorder is chosen first and deterministically. Neither an expectation
over \(A\), nor \(\max_\sigma\min_A\), is the same object.

For iid mean-zero variance-one disorder, the exact covariance in this
\(i<j\) convention is

\[
 \mathbb E\,\mathcal H_{n,A}(\sigma)\mathcal H_{n,A}(\tau)
 =\frac n2R(\sigma,\tau)^2-\frac12,
 \qquad
 R(\sigma,\tau)=\frac1n\sum_i\sigma_i\tau_i.
\]

Hence the limiting SK covariance is \(\xi(r)=r^2/2\), not \(r^2\).
Quenched Rademacher universality concerns a typical iid \(A\) under this
covariance and does not commute with the outer minimum.

Thermalize only the inner antipodal maximum. Define the normalized
antipodal partition function

\[
 \widehat Z_n(A;\beta)
 =2^{-n}\sum_{\sigma\in\{\pm1\}^n}
 \cosh\!\left(\frac{\beta h_A(\sigma)}{\sqrt n}\right),
 \qquad \beta>0,
\]

and the adversarial pressure numerator and density

\[
 a_n(\beta)=\min_{A\in\mathcal A_n}\log\widehat Z_n(A;\beta),
 \qquad
 \phi_n(\beta)=\frac{a_n(\beta)}n.
\]

The minimum over \(A\) remains a zero-temperature deterministic minimum;
\(\beta\) softens only \(\max_{s,\sigma}\).

If \(E_n=\binom n2\), \(a\) is the sign vector of \(A\), and
\(\mathcal C_n^+\) is the antipodal cut code, then every codeword has two
preimages among \((s,\sigma)\), and

\[
 \widehat Z_n(A;\beta)
 =2^{-n}\sum_{c\in\mathcal C_n^+}
 \exp\!\left(\frac\beta{\sqrt n}\langle a,c\rangle\right)
 =2^{-n}\sum_{c\in\mathcal C_n^+}
 \exp\!\left(\frac\beta{\sqrt n}(E_n-2d(a,c))\right).
\]

Thus \(\widehat Z_n\) is one Laplace transform of the coset-distance
enumerator, while

\[
 Q(A)=E_n-2d(a,\mathcal C_n^+),
 \qquad
 M_n=E_n-2\rho(\mathcal C_n^+).
\]

Finally, for every deterministic \(A\),

\[
 2^{-n}e^{\beta Q(A)/\sqrt n}
 \le \widehat Z_n(A;\beta)
 \le e^{\beta Q(A)/\sqrt n}.
\]

After minimizing over \(A\), this is the exact sandwich

\[
 \boxed{
 g_n-\frac{\log2}{\beta}
 \le \frac{\phi_n(\beta)}\beta
 \le g_n.}
 \tag{1}
\]

## 2. The one frozen architecture

Use a deterministic analogue of Guerra--Toninelli splitting, but formulate
it as Laplace-transform order rather than as a pointwise energy comparison.
Given \(A\in\mathcal A_m\) and \(B\in\mathcal A_k\), retain \(A\) on the
first principal block, retain either \(B\) or \(-B\) on the second block,
and choose the \(mk\) cross signs. The relative sign is legitimate because
\(\widehat Z_k(B;\beta)=\widehat Z_k(-B;\beta)\). The desired completion
has antipodal Laplace transform at most the product of the two block
transforms, up to a summably sublinear defect.

This differs from a pointwise discrepancy gluing: it need not control every
spin configuration. It also differs from a quenched interpolation: the
principal blocks are arbitrary prescribed deterministic arrays and the
endpoint cross couplings must be binary.

There is an unavoidable quantifier boundary. If the lemma is required only
for thermal minimizers \(A_m(\beta),A_k(\beta)\), it is just the desired
near-subadditive recurrence for \(a_n(\beta)\) with the witnesses written
out. The noncircular version below is therefore uniform over all prescribed
block arrays. This is much stronger, and is the main reason the proposal is
currently low-confidence.

## 3. Best exact implication

Let the imported nearly-subadditive theorem be Theorem 5 of Zoltan
Furedi and Imre Ruzsa, [*Nearly subadditive
sequences*](https://arxiv.org/abs/1810.11723): if

\[
 b_{m+k}\le b_m+b_k+f(m+k)
\]

for all sufficiently large comparable pairs \(m\le k\le\mu m\), where
\(f\) is nonnegative, nondecreasing, and
\(\sum_{N\ge1}f(N)/N^2<\infty\), then \(b_n/n\) converges. This is the
de Bruijn--Erdos extension of Fekete's lemma. Panchenko's
[*A Note on the Free Energy of the Coupled System in the
Sherrington--Kirkpatrick Model*](https://arxiv.org/abs/math/0405359)
is the spin-glass precedent for paying a sublinear split defect and then
using this type of restricted Fekete argument.

The proposed implication is

\[
 \boxed{
 \text{de Bruijn--Erdos/Furedi--Ruzsa near-subadditivity}
 + L_{\mathrm{Lap}}
 + \text{the exact sandwich (1)}
 \Longrightarrow
 \lim_{n\to\infty}\frac{M_n}{n^{3/2}}\ \text{exists}.}
\]

Indeed, \(L_{\mathrm{Lap}}\) applied to minimizers of the two block
pressures gives

\[
 a_{m+k}(\beta)
 \le a_m(\beta)+a_k(\beta)+C_\beta(m+k)^{1-\delta_\beta}.
 \tag{2}
\]

Since
\(\sum_N N^{1-\delta_\beta}/N^2
=\sum_NN^{-1-\delta_\beta}<\infty\), the imported theorem yields

\[
 \ell(\beta):=\lim_{n\to\infty}\phi_n(\beta)
\]

for every fixed \(\beta>0\). From (1),

\[
 \liminf_ng_n\ge\frac{\ell(\beta)}\beta,
 \qquad
 \limsup_ng_n\le\frac{\ell(\beta)}\beta+\frac{\log2}{\beta}.
\]

Therefore

\[
 0\le\limsup_ng_n-\liminf_ng_n\le\frac{\log2}{\beta}
 \quad\text{for every }\beta>0.
\]

Sending \(\beta\to\infty\) proves convergence. Notice the order:
first \(n\to\infty\) at each fixed \(\beta\), then \(\beta\to\infty\).
No uniform bound on \(C_\beta\) or \(\delta_\beta\) as
\(\beta\to\infty\) is needed. Conversely, a result only for
\(0<\beta\le\beta_0<\infty\) leaves the nonzero gap
\((\log2)/\beta_0\) and is insufficient.

## 4. Exact missing lemma \(L_{\mathrm{Lap}}\)

For \(A\in\mathcal A_m\), \(B\in\mathcal A_k\), a relative orientation
\(\eta\in\{\pm1\}\), and a rectangular sign matrix
\(D\in\{\pm1\}^{m\times k}\), write

\[
 A\oplus_D(\eta B)
 =\begin{pmatrix}A&D\\D^{\mathsf T}&\eta B\end{pmatrix}
 \in\mathcal A_{m+k}.
\]

The fully quantified lemma is:

\[
\boxed{
\begin{minipage}{0.91\linewidth}
For every fixed \(\beta\in(0,\infty)\), there exist constants
\(C_\beta<\infty\), \(\delta_\beta\in(0,1]\), and
\(N_\beta<\infty\) such that, for every pair of integers
\(N_\beta\le m\le k\le2m\), and every pair of deterministic arrays
\(A\in\mathcal A_m\), \(B\in\mathcal A_k\), there exist
\(\eta=\eta(A,B,\beta)\in\{\pm1\}\) and
\(D=D(A,B,\beta)\in\{\pm1\}^{m\times k}\) for which
\[
 \log\widehat Z_{m+k}(A\oplus_D(\eta B);\beta)
 \le
 \log\widehat Z_m(A;\beta)+\log\widehat Z_k(B;\beta)
 +C_\beta(m+k)^{1-\delta_\beta}.
\]
The choices of \(\eta,D\) may depend on \(\beta\); no common completion
for two different temperatures is asserted.
\end{minipage}}
\tag{L_{\mathrm{Lap}}}
\]

The range \(m\le k\le2m\) is exactly enough for the restricted
near-subadditivity theorem. Replacing the displayed defect by an arbitrary
\(o(m+k)\) is not enough: the de Bruijn--Erdos summability condition is
essential in general. The power saving supplies it explicitly.

The lemma is required at every fixed temperature, including the entire
low-temperature/replica-symmetry-breaking regime. It already contains the
joint \(A,-A\) channel through \(\cosh\); two separately chosen one-sided
completions do not imply it.

## 5. Why the lemma is strictly weaker than full deterministic or coset optimization

The data used by \(L_{\mathrm{Lap}}\) at a fixed temperature are a
one-sided inequality for one exponential moment and the existence of one
completion. Full coset optimization asks, for every \(A\), for the extreme
coefficient \(d(a,\mathcal C_n^+)\); the full histogram asks for every
coefficient of the distance enumerator. These are not recoverable from the
fixed-temperature datum.

Here is an exact non-recovery witness inside the present coupling class,
not a description-length argument. Order the edges of \(K_6\)
lexicographically. Let \(A_9\) have all edges negative except
\(45,46,56\), and let \(A_{11}\) have all edges negative except
\(36,45\). Direct enumeration of the 64 spins gives the absolute-energy
histograms

\[
 \begin{array}{c|ccc}
 A_9:&|h|=1&7&9\\
     &48&12&4
 \end{array},
 \qquad
 \begin{array}{c|ccccc}
 A_{11}:&|h|=1&3&5&7&11\\
        &28&16&14&4&2.
 \end{array}
\]

Consequently \(Q(A_9)=9\) and \(Q(A_{11})=11\), and their augmented
coset histograms differ. Nevertheless, at

\[
 \beta_*=\sqrt6\,
 \operatorname{arcosh}\sqrt{\frac{5+\sqrt{17}}8},
\]

they have exactly the same \(\widehat Z_6\). To verify this, set
\(x=\beta/\sqrt6\) and \(c=\cosh x\). Sixty-four times the difference of
their normalized antipodal partition functions is

\[
 20\cosh x-16\cosh3x-14\cosh5x
 +8\cosh7x+4\cosh9x-2\cosh11x
 =-256c^3(c^2-1)^2(8c^4-10c^2+1),
\]

which vanishes at \(c^2=(5+\sqrt{17})/8\). Hence even the exact value of
the fixed-temperature Laplace statistic does not determine the maximum,
much less the whole histogram.

Knowing \(\widehat Z_A(\beta)\) as an analytic function of \(\beta\) for
the same fixed \(A\) would determine its finite energy histogram. The lemma
does not give that information: its completion may change with \(\beta\),
and it gives only an upper inequality. This temperature-dependent
existential quantifier is precisely what makes the lemma strictly weaker.
It also yields no exact finite value of \(a_n(\beta)\) or \(M_n\), only a
near-subadditive upper comparison sufficient for existence of a limit.

## 6. Decisive falsifier

For fixed \(\beta,m,k,A,B\), define the optimal completion defect

\[
 \Delta_{m,k}^{\beta}(A,B)
 =\min_{\eta\in\{\pm1\},\,D\in\{\pm1\}^{m\times k}}
 \left[
 \log\widehat Z_{m+k}(A\oplus_D(\eta B);\beta)
 -\log\widehat Z_m(A;\beta)
 -\log\widehat Z_k(B;\beta)
 \right]
\]

and

\[
 \Gamma_{m,k}(\beta)=\max_{A\in\mathcal A_m,\,B\in\mathcal A_k}
 \Delta_{m,k}^{\beta}(A,B).
\]

Then \(L_{\mathrm{Lap}}\) is falsified decisively by any
\(\beta_0>0\), \(c>0\), and balanced sequence
\(m_j\le k_j\le2m_j\), \(m_j\to\infty\), for which

\[
 \Gamma_{m_j,k_j}(\beta_0)\ge c(m_j+k_j).
 \tag{3}
\]

This is an exact structural obstruction: a fixed linear completion cost
cannot satisfy any power-saving defect. At finite orders, all quantities in
\(\Delta\) are finite sums of exponentials of integers, so exhaustive
enumeration (with switching and permutation quotients if desired) gives a
certificate against any proposed numerical triple
\((C_{\beta_0},\delta_{\beta_0},N_{\beta_0})\). A single finite violation
cannot refute the existential constants; a family proving (3) can.

A useful warning for this test is the exact random-cross-sign identity. If
\(D\) is uniform and the principal blocks are fixed, then

\[
 \mathbb E_D\widehat Z_{m+k}(A\oplus_D(\eta B);\beta)
 =\cosh\!\left(\frac\beta{\sqrt{m+k}}\right)^{mk}
 2^{-(m+k)}\sum_{\sigma,\tau}
 \cosh\!\left(
 \frac\beta{\sqrt{m+k}}(h_A(\sigma)+\eta h_B(\tau))
 \right).
\]

For balanced blocks the displayed cross factor has logarithm
\(\beta^2mk/(2(m+k))+O(1)=\Theta(m+k)\). A proof by averaging must cancel
this linear annealed cost against the dilution of the two block
Hamiltonians. Retaining any positive linear residue is fatal; it does not
become a Fekete-negligible error.

## 7. Closest imported theorem and unmatched hypotheses

The closest spin-glass theorem is Guerra and Toninelli,
[*The Thermodynamic Limit in Mean Field Spin Glass
Models*](https://arxiv.org/abs/cond-mat/0204280). Their split interpolation
compares a full Gaussian SK sample with two independently resampled block
samples. Gaussian integration by parts turns the derivative into the
convex-overlap remainder and proves exact superadditivity of the quenched
pressure. Their non-Gaussian extension allows a sublinear averaged error.

Every unmatched hypothesis is material here:

1. Guerra--Toninelli average over compatible random ensembles;
   \(L_{\mathrm{Lap}}\) is uniform over every prescribed deterministic
   pair \(A,B\).
2. Their full disorder and two block disorders are freely and independently
   resampled. Here the two principal blocks must be retained and only the
   cross signs may be selected.
3. Gaussian integration by parts, or its averaged symmetric-disorder
   replacement, creates a covariance expression. A deterministic sign
   selection has no such identity.
4. Convexity signs an expected replica-overlap remainder. It does not sign
   the derivative for a single adversarial array.
5. The endpoint of the proposed interpolation must lie in
   \(\{\pm1\}^{\binom{m+k}{2}}\); a Gaussian or continuously weighted
   endpoint is not admissible.
6. The present partition function contains the common antipodal channel
   \(\cosh\). Separate one-sided comparisons for \(A\) and \(-A\) may
   choose different cross completions and do not prove the joint bound.
7. The conclusion is needed for every fixed \(\beta>0\), not merely in the
   annealed/high-temperature phase.
8. The deterministic defect must satisfy the de Bruijn--Erdos summability
   threshold. A fixed \(c_\beta(m+k)\) loss, however small its coefficient,
   is insufficient.

Panchenko's coupled-system theorem is the closest precedent for a
sublinear-defect restricted-Fekete step, but its constraint is on replica
overlap; Gaussian disorder remains independently sampled. Carmona--Hu and
Chatterjee universality transfer the typical iid Rademacher pressure and
ground state, not the uniform completion statement. None supplies a
missing hypothesis above.

## 8. Circularity audit

The following strengthenings or weakenings would be circular or invalid.

- Requiring the completion lemma only when \(A,B\) minimize their block
  pressures is equivalent to assuming (2), the recurrence one needs. The
  uniform all-block quantifier in \(L_{\mathrm{Lap}}\) must not be silently
  dropped.
- Selecting a different completion for the \(+A\) and \(-A\) channels
  replaces the absolute maximum by two unrelated one-sided problems.
- Assuming that a typical iid cross block works imports annealed/quenched
  information at the wrong quantifier. The deterministic minimum may live
  among exponentially rare completions.
- Assuming the same completion works for all \(\beta\), then taking
  \(\beta\to\infty\), would amount to a much stronger zero-temperature
  gluing recurrence. It is not part of the lemma.
- Conversely, allowing the completion to depend on \(\beta\) is harmless
  for the implication because the limit is first taken at each fixed
  temperature. One may not use those temperature-dependent witnesses to
  claim a coherent optimizer.
- Replacing the binary completion by an interpolated Gaussian or real
  matrix does not return an admissible coupling array without a new
  rounding theorem at \(o(n)\) pressure cost.
- An \(o(n)\) defect without the summability condition is not enough in
  general; nearly subadditive sequences can still have nonconvergent
  slopes.
- A strict-high-temperature proof cannot be extrapolated through the phase
  transition. The zero-temperature sandwich requires all fixed
  temperatures.
- Computing the entire coset enumerator to verify the lemma would defeat
  its advertised strictness. A viable proof must establish the Laplace
  inequality directly, for example through a deterministic signing or
  convex-order argument.

The audit leaves a genuine, noncircular lemma, but no currently imported
spin-glass theorem proves it. In particular, this report does not claim a
viable proof of convergence.

## 9. Confidence and candidate card

| Field | Entry |
|---|---|
| Domain | Mean-field spin glasses; deterministic antipodal SK pressure |
| Imported theorems | Guerra--Toninelli split interpolation (mechanism); de Bruijn--Erdos/Furedi--Ruzsa nearly-subadditive theorem (implication); Panchenko restricted-Fekete precedent |
| Problem translation | \(g_n=\min_A\max_{s,\sigma}s\mathcal H_{n,A}(\sigma)/n\), with \(\xi(r)=r^2/2\) for random disorder and exact sandwich (1) |
| Proposed mechanism | All-temperature deterministic Laplace-order block completion with summably sublinear pressure defect |
| Exact missing lemma | \(L_{\mathrm{Lap}}\) in Section 4 |
| Why strictly weaker | One fixed-temperature moment does not recover \(Q\): the explicit order-six arrays above have equal moment at \(\beta_*\) but maxima 9 and 11; witnesses may change with \(\beta\) |
| Falsification test | Prove \(\Gamma_{m_j,k_j}(\beta_0)\ge c(m_j+k_j)\) on a balanced sequence |
| Specialist confidence | Truth of \(L_{\mathrm{Lap}}\): 0.08; tractability: 0.03; correctness of the conditional implication: 0.99 |

Overall judgment: **hold, presently non-viable**. The implication is exact,
and the lemma is genuinely less informative than full optimization, but the
only clearly noncircular quantifier is a worst-case deterministic completion
statement far beyond existing quenched interpolation or universality.
