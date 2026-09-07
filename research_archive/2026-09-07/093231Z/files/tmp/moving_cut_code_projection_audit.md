# Moving-projection certificates for augmented cut-code cosets

Date: 2026-08-15. This is an agent-authored research draft. It does not
modify the project ledger or steering file.

## 1. Verdict

The moving-projection theorem in Chapter 2 of OpenAI's
[*Ten Advances in Mathematics and Theoretical Computer Science*](https://cdn.openai.com/pdf/ten-proofs-oai.pdf)
is a packing theorem. Applied directly to the augmented cut code, it sees
only the code's minimum distance, which is \(n-1\), and is invariant under
translation to a coset. It therefore does **not** itself bound the covering
radius or the extreme weight in a coset.

There is, however, an exact coset inequality that genuinely retains the two
mixed channels in the paper's complete Gram remainder. For every moving
projection kernel it gives

\[
  (\lambda-\mu(a))_+T_a\le J_{\mathcal C_n^+},
  \tag{1.1}
\]

where \(\mu(a)=Q(a)/E_n\), \(T_a\) is one weighted coset enumerator, and
\(J_{\mathcal C_n^+}\) is an explicit quantity determined solely by the
known weight distribution of the cut code. This is a rigorous inequality at
the correct possible leading scale: choosing representation degree
\(L=\Theta(n)\) makes \(\lambda=\Theta(n^{-1/2})\), hence
\(E_n\lambda=\Theta(n^{3/2})\).

The obstruction is now precise. At that degree, \(T_a\) is a signed weighted
sum over even Eulerian subgraphs with up to \(2L=\Theta(n)\) edges. A uniform
lower bound for this one joint sum is not supplied by the packing theorem.
The first finite test of the paper's rank-one moving family finds that every
tested higher-level kernel is worse than the level-one Parseval bound on the saved
exact minimizers through order 14 and on the order-18 conference signing.
Thus the new theorem demonstrates a viable *architecture*, but the tested
published kernels do not yet remove the project's growing-hierarchy obligation.

## 2. Exact cut-code mapping

Put

\[
  E=E_n=\binom n2,
  \qquad \Omega_E=\{\pm1\}^{E},
  \qquad \tau(z)=\frac1E\sum_{e\in E(K_n)}z_e.
\]

The augmented cut code in sign coordinates is

\[
  \mathcal C_n^+
  =\{(\sigma x_ix_j)_{i<j}:\sigma,x_i\in\{\pm1\}\}
  \le \Omega_E.
  \tag{2.1}
\]

It has \(2^n\) words. For a signing \(a\in\Omega_E\), define

\[
  \mu(a)=\max_{c\in\mathcal C_n^+}\tau(ac).
\]

Because the all-minus word belongs to \(\mathcal C_n^+\), this is exactly

\[
  \mu(a)
  =\frac1E\max_{x\in\{\pm1\}^n}
       \left|\sum_{i<j}a_{ij}x_ix_j\right|
  =\frac{Q(a)}E.
  \tag{2.2}
\]

Equivalently,

\[
  d(a,\mathcal C_n^+)=\frac{E-Q(a)}2,
  \qquad
  M_n=E-2\rho(\mathcal C_n^+).
  \tag{2.3}
\]

For \(n\ge5\), the minimum distance of the code is \(n-1\). Indeed, a
nonconstant cut has weight \(r(n-r)\ge n-1\), while a complemented cut has
weight \(E-r(n-r)\ge E-\lfloor n^2/4\rfloor\ge n-1\). A singleton cut
attains \(n-1\).

Consequently, a direct use of Theorem 2.1 of the paper has ambient length
\(E\), distance \(n-1\), and

\[
  s=1-\frac{2(n-1)}E=1-\frac4n.
  \tag{2.4}
\]

It produces an upper bound on the size \(2^n\) from these *internal*
distances. Translating \(\mathcal C_n^+\) by \(a\) preserves every internal
distance, so that conclusion is identical for every coset. The quantity in
(2.2), by contrast, is the distance from the external origin to that coset.
This is the packing-versus-covering mismatch.

## 3. What the new joint Gram identity really supplies

Use the paper's construction at arbitrary ambient Hamming length \(E\), with
parameters \(0\le k<L\le E-k\). It attaches to each word \(z\in\Omega_E\)
an equal-rank orthogonal projection \(P_z\). Its scalar overlap is

\[
  K(z,w)=\operatorname{tr}(P_zP_w)=\|P_zP_w\|_{HS}^2\ge0.
\]

Translation and permutation equivariance make this a radial kernel, so write

\[
  K(z,w)=\kappa(zw).
\]

Most importantly, equations (25)--(27) of the paper construct maps
\(\Theta_z\) for which

\[
  \langle\Theta_z,\Theta_w\rangle_{HS}
  =\bigl(\tau(zw)-\lambda\bigr)\kappa(zw),
  \tag{3.1}
\]

where \(\lambda\) is the top eigenvalue of the corrected transition matrix.
Thus both

\[
  \kappa(z)\quad\text{and}\quad
  F(z):=(\tau(z)-\lambda)\kappa(z)
  \tag{3.2}
\]

are positive-definite functions on \(\Omega_E\), and \(\kappa(z)\ge0\)
pointwise. The
[*Mathematical Discovery Notes*](https://cdn.openai.com/pdf/reasoning-walkthroughs.pdf)
emphasize that (3.1) is obtained only after
retaining and cancelling both mixed terms; individual harmonic channels are
not being declared positive. The argument below uses (3.1) as one complete
Gram remainder and never pays its channels separately.

## 4. A rigorous moving-projection coset inequality

Define

\[
  T_a=\sum_{c\in\mathcal C_n^+}\kappa(ac),
  \qquad
  J_{\mathcal C_n^+}
  =\sum_{c\in\mathcal C_n^+}F(c).
  \tag{4.1}
\]

> **Theorem 4.1 (joint-Gram coset inequality).** For every signing
> \(a\in\Omega_E\),
> \[
>   (\lambda-\mu(a))_+T_a\le J_{\mathcal C_n^+}.
>   \tag{4.2}
> \]
> In particular, if \(T_a>0\), then
> \[
>   \boxed{\mu(a)\ge
>   \lambda-\frac{J_{\mathcal C_n^+}}{T_a}.}
>   \tag{4.3}
> \]

**Proof.** Use the normalized Fourier expansion

\[
  F(z)=\sum_{\xi\in\widehat\Omega_E}\widehat F(\xi)\chi_\xi(z).
\]

Because \(F\) is the complete Gram kernel in (3.1), Bochner's theorem on the
finite abelian group gives \(\widehat F(\xi)\ge0\) for every \(\xi\). If
\((\mathcal C_n^+)^\perp\) is the annihilator of the code, then

\[
\begin{aligned}
  S_a:=\sum_{c\in\mathcal C_n^+}F(ac)
  &=|\mathcal C_n^+|
    \sum_{\xi\in(\mathcal C_n^+)^\perp}
      \widehat F(\xi)\chi_\xi(a),\\
  J_{\mathcal C_n^+}
  &=|\mathcal C_n^+|
    \sum_{\xi\in(\mathcal C_n^+)^\perp}\widehat F(\xi).
\end{aligned}
\tag{4.4}
\]

Hence \(-J_{\mathcal C_n^+}\le S_a\le J_{\mathcal C_n^+}\), and in
particular \(J_{\mathcal C_n^+}\ge0\). On the other hand,
\(\tau(ac)\le\mu(a)\) and \(\kappa(ac)\ge0\), so

\[
  S_a
  =\sum_c(\tau(ac)-\lambda)\kappa(ac)
  \le(\mu(a)-\lambda)T_a.
  \tag{4.5}
\]

If \(\mu(a)<\lambda\), combining (4.4) and (4.5) proves (4.2). If
\(\mu(a)\ge\lambda\), (4.2) is automatic. This also proves (4.3). \(\square\)

This inequality is not obtained by summing the packing theorem over the
coset. Indeed,

\[
  \sum_{c,d\in a\mathcal C_n^+}F(cd)
  =|\mathcal C_n^+|J_{\mathcal C_n^+}
  \tag{4.6}
\]

is independent of \(a\). The location-sensitive information is exactly the
cross term \(S_a\), and Fourier positivity bounds it by the internally
computable budget \(J_{\mathcal C_n^+}\).

### 4.1 Comparison with the rooted annular full-Gram inequality

There is a generic rooted inequality using the same unsplit Gram remainder.
Let \(r=\kappa(1)=\operatorname{rank}P_z\) and

\[
  m_s=\min\{\kappa(z):|\tau(z)|\le s\}.
\]

If \(\mu(a)\le s<\lambda\), then every root-to-coset Gram entry is at most
\(-(\lambda-s)m_s\). Cauchy--Schwarz applied to
\(\Theta_1\) and \(\sum_c\Theta_{ac}\) gives

\[
  (\lambda-s)^2|\mathcal C_n^+|m_s^2
  \le(1-\lambda)rJ_{\mathcal C_n^+}.
  \tag{4.7}
\]

The cut-specific Fourier argument above is stronger in the relevant regime.
Indeed \(T_a\ge|\mathcal C_n^+|m_s\), so (4.2) gives

\[
  (\lambda-s)|\mathcal C_n^+|m_s
  \le J_{\mathcal C_n^+}.
  \tag{4.8}
\]

For \(\lambda<1/2\), one has
\((\lambda-s)m_s\le(1-\lambda)r\); hence the lower bound on
\(J/|\mathcal C_n^+|\) in (4.8) dominates the one obtained by squaring in
(4.7). The improvement comes from subgroup Fourier positivity:
\(|S_a|\le J\), rather than generic Gram Cauchy--Schwarz. Thus the rooted
annular inequality and (4.2) are not competing mechanisms; (4.2) is its
coset/subgroup strengthening. The same denominator collapse found below
therefore obstructs the rooted version at least as strongly for these
kernels.

## 5. Every numerator is explicit

For a vertex sign vector with \(r\) negative coordinates, the corresponding
unaugmented cut word has normalized edge correlation

\[
  q_r=\frac{(n-2r)^2-n}{n(n-1)}.
  \tag{5.1}
\]

Adding both global signs and correcting the \(x\sim-x\) double count gives

\[
\begin{aligned}
J_{\mathcal C_n^+}
=\frac12\sum_{r=0}^n\binom nr\bigl[&
  (q_r-\lambda)\kappa(q_r)\\
  &+(-q_r-\lambda)\kappa(-q_r)\bigr].
\end{aligned}
\tag{5.2}
\]

Thus the numerator of (4.3) depends only on \(n\) and the chosen transition
graph, not on an unknown optimizer or coset.

The denominator is also an exact single scalar:

\[
  T_a=\frac12\sum_{x\in\{\pm1\}^n}
  \left[
    \kappa\!\left(\frac{H_a(x)}E\right)
    +\kappa\!\left(-\frac{H_a(x)}E\right)
  \right].
  \tag{5.3}
\]

It is not the full maximum. Nevertheless, controlling (5.3) uniformly is
the remaining mathematical obligation; merely renaming it is not progress.

## 6. Exact relation to the Eulerian hierarchy

The dual code is

\[
  (\mathcal C_n^+)^\perp
  =\{F\subseteq E(K_n):\deg_F(v)\text{ is even for every }v,
                         \ |F|\text{ is even}\}.
  \tag{6.1}
\]

Since \(\kappa\) is positive definite, all its Walsh coefficients are
nonnegative. Therefore

\[
  \frac{T_a}{|\mathcal C_n^+|}
  =\sum_{F\in(\mathcal C_n^+)^\perp}
       \widehat\kappa(F)\prod_{e\in F}a_e.
  \tag{6.2}
\]

For the whole-cube construction using Fourier levels \(k,\ldots,L\), the
matrix entries of the moving embedding have degrees at most \(L\), and
\(K=\|\Psi_z^*\Psi_w\|_{HS}^2\) has Walsh degree at most \(2L\). Hence
(6.2) is a nonnegative radial weighting of the signed even-Eulerian
coefficients through size \(2L\). This is exactly the hierarchy exposed in
artifacts/eulerian_free_energy_identity.md, but with weights generated by
the moving representation rather than the product kernel \(\rho^{|F|}\).

The corrected transition entries satisfy

\[
  c_i^{(k)}
  =\frac{(i-k+1)(E-i-k)}{E\sqrt{(i+1)(E-i)}}
  \le\sqrt{\frac{L+1}{E}}
  \qquad(k\le i<L).
\]

The spectral radius of the path is at most twice its largest edge weight, so

\[
  \lambda\le2\sqrt{\frac{L+1}{E}}.
  \tag{6.3}
\]

Conversely, for the rank-one case \(k=0\), the final \(2\times2\) principal
block gives
\(\lambda\ge c_{L-1}^{(0)}=\Theta(\sqrt{L/E})\) whenever
\(L=o(E)\). Thus \(L=\Theta(n)\) genuinely attains, rather than merely permits,
the scale \(\lambda=\Theta(n^{-1/2})\).

Therefore any certificate with \(\lambda\ge c/\sqrt n\) must use

\[
  L+1\ge \frac{c^2E}{4n}=\Theta(n).
  \tag{6.4}
\]

This proves that the correct \(n^{3/2}\) scale necessarily opens a
linear-size Eulerian hierarchy. The moving construction compresses that
hierarchy to the one scalar (6.2), which is a genuine improvement in
description, but no published inequality gives the uniform lower bound on
that scalar needed in (4.3).

## 7. Finite test of the paper's rank-one family

For \(k=0\), let \(v\) be the positive Perron vector of the path with edges

\[
  c_i=\frac{\sqrt{(i+1)(E-i)}}E,
  \qquad 0\le i<L.
\]

Set

\[
  p_i=\frac{\sqrt{\binom Ei}\,v_i}
            {\sum_{j=0}^L\sqrt{\binom Ej}\,v_j},
  \qquad
  g_L(t)=\sum_{i=0}^Lp_i
  \frac{K_i(E(1-t)/2)}{\binom Ei}.
  \tag{7.1}
\]

Then the rank-one projection overlap is \(\kappa_L(t)=g_L(t)^2\).
The accompanying script evaluated (5.2)--(5.3) for all
\(1\le L\le\min(E-1,3n)\). The inputs were the saved exact minimizers for
\(n=6,8,10,12,14\), including the conference-optimal orders, and the saved
order-18 conference energy histogram.

| \(n\) | signing cap \(Q(a)\) | best bound \(E(\lambda-J/T)\) | best \(L\) |
|---:|---:|---:|---:|
| 6 | 5 | 3.872983 | 1 |
| 8 | 10 | 5.291503 | 1 |
| 10 | 13 | 6.708204 | 1 |
| 12 | 18 | 8.124038 | 1 |
| 14 | 21 | 9.539392 | 1 |
| 18 | 33 | 12.369317 | 1 |

In every row the best result is exactly \(\sqrt E\), the elementary
Parseval/RMS bound. At \(L=1\), \(\lambda=E^{-1/2}\), the joint remainder
has Walsh degree three, and the dual distance four makes \(J=0\). Equation
(4.3) then reduces to \(\mu(a)\ge E^{-1/2}\).

The first level whose \(\lambda\) exceeds the *actual* normalized cap is
already vacuous because the coset denominator collapses. Two examples are

\[
\begin{array}{c|c|c|c|c|c}
n&L&\lambda&J/|\mathcal C|&T_a/|\mathcal C|&\lambda-J/T_a\\ \hline
14&3&0.242814&6.3352\cdot10^{-5}&1.0440\cdot10^{-5}&-5.82536\\
18&4&0.229312&5.9675\cdot10^{-6}&3.9254\cdot10^{-8}&-151.793
\end{array}
\]

This is finite numerical evidence, not a no-go theorem for \(k>0\) or for a
new \(S_n\)-adapted transition graph. It does rigorously identify the ratio
that any such construction must improve.

Reproduction:

    .venv/bin/python tmp/moving_cut_code_kernel_audit.py

The script passes py_compile; its SHA-256 after this audit is recorded in the
handoff.

## 8. Exact sufficient statement and falsification criterion

A moving-representation route would prove the conference constant, and hence
convergence, if it supplied kernels \((\kappa_n,F_n)\) satisfying the joint
Gram identity (3.1) and

\[
  \lambda_n=\frac{1-o(1)}{\sqrt n},
  \qquad
  \sup_{a\in\Omega_E}
  \frac{J_{\mathcal C_n^+}}{T_a}=o(n^{-1/2}).
  \tag{8.1}
\]

Indeed, (4.3) would give uniformly

\[
  Q(a)\ge\frac{1-o(1)}{\sqrt n}E
       =\left(\frac12-o(1)\right)n^{3/2}.
\]

The conference upper bound would then force
\(M_n/n^{3/2}\to1/2\). A weaker constant in (8.1) would give a rigorous
lower-bound improvement but would not by itself prove convergence.

This formulation is falsifiable without solving the original optimization:

1. Specify the transition graph and its projections independently of
   \(M_n\).
2. Compute the explicit numerator (5.2).
3. Test (5.3) on held-out exact minimizers, conference signings, and random
   nonconference signings.
4. Reject the family if \(J/T_a\) exceeds \(\lambda\) by a fixed
   \(n^{-1/2}\)-scale amount on a scalable family, or if \(T_a=0\).
5. Continue only if a representation-theoretic identity gives a uniform
   lower bound for (6.2) that uses the complete weighted Eulerian remainder.

The required new ingredient is therefore not another packing estimate. It
is an \(S_n\)-adapted lower bound for a moving kernel's coset convolution.
That statement is genuinely covering-specific and is the exact point at
which the August code theorem stops.
