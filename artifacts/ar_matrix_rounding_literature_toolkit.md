# Matrix sign-rounding literature toolkit

Audit date: 2026-08-16. This note records theorem-level facts relevant to
rounding a symmetric hollow matrix \(W=(w_{ij})\in[-1,1]^{n\times n}\) to a
symmetric hollow sign matrix \(A\in\{\pm1\}^{n\times n}\). It is a literature
audit, not a proposed recovery argument.

## 1. Verdict at the required little-o scale

Put \(E=A-W\), and define

\[
 q(E):=\max_{x\in\{\pm1\}^n}|x^{\mathsf T}Ex|,
 \qquad
 B(E):=\|E\|_{\infty\to1}
      =\max_{x,y\in\{\pm1\}^n}|x^{\mathsf T}Ey|.
\]

The same-spin defect written as an unordered-edge sum is \(q(E)/2\).

There is **no theorem, and in fact no possible theorem**, that rounds every
such \(W\) with either

\[
 \|A-W\|_{\rm op}=o(\sqrt n)
 \quad\text{or}\quad
 B(A-W)=o(n^{3/2})
\]

uniformly over \(W\). The single legal input \(W=0\) rules out both
conclusions for every \(A\). Thus results stated as \(O(\sqrt n)\),
\((2+o(1))\sqrt n\), or \(O(n^{3/2})\) are at the natural scale but do not meet
the required little-\(o\) threshold.

The sharpest direct deterministic supported-outcome theorem is
Wang--Lau--Zhou (STOC 2026): for the natural biased edgewise rounding it
selects one supported outcome with operator norm at most the free-probability
leading term \(2\sqrt{v_{\max}}\), plus lower-order covariance terms. Here

\[
 v_{\max}(W):=\max_i\sum_{j\ne i}(1-w_{ij}^2).
\]

Consequently it reaches \(o(\sqrt n)\) only in a non-worst-case regime such as
\(v_{\max}=o(n)\), together with the explicit logarithmic remainder condition
given in Section 3. For existence without deterministic selection,
Bandeira--van Handel's bounded independent-entry estimate already gives

\[
 \mathbb E\|A-W\|_{\rm op}
 \le C\bigl(\sqrt{v_{\max}}+\sqrt{\log n}\bigr)
\]

for the same biased rounding, and hence \(v_{\max}=o(n)\) alone is sufficient.
At \(W=0\), \(v_{\max}=n-1\), and both bounds remain of order \(\sqrt n\), as
they must.

For scalar same-spin or bilinear discrepancy, ordinary biased independent
rounding plus scalar Bernstein gives an outcome with error

\[
 O\!\left(\sqrt{nV}+n\right),
 \qquad
 V(W):=\sum_{i<j}(1-w_{ij}^2).
\]

This is \(o(n^{3/2})\) when \(V=o(n^2)\), but only \(O(n^{3/2})\) in the dense
worst case. Classical shifted/linear discrepancy gives the analogous bound
\(O(\sqrt{rn})\) if only \(r\) edge coordinates are genuinely fractional.

## 2. Exact translations and unavoidable lower bounds

### 2.1 Same-spin and bilinear discrepancy are equivalent up to factor two

For every real symmetric hollow \(E\),

\[
 \boxed{\quad \tfrac12 B(E)\le q(E)\le B(E).\quad}
\]

The upper inequality follows by taking \(x=y\). Conversely, for sign vectors
\(x,y\), let \(u=(x+y)/2\) and \(v=(x-y)/2\). Symmetry gives

\[
 x^{\mathsf T}Ey=u^{\mathsf T}Eu-v^{\mathsf T}Ev.
\]

Each of \(u,v\) is a ternary vector. Complete its zero coordinates by
independent uniform signs. Because \(E\) is hollow, its ternary quadratic
form is the expectation of the completed sign quadratic forms, and hence has
absolute value at most \(q(E)\). Therefore \(B(E)\le2q(E)\). In particular,
the two scalar targets have exactly the same little-\(o(n^{3/2})\) threshold.

### 2.2 Entrywise distance already forces spectral error

Let \(d_{ij}=1-|w_{ij}|\). Every sign choice obeys
\(|a_{ij}-w_{ij}|\ge d_{ij}\), so

\[
 \boxed{\quad
 \|A-W\|_{\rm op}
 \ge {\|A-W\|_F\over\sqrt n}
 \ge\left({2\over n}\sum_{i<j}d_{ij}^2\right)^{1/2}.
 \quad}
\]

Thus a necessary condition for spectral \(o(\sqrt n)\) is
\(\sum_{i<j}d_{ij}^2=o(n^2)\).

There is also a direct scalar lower bound. The sharp \(p=1\) Khintchine
constant \(1/\sqrt2\) gives, for any fixed \(A\),

\[
 \begin{aligned}
 B(E)&=\max_y\|Ey\|_1
      \ge \mathbb E_y\|Ey\|_1\\
 &\ge {1\over\sqrt2}\sum_i
       \left(\sum_{j\ne i}e_{ij}^2\right)^{1/2}
 \ge {1\over\sqrt2}\sum_i
       \left(\sum_{j\ne i}d_{ij}^2\right)^{1/2}.
 \end{aligned}
\]

Together with \(q(E)\ge B(E)/2\), this is a deterministic obstruction for
both scalar norms. The sharp Khintchine constant is due to Haagerup,
[The best constants in the Khintchine inequality](https://doi.org/10.4064/sm-70-3-231-283),
*Studia Mathematica* 70 (1981), 231--283.

### 2.3 The decisive input \(W=0\)

If \(W=0\), every admissible \(A\) has \(n(n-1)\) unit off-diagonal entries.
Consequently

\[
 \boxed{\|A\|_{\rm op}\ge\sqrt{n-1}.}
\]

The scalar obstruction has an exact finite-\(n\) form. For a uniform sign
vector \(y\), every row of \(Ay\) is distributed as a sum \(S_{n-1}\) of
\(n-1\) independent Rademacher variables. Hence

\[
 B(A)\ge n\mu_{n-1},\qquad
 q(A)\ge {n\mu_{n-1}\over2},\qquad
 \max_x\left|\sum_{i<j}a_{ij}x_ix_j\right|
 \ge {n\mu_{n-1}\over4},
\]

where \(\mu_m=\mathbb E|S_m|\), with

\[
 \mu_{2k}={2k\binom{2k}{k}\over2^{2k}},\qquad
 \mu_{2k+1}={(2k+1)\binom{2k}{k}\over2^{2k}},\qquad
 \mu_m\sim\sqrt{2m/\pi}.
\]

Thus \(B(A)\) is always at least
\((\sqrt{2/\pi}+o(1))n^{3/2}\), \(q(A)\) at least
\((1/\sqrt{2\pi}+o(1))n^{3/2}\), and the unordered same-spin defect at least
\((1/(2\sqrt{2\pi})+o(1))n^{3/2}\).

For context, Erdős--Spencer's primary result
[Imbalances in k-colorations](https://doi.org/10.1002/net.3230010407),
*Networks* 1 (1971/72), 379--385, proves the related induced-subgraph
\(k\)-color discrepancy scale \(\Theta(n^{(k+1)/2})\); at \(k=2\) this is also
\(\Theta(n^{3/2})\). The elementary norm argument above is the exact one for
the present matrix target.

## 3. Direct entrywise rounding: the sharp 2026 theorem

For every unordered edge \(e=\{i,j\}\), set
\(S_e=e_ie_j^{\mathsf T}+e_je_i^{\mathsf T}\), and independently choose

\[
 A_e=1\text{ with probability }{1+w_e\over2},\qquad
 A_e=-1\text{ with probability }{1-w_e\over2}.
\]

Then

\[
 E=A-W=\sum_e Z_e,\qquad Z_e=(A_e-w_e)S_e.
\]

The summands are centered, self-adjoint, and have two-point support. The exact
hypotheses of Wang--Lau--Zhou's simplified Theorems 1.4--1.5 are: \(A_0\) is
deterministic Hermitian; \(Z_1,\ldots,Z_N\) are independent self-adjoint
\(d\times d\) random matrices; \(\mathbb EZ_e=0\);
\(\|Z_e\|\le\rho\) almost surely; and every support has polynomially bounded
size (with the supports and probabilities explicitly enumerable for the
polynomial-time conclusion). The laws need not be identical, symmetric, or
unbiased Bernoulli laws. Thus the two-point variables above, with arbitrary
edge-dependent probabilities \((1+w_e)/2\), are covered exactly; an edge with
\(w_e=\pm1\) is simply a degenerate one-point centered summand.

In the notation of that theorem,

\[
 \sigma^2:=\left\|\sum_e\mathbb E Z_e^2\right\|=v_{\max},
 \qquad
 \nu^2:=\|\operatorname{Cov}(\operatorname{vec}E)\|
       =2\max_e(1-w_e^2)\le2,
 \qquad
 \rho:=\max_{e,z\in\operatorname{supp}Z_e}\|z\|\le2.
\]

The covariance identity uses the Frobenius-orthogonality of the \(S_e\)'s.

Theorem 1.5 and its corollary in Robert Wang, Lap Chi Lau, and Hong Zhou,
[Derandomizing Matrix Concentration Inequalities from Free Probability](https://arxiv.org/abs/2601.08111),
STOC 2026, give a deterministic supported outcome \(Z'\) satisfying

\[
 \begin{aligned}
 \|Z'\|\le \|X_{\rm free}\|
 +O(L)\big(&\sigma^{1/2}\nu^{1/2}\log^{3/4}d
 +\sigma^{2/3}\rho^{1/3}\log^{2/3}d
 +\rho\log d\big),\\
 L&=1+\log(N\rho^2/\sigma^2),
 \end{aligned}
\]

for \(X=\sum_{e=1}^N Z_e\) in dimension \(d\). Pisier's free bound gives
\(\|X_{\rm free}\|\le2\sigma\) in this centered case. Specializing
\(N=\binom n2\), \(d=n\), and absorbing the absolute bounds on \(\nu,\rho\),
one obtains, when \(v_{\max}>0\),

\[
 \boxed{
 \begin{aligned}
 \|A-W\|_{\rm op}\le 2\sqrt{v_{\max}}
 +O(L)\big(&v_{\max}^{1/4}\log^{3/4}n\\
            &+v_{\max}^{1/3}\log^{2/3}n+\log n\big).
 \end{aligned}}
\]

The displayed formula, rather than a hidden-tilde abbreviation, records the
conditioning logarithm. If \(v_{\max}=0\), then \(W\) is already a sign
matrix and \(A=W\). A sufficient little-\(o\) hypothesis is

\[
 v_{\max}=o(n)
 \quad\text{and}\quad
 L\big(v_{\max}^{1/4}\log^{3/4}n
       +v_{\max}^{1/3}\log^{2/3}n+\log n\big)=o(\sqrt n).
\]

Under ordinary polynomial conditioning, \(L=O(\log n)\), so the second
condition is lower order throughout the usual dense asymptotic regimes. For
arbitrary \(W\), the theorem says only \(2\sqrt n+\widetilde O(n^{1/3})\).

The same paper's graph-signing corollary says that a \(k\)-regular graph, for
\(k\) at least polylogarithmic in its number of vertices, has a deterministically
found signing whose signed adjacency norm is

\[
 2\sqrt{k}\bigl(1+\widetilde O(k^{-1/6})\bigr).
\]

For \(K_n\), this is an every-order \(O(\sqrt n)\) statement, not an
\(o(\sqrt n)\) statement.

The probabilistic inputs behind this result are the sharp free-probability
matrix concentration theory of Bandeira--Boedihardjo--van Handel,
*Inventiones Mathematicae* 234 (2023), 419--487, and
Bandeira--Cipolloni--Schröder--van Handel,
[Matrix Concentration Inequalities and Free Probability II](https://arxiv.org/abs/2406.11453)
(2024). These locate the correct leading natural scale; they do not remove
that leading term in the worst case.

## 4. Scalar biased rounding and exact sufficient variance condition

The same independent biased edge rounding yields a concise existence theorem
for the scalar norms. For fixed \(x\), write

\[
 H_E(x)=\sum_{i<j}(A_{ij}-w_{ij})x_ix_j.
\]

Its summands have total variance \(V=\sum_{i<j}(1-w_{ij}^2)\) and absolute
value at most \(2\). Scalar Bernstein and a union bound over \(2^n\) spins
give the explicit tail

\[
 \Pr\!\left\{\max_x|H_E(x)|\ge t\right\}
 \le 2^{n+1}\exp\!\left[-{t^2\over2(V+2t/3)}\right].
\]

For the bilinear form,

\[
 x^{\mathsf T}Ey
 =\sum_{i<j}(A_{ij}-w_{ij})(x_iy_j+x_jy_i).
\]

Each summand has absolute value at most \(4\) and the total variance is at
most \(4V\). A union bound over \(4^n\) pairs gives

\[
 \Pr\!\left\{B(E)\ge t\right\}
 \le 2\,4^n\exp\!\left[-{t^2\over2(4V+4t/3)}\right].
\]

In particular, for an absolute constant \(C\), some supported sign outcome
obeys simultaneously

\[
 \max_x|H_E(x)|\le C(\sqrt{nV}+n),
 \qquad
 B(E)\le C(\sqrt{nV}+n).
\]

Therefore \(V=o(n^2)\) is a clean sufficient condition for both requested
scalar errors to be \(o(n^{3/2})\). For \(W=0\), \(V=\binom n2\), and the
conclusion stops at \(O(n^{3/2})\), matching the lower-bound scale.

At operator level, the classical self-adjoint matrix Bernstein inequality
(Tropp, [User-Friendly Tail Bounds for Sums of Random Matrices](https://arxiv.org/abs/1004.4389),
*Foundations of Computational Mathematics* 12 (2012), 389--434) specializes to

\[
 \Pr\{\|E\|_{\rm op}\ge t\}
 \le2n\exp\!\left[-{t^2\over2v_{\max}+(4/3)t}\right].
\]

Its generic dense guarantee is \(O(\sqrt{n\log n})\). The sharper
independent-entry theorem of Bandeira--van Handel,
[Sharp nonasymptotic bounds on the norm of random matrices with independent entries](https://arxiv.org/abs/1408.6185),
*Annals of Probability* 44 (2016), 2479--2506, has the following exact form
useful here. Their Corollary 3.6 assumes that \(X\) is symmetric,
\(X_{ij}=\xi_{ij}b_{ij}\), and that
\(\{\xi_{ij}:i\ge j\}\) are independent, symmetric, unit-variance random
variables. With

\[
 \sigma=\max_i\left(\sum_jb_{ij}^2\right)^{1/2},
\]

it states, for every \(\alpha\ge3\),

\[
 \mathbb E\|X\|
 \le e^{2/\alpha}\left[
       2\sigma+
       14\alpha\max_{ij}\|X_{ij}\|_{2\lceil\alpha\log n\rceil}
       \sqrt{\log n}\right].
\]

The centered biased Bernoulli residual \(\zeta_{ij}=A_{ij}-w_{ij}\) is not
symmetric in law unless \(w_{ij}=0\), but it is still covered by the standard
symmetrization step in the proof. If \(E'\) is an independent copy, convexity
gives \(\mathbb E\|E\|\le\mathbb E\|E-E'\|\). The upper-triangular entries of
\(E-E'\) are independent and symmetric,

\[
 \mathbb E(E_{ij}-E'_{ij})^2=2(1-w_{ij}^2),
 \qquad |E_{ij}-E'_{ij}|\le2.
\]

Taking \(\alpha=3\) in the displayed corollary therefore gives the fully
explicit bound

\[
 \boxed{\quad
 \mathbb E\|A-W\|_{\rm op}
 \le e^{2/3}\left(2\sqrt{2v_{\max}}+84\sqrt{\log n}\right).
 \quad}
\]

No lower bound on the biases, no identical-distribution assumption, and no
uniform subgaussian assumption on the variance-normalized entries is needed.
In particular, some supported sign outcome satisfies the same bound. This
proves operator \(o(\sqrt n)\) whenever \(v_{\max}=o(n)\), with no numerical
conditioning logarithm; unlike Wang--Lau--Zhou, it is an expectation/existence
argument and does not retain the sharp leading constant two or supply the
same deterministic outcome-selection algorithm.

The inhomogeneous random graph result of Lu--Peng,
[Spectra of edge-independent random graphs](https://doi.org/10.37236/3576),
*Electronic Journal of Combinatorics* 20(4) (2013), P27, recover sharp
\((2+o(1))\sqrt\Delta\)-type behavior in dense regimes. That is still the
natural \(O(\sqrt n)\) scale rather than the required little-\(o\) scale.

## 5. Shifted discrepancy and partial coloring

Let the columns be the \(r\) genuinely fractional edges and let a row encode
either \(x_ix_j\) for a spin \(x\), or \(x_iy_j+x_jy_i\) for a sign pair
\((x,y)\). There are respectively \(m=2^n\) or \(m=4^n\) rows.

Spencer's general discrepancy bound, from
[Six standard deviations suffice](https://www.ams.org/journals/tran/1985-289-02/S0002-9947-1985-0784009-0/),
*Transactions of the AMS* 289 (1985), 679--706, gives for every restriction
to \(s\) columns

\[
 \operatorname{disc}=O\!\left(\sqrt{s\log(2m/s)}\right).
\]

Lovász--Spencer--Vesztergombi,
[Discrepancy of set-systems and matrices](https://doi.org/10.1016/S0195-6698(86)80041-5),
*European Journal of Combinatorics* 7 (1986), 151--160, compare linear
(shifted) discrepancy with hereditary discrepancy within an absolute factor
(at most two under the standard sign normalizations). Hence an arbitrary
fractional vector on those \(r\) edges can be rounded with

\[
 O\!\left(\sqrt{rn}\right)
\]

error in all same-spin rows, and likewise in all bilinear rows. This reaches
little-\(o(n^{3/2})\) when \(r=o(n^2)\), but only the natural
\(O(n^{3/2})\) scale when \(r=\Theta(n^2)\).

Modern partial-coloring results do not change that worst-case exponent for
this encoding. Dadush--Jiang--Reis,
[A new framework for matrix discrepancy](https://arxiv.org/abs/2111.03171)
(2022), prove matrix partial-coloring/full-coloring estimates governed by
rank, blocks, and Schatten norms. With \(N\asymp n^2\) edge variables, their
general matrix-Spencer scale is \(O(\sqrt N)=O(n)\), not \(o(\sqrt n)\).
Lau--Wang--Zhou,
[Spectral Sparsification by Deterministic Discrepancy Walk](https://doi.org/10.1137/1.9781611978315.24),
SOSA 2025, give a deterministic proof of the Reis--Rothvoss matrix partial
coloring theorem and sparsification applications, but no stated
little-\(o(\sqrt n)\) completion theorem for linked signed adjacency entries.

## 6. Matrix Spencer theorems: why their parameter is too large here

Treating each edge matrix \(S_e\) as one arbitrary contraction invokes a
number of variables \(N=\binom n2\). The matrix Spencer scale \(O(\sqrt N)\)
is therefore \(O(n)\).

* Bansal--Jiang--Meka,
  [Resolving Matrix Spencer Conjecture Up to Poly-logarithmic Rank](https://doi.org/10.1137/23M1592201),
  *SIAM Journal on Computing* (2024), prove the conjectured \(O(\sqrt N)\)
  bound when the summands have rank at most \(N/\log^3N\). Edge matrices have
  rank two and satisfy the rank hypothesis, but the resulting \(O(n)\) bound
  is much weaker than \(o(\sqrt n)\).

* Akbas--Sra,
  [An Algebraic Matrix Spencer Theorem](https://arxiv.org/abs/2606.16005)
  (2026), prove an \(O(\sqrt N)\) signing theorem for contractions in a
  finite-dimensional \(C^*\)-algebra whose algebraic dimension is comparable
  to \(N\). The edge matrices generate \(M_n\), of dimension \(n^2\asymp N\),
  so again the direct scale is \(O(n)\).

The latter paper also proves a useful no-go statement: there are \(N\)
diagonal \(N\times N\) contractions whose signing discrepancy is
\(\Theta(\sqrt{\log N})\), although

\[
 \left\|\sum_i A_i^2\right\|^{1/2}=\Theta(1).
\]

Thus a blanket matrix-Spencer bound of the form
\(C\|\sum_iA_i^2\|^{1/2}\) is false. The free norm and covariance structure
in Wang--Lau--Zhou cannot simply be replaced by the row-variance parameter in
an arbitrary deterministic signing theorem.

A genuine shifted, variance-sensitive theorem exists for rank-one positive
semidefinite summands. Kyng--Luh--Song,
[Four Deviations Suffice for Rank 1 Matrices](https://arxiv.org/abs/1901.06731),
*Advances in Mathematics* 375 (2020), 107366, show that for independent
finite-support scalars \(\xi_i\), one can select supported values \(\varepsilon_i\)
so that

\[
 \left\|\sum_i(\mathbb E\xi_i-\varepsilon_i)u_iu_i^*\right\|
 \le4\left\|\sum_i\operatorname{Var}(\xi_i)(u_iu_i^*)^2\right\|^{1/2}.
\]

An edge matrix is indefinite rank two:

\[
 S_{ij}=\tfrac12(e_i+e_j)(e_i+e_j)^{\mathsf T}
       -\tfrac12(e_i-e_j)(e_i-e_j)^{\mathsf T}.
\]

The two rank-one coefficients must use the same edge choice with opposite
signs. The rank-one theorem does not preserve this linkage, so it is not a
rounding theorem for \(W\mapsto A\). The same linkage issue separates the
MSS/Kadison--Singer rank-one machinery from arbitrary weighted adjacency
rounding.

## 7. Pipage rounding and interlacing/signing results

### Pipage rounding

Harvey--Olver,
[Pipage Rounding, Pessimistic Estimators and Matrix Concentration](https://arxiv.org/abs/1307.2274),
SODA 2014, prove randomized and deterministic pipage rounding over matroid
base polytopes. For \(0\preceq M_i\preceq RI\), it preserves the matrix
Chernoff upper tail

\[
 \Pr\!\left\{\lambda_{\max}\!\left(\sum_i X_iM_i\right)
      \ge(1+\delta)\mu\right\}
 \le d\left({e^\delta\over(1+\delta)^{1+\delta}}\right)^{\mu/R}.
\]

This is valuable when hard matroid constraints accompany PSD selection.
Here the residual edge contribution is indefinite. PSD decompositions add
degree-diagonal terms that cancel only if the relevant degree constraints are
also preserved. Even in a legal PSD encoding, the theorem supplies a
Chernoff/natural-scale bound, not a universal little-\(o\) improvement.

### Interlacing and graph signings

Marcus--Spielman--Srivastava,
[Interlacing Families I: Bipartite Ramanujan Graphs of All Degrees](https://arxiv.org/abs/1304.4132),
*Annals of Mathematics* 182 (2015), prove that every graph has a signing for
which the **largest** eigenvalue of the signed adjacency matrix is at most the
spectral radius of its universal cover. For a bipartite graph the spectrum is
symmetric, so this is a two-sided operator bound; for a \(d\)-regular
bipartite graph it is at most \(2\sqrt{d-1}\). For nonbipartite \(K_n\), the
original theorem is one-sided and cannot be quoted directly as an operator
norm bound.

Xu--Zhang,
[An Improved Upper Bound for the Bilu--Linial Conjecture via Interlacing Families](https://arxiv.org/abs/2606.28797)
(2026), prove the all-graph, two-sided statement

\[
 \|A_s\|_{\rm op}\le2\sqrt{3(d-1)}
\]

for some signing of every graph of maximum degree \(d\). Wang--Lau--Zhou's
2026 result improves the dense regular bound to asymptotic leading constant
two. Both give \(O(\sqrt n)\) for \(K_n\), and neither gives
\(o(\sqrt n)\).

For comparison, the \(W=0\) spectral lower bound is attained on infinite
subsequences. A symmetric conference matrix \(C\) has diagonal zero,
off-diagonal signs, and \(C^2=(n-1)I\), so
\(\|C\|_{\rm op}=\sqrt{n-1}\). Paley's construction gives such matrices of
order \(q+1\) for prime powers \(q\equiv1\pmod4\); see
[On Orthogonal Matrices](https://doi.org/10.1002/sapm1933121311) (1933).

There is also an elementary every-order upper bound at \(W=0\). Let \(m\) be
the least power of two at least \(n\), take an \(n\times n\) principal
submatrix \(B\) of the symmetric Sylvester Hadamard matrix \(H_m\), and erase
its diagonal. The result \(A\) is symmetric, hollow, and signed, with

\[
 \|A\|_{\rm op}\le\|B\|_{\rm op}+1\le\sqrt m+1
 <\sqrt{2n}+1.
\]

Thus the correct all-order spectral scale at the maximally fractional input
is exactly \(\Theta(\sqrt n)\), while the scalar scale is exactly
\(\Theta(n^{3/2})\).

## 8. Applicability ledger

| Primary theorem | What it actually rounds/signs | Direct scale here | Meets the required little-\(o\)? |
|---|---|---:|---|
| Wang--Lau--Zhou, STOC 2026 | Exact supported outcomes of centered self-adjoint finite-support summands | \(2\sqrt{v_{\max}}\) plus explicit lower-order terms | Yes under the displayed sublinear-variance/remainder condition; no uniformly |
| Bandeira--van Handel, 2016, Cor. 3.6 + symmetrization | Independent centered bounded upper-triangular entries; arbitrary non-identical biases are allowed | \(C(\sqrt{v_{\max}}+\sqrt{\log n})\) in expectation | Yes if \(v_{\max}=o(n)\); no uniformly |
| Scalar Bernstein + union bound | Biased edgewise sign rounding | \(O(\sqrt{nV}+n)\) for same-spin and \(\infty\to1\) | Yes if \(V=o(n^2)\); no uniformly |
| Spencer + LSV linear discrepancy | Arbitrary fractional vector in exponentially many scalar constraints | \(O(\sqrt{rn})\) | Yes if \(r=o(n^2)\); no for \(r\asymp n^2\) |
| Bansal--Jiang--Meka / Akbas--Sra matrix Spencer | Signs on \(N\asymp n^2\) arbitrary matrix contractions | \(O(\sqrt N)=O(n)\) | No |
| Kyng--Luh--Song | Shifted finite-support rank-one PSD sums | \(4\sigma\) | Not directly applicable: linked indefinite edge pair |
| Harvey--Olver pipage | PSD selection sums under matroid constraints | Matrix-Chernoff scale | Not directly applicable to indefinite edge residuals; no universal little-\(o\) |
| MSS interlacing | Graph signing; two-sided automatically for bipartite graphs | \(2\sqrt{d-1}\) | No at \(d\asymp n\); original all-graph statement is one-sided |
| Xu--Zhang 2026 | Two-sided signing of every max-degree-\(d\) graph | \(2\sqrt{3(d-1)}\) | No at \(d\asymp n\) |

The strict conclusion of this audit is therefore twofold: arbitrary-\(W\)
little-\(o\) rounding is mathematically impossible, while the literature does
contain exact every-order rounding theorems at the required scale once the
input supplies a genuinely subcritical variance or fractional-support
parameter. Merely citing a natural-scale \(O\)-bound does not cross the
threshold.
