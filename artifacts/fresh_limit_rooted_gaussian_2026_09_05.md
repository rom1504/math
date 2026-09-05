# Gaussian limit of the rooted response: contraction proof

Date: 2026-09-05. Status: independently audited analytic theorem.
The director reconstructed the argument; two independent agents checked all
four contraction cases, symmetrization, transfer, limiting quantifiers, and
primary theorem hypotheses. The exact rational numerical certificate was
also independently rerun. See `fresh_limit_rooted_independent_audit_2026_09_05.md`
and `fresh_rooted_literature_adversarial_audit_2026_09_05.md`.

Use the normalization and lemmas of
`fresh_limit_rooted_response_2026_09_05.md`: `m=n-1`, `B=A/sqrt(m)`, `Q=B²`,
and `q(A)=O(n^(3/2))`. In particular

\[
 Q_{jj}=1,\quad |B_{jk}|\le m^{-1/2},\quad
 \|Q\|_{\rm op}=O(\sqrt n),\quad
 \operatorname{Tr}Q^2=O(n^{3/2}).                       \tag{1}
\]

The stronger claim is that, for every fixed even smooth `h` with suitable
polynomial-growth derivatives and every coordinate `i`,

\[
 \sum_jB_{ij}S_jh((BS)_j)
       \ \Longrightarrow\ N(0,\mathbb Eh(Z)^2),       \tag{2}
\]

uniformly over coordinates and low-cap sequences. It is only a theorem about
this rooted field. It asserts nothing about the entire AMP response vector or
the nonrooted field `B f(BS)`.

## 1. The rooted Wick tensor

Fix an even integer `r≥2`. Write `b_j=B_j` for row vectors, `a_j=B_ij`, and

\[
 T=\sum_j a_j e_j\otimes b_j^{\otimes r},\qquad
 K=\operatorname{Sym}(T).
\]

For an isonormal standard Gaussian vector `Z`, the multiple integral of `K`
is exactly

\[
 I_{r+1}(K)=\sum_j a_j Z_jH_r(b_jZ),                    \tag{3}
\]

because `ej·bj=Bjj=0`. Its variance tends to `r!`, by the exact covariance
formula in the weak-decoupling note. We prove that every nontrivial
self-contraction of `K` tends to zero.

Expand `K⊗_ℓK`, `1≤ℓ≤r`, as the average of contractions between two permuted
copies of `T`. Up to a permutation of the remaining tensor coordinates, every
term has one of the following four forms. There are only finitely many terms,
depending on `r`, so bounding the norm of each form suffices.

### A. The two roots are contracted with each other

The root inner product forces the two summation indices to agree. The
remaining tensor is

\[
 \sum_j a_j^2
       b_j^{\otimes(r-\ell+1)}\otimes b_j^{\otimes(r-\ell+1)}.
\]

Its squared norm is at most

\[
 \sum_{j,k}a_j^2a_k^2|Q_{jk}|^{2(r-\ell+1)}
       \le{\operatorname{Tr}Q^2\over m^2}=o(1).        \tag{4}
\]

The exponent is at least two because `ℓ≤r`.

### B. Neither root is contracted

All `ℓ` contractions are branch-to-branch. The remaining roots force equal
summation indices when computing the squared norm, so that norm equals

\[
 \sum_{j,k}a_j^2a_k^2Q_{jk}^{2\ell}
       \le{\operatorname{Tr}Q^2\over m^2}=o(1).        \tag{5}
\]

### C. Exactly one root is contracted to an opposite branch

Choose the orientation in which root `j` is contracted to a branch at `k`.
There are `ℓ-1` further branch-to-branch contractions. The coefficient is
`a_j a_k B_kj Q_jk^(ℓ-1)`, and the remaining tensor is

\[
 b_j^{\otimes(r-\ell+1)}\otimes e_k
           \otimes b_k^{\otimes(r-\ell)}.
\]

For fixed `k` put `(v_k)_j=a_j B_kj Q_jk^(ℓ-1)` and
`P=Q^{∘(r-ℓ+1)}`. Its squared norm is exactly

\[
 \sum_k a_k^2 v_k^\top P v_k.
\]

Schur contraction gives `||P||op≤||Q||op`; also
`||v_k||²≤m^-1 Σ_j a_j²≤1/m`. Thus

\[
 \|\text{term}\|^2\le{\|Q\|_{\rm op}\over m}=o(1). \tag{6}
\]

### D. Each root is contracted to an opposite branch

This requires `ℓ≥2`. The coefficient matrix for the remaining two groups of
branches is

\[
 C_{jk}=a_ja_kB_{jk}^2Q_{jk}^{\ell-2},\qquad C_{jj}=0.
\]

Again put `P=Q^{∘(r-ℓ+1)}`. The squared norm is
`Tr(C P Cᵀ P)`, and hence is at most

\[
 \|P\|_{\rm op}^2\|C\|_F^2
 \le{\|Q\|_{\rm op}^2\over m^2}
           \sum_{j,k}a_j^2a_k^2
 \le{\|Q\|_{\rm op}^2\over m^2}=o(1).                 \tag{7}
\]

These cases exhaust the possibilities: a root is either left uncontracted,
paired to the other root, or paired to a branch. The estimates hold uniformly
in `i`. Consequently every self-contraction of `K` vanishes.

## 2. Applying the fixed-chaos Gaussian limit theorem

Theorem 1 of [Nualart and Peccati, *Central limit theorems for sequences of
multiple stochastic integrals*](https://arxiv.org/pdf/math/0503598) states that,
at fixed chaos degree with variance tending to one, vanishing of all
nontrivial self-contractions is equivalent to convergence to a standard
normal law. Applying it after normalization by `sqrt(r!)` proves that (3)
converges to `N(0,r!)`.

For a finite list of even `r`, the chaos degrees `r+1` are distinct, so their
covariances are exactly zero. Include also the degree-one Gaussian `b_iZ`.
The multivariate theorem gives joint convergence to independent normal
variables. A full primary proof and precise statement are available as
Theorem 7 of [Nualart and Ortiz-Latorre, *Central limit theorems for multiple
stochastic integrals and Malliavin calculus*](https://arxiv.org/pdf/math/0703240).
Its hypotheses are converging covariance and componentwise contraction
vanishing; the mixed contractions need not be imposed separately.

For every fixed polynomial truncation of `h`, the explicit distinct-index
comparison and low-influence Lindeberg argument in the weak-decoupling note
transfer this Gaussian conclusion to Rademacher inputs. The tail-transport
estimate there gives

\[
 \limsup_n\sup_i\mathbb E\left|
 \sum_jB_{ij}S_j(h-P)((BS)_j)\right|^2
                   \le\mathbb E(h(Z)-P(Z))^2.
\]

First take `n→∞`, then let the Hermite truncation degree tend to infinity.
This proves (2). The same uniform L² control gives convergence of expectations
against the linearly growing, Lipschitz function `z↦max(a,|z|)`.

## 3. Consequence for the original Boolean lower bound

For the smoothed two-orientation threshold functions `f,h`, the exact
conditional-dither convexity and the own-spin removal argument give

\[
 \liminf_n{\ell\over n}\ge
       L_\tau(t):=\mathbb E\max\{a_\tau(t),
                    \sqrt{\mathbb Eh_\tau(Z)^2}\,|Z|\}.
\]

Let `τ↓0` only after the dimension limit. Write

\[
 a(t)=2\phi(t),\qquad b(t)=2\Phi(t)-1.
\]

Then `Ehτ²→b(t)`, and the limiting field lower bound is

\[
 L(t)=a(t)\left[2\Phi\!\left({a(t)\over\sqrt{b(t)}}\right)-1\right]
       +2\sqrt{b(t)}\,\phi\!\left({a(t)\over\sqrt{b(t)}}\right).
                                                               \tag{8}
\]

The one-probe energy baseline is `a(t)b(t)`. For every fixed `p∈[0,1]`, the
exact partial-response inequality therefore implies

\[
 \liminf_n{M_n\over n^{3/2}}\ge
 { (1-p)^2 a(t)b(t)+p(1-p)L(t)\over1+p^2}.              \tag{9}
\]

The restriction to low-cap sequences entails no loss: if a competing
subsequence has unbounded normalized cap there is nothing to prove, and the
known universal upper construction provides a bounded-cap minimizing
sequence. Also `sqrt(n-1)/sqrt(n)→1` converts the present normalization.

Numerical optimization diagnostics (the last fixed choice is certified):

- At the old threshold `t≈0.876900985553`, `L(t)≈0.7718805447880719`, and
  optimizing `p` gives approximately `0.3396480554627936`.
- Joint optimization gives `t≈0.87482575558124`, `p≈0.06383197972828`, and
  approximately `0.3396496550505262`.
- The rational choices `t=7/8`, `p=8/125` give an exact rational interval
  certificate

\[
 0.339649621211865756397462117164309583060688819678602946169260
 \ \le\ \text{right side of (9)}\ \le
 0.339649621211865756397462117164309583060688819678602946169325.
\]

The reproducible source is `computations/fresh_limit_rooted_lower_certificate.py`,
with stdout saved as
`computations/results/fresh_limit_rooted_lower_certificate.json`. It uses only
exact `Fraction` arithmetic and outward rounding to a rational `10^-60` grid.
Machin's identity supplies `π` through alternating arctangent series;
exponential and Gaussian-integral series have explicit rational remainder
bounds. In particular, the lower endpoint is rigorously greater than `0.3396`.

## 4. A scoped optimality statement for the initial two-channel probe

This section optimizes only the displayed rooted-only lower-bound mechanism,
not arbitrary algorithms or the original Boolean optimum.

Suppose the initial probe means have `f` odd, `h` even, and `|f|+|h|≤1`.
Signs of the channels can be chosen so that `a=E[Z f(Z)]≥0` and `b=Eh(Z)≥0`.
Replacing `h` by `|h|` preserves `v=Eh(Z)²`, increases `b`, and preserves
the cube constraint. Thus it suffices to consider `0≤h≤1`. The same
argument yields the right side of (9) with
`L=E max(a,sqrt(v)|Z|)`. It is coordinatewise nondecreasing in `a,b,v` for
`p∈[0,1]`. At fixed `b`, one has `v≤b`. Moreover

\[
 a\le\mathbb E|Z|(1-h(Z)).
\]

Among measurable `h∈[0,1]` of mean `b`, the last expression is maximized by
placing `h=1` on the central interval of Gaussian measure `b`, and zero
outside. This is the elementary rearrangement inequality: moving any positive
mass of `h` from a larger `|Z|` to a smaller one cannot increase
`E|Z|h(Z)`. The same central indicator attains `v=b`, and choosing
`f(z)=sign(z)(1-h(z))` attains its maximal `a=2φ(t)`.

Thus central hard thresholds are optimal within this entire odd/even
two-channel certificate family. A different conversion of the updated energy
or additional rounds is outside this statement.

The larger conference-specific field value `0.805416...` is not claimed.
It retains extra noise in the nonrooted channel that the own-spin Jensen
step discards. Likewise, (9) does not address convergence or the upper bound
`1/2`.
