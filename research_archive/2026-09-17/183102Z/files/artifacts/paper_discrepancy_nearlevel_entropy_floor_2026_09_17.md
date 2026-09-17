# An unavoidable entropy floor for complete quadratic near-level codes

2026-09-17 campaign. Derived by the discrepancy agent; independently
reconstructed **PASS** by the root, BH/mechanism, and localization agents.
This result corrects the applicability of earlier full-code entropy
criteria: a finite entropy slope is impossible for every bounded,
positive-cap quadratic sequence, not merely unproved for minimizers.

## 1. Exact finite theorem

Let H(x)=Σ_(i<j) a_ij x_i x_j be any homogeneous quadratic form on
{−1,1}^n, n≥2. No sign-entry assumption is needed. Put
Q=max_x |H(x)|>0 and E(t)={y:Q−|H(y)|≤t}.
For every integer 0≤r≤n and every t>0,

$$
|E(t)|\ \ge\
\left(1-\frac{4r(n-r)Q}{n(n-1)t}\right)_+
          \binom nr.                                  \tag{1}
$$

The positive-part notation makes the estimate harmless when the
displayed factor is nonpositive.

Choose x and σ∈{−1,1} with σH(x)=Q, and flip a uniformly chosen
r-element coordinate set to obtain Y. For each i≠j,

$$
\mathbb E\frac{Y_iY_j}{x_ix_j}
=1-\frac{4r(n-r)}{n(n-1)}.
$$

Indeed, exactly one of i,j is flipped with probability
2r(n−r)/(n(n−1)). Consequently the nonnegative signed deficit
D=Q−σH(Y) satisfies

$$
\mathbb E D=\frac{4r(n-r)}{n(n-1)}Q.                    \tag{2}
$$

Markov's inequality implies that at least the fraction in (1) has
D≤t. Such a word also has Q−|H(Y)|≤t. The law of Y is uniform on a
Hamming sphere of size binom(n,r), proving (1).
No independent-coordinate approximation or concentration hypothesis
is used.

## 2. Exact asymptotic entropy envelope

Suppose Q_n/n^(3/2)→c∈(0,∞), and write

$$
E_n(\eta)=\{x:Q_n-|H_n(x)|\le\eta n^{3/2}\},\qquad
s_-(\eta)=\liminf_n\frac1n\log|E_n(\eta)|.
$$

For every 0<η<c,

$$
\boxed{\displaystyle
s_-(\eta)\ge
h\left(\frac{1-\sqrt{1-\eta/c}}2\right).}              \tag{3}
$$

For η≥c, s_−(η)=log2. To prove (3), choose any fixed
0<ρ<1/2 with 4cρ(1−ρ)<η and take r=⌊ρn⌋ in (1).
The fraction is bounded away from zero, so Stirling's formula gives
s_−(η)≥h(ρ). Increase ρ to the root in (3). For η≥c,
every ρ<1/2 satisfies the strict inequality (including η=c),
so s_−≥log2, the trivial upper bound.

In particular,

$$
\boxed{\displaystyle
\liminf_{\eta\downarrow0}
\frac{s_-(\eta)}{\eta\log(1/\eta)}
\ge\frac1{4c}.}                                       \tag{4}
$$

Thus for the usual upper entropy profile s(η)=limsup_n n^(−1)log|E_n(η)|,

$$
\lim_{\eta\downarrow0}\frac{s(\eta)}{\eta}=+\infty.       \tag{5}
$$

If the cap sequence does not converge but has finite positive limsup
C, replace c by C: the same lower bound follows from the eventual
upper cap bound. The result is not restricted to full-sign matrices.

## 3. Consequences and the corrected scope

Any theorem assuming a finite full-near-level entropy slope
limsup_(η↓0) s(η)/η has **no examples** among bounded positive-cap
quadratic sequences. Its finite bridge estimates may remain correct,
but its advertised asymptotic hypothesis is universally impossible.

Consequently:

- The direct full-code entropy-slope augmentation criterion is a
  formally valid implication with an impossible hypothesis.
- If a factorization-weighted entropy criterion uses e(η)s(η)/η
  and e(η)→e₀>0, its finite-slope regime is also impossible.
- The weighted case e₀=0 is not excluded by this argument.
- A low-response **center code plus a Hamming residual cover** is
  not excluded. In particular, cover costs ρh(ρ)/η may vanish when
  ρ=O(η), despite h(ρ)/η diverging.
- None of these statements proves convergence or changes the known
  unrestricted quadratic-cap constants.

Combining (4) with the flat-diagonal response-to-entropy theorem gives
an additional necessary response scale. If a full near-level code in
the regime k²=o(n) has uniform matrix response at most μ(η)→0, then

$$
\liminf_{\eta\downarrow0}\frac{\mu(\eta)^2}{\eta}
\ge\frac1{4c}.                                        \tag{6}
$$

Here the entropy upper bound is
(1+o(1))μ²log(1/μ²). Monotonicity of ulog(1/u) near zero
gives (6). The director subsequently proved a stronger, aspect-free
mean-response rigidity inequality by convex concentration and biased
flips; that is a separate result and does not enter (1)–(5).

## 4. Reproducibility

The finite identity and counting lower bound are tested exhaustively
for all full-sign matrices through order five in

computations/paper_discrepancy_2026_09_17_nearlevel_entropy_floor.py

The test enumerates every spin word, every flip radius, and every
positive integer deficit threshold through 2Q. Its arithmetic is
integral/rational; no floating-point conclusion is used.

