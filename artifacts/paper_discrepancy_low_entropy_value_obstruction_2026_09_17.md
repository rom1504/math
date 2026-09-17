# Low entropy alone does not make physical bridge responses cheap

2026-09-17 campaign. This is a finite probabilistic code construction,
NOT a construction of actual quadratic near-ground codes. The final
section quantifies precisely why the random example cannot be such a
near-ground code at positive normalized cap.

## 1. A polynomial code can force the independent-sign response constant

Put d=binom(n,2), and let μ_n=E|Σ₁ⁿ ε_i|. Sample M independent
uniform Boolean words X¹,…,Xᴹ and let π be their empirical law.
For every u>0, with failure probability at most exp(−u),

$$
\inf_{h\in\{\pm1\}^n}\mathbb E_{\pi}|h\cdot X|
\ge\mu_n-\sqrt{\frac{2n(n\log2+u)}M}.                \tag{1}
$$

Indeed the centered absolute overlap for a fixed h is n-subGaussian
by bounded differences. The empirical average is n/M-subGaussian.
Use its lower tail and union over all 2ⁿ physical sign queries.
Repeated sampled words are harmless: π is an empirical probability
law, and its support C has size at most M. Adding all antipodes to C
does not alter the conclusion or more than double its size.

By averaging (1), EVERY probability law ν on physical sign columns,
without any isotropy requirement, therefore has

$$
\max_{x\in C}\mathbb E_{h\sim\nu}|h\cdot x|
\ge\mu_n-\sqrt{\frac{2n(n\log2+u)}M}.                \tag{2}
$$

For example M=n³ and u=n give a code with log|C|≤3log n and
lower response `(sqrt(2/pi)−o(1)) sqrt(n)` against ALL sign laws.
The asymptotic value of μ_n follows either from the elementary central
binomial formula or the ordinary one-dimensional CLT with uniform
integrability. A direct calculation gives μ_n= n binom(n−1,floor((n−1)/2))
/2^(n−1).

Consequently small entropy by itself cannot guarantee a physical-column
mean response constant below sqrt(2/pi). For comparison, the reported
upper cap constant below 0.494 has linear dilution budget below
`(3/2)·0.494=0.741`, whereas sqrt(2/pi)>0.797. This comparison is
only a warning about entropy-only value arguments, not a recurrence
obstruction for actual minimizers.

## 2. The same random code is excluded by actual quadratic energy geometry

For x∈{±1}ⁿ write z(x)=(x_i x_j)_(i<j)∈{±1}ᵈ. For a uniform
word X,

$$
\mathbb E z(X)z(X)^T=I_d,\qquad\|z(X)\|^2=d.
$$

The elementary rank-one matrix Chernoff bound gives

$$
\Pr\!\left\{\frac1M\sum_{a=1}^M
             z(X^a)z(X^a)^T\not\preceq2I_d\right\}
\le d\exp[-M/(3d)].                                 \tag{3}
$$

For instance M=n³ makes this failure exponentially small. Hence, for
all sufficiently large n, a single polynomial-size code and empirical
law simultaneously satisfy (1) with u=n and

$$
\mathbb E_\pi H_A(X)^2\le2\sum_{i<j}A_{ij}^2          \tag{4}
$$

for EVERY real hollow symmetric matrix A, selected after the code.
For a full signing the right side is 2d=n(n−1).

If the entire support C were contained in the absolute near-level
set of A at deficit T<Q(A), then (4) would force

$$
\boxed{\quad Q(A)-T\le\sqrt{n(n-1)}.\quad}           \tag{5}
$$

Thus when Q(A)/n^(3/2)→c>0, this random code cannot be a full
near-ground support at any T=o(n^(3/2)). Its bad physical response
does NOT provide a counterexample inside the actual regularized
near-minimizer class.

Conversely, any probability law π supported on |H_A|≥Q(A)−T
for a full signing must satisfy the useful quadratic-feature condition

$$
\left\|\mathbb E_\pi z(X)z(X)^T\right\|_{\rm op}
\ge\frac{[Q(A)-T]_+^2}{d}.                           \tag{6}
$$

At positive normalized cap and a vanishing relative deficit, this
lower bound is of order n, not order 1. Equation (6) is not itself a
cheap-column theorem. It isolates a structural piece of actual energy
geometry that an entropy-only argument throws away, and that a future
Gaussian-parent value bound would have to exploit.

The matrix concentration ingredient in (3) is the same finite Chernoff
lemma reconstructed in the localization track's critical-realization
artifact. Its role here is only to demonstrate the simultaneous scope
separation; no conditional law of a selected optimizer is assumed.

## 3. The exact first-order value obligation for a small old code

Let the old word be restricted to a fixed finite code C, keep the
complete child B of order q, all new spins, and both polarities. For
each fixed old x suppose the bridge fields v_j(x) are independent
and symmetric over j. Define

$$
m(x)=\sum_{j=1}^q\mathbb E|v_j(x)|,\qquad
L_C=\max_{x\in C}\{|H_A(x)|+m(x)\}.
$$

Assume the centered variable Σ_j|v_j(x)| is V_*-subGaussian uniformly
over x. Then the expected restricted absolute cap obeys the sandwich

$$
\boxed{\quad
L_C\le\mathbb E Q_C(\mathrm{parent})
\le L_C+Q(B)+\sqrt{2V_*\log|C|}.\quad}              \tag{7}
$$

For the lower bound fix x, choose its favorable old-energy polarity s,
and take y_j=s sign(v_j(x)), breaking zero ties by independent fair
signs. These y_j are independent and fair. The expected child quadratic
energy is zero, while the bridge contributes exactly m(x) in signed
expectation. The true maximum is at least this randomized choice.
Then maximize the deterministic lower expectation over x.

For the upper bound, deterministically bound the child term by Q(B)
and maximize each bridge field in absolute value. The maximum over
the fixed old code of the centered absolute-field sum costs at most
the displayed subGaussian maximum bound. No independence across
different old words is used, and |C|=1 is covered with zero overhead.

In the exact-sign frame and its Gaussian comparison from the director's
codewise universality theorem, write

$$
V_j(x)=\sum_t(h_j\cdot x_{:,t})^2+\ell.
$$

For the sign law,
`m_R(x)=Σ_j E|Σ_t(h_j·x_col) ε_(j,t)+Σ_left x_i ξ_(j,i)|`.
For the Gaussian law,

$$
m_G(x)=\sqrt{\frac2\pi}\sum_j\sqrt{V_j(x)}.          \tag{8}
$$

Both laws meet the assumption of (7) with
`V_*=sup_x Σ_j V_j(x)`: bounded differences gives it for signs,
and the Euclidean Lipschitz bound gives it for Gaussians.

Consequently, when q≈εn, V_*=O(qn), log|C|=o(n), and
Q(B)=O(q^(3/2)), the leading linear-in-ε cap cost is the actual
penalized column-mean profile L_C. The new-child geometry affects
only the O(ε^(3/2)) term in this sandwich. Small code entropy controls
selection fluctuations but does not bound those means.

The aggregate covariance certificate alone yields
`m_G(x)≤sqrt(2/pi) q sqrt((1+o(1))n)`. Its leading constant
sqrt(2/pi) exceeds the dilution budget `(3/2)c` at the reported
near-optimal c. Thus universality plus this covariance bound does not
close the value step. Conversely, the random code in Section 1 shows
that an entropy-only theorem cannot always reduce the physical mean
below that constant. Actual quadratic near-energy geometry remains
the essential additional input.
