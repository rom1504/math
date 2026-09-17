# A nonvacuous full-sign application: spectrally pinned Paley perturbations

2026-09-17 campaign. This note supplies an actual infinite family, not
an assumed geometry of minimizers. Its cap coefficient can be any
c>1/2. A full-sign bridge preserves the old matrix and improves its
normalized cap with zero first-order unnormalized extension cost.
This does **not** improve the known unrestricted upper bound below
1/2. The last section proves an obstruction in this certificate itself.

The BH/mechanism researcher independently read and reconstructed the
entire frozen proof and returned **PASS**, including the Fourier
sampling, both energy sectors, paired-bridge constants, and certificate
floor. The finite replay also passed.

## 1. Full-sign matrices with a pinned Boolean maximizer

Fix λ>1. Along N=5^a tending to infinity there are hollow symmetric
full-sign matrices A_N such that

$$
A_N{\bf1}=d_N{\bf1},\qquad
\frac{d_N}{\sqrt N}\longrightarrow\lambda,\qquad
\left\|A_N\big|_{{\bf1}^\perp}\right\|_{\rm op}
\le(1+o(1))\sqrt N.                                   \tag{1}
$$

Consequently, for all sufficiently large N,

$$
Q(A_N)=\frac{Nd_N}{2},\qquad
\frac{Q(A_N)}{N^{3/2}}\longrightarrow\frac\lambda2.       \tag{2}
$$

The only absolute maximizers are ±1. The construction and its spectral
estimate are elementary and are given next.

### 1.1 A sparse random change of Paley generators

Work in the finite field F_N. Let χ be its quadratic character, extended
by χ(0)=0. Since N≡1 mod4, χ(−1)=1. The Paley sign matrix

$$
P_{uv}=\chi(u-v)
$$

is real symmetric, hollow and full sign, with P1=0. For any nontrivial
additive character ψ, put G=Σ_t χ(t)ψ(t). Then G is real and |G|=√N.
For completeness, conjugation gives Ḡ=χ(−1)G=G, while changing
variables t=as in |G|² gives

$$
|G|^2=\sum_a\chi(a)\sum_{s\ne0}\psi((a-1)s)=N.
$$

Thus the additive-character eigenvalues of P off 1 are χ(b)G,
all of magnitude √N.

There are M=(N−1)/4 unordered negative-generator pairs {t,−t}.
Choose m=⌊λ√N/4⌋ of them uniformly without replacement, and reverse
the signs of every corresponding Cayley edge. The resulting A has
every off-diagonal entry still ±1, and has exact row sum d_N=4m.

For nonzero frequency b, the eigenvalue perturbation is a sum of m
sampled real numbers

$$
w_b(\{t,-t\})=2(\psi(bt)+\psi(-bt))
            =4\operatorname{Re}\psi(bt)\in[-4,4].
$$

The sum over the whole population is −1−χ(b)G. Hence its sampled
mean is −4m(1+χ(b)G)/(N−1), of absolute value at most
4m(1+√N)/(N−1).

Reveal the sampled pairs in random order. The Doob martingale for the
final sum has increments bounded by eight. Indeed, after i−1 draws,
revealing the next population value changes the conditional final mean
by (M−m)/(M−i) times its deviation from the remaining population mean,
whose absolute value is at most eight. For large N, m<M.
The elementary bounded-increment exponential inequality gives

$$
\Pr\{|Z_b-\mathbb E Z_b|\ge t\}
\le2\exp[-t^2/(128m)].
$$

Use t=√(512m log N) and union bound over the N−1 nontrivial
frequencies. The failure probability is less than 2N^(−3)<1. Fix a
successful subset. It obeys the explicit estimate

$$
\left\|A\big|_{{\bf1}^\perp}\right\|_{\rm op}
\le \sqrt N+\frac{4m(1+\sqrt N)}{N-1}
                +\sqrt{512m\log N}
=(1+o(1))\sqrt N.                                    \tag{3}
$$

This proves (1). The randomization is only a proof of existence of
a deterministic generator subset; no matrix entries are fractional.
Unlike independent edge changes, this Cayley modification preserves
exact regularity.

### 1.2 The full near-level code is a narrow two-center tube

More generally let a full-sign matrix A obey A1=d1, and let
s=||A|_(1⊥)||op<d. For x∈{−1,1}^N, put
ρ=N^(−1)min(d_H(x,1),d_H(x,−1))≤1/2. The squared normalized
projection onto 1 is (1−2ρ)². Therefore

$$
H_A(x)\le Q(A)-2(d-s)N\rho(1-\rho),\qquad
-H_A(x)\le\frac{sN}{2}.                              \tag{4}
$$

These inequalities prove (2) and uniqueness of its two maximizers.
For the family (1), write g=λ−1. For any fixed 0<η<g/2, every word
of the full absolute near-level code lies within radius
(r(η)+o(1))N of ±1, where

$$
r(\eta)=\frac{1-\sqrt{1-2\eta/g}}2
       =\frac{\eta}{2g}+O(\eta^2).                    \tag{5}
$$

This is a genuine complete-code statement: the negative-energy
sector is excluded by its fixed gap, not dropped by assumption.
It is consistent with the universal ηlog(1/η) entropy floor,
because the tube itself has that much local entropy.

## 2. An exact constraint-preserving full-sign bridge

Let q=⌊εN⌋. Pair the first 2⌊N/2⌋ old vertices. Independently in
each new column take

$$
C_{2i-1,j}=\xi_{ij},\qquad C_{2i,j}=-\xi_{ij},
$$

with all ξ fair independent signs. If N is odd, fill the last
entry of each column by an independent fair sign. Every bridge entry
is an exact sign. For even N, Cᵀ1=0 exactly; for odd N its coordinates
have magnitude one. Thus the two old maximizing responses are
preserved up to the negligible odd-order error q.

For a word x at distance r from its nearer uniform word, at most r
pairs have a nonzero coefficient x_(2i−1)−x_(2i)=±2.
Writing B(x)=Σ_j|xᵀC_j|, the paired part has

$$
\mathbb E B_{\rm pair}(x)\le2q\sqrt r.
$$

Changing one active scalar driver changes B_pair by at most four.
The bounded-difference MGF bound therefore has variance proxy at most
4qr, and its upper tail at logarithmic confidence u is
√(8qru). For each r there are at most 2 binom(N,r) words. Union over
all radii, using u_r=log(2binom(N,r))+3log(N+1), proves existence of a
single bridge with, simultaneously for every old word,

$$
B(x)\le q+2q\sqrt r+
 \sqrt{8qr[\log(2\binom Nr)+3\log(N+1)]}.               \tag{6}
$$

At r=0 the paired part is identically zero. Consequently, uniformly
in 0≤ρ≤1/2 and after N→∞ at fixed ε,

$$
\frac{B(x)}{N^{3/2}}
\le b_\epsilon(\rho)+o(1),\qquad
b_\epsilon(\rho)=2\epsilon\sqrt\rho+
                  \sqrt{8\epsilon\rho h(\rho)}.       \tag{7}
$$

The maximization over all new signs is exact in B. No deletion of
old words or physical edge support is used.

## 3. Actual parent bound and zero linear extension cost

Let D_q be any sequence of hollow full-sign new blocks with
limsup Q(D_q)/q^(3/2)≤c_D. For every fixed ε>0 the construction gives

$$
\begin{aligned}
\limsup\frac{Q\left(\begin{smallmatrix}A_N&C\\C^T&D_q\end{smallmatrix}\right)}
                 {N^{3/2}}
\le c_D\epsilon^{3/2}
+\max\Big\{&
 \frac\lambda2+
 \sup_{0\le\rho\le1/2}
 [b_\epsilon(\rho)-2(\lambda-1)\rho(1-\rho)],\\
&\sup_{0\le\rho\le1/2}
 [\tfrac12-\tfrac{\lambda+1}{2}(1-2\rho)^2
                                    +b_\epsilon(\rho)]
                         \Big\}.                    \tag{8}
\end{aligned}
$$

Both sectors are retained. Formula (8) follows from (4), its lower
Rayleigh bound, (7), and the pointwise inequality

$$
Q(\text{parent})\le Q(D_q)+\max_x[|H_A(x)|+B(x)].
$$

For fixed g=λ−1>0 and sufficiently small ε, the negative sector in
(8) is below λ/2, because its deficit is at least g/2 and
max_ρ b_ε(ρ)≤√2 ε+2√(εlog2)→0.

The positive excess has the quantitative bound

$$
F_g(\epsilon):=
\sup_{0\le\rho\le1/2}
[b_\epsilon(\rho)-2g\rho(1-\rho)]
\le \frac{2\epsilon^2}{g}+o_g(\epsilon^2).             \tag{9}
$$

To check it, use h(ρ)≤ρlog(e/ρ) and
2gρ(1−ρ)≥gρ. Put ρ₀=exp(1−g²/(32ε)).
For ρ≥ρ₀, the entropy term is at most gρ/2, and

$$
2\epsilon\sqrt\rho-\frac g2\rho\le\frac{2\epsilon^2}{g}.
$$

For ρ≤ρ₀, the remaining positive terms are bounded by
2ε√ρ₀+(g/2)ρ₀=o_g(ε²). This proves (9).

In particular, genuine full-sign parents exist with

$$
Q(\text{parent})\le Q(A_N)+
  \left[c_D\epsilon^{3/2}+O_g(\epsilon^2)+o_N(1)\right]N^{3/2}.
                                                               \tag{10}
$$

One may take c_D=1 by the elementary random-sign union bound, so the
application does not depend on a conjectural or optimized child family.
The order is N→∞ first, then ε↓0 with fixed λ.
After normalization by (N+q)^(3/2), (10) is strictly below λ/2
for every sufficiently small ε>0: its leading change is
−3λε/4. The old principal block is unchanged.

This supplies an unconditional positive model for zero-linear-slope
augmentation. It is not an iteration theorem: the new parent need not
retain the spectral pinning hypothesis.

## 4. A rigorous barrier in this particular certificate

The positive-sector supremum in (8), evaluated at ρ=1/2, is at least

$$
\frac12+\sqrt2\,\epsilon+2\sqrt{\epsilon\log2}.
$$

Thus the numerical right side of the normalized certificate (8) is
always at least

$$
\frac{\frac12+\sqrt2\,\epsilon
      +2\sqrt{\epsilon\log2}+c_D\epsilon^{3/2}}
     {(1+\epsilon)^{3/2}}.                            \tag{11}
$$

For any 0≤c_D≤1/2 and any ε>0, (11) is strictly larger than c_D.
Indeed,

$$
(1+\epsilon)^{3/2}-\epsilon^{3/2}
=\frac32\int_0^1\sqrt{\epsilon+t}\,dt
\le1+\frac32\sqrt\epsilon,
$$

whereas the first three numerator terms in (11) strictly exceed
c_D(1+(3/2)√ε). Therefore optimizing λ and ε cannot improve
the supplied child constant. With c_D=1/2, this certificate cannot
cross 1/2. Using a known child bound below 1/2 can inherit that
improvement, but cannot improve it.

This is a **certificate limitation**, not a lower bound on the actual
caps of these parents. It separates a real nonvacuous application from
an unsupported claim about asymptotic minimizers.

## 5. Finite replay and declared verification limits

`computations/paper_discrepancy_2026_09_17_pinned_signings.py` constructs
sampled Paley perturbations at prime orders 13, 17, 29, 101 and 257,
checks every Fourier mean and the displayed simultaneous-error bound,
and verifies exact regularity and full sign support. At orders 13 and
17 it enumerates every spin pair, the full near-level spectral tube,
and every new spin for a three-vertex paired-sign bridge. The exact
caps are respectively 52→56 at orders 13→16 and 102→106 at
orders 17→20, lowering their normalized caps. These are demonstrations,
not optimal signings.

The same replay checks the certificate floor on 4,004 parameter pairs.
Fourier norms in this diagnostic use floating-point arithmetic; the
asymptotic existence and certificate-floor proofs above are analytic
and do not rely on those numerical values. The cap enumerations and
matrix-entry checks use integer arithmetic.

## 6. Archive collision check

The archive already contains
`paley_correlated_shell_optimization_audit.md`, auditing
`perturbed_conference_stratified_entropy.md`. Its perturbations are
independent edge-sign flips, and its positive conclusion is conditional
on an unproved deterministic Paley shell profile. It also correctly
distinguishes a scalar-union certificate barrier from a family barrier.

The construction here instead flips entire negative Cayley-generator
pairs: it preserves exact regularity, creates a dominant uniform
eigenvector, and proves the needed full near-level tube by a spectral
gap. Its purpose is an actual nonoptimal bounded-cap application,
not the archive's conditional sub-half optimization. The present
paired-column bridge and its certificate floor are separately derived
above. No external literature priority is claimed for Paley perturbation
as a broad idea.
