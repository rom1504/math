# Universal entropy from a small isotropic sign response

2026-09-17 paper-combination campaign. Developed by the discrepancy
agent; the root and BH/mechanism agents independently reconstructed the
proof and reported **PASS**. All logarithms are natural. This is a
self-contained theorem, not a claim of literature priority or an
unconditional improvement to the full-sign quadratic-cap problem.

## 1. The finite theorem

Let ν be any probability law on {−1,1}^N satisfying E hhᵀ=I_N.
Its mean need not be zero. Suppose a nonempty code C⊆{−1,1}^N satisfies

$$
\max_{x\in C}\mathbb E_\nu|h\cdot x|
\le A=\delta\sqrt N.
$$

Then A≥1. For every integer m≥2A² log 2,

$$
\boxed{\log|C|\le m\log(1+2/\delta^2)
       +N h\left(e^{-m/(2A^2)}\right),}                 \tag{1}
$$

where h is binary entropy. In particular, for 0<δ²≤1/2,

$$
\boxed{\frac{\log|C|}{N}
\le\delta^2\bigl(2\log(1/\delta^2)+1\bigr)
                       \log(1+2/\delta^2)+h(\delta^2)
=\bigl(8+o(1)\bigr)\delta^2\log^2(1/\delta).}            \tag{2}
$$

The bound is uniform over N, the sign law, its support size, and the
code. Thus δ_N→0 implies log|C_N|=o(N) with **no aspect-ratio
or atom-count hypothesis**.

### Proof: an importance-weighted information channel

Write a_x=E|h·x|. Isotropy and |h·x|≤N imply
N=E(h·x)²≤Na_x, so A≥a_x≥1.

Take a common reference law R=ν×Uniform{−1,1} on (h,σ).
Define Q_x by the nonnegative density

$$
q_x(h,\sigma)
=1-\frac{a_x}{A}+\frac{2(\sigma\,h\cdot x)_+}{A}.        \tag{3}
$$

Averaging over σ proves that it integrates to one. Furthermore,

$$
\mathbb E_{Q_x}\sigma h
=\frac1A\mathbb E_\nu h(h\cdot x)=\frac{x}{A}.           \tag{4}
$$

For b=2(σh·x)_+, one has E_R b=a_x and E_R b²=2N.
Consequently

$$
\mathbb E_Rq_x^2
=1+\frac{2N-a_x^2}{A^2}
\le1+\frac2{\delta^2},
\qquad
D(Q_x\Vert R)\le\log\mathbb E_Rq_x^2
\le\log(1+2/\delta^2).                                 \tag{5}
$$

The relative-entropy inequality is Jensen applied to log q_x under
Q_x; zero-density values contribute zero.

Let X be uniform on C. Conditional on X=x, draw m independent
pairs from Q_x, and denote the resulting observation by Z. The
common product reference gives

$$
I(X;Z)\le\mathbb E_X D(Q_X^m\Vert R^m)
       \le m\log(1+2/\delta^2).                         \tag{6}
$$

For each coordinate i, predict X_i by the majority of the sampled
signs σ_j h_(j,i). Conditional on X=x, their mean is x_i/A
and they are independent. Hoeffding's inequality bounds the prediction
error by ε=exp(−m/(2A²)). Since ε≤1/2, the binary
coordinate entropy bound gives

$$
H(X\mid Z)\le\sum_i h(\Pr\{\widehat X_i\ne X_i\})
            \le N h(\epsilon).
$$

Adding to (6) proves (1). Choose m=⌈2A² log(1/δ²)⌉.
Its error is at most δ² and, because A²≥1, it obeys
m≤A²(2log(1/δ²)+1). This proves (2), including integer rounding.

The encoding law depends on the codeword. It is used only to prove
entropy; it is **not** a physical bridge construction that may depend
on a query.

## 2. A sharp-scale Walsh lower example

Let N=4^r, with coordinates indexed by F₂^r×F₂^r, and let
ν be uniform on the antipodal Walsh characters. This is a centered
isotropic sign law. For each binary r×r matrix T, put

$$
S_T=\{(u,Tu):u\in\mathbb F_2^r\},\qquad
f_T=1-2\,1_{S_T}.
$$

There are 2^(r²) distinct words. The normalized Fourier transform of
1_(S_T) is 2^(−r) on its annihilator and zero elsewhere. Thus

$$
\mathbb E_\nu|h\cdot f_T|
=\|\widehat f_T\|_1=3-2^{2-r}\le3,
\qquad
\log|\{f_T\}|=\frac{(\log N)^2}{4\log2}.                 \tag{7}
$$

For δ=3/√N, the upper bound (2) is also of order
(log N)². Hence its double-logarithmic dependence is sharp in this
joint regime, up to constants. In particular, a putative universal
Cδ²Nlog(1/δ)+O(log N) bound is false.

Even replacing that additive remainder by O((log N)²) does not
repair the single-log main term: take t=2^s independent graph choices
in disjoint slabs of F₂^s×F₂^r×F₂^r. There are
2^(tr²) words and each has response at most 1+2t, by the Fourier
triangle inequality. With t the largest power of two below √r,
N=t4^r and δ=(1+2t)/√N, the claimed repaired upper bound
would be O(r²) while the logarithmic cardinality is Ω(r^(5/2)).

Affine-subspace indicators as small-spectral-norm Boolean functions
are classical; see Green–Sanders,
[*Boolean functions with small spectral norm*](https://arxiv.org/abs/math/0605524),
and Shpilka–Tal–Volk,
[*On the Structure of Boolean Functions with Small Spectral Norm*](https://arxiv.org/abs/1304.0371).
The exact count and all constants in (7) are proved above, without
importing their structural theorems. A preliminary primary-source
search located those Walsh-specific results but did not establish a
literature novelty claim for the support-free isotropic theorem (1).

## 3. Matrix-array corollary without an aspect-ratio restriction

Let ν now be isotropic on {−1,1}^k, let N=kp, and identify a sign
word with a k×p sign matrix X. Define

$$
\Phi_\nu(X)=\frac1{\sqrt N}
            \mathbb E_{h\sim\nu}\|h^T X\|_2.
$$

Take g uniform on {−1,1}^p, independent of h. The physical sign
law V=h⊗g is centered and isotropic in dimension N, and

$$
\mathbb E|\langle V,X\rangle|
\le\mathbb E_h\|h^TX\|_2=\sqrt N\,\Phi_\nu(X).
$$

Therefore the **entire** Boolean code with Φ_ν(X)≤μ satisfies
(2) with δ=μ. In particular its entropy is
O(μ²log²(1/μ))N, for all k,p, with no condition such as
k²=o(N). This is distinct from the sharper flat-diagonal Gaussian
information argument, which gives asymptotic
(1+o(1))μ²log(1/μ²) after its covariance-label entropy vanishes
under k²=o(N). Both statements are support-free.

If fewer than k coordinates are left over, independent fair signs on
those coordinates preserve isotropy, and their response contributes
at most the square root of their number. Alternatively append their
at-most-klog2 entropy to a main-array estimate when k=o(N).

## 4. Scope for the quadratic-cap campaign

The theorem shows that an arbitrary isotropic sign law cannot have
vanishing normalized mean response on an exponentially large Boolean
code. It is not by itself a full-sign augmentation theorem: entropy
control and a realizable bridge with sharp concentration remain
different requirements. For the campaign's matrix-array laws,
support-free sampled-mode realization supplies the latter in the
k²=o(N) regime; there the stronger flat-diagonal bound is preferable.

There is now a stronger scope correction: the full near-level entropy
slope is universally infinite for bounded positive-cap quadratic
sequences. Thus the earlier finite full-code entropy-slope deployment
has no examples. Low-response **center codes with a Hamming residual
cover** remain viable, including critical k=O(√N) sampled-mode
realization. See [the exact entropy-floor correction](paper_discrepancy_nearlevel_entropy_floor_2026_09_17.md).
Neither convergence nor a new unrestricted cap constant follows.

The finite proofs and lower examples are also recorded in Sections
20–21 of [the discrepancy campaign artifact](paper_discrepancy_2026_09_17.md).
Reproducible exact checks:

- computations/paper_discrepancy_2026_09_17_scalar_entropy_counterexample.py:
  all graph maps through r=3, plus every two-slab map at r=2.
- computations/paper_discrepancy_2026_09_17_scalar_entropy_encoding.py:
  exact density normalization, all signed means, second moments, and
  numerical relative-entropy bounds for all 530 tested graph words.
