# An all-law near-level obstruction survives small principal extensions

2026-09-17 campaign. Discrepancy-track deduction from the localization
track's quadratic response dual. The localization researcher independently
reconstructed the finite lifting and covariance argument and returned
**PASS**. This is a limitation of generic regularization on an explicit
near-half family, not a statement about unknown asymptotic minimizers.
The same researcher subsequently audited Sections 4--5, including the
covariance-free response and arbitrary-future-bridge bound, with PASS.

## 1. A finite principal-extension lemma

Let A be a hollow full-sign matrix of order n. Suppose a probability
law π is supported on Boolean words x with

$$
Q(A)-|H_A(x)|\le a,
$$

and suppose, for every isotropic sign law ν on {−1,1}^n,

$$
\mathbb E_{x\sim\pi}\mathbb E_{h\sim\nu}|h\cdot x|
\ge b\sqrt n.                                         \tag{1}
$$

Let P be ANY hollow full-sign principal extension of A by q vertices,

$$
P=\begin{pmatrix}A&C\\C^T&D\end{pmatrix},\qquad
\Delta=Q(P)-Q(A).
$$

Then Δ≥0. There is a deterministic lift x↦z_x=(x,y_x), chosen
from P and x without reference to a queried h, such that

$$
Q(P)-|H_P(z_x)|\le\Delta+a.                            \tag{2}
$$

For EVERY isotropic law ν_P on {−1,1}^(n+q), these lifts obey

$$
\mathbb E_{x\sim\pi}\mathbb E_{h\sim\nu_P}|h\cdot z_x|
\ge b\sqrt n-\sqrt q.                                 \tag{3}
$$

The right side may be replaced by its positive part.

### Proof

For each x choose σ_x with σ_xH_A(x)=|H_A(x)|.
Average the signed parent energy over uniform y. Both its new
quadratic part and bridge part have average zero, hence

$$
\mathbb E_y\,\sigma_x H_P(x,y)=|H_A(x)|.
$$

Some y_x attains at least this average. This proves (2), and applying
the same observation to an old maximizer proves Δ≥0.

The old marginal of an isotropic parent sign law is isotropic. The new
marginal is also isotropic, irrespective of cross-block dependence.
For each fixed lift, triangle inequality and Cauchy–Schwarz give

$$
\mathbb E|h\cdot z_x|
\ge\mathbb E|h_{\rm old}\cdot x|
    -\mathbb E|h_{\rm new}\cdot y_x|
\ge\mathbb E|h_{\rm old}\cdot x|-\sqrt q.
$$

Average over π and apply (1). This proves (3).
No independence between the old and new coordinates of h is used.
No estimate on Q(D), nor any particular bridge construction, is needed.

## 2. Exact consequence for low-response center covers

Put N=n+q. Let B⊆{−1,1}^N be any nonempty center family and suppose

$$
\sup_{f\in B}\mathbb E_{\nu_P}|h\cdot f|\le\mu\sqrt N.
$$

If every word in the parent near-level set of deficit Δ+a is within
Hamming distance rN of B, then every lift is. Isotropy gives

$$
\mathbb E|h\cdot z_x|
\le(\mu+2\sqrt r)\sqrt N.
$$

Comparing with (3) proves the necessary radius bound

$$
\boxed{\displaystyle
r\ge\frac14\left[
 b\sqrt{\frac nN}-\sqrt{\frac qN}-\mu
                         \right]_+^2.}               \tag{4}
$$

The center family may be exponentially large. The sign law may be
selected after seeing the entire parent, code, and centers.

## 3. The explicit Hadamard-family corollary

The localization track's
[all-law Hadamard obstruction](paper_localization_all_law_hadamard_obstruction_2026_09_17.md)
constructs, for n=p² and p a power of two, a hollow full signing
A=F−I with

$$
Q(A)=\frac{n(p+1)}2,\qquad
a=n,\qquad
b=b_p=\frac1{\sqrt{6[\,2+p/(p-3)\,]}}.
$$

Its law π is the explicitly defined random-involution mixture of
both global eigensectors. The complete proof of (1) there uses a
quadratic correction, removed exactly by isotropy; it does not assume
π itself is isotropic.

Consider ANY sequence of principal extensions P_N of this family with

$$
q=o(n),\qquad Q(P_N)-Q(A_n)=o(n^{3/2}).
$$

For every fixed η>0, the threshold ηN^(3/2) eventually contains all
the lifts (2), because Δ+n=o(n^(3/2)). Hence if μ_N→0, the Hamming
radius needed to cover the full η-near-level code by centers of
response at most μ_N√N satisfies

$$
\liminf_N r_N(\eta)\ge\frac1{72}.                      \tag{5}
$$

This holds for every isotropic physical sign law, after arbitrary
switching, permutation, mode grouping, or adaptation to the parent.
It survives every full-sign regularization satisfying the two
displayed small-extension conditions.

The order of limits is essential: η is fixed before N→∞.
A regularization theorem can validly control a window η_N→0 and
give sublinear Gaussian width there without contradicting (5).
Neither small entropy nor sublinear Gaussian width alone implies
a low-response isotropic center law.

Finally this family's cap tends to 1/2, above the reported unrestricted
upper bound below 0.494. Thus (5) rules out a universal preprocessing
argument, but does not prove the same obstruction for actual asymptotic
minimizers or settle the original convergence question.

## 4. A covariance-free strengthening for smaller extensions

Section 7 of the linked all-law obstruction strengthens its old ground
law: there is a law π_mix on words of deficit at most n such that, for
EVERY deterministic physical sign h,

$$
\mathbb E_{x\sim\pi_{\rm mix}}|h\cdot x|\ge c_p,
\qquad c_p\sim\sqrt n/\sqrt{18}.                      \tag{6}
$$

For the precise finite constant, put

$$
C_p=2+p/(p-3),\quad
a_p=\frac{p^2}{(p-1)\sqrt{6C_p}},\quad
d_p=\frac{p^{3/2}}{\sqrt3},\quad
c_p=\frac{a_pd_p}{a_p+d_p}.
$$

This mixes the identity-involution positive ground law with the
fixed-point-free two-eigensector law. Neither π_mix nor a queried law
is required to be isotropic in (6).

Use the deterministic lift from Section 1 for this law. For every
parent physical sign h, the pointwise triangle inequality gives

$$
\mathbb E_x|h\cdot z_x|\ge c_p-q,                    \tag{7}
$$

because |h_new·y_x|≤q. Thus if q=o(√n) and Δ=o(n^(3/2)), every
full fixed-η parent nearcode has maximum response at least
(1/√18−o(1))√N against EVERY sign law, with no covariance condition.

This does not give a covariance-free Hamming-radius lower bound:
the `2 sqrt(rN)` continuity estimate in Section 2 still needs isotropy
(or a specified covariance bound). Equation (7) itself is unconditional.

## 5. Actual microscopic-regularity examples still have a bridge cost

Here is a concrete application to the director's
[sparse random-star regularizer](paper_director_sparse_random_regularization_2026_09_17.md).
Take the Hadamard A_n above as the retained old block, and adjoin
q₀→∞ vertices with q₀³≤n using that theorem's independent sign bridge
and a low-cap child. Its proof applies to an arbitrary retained block;
no pre-existing full matrix on n+q₀ vertices is needed for this
extension formulation. Some realized full signing W_N, N=n+q₀, has

$$
\Delta=Q(W_N)-Q(A_n)=O(n\sqrt{q_0}),\qquad
w\bigl(E_{W_N}(N/\sqrt{q_0})\bigr)=O(N/\sqrt{q_0}),   \tag{8}
$$

and the corresponding subexponential cardinality bound. In particular
q₀≈n^(1/3) gives a genuine full-sign sequence whose window and Gaussian
width are both of order N^(5/6), with cap tending to 1/2 after
normalization. One can also take q₀ of logarithmic order.

Now form ANY further full-sign principal extension P of W_N by r
vertices, with arbitrary bridge C' and arbitrary child D'. Then

$$
\boxed{\quad
Q(P)\ge Q(W_N)-(\Delta+n)+r(c_p-q_0)-Q(D').\quad}      \tag{9}
$$

To prove (9), use the mixed-law lifts z_x. Their W_N energy is within
Δ+n of Q(W_N), and (7) holds for every column of C'. Given z_x,
orient each new spin so that its bridge contribution has the same
sign as H_W(z_x). The resulting absolute parent energy is at least
|H_W(z_x)|+Σ_j|C'_j·z_x|−Q(D'). Average over x. This proves (9)
without independence, isotropy, or a probability law for C'.

Consequently, fix ε>0, take r/N→ε, and suppose
Q(D')≤L r^(3/2)+o(N^(3/2)) for a fixed finite L. Since
q₀=o(√n) and Δ=o(n^(3/2)),

$$
\liminf\frac{Q(P)-Q(W_N)}{N^{3/2}}
\ge\frac{\varepsilon}{\sqrt{18}}-L\varepsilon^{3/2}.  \tag{10}
$$

Thus microscopic small width and entropy do not universally imply
an actual low-cap extension with zero linear cost. The obstruction
holds against ALL choices of the future bridge, not only a specified
rounding algorithm or isotropic mode family. It does not exclude a
useful nonzero slope, and it is not a statement about the unknown
asymptotic minimizers: these examples still have normalized cap 1/2.

The relevant lifted code lies in a window Δ+n, generally much larger
than the controlled microscopic window N/√q₀. This explicit scale
separation is why the valid regularizer and (10) coexist.
