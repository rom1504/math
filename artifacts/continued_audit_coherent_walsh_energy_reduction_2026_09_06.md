# Exact coherent Walsh structure and first-degree energy reduction

Date: 2026-09-06. Independently derived by the audit agent and cross-checked
with the feedback agent. This sharpens the retained-coherent energy term
without substituting a Gaussian law for Boolean seeds.

## 1. Polynomial theorem

Let S be an independent Boolean seed vector and let M_1,...,M_d be finitely
many matrices with n rows and uniformly bounded operator norms. The input
dimension can be n, dn for fixed colors, or any common dimension. Let
`D_i=P((M_1S)_i,...,(M_dS)_i)` for a fixed polynomial P. Deterministic
bounded row-dependent coefficients from a fixed finite polynomial family
are also allowed. Write D^(q) for its exact Boolean Walsh degree-q part.

For every fixed q>=2,

```math
\max_i\sum_j|\mathbb E D_i^{(q)}D_j^{(q)}|
+\max_j\sum_i|\mathbb E D_i^{(q)}D_j^{(q)}|=O(1). \tag{1}
```

The same bound holds for cross-covariances of two fixed coherent polynomial
families at common Walsh degree q>=2. This is stronger than a covariance
operator bound and uses the exact Boolean expansion.

Let B be an n-by-n symmetric hollow signing normalized by sqrt(n-1), with
fixed bounded operator norm. Put

```math
\mu_i=\mathbb E D_i,\qquad K_{ia}=\mathbb E[S_aD_i].
```

Then

```math
\frac{\mathbb E D^{\mathsf T}BD}{2n}
=\frac{\mu^{\mathsf T}B\mu}{2n}
 +\frac{\operatorname{Tr}(BKK^{\mathsf T})}{2n}+o(1). \tag{2}
```

For an odd coherent response, mu=0. The projection in (2) is the FIRST
BOOLEAN WALSH projection, not the first Gaussian coefficient of P under
a fictitious Gaussian coherent law.

## 2. Exact covariance proof

The exact Walsh-kernel argument in the retained-coherent audit expresses
each positive degree q as a finite sum of distinct-marked-slot products
of row matrices U_1,...,U_q, with bounded row scalars. Every U_l has bounded
operator norm. The same holds for a second kernel V_1,...,V_q. Match their
marked sets bijectively; at each fixed matching permutation, covariance
is a distinct-label sum

```math
\sum_{a_1,\ldots,a_q\ {m distinct}}
\prod_{\ell=1}^q U_{\ell,ia_\ell}V_{\ell,ja_\ell}. \tag{3}
```

Expand distinctness by inclusion-exclusion over set partitions pi of the
q slots. The resulting term is

```math
\prod_{C\in\pi}\Gamma_C(i,j),\quad
\Gamma_C=U_CV_C^{\mathsf T},\quad
U_C=\mathop{\circ}_{\ell\in C}U_\ell,\quad
V_C=\mathop{\circ}_{\ell\in C}V_\ell. \tag{4}
```

Each U_C and V_C has bounded operator norm. For |C|>=2 it also has bounded
absolute row and column sums, by Cauchy--Schwarz on two entrywise factors.
The bounded row scalars do not alter any conclusion.

If pi has at least two blocks, the matrices Gamma_C have bounded row and
column Euclidean norms. Retain two factors in Cauchy--Schwarz and bound
the others entrywise; the product in (4) has bounded absolute row and
column sums. If pi has one block, its size is q>=2. Both U_C,V_C then
have bounded absolute row AND column sums; their product U_C V_C^T has
the same bounded-sum property. This proves (1), including all internal
coherent collision corrections.

Different Walsh degrees are exactly orthogonal across all root pairs.
For q>=2, flatness and hollowness give

```math
\left|\frac{\operatorname{Tr}(B\operatorname{Cov}(D^{(q)}))}{2n}\right|
\le\frac1{2n\sqrt{n-1}}\sum_{ij}
 |\mathbb E D_i^{(q)}D_j^{(q)}|=O(n^{-1/2}).
```

The degree-zero and degree-one covariances are exactly mu mu^T and KK^T.
This proves (2).

## 3. Bounded coherent closure

The energy conclusion extends to fixed bounded continuous functions of
the finite coherent field family. The fields are uniformly subgaussian.
A uniform L2 polynomial approximation can be obtained as follows: choose
a compact box whose complement has uniformly small probability; extend
the bounded response periodically on a larger box while preserving its
bound; approximate it by bounded Fejer trigonometric polynomials; and
then approximate their finitely many Fourier exponentials by Taylor
polynomials in L2 using uniform subgaussian moment bounds. The bounded
trigonometric stage is useful because its norm on the tail does not grow
uncontrollably with the approximation degree.

For every fixed polynomial, take the matrix limit first. Then remove the
L2 approximation error. Both the degree-zero and first-Walsh projections
are contractions in rowwise L2. Matrix Cauchy--Schwarz and the fixed
operator cap of B therefore transfer (2), even without a covariance-op
bound for the limiting nonpolynomial response.

In the intended application

```math
c_i^0=H((BS)_i)\mathbb E_N\psi(b(QS)_i+\sigma_iN),
```

the deterministic variance profile stays in a bounded interval; finite
variance bins supply the uniform row-dependent approximation. H may be
bounded even Gaussian-a.e.-continuous: approximate H separately under the
standard Gaussian law, using the common standardized Rademacher-row-sum
law of (BS)_i. The other factor is bounded and Lipschitz. Consequently

```math
\frac{\mathbb E(c^0)^{\mathsf T}Bc^0}{2n}
=\frac{\operatorname{Tr}(B K_0 K_0^{\mathsf T})}{2n}+o(1),
\qquad (K_0)_{ia}=\mathbb E[S_a c_i^0]. \tag{5}
```

This makes one part of the retained energy a deterministic first-Walsh
kernel trace. It does not dispose of the separate c0--Z bare-star cross,
and K_0 still uses the actual coherent Boolean law.

## 4. A useful transported flattening, with a critical scope warning

For a fixed positive-degree coherent Walsh kernel A_(q,a), transport it
by the flat B: `L_(q,i)=sum_a B_ia A_(q,a)`. The same exact partition
expansion verifies small proper flattenings. If no equality block
straddles a cut, the two group root maps are bounded and the middle
diagonal b_i supplies O(n^-1/2). If a block straddles the cut, the
flattening is block-diagonal in its shared labels. In each block pull
out the flat B entry; one straddling product has at least two original
row factors, so its absolute COLUMN sum is bounded. Remaining shared
factors and uncontracted row Hilbert norms are bounded. The resulting
block operator norm is again O(n^-1/2).

This estimate is promising for a generalized transported-Walsh theory.
It does NOT justify restoring pre-Walsh tensor diagonals at small Hilbert
cost. An actual-signing counterexample is any normalized conference matrix
B, with Q=B²=I. The coherent input `(QS)_i³=S_i` has ZERO degree-three
Walsh part. But the unrestricted Gaussian transported tensor

```math
\sum_a B_{ia} e_a^{\otimes3}
```

has row Hilbert norm one, supported entirely on repeated marked labels.
Deleting those repeats therefore costs one, despite its small proper
flattenings. Any generalized star proof must retain the exact source
Walsh projection and its diagonal-partition terms. It cannot import
the scalar raw-Hermite collision restoration without a new argument.

Only (1)–(5) are used here as energy-reduction theorems. A general
next-stage feedback law is not claimed by this artifact.

## 5. Reproducible finite identity check

The script `computations/continued_audit_coherent_walsh_identity_2026_09_06.py`
checks the exact cube decomposition for V=QS. Its first-Walsh matrix is

```math
K_{ia}=3(QQ^{\mathsf T})_{ii}Q_{ia}-2Q_{ia}^3,
```

and the degree-three covariance is exactly

```math
6\left[(QQ^{\mathsf T})^{\circ3}
-3(QQ^{\mathsf T})\circ
       ((Q^{\circ2})(Q^{\circ2})^{\mathsf T})
+2(Q^{\circ3})(Q^{\circ3})^{\mathsf T}\right].
```

All 64 Boolean seeds of the order-six Steiner signing verify the two
formulas, evaluated to floating roundoff. The script also checks the
identity-matrix cube has zero third-Walsh part, and displays finite
Steiner scaling through n=496. The unbounded polynomial cube's numerical
energies are not Boolean-cap certificates; they only test the exact
Walsh decomposition and the decay mechanism.
