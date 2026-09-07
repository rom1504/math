# Independent audit: affine Gram laws miss the width of low-cap signings

Date: 2026-09-07. Status: proof and constants PASS.

This audits the director's proposed scalar-affine Gaussian witness
exclusion. The scope is signed expected energy, not expected absolute
energy and not arbitrary Gaussian covariances.

## 1. A finite statement with the inequality directions explicit

Let A be a hollow symmetric full-sign matrix of order n. Define
H_A(x)=x^T A x/2, P=max H_A, R=max(-H_A), w=(P+R)/2, and Q=max(P,R).
Suppose known bounds w>=w0 and Q<=q0 give L=2w0-q0>n/2. Then

```math
P\ge 2w-Q\ge L,\qquad R\ge 2w-Q\ge L.                 \tag{1}
```

For s>=0 such that G=I+sA is PSD, let X be the signs of a centered
Gaussian with covariance G. Every coordinate variance is one. With
mu=-lambda_min(A), the Rayleigh bound and the PSD restriction imply

```math
R\le n\mu/2\quad\Longrightarrow\quad
\mu\ge 2R/n\quad\Longrightarrow\quad
s\le1/\mu\le n/(2R)\le n/(2L)<1.                     \tag{2}
```

Thus taking reciprocals does not reverse the desired conclusion. The
exact bivariate Gaussian angle identity, including singular full G,
gives E X_i X_j=(2/pi)arcsin(s A_ij). Since A_ij is +/-1,

```math
0\le\mathbb E H_A(X)
 ={n(n-1)\over\pi}\arcsin(s)
 \le {n(n-1)\over\pi}\arcsin\!\left({n\over2L}\right). \tag{3}
```

Apply the same argument to -A, whose negative cap is P, to obtain the
identical bound on -E H_A for covariance I-sA. This proves an absolute
bound on the *signed expectation* for either affine covariance orientation.

## 2. Uniform asymptotics and the numerical threshold

The freshly audited universal width lower bound gives
W_n/n^(3/2)>=c-o(1), with c=.4333221116640807 as a rigorous downward
rounding of the certificate. Suppose a full-sign sequence satisfies
Q(A_n)/n^(3/2)<=u+o(1), u=.494515125. This includes any sequence of
asymptotic original minimizers, but does not require optimality.

Equation (1) gives both P and R at least
(2c-u-o(1))n^(3/2). Every admissible affine parameter is therefore
O(n^(-1/2)), uniformly over the chosen signs and Gaussian laws. Using
arcsin(s)=s+O(s^3) in (3) gives

```math
\limsup_n {\left|\mathbb E H_{A_n}(X_n)\right|\over n^{3/2}}
\le {1\over2\pi(2c-u)}
< .427687444510305
< c-.005634667153775.                                 \tag{4}
```

The normalized cubic remainder is O(n^(-1)): the prefactor after
normalization is (n-1)/(pi sqrt(n)), while s^3=O(n^(-3/2)). Thus there is
no nonuniform endpoint or limit-order issue with an n-dependent s.

For general positive c and u<2c, this ceiling is below c exactly when

```math
u<2c-{1\over2\pi c}.                                  \tag{5}
```

At the certified c, the threshold is approximately .499354074019383.
The new upper .494515125 satisfies (5); the previous upper
.499432220485404 does not satisfy this sufficient criterion.

The source width theorem is
`artifacts/decisive_audit_certified_minimum_width_lower_2026_09_07.md`.
No spectral deletion, near-ground entropy, or additional minimizer
property is used in the present deduction.

## 3. Mixtures and the precise obstruction scope

The bounds are uniform over all admissible parameters of both signs.
Consequently any convex mixture of these centered Gaussian-sign laws
also has |E H_A| bounded by (3) and asymptotically by (4). One may
equivalently choose an orientation before sampling the corresponding
law and bound the expectation of the oriented energy. A common positive
scaling of the Gaussian covariance does not change its sign law.

These statements do not bound E|H_A(X)|, and do not cover choosing
the orientation after seeing X, conditioning on an energy event,
nonzero Gaussian means, or nonlinear/anisotropic Gram matrices.

The Steiner trace obstruction used an affine Gram Gaussian-sign law
whose signed expected energy exceeded w by a fixed multiple of
n^(3/2). Equation (4) proves that *this particular separating-witness
mechanism* cannot occur for the above low-cap sequences. It does not
prove that a small-trace PSD ramp exists there: other Boolean laws or
other obstructions could still rule it out.

## 4. Exact numerical audit

`computations/decisive_audit_nearmin_affine_gram_constants_2026_09_07.py`
uses only rational arithmetic for its assertions. It encloses pi through
Machin's identity and alternating arctangent series, uses the rounded-down
c above, and verifies (4), the strict positive width margin, and both
comparisons with the threshold in (5). Printed decimal values are
diagnostic only. All exact assertions passed.

## 5. Robustness to normalized-Frobenius-small covariance changes

Let G0=I+sA be an admissible correlation matrix with |s|<=1/6, and
let G be any other correlation matrix. Set E=G-G0. Both diagonal
vectors are one, so E has zero diagonal. For any |u|<=1/6 and
v in [-1,1],

```math
|\arcsin(v)-\arcsin(u)|\le2|v-u|.                      \tag{6}
```

Here is an elementary global check, including v at an endpoint. If u
and v have opposite signs, the secant is a weighted average of the two
secants from zero, each at most pi/2<2. By symmetry it remains to
consider u,v>=0. When v<=u, both are at most 1/6, and the derivative
is at most 6/sqrt(35)<2. When v>=u, convexity makes the secant largest
at v=1, and the resulting secant is increasing in u. Its maximum is
therefore

```math
{\pi/2-\arcsin(1/6)\over1-1/6}\le {3\pi\over5}<2.
```

Applying the Gaussian angle identity entrywise and then Cauchy--Schwarz
over unordered edges gives the exact finite estimate

```math
\begin{split}
|\mathbb E H_A(\operatorname{sign}N(0,G))
 -\mathbb E H_A(\operatorname{sign}N(0,G0))|
&\le {4\over\pi}\sum_{i<j}|E_{ij}|\\
&\le {2\over\pi}\sqrt{n(n-1)}\,\|E\|_F
\le {2n\over\pi}\|E\|_F.                              \tag{7}
\end{split}
```

Consequently ||E||F=o(sqrt(n)) changes the normalized signed expected
energy by o(1). In particular this holds when the *difference* E has
uniformly bounded operator norm and rank o(n), since
||E||F<=sqrt(rank(E))||E||op. The affine reference must itself be PSD;
the numerical condition |s|<=1/6 alone does not imply admissibility.
The near-minimizer application already supplies admissibility and
s=O(n^(-1/2)), so the scalar cutoff holds eventually.

Thus the signed-mean exclusion in (4) survives these perturbations.
No bound on expected absolute energy or on arbitrary anisotropic
covariances follows from this corollary.
