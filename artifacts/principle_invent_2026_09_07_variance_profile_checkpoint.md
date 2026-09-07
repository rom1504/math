# Bounded variance-profile checkpoint

2026-09-07. Status: exact normalizations, primary-source inspection, and
explicit finite tests. Favorable flatification remains open.

## 1. Square bipartite comparison

For a full sign n-by-n C and N=2n, set

```math
W=\sqrt{(2n-1)/n}\begin{pmatrix}0&C\\C^T&0\end{pmatrix}.
```

Every row has squared norm N-1, maxentry<sqrt(2), and

```math
\frac{Q(W)}{N^{3/2}}
=\frac{\sqrt{1-1/(2n)}}2\frac{\beta(C)}{n^{3/2}}.
```

This calculation was already in the archive's favorable-flatification
scope note. The separately audited marked-response proof gives the square
full-sign lower `liminf beta(C)/n^(3/2)>=2c_*=.8666442233281615...`.
Therefore an actual full-sign square bipartite profile cannot provide a
counterexample below the current original lower endpoint c_*.

Primary sources inspected directly:

- Pellegrino--Raposo, *Upper bounds for the constants of Bennett's inequality
  and the Gale--Berlekamp switching game*, arXiv:2111.00445v3,
  https://arxiv.org/html/2111.00445 . Its introduction establishes an
  all-large-order upper 1+epsilon through Hadamard orders, not a sub-.8666
  asymptotic family.
- Dumitrescu, *Geometric Variants of the Gale--Berlekamp Switching Game*,
  arXiv:2412.16994v1, https://arxiv.org/html/2412.16994v1 . It records the
  classical sqrt(2/pi) asymptotic lower bound. This does not supersede the
  stronger archived marked-response theorem.

No literature theorem identifying the optimized bilinear constant, or
giving the desired weighted counterexample, was found in this search.

## 2. Concrete bounded-amplitude tensor-rotation candidate

Let H_s be the Sylvester Hadamard matrix of order s, n=2s, and

```math
C_\theta=\sqrt2\begin{pmatrix}
\cos\theta\,H_s&\sin\theta\,H_s\\
-\sin\theta\,H_s&\cos\theta\,H_s
\end{pmatrix}.
```

Every row and column has squared norm n; entries are bounded by sqrt(2).
Bipartite dilation therefore produces an actual hollow row-regular weighted
profile with bounded amplitudes (at most 2 after the dilation scaling).

For theta=pi/8, exhaustive enumeration at n=2,4,8,16 gives the normalized
bilinear value cos(pi/8)=.9238795325112867. This suggested a possible
constant-gap orthogonal family. It is NOT an asymptotic upper certificate.
An explicit order-64 witness falsifies that tempting extrapolation.

Writing x=(x_1,x_2), y=(y_1,y_2), define the exact integers

```math
A=x_1^TH_s y_1+x_2^TH_s y_2,\qquad
B=x_1^TH_s y_2-x_2^TH_s y_1.
```

The preserved order-64 witness has A=320 and B=116, hence value

```math
\frac{\sqrt2[320\cos(\pi/8)+116\sin(\pi/8)]}{64^{3/2}}
=.9392165944871081\ldots>\cos(\pi/8).
```

The witness and its y vector are printed by the reproducible script
`computations/principle_invent_2026_09_07_rotation_profile_check.py`.
It starts from the fixed NumPy seed 20260907 and retains the maximizing
iterate; multiplication of the returned vectors checks the integer pair
without trusting the optimizer's global quality. Other returned lower
witnesses were `(n,A,B)=(32,128,0),(128,928,256),(256,2456,1072)`.

Tensoring this explicit witness with a regular Hadamard H_(4^r) preserves
its normalized value, so the violation of the cos(pi/8) conjecture is
scalable on the corresponding infinite order subsequence. This is a scoped
failure of one candidate, not a weighted-profile universality theorem.

## 3. Positive theorem target selected instead

The new substantive target is the fixed-finite-species marked lower
extension in `principle_invent_2026_09_07_finite_species_marked_lower.md`.
Its proposed new mechanism is replacement of the exceptional kernels by
finite Schur products of the PSD subcorrelations `Q_s=B P_s B`.
It would establish c_* uniformly for every fixed block-constant row-regular
profile, including all two-block profiles. It still would not identify
their optimized values with the original full-sign optimum.
