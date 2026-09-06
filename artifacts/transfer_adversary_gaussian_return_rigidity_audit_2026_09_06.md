# Independent audit: inhomogeneous Gaussian nonlinear return

Date: 2026-09-06. Audited source:
`transfer_director_inhomogeneous_gaussian_return_rigidity_2026_09_06.md`.
Verdict: PASS as a finite Gaussian-kernel theorem. The actual Boolean
mixed-covariance comparison is a separate unresolved input at this audit.

## Signed coefficients and the exact n-dimensional rank bound

For `T_ij=<v_i,v_j>`, odd `k=2r-1`, and ANY real diagonal entries d_i,
form the symmetric operator

```math
S=sum_i d_i (v_i^(tensor r))(v_i^(tensor r))^T.
```

Even with mixed-sign d_i, its rank is at most n: its range is contained
in the span of the n displayed vectors. Therefore

```math
tr(T D T^(circ k) D)=||S||_F^2
>= (tr S)^2/n = (sum_i d_i T_ii^r)^2/n.
```

The inequality follows from Cauchy--Schwarz on the nonzero real
eigenvalues; S need NOT be PSD. This removes the need for a common
sign of higher row coefficients. With
`d_i=alpha_ik/q_i^(k/2)`, the trace on the right is exactly
`sum_i sqrt(q_i) alpha_ik`, giving n times the squared kth Hermite
coefficient of the averaged response. No ambient tensor dimension
replaces the crucial rank n.

The Schur matrix `D T^(circ k) D` IS PSD, despite signed D, by
congruence of a Gram matrix. Thus `Q>=T/C` permits the proposed trace
comparison, and the operator bound gives
`tr(Q D T^(circ k) D)<=||Q||op sum_(q_i>0) alpha_ik^2<=||Q||op n`.
Coordinates q_i=0 have zero Gram rows. One should either omit their
coefficients from this trace sum or define their ignored coefficients
as zero; including an arbitrary nonzero alpha from an ignored row would
not be an equality. Parseval gives a finite total trace, legitimizing
the infinite sum before the subsequent fixed-degree truncation.

## Nonlinear mass and a uniform finite degree

Let the averaged response g be odd, nondecreasing, bounded in absolute
value by M0, and have first orthonormal Hermite coefficient a>=c>0.
Pointwise,

```math
|g(t)-at| >= (c|t|-M0)_+.
```

Since a is the actual orthogonal first coefficient and g is odd,
Parseval identifies its squared distance from a times t with the
sum of its odd Hermite masses of degrees at least three. This gives
the claimed strictly positive
`delta=E(c|G|-M0)_+^2`. Its displayed closed form follows by the first
two Gaussian tail-moment identities and has the correct factor two.

For completeness, the monotone-tail estimate has two elementary
independent reconstructions. A bounded increasing function has a
positive Stieltjes increment measure of mass at most 2M0. Weighted
Cauchy--Schwarz bounds its squared Gaussian noise difference by this
mass squared times maximum threshold disagreement. For Gaussian
correlation rho>=0, that disagreement is at most `arccos(rho)/pi`.
The source's derivative test for the two-coordinate Gaussian CDF is
correct. Alternatively, realize the two coordinates as projections of
an isotropic planar Gaussian onto axes separated by angle theta.
Conditional on radius, each threshold event is a circular arc; its
symmetric difference with a theta-rotation has normalized length at
most theta/pi. This also handles all thresholds uniformly.

Hence `E(g(X)-g(Y))^2<=4M0^2 arccos(rho)/pi`. Mehler's identity
has the factor two, giving exactly the source's tail inequality for
`rho=1-1/K`:

```math
sum_(j>K) a_j^2
<=2M0^2 arccos(1-1/K)/[pi(1-(1-1/K)^(K+1))]
<=3M0^2/sqrt(K).
```

The last estimate uses `arccos(1-1/K)<=pi/sqrt(2K)` and denominator
greater than one half, leaving constant `2sqrt(2)<3`. Thus an integer
`K0>=max(2,36M0^4/delta^2)` leaves at least delta/2 mass in a fixed
finite set of odd degrees at least three. Pigeonhole gives one trace
at least `n delta/(2K0)`, and comparison with Q gives the additional
factor 1/C. All constants and normalization factors check.

## Actual-feedback boundary

The retained degree may depend on the finite input but lies in a set
fixed before the order limit. A Boolean application needs its mixed
covariance error and probe second moment uniformly over that finite
set. If the displayed mixed identity and PSD comparison hold, the
claimed full-row correlation `n delta/(2 C0 K0)+o(n)` follows.
For a general nonsymmetric linear map use `Q=B^T B`; writing B squared
requires symmetric B, as in the intended signing setting.

This algebra does not establish that identity for the literal old
coherent fields, nor the probe's joint Gaussian law with those fields.
Marginal Gaussian convergence is insufficient. A Frobenius-small
covariance replacement cannot be applied through an arbitrary Boolean
sign operation without a separate comparison. Actual functions must
also satisfy the required oddness, monotonicity and uniform first-mass
and amplitude hypotheses; the finite lemma does not infer them from
names of response variables. Finally, a positive constant at each fixed
operator cap is not automatically large enough to pay an unrestricted
spectral-deletion loss.

## Finite algebra replay

`computations/transfer_adversary_gaussian_return_rank_2026_09_06.py`
passed 640 exact integer tensor-rank and PSD-trace tests, including
111 mixed-sign diagonal inputs and 40 inputs with a zero Gram diagonal.
The result JSON has the same basename under `computations/results/`.
These are exact identity regression tests, not numerical evidence for
the still separate Boolean mixed-field closure.
