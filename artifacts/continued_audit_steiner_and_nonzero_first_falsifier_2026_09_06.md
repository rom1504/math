# Steiner feedback audit and a bounded nonzero-first covariance falsifier

Date: 2026-09-06. The director's construction and executable are
`continued_director_steiner_feedback_test_2026_09_06.md` and
`computations/continued_director_steiner_signing_2026_09_06.py`.
This note independently verifies their integer construction, closes the
hard-threshold consequence using the audited zero-first theorem, and gives
an actual-signing falsifier of its tempting nonzero-first extension.

## 1. Exact construction audit

The nonzero vectors of F_2^k form v=2^k-1 points; the triples
`{a,b,a+b}` partition the unordered pairs. Every point lies in
`r=(v-1)/2` triples, and there are `b=vr/3` triples. A Sylvester matrix
of order r+1 has first row all ones. Insert its r other rows into the
incident triples of each point, using each row once, to form the b-by-N
matrix E, where `N=v(r+1)`.

Two different columns associated with the same point have inner product
-1, because deleting the all-ones row from a Hadamard leaves off-diagonal
Gram entries -1. Columns associated with different points overlap in
exactly one triple, so their inner product is a sign. Every column has
squared norm r. Distinct rows of E either have disjoint point support or
share one point, where two distinct Hadamard rows are orthogonal. Each
row has squared norm `3(r+1)`. Thus

```math
EE^{\mathsf T}=3(r+1)I,\qquad A=E^{\mathsf T}E-rI
```

is an exact hollow signing. Its Gram identity gives

```math
A^2=(r+3)A+r(2r+3)I=(r+3)A+(N-1)I. \tag{1}
```

The nonzero eigenvalue of `E^T E` is `3(r+1)` with multiplicity b;
the rest are zero. Hence A has eigenvalues `2r+3` and `-r`, with the
claimed multiplicities. This also checks `N-1=r(2r+3)` algebraically.

The executable was read completely and independently replayed at
`N=6,28,120,496`; every incidence, row-Gram, sign-entry, and quadratic
identity assertion passed using integer arithmetic. No output option was
provided, so no stored report was overwritten.

Set `B=A/sqrt(N-1)`, `Q=B²`, and `gamma=(r+3)/sqrt(N-1)`. Then

```math
Q=I+\gamma B,\quad B^3=\gamma I+(1+\gamma^2)B,
\quad\gamma\to1/\sqrt2,\quad\|B\|_{\rm op}\to\sqrt2. \tag{2}
```

In particular `Tr(B³)/N=gamma` exactly. This is not a near-minimizer
assertion.

## 2. The zero-first hard-threshold consequence is now unconditional

Let f be any fixed nonzero bounded odd Gaussian-a.e.-continuous scalar
function with zero first Gaussian coefficient. Write
`tau²=sum_(p>=3 odd) f_p²>0`. Entrywise sign flatness and (2) give

```math
Q^{\circ p}=I+
 \frac{\gamma^p}{(N-1)^{(p-1)/2}}B\qquad(p\ge3\text{ odd}).
```

Therefore the entire INFINITE Hermite sum is exactly

```math
R_f=\tau^2I+\beta_N B,\qquad
\beta_N=\sum_{p\ge3\ {m odd}}
 f_p^2\frac{\gamma^p}{(N-1)^{(p-1)/2}}=O(\tau^2/N). \tag{3}
```

For all sufficiently large N, the squared ratio
`gamma²/(N-1)<1`; a geometric domination by the p=3 factor justifies
the bound and the infinite sum. All coefficients are nonnegative.
Consequently

```math
T=BR_fB=\tau^2Q+\beta_NB^3,\qquad
T_{ii}=\tau^2+\beta_N\gamma\ge\tau^2>0. \tag{4}
```

Thus the positive variance floor required for hard thresholding holds
uniformly, rather than being an assumption left over from the general
theorem. The normalized Gaussian correlations have the form

```math
\frac{T_{ij}}{T_{ii}}=\alpha_N B_{ij}\quad(i\ne j),\qquad
\alpha_N=\frac{\tau^2\gamma+\beta_N(1+\gamma^2)}
                  {\tau^2+\beta_N\gamma}\to1/\sqrt2. \tag{5}
```

The audited normalized-nuclear covariance theorem, followed by fixed-width
smooth threshold approximation and the floor (4), now proves for the
ACTUAL Boolean response `C=sign(Bf(BS))`

```math
\frac{\mathbb E C^{\mathsf T}BC}{2N}
\longrightarrow\frac1{\pi\sqrt2}. \tag{6}
```

Indeed the comparison Gaussian sign energy is EXACTLY
`sqrt(N-1) arcsin(alpha_N/sqrt(N-1))/pi`.
For an old bounded even mask H with mean mu, its Gaussian covariance has
off-diagonal expansion `mu²+O(Q_ij²)`. The extra terms vanish in the
normalized energy trace, yielding `mu²/(pi sqrt(2))`.

This is an actual expected self-energy, not merely a lower bound or a
one-root prediction. It does not require `|f|+H<=1` unless the sum of
the old and new responses is also being used as a cube endpoint.

The proposed example
`f(z)=c[sin(z)-(exp(3/2)/2)sin(2z)]`, scaled to be bounded by one,
has zero first chaos because `E[N sin(tN)]=t exp(-t²/2)`.
It is not identically zero and therefore has positive Gaussian variance.

## 3. Exact return law already defeats Gaussianization

Use the same actual signings and a fresh Boolean seed S. Put `G=BS`.
Hollowness makes G_i independent of S_i exactly. Equation (2) gives

```math
(QS)_i=S_i+\gamma G_i. \tag{7}
```

Its variance is exactly `1+gamma²`, but it is not a centered Gaussian
with that variance. In fact, for `t=pi/4`, own-spin averaging gives

```math
\mathbb E\sin^2(t(QS)_i)=\frac12 \tag{8}
```

for every coordinate and EVERY matrix order: the expectation of
`cos((pi/2)(S_i+gamma G_i))` is zero. A variance-matched Gaussian W_i
instead satisfies

```math
\mathbb E\sin^2(tW_i)
=\frac12\left(1-e^{-\pi^2(1+\gamma^2)/8}\right).
```

The limiting diagonal covariance gap is

```math
g_*:=\frac12e^{-3\pi^2/16}=0.078575289984268\ldots>0. \tag{9}
```

No unidentified cross-root covariance is needed to expose this failure:
it is already on the diagonal.

## 4. The failure survives bounded functions and pointwise feasibility

The linear input in Section 3 is only an intermediate comparison. Define
the fixed continuous odd clipping function

```math
h(z)=\max(-4,\min(z,4)),\qquad f(z)=h(z)/8,
\qquad H(z)=1/2,\qquad\psi(y)=\sin(2\pi y). \tag{10}
```

Then f is bounded by 1/2, H is bounded even, psi is smooth bounded odd,
and `|f|+H<=1` pointwise. The actual feasible mixed response is

```math
C_i=\tfrac12\sin\left(\tfrac\pi4[B h(G)]_i\right).
```

Let `delta=||h(N)-N||2`. Direct integration gives

```math
\delta^2=2[(1+16)\overline\Phi(4)-4\phi(4)],
\qquad\delta=0.0024860442890247\ldots. \tag{11}
```

For a proof without decimal reliance, replacing `exp(-u²/2)` by one in
the shifted tail integral gives `delta²<=phi(4)/16`, and hence
`delta<1/256`. The latter follows already from `e²>7` and `sqrt(2pi)>2`.

For the actual sign fields, all signed row sums have the same normalized
Rademacher-sum law. Consequently

```math
\limsup\left(\frac1N\mathbb E\|B[h(G)-G]\|^2\right)^{1/2}
\le L\delta, \tag{12}
```

with any eventual operator cap L; use L=2 for this family. The bounded
function `sin²(ty)` is t-Lipschitz, so (12) changes its averaged expectation
by at most `tL delta`.

Now consider the tempting generalized covariance prescription obtained by
adding the FIRST Hermite level to the zero-first theorem:

```math
T_h=B\left[\sum_{p\ge1\ {m odd}}h_p^2Q^{\circ p}\right]B,
\qquad Z_h\sim N(0,T_h). \tag{13}
```

Couple Z_h with the linear comparison Gaussian `W~N(0,Q²)` using
independent Gaussian channels of covariance `B Q^(circ p)B`. Their
normalized mean-square difference is at most `L² delta²`, because the
coefficient difference is the Gaussian L2 error `h-id`. Thus replacing
W by Z_h costs at most another `tL delta` in the averaged sine-square test.

Combining (8)–(13), the actual and prescribed Gaussian diagonal
covariances for the bounded feasible C differ by at least

```math
\liminf\frac1N\sum_i
\left(\mathbb E C_i^2-\tfrac14\mathbb E\sin^2(tZ_{h,i})\right)
\ge\frac14(g_*-\pi\delta)>0. \tag{14}
```

Numerically the right side is greater than `.01769`. A coarse rational
certificate is enough: `pi²<10` and `e²<8` give `g_*>1/16`, while
`pi delta<1/64`; hence (14) is strictly greater than `3/256`.

For a symmetric matrix, its nuclear norm is at least the absolute trace.
Therefore (14) is a positive normalized-nuclear covariance discrepancy.
This falsifies the nonzero-first extension even for fixed bounded input,
fixed smooth bounded response, an even bounded mask, feasible endpoints,
and an explicit bounded-op dense hollow-signing family.

The failure is a failure of the proposed transfer law, not a contradiction
to the zero-first theorem or a nonconvergence result for the minimizers.

The executable
`computations/continued_audit_nonzero_first_gaussian_failure_2026_09_06.py`
replays the exact signing identities, the finite Rademacher characteristic
function, binomial clipping errors, and the Gaussian comparison bounds.
It uses no random sampling or solver and writes no stored report.

## 5. Optional sign-preserving response variant

The falsifier can also satisfy the gain-side condition `y psi(y)>=0`.
Define `rho(y)=sin(pi y/4)` for `|y|<=4` and zero outside, and replace
psi in (10) by `rho(8y)`. This function is odd, bounded, continuous,
Lipschitz, and sign-preserving. Its square is still pi/4-Lipschitz before
the rescaling. Thus the two clipping comparison errors above are unchanged.

For the exact linear return (7), the only loss from using rho instead of
the full sine occurs when `|S_i+gamma G_i|>4`. The flat Rademacher row sum
G_i is unit-subgaussian, so this probability is at most
`2exp(-9/(2gamma²))`, tending to `2exp(-9)`. On the Gaussian comparison
side truncating the squared sine can only decrease its expectation.
Therefore the final normalized trace gap is at least

```math
\frac14\left(g_*-2e^{-9}-\pi\delta\right)>\frac{95}{8192}>0.
```

The rational inequality uses the bounds already given together with
`2e^-9<1/2048`. Hence even all the sign/feasibility conditions of the
endpoint-gain assertion do not repair the NONZERO-FIRST covariance claim.
