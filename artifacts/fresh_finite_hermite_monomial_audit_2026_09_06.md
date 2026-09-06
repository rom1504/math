# A finite-degree nonlinear Hermite mass certificate from one monomial

Date: 2026-09-06. Independent audit of the director's short monomial
argument. The argument passes and needs neither Gaussian regularity of
the response nor a Hermite-tail estimate.

## 1. Statement

Let `G` be standard Gaussian and let `F` be any jointly defined random
variable with `|F|<=1`. It need not be a function of `G` alone. Suppose

```
a=E[G F]>=a0>0.
```

Write `h_r=He_r/sqrt(r!)` for the orthonormal probabilists' Hermites.
For every odd integer `D>=3`, put

```
s_D^2 = sum_(3<=r<=D, r odd) |E[F h_r(G)]|^2,
b_D=D!!,
V_D=(2D-1)!!-(D!!)^2 >0.
```

Then

```
s_D >= [a0 D!! - E|G|^D]_+ / sqrt(V_D).                (1)
```

In particular,

```
E|G|^D / D!! <= 1/sqrt(D),                             (2)
```

so any odd `D>a0^(-2)` gives a strictly positive uniform lower bound.
No oddness assumption on `F` is needed; only its odd Hermite
coefficients appear in the argument. An odd response is of course a
special case.

## 2. Proof

The polynomial `p_D(G)=G^D-b_D G` is orthogonal to `G`, because
`E G^(D+1)=D!!`. It is odd and of degree `D`, so it lies exactly in
the span of `h_3,h_5,...,h_D`. Its squared norm is

```
||p_D||_2^2 = E G^(2D)-b_D^2 = V_D.
```

Boundedness of `F` gives

```
E[F p_D] = E[F G^D]-b_D a
          <= E|G|^D-b_D a0.
```

When the right side is negative, Cauchy--Schwarz in the indicated
finite Hermite subspace gives
`b_D a0-E|G|^D <= |E[F p_D]| <= s_D sqrt(V_D)`.
If it is nonnegative, the positive-part assertion is trivial.
This proves (1), also for arbitrary jointly defined bounded `F` by
conditional expectation onto `G` if desired.

The requested odd-moment ratio needs only another Cauchy--Schwarz
application:

```
(E|G|^D)^2
 <= E|G|^(D-1) E|G|^(D+1)
 = (D-2)!! D!! = (D!!)^2/D.
```

This proves (2) without importing a Gamma-function estimate.

## 3. A completely elementary explicit floor

For odd `D>=1` let `R_D=(2D-1)!!/(D!!)^2`. Direct cancellation gives

```
R_1=1,
R_(D+2)/R_D = [(2D+3)(2D+1)]/(D+2)^2 <4.
```

Thus `R_D<=2^(D-1)`. Equations (1)--(2) imply

```
s_D >= (a0-1/sqrt(D))_+ / sqrt(R_D-1)
     >= (a0-1/sqrt(D))_+ * 2^(-(D-1)/2).
```

Choosing an odd `D>=4/a0^2` gives the convenient uniform certificate

```
s_D >= a0 * 2^(-(D+1)/2) >0.                           (3)
```

For the campaign's `a0=1/200`, `D=160001` is a safe explicit choice
in (3); `D=40001` already yields positivity in (1)--(2). The constants
are intentionally loose. The improvement is the polynomial degree
bound `D=O(a0^(-2))`, not a claim that this very large degree is
computationally practical.

## 4. Exact role in the joint response theorem

The lower bound is on the entire finite coefficient mass `s_D`.
The existing joint direction chooses
`h=sum_(r>=3 odd) (E[F h_r(G)]/s_D) h_r`, through degree `D`.
Then its weighted energy coefficients are nonnegative squares, giving
the convex Schur mixture needed for variance normalization.

One must not instead use `p_D` as the response direction merely
because its unweighted correlation has a known sign: the different
Hermite transport variances need not be equal, and mixed signed
coefficients would lose the Schur-mixture lower bound. The monomial
is a separating test proving finite mass; the actual jointly
normalized response remains coefficient-aligned as before.

This lemma removes the tail-regularity/compactness step for finite
degree selection. It does not supply local cube slack, force a
nonzero edge coefficient in a new response family, or prove original
minimum convergence.

## 5. Independent general-mask projection check

For the particular paired family with `0<=H<=1`, `|F|<=1-H`,
`W=U H`, and `G0=U1`, where `U` is the old Hilbert isometry, let
`p=EH`, `u=EH^2`, `J=E[F W]`, and `a=E[F G0]`. The isometry gives
`<W,G0>=p`, `||W||_2^2=u`; feasibility gives
`||F||_2^2<=1-2p+u`. Orthogonal projection off `G0` therefore yields

```
J <= p a + sqrt(u-p^2) sqrt(EF^2-a^2).
```

Since `p^2<=u<=p`, the remaining product is at most
`sqrt(p)(1-p)<=2/(3sqrt(3))`. If `J>=43/100`, this first proves
`a>0`; then `p a<=a` proves

```
a >= 43/100-2/(3sqrt(3)) > 9/200.
```

For the last rational comparison,
`(77/200)^2-4/27=83/1080000>0`. Thus the general paired family,
whenever it satisfies these hypotheses, admits the safe monomial
degree `D=1977` in (3); `D=495` already gives positivity. This is an
independent check of the director/algebra projection branch, not a
solution of its remaining local-slack problem.

The sharper director/algebra projection in
`fresh_general_mask_nonlinear_direction_2026_09_06.md` gives
`a>=9923/202500>49/1000`, reducing the common positive-mass degree
to `417`. That entire note, including its denominator bound and its
low-value creation-space slack counterexample, was independently read
and checked during this audit. The simpler degree choices above are
valid general estimates, not the best campaign degree bound.
