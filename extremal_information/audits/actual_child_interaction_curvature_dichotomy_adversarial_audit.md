# Adversarial audit of the interaction-curvature dichotomy

**Object audited.**
[`../drafts/actual_child_interaction_curvature_dichotomy.md`](../drafts/actual_child_interaction_curvature_dichotomy.md),
Theorems IC.1--IC.2 and the proposed smallest missing lemma.

**Verdict: PASS after a scope correction.**  The escort interpolation,
normalizations, direction of KL, cumulant signs, variance identity, constants,
bit-flip bound, and conditional Renyi-two estimate are all correct.  The
original draft over-described a linear canonical mismatch as quantitatively
irreducible row dependence and called the uniform-curvature target weaker in
data requirements without separating syntactic information from mathematical
strength.  The source now states the exact scope: the weighted variance
integral is equivalent to the canonical error `J=o(N)`; the supremum bound is
stronger; and the canonical error only upper-bounds the distance to the best
row product.

## 1. Endpoint and normalization reconstruction

Write

```math
h=\log p-\sum_i\log p_i,
\qquad
{dr\over dU}={\prod_i p_i^{-\lambda}\over\prod_i Z_i},
\quad Z_i=E_{U_i}p_i^{-\lambda}.
```

Multiplication by `exp(-s h)` gives

```math
{dq_s\over dU}
\propto
p^{-s}\prod_i p_i^{-(\lambda-s)}.                  \tag{AIC.1}
```

Thus `q_0=r`.  At `s=lambda`, all row-marginal factors cancel and
`dq_lambda/dU` is proportional to `p^{-lambda}`, so `q_lambda=q`.  There is
no omitted factor `prod_i Z_i`: it is absorbed into the normalizer in
(AIC.1).

The centered and uncentered definitions of the path agree because
`exp(s E_r h)` is constant in the bridge.  Hence the path is well defined
whenever the channel amplitude is finite; all likelihoods are then strictly
positive.

## 2. Direction of KL and cumulant signs

Let

```math
\psi(s)=\log E_r e^{-s h},
\qquad m=E_rh,
\qquad K(s)=\psi(s)+sm.
```

Finite differentiation gives

```math
\psi'(s)=-E_{q_s}h,
\qquad
\psi''(s)=\operatorname{Var}_{q_s}(h),              \tag{AIC.2}
```

so `K(0)=K'(0)=0` and `K''=Var_(q_s)(h)`.  More
decisively, for every `s`, not only at the endpoint,

```math
\begin{aligned}
D(r\Vert q_s)
 &=E_r\log {dr\over dq_s}\\
 &=sE_rh+\log E_r e^{-sh}=K(s).                     \tag{AIC.3}
\end{aligned}
```

Thus the direction in IC.7 is `D(r||q)`, and both signs in the centered
lower-tail cumulant are correct.  Replacing it by `D(q||r)` would give a
different expression and is not licensed.

Twice integrating (AIC.2) yields

```math
K(\lambda)
=\int_0^\lambda(\lambda-s)
  \operatorname{Var}_{q_s}(h)\,ds.                 \tag{AIC.4}
```

The weight has mass `lambda^2/2`.  Therefore

```math
K(\lambda)\ge\eta N
\quad\Longrightarrow\quad
\max_{0\le s\le\lambda}\operatorname{Var}_{q_s}(h)
\ge {2\eta N\over\lambda^2},                       \tag{AIC.5}
```

with the constant in IC.9 exactly correct.  For fixed `lambda`, (AIC.4)
also proves the exact equivalence between `J=o(N)` and an `o(N)` weighted
integral.  It does **not** make `sup_s Var_(q_s)(h)=o(N)` necessary: a narrow
curvature spike can have a sublinear weighted integral.

## 3. One-bit oscillation

For `rho=tanh u`, the forward channel likelihood has the planted form

```math
p(B)=E_Q\prod_e(1+\rho B_eQ_e).
```

For a fixed latent word, flipping coordinate `e` multiplies its channel
factor by

```math
{1-\rho B_eQ_e\over1+\rho B_eQ_e}
=e^{-2uB_eQ_e}\in[e^{-2u},e^{2u}].                 \tag{AIC.6}
```

A positive mixture preserves these pointwise ratio bounds, proving
`osc_e(log p)<=2u`.  The same comparison survives averaging all rows except
row `i`, so `osc_e(log p_i)<=2u` for `e` in row `i`; the other row marginals
do not depend on `e`.

By (AIC.1), and crucially because `0<=s<=lambda`, the unnormalized log
density of `q_s` therefore has coordinate oscillation at most

```math
2us+2u(\lambda-s)=2\lambda u.                       \tag{AIC.7}
```

No cancellation was assumed.  Conditioning other coordinates and summing
over hidden coordinates preserve a likelihood-ratio interval, hence preserve
the bound in (AIC.7).

## 4. Conditional Renyi-two constant

Apply the audited cube lemma AC.1 with conditional half-log-odds at most
`lambda u`.  For any retained set `S`, after arbitrary disjoint conditioning
and marginalization,

```math
D_2((q_s)_S\Vert U_S)
\le |S|\log(1+\tanh^2(\lambda u))
\le |S|\lambda^2u^2.                                \tag{AIC.8}
```

Taking `S` to be a row of length `n`, conditioning the row prefix, and
marginalizing later rows gives IC.10.  At `u=beta/sqrt(N)`, its right side
is `lambda^2 beta^2 n/N<=lambda^2 beta^2`, even without a comparable-split
assumption.  This verifies that extensive curvature, if present, cannot be
attributed to an escaping conditional row Renyi-two constant along this
specific path.

## 5. What the dichotomy does and does not prove

The exact consequence is

```text
canonical row-erased product is accurate at o(N),
or some regular hybrid escort has Omega(N) interaction variance.
```

The second alternative is genuine nonproduct structure, but its magnitude
is measured relative to the particular product `r`.  Since

```math
\inf_{a=\otimes_i a_i}D(a\Vert q)\le D(r\Vert q)=\mathcal J,
                                                               \tag{AIC.9}
```

a linear `J` does not imply a linear best-product distance.  Consequently
IC.1 neither decides the requested basin/dependence dichotomy nor proves
that the no-gain phase fails.  The implication in the favorable direction
is valid: `J=o(N)` certifies that the optimal directed product cost is
`o(N)`.

The proposed missing lemma has two distinct readings:

1. `int_0^lambda (lambda-s) Var_(q_s)(log G) ds=o(N)` is an exact scalar
   reformulation of `J=o(N)`.  It removes the variational oracle and declares
   less output, but is not a strict implication-theoretic reduction.
2. `sup_s Var_(q_s)(log G)=o(N)` is a stronger sufficient lemma.  It may be
   more localizable, but cannot be advertised as weaker than the canonical
   problem.

Both formulations use only the child-derived collision--cavity partition
`G`, not a target-order optimizer.  That is a real architectural narrowing.
It becomes a strict mathematical reduction only if optimizer-specific
structure controls the curvature from a state demonstrably coarser than the
full bridge likelihood or child energy landscape.  The scalable rank-one
example CC.18--CC.24 correctly prevents any such conclusion from local
oscillation and bounded row Renyi complexity alone.

## 6. Audit disposition

After the source corrections, IC.1--IC.2 are rigorous and useful as an exact
actual-law diagnostic.  They should be classified as a **RESET in explicit
state/observable**, but not yet as a solved structural branch: the new SML
is a scalar canonical-certificate criterion, not evidence that its variance
is simpler to bound for optimizing children.

## 7. Addendum: tilted average-influence criterion IC.3

**Verdict on IC.3: PASS after a comparability correction.**  The sign of
`tau`, the normalization of `A_s`, the integration constant, and the
converse scale are correct.  The original introductory sentence called the
criterion strictly weaker than the projective-diameter statistic.  That
logical comparison does not hold without a uniform bound on individual row
increments, because `tau` grows exponentially.  The source now calls it a
different, potentially weaker, lower-information certificate and states an
exact sufficient comparison.

### 7.1 Modified log-Sobolev sign

For a product law, entropy tensorization followed by the elementary
comparison form of the modified log-Sobolev inequality says

```math
\operatorname{Ent}(e^{aZ})
\le E e^{aZ}\sum_i
 \tau\{-a(Z-Z_i)\},
\qquad \tau(v)=e^v-v-1,                              \tag{AIC.10}
```

whenever `Z_i` is independent of coordinate `i`.  The sign can be checked
without citation: condition on all other rows, use the entropy variational
bound with comparison constant `exp(aZ_i)`, and write

```math
e^{aZ}\{a(Z-Z_i)-1+e^{-a(Z-Z_i)}\}
=e^{aZ}\tau\{-a(Z-Z_i)\}.                           \tag{AIC.11}
```

In IC.3, `Z=-h`, `Z_i=-bar h_i`, and `a=s`.  Therefore

```math
-s(Z-Z_i)=s(h-\bar h_i),                            \tag{AIC.12}
```

which verifies the positive sign inside `tau` in IC.17.  Replacing it by
`tau(-s(h-bar h_i))` would be an error because `tau` is asymmetric.

### 7.2 Normalization and integration

Let `M(s)=E_r e^{-sh}` and `psi=log M`.  Direct calculation gives

```math
{\operatorname{Ent}_r(e^{-sh})\over M(s)}
=-sE_{q_s}h-\psi(s)
=sK'(s)-K(s).                                       \tag{AIC.13}
```

Dividing IC.17 by `M(s)` changes its right side to the `q_s` expectation.
The definition

```math
\mathcal A_s={1\over s^2}E_{q_s}\sum_i
 \tau\{s(h-\bar h_i)\}                              \tag{AIC.14}
```

therefore gives exactly

```math
sK'(s)-K(s)\le s^2\mathcal A_s.                     \tag{AIC.15}
```

Moreover `tau(sx)/s^2 -> x^2/2` and `q_s -> r` on the finite cube, so IC.14
is the correct continuous endpoint.  Since

```math
\left({K(s)\over s}\right)'
={sK'(s)-K(s)\over s^2}
\le\mathcal A_s
```

and `K(s)/s -> K'(0)=0`, integration from zero to `lambda` yields

```math
K(\lambda)\le\lambda\int_0^\lambda\mathcal A_s\,ds. \tag{AIC.16}
```

Thus the factor in IC.15 is `lambda`, not `1`, `lambda^2`, or
`1/lambda`.

If `K(lambda)>=eta N`, (AIC.16) forces

```math
\int_0^\lambda\mathcal A_s\,ds\ge {\eta N\over\lambda},
\qquad
\sup_{0\le s\le\lambda}\mathcal A_s
\ge {\eta N\over\lambda^2}.                        \tag{AIC.17}
```

The claimed converse scale is therefore exact for this averaging argument
when `lambda` is fixed.

### 7.3 Comparison with projective diameter

Let `delta_i(B_-i)` be the row oscillation from CC.14.  Because `bar h_i`
is an average of the row section under `r_i`, it lies between that section's
minimum and maximum.  Hence

```math
|h-\bar h_i|\le\delta_i(B_{-i}).                    \tag{AIC.18}
```

For `d>=0`, `max_{|x|<=d}tau(sx)=tau(sd)`, which proves IC.19.  If all row
ranges are bounded by a fixed `C`, then uniformly for `0<=s<=lambda`,

```math
{\tau(s\delta_i)\over s^2}
\le c_{C,\lambda}\delta_i^2,                        \tag{AIC.19}
```

so `sum_i sup delta_i^2=o(N)` implies `sup_s A_s=o(N)`.  Without the bound
`C`, a sparse `delta_i` tending to infinity can have sublinear squared sum
but exponentially large `tau(s delta_i)`.  Therefore the average-influence
and projective-diameter hypotheses are not intrinsically ordered, even
though `A_s` uses an average rather than a worst-case range.

IC.3 is nevertheless a legitimate new proof interface: it preserves the
signed complete row replacement under the relevant tilted law and may be
small through probability-weighted cancellation where a supremum range is
large.  Establishing that behavior for optimized children remains the new
mathematical work.
