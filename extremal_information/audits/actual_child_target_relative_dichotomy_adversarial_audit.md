# Adversarial audit: actual-child target-relative dichotomy

**Object audited:**
[`drafts/actual_child_target_relative_dichotomy.md`](../drafts/actual_child_target_relative_dichotomy.md)

**Verdict:** **PASS after one asymptotic repair.**  TR.1--TR.10 have the
correct signs, constants, and orientation normalization.  The product-basin
versus reverse-projection dichotomy is rigorous conditional on fixed-tilt
target reachability.  I repaired the choice of the basin slack when the
target is reached only up to `o(N)`: one must also absorb a possibly negative
`h=T-L_0`.  The theorem is a strict target-relative diagnostic, but target
reachability and a low-information way to decide `Delta` remain separate
obligations.

## 1. Fixed-orientation zero-bridge floor

For a fixed signing orientation `epsilon`, division by the zero-bridge
partition function gives TR.6.  The involution `x->-x` preserves both
quadratic child Hamiltonians and the fixed orientation, while sending
`Q=tau xy^T` to `-Q`.  Therefore the law of `Q` under the zero-bridge Gibbs
measure is central and

```math
E e^{t\langle B,Q\rangle}=E\cosh(t\langle B,Q\rangle)\ge1.
```

No average over the two orientations and no factor of two is missing:
`L_epsilon(B)` is the pressure of the actual block signing
`(A,epsilon D,B)`.  If one later packages both orientations into a joint
disorder law, its prior contributes at most `log 2`, but that is a different
object from this fixed-sector theorem.

Since `L>=L_0`, also `V_lambda>=L_0`.  Thus exact target reach
`V_lambda<=T` implies `h>=0`, while approximate reach
`V_lambda<=T+e_N`, `e_N=o(N)`, implies only `(-h)_+=o(N)`.

## 2. Product Renyi transfer and basin mass

Every global row-product minimizer is a coordinate minimizer.  Entropy makes
the coordinate minimizer the positive Gibbs row in AC.17, and its one-bit
oscillation is `2lambda t`.  Hence TR.7 is exactly AC.18.  Renyi-two
divergence is additive under products, giving

```math
D_2(p^*||U_B)\le m n\lambda^2t^2
=\lambda^2\beta^2mn/N.
```

There is no hidden extra factor of two.

The Markov step is valid because `L-L_0>=0` and

```math
E_{p^*}L\le V_lambda^row=T+Delta.
```

For `a>Delta` and `h+a>0`, it yields exactly
`p^*(G_a)>=(a-Delta)/(h+a)`.  Cauchy--Schwarz with
`f=dp^*/dU` then gives

```math
U(G_a)\ge p^*(G_a)^2/E_Uf^2,
```

which is TR.9.  Thus the lower mass is genuinely `exp[-O(N)]` at comparable
splits and fixed `beta,lambda`; it is not a claim of polynomial mass.

## 3. Directed projection branch

TR.10 is an immediate exact consequence of

```math
I_lambda^leftarrow
=lambda(V_lambda^row-V_lambda).
```

If `V_lambda<=T+e`, then
`I_lambda^leftarrow>=lambda(Delta-e)`.  Forward total correlation cannot be
substituted here.  No target or orientation constant is lost.

## 4. The repaired asymptotic choice

The original example `a_N=2(Delta_N)_++sqrt(N)` need not satisfy
`h_N+a_N>0` under approximate reach: one may have
`h_N=-N/log N` while still `h_N=-o(N)`.  The draft now uses

```math
a_N=2(Delta_N)_++2(-h_N)_++sqrt(N).
```

Because `V_lambda>=L_0`, TR.11 gives `(-h_N)_+=o(N)`, so this `a_N` is
`o(N)`, has positive denominator, and leaves numerator at least `sqrt(N)`.
With `h_N=O_beta(N)`, the prefactor in TR.9 is at least polynomial, hence
does not alter the `exp[-O(N)]` scale.

For a Hammersley-summable defect `O(N^(1-delta))`, both
`(Delta_N)_+` **and** `(-h_N)_+` need that rate.  Exact reachability makes
the second condition automatic.  The draft now states this qualification.

## 5. How much the SML narrows

Conditional on fixed-`lambda` target reachability, the theorem is a real
strict narrowing: the entire product side is reduced to the scalar target
excess

```math
Delta_N=V_lambda^row-T_N.
```

If its positive part is sublinear, an explicit row-product witness plus
TR.9 supplies a sublinear-pressure, linearly rare basin.  If it is linear,
TR.10 certifies the exact extensive directed resource.  This is more than a
renaming of the uniform-to-escort gain split.

It is not yet an information-theoretic or Level-5-to-6 reduction.  Three
obligations remain distinct:

1. prove that one fixed `lambda` reaches the target to the required rate;
2. decide the scale of `Delta_N` from a statistic smaller than the full
   bridge response table (the AC.17 oracle reconstructs that table on point
   masses);
3. for convergence-scale use of branch 1, strengthen bare `o(N)` to the
   summable defect required by the archived recurrence criterion.

Branch 2 is a genuine structural theorem about the actual escort, but by
itself it is a diagnosis rather than a composition mechanism.
