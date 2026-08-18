# Adversarial audit: low-transport feature ceiling

**Object audited:**
[`../drafts/actual_child_low_transport_feature_ceiling.md`](../drafts/actual_child_low_transport_feature_ceiling.md).

**Verdict:** **PASS.**  LT.1's fibre coupling, LT.2's coordinate marginal,
LT.3's Walsh-pivot argument, all constants, and the three finite comparison
values are correct.  LT.4 is also a valid generic sharpness example; it is
properly separated from the actual-child theorem.

## LT.1

Conditioned on a feature vector, the fair bridge law is the product of the
uniform row fibres.  A `W_infinity` coupling of two `i`th fibres with Hamming
distance at most `delta_i`, tensored with identity couplings on the other
rows, is therefore a valid coupling of the full conditional laws.  The
pointwise pressure oscillation is at most `2t delta_i`, so the inverse weight
ratio lies in `exp(+-2 lambda t delta_i)`.  Taking expectations preserves
this ratio and gives coordinate range `2 lambda t delta_i` for `g`.

Hoeffding's bounded-difference constant is consequently

```math
 {1\over8}\sum_i(2\lambda t\delta_i)^2
 ={\lambda^2t^2\over2}\sum_i\delta_i^2.
```

Since `E_nu exp(g)=1`, the centered log-MGF is `-E_nu g=D(nu||Q)`.
Finally `nu` is a feasible row product, so both inequalities in LT.5 have
the stated direction.

## LT.2

For the selected bit set `S`, marginalizing the positive inverse weight
preserves a `2 lambda t` log-density oscillation per retained coordinate:
the uniform fibres before and after flipping a selected bit are related by
the flip bijection.  Applying the same product-cube inequality gives
`D(U_S||q_S)<=lambda^2t^2|S|/2`.  The declared feature is a rowwise function
of `B_S`, so KL data processing sends this to `D(nu||Q)` in the required
direction.  No independence under `q_S` is assumed.

## LT.3

A rank-`r_i` linear parity map has `r_i` pivot columns spanning its image,
with `r_i<=k_i`.  Every attainable syndrome difference is a sum of at most
`r_i` pivots.  Translation by the corresponding Hamming vector is a
measure-preserving bijection between the two uniform affine fibres, giving
an explicit `W_infinity` coupling of radius at most `k_i`.  This proves
LT.12 even when the recorded parities overlap or collectively inspect all
row coordinates.

## Constants and finite table

Substitution in LT.10 gives

```text
N=8:  lambda^2 * (16/8) * 4 / 2 = 115.8681822962743
N=9:  1^2       * ( 4/9) * 4 / 2 = 8/9
N=9:  1^2       * (16/9) * 4 / 2 = 32/9.
```

These match the table.  The selected masks are one-coordinate Walsh
characters, hence `K=4` in each case.

For LT.4, a bit flip changes `F_r` by at most `2/sqrt(r)`, hence the pair
log-density by at most `2 gamma/sqrt(r)`.  The bounded-potential row density
estimate stated in the proof is conservative but valid.  The majority-sign
pushforward has correlation (LT.17), and CA.3 applies because
`4 atanh(rho_r)<4`.  Its one-pair reverse product value is exactly
`-(1/2)log(1-rho_r^2)`; additivity and data processing give the extensive
claim.  The all-plus word proves the asserted essential-supremum transport
lower bound.
