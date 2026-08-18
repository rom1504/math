# Independent audit: linear-entropy conference conditioning theorem

**Frozen source:**
`extremal_information/drafts/conference_linear_entropy_conditioning_theorem.md`

**SHA-256:**
`486916f58c9d6b714cfc32d257f68283e16a292a75119e4e6bfc70bf4797a163`

**Verdict:** **PASS.**  The two-sided regular input and its hypotheses, net
constants, arbitrary-event probability and `L^1` conclusions, finite-tilt
condition, one-sided KL transport, liminf argument, and the `3/128`
small-temperature limit are correct.  No repair is required.

## 1. Two-sided regular-sector input

The qualifier in (LC.3) is essential and sufficient.  Since
`kappa>beta(3+delta)/sqrt(2)`, the operator ball contains the overwhelming
uniform-bridge sector.  The convex dimension-free Frobenius-Lipschitz
extension therefore has mean `h_beta r+o(r)`.  Two-sided convex
concentration applied to that extension proves (LC.5).  The source does not
make the false claim for an arbitrarily small `kappa`.

## 2. Net and operator constants

For fixed unit `u,v`, the coefficient squares in `u^TBv` sum to one, so
Hoeffding gives `2exp(-z^2/2)`.  A `1/4`-net has at most `9^r` points, two
nets cost `exp(2r log 9)`, and the bilinear comparison loses a factor at
most two.  Consequently `||B||>Lsqrt(r)` forces a net value exceeding
`Lsqrt(r)/2`, giving exactly

```math
2\exp\{-(L^2/8-2\log9)r\}.
```

Block triangle inequality costs one additional `sqrt(r)` from the
conference diagonal blocks.  Thus the strict endpoint

```math
L<{1\over\sqrt2\beta}-1
```

and the formula for `I_net(beta)` are correctly normalized.  Under
`C<I_net(beta)`, the interval in (LC.17) is nonempty; its lower endpoint is
automatically greater than two, so the typical-sector centering hypothesis
in (LC.3) is available.  The explicit beta choice in (LC.20) leaves unit
exponential margin as claimed.

## 3. Arbitrary events and `L^1`

Dividing the norm-tail probability by `U(F)>=exp(-Cr)` produces
`2exp(-ar)`, while division of the regular pressure tail produces
`2exp(-c r^2+Cr)`.  These bounds are uniform even when `F` is selected from
the complete pressure landscape.

On the regular event, a `2r`-spin interaction of norm at most `kappa` has
pressure at most `kappa r`.  Globally the `2r^2-r` signed edges give

```math
f\le {\beta\over\sqrt{2r}}(2r^2-r)
\le\sqrt2\beta r^{3/2}.
```

The latter, multiplied by `exp(-ar)`, is negligible.  On the regular sector
the normalized variables are uniformly bounded, so uniform convergence in
probability implies the asserted conditional `L^1` convergence.  The same
domination works for max-density laws.

Letting `L` increase to the strict endpoint after taking the `r` limsup
justifies (LC.26); no endpoint operator theorem is being assumed.

## 4. Finite negative-disorder tilt

After conditioning by an event of speed `C`, every fixed lower deviation
has any certified exponential rate `d<I_net(beta)-C`.  Splitting
`E exp(-lambda f)` at `(h_beta-eta)r` compares the exponents
`lambda(h_beta-eta)` and `d`.  Such a `d` with the desired ordering exists
exactly when

```math
C+\lambda h_\beta<I_{net}(\beta).
```

The first term is then dominant, giving the lower bound on the soft minimum;
Jensen and conditional `L^1` convergence give the matching upper bound.
All choices have fixed positive margins, so the convergence is uniform over
admissible `F`.  The source correctly stops at this finite interval and does
not infer all-fixed-tilt annealing.

## 5. KL conditional transport

For `p=U(K^c)` and `theta=q(K^c)`, binary data processing gives

```math
D(q\|U)\ge\theta\log(1/p)-\log2,
```

so `limsup theta<=C/a`.  When `C<a`, the exact KL chain rule leaves positive
regular mass and gives `D(q(.|K)||U(.|K))=O(r)`.

The transport step is correctly centered at the exact conditional/extension
mean.  Conditioning the all-`s>=0` negative MGF changes it only by
`log(1/U(K))` and an exponentially small centering error.  Optimizing at
`s=Theta(sqrt(r))` therefore costs only `O(sqrt(r))+o(r)`, proving regular
conditional mean `h_beta r-o(r)`.  Positivity on `K^c` yields
`(1-C/a)_+h_beta`; sending `a` upward proves (LC.36).

The subsequence argument in LC.5 is sound: entropy liminf below the claimed
threshold supplies a subsequence with a slightly larger uniform budget,
and LC.4 then contradicts the assumed pressure limsup.

## 6. Small-beta constant

The displayed expansions give

```math
{\gamma(\beta)\over h_\beta}
={3\beta^2\over8}+O(\beta^4),
\qquad
I_{net}(\beta)={1\over16\beta^2}+O(\beta^{-1}).
```

Their product tends to

```math
{1\over16}\,{3\over8}={3\over128},
```

so (LC.49) is correct.

## Scope

The certified entropy exponent is deliberately crude and finite.  The
theorem gives arbitrary-event typicality below that exponent, a finite tilt
window, and a one-sided KL mean floor.  It does not establish a
superexponential full lower tail, two-sided typicality under KL alone, or an
all-fixed-tilt result.

## Corrections

None required.
