# Wave 23 Route 3: nonuniform flip weights and the exact dual-sign wall

This note asks whether optimized nonuniform weights or LP duality can turn the
edge-cube cover (10.733) into a Gibbs-weighted congestion/common-coordinate
bound, improving the uniform wrong-way transform (10.735).  The answer for
this mechanism is negative and sharp.  Uniform perturbation weighting is the
minimax-optimal nonnegative transform for a coordinate star; the full
fractional-cover dual still lower-bounds affinity; and an exact two-center
abstract star model satisfies every nonuniform perturbation inequality and
the coordinate pair law while giving only the parity-scale cavity reward.

This is a scoped certificate/row-marginal wall.  It neither constructs an
asymptotic complete-signing minimizer nor falsifies (10.617).  All finite LPs
and enumerations below are verified by `nonuniform_cover_dual_r23.py`.

## 1. The arbitrary-weight transform

Keep the notation of (10.733).  For a block `H\subseteq E(K_r)`, let
`h=|H|`, and for every `S\subseteq H` set

```math
D_S
=D_\beta(B^S;q)
=\sum_\omega
\exp\left\{-\beta\left[
\Delta_\omega+4\sum_{e\in S}s_e(\omega)
\right]\right\}.
```

Exact global minimality gives `D_S\ge1` separately for every `S`.  For an
arbitrary probability law `p` on the perturbations define

```math
F_p(s_H)
=\mathbb E_{S\sim p}
\exp\left(-4\beta\sum_{e\in S}s_e\right).
```

Then the strongest direct nonnegative aggregation is the identity

```math
\boxed{
\mathbb E_{S\sim p}D_S
=D_\beta(B)\,
\mathbb E_{\nu_{\beta,B}}F_p(s_H)
\ge1.
}
\tag{R23.1}
```

More generally, for arbitrary coefficients `\lambda_S\ge0`, the right side
is at least `\sum_S\lambda_S`.  Nonnegativity is essential: it is the only
reason the inequalities `D_S\ge1` may be summed.

For a coordinate star `H=\delta(i)`, put `m=r-1` and

```math
g(s)=\exp\left(-2\beta\sum_{e\in H}s_e\right).
```

Then `\mathbb E_\nu g=\operatorname{BC}_i`.  Any pointwise comparison

```math
F_p(s)\le C_p g(s)
```

used in (R23.1) gives only

```math
\operatorname{BC}_i\ge\frac1{C_pD_\beta(B)}.
\tag{R23.2}
```

Thus arbitrary weighting has the same sign problem as (10.735).

## 2. Uniform weighting is exactly minimax

The sign obstruction cannot be improved even within its wrong direction.
Encode `S` by `a_S\in\{-1,1\}^m`, with coordinate `-1` on `S` and `1`
off `S`.  Directly,

```math
\frac{F_p(s)}{g(s)}
=\mathbb E_{S\sim p}e^{2\beta\langle a_S,s\rangle}.
```

Averaging this ratio uniformly over `s\in\{-1,1\}^m` gives, for every `p`,

```math
2^{-m}\sum_s\frac{F_p(s)}{g(s)}
=\cosh(2\beta)^m.
```

Consequently

```math
\boxed{
\min_{p\in\mathcal P(2^H)}
\max_{s\in\{-1,1\}^m}\frac{F_p(s)}{g(s)}
=\cosh(2\beta)^m.
}
\tag{R23.3}
```

Uniform `p` attains equality at every `s`, reproducing (10.735).  For
`\beta>0`, the tensor exponential kernel is invertible, so uniform `p` is the
unique minimizer.  Hence no nonuniform perturbation law gives a stronger
universal pointwise extraction of affinity than the existing uniform
identity.

## 3. The exact fractional-cover dual also has the wrong direction

For the star at `i`, let

```math
A_i(S,\omega)
=\mathbf1\left\{
\Delta_\omega+4\sum_{e\in S}s_e(\omega)\le0
\right\}
```

be the exact witness incidence matrix, and set

```math
a_i(\omega)=e^{-\beta(\Delta_\omega+2h_i(\omega))},
\qquad
G_i=\sum_\omega a_i(\omega)
=D_\beta(B)\operatorname{BC}_i.
```

All row perturbations are covered, so the all-one state vector is feasible
in the weighted fractional-cover primal

```math
\tau_i
=\min\left\{
\sum_\omega a_i(\omega)x_\omega:
x_\omega\ge0,
\ \sum_\omega A_i(S,\omega)x_\omega\ge1\quad( S\subseteq H)
\right\}.
```

Its exact dual is

```math
\boxed{
\tau_i
=\max\left\{
\sum_S\lambda_S:
\lambda_S\ge0,
\ \sum_S A_i(S,\omega)\lambda_S\le a_i(\omega)quad(\omega)
\right\}
\le G_i.
}
\tag{R23.4}
```

Thus optimized witness congestion supplies a lower bound on `G_i`, hence a
lower bound on `\operatorname{BC}_i` and an upper bound on the reward.  This
LP is the strongest linear consequence of the exact cover with these state
costs.  Optimizing jointly over coordinates can only lower-bound a
nonnegative weighted average of the `G_i`; it can force a coordinate with
*large* affinity, not a common coordinate with small affinity.

This is not an artifact of choosing the wrong primal.  Cover feasibility is
monotone upward in the available state weights.  The desired conclusion
`G_i` small is in the opposite, downward direction.

## 4. Sharp upper-moment LP and a two-center wall

One can add more than incidence: retain every partition inequality from
(R23.1), complement symmetry, and the exact coordinate-pair Gibbs law.  Let
`\nu(y)` be the star-pattern marginal, `h(y)=\sum_e y_e`, and impose

```math
\nu(-y)=e^{-4\beta h(y)}\nu(y)
\qquad(h(y)>0),
```

together with

```math
\mathbb E_\nu e^{-4\beta\sum_{e\in S}y_e}
\ge\frac1{D_\beta(B)}
\qquad(S\subseteq H).
\tag{R23.5}
```

The linear program maximizing
`\mathbb E_\nu e^{-2\beta h(y)}` under all these constraints has the exact
value

```math
\boxed{
\sup\operatorname{BC}
=\operatorname{sech}(2\beta\varepsilon),
\qquad
\varepsilon=m\bmod2.
}
\tag{R23.6}
```

The upper bound follows immediately from the pair representation (10.710):
every unordered pair contributes `\operatorname{sech}(2\beta u)`, where
`u\ge\varepsilon` by parity.  More importantly, none of the simultaneous-
flip moment constraints lowers this ceiling.

To see sharpness, choose a pattern `y` with `h(y)=\varepsilon`.  Give `y`
deficit zero and `-y` deficit `4\varepsilon`.  Their exact restricted radii
from (10.734) are both

```math
\rho=\frac{m-\varepsilon}{2}=\left\lfloor\frac m2\right\rfloor.
```

The two antipodal radius-`\rho` balls cover the entire `m`-cube.  The Gibbs
weights are `1` and `e^{-4\beta\varepsilon}`, so their affinity is precisely
`\operatorname{sech}(2\beta\varepsilon)`.  For any perturbation `S`, write
`a=\sum_{e\in S}y_e`.  Its normalized partition moment is

```math
\boxed{
\frac{D_S}{D}
=\frac{\cosh(2\beta(2a-\varepsilon))}
       {\cosh(2\beta\varepsilon)}
\ge1.
}
\tag{R23.7}
```

The last inequality holds because `2a-\varepsilon` has the same parity as
`\varepsilon`, hence absolute value at least `\varepsilon`.  This two-center
model therefore satisfies, simultaneously:

- the complete restricted star cover;
- `\Delta_{-y}=\Delta_y+4h(y)`;
- the radius/exponent identity (10.734);
- complement symmetry `D_S=D_{H\setminus S}`;
- every nonuniform inequality (R23.1), in fact with `D_S\ge D`;
- and the maximal affinity in (R23.6).

For odd `r` (`m` even), this permits `\operatorname{BC}=1` and zero reward.
For even `r` (`m` odd), it gives only the universal parity reward

```math
\kappa\ge\frac1\beta\log\cosh(2\beta),
```

which is constant.  The model is an abstract star marginal, not the full
augmented-cut state space of a complete signing.  It is therefore a precise
wall for the weighted-cover/moment axioms, not a counterexample to
(10.617).

## 5. Where an apparent upper bound becomes invalid or circular

There are only two ways to try to reverse (R23.1), and both expose the
missing input.

1. **Signed perturbation coefficients.**  Negative coefficients cannot be
   multiplied into `D_S\ge1` while retaining the inequality direction.
   The exponential functions indexed by all star subsets form a basis, and
   the unique exact expansion of `g` has the positive uniform coefficients
   in (10.735).  An upper-affinity dual would require coefficients in the
   unavailable sign cone.
2. **Insert the actual competitor sums.**  If one keeps
   `\sum_S\lambda_SD_S` rather than replacing it by its lower bound, a
   pointwise lower domination of `g` can indeed yield an upper bound.  But
   the right side now contains the unknown competitor partition functions.
   For uniform weights it is exactly
   `D\cosh(2\beta)^m\operatorname{BC}_i`; the proposed bound is an identity.
   Controlling the weighted competitor sum from above is therefore the
   desired Gibbs-overlap statement in another form, not a consequence of
   minimality.  The generic estimate
   `D_S\le2^r e^{4\beta\min(|S|,m-|S|)}` is exponentially too large and gives
   no nontrivial affinity bound.

Thus a proposed “optimized dual” must display an independent upper bound on
competitor partition sums or a Gibbs-weighted Hall/congestion theorem.  If
its weights are chosen using the unknown witness Gibbs masses or the target
affinity itself, that is precisely the circular step.

## 6. Exhaustive `A_8/A_9` audit

The checker enumerates every oriented state, every coordinate, and every one
of the `2^{r-1}` star perturbations.  At `\beta=0.5,1` it solves (R23.3), the
weighted cover dual (R23.4), and the all-moment upper LP (R23.5).

- The minimax LP always returns `\cosh(2\beta)^{r-1}` and the uniform law to
  numerical precision, for both minimizers and every audited temperature.
- For `A_8`, the optimized cover lower bound `\tau_i/G_i` is about `0.05996`
  at `\beta=0.5` and `0.00925` at `\beta=1`; all coordinates are equivalent.
- For `A_9`, `\tau_i/G_i` ranges from `0.04251` to `0.12247` at
  `\beta=0.5`, and from `0.001247` to `0.01827` at `\beta=1`.  Even the fully
  optimized exact-cover LP captures only a small fraction of the actual
  geometric mass, always from below.
- The all-moment upper LP returns exactly `\operatorname{sech}(2\beta)` for
  `A_8` (`m=7`) and `1` for `A_9` (`m=8`), matching the two-center wall.  It
  does not recover the much smaller actual affinity ranges
  `0.323374`/`0.104231` on `A_8` and `0.263762--0.430737`/
  `0.172025--0.334714` on `A_9`.

These are exhaustive star audits; arbitrary full edge-set coverage remains
guaranteed algebraically by (10.733), rather than by enumerating `2^N`
perturbations.

## 7. Scoped conclusion

- **Verified strongest transform:** (R23.1) is the arbitrary nonnegative
  perturbation-weight theorem, and uniform weighting is its sharp minimax
  specialization for a star, (R23.3).
- **Verified exact LP direction:** fractional weighted congestion gives
  `\tau_i\le D\operatorname{BC}_i`, never the upper affinity needed for a
  large cavity reward.
- **Sharp abstract wall:** even all perturbation moments plus exact pair
  balance and the restricted cover allow affinity `1` or
  `\operatorname{sech}(2\beta)` according only to parity.
- **No route to (10.617):** nonuniform weighting does not improve
  `\max_i\kappa_i` beyond a constant parity bound.  A continuation needs new
  complete-signing structure that upper-bounds a competitor partition sum or
  forces Gibbs-weighted witness congestion on one common coordinate; exact
  cover LP duality alone cannot provide it.

