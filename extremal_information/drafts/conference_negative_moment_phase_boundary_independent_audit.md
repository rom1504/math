# Independent audit: conference negative-moment phase boundary

**Verdict: PASS.**  I found no mathematical, normalization, quantifier, or
source-mapping defect requiring a repair to the frozen task-local note.

This audit applies only to
`conference_negative_moment_phase_boundary.md` at SHA-256

```text
425babaa7352a9997ee130eaa8c92f25da82e41bdf24f97068d28456964a7d87
```

It is independent of the source note and makes no change to it.  In
particular, the verdict distinguishes the proved transfer/equivalence
statements from the open conference lower-LDP lemma.

## 1. Evidentiary-status table

The opening table is accurate.

1. The convergence `R^E_(lambda,r)/r -> h_beta` below the frozen
   `lambda_*(beta)` is explicitly attributed to the earlier CT theorem and
   is not claimed as a new proof.
2. NP.1 proves only a lower bound for every fixed positive `lambda`; it
   gives a positive defect precisely on the displayed interval
   `lambda < 4 gamma(beta)/beta^2`.
3. The all-fixed-tilt assertion is explicitly open.  NP.2 is an exact
   reduction of it to superexponential lower tails, not a proof of those
   tails.

The note also correctly avoids comparing the frozen, unpublished
`lambda_*(beta)` with the explicit NP.1 endpoint.

## 2. Normalization and conference input

The conventions are internally consistent and agree with the archived
conference calculation.

* With `H_S(z)=z^T S z/2`, the block matrix in NP.1 has cross contribution
  exactly `x^T B y`, rather than twice that quantity.
* The parent inverse-temperature scale is correctly
  `t=beta/sqrt(2r)`.
* The archived limits map to
  `tau_beta=2 psi(beta)` and
  `h_beta=2 psi(beta/sqrt(2))+beta^2/4`.
* In the stated range `0<beta<sqrt(2)/6`, the archived result supplies
  `gamma(beta)=h_beta-tau_beta>0`, convergence of both orientation means,
  and `L/r -> h_beta` in probability.
* The thin-tail event has exponentially vanishing complement, which is
  more than enough for all conditioning steps used here.

No hidden factor of two enters the bridge Lipschitz estimate or the parent
subtraction.

## 3. Audit of Theorem NP.1

The constants and conditioning are correct.

Flipping one bridge entry changes the cross Hamiltonian by `2` for every
spin configuration.  Since `u -> log cosh(u)` is 1-Lipschitz, it changes
`L` by at most

```math
c_e=2t={2\beta\over\sqrt{2r}}.
```

There are `r^2` independent bridge entries, and hence

```math
\sum_e c_e^2=r^2{4\beta^2\over2r}=2\beta^2r.
```

The bounded-difference exponential estimate with constant `1/8` therefore
gives

```math
\log\mathbb E\exp\{-\lambda(L-\mathbb EL)\}
\le {\lambda^2\beta^2r\over4}.
```

Since the two orientation means are both `h_beta r+o(r)`, averaging their
negative moments preserves this exponent and yields NP.13.  For the
conditioning event,

```math
\mathbb E_{U^E}e^{-\lambda L}
\le {\mathbb E_Ue^{-\lambda L}\over U(E)},
```

so the soft minimum changes by at worst
`lambda^(-1) log U(E)=o(1)`.  Subtracting
`T_r=tau_beta r+o(r)` gives exactly

```math
\liminf {R^E_{\lambda,r}-T_r\over r}
\ge\gamma(\beta)-{\lambda\beta^2\over4}.
```

The endpoint and strict positivity assertion follow with no missing
uniformity: `lambda` is fixed throughout.

The entropy-transport reformulation is also correct.  The
Donsker--Varadhan variational inequality combined with the same MGF gives

```math
\mathbb E_U L-\mathbb E_qL
\le\beta\sqrt{rD(q\Vert U)}.
```

For the joint orientation law, entropy decomposes into the orientation
part plus conditional bridge entropies, while the two reference means
differ by `o(r)`; Jensen's inequality gives the same leading bound.  The
optimization

```math
\inf_{d\ge0}\{-\beta\sqrt d+d/\lambda\}
=-\lambda\beta^2/4
```

is exact.  Thus the stated method-optimality is appropriately limited to
the global bounded-difference/transport input; it is not advertised as an
optimal lower-tail theorem.

## 4. Audit of Theorem NP.2 and its quantifiers

The equivalence is correct, including the order of quantifiers.

Assume the lower tails are superexponential at speed `r`.  The event
`X_r <= (h+epsilon)r` has probability tending to one and supplies the lower
bound on the negative moment.  Splitting at `(h-epsilon)r` gives the upper
bound; the lower-tail term is eventually smaller than `exp(-Kr)` for every
fixed `K`.  Taking logarithms and then sending `epsilon` down to zero proves
NP.19 for each fixed `lambda`.

Conversely, Markov's elementary lower contribution from
`{X_r <= ar}` yields

```math
\limsup_r r^{-1}\log\mu_r\{X_r\le ar\}
\le-\lambda(h-a)
```

for every fixed finite `lambda`.  Given any desired finite `K`, one first
chooses the fixed number `lambda>K/(h-a)` and only then takes the limit in
`r`.  Sending this post-limit bound to `-infinity` proves NP.20.  The proof
does **not** assert control for a growing sequence `lambda=lambda_r`, and
the theorem does not need such control.

Applying the result to the joint output and to its conditioned version is
valid: the archived conference theorem supplies convergence in probability,
and conditioning on an event of probability tending to one preserves it.
The notation `e^{-omega_delta(r)}` with
`omega_delta(r)/r -> infinity` is precisely equivalent to NP.20.  A
speed-`r^2` lower bound is sufficient, while a finite speed-`r` lower rate
can create a finite transition.

## 5. Reverse-Renyi identity

The identity and its equivalence to the tail statement are exact.  From

```math
\Pi(o)={e^{L_o}\over |\Omega|a},
\qquad a=\mathbb E_Ue^L,
```

the order-`1+lambda` Renyi divergence is

```math
D_{1+\lambda}(U\Vert\Pi)
=\log a+{1\over\lambda}\log\mathbb E_Ue^{-\lambda L}.
```

Rearrangement gives NP.27.  Since `log a_r/r -> h_beta`, NP.2 says exactly
that `R_(lambda,r)/r -> h_beta` for every fixed positive `lambda` iff every
corresponding Renyi divergence is `o(r)`, iff the unconditioned fixed lower
deviations are superexponential at speed `r`.

The warning about ordinary reverse KL is valid: the `lambda downarrow 0`
limit alone does not control every fixed higher Renyi order.

## 6. LDP phase formula

Under the additionally stated good speed-`r` LDP (and hence its applicable
Laplace principle), NP.29 is correct:

```math
{R_{\lambda,r}\over r}
\longrightarrow
\inf_{a\ge0}\left\{a+{I_\beta(a)\over\lambda}\right\}.
```

The typical value `h_beta` minimizes this expression exactly when
`I_beta(a)>=lambda(h_beta-a)` for all `a<h_beta`; values above `h_beta`
cannot beat `h_beta`.  Therefore the critical typical-branch tilt is the
slope infimum in NP.31.  NP.32 is separately and correctly stated as the
criterion for crossing the smaller same-temperature target `tau_beta`.
Only `a<=tau_beta` can realize that crossing.

One scope condition is important but already present: a one-sided upper
bound is not by itself enough for NP.29; the claimed formula uses the
additional full lower-tail LDP/Laplace principle.  The note makes that
assumption explicitly.  An isolated bridge of probability
`exp(-Theta(r^2))` has infinite speed-`r` rate and therefore cannot locate
the finite slope boundary, as stated.

## 7. Primary-source audit of arXiv:2603.06368

The cited paper is correctly scoped.  Chen, Guionnet, Ko,
Lacroix-A-Chez-Toine, and Mourrat prove a speed-`N` **upper** large-deviation
principle for the ground-state maximum of Gaussian mixed `p`-spin Ising
models.  Their route uses positive fractional moments `E[Z_N^s]` for
`0<s<1` and then a positive Laplace transform of the maximum.

The paper itself distinguishes the missing opposite tail: at zero external
field, lower deviations are expected at speed `N^2`, with the cited proof
being for spherical rather than Ising models.  Consequently none of the
following gaps is cosmetic:

1. full Gaussian mixed-spin disorder versus a random Rademacher bipartite
   bridge with fixed conference diagonal blocks;
2. zero-temperature maximum versus finite-temperature normalized cosh
   pressure;
3. upper disorder deviations versus the required lower pressure deviation;
4. positive fractional powers versus negative moments/reverse Renyi.

The paper is therefore legitimate motivation for the conjectural scale of
NP.33, but supplies no theorem applicable to NP.33.  The note does not
overclaim the import.

## 8. Archive and novelty judgment

The archive comparison is accurate and suitably conservative.

* The raw MGF estimate predates this note; the new NP.1 increment is its
  transfer through the later uniform thin-tail conditioning with the full
  explicit interval.
* The reverse-Renyi algebraic identity is archived; the new useful
  synthesis is the exact all-fixed-order equivalence with superexponential
  lower pressure deviations.
* The LDP variational formula is a standard Laplace-principle consequence,
  used here as a precise diagnostic rather than claimed as a new theorem.
* The external paper motivates a direction and scale only.

I found no archived theorem that already supplies NP.33 and no primary
result in the cited source that can be imported to do so.  The proposed
conference bridge lower-LDP lemma is strictly about one explicit random
model and is not a disguised statement of the original signing
optimization.

## Final judgment

**PASS, with no source repair requested.**  NP.1 is a valid explicit
fixed-tilt wall; NP.2 is a valid iff theorem for all *fixed* tilts; the
Renyi and LDP translations are exact under their stated hypotheses; and the
March 2026 paper is correctly treated as motivation rather than imported
progress.  The remaining NP.33 lower-deviation estimate is open and is
clearly marked as such.
