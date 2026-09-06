# Independent audit: the optimized bare certificate is not the original optimum

Date: 2026-09-06. Status: PASS for
`decisive_independent_bare_scalar_optimum_not_original_2026_09_06.md`.
I independently reconstructed the parameter compactness, continuity,
partition normalization, and strict finite-temperature extraction.

## 1. Exact pressure normalization from the signing construction

Use the notation of `transfer_reconstruction_standalone_2026_09_06.md`,
Sections 5--6: ambient fibre size `m`, retained size `k`, `N=mk`,
`p=k/m`, and full retained matrix `W_T`. Its exact defect satisfies

```math
D_sigma=2(m²k-sigma x^T W_T x).
```

Thus, before the Markov/union-bound extraction,

```math
exp[-tD_sigma/(2k)]
=exp[-tm²]exp[(t/k)sigma x^T W_T x].
```

The PSD/Finner bound used for (34), summed over all spins and orientations,
therefore bounds the expectation of the actual partition function, not
merely the event probability at one energy level. For the hollow matrix
`A=W_T-diag(W_T)`, its energy is `H_A=x^T A x/2`. Consequently the correct
inverse temperature is

```math
beta=2t/k.
```

Removing the diagonal multiplies the partition upper bound by at most
`exp(t |Tr W_T|/k)<=exp(tm)`, whose logarithm is `o(m²)`. The factor 2
for the two orientations also has negligible logarithm. Independent
fibre bases and the fixed-depth type theorem hence give a deterministic
signing with

```math
log(Z_++Z_-)/m² <= t+p log2+(B^r Phi_t)(nu_p)+o(1).
```

The new exact deep value `H=E` allows the finite-depth term to approach
`E_t(nu_p)` from above, with depth fixed before taking the order limit.

## 2. Entropy extraction and strict improvement

Around either signed absolute extremizer, independent spin flips of
probability `delta` give entropy `N h(delta)` and expected signed energy
`(1-2delta)² Q(A)`. The finite Gibbs variational bound is therefore

```math
log(Z_++Z_-) >= N h(delta)+beta(1-2delta)² Q(A).
```

Here `N/m²=p` and `beta N^(3/2)/m²=2t sqrt(p)`, exactly producing

```math
C_delta(p,t)=
[t+p log2+E_t(nu_p)-p h(delta)]/[2t sqrt(p)(1-2delta)²].
```

The condition `C_delta<C` is precisely
`h(delta)/delta > (8tC/sqrt(p))(1-delta)`. It holds at every fixed
finite positive `p,t,C` for sufficiently small positive `delta`.
This is an actual-spin extraction, not a Gaussian or fractional-spin
argument, and it needs no equality in the Finner pressure bound.

## 3. Why optimizing cannot escape to a boundary

The trivial channel and full-revelation channel give all three lower
bounds used in the source. For `t->infinity`,
`C>=(1-log(2)/t)/2` eventually, uniformly in `p`. For `p->1` with bounded
positive `t`, the full-revelation bound tends to `1/2`. The other separate
boundaries diverge as stated.

The remaining simultaneous boundary `p->0,t->0` is important: the exact
small-`t` expansion is `t+g_t(1)=t²+O(t⁴)`. Consequently, with
`r=t/sqrt(p)`, the lower bound is
`[(1+o(1))r+log(2)/r]/2`, at least
`sqrt((1+o(1))log(2))`. This is uniformly above `1/2`; no restriction on
the limiting behavior of `r` is needed. Hence a sublevel below `1/2`
is contained in a compact interior parameter rectangle.

For continuity, transport a source label through a coupling at `W2`
distance `d`. Information decreases by data processing. The new
mean squared reconstruction error is at most the old one plus
`2 sqrt(D_old)d+d²<=2d+d²`; `g_t` is decreasing and `t`-Lipschitz in
variance. Reverse the coupling argument to obtain the absolute
continuity bound. Its temperature derivative has magnitude at most
the residual variance, at most one here. This verifies joint continuity
on the compact rectangle and attainment of `C_*` in its interior.

## 4. Order of choices and resulting scope

First select the attained parameters `(p_*,t_*)`, then a fixed positive
`delta_*` giving a strict extraction improvement, then a finite depth
whose approximation error is smaller than part of that fixed gap.
Finally take the constructed-order limit and use the already proved
relatively dense Hadamard/all-order interpolation. This gives
`limsup M_n/n^(3/2)<C_*`, not merely a sequence of pointwise improvements
whose gains could disappear after optimization.

The theorem is conditional only on the explicitly named, independently
proved construction-pressure realization and one bare value below `1/2`.
It disproves the proposed identification of the UNCORRECTED scalar
certificate optimum with the original minimax limit. It does not establish
convergence or exclude an appropriately entropy-corrected construction
from ultimately being optimal.
