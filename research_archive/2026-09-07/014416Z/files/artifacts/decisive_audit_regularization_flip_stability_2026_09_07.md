# Independent audit: optimized Gaussian regularization and block comparison

Date: 2026-09-07. Verdict: PASS for the corrected version of
`decisive_bridge_gaussian_regularization_flip_stability_2026_09_07.md`.
I checked the derivative and block signs independently. The result remains
an exact reduction, not a proof of convergence.

## 1. Flip identity and endpoint

For edge `v=beta sqrt(u) A_e/sqrt(n)`, write the edge cavity partition
as a constant times `cosh(h+v)`. The cavity h can be random. Taylor
expansion of log cosh from h+v to h-v yields

`log cosh(h+v)-log cosh(h-v)=2v m-2v^2 m'+error`,

with `|error|<=(8/3)|v|^3`, because the third derivative has absolute
value at most 2. Dividing by `4u beta n` and summing fewer than n^2/2
edges gives exactly the stated `beta^2 sqrt(u)/(3 sqrt(n))` remainder.
The deterministic derivative and the negative Gaussian-variance
derivative have respectively the coefficients
`A_e/(2 sqrt(u) n^(3/2))` and `-beta/(2n^2)`.

For the ABSOLUTE partition, the expected pressure is even and convex
in a at every fixed s. Thus the nonnegative mean derivative leaves
the branch derivative bounded below by `-beta/4`. Together with
nonnegative optimizer flip costs this gives the uniform bound on S.
The minimum of finitely many smooth branches has the stated almost
everywhere envelope derivative. Integration down to zero is valid;
the bound controls the endpoint and rules out a concentrating atom.
The integrated error is `2 beta^2/(9 sqrt(n))`.

One correction found in the original draft is now incorporated:
general-weight flip stationarity requires `w_e>0`. At `w_e=0`, an
unused sign A_e can be chosen arbitrarily, flipping it has zero cost,
and division by its edge field is invalid. The main u>0 theorem is
unaffected. Zero-weight block endpoints require their own symmetry
or a justified positive-weight limit, as the corrected note states.

## 2. Block comparison and exact covariance normalization

Global spin reversal on one block, with the Gaussian cross entries
reversed too, proves evenness in the common deterministic cross-edge
coefficient. Convexity therefore permits deleting deterministic cross
means for a LOWER bound on Gaussian-averaged pressure.

For one-sided pressure, compare the remaining full Gaussian noise
with independent block noises. The full-minus-block covariance is

`(nt^2/2)[q^2-sum_i r_i q_i^2]+t^2(k-1)/2`,

for k blocks. The last scalar is exactly the self-covariance
difference and cancels from the interpolation derivative. Jensen
gives the required sign on the bracket. The child deterministic
amplitude is `sqrt(r_i)`, not 1, since its exponent must still be
`beta A_e/sqrt(n)`. These checks give the displayed LOWER comparison.

Absolute pressures differ uniformly by O(n^(-1/2)) from the maximum
of the two one-sided expected pressures, by Gaussian concentration.
The same polarity has to be used on all blocks. Moving the maximum
inside the block sum is invalid.

## 3. The natural random-cross variance repair does not chain

One might choose cross signs independently and replace them, by
Lindeberg comparison, by added cross Gaussian variance 1. This is
legitimate at fixed beta, with vanishing normalized error. It seems
to supply the missing child variance, but its comparison direction
is unhelpful.

The hybrid noise has variance t^2 internally and t^2+1 across blocks.
Compare it with independent child noises of variance
`t^2+1-r_i`. Their covariance difference, modulo the cancelling
scalar convention, is

`n(t^2+1)/2 [q^2-sum_i r_i q_i^2] <= 0`.

Guerra interpolation therefore makes the RANDOMIZED parent pressure
at least the weighted child pressure. Minimization over cross signs
only makes the optimized parent at most the randomized one. These
oppositely placed inequalities do not imply an original comparison.
Thus this concrete missing-variance repair does not settle the
uniqueness of the integrated stability density.

## 4. Independent check of the one-sided-collapse addendum

This new addendum is especially useful for scope. For a deterministic
signing A, let D_z be a uniform spin gauge. Gaussian gauge invariance
gives `P^+(D_z A D_z;a,s)=P^+(A;a,s)` after expectation. Since the
mean over z of the hollow gauged matrix is zero, convexity gives

`P^+(A;a,s)>=P^+(0;0,s)`.

Conversely the all-negative signing has
`H_A(x)=[n-(sum_i x_i)^2]/2<=n/2`. Hence

`min_A P^+(A;a,s)<=P^+(0;0,s)+a/(2 sqrt(n))`.

This proves that the optimized one-sided pressure has exactly the
Gaussian limit. The two-orientation constraint is therefore not a
small technical correction: dropping it trivializes the designer
problem. An ordinary one-sided SK limit, by itself, supplies no
absolute-signing limit theorem.

The uniformly bounded nonnegative S densities supply subsequential
compactness, not uniqueness of their total mass. Neither edge
stationarity nor ordinary Gaussian overlap identities currently
give the missing comparison of optimizer densities under block
variance redistribution.

## 5. A finite-temperature refinement of the EXISTING spectral refill

The cap-level same-order spectral refill was already proved in
`fresh_range_and_spectral_regularization_2026_09_05.md`; it is not new.
The following pressure estimate has a smaller deletion loss and follows
directly without a replacement theorem.

Let `P_n(A;a,s)` be the absolute normalized pressure, fix a vertex core
I, and let A0 retain A only on I x I, with all other entries zero. For
independent signs R on the missing edge set E0, put B=A0+R. Then

```math
P_n(A_0;a,s)\le P_n(A;a,s),\qquad
0\le \mathbb E_R P_n(B;a,s)-P_n(A_0;a,s)
\le\frac{|E_0|}{\beta n}\log\cosh\frac{\beta a}{\sqrt n}
\le\frac{\beta a^2|E_0|}{2n^2}.                    (5)
```

For the first inequality average all spin gauges on I-complement.
Each gauged signing has the same Gaussian-averaged pressure; its mean
is A0, so convexity proves the assertion. The same argument gives
`P_n(B)>=P_n(A0)` for EVERY deterministic refill. For the upper
inequality condition on all Gaussian entries and average the partition
sum over R. Each missing edge multiplies that average exactly by
`cosh(beta a/sqrt(n))`, since its signed spin coefficient has square
one. Jensen for the outer logarithm proves (5). This preserves the
absolute objective and its shared polarity throughout.

For completeness the simultaneous operator control can be imposed
without sacrificing the mean bound. A symmetric hollow independent-sign
matrix supported on ANY fixed edge mask satisfies

`Pr(||R||op>8 sqrt(n)) <= q_n:=2 exp[-(4-log 9)n]`.

For a unit vector v, Hoeffding gives
`Pr(|v^T R v|>t)<=2 exp(-t^2/4)`, because
`sum_(i<j)(2v_i v_j)^2<=2`. A 1/4-net of size at most 9^n,
with the usual factor 2 for a symmetric quadratic form, proves the
claim. The nonnegative variable in (5), conditioned on this good
operator event, has mean at most its unconditioned mean divided by
`1-q_n`. Thus at least one refill has BOTH the operator bound and
pressure increase at most the right side of (5) divided by `1-q_n`.

If `Q(A)<=C n^(3/2)`, the independently reconstructed simultaneous
diagonal majorant supplies I with fewer than delta n deleted vertices
and `||A_I||op<=4 K_G C sqrt(n)/delta`. There are at most delta n^2
missing edges. The resulting actual SAME-ORDER signing B satisfies

```math
\|B\|_{op}\le(4K_G C/\delta+8)\sqrt n,\qquad
P_n(B;a,s)\le P_n(A;a,s)
 +\frac{\beta a^2\delta}{2(1-q_n)}.                 (6)
```

In particular, at fixed beta,a>0,s, the operator-constrained optimized
pressures approximate the unrestricted optimized pressure uniformly
over orders at rate O(1/L). To see that the requisite C is uniform,
Jensen gives `P_n(A;a,s)>=a Q(A)/n^(3/2)`. Averaging a fully random
signing and using the Gaussian exponential moment gives the uniform
upper bound
`min_A P_n <=2 log(2)/beta+beta(a^2+s^2)/4`.
Thus an optimizer has bounded C depending only on the fixed parameters.
Set `delta=4 K_G C/(L-8)` in (6).

This strengthens the old cap refill's O(L^(-1/2)) loss at FINITE
temperature, but it does not itself prove convergence for a fixed
operator class. Nor does it give an operator bound uniform as the
ground-state approximation error goes to zero. It is a quantitative
comparison available to a future bounded-operator pressure argument,
not a recovered original-limit proof.

Finite diagnostic: `computations/decisive_audit_pressure_refill_2026_09_07.py`
checks 15 zero-noise parameter/core cases and all 276 associated sign
refills, with exact spin enumeration; every inequality passes. These
floating log-sum-exp checks are not a substitute for the proof above.
