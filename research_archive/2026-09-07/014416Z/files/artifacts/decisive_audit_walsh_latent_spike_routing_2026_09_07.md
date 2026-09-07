# Independent audit: finite latent Walsh spike routing

Date: 2026-09-07. PASS for
`decisive_bridge_walsh_latent_spike_routing_2026_09_07.md`, including its
explicit residual-block-formation and missing-edge repair clarifications.

## 1. Common-mask source information

For `Y=(M sigma_1,...,M sigma_k)`, the source entropy is
`h(p)+pk log2`. Independent coordinate channels conditional on `Y` give
`I(Y;L_vector)<=sum_i I(Y_i;L_i)`. Refinement by the entire label vector
can only reduce each conditional variance. Coordinatewise sign symmetry
of the common-mask source and scalar channel makes the label orientations
independent fair signs, conditionally on the absolute-label tuple; fixed
labels may receive redundant independent signs. Thus
`E[Y_i|T,epsilon]=epsilon_i a_i(T)` and the average residual covariance
is diagonal. The shared mask does not invalidate this independent sign
action.

Exact conditional physical types on label cells count actual common
selectors and actual `k`-spin tuples. Dividing by the SINGLE selector
normalizer `binom(m,ell)` subtracts `m h(p)`, leaving precisely
`m[pk log2-I(Y;L_vector)]`. Subtracting `k h(p)` would be wrong, but is
not done in the source proof. Rational/dyadic approximations are fixed
before order limits. A subsequent small change in density requires whole
physical-row repairs, with their entropy preimage cost tending to zero;
the global quadratic L2 comparison controls the associated energy change.

## 2. Literal Walsh algebra and routing

The mean function of row `i` uses exactly one orientation character
`epsilon_i`, no other orientation characters, and any shared amplitude
characters. Its support lies in the coset `epsilon_i span(T)`.
The `k` such cosets are disjoint. This is exact Walsh multiplication,
not a conclusion from pairwise orthogonality of an arbitrary Hadamard.

For each fixed amplitude frequency, a separate perfect matching of the
group graph prescribes that frequency at column `i` of its partner block.
At each group this fixes only `kS` columns. The probability cost per group
is at most `kS log m`; after `n` groups it remains `o(n^2)`.
Exact conditional types fix every coefficient in the finite label span.
Each prescribed block is diagonal, with identical diagonal spike values
at its two matching endpoints and zero offdiagonal entries. Since
`A_ii=1`, the seed reflection preserves this diagonal matrix. The folded
kernel is at least `1/2` on each matched spike edge, independently of its
possibly growing energy. Thus these are matched spikes, not spikes
discarded under an unjustified two-sided truncation.

## 3. Residual profile and exact type repair

Before repair, the centered physical vectors are independent over the
`m` rows and uniformly bounded. Their covariance is a function of the
fixed finite label bits. Distinct frequency covariance therefore vanishes
unless the frequency difference lies in the label span. Uniform column
grouping makes such exceptional one-/two-block choices an `O(1/m)`
fraction. For generic blocks, the fixed-dimensional bounded-array
characteristic-function expansion gives a joint Gaussian limit with
covariance `v I`, where `v=D_joint/p`. The two-block version gives the
empirical weak-law statement; it does not presume independent spectra.

Whole-row conditional-type repair changes only `O(sqrt(m log m))` rows
per fixed cell, has `exp(o(m))` preimages, and preserves the common-mask
source alphabet. Parseval bounds its mean squared spectral perturbation
by `o(1)`. After repair the finite label-span spikes are exact, so a second
Parseval identity fixes residual second moment to `k^2v+o(1)`. This gives
the required W2 profile convergence, including the degenerate case `v=0`.

The source now explicitly separates random formation of residual ordered
column blocks from a subsequent uniform permutation of their destinations.
The good empirical-profile event concerns the formed multiset, so conditioning
on it preserves that latter uniform row-permutation mechanism. This is
necessary for applying the positive transport theorem.

## 4. Pressure comparison and normalization

The global gradient bound for the folded kernel gives
`|log weight(X)-log weight(Z)|<=2t ||X-Z|| (||X||+||Z||)`.
Optimal empirical matching of row profiles therefore yields a genuine
two-sided W2 pressure comparison under a deterministic energy bound.
Bounded finite approximations can be used before the order limit and then
removed, without deleting uncontrolled energetic residuals.

For the finitely many prescribed missing matching edges, the source now
uses the needed row-law argument: insert fixed dummy colors, sample a
uniform complete-row arrangement, remove the missing slots, and minimally
repair at most `S` remaining colors per row to the desired counts. The
repair is permutation-equivariant, hence the target row law is uniform.
Bounded log kernels cost only `O(Sn)`. Simply deleting factors without
this row-law comparison would not by itself justify the lower bound.

The symmetric folded Gaussian transport equals the unfolded Gaussian
transport: introduce the sign as a variational variable and average the
sign/endpoint-reflection group. This fixes both conditional marginals in
each sign sector. The isotropic quadratic Gaussian calculation then gives
`k^2 g_t(v)`, independently of `A`.

Finally `m=kn`. Summing the source rate over the `n` groups contributes
`k[pk log2-I(Y;L_vector)] n^2`, while the residual pressure contributes
`k^2 g_t(D_joint/p) n^2`. Applying the information and variance bounds,
and monotonicity of `g_t`, gives exactly
`k^2[p log2-I_scalar+g_t(D_scalar)]`. Supremizing finite scalar channels
is the stated `k^2[p log2+E_t(nu_p)]` lower bound. All channel and support
parameters remain fixed during the order limit.

## 5. Audited exact self-block extension

The amended source proof restores all within-group factors exactly. The
full finite label span has `2^(r+k)` columns, while the union of the row
mean supports has at most `k 2^r`. Each remaining label-span column is
zero in every row under the exact conditional types. There are at least
`2^r(2^k-k)>=k` such columns for k>=2. Routing k distinct common zero
columns to the self block costs only another `O(kn log m)` in the log
probability, and every within-group defect is then zero. This verifies
the FULL positive partition extension, with the source's explicit type
approximation and parity qualifications retained.

The audit applies to the specified Walsh/full-column-permutation ensemble.
It does not extend to arbitrary recursively randomized bases or assert a
lower bound on the original minimizing cap. Those exclusions remain real.
