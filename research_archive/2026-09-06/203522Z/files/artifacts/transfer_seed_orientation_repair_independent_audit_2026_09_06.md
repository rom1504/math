# Independent audit of all-order orientation-balanced near-minimizers

Date: 2026-09-06. Audit of
`transfer_director_balanced_orientation_repair_2026_09_06.md`.
Verdict: **PASS**. No extension of its scope is claimed.

For a fixed spin, an independently signed e-edge remainder has an
energy equal in law to a sum of e independent signs. Chernoff's bound
and the `2^N`-spin union bound give
`Pr(Q(V)>t)<=2^(N+1) exp(-t^2/(2e))`.
At `t=sqrt(2e(N+2)log2)` the right side is one half, so a deterministic
remainder with the asserted cap exists. The assertion is existential;
no efficient search is presumed. For e zero the zero remainder works.
For either orientation, maxima change by at most the actual cap of V.

Orient the seed so `P=Q>=R` and put `Delta=P-R`.
Hollow energies have mean zero under uniform spins, so `R>=0` and
`Delta<=Q`. The least k with `binom(k,2)>=Delta` obeys
`k<=sqrt(2Delta)+2`; when Delta is zero choose k zero.
For positive Delta the overshoot lies in `[0,k)`.

The negative k-clique energy is `(k-(sum x)^2)/2`. Its positive
maximum is `floor(k/2)` for BOTH parities of k, and its negative
maximum is `binom(k,2)`. Thus for the block-diagonal seed plus clique,
`P(W)=P+floor(k/2)` and `R(W)=R+binom(k,2)`.
Both lie in `[Q,Q+k]`. Filling the nk bridge edges therefore gives
`Q(seed)<=Q(B)<=Q(seed)+k+L` and oriented gap at most `k+2L`.
The lower cap bound is exact principal monotonicity, not a perturbation
estimate. If `Q(seed)<=C n^(3/2)`, then k is `O_C(n^(3/4))` and
`L=O_C(sqrt(n*k*(n+k)))=O_C(n^(11/8))`.

For the all-order version take
`r=ceil(K N^(3/4))`, `m=N-r`, with `K>sqrt(2C)+3` fixed.
At sufficiently large N, m is positive and the k obtained from an
EXACT m-minimizer is at most r. The unused `r-k` vertices can remain
isolated in the intermediate weighted matrix. The number of subsequently
filled edges is exactly
`m*r+[r(r-1)-k(k-1)]/2<=N*r`.
The same cap and orientation argument applies to that whole remainder.
Finally `M_m<=M_N` by principal monotonicity, and hence

```math
M_N\le Q(B_N)\le M_N+O_C(N^{11/8}),\qquad
|P(B_N)-R(B_N)|=O_C(N^{11/8}).
```

This final comparison assumes no continuity of normalized M and no
limit for it. Dividing by `N^(3/2)` gives exactly `O_C(N^(-1/8))`
for both error quantities. For a liminf-realizing seed subsequence,
the nearby order satisfies `(n+k)/n->1`, so the lower principal cap
and the upper construction cap preserve the same normalized limit.

`computations/transfer_seed_orientation_repair_audit_2026_09_06.py`
independently checks the finite clique formulas and actual cap inequalities
by integer spin enumeration, allowing unused new vertices as in the
all-order construction. This finite replay supplements the proof;
it is not evidence used to infer the probabilistic or asymptotic claim.

The result is existence of a selectable, normalized-orientation-balanced
near-minimizer at every large order. Its raw errors grow as `N^(11/8)`;
it does not provide `o(sqrt N)` ground slack, an isotropic ground law,
or exact minimizers. The new clique may have norm on the larger
`N^(3/4)` scale. No bounded normalized operator norm follows.
Since the construction starts from a minimizer at a nearby growing
order, it is not a fixed-seed all-order transfer and does not prove
convergence or nonconvergence.
