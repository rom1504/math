# Independent audit of the projective-shell local-field roof

Date: 2026-08-17.

**Verdict: PASS after two scope/terminology repairs.**  The final audited
draft has SHA-256
`ef83fc221171fd31d648f72b35a80346d21c1f44ac91d00f70716f9ebb57f7f2`.
Lemma PR.1, the exact
two-sided error bound in Theorem PR.2, the sorting evaluator, the bit count,
and the comparison with AO.2 are correct.  The supplied finite verifier also
passes.  The two repairs do not change the theorem:

The repaired points are:

1. an ordinary small radius-`R` covering number does not by itself bound the
   size of the inclusion-maximal `R`-separated net used by the theorem;
   Section 5 should instead ask for a small maximal internal net, a bound on
   the `R`-packing number, or (sufficiently) a small radius-`R/2` cover;
2. when `G=o(n^(3/2))`, target-scale approximation of the declared query
   ball is already trivial from the single scalar `Q(A)`.  The atlas is
   genuinely finer only when its error is compared with `G` (for example
   `k_R^2=o(G)`), or when `G` itself is macroscopic.  Thus the valid theorem
   is a fine-scale static response result, not on its own a strict reduction
   of the near-minimizer/convergence frontier.

The final draft now distinguishes the net/cover notions and states the
scalar baseline and fine-scale qualification explicitly.  No canonical file
was changed in this audit.

## 1. Projective geometry and orientation bookkeeping

For signed augmented cuts `z=sigma c(x)` and
`z_0=sigma_0 c(u)`,

```math
d_E(z,z_0)={E-\langle z,z_0\rangle\over2},
\qquad
d_{\rm P}(z,z_0)={E-|\langle z,z_0\rangle|\over2}.
```

Choose `eta` so that `d_E(z,eta z_0)=d_P(z,z_0)`, and let

```math
d=\min\{d_H(x,u),n-d_H(x,u)\}\le\lfloor n/2\rfloor.
```

If `sigma=eta sigma_0`, this oriented distance is `d(n-d)`.  If the
orientations disagree, it is instead

```math
E-d(n-d)\ge E-\lfloor n^2/4\rfloor.
```

The strict radius condition in PR.4 therefore rules out the second branch,
and monotonicity of `d(n-d)` on `[0,n/2]` gives `d<=k_R`.  After the harmless
global-spin choice `tau`, this is exactly

```math
sigma=eta sigma_0,\qquad x=tau u^S,\qquad |S|=d.
```

Thus Lemma PR.1 is correct.  The hypothesis `n>=3` is useful: an augmented
cut then has a unique orientation (and a vertex representative unique up to
global spin).  The displayed `d_P` is a metric on the projective quotient;
on `mathcal Z_n` itself it is a pseudometric because `z` and `-z` have
distance zero.  This harmless convention should be kept in mind when using
the word “family.”

## 2. The shell restriction and local Taylor identity

For a signed ground state `(sigma_*,x_*)`, global negation of `x_*` permits
the favorable field sign, so

```math
\mathcal B_A(g)\ge Q+|g\mathbin\cdot x_*|\ge Q.
```

At an optimizer `(sigma,x)`, consequently,

```math
sigma H_A(x)
=\mathcal B_A(g)-g\mathbin\cdot x
\ge Q-\|g\|_1\ge Q-G.
```

The positive shell really has width `G`, not `2G`; no absolute-value or
global-spin loss is missing.

For a centre, write

```math
d^r_{ij}=sigma_r a_{ij}u_i^ru_j^r,
\qquad
ell_{r,i}=\sum_{j\ne i}d^r_{ij}.
```

If `S` is flipped, direct edge bookkeeping gives

```math
sigma_rH_A((u^r)^S)
=h_r-2\sum_{i\in S}ell_{r,i}
  +4\sum_{\{i,j\}\subseteq S}d^r_{ij}.
```

After changing the quadratic orientation by `eta`, the omitted term is
multiplied by `eta` and has absolute value at most

```math
4{|S|\choose2}\le2k_R(k_R-1).
```

The field term for the genuine competitor
`(eta sigma_r,tau (u^r)^S)` is exactly

```math
tau g\mathbin\cdot u^r-2tau\sum_{i\in S}g_i u_i^r.
```

Hence every atlas tuple is within the displayed error of an actual response
competitor.  Conversely, maximality of an inclusion-maximal separated net
covers the optimizer's shell atom, and Lemma PR.1 supplies the required
tuple.  Both directions of PR.11 pass; there is no separately paid
orientation or field channel.

## 3. Net size and the precise packing/covering statement

In a finite metric space, every inclusion-maximal family with mutual
distance `>R` is an `R`-cover.  Thus the particular `mathcal C` used in the
proof simultaneously is a separated packing and a cover, exactly as
claimed.

There is, however, a standard same-radius distinction relevant to the scope
paragraph.  A small *minimum* `R`-cover need not imply that an arbitrary (or
even maximum) `R`-separated family is small: one radius-`R` ball may contain
many points mutually more than `R` apart.  What suffices for PR.21 is any of
the following:

* a directly constructed inclusion-maximal internal `R`-net of size
  `o(n/log n)`;
* an `o(n/log n)` upper bound on the `R`-packing number;
* an `o(n/log n)` cover by radius-`R/2` balls (by the triangle inequality).

Accordingly, “any bound on the projective covering number” in the current
scope statement is too loose if it means the ordinary same-radius covering
number.  This is a terminology repair, not a defect in PR.2.

For `R=floor(gamma E)`, integer separation gives

```math
d_P(z^r,z^s)>gamma E,
\qquad
{|\langle z^r,z^s\rangle|\over E}<1-2gamma.
```

The floor therefore causes no loss.  The radius restriction implies
`2gamma<1` at finite `n`, so AO.2 may indeed be invoked with its projective
gap parameter `Gamma=2gamma`.

## 4. Evaluation and retained information

For fixed `(r,eta,tau)`, the only dependence on `S` is

```math
\sum_{i\in S}
(-2eta ell_{r,i}-2tau g_i u_i^r),
\qquad |S|\le k_R.
```

Taking the largest at most `k_R` positive increments is exact.  Four sorts
per centre therefore evaluate the formula in `O(Ln log n)` comparisons.

All stored values are integral.  A centre and its orientation cost `n+1`
bits, `h_r` costs `O(log n)` bits, and the `n` local fields in
`[-(n-1),n-1]` cost `O(n log n)` bits.  Thus PR.17 is correct.  Also

```math
\sum_i ell_{r,i}=2h_r,
```

so the stated redundancy of the baseline is exact.  If
`L=o(n/log n)`, this is genuinely `o(n^2)` retained data and need not store
the edge matrix.  The theorem correctly does not claim an efficient encoder:
finding the shell/net can still require the original hard search.

## 5. Scaling and the AO.2 comparison

Because `k_R<=n/2`,

```math
k_R(n-k_R)\le gamma E,
\qquad n-k_R\ge n/2,
```

and hence

```math
k_R\le {2gamma E\over n}=gamma(n-1).
```

The response error is therefore at most
`2gamma^2(n-1)^2`, and is `o(n^(3/2))` under
`gamma=o(n^(-1/4))`.  These constants and the critical exponent pass.

On the separated branch AO.2 receives `Gamma=2gamma`.  Its exact leading
coefficient (before lower-order terms) is

```math
{gamma(1-1/n)\over2}\min\{alpha,h/\sqrt n\},
```

so the draft's order notation
`asymp gamma min(alpha,lambda)n^(3/2)` is correct.  With fixed positive
`alpha,lambda`, positivity uniformly over a vanishing scale requires

```math
G=o(gamma n^(3/2)),
\qquad
gamma\gg n^{-1/4},
```

the latter because the simultaneous sparse-flip fluctuation is
`O(sqrt(alpha)n^(5/4)+n)`.  Optimizing `alpha<=1` cannot improve that
exponent.  Thus the stated mismatch between the two existing tools is real.

One wording nuance is useful: if `gamma->0`, AO.2 gives a physical gap of
order `gamma n^(3/2)`, not a fixed macroscopic gap.  “Physical branch above
the critical scale” is correct only at that scale of accuracy; a fixed
target-scale physical packing follows when `gamma` is bounded away from
zero.

## 6. Nontriviality audit: the subleading-query caveat

For every query in the declared ball there is a universal one-scalar
sandwich, independent of shell geometry:

```math
Q\le\mathcal B_A(g)\le Q+\|g\|_1\le Q+G.           \tag{A.PR.1}
```

Thus retaining only `Q` already approximates the whole ball to error `G`
(and the estimator `Q+G/2` has worst-case error at most `G/2`).  It follows
that if `G=o(n^(3/2))`, an `o(n^(3/2))` response approximation is automatic
and contains no near-minimizer structural information.

PR.2 is nevertheless nontrivial at finer resolution.  A simple sufficient
regime in which it beats the scalar baseline asymptotically is

```math
k_R^2=o(G),
```

or, using the displayed coarse radius estimate,

```math
gamma=o(\sqrt G/n).
```

It is also target-scale nontrivial when `G` itself is of order
`n^(3/2)` and the atlas error is subleading.  These are the regimes in which
the local fields genuinely answer information that `Q` alone does not.

Consequently the theorem should be classified as:

* a rigorous universal local-field response approximation;
* a strict subquadratic *representation* conditional on a small maximal
  net;
* a useful fine-scale packing-versus-local-chart dichotomy;
* **not**, without an additional nontrivial accuracy comparison or update
  law, a strict reduction of the original near-minimizer frontier.

In particular, saying that PR.2 by itself “sharpens `L_projective`” should
be understood as organizing its multiscale geometry, not as completing a
new implication arrow toward convergence or physical low-cap compression.

## 7. Finite verifier

I independently ran

```text
python3 extremal_information/experiments/verify_nearmin_projective_shell_roof.py
```

and obtained

```text
PASS: projective shell roof; 88696 matrix/query/radius checks; max centres 51
```

The script checks all admissible radii for all hollow signings through order
five (with exhaustive/broad deterministic query banks as documented), then
random real-valued fields through order nine.  It verifies cover and
separation, the signed-chart representation, sorting versus subset
enumeration on a sample, and the exact error constant.  It is a strong
regression check, though naturally it does not verify the information-cost
interpretation, the AO.2 asymptotics, or the subleading-query caveat; those
were checked analytically above.

## Final recommendation

**PASS PR.1--PR.3 without formula changes.**  Before canonicalization, add
the same-radius packing/covering qualification and the trivial baseline
(A.PR.1).  With those scope repairs, the exact frontier is clear: this is a
valid one-block, formal-field, fine-resolution compression theorem with no
transition law, no all-order realization, and no theorem that exact
minimizers have the required small net or large resolvable packing.
