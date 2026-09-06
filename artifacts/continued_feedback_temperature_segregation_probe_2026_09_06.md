# Bounded temperature-segregation diagnostic

Date: 2026-09-06. Numerical necessary-condition test only; no certified
counterexample and no global temperature-alignment theorem is claimed.

The symmetric input has magnitudes a and 1 with masses q and 1-q.
Use an admissible deterministic pair: at magnitude a take B=A, while at
magnitude 1 take B=-A. The pair information cost is exactly
`H(input)=h(q)+log2`. The two symmetrized output laws are ternary:

- plus: mass q at magnitudes sqrt(2)a, and mass 1-q at zero;
- minus: mass 1-q at magnitude sqrt(2), and mass q at zero.

Thus both child latent envelopes are covered by the exact analytic
ternary concave-hull reduction. Their floating outer evaluations give
the selected Bellman lower diagnostic

`[E(plus)+E(minus)-h(q)-log2]/2`.

The source envelope has a separate finite-reproduction-grid RD dual
upper calculation, with its explicit unbiased-rounding correction.
Floating arithmetic was used throughout this diagnostic, so even these
enclosures are not asserted as interval certificates.

The grid was q in {.1,.25,.5,.75,.9}, a in {.1,.4,.8}, and
t in {.25,1,4}: 45 cases, all completed without errors. The closest
reported lower-policy/source-upper difference was -.22414558988 at
q=.1,a=.1,t=4. This family therefore gave no promising numerical
falsifier, and the search was stopped at the declared finite scope.
Failure to find a violation is not a proof of `B E<=E` even within
this policy family.

Files:

- `computations/continued_feedback_temperature_segregation_probe_2026_09_06.py`
- `computations/results/continued_feedback_temperature_segregation_probe_2026_09_06.json`

The small-temperature child region is treated analytically: if the
normalized child temperature divided by its retention is at most 1/2,
the whole ternary inner maximizer has zero posterior magnetization and
the child envelope is the Gaussian branch exactly. This avoids using
a reversed numerical outer interval. The saved grid also keeps the
ordered coexistence fields in a numerically stable finite range.
