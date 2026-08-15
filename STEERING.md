# Strategic steering

Evidence cutoff: Section 10.139 (2026-08-15). Status: **consolidation
diagnostic active; automatic neighboring-route generation paused**.

Next mandatory blank-slate refresh: this consolidation diagnostic is the
current refresh. If ordinary numbered waves later resume, keep Wave 61 as the
next scheduled boundary; refresh earlier after a decisive result or route
change.

## User-stated objective and workflow directives

The objective is to determine whether `M_n/n^(3/2)` converges. The value
`1/2` is a conjectural mathematical assessment, not a separate user
objective. Reproducibility, verification, Git checkpoints, the investigation
of both convergence and nonconvergence, and the README stopping rule remain
user directives.

The user authorized a second sustained autonomous phase. Their latest input
supplied external research suggestions, not mathematical directives. Those
suggestions have been independently tested. Route selection and all
conjectures below remain agent-authored. `README.md` is unchanged.

## Agent-authored frontier

```math
0.336493364431\ldots
\le \liminf_{n\to\infty}{M_n\over n^{3/2}}
\le \limsup_{n\to\infty}{M_n\over n^{3/2}}
\le {1\over2}.                                      \tag{S1}
```

No exact value, asymptotic bound, recurrence step, convergence theorem, or
genuine nonconvergence mechanism improved at the last two substantive
checkpoints.

## What the moving-representation audit proved

For an operator kernel `K` on the edge cube whose pointwise values, Fourier
coefficients, and complete remainder coefficients have the required Loewner
positivity, the full same-switch theorem is

```math
\boxed{(\lambda-\mu(a))_+T_a\preceq J.}             \tag{S2}
```

It preserves joint matrix cancellation. If the root is constant, however, a
least generalized eigenvector scalarizes the certificate with exactly the
same bound. Correct scale requires Fourier support
`exp(Omega(n log n))`; rank alone supplies no strength.

Partial matchings have that support size and internal transition order, but
the two direct uses are closed:

```math
\lambda-{J\over T_a}\le {2\over\sqrt{\binom n2}}=O(n^{-1})
```

for the radial matching square, and

```math
\lambda\le {1\over n-1}
```

for every arbitrary-rank kernel whose Fourier support is directly on
matchings. The second result uses an operator star marginal and assumes no
radiality, equivariance, commutativity, or scalar decomposition.

Every pure add/delete partial-matching representation also satisfies

```math
\lambda\le
\left({4\over3\sqrt3}+o(1)\right)n^{-1/2},           \tag{S3}
```

strictly below the coefficient `1` required by the conference-scale
convergence criterion. Broader difference-of-matchings blocks remain open,
but their signed-hafnian root is not demonstrably simpler than the original
tail and the pure matching transition ceiling still applies.

## Exact sufficient convergence lemma

A moving-kernel family defined independently of `M_n` would prove convergence
to `1/2` if

```math
\lambda_n={1-o(1)\over\sqrt n},
\qquad
\sup_a\lambda_{\min}(T_a^{-1/2}J_nT_a^{-1/2})
=o(n^{-1/2}).                                      \tag{S4}
```

Then (S2) gives `Q(a)>=(1/2-o(1))n^(3/2)` uniformly, and conference orders
supply the matching upper subsequence. The remaining root estimate can encode
the full signed coset histogram, so (S4) is not presently a strict reduction.

## Known obstructions

- The imported packing bound is translation-blind; only the rooted cross term
  in (S2) sees a cut-code coset.
- Canonical radial moving kernels collapse to RMS on tested signings.
- Scalar partial transversals have too little support by the
  Bollobas--Lee--Letzter theorem.
- Constant-root operator rank scalarizes exactly.
- Direct matching support, radial matching squares, and pure matching
  add/delete transitions are closed by the new scalable theorems.
- Fixed small disorder tilts retain a positive linear conference defect; the
  known algebraic bridge costs `Theta(n^2)` entropy.
- Selected-prior/common-active-face, scalar atoms, separately paid channels,
  ordinary polarization, same-map Gaussian response, fixed-level SOS, Walsh
  basin classification, and eigenvalue-only polar variants remain inactive.

## Ranked formulations entering the diagnostic

No route is currently entitled to automatic continuation. The bounded
blank-slate audit must first formulate at most three mechanisms from the
original problem and only afterward compare them with the archive.

Provisional alternatives to compare fairly are:

1. a graph-orbit moving representation outside partial matchings, but only if
   it has coefficient one and an algebraic root identity satisfying (S4);
2. a diffuse entropy-tilted bridge law with `O(n)` relative entropy and a
   linear pressure gain, yielding a summable composition defect;
3. a genuinely different compactness, interpolation, exchange, or
   nonconvergence mechanism with an exact theorem and falsifier.

The diagnostic must reject vocabulary-only reformulations and explicitly
test whether each obligation is weaker than bare parent/coset optimization.

## Genuine nonconvergence standard

A route falsifier is not evidence of nonconvergence. A genuine result must
produce some fixed `epsilon>0` and two infinite subsequences whose normalized
values differ by at least `epsilon`, or equivalently prove
`liminf<limsup`. No current construction does this.

## Consolidation decision

This is the second consecutive substantive checkpoint without primary
progress. The moving-representation analogy yielded rigorous mathematics but
not a strict reduction of the original problem. The automatic loop is paused
under the README rule while one bounded blank-slate diagnostic is run. After
its synthesis, make and commit one evidence-based choice: resume with one
named theorem target, change architecture, seek external mathematical review,
or stop. Do not open another neighboring kernel variant by default.
