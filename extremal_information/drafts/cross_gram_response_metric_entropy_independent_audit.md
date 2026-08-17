# Independent audit: cross-Gram response metric entropy

**Audited files.** `cross_gram_response_metric_entropy.md` and
`verify_cross_gram_response_metric_entropy.py`.

**Verdict.**  GE.1's covering exponent and constants are correct; its net
centres remain inside the declared relaxed PSD trace ball.  GE.2--GE.4 have
the correct metric factors, trust-coordinate normalization, hard-edge
modulus, and trust-margin criterion.  Most importantly, the draft correctly
separates the compressible scaling `m=Theta(r/p)` from the anti-pin scaling
`m=r`.  The only recurring terminology qualification is that `d_q` is a
pseudometric on raw matrices and a metric on the induced response quotient.
Also, replacing ambient centres by physically realizable centres doubles
the metric radius and can lose the low-rank property; this does not affect
the stated relaxed theorem but must be budgeted in any later realizability
corollary.

## 1. State space, trace, and metric factors

For actual ports and `J=H/r`, the orthogonal projectors
`P^+-=(I+-J)/2` give

```math
K^\pm_{ij}={w_i^TP^\pm w_j\over n}.                    \tag{AGE.1}
```

Thus both sectors are PSD.  Since every port has squared norm `n`,

```math
\operatorname{tr}K^++\operatorname{tr}K^-
=\sum_i{\|w_i\|^2\over n}=p.                         \tag{AGE.2}
```

In particular each trace is at most `p`, so the actual state space really
is contained in GE.1.  The relaxed space drops the fixed diagonal and the
rank constraints imposed by the eigenspace dimensions, which is legitimate
for an upper covering bound.

Because

```math
\Delta G+\Delta R=2\Delta K^+,
\qquad
\Delta G-\Delta R=2\Delta K^-,                       \tag{AGE.3}
```

GE.4 follows with exactly the displayed factor two.  As in the collective
packing note, `d_q` is only a pseudometric before quotienting: a traceless
diagonal perturbation has zero Boolean quadratic response.  It is a metric
on the response equivalence classes, and on a fixed-sector-diagonal slice.
All later response inequalities depend only on this quotient, so no proof
uses false positive definiteness.

## 2. GE.1 net cardinality and internal centres

For one sector, retain eigenvalues greater than `eta p/4`.  The trace bound
implies

```math
\operatorname{rank}K_{hi}< {4\over\eta},
\qquad
\|K_{lo}\|_{op}\le {\eta p\over4}.                   \tag{AGE.4}
```

The weaker integer statement `rank<=floor(4/eta)` in GE.7 is therefore
valid, including when `4/eta` is an integer.  Since a Boolean query has norm
`sqrt(p)`,

```math
q_p(K_{lo},0)le {p(\eta p/4)\over p^2}={\eta\over4}. \tag{AGE.5}
```

Put `r_eta=ceil(4/eta)` and factor `K_hi=BB^T` after zero padding.  The
factor lies in the Euclidean ball of radius `sqrt(p)` in dimension
`d=pr_eta`.  A maximal internal net of radius
`rho=eta sqrt(p)/8` has cardinality at most

```math
\left(1+{2\sqrt p\over\rho}\right)^d
=\left(1+{16\over\eta}\right)^{pr_\eta}.             \tag{AGE.6}
```

“Internal” matters: every selected `C` still obeys `||C||_F<=sqrt(p)`.
Consequently `CC^T` is PSD, has trace at most `p`, and rank at most
`r_eta`.  Moreover,

```math
\|BB^T-CC^T\|_{op}
\le(\|B\|_{op}+\|C\|_{op})\|B-C\|_{op}
\le {\eta p\over4}.                                  \tag{AGE.7}
```

The factor approximation and discarded tail each cost `eta/4` in `q_p`,
so one sector is covered to `eta/2`.  The Cartesian product of the two nets
has cardinality

```math
\left(1+{16\over\eta}\right)^{2pr_\eta},             \tag{AGE.8}
```

and GE.4 turns its sector radius into `d_q` radius `eta`.  Taking logs gives
GE.5 exactly.

The trace hypothesis is essential to this argument twice: it bounds both
the retained rank and the factor-ball radius.  Without a trace or comparable
nuclear bound, the stated dimension-free truncation is false.

If one wants centres in a smaller realizable subset, selecting one subset
point from every ambient ball that meets the subset gives pairwise distance
at most `2eta` by the triangle inequality.  This standard observation is
correct.  The selected point need not retain the ambient centre's low rank,
and a response theorem demanding radius `eta` must construct the ambient
cover at radius `eta/2` first.  GE.1 itself makes neither stronger claim.

## 3. Trust-coordinate normalization

For a fixed channel, set

```math
a=g+\sigma h,
\qquad b=g-\sigma h.                                  \tag{AGE.9}
```

These are twice the dangerous and safe sector queries and are nonnegative.
With `t=2alpha-1`,

```math
2\alpha g+\sigma h
={(t+2)a+tb\over2},
\qquad
4\alpha^2-1=t(t+2).                                   \tag{AGE.10}
```

Substitution in SA.19 produces

```math
{1\over2}+{t\over2}+{\mu^2a\over4t}
 +{\mu^2b\over4(t+2)},                                \tag{AGE.11}
```

so GE.14 has no missing factor.  The canonical verifier compares the
`alpha` and `t` optimizations on 100 random channels and passes.

The draft explicitly repairs the formal hypothesis issue found in the
collective-packing note: the trust dual identity uses no assumption
`2m>r` and is valid for every `m>=0`; only the original anti-pin compiler
needed that inequality.

## 4. Global modulus and hard example

For GE.2, changing `B` costs at most `|Delta B|/8`.  If
`A'=A+delta`, evaluate the primed objective at `t+s`, where `t` is an
arbitrarily accurate old minimizer.  The old `B` term decreases and

```math
{A+\delta\over4(t+s)}-{A\over4t}
\le {\delta\over4s}.                                  \tag{AGE.12}
```

The total additional cost is `s/2+delta/(4s)`, minimized at
`s=sqrt(delta/2)` with value `sqrt(delta/2)`.  Swapping the two states gives
the absolute estimate.  At `B=0`, direct minimization gives
`F(A,0)-F(0,0)=sqrt(A/2)`, so the exponent and constant are sharp on the
relaxed cone.

If `d_q<=eta`, both channel coordinates in AGE.9 change by at most
`eta p^2`.  Multiplication by `mu^2` makes the dangerous difference
`c^2eta`, where `c=mu p`, and GE.2 yields exactly

```math
c\sqrt{\eta/2}+c^2\eta/8.                             \tag{AGE.13}
```

The hard state in GE.22 has trace `eta p/2<=p`, query distance `eta`, and
all-positive dangerous coordinate `mu^2a=c^2eta`; hence it attains the
square-root term.  It is a relaxed PSD state rather than an actual
fixed-diagonal port state, as the draft correctly says.

## 5. Trust-margin criterion and constants

The derivative of the objective in GE.16 is

```math
{1\over2}-{A\over4t^2}-{B\over4(t+2)^2}.              \tag{AGE.14}
```

It is increasing because the objective is convex.  Therefore the minimizer
is at least `tau` exactly when the derivative at `tau` is nonpositive, which
is GE.23:

```math
{A\over\tau^2}+{B\over(\tau+2)^2}\ge2.               \tag{AGE.15}
```

The sufficient condition `A>=2tau^2` is correct.  If both compared channel
minimizers are at least `tau`, testing each objective at the other's
minimizer bounds its two coordinate sensitivities by
`1/(4tau)` and `1/[4(tau+2)]`.  Coordinate changes at most `c^2eta` give
GE.24.  Substituting `tau=2gamma` gives GE.25 with the displayed factor
`1/8`.  The canonical randomized derivative/criterion checks pass.

The hypothesis must cover every channel that can become an outer maximizer
in either state; a margin only on the currently winning channel of one state
is not stable under comparison.  The draft's all-channel formulation safely
ensures this.

## 6. Critical scaling and scope

At fixed total repeated-port mass,

```math
c={mp\over r}=O(1)
\quad\Longleftrightarrow\quad
m=O(r/p).                                             \tag{AGE.16}
```

The `d_q` cover at radius `epsilon^2` then gives response error below
`epsilon` and logarithmic state count

```math
O\left({p\over\epsilon^2}\log{1\over\epsilon}\right).\tag{AGE.17}
```

For the integral choice `m=floor(r/p)` and `p<=r`, the auxiliary order
`pm` is at most `r`; completing it costs `O(r^2)=O(n)`, which is negligible
relative to `rn=n^(3/2)`.  These normalizations are correct.

This is not the anti-pin regime.  With the original choice `m=r`, one has
`c=p`, so a fixed response error near the hard boundary requires roughly
`eta=O(epsilon^2/p^2)`, not fixed metric accuracy.  Feeding this scale into
GE.5 gives a much larger bound (of order
`p^3 epsilon^(-2) log(p/epsilon)` in logarithmic size), not `O(p)`.
Accordingly GE.1 alone does not compress the growing-port anti-pin response.
The draft states this limitation clearly and does not claim a Boolean
integrality-gap or congruence theorem.

The response-cover centres in GE.5 are internal to the relaxed PSD state
space, but not necessarily actual Boolean-port Gram pairs.  A physically
realizable response net requires the doubled-radius selection argument (and
retuned constants) or a separate rounding theorem.

## 7. Verifier assessment

The canonical verifier passes all 727 diagnostics.  It correctly checks:

- GE.3/GE.4's factor two;
- truncation rank, tail operator norm, and query error;
- the `alpha`-to-`t` trust formula;
- GE.2 and its exact hard edge;
- GE.3 on random PSD sector pairs;
- the derivative form of GE.23;
- the sharp relaxed state GE.22.

It does not construct the Euclidean factor net itself; the volume estimate
and internal-centre argument are deterministic proofs audited in Section 2.
No experimental claim substitutes for those steps.

Run:

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_cross_gram_response_metric_entropy.py
```

It reports `cross-Gram response metric-entropy checks passed: 727`.
