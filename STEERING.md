# Strategic steering

Evidence cutoff: row-sign augmented-greedy/local-stability checkpoint (2026-08-13).
Status: **joint same-switch campaign active; diffuse unmatched core is leading**.

## User-stated objective and workflow directives

The research objective is to determine whether `M_n/n^(3/2)` converges. The
conjectural value `1/2` is not an additional user objective.

The user authorized sustained computational--composition research and a
second three-agent phase from commit `eec5aed`. Reproducible artifacts,
claim classification, regular integration, Git checkpoints, and the existing
consolidation/stopping discipline are workflow directives. They do not make
any mathematical route or ranking a user directive.

The user explicitly authorized a focused campaign seeking one rigorous joint
same-switch inequality at the correct leading scale, with its minimal theorem
formulated before computation. Scalar atom decompositions, separately paid
channels, ordinary polarization, and canonical same-map Gaussian responses
are excluded. Every candidate must be rejected if it reconstructs the full
parent optimization, has a fixed leading loss, or fails the conference and
exact-small-minimizer tests. These are workflow and falsification directives;
mathematical statements and route judgments remain agent-authored.

## Agent-authored campaign assessment

The rigorous frontier is unchanged:

```math
0.336493364431\ldots
\le \liminf_{n\to\infty}\frac{M_n}{n^{3/2}}
\le \limsup_{n\to\infty}\frac{M_n}{n^{3/2}}
\le \frac12.                                      \tag{S1}
```

The requested escape now exists. Take `X` uniform and
`Y=sign(AX)`. Its doubled asymmetric response is signing-independent:

```math
\mathbb E[X^{\mathsf T}AY]
=\mathbb E\|AX\|_1
=\left(\sqrt{2/\pi}+o(1)\right)n^{3/2}.          \tag{S2}
```

Agreement/disagreement recoupling can retain the whole cross field by
collapsing the anchored shore to one weighted vertex. On that weighted shore,
deterministic best-improvement one-spin ascent gives a rigorous samplewise
certificate in `O(n^3)` time. If `Delta_gr(A,X)` is its clipped loss, then

```math
Q(A)\ge\mathbb E\|AX\|_1-\mathbb E\Delta_{\rm gr}(A,X). \tag{S3}
```

This escapes both earlier no-go theorems: it is asymmetric, uses one common
parent spin, preserves all shore/cross-field cancellation before absolute
values, and neither decomposes into paid scalar channels nor evaluates the
full Boolean maximum.

All tested exact cases from order 6 through 14 give a doubled normalized
certificate above `0.672986...` (orders 3 and 5 are small-order exceptions).
Conference samples give defects `0.001278` at order 30, `0.000124` at 54,
`0.000009` at 90, and zero at 98. These are a computational scaling law, not
an asymptotic proof. The spectral and one-shot implementations miss the
coefficient on tested larger conference orders and are inactive.

## Exact sufficient successor lemma

Let `c_*=0.672986728863...` be the current doubled lower coefficient. A new
rigorous bound follows from the exact uniform statement

```math
\mathbb E_X\Delta_{\rm gr}(A,X)
\le\left(\sqrt{2/\pi}-c_*-\eta\right)n^{3/2}
+o(n^{3/2})                                      \tag{S4}
```

for some `eta>0`, uniformly over signings with project-scale cap. The stronger
statement `E Delta_gr=o(n^(3/2))` would improve the project lower constant to
`1/sqrt(2pi)=0.398942...`. It would improve the interval, not alone prove
convergence.

The defect now has a path-independent potential reduction. At a terminal set
on either shore it is twice the positive part of unmatched original row-field
mass minus terminal removal margin. Let `kappa_I,kappa_J` be the least numbers
of largest outside fields needed to carry each unmatched mass, and put
`kappa_*=min(kappa_I,kappa_J)` in the hard branch and zero otherwise. For
`Q(A)<=K n^(3/2)`,

```math
{\mathbb E\Delta_{\rm gr}\over n^{3/2}}
\le C\alpha\sqrt{\log(2e/\alpha)}
+K\Pr\{\kappa_*>\alpha n\}.                     \tag{S5}
```

Thus the exact leading target is `kappa_*/n -> 0` in probability. This is
strictly more specific than the original cap problem: it concerns the output
of two deterministic polynomial local searches. Its exact falsifier is a
nonvanishing probability of simultaneous linear-sized diffuse unmatched
cores on both shores; isolated large fields or one bad shore cannot defeat it.
An explicit regular-Hadamard/Seidel family has a bad stable core at one
row-sign realization, so project-scale norm plus terminal stability alone is
insufficient. Its possible normalized defect tends to `1/(4sqrt(2))`, but
nonvanishing basin probability under random `X` is not known.

The cleanest convergence interface remains a family of explicit constructors,
defined without `M_n`, such that for every child
signing `A` of order `n` and every `k>=2`,

```math
\operatorname{cap}(\mathcal T_k(A))
\le (1+\eta_k)k^{3/2}\operatorname{cap}(A)
   +Cnk^{3/2},
\qquad \eta_k\longrightarrow0.                    \tag{S6}
```

with `C` uniform and `T_k(A)` a signing of order `kn`. Applying (S6) to a
large near-liminf seed, sending `k` to infinity, and filling a fixed remainder
with `O(nN)` new-edge cost would give `limsup<=liminf`, hence convergence.

The signed elliptope has the right factor but the sharp conference floor;
fixed-level SOS inherits its child obstruction, while scalar/separately paid
channels have a fixed gap. These are inactive without a new global separator
or an inexact quotient with summable error.

## Ranked alternatives and falsification criteria

1. **Diffuse-core theorem.** Prove (S5)'s probability term vanishes, or the
   weaker coefficient bound (S4). Stop only upon a scalable simultaneous-core
   construction or a proof that the terminal condition hides full cap
   optimization; current evaluation is explicitly polynomial.
2. **Boolean joint composition.** Find a global separator or inexact quotient
   and prove (S6). Signed elliptope, fixed-level duplication-closed SDP/SOS,
   bounded-local facets, and scalar atoms are stopped.
3. **Genuine nonconvergence.** It requires a fixed `epsilon>0` and two
   asymptotically separated subsequences, or strict `liminf<limsup`. No
   candidate subsequences are known.

## Checkpoint decision

The user-requested milestone is met by (S3), now with an explicit polynomial
certificate and the correct raw leading coefficient. The scaling experiment
also isolates the precise uniform lemma (S5). No rigorous asymptotic bound,
recurrence, or convergence mechanism has yet improved, so this is material
campaign progress but does not reset the README's proof-level stopping count.

Continue only on the diffuse-core theorem and its falsification criterion; do
not reopen the deprioritized spectral/one-shot or proved scalar-atom/SOS
implementations. Composition remains secondary because it could prove convergence.
Wave 61 is the next blank-slate boundary if ordinary waves resume.
