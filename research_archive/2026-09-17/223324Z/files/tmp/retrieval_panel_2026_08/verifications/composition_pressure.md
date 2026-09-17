# Archive verification: composition and pressure proposals

## Verdict

I find **no A or B candidate** in the four assigned reports.  There are eight
`C` reformulations of archived obligations and two `D` architectures with a
proved wrong-scale/fixed-loss obstruction.

| Frozen architecture | Class | Exact archive alignment |
|---|:---:|---|
| Extremal `L_IC` | C | Ledger 10.115.2--10.115.3 and 10.140.1: summable `2/3`-power composition plus the exact full bridge objective |
| Extremal hereditary thinning (3) | C | Ledger 10.113.2: fixed-size thinning is induced-pattern universality followed by the desired scalar inequality; proportional thinning is the archived adaptive-selector game |
| Extremal finite-template random switching | D | Fixed templates are variance-driven and forget `M_k`; the growing-template regime needs the full local Laplace/coset histogram and lacks a uniform absorber |
| Discrepancy `L_advGT` | C | The opposite-direction Hammersley recurrence for the already archived minimized pressure-limit obligation; no deterministic interpolation mechanism is supplied |
| Discrepancy star statement (S) | C | Ledger 10.44, especially (10.254)--(10.256): the same derivative-scale one-vertex insertion target, whose exact state is the full near-cap deficit hierarchy |
| Discrepancy vector/factorization route | D | Exact `sqrt(n)` integrality gap already on `W_n`, with any same-spin return requiring the archived fixed-loss polarization/recoupling step |
| Spin-glass `L_Lap` | C | The exact fixed-temperature bridge minimization already isolated in 10.131.3--10.131.4, 10.135.3, and `two_temperature_bridge_audit.md` |
| Contrarian `FT` | C | Exactly the archived statement that minimized antipodal pressure converges at every fixed inverse temperature |
| Contrarian `PG` | C | A different homogeneous scalarization of the same exact block-completion obligation; not falsified by the archived Pythagorean *lower*-bound counterexample |
| Contrarian `CBL` | C | `dependent_profile_recovery.md`, (14)--(15): the same arbitrary dependent compressed-lift absorption/equality obligation, now rounded to every fibre size |

No item warrants specialist revision under the A/B rule.

## Common normalization and implication audit

Put

```math
\overline Z_n(A,\beta)=2^{-n}\sum_x
 \cosh\!\left({\beta H_A(x)\over\sqrt n}\right),
\qquad
\mathsf F_n(\beta)=\min_A\log\overline Z_n(A,\beta).
```

The discrepancy report's average over `(sigma,x)`, the spin report's
`widehat Z`, and the contrarian's codeword average are this same object.
The augmented code has `2^n` words (for the relevant `n>=3`), each with two
`(sigma,x)` preimages.  Hence there is no missing factor two, and

```math
\beta {M_n\over n^{3/2}}-\log2
\le {\mathsf F_n(\beta)\over n}
\le \beta {M_n\over n^{3/2}}.                 \tag{V.1}
```

Thus convergence of `mathsf F_n(beta)/n` for every fixed `beta>0`, followed
by `beta -> infinity`, does imply convergence of `M_n/n^(3/2)`.  No
uniformity as `beta -> infinity` is needed.  The pressure reductions in all
three reports have the correct normalization and order of limits.

I checked the current primary text of Füredi--Ruzsa, Theorem 5.  The needed
hypotheses are a real all-integer sequence, fixed `mu>1` and threshold,
nonnegative nondecreasing `f`, near-subadditivity on
`N<=m<=n<=mu m`, and `sum f(r)/r^2<infinity`.  Both
`f(r)=Kr/log^2(e+r)` (eventually monotone) and
`f(r)=C_beta r^(1-delta_beta)` satisfy them.  Applying the theorem to
`-mathsf F_n(beta)` also validates the direction used in `L_advGT`.
Hammersley's discrete consequence for a power defect is therefore valid
even if that older citation is stated loosely.  Bare `o(r)` is not enough.

The pressure directions are different but consistent:

- `L_advGT` asks for a **lower** bound
  `mathsf F_(m+n)>=mathsf F_m+mathsf F_n-o(m+n)` and applies the theorem to
  `-mathsf F`.
- `L_Lap` asks for an **upper** bound supplied by one binary completion.
  Random bridge averaging naturally gives this direction, but at contracted
  child temperatures and with a linear annealed bridge term.

Neither direction is inherited from Guerra--Toninelli: its disorder average
and Gaussian integration-by-parts identity do not survive the outer
deterministic minimum.  The reports correctly acknowledge this unmatched
hypothesis.

## Exact architecture checks

### 1. Extremal proposals

**`L_IC` (C).**  Its block identity

```math
Q\!\begin{pmatrix}A&R\\R^T&B\end{pmatrix}
=\max_{x,y}\{|H_A(x)+H_B(y)|+|x^TRy|\}
```

is exact, as is the Füredi--Ruzsa implication for
`b_n=M_n^(2/3)`.  But the quantifiers align with the archived obligation:
for **every** comparable low-cap child pair, choose one of `2^(mn)` bridges
that satisfies all `2^(m+n)` joint state constraints with a summable
transformed defect.  Ledger 10.115.2--10.115.3 already states precisely this
`2/3`-power target and exact bridge state; 10.140.1 rejects the universal
completion reformulation as non-reducing.  Excluding high-cap inputs does
not compress the low-cap response.  The claimed non-recovery of high-cap
histograms is therefore not a strict reduction of the operative bridge
state.

**Hereditary thinning (C).**  Align the quantifiers by taking a liminf
sequence of competitive order-`N` parents and then fixing `n`.  The archived
rectangle estimate makes every such parent cut-quasirandom, so every fixed
signed `n`-vertex pattern, including an order-`n` minimizer, occurs for all
large `N`.  Consequently

```math
\min_{|S|=n}Q(A_N[S])=M_n,
```

and (3) becomes exactly
`M_n/n^(3/2)<=liminf M_N/N^(3/2)+r(n)`.  If `n/N` stays positive instead,
the selection retains the full revealed-versus-hidden restriction profile.
This is the collision proved in ledger 10.113.2, not a new hereditary
mechanism.

**Finite-template random switching (D as stated).**  The report's variance
calculation is correct.  With `K_k` blocks in an edge decomposition, the
total variance for each global spin is `E_N=Theta(N^2)`.  A Bernstein union
bound at `Theta(N^(3/2))` can pay `2^N` tests from cap-only information only
when `k=O(N^(1/3))`; for fixed `k` its leading value is universal and forgets
`M_k`.  Once `k` is large enough for the local optimum to matter, the needed
log-mgf is the complete local energy/coset histogram, and the design theorem
is not uniform in that growing template.  Fixed-template absorption itself
costs only `O_k(N)` and is harmless, but it does not cure this fixed leading
probabilistic loss.  (A genuinely lossless all-order propagation of selected
liminf templates would in fact suffice for convergence; the report's claim
that an upper construction can never suffice is too strong.)

### 2. Discrepancy proposals

**`L_advGT` (C).**  The implication and the superadditive sign are correct.
Nevertheless the lemma is simply a Hammersley-summable inequality for the
same scalar minimized pressure whose fixed-temperature convergence is the
archived target.  The imported quenched theorem has the min/expectation order
wrong and supplies none of the inequality.  The abstract scalar-pressure
countermodel in ledger 10.135.2 shows that the known centering, convexity,
restriction, and entropy axioms cannot imply it.  This is a desired pressure
recurrence, not a new interpolation state.

**Star statement (S) (C).**  The exact identity and its implication are
correct.  It is the archived insertion functional

```math
\Delta_A(b)=\max_x\{|b\cdot x|-[M(A)-|H_A(x)|]\}.
```

Statement (S) is a slightly stronger version of (10.256), with summable
remainder `Cn^(1/2-delta)`.  Controlling it requires the entire exponential
near-cap slack profile; endpoint-only balancing is defeated by the thick-cap
example.  That example is not an asymptotic minimizer falsifier, so the
correct label is C rather than D.

**Vector/factorization route (D).**  Orthogonal edge vectors give
`vdisc(W_n)<=sqrt(E_n)=Theta(n)`, while the proved lower frontier gives
`disc(W_n)=Omega(n^(3/2))`.  This is an exact `Omega(sqrt n)` gap on the
target matrix, before rounding.  Independent left/right or bilinear
factorizations additionally encounter `Q<=B<=2Q` and the archived
same-spin recoupling obligation.  The proposed surrogate is rigorously at
the wrong scale.

### 3. Spin-glass proposal

**`L_Lap` (C, not B).**  For fixed `beta,A,B`, define

```math
\Delta_\beta(A,B)=\min_{\eta,D}
 [\log\overline Z_{m+n}(A\oplus_D\eta B,\beta)
  -\log\overline Z_m(A,\beta)-\log\overline Z_n(B,\beta)].
```

Then `L_Lap` is exactly `sup_(A,B) Delta_beta(A,B)=O_beta(N^(1-delta))`.
This is the full binary Gibbs bridge completion: the minimization ranges over
`2^(mn+1)` outputs and each value sums all parent spin states.  The archive's
noisy-code, reverse-KL, and Rényi reveal formulas are exact rewritings of
this same state; without a closure theorem, sequential evaluation is full
backward dynamic programming.

The order-six equal-Laplace-moment example in the report is correct (I
independently reproduced the histograms and `beta_*=0.8975636801...`).  It
proves that **one scalar moment** does not determine a histogram.  It does
not make the universally quantified completion lemma a compressed bridge
mechanism.  Allowing `(eta,D)` to depend on `beta` also does not create a new
strict reduction: every archived fixed-temperature bridge criterion already
allows a new witness at each `beta`.

Conference results do not justify D.  They prove a positive linear defect
for uniform-output reverse KL, fixed small Rényi tilts, fixed quantiles, and
polynomial sampling in a strict high-temperature interval.  They do **not**
exclude an exponentially rare best bridge, and `L_Lap` is precisely allowed
to choose that bridge.  Conversely, those results explain why no cited
spin-glass theorem proves `L_Lap`, especially at every fixed temperature.

### 4. Contrarian-short proposals

**`FT` (C).**  This is exactly the conclusion of the archived
finite-temperature program, with optional local uniformity.  Pointwise
convergence at every fixed `beta` already suffices by (V.1).  Naming the
thermodynamic limit is not a comparison theorem or information compression.

**`PG` (C).**  Its normalization is correct:
`lambda_n=M_n^2/n^2=n(M_n/n^(3/2))^2`, so a summably almost-subadditive
`lambda_n` would force convergence.  For children with normalized caps
`c_A,c_B` and `theta=n/(n+m)`, its zero-defect target is the power mean

```math
{Q(parent)\over(n+m)^{3/2}}
\le [\theta c_A^2+(1-\theta)c_B^2]^{1/2}.
```

This differs algebraically from the archived `2/3`-power and linearized
targets for unequal child constants, but it asks for the identical hard
object: one bridge controlling the exact joint parent maximum for every
low-cap pair.  It coincides with them when the child normalized constants
agree, the asymptotically decisive regime.  The archived
`pythagorean_centered_width_block.md` falsifies a Pythagorean **lower** gain
from scalar cross norm; it does not falsify this existential upper
completion.  Thus PG is C, not D, but the squared normalization supplies no
new mechanism and is not B.

**`CBL` (C).**  The parity rounding `r_t=t^(3/2)+O(1)`, coarse-state lower
bound, fixed-seed order of limits, and padding argument are all correct.
However the lemma is the existing arbitrary dependent compressed-lift
target (after cancelling the harmless factor two in the older artifact's
`x^T A x` convention): prescribe critical block sums and minimize `Q` over
the entire affine fibre, demanding

```math
G_t(A)\le t^{3/2}Q(A)+o((kt)^{3/2}).
```

The rounding extends the exact-square formulation to every `t` but does not
reduce its state.  The fibre contains exponentially many signings and the
cap tests exponentially many microscopic spins.  The Frobenius/ANOVA,
independent random-lift, Hadamard, and one-step Onsager no-gos explain why
separately paying the forced microscopic residual has a leading loss; they
do not falsify arbitrary seed-dependent cancellation.  The archive already
records that exact dependent lifts can absorb finitely and that the only
surviving theorem is normalized equality on centered/chiral profiles.
Accordingly CBL is an unresolved archived obligation (C), not a rigorously
false architecture (D).

## Final pressure/composition boundary

The proposals correctly avoid factor-of-two mistakes and, except for the
noted overstatement about selected upper constructions, their conditional
implications are sound.  What none supplies is a new theorem between the
input data and the needed recurrence.  At the target scale the surviving
statements require one of: the full state-dependent bridge maximum, the full
Gibbs bridge output, the exponential near-cap/coset profile, or the full
dependent affine lift fibre.  Random, separately paid, bounded-profile, and
strict-high-temperature substitutes have the archived fixed-loss or
conference obstruction.  Therefore the exact archive classification is
`C x 8, D x 2`, with no candidate eligible for an A/B revision.
