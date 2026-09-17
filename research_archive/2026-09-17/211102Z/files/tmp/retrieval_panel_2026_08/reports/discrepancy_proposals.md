# Discrepancy / vector-balancing domain report

**Proposal freeze:** 2026-08-15, before any archive exposure. I used only the
naked problem, the discrepancy retrieval packet, the report schema, and the
primary literature cited below.

## Executive verdict

There is one clean candidate that preserves the problem exactly: prove an
almost-superadditive thermodynamic limit for the **minimum two-sided soft
discrepancy** of the degree-two Walsh matrix.  The zero-temperature transfer
then loses only `log(2)/beta`, which can be sent to zero *after* the size limit.
This is the best proposal below.

There is also an exact one-vertex recurrence that turns the issue into a
weighted linear discrepancy problem for the near-maximal level sets of one
optimal quadratic form.  It is more recognizably a partial-coloring problem,
but the increment estimate needed for convergence is extremely sharp.

Generic vector discrepancy, `gamma_2`, Hadamard factorization, and independent
bipartite switching do **not** presently give a viable third route.  The most
natural vector relaxation is wrong by a factor of order `sqrt(n)` on this very
matrix, before any rounding.  No abundance theorem in the packet survives the
outer deterministic minimization and all `2^n` tests with zero leading loss.

## 1. Exact native translation

Put `E_n = binom(n,2)` and

\[
  X_n=\{-1,1\}^n/(x\sim -x).
\]

Let `W_n` be the `2^(n-1) by E_n` degree-two Walsh matrix

\[
  W_n(x,\{i,j\})=x_i x_j,\qquad x\in X_n,\quad i<j.
\]

For an edge coloring `a in {+1,-1}^{E_n}`,

\[
  (W_na)(x)=\sum_{i<j}a_{ij}x_ix_j=H_a(x).
\]

Consequently, in standard matrix-discrepancy notation,

\[
  \boxed{M_n=\operatorname{disc}(W_n)
  :=\min_{a\in\{\pm1\}^{E_n}}\|W_na\|_\infty.}
\]

Thus the quadratic appearance is in the **row indexing**; the outer variable
is an ordinary coloring of the `E_n` columns.  The rows are highly structured,
not an arbitrary family of `2^(n-1)` constraints.  Character orthogonality
also gives the exact identity

\[
  2^{-(n-1)}W_n^{\mathsf T}W_n=I_{E_n}.
\]

The code translation checks as follows.  Write

\[
  w_x=(x_ix_j)_{i<j},\qquad
  \mathcal C_n^+=\{\sigma w_x:\sigma=\pm1, x\in X_n\}.
\]

For sign vectors, Hamming distance satisfies

\[
  d(a,\sigma w_x)=\frac{E_n-\sigma\langle a,w_x\rangle}{2}.
\]

Since `|C_n^+|=2^n` for `n>=3`,

\[
  d(a,\mathcal C_n^+)=\frac{E_n-Q(a)}2,\qquad
  Q(a)=E_n-2d(a,\mathcal C_n^+),
\]

and hence

\[
  M_n=E_n-2\rho(\mathcal C_n^+).
\]

Three nearby translations are not interchangeable:

* If `ell(a)` is the minimum number of negative edges after vertex
  switching, then

  \[
    \ell(a)=\frac{E_n-\max_x H_a(x)}2,\qquad
    Q(a)=E_n-2\min\{\ell(a),\ell(-a)\}.
  \]

  Covering radius for the unaugmented cut code controls only `ell(a)`.

* If `g_a(S)=sum_{ij in binom(S,2)} a_ij`, `T=sum_{i<j}a_ij`, and
  `d_i=sum_{j ne i}a_ij`, then for `S={i:x_i=1}`,

  \[
    4g_a(S)=T+\sum_i d_i x_i+H_a(x).
  \]

  The constant and degree channels are of the same potentially relevant
  order, so induced-set discrepancy has no sharp-constant transfer for free.

* For the symmetric hollow matrix `A`,

  \[
    x^{\mathsf T}Ax=2H_a(x),\qquad
    \max_{r,c}|r^{\mathsf T}Ac|\ge 2Q(a).
  \]

  Independent `r,c` are a relaxation.  Polarization or Grothendieck-type
  comparison with any fixed factor greater than one cannot identify a limiting
  constant for `M_n/n^(3/2)`.

## 2. Frozen architectures (at most three)

### A. Two-sided soft discrepancy and adversarial interpolation — best

Replace the maximum over all antipodal cut words by its log-sum-exp, but keep
the minimum over edge colorings outside.  Seek a deterministic analogue of
Guerra--Toninelli interpolation with a sublinear additivity defect.  This
handles all `2^n` tests in the partition function itself and retains one
global orientation sign; it does not union-bound tests or replace same-spin
configurations by independent bipartite signs.

The exact missing lemma is in Section 4.

### B. Near-extremal-profile star balancing

If `a` colors `K_n` and `b_i` colors the edge from vertex `i` to a new vertex,
then for `t in {+1,-1}`

\[
  H_{a\oplus b}(x,t)=H_a(x)+t\langle b,x\rangle.
\]

The elementary identity

\[
  \max_{t=\pm1}|u+tv|=|u|+|v|
\]

therefore gives the exact recurrence

\[
  \boxed{M_{n+1}=\min_{a\in\{\pm1\}^{E_n}}
  \min_{b\in\{\pm1\}^n}\max_{x\in X_n}
  \bigl(|H_a(x)|+|\langle b,x\rangle|\bigr).}
\]

A sufficient missing statement for convergence would be: there are
`delta>0`, `C<infinity`, and `n_0` such that for every `n>=n_0` some exact
minimizer `a` with `Q(a)=M_n` and some `b in {+1,-1}^n` satisfy

\[
  \tag{S}
  \max_x\bigl(|H_a(x)|+|\langle b,x\rangle|\bigr)
  \le M_n+\frac{3M_n}{2n}+C n^{1/2-\delta}.
\]

Indeed, convexity gives `(n+1)^(3/2) >= n^(3/2)(1+3/(2n))`, so (S) implies

\[
  \frac{M_{n+1}}{(n+1)^{3/2}}
  \le \frac{M_n}{n^{3/2}}+\frac{C}{n^{1+\delta}}.
\]

The total upward variation is summable; boundedness below then forces the
normalized sequence to converge.

For a proposed partial-coloring proof, the constraints are exactly

\[
  |\langle b,x\rangle|\le R-|H_a(x)|\quad\text{for every }x\in X_n,
  \qquad R=M_n+\frac{3M_n}{2n}+Cn^{1/2-\delta}.
\]

At a Lovett--Meka stage with uncolored coordinate set `U`, `|U|=u`, one must
allocate thresholds `s_{U,x}` satisfying

\[
  \sum_{x\in X_n}\exp\!\left(-\frac{s_{U,x}^2}{16u}\right)\le\frac{u}{16}
\]

and ensure that the thresholds accumulated over all adaptive stages are at
most `R-|H_a(x)|` for every row.  This is how entropy would have to pay for
the complete exponential family.  Merely asserting that most states have
large slack, or that there are many good edge colorings, is insufficient.
The deterministic outer minimum survives only if the argument jointly
selects the exact (or adequately near-) optimal `a` and the star coloring
`b`; an average over random `a` does not do this.

This architecture is exact but less credible than A.  The supplied anchors
already show why an asymptotic theorem is essential: `M_6=5` and `M_7=9`, so
even the unrestricted one-step increment at `n=6` exceeds `3M_6/(2n)` by
`2.75`, an amount of order `sqrt(n)`.

### C. Tensor/factorization/vector relaxation — stopped

The natural vector-discrepancy relaxation of the exact Walsh matrix is

\[
  \operatorname{vdisc}(W_n)=
  \min_{\|u_e\|_2=1}\max_{x\in X_n}
  \left\|\sum_e W_n(x,e)u_e\right\|_2.
\]

Taking the `u_e` orthonormal gives

\[
  \operatorname{vdisc}(W_n)\le\sqrt{E_n}=\Theta(n),
\]

whereas the supplied rigorous lower frontier gives
`M_n=Omega(n^(3/2))`.  Thus this relaxation has an integrality gap of order
at least `sqrt(n)` on the target family itself.  Exact tensorization of
`gamma_2` cannot repair a wrong-scale relaxation, and the general
`gamma_2`/hereditary-discrepancy comparisons have logarithmic losses that are
known to be real.  A Hadamard block construction also changes scalar same-spin
variables into vector-valued or independent left/right variables; returning
to one scalar spin per vertex needs an exact transfer theorem, not a constant
factor comparison.  I therefore do not regard this as a viable architecture.

Nor is there a viable nonconvergence proposal here: special matrix orders can
give an upper subsequence, but genuine nonconvergence additionally requires a
strictly larger lower bound on another subsequence.  None of switching,
factorization, or Hadamard availability supplies that second half.

## 3. Best implication

For `beta>0`, define the normalized, explicitly two-sided partition function

\[
  Z_n(a,\beta)=2^{-(n+1)}
  \sum_{x\in\{\pm1\}^n}\sum_{\sigma=\pm1}
  \exp\!\left(\frac{\beta\sigma H_a(x)}{\sqrt n}\right).
\]

Each element of `C_n^+` occurs twice in the displayed sum, so this is exactly
the uniform exponential moment over the `2^n` distinct antipodal codewords.
Set

\[
  F_n(\beta)=\min_{a\in\{\pm1\}^{E_n}}\log Z_n(a,\beta),
  \qquad p_n(\beta)=\frac{F_n(\beta)}n.
\]

The best implication is

\[
\boxed{\text{Hammersley's approximate superadditivity theorem}
\ +\ \boxed{L_{\rm advGT}}
\ +\ \text{the log-sum-exp inequality}
\ \Longrightarrow\ M_n/n^{3/2}\text{ converges}.}
\]

The approximate subadditivity theorem used here is J. M. Hammersley,
“Generalization of the Fundamental Theorem on Subadditive Functions,”
*Proc. Cambridge Philos. Soc.* 58 (1962), 235--238
([primary DOI](https://doi.org/10.1017/S030500410003646X)).  Its summability
condition applies to an error `O(N^alpha)` with `alpha<1`; equivalently one
may apply it to `-F_n`.

## 4. Exact missing lemma `L_advGT`

> **Adversarial two-sided Guerra--Toninelli lemma.** For every finite
> `B>0`, there exist constants `C_B<infinity` and `alpha_B<1` such that, for
> all integers `n,m>=3` and every `beta in (0,B]`,
> \[
>   \boxed{
>   F_{n+m}(\beta)\ge F_n(\beta)+F_m(\beta)
>   -C_B(n+m)^{\alpha_B}.}
> \]
> The same `C_B,alpha_B` must work uniformly for all splits `n+m` and all
> temperatures in `(0,B]`.  The `F_k` in all three terms are the deterministic
> minima over `{+1,-1}` edge colorings and use the single global `sigma`
> channel in the definition above.

The scale is exact.  The supplied discrepancy upper bound gives
`F_n(beta)=O_beta(n)`, and the zero-temperature scale of `F_n` is linear in
`n`; hence an `o(n)` defect is what approximate superadditivity requires.
An `O(n)` defect is useless.  In zero-temperature energy units the allowed
error is `O(n^(alpha_B+1/2))=o(n^(3/2))` for fixed `beta`.

### Proof that the lemma implies convergence

For every fixed `a`, the largest exponential term is
`exp(beta Q(a)/sqrt(n))`.  There are `2^n` distinct codewords.  Therefore

\[
  \frac{\beta Q(a)}{\sqrt n}-n\log2
  \le \log Z_n(a,\beta)
  \le \frac{\beta Q(a)}{\sqrt n}.
\]

Minimizing over `a` and writing `mu_n=M_n/n^(3/2)` gives the factor-checked
squeeze

\[
  \tag{1}
  \frac{p_n(\beta)}\beta
  \le \mu_n
  \le \frac{p_n(\beta)}\beta+\frac{\log2}\beta.
\]

The supplied `O(n^(3/2))` upper bound on `M_n` makes `F_n(beta)=O_beta(n)`;
also `F_n(beta)>=0` because averaging over `sigma` gives a `cosh`, at least
one.  Applying Hammersley's theorem to `-F_n(beta)` and the lemma shows that

\[
  p(\beta)=\lim_{n\to\infty}p_n(\beta)
\]

exists for every fixed `beta`.  Taking `liminf` and `limsup` in (1),

\[
  \frac{p(\beta)}\beta
  \le\liminf_n\mu_n\le\limsup_n\mu_n
  \le\frac{p(\beta)}\beta+\frac{\log2}\beta.
\]

The gap is at most `log(2)/beta` for every `beta`; sending `beta` to infinity
proves equality of `liminf` and `limsup`.  No value of the limiting constant
is required.  Notice the safe order of limits: first `n to infinity` at fixed
temperature, then `beta to infinity`.

## 5. Why `L_advGT` is strictly weaker than full Boolean/coset optimization

The full data at order `n` are the pointwise energy/coset-weight histograms

\[
  N_{a,n}(t)=\#\{c\in\mathcal C_n^+:\langle a,c\rangle=t\}
  \quad\text{for every }a\in\{\pm1\}^{E_n}.
\]

Those data determine every `Z_n(a,beta)`, every optimizer, and `M_n`.
The converse fails for three concrete reasons.

1. For a fixed temperature, `F_n(beta)` is only the lower envelope, over
   `a`, of one scalar Laplace transform
   \[
     \log\left(2^{-n}\sum_tN_{a,n}(t)e^{\beta t/\sqrt n}\right).
   \]
   Any coloring whose transform is never exposed on that lower envelope is
   completely invisible; its histogram can change without changing `F_n`.
   This is an algebraic loss of information, not a description-length claim.

2. The lemma does not even provide the exact lower envelope.  It asserts one
   family of inequalities with an unspecified sublinear defect.  Arbitrary
   changes at finitely many orders are absorbed by enlarging `C_B`, so the
   lemma cannot recover any prescribed finite value `M_k`, let alone its
   minimizing cosets or weight distributions.

3. A proof of the lemma may use only the soft moment at the optimizing
   temperature; it never asks which Boolean state attains `Q(a)`.  The
   zero-temperature squeeze recovers only equality of two asymptotic limits,
   with an error `log(2)/beta`, not any finite optimum or histogram.

Thus full Boolean/coset optimization implies enough information to evaluate
the quantities in the lemma, whereas the lemma is invariant under changes
that full optimization detects.  It is strictly weaker even though it is
strong enough, through (1), to settle convergence.

For architecture B the strictness is sharper still: (S) needs the geometry of
the near-top level sets of only **one** optimal coloring and one feasible star
coloring.  It neither determines `M_{n+1}` (the true optimum may be smaller)
nor accesses the histograms of any other coset.

## 6. Decisive falsifiers

### For the best lemma

Define the exact superadditivity defect

\[
  \Delta_B(n,m)=\sup_{0<\beta\le B}
  \bigl[F_n(\beta)+F_m(\beta)-F_{n+m}(\beta)\bigr].
\]

The route is dead if there are `B,c>0` and a sequence of splits with
`min(n_k,m_k) to infinity` such that

\[
  \Delta_B(n_k,m_k)\ge c(n_k+m_k).
\]

This is a decisive structural falsifier: it contradicts every possible
choice `alpha_B<1`.  At finite order, the test is exact—enumerate edge-signing
orbits, compute each finite energy histogram, take its exponential moment,
and minimize it with certified interval arithmetic in `beta`.  Once a proof
claims concrete `C_B,alpha_B`, a single pair with
`Delta_B>C_B(n+m)^(alpha_B)` refutes that claim.  The most important diagnostic
is to compute the two-sided `F_n` itself, not the one-sided pressure; good
behavior of the latter with a linear defect in the former falsifies the
proposed transfer.

### For star balancing

Let

\[
  R_n^*=\min_{a:Q(a)=M_n}\min_b
  \max_x\bigl(|H_a(x)|+|\langle b,x\rangle|\bigr)
\]

and

\[
  \eta_n=\frac1{\sqrt n}
  \left(R_n^*-M_n-\frac{3M_n}{2n}\right).
\]

Statement (S) says `eta_n=O(n^(-delta))` for some `delta>0`.  Hence a proof
that `limsup eta_n>0` is decisive.  Equivalently, for the danger set

\[
  D_a(s)=\{x:M_n-|H_a(x)|\le s\},
\]

it suffices to show along infinitely many orders that every optimal `a` and
every `b` have an `x in D_a(s_n)` with `|<b,x>|` exceeding the remaining
budget in (S).  This is a finite Hamming-covering computation at each order
and tests the exact same-spin constraints.

### For factorization

The orthonormal-vector construction
`vdisc(W_n)<=sqrt(E_n)` together with `M_n=Omega(n^(3/2))` already falsifies
lossless use of ordinary vector discrepancy.  For any proposed substitute,
compute its value on `W_n`; a ratio bounded away from one (or of the wrong
order) rules out a leading-constant transfer before tensorization.

## 7. Closest imported theorem and every unmatched hypothesis

The closest theorem is Francesco Guerra and Fabio L. Toninelli,
“The Thermodynamic Limit in Mean Field Spin Glass Models,” *Communications in
Mathematical Physics* 230 (2002), 71--79
([primary arXiv](https://arxiv.org/abs/cond-mat/0204280),
[DOI](https://doi.org/10.1007/s00220-002-0699-y)).  For the SK Hamiltonian
with independent centered unit Gaussian couplings, their Theorem 1 proves
superadditivity of the quenched expected log partition function by Gaussian
interpolation and convexity of the squared replica overlap; Theorem 2 gives
the thermodynamic limit.  Their Section 4.2 extends existence to symmetric
non-Gaussian disorder with finite fourth moment, explicitly including
Rademacher couplings, with interpolation error of sublinear total order.

The unmatched hypotheses are all substantive:

1. **Expectation versus adversarial minimum.** They control
   `E_J log Z(J)`.  We need `min_a log Z(a)`.  The elementary inequality
   `min <= expectation` has the wrong direction for the lower bound in
   `L_advGT`.

2. **Independent block disorders.** Their interpolation replaces one random
   system by two independently resampled block systems.  Restrictions of one
   deterministic minimizing signing are neither independent nor guaranteed
   to minimize their block pressures.

3. **Integration by parts.** Gaussian integration by parts (or its averaged
   non-Gaussian replacement) supplies the signed derivative.  A discrete
   minimizer over `{+1,-1}` couplings has no corresponding identity; local
   edge-flip optimality alone does not reproduce the replica-overlap formula.

4. **Convex covariance.** Standard SK covariance is a convex function
   `q^2` of one overlap.  The exact two-sided augmentation has configurations
   `(sigma,x)` and covariance proportional to `sigma*tau*q(x,y)^2`, which is
   not the same convex-overlap structure.

5. **One channel versus two.** The standard theorem treats one partition
   function.  Replacing our global `sigma` by independent orientation signs
   in the two blocks makes the partition function factor, but changes the
   problem.  The two orientations can have exponentially different weights,
   so this cannot be dismissed without an exact comparison.

6. **Uniform deterministic error.** We need an `o(n+m)` defect uniformly over
   every split and bounded temperature interval.  Concentration around a
   quenched mean does not provide a single signing that simultaneously
   satisfies this family of comparisons.

7. **Outer-minimum stability.** Any abundance or entropy argument must
   produce a family surviving the minimization at all three sizes.  Existing
   high-entropy coloring samplers concern a prescribed linear system and do
   not give this stability.

For architecture B, the closest imported result is the Lovett--Meka partial
coloring theorem
([primary arXiv](https://arxiv.org/abs/1203.5747),
[DOI](https://doi.org/10.1137/130929400)).  Missing are: a near-extremal-level
entropy bound for some optimal `a`; a multistage budget valid for the adaptive
uncolored coordinate sets; no accumulated leading loss; and the sharp
increment `3M_n/(2n)+o(sqrt(n))` rather than a generic `O(sqrt(n))` with an
uncontrolled constant.

For architecture C, the closest exact tensor theorem is
Matousek--Nikolov--Talwar's multiplicativity of `gamma_2`
([primary arXiv](https://arxiv.org/abs/1408.1376)).  The missing hypothesis is
a lossless identification of `M_n` with that surrogate.  Li--Nikolov's sharp
separation results ([primary arXiv](https://arxiv.org/abs/2303.08167)) and the
explicit orthonormal relaxation above show that no generic identification is
available.

## 8. Circularity audit

Any proof of the best lemma must avoid the following hidden reconstructions of
the desired result.

1. Do not assume that minimizing signings are random-like, self-averaging, or
   locally distributed as independent Rademachers.  That is precisely the
   missing bridge from quenched expectation to deterministic minimum.

2. Do not commute `min_a` with disorder expectation, logarithm, a random
   partition, or a block restriction.  None of these commutations is valid by
   minimax alone.

3. Do not assume that restrictions of a global near-minimizer are near
   minimizers at their own orders.  Such hereditary stability would itself be
   a strong scale-transfer theorem.

4. Keep one global orientation `sigma`.  Giving blocks independent
   orientations proves a statement for `|H_1|+|H_2|`, not for
   `|H_1+H_2|` with the actual cross term.

5. The cross edges number `nm`.  Deleting them, filling them arbitrarily, or
   paying their absolute contribution produces an error much larger than
   `o((n+m)^(3/2))` in balanced splits.

6. The lemma is only at fixed bounded `beta`.  One must take `n to infinity`
   before `beta to infinity`; choosing `beta=beta_n` without a uniform theorem
   would be circular.

7. A Gaussian/Bernoulli-process expectation over coefficients is not the
   minimum coefficient signing.  Majorizing measures may handle the `2^n`
   correlated tests for one random process, but do not supply the adversarial
   outer minimum.

8. Do not replace the two-sided statistic by frustration, induced-set
   discrepancy, or the independent bipartite norm.  Every such move needs an
   equality or a `1+o(1)` transfer; a universal constant factor is a fixed
   leading loss.

9. A finite witness family must approximate the **deterministic** maximum for
   every candidate signing relevant to the outer minimum.  Distributional
   Gaussian sparsification with offsets does not provide that statement.

10. Proving `L_advGT` by first assuming convergence (or a common limiting
    optimizer/profile) merely repackages the target.  The proof must derive
    the sublinear defect directly from finite-size structure.

## 9. Confidence and candidate card

| Item | Truth | Tractability | Conditional implication |
|---|---:|---:|---:|
| `L_advGT` | 45% | 12% | 99% |
| Star-profile statement (S) | 25% | 8% | 98% |
| Lossless factorization route | <5% | <5% | not viable as stated |

My confidence that the best route can currently be completed from the cited
literature without a new adversarial interpolation idea is about 10%.  The
conditional implication itself is essentially airtight: its only inputs are
the quantified lemma, approximate superadditivity, and the exact soft-max
squeeze.

| Field | Frozen entry |
|---|---|
| Domain | Matrix discrepancy of the degree-two Walsh system; antipodal cut-code covering |
| Imported theorem(s) | Guerra--Toninelli interpolation; Hammersley approximate subadditivity |
| Problem translation | `M_n=disc(W_n)=E_n-2 rho(C_n^+)`, with `|C_n^+|=2^n` |
| Proposed mechanism | Deterministic two-sided soft pressure with adversarial near-superadditivity |
| Exact missing lemma | `L_advGT` in Section 4 |
| Why strictly weaker | Only a lower-envelope pressure inequality; finite orders, optimizers, and nonexposed coset histograms are unrecoverable |
| Ledger collisions | Not inspected by design; proposal frozen before archive exposure |
| Falsification test | Linear superadditivity defect `Delta_B(n,m)` along a growing sequence |
| Specialist confidence | 45% truth / 12% tractability / 99% implication |
| Director judgment | Hold for an adversarial-interpolation attack; reject generic factorization |
