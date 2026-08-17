# Independent audit: mesoscopic pair-query visibility

**Verdict: PASS, with a strict rooted-versus-unrooted frontier
qualification.**

All inequalities and constants in MQ.1--MQ.3 are correct under the
repository's edge-sum normalization.  The exact-sign completion in MQ.3 is
also correct for both augmented orientations.  The result is not subsumed by
the closest archived compiler theorems: its linear public shore uses joint
column cancellation to select one mesoscopically separated **root** at
leading scale while keeping the complete parent cap at the target scale.

The scope limitation in Sections 5--6 is mathematically essential.  The
theorem gives a packing of rooted slices `(A,z)`, not a packing of unrooted
exact-sign children.  It therefore does not yet amplify the FR.5 packing into
separated scalar parent responses, does not supply a reusable contextual
state, and has no cross-order or convergence implication.

No repair to the theorem statement or proof is required.  Two cautions
should remain visible in any canonical summary:

1. FR.5 proves a lower separation of order `M_n`, not an upper
   `O(M_n)` diameter.  MQ.2--MQ.3 apply to those pairs which also satisfy the
   upper bound in MQ.27; MQ.1 has no such upper-distance restriction.
2. The finite verifier checks the exact selector/completion algebra at
   `d=0`.  It does not computationally certify MQ.1, the probabilistic
   arbitrary-order existence theorem, or the asymptotic projective
   conversion.  Those parts are proved analytically.

## 1. Frozen sources and rerun

This audit concerns the following frozen files:

```text
extremal_information/drafts/mesoscopic_pair_query_visibility.md
sha256 3db483e770bfcb5df9e5adde53d047c8c230edb4e05f4142a6ef66c0dab195cc

extremal_information/experiments/verify_mesoscopic_pair_query_selector.py
sha256 f5917b56f5656fae8e0f54f9ca2800bee7560513299eebcafdc0349710448844

extremal_information/experiments/mesoscopic_pair_query_selector_results.json
sha256 49a989fd4bfaeb44ea402066dd1ef2d2a87d46b1c5b10ca51bfa70bb44e8c9e4
```

The verifier was rerun without using `/tmp`:

```text
.venv/bin/python \
  extremal_information/experiments/verify_mesoscopic_pair_query_selector.py \
  --output /home/math/quadra/tmp/mesoscopic_pair_query_selector_audit.json
```

The rerun output has SHA-256
`49a989fd4bfaeb44ea402066dd1ef2d2a87d46b1c5b10ca51bfa70bb44e8c9e4`,
exactly matching the frozen result.  All twelve instances at orders
`4<=n<=7` pass.  The script also compiles under the project virtual
environment, and `git diff --check` reports no whitespace error in the two
task sources.

## 2. Normalization and augmented cuts: PASS

The draft uses

```math
H_A(x)=\sum_{i<j}A_{ij}x_ix_j=\langle A,c(x)\rangle.
```

There is no hidden factor two.  For an augmented cut
`w=tau c(x)`,

```math
\max_w\langle b,w\rangle
=\max_x|H_b(x)|=Q(b).
```

There are at most `2^n` augmented cut words: `c(x)=c(-x)` gives at most
`2^(n-1)` ordinary cut words, and the extra orientation doubles this.  Thus
the union bound in MQ.1 is over the correct query class.

For two augmented words `z=sigma c(u)` and `z'=sigma'c(v)`, let

```math
k=\min\{d_H(u,v),n-d_H(u,v)\},\qquad q=k(n-k).
```

Their raw edge Hamming distance is either `q` or `E-q`, depending on the
relative augmented orientation.  Hence in both cases

```math
d_P(z,z')=\min\{q,E-q\}.                         \tag{A.1}
```

Since `q<=floor(n^2/4)=E/2+O(n)`, the alternative
`d_P=E-q=o(E)` is impossible for large `n`.  Therefore in the mesoscopic
regime of MQ.27 one has exactly

```math
d_P(z,z')=k(n-k).                                \tag{A.2}
```

If `aM_n<=d_P<=bM_n` and `M_n` is between two fixed positive multiples of
`n^(3/2)`, then

```math
{d_P\over n}\le k\le {2d_P\over n},
```

so `k` is between fixed positive multiples of `sqrt n`.  This proves the
conversion to MQ.9, including the case `sigma'=-sigma`.

## 3. MQ.1: PASS

Let `D={e:z_e!=z'_e}`, `|D|=h`, and `mu=t/h`.  The proposed coordinate law
exists for every `0<=mu<=1`.  For an augmented cut `w`,

```math
|\mathbb E\langle b,w\rangle|
=\mu\left|\sum_{e\in D}z_ew_e\right|\le\mu h=t.
```

For a sum of `E` independent variables in intervals of length two,
Hoeffding gives

```math
\Pr\{|S-\mathbb ES|\ge s\}\le2e^{-s^2/(2E)}.
```

With

```math
s=\sqrt{2E(n+3)\log2},
```

the union failure probability over at most `2^n` augmented cuts is

```math
2^n\,2e^{-(n+3)\log2}=1/4.                       \tag{A.3}
```

On the disagreement set,

```math
{\langle b,z\rangle-\langle b,z'\rangle\over2}
=\sum_{e\in D}b_ez_e
```

has mean `t`.  The one-sided Hoeffding bound at
`sqrt(2h log4)` is `1/4`.  The two good events therefore intersect with
positive probability.  This proves MQ.4--MQ.5 with the displayed constants.

Taking `t=delta M/2` is allowed because `h>=delta M`; it gives signal
`delta M-O(sqrt h)=delta M-O(n)` and cap `O(n^(3/2))`.  Projective
separation is enough for this application because raw Hamming distance is
always at least projective Hamming distance.  The query is orientation
sensitive, exactly as a labelled augmented-cut evaluation should be.

The same-order non-overlay observation MQ.6a is also correct.  Uniqueness of
Boolean Fourier coefficients forces `r=0` and

```math
C_e=A_e+lambda b_e.
```

If both agreement classes occur, exact signs require simultaneously
`|1+lambda|=1` and `|1-lambda|=1`, whose only common real solution is
`lambda=0`.  If only one class occurs, `b=+-A`, so the putative query is
only a scalar copy of the child.

## 4. MQ.8 and the public-block scales: PASS

For `c,k>=0`,

```math
\min\{k,c\}\ge {kc\over k+c}.
```

Titu Andreescu/Cauchy--Schwarz gives

```math
\sum_j{c_j\over k+c_j}
=\sum_j{c_j^2\over kc_j+c_j^2}
\ge {L^2\over kL+V},
```

which is MQ.8.

A rectangular Rademacher matrix `C` of size `m times n` exists with
`||C||<=L_0(sqrt m+sqrt n)` for an absolute `L_0`.  Row and column gauges
preserve this norm and its bipartite cap.  If `x_0` maximizes
`||C^Tx||_1`, then the sharp `p=1` Khintchine lower bound gives

```math
L\ge\mathbb E_x\|C^Tx\|_1
 \ge n\sqrt{m/2}.
```

Under `k<=Ksqrt n`, one has `m=n-k>=n/2` eventually, so
`L>=n^(3/2)/2`.  Conversely,

```math
L\le\sqrt n\,\|C\|\sqrt m\le2L_0n^{3/2},
\qquad
V=\|C^Tx_0\|_2^2\le\|C\|^2m\le4L_0^2n^2.
```

All constants in MQ.19--MQ.20 are therefore valid.

## 5. MQ.2: PASS

After replacing `v` by `-v` if needed, its disagreement set with `u` has
size `k`.  This replacement changes neither a quadratic child energy nor
the free-shore response `F_B(v)`.

The row and column gauges in MQ.15 are feasible.  Explicitly, if
`C^(0)` is the original public block and
`(C^(0))^Tx_0=(eta_jc_j)_j`, multiplying outside row `i` by
`x_(0,i)u_i` and column `j` by `eta_jy_(0,j)` gives

```math
C^Tu_{S^c}=(c_jy_(0,j))_j.
```

For arbitrary `x`, the synchronized rows contribute
`a y_0`, where `a=sum_(i in S)u_ix_i`, and the public rows contribute
`C^Tx_(S^c)`.  Hence

```math
F_B(x)\le n|a|+\|C^Tx_{S^c}\|_1\le nk+L.
```

At `u`, every column has sign `y_(0,j)` and magnitude `k+c_j`, so equality
holds and `y_0` is a maximizing shore state.  At `v`, the corresponding
magnitude is `|c_j-k|`.  Therefore

```math
F_B(u)-F_B(v)
=2\sum_j\min\{k,c_j\}
\ge {2kL^2\over kL+V}.                            \tag{A.4}
```

Using `k>=kappa sqrt n`, `k<=Ksqrt n`, and the bounds above, the numerator
is at least `(kappa/2)n^(7/2)` and the denominator is at most
`(2KL_0+4L_0^2)n^2`.  Thus

```math
F_B(u)-F_B(v)
\ge {kappa\over4KL_0+8L_0^2}n^{3/2},
```

exactly MQ.13.  Also

```math
\|B\|_{infinity to1}=F_B(u)=nk+L
\le(K+2L_0)n^{3/2}.
```

This establishes MQ.10--MQ.11.  The construction is existential because it
chooses a bipartite Boolean maximizer `x_0`, but this optimizes only the
public block `C`, not the child `A`.  It therefore does not reconstruct the
child landscape or use the unknown optimum at the target order.

## 6. MQ.3 and the orientation `sigma`: PASS

From `z=sigma c(u)` and MQ.21,

```math
sigma H_A(u)=\langle A,z\rangle\ge Q(A)-d.
```

The column gauge gives `u^TBy_0=F_B(u)`, while
`H_D(y_0)=Q(D)`.  In the exact-sign parent MQ.22, the three target
contributions consequently have the common sign `sigma`:

```math
H_A(u)+sigma u^TBy_0+sigma H_D(y_0)
=sigma\{\langle A,z\rangle+F_B(u)+Q(D)\}.        \tag{A.5}
```

For every `(x,y)`, triangle inequality and MQ.10 give

```math
|H_A(x)+sigma x^TBy+sigma H_D(y)|
\le Q(A)+F_B(u)+Q(D).                             \tag{A.6}
```

Equations (A.5)--(A.6) prove MQ.23.  When `d=0`, equality holds, so the
target is a global absolute maximizer and MQ.24 follows.  This argument is
unchanged for `sigma=-1`; then all three target contributions are negative.

For a rival old spin `v`, its rooted response is at most
`Q(A)+F_B(v)+Q(D)`.  The target rooted response is at least
`Q(A)+F_B(u)+Q(D)-d`, which proves MQ.25.  If MQ.2 replaced the original
`v` by `-v`, the rooted response is still unchanged: substitute
`y\mapsto-y` and use the evenness of `H_D`.

The displayed block matrix is a genuine hollow exact signing of order
`2n`: its two diagonal blocks are hollow exact signings and every cross edge
is the sign `sigma B_(ij)`.  With the edge-sum convention, each cross edge
appears once as `x^TBy`; no factor two is missing.  MQ.10 and the assumed
caps of `A,D` give total cap `O(n^(3/2))=O((2n)^(3/2))`.

Finally, a suitable `D` exists at every order.  A random hollow signing has
cap `O(n^(3/2))` by a direct Hoeffding union bound; global edge negation and
vertex switching can make any selected absolute ground positive and move it
to the prescribed `y_0` without changing the cap.

## 7. Verifier scope: PASS

The script constructs independent random `A,D,C`, performs the exact row
and column gauges, exhausts every old spin for `F_B`, and then exhausts all
`2^(2n)` parent states.  It checks:

- `F_B(u)=nk+L=max_xF_B(x)`;
- `F_B(v)=sum_j|c_j-k|`;
- the MQ.8 lower bound for the pair gap;
- both positive and negative values of the augmented orientation `sigma`;
- exact target-parent equality at `d=0`; and
- the rooted rival exclusion inequality.

The seeded corpus includes both signs of `sigma` (as determined by the
random ground state).  The frozen output's twelve records all pass.

This is an appropriate finite algebra verifier, not a certificate for the
asymptotic theorem.  It does not implement the probabilistic selection in
MQ.1, test a nonzero deficit `d`, or prove the uniform rectangular random
matrix norm bound.  Those omissions do not affect the proof, but the
canonical description should continue to call the experiment an algebra
and sign check.

## 8. Archive collision audit

The mechanism is close to, but not contained in, four archived results.

1. **Theorem 36.3 (pinned-response entropy).**  It compares many distinct
   sparse-flipped children under amplitude-`n` labelled fields and incurs a
   quadratic physical calibration.  MQ.1 compares two labelled roots; MQ.3
   gives a target-scale physical completion but still has only one child.
2. **Theorem 36.11/AO.2.**  Its `Theta(sqrt n)` rank-one shore compiles a
   fixed *fractional* projective gap.  At edge distance `Theta(M_n)`, its
   scalar response is only `Theta(n)`.  MQ.2 instead uses `Theta(n)` public
   columns with a spectrally flat diffuse block, synchronizes only the
   `Theta(sqrt n)` differing rows, and retains cancellation inside each
   column before the absolute values are summed.  This is a different
   leading-scale rooted selector.
3. **Theorems 21.29 and 21.31 (coordinate compiler and universal-pin
   barrier).**  The exact coordinate compiler has a quadratic baseline.
   The barrier concerns one future which robustly pins a target across every
   possible child.  MQ.3 is pair/root conditioned and assumes that `u` is
   already near-ground for a bounded-cap child, so it neither invokes nor
   contradicts universal locking.
4. **Theorem 21.66 (rowwise microcanonical compilation).**  That theorem
   realizes a prescribed field, but scalar cap control alone does not
   preserve a response gap; its successful applications require a separate
   affine endpoint response law.  MQ.2 directly constructs one joint
   free-shore response with the needed root gap.  It is not a separately
   paid scalar-channel decomposition.

The general public-continuation nonexpansiveness principle does not subsume
MQ.3 either: `(A,z)` and `(A,z')` have the same unrooted child landscape, so
that principle predicts exactly the remaining failure to turn this rooted
selector into a child response packing.

The word "bank" must be read literally as an externally indexed collection
of pair-specific contexts.  The MQ.1 signing `b` depends on `(z,z')`; the
MQ.2 block `B` depends on `(u,v)` and on the selected shore ground `y_0`.
When these roots are obtained from the shell of `A`, the bank is therefore
root-dependent and indirectly signing-dependent.  It is not one universal
compiler derived from `A` alone, and it cannot be applied to a child whose
root label has been forgotten.  This dependence is allowed by the theorem
and is exactly why its output is a rooted-slice response table rather than a
scalar contextual packing.

## 9. Precise frontier classification

The correct classification is:

```text
PROVED, universally for each eligible pair:
  * an exact-sign quadratic labelled query detects Theta(M_n) raw/projective
    separation with O(n^(3/2)) cap;
  * if the corresponding vertex roots are Theta(sqrt n) apart, a linear
    diffuse shore gives a complete exact-sign order-2n parent of
    O(n^(3/2)) cap which selects one root over the other by
    Theta(n^(3/2));
  * the same completion aligns either augmented energy orientation.

NOT PROVED:
  * that every pair in the FR.5 packing has the upper mesoscopic distance
    needed by MQ.2;
  * that root labels can be encoded into distinct exact-sign children at
    leading uniform distance;
  * any pairwise-separated scalar parent-response family;
  * a reusable query congruence, cross-order transfer, or recurrence.
```

In the campaign's scale this is a **Level-5 rooted-slice theorem**, not a
Level-6 transfer theorem.  It materially removes the physical-realization
objection for a *declared mesoscopic root pair*: a linear shore need not pay
quadratic cap, and joint same-switch cancellation can preserve a leading
selector margin.  It does **not** close the selected energy-to-ambient
amplification arrow because the child-state encoding is absent.

The next non-equivalent lemma is therefore exactly the one identified by the
draft: convert a rooted optimizer-exclusion table into separated scalar
responses of genuinely distinct exact-sign child systems while retaining
less information than one root-dependent near-minimizer child per word.
More pair-conditioned shores alone cannot supply that missing state map.
