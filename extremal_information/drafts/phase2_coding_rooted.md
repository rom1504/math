# Phase 2 coding report: a syndrome-rooted feature algebra

**Status.**  The theorem and examples below are proved.  The syndrome
interpretation is classical coding theory; the operational minimality,
response-complexity statement, and outer-spectrum falsifier are the additions
relevant to this program.  Nothing here concerns the original signing
convergence problem.

## 1. Director verdict

There is a nontrivial code model in which the smallest reusable extremal state
is much smaller than the full rooted distance landscape and closes exactly
under repeated composition.

Fix a labeled syndrome interface

```math
G=\mathbb F_2^w.
```

For a full-row-rank parity-check fragment `H`, retain the coset-leader profile

```math
\lambda_H(s)=\min\{\operatorname{wt}(e):He=s\},
\qquad s\in G.                                      \tag{CR.1}
```

If two fragments are joined by concatenating their columns, their code is a
fibre code, not a Cartesian product: corrections in the two blocks may cancel
their syndromes.  Nevertheless,

```math
\lambda_{[H_1\ H_2]}(s)
=\min_{u\in G}\{\lambda_{H_1}(u)+\lambda_{H_2}(s+u)\}. \tag{CR.2}
```

For binary Hamming correction this min-plus algebra synchronizes further:
`lambda_H` is equivalent to the set of distinct nonzero column types of `H`,
and composition is set union.  Thus the exact feature algebra has
`Theta(2^w)` bits, independent of composition depth.  This is exponentially
smaller than a complete labeled root table, and is polynomial size whenever
`w=O(log n)`.

The compression is sharp for the declared experiment.  If a summary must
answer the covering radius after **every future parity-check fragment** is
appended, then special code environments expose every column-support bit.
Uniform response error below `1/2` requires `2^w-O(w)` bits in the worst case.
This is an operational lower bound inside a structured linear-code class, not
the universal arbitrary-kernel packing of the feature-growth draft.

The root-averaged outer spectrum is strictly too coarse for this operation.
An exact pair below has the same complete outer pressure at every temperature,
but the same appended fragment gives covering radii two and one.  What is
missing is not the full labeled root map; it is the alignment of the distance
layers with addition in the syndrome group.

## 2. Model and exact state

Let `H` be a `w times n` binary matrix of rank `w`, with nonzero columns
`h_1,...,h_n in G`.  Put

```math
C_H=\ker H\subseteq\mathbb F_2^n,
\qquad
S_H=\{h_i:1\le i\le n\}\subseteq G\setminus\{0\}.  \tag{CR.3}
```

Repeated columns are permitted.  The length and the fixed labeled syndrome
space are part of the apparatus.  Zero columns can be deleted for the
covering-radius experiment and are excluded only to avoid irrelevant
bookkeeping.

For another fragment `E` with `w` rows define the future-environment response

```math
\mathcal R_H(E)=\rho\!\left(\ker[H\ E]\right),       \tag{CR.4}
```

where `rho` is Hamming covering radius.  This operation is richer than a
Cartesian product:

```math
\ker[H\ E]
=\{(x,y):Hx+Ey=0\},                                  \tag{CR.5}
```

so a nonzero syndrome in the first block may be cancelled in the second.

### Theorem CR.1 (syndrome-rooted response algebra)

For binary full-rank fragments over a fixed labeled `G=F_2^w`:

1. For every root `x in F_2^n`,

   ```math
   d(x,C_H)=\lambda_H(Hx),
   \qquad
   \rho(C_H)=\max_{s\in G}\lambda_H(s).              \tag{CR.6}
   ```

2. The function `lambda_H` depends only on `S_H`, and determines it:

   ```math
   s\in S_H\quad\Longleftrightarrow\quad\lambda_H(s)=1
   \qquad(s\ne0).                                    \tag{CR.7}
   ```

3. Concatenation obeys the associative min-plus group convolution (CR.2).
   Equivalently,

   ```math
   S_{[H_1\ H_2]}=S_{H_1}\cup S_{H_2}.               \tag{CR.8}
   ```

   The numerical convolution is nonexpansive:

   ```math
   \|f\star g-f'\star g'\|_\infty
   \le \|f-f'\|_\infty+\|g-g'\|_\infty.             \tag{CR.9}
   ```

4. Assume `w>=2`.  Up to one-to-one recoding, `lambda_H` is the coarsest
   exact deterministic quotient for the complete response experiment over
   unrestricted appended fragments on the same labeled syndrome group,

   ```math
   \bigl(\mathcal R_H(E)\bigr)_{E}.                  \tag{CR.10}
   ```

   It therefore restores exact compositional sufficiency without retaining
   the labeled map `x -> d(x,C_H)` or the code itself.

#### Proof

For (CR.6), a correction `e` sends `x` into `C_H` exactly when
`H(x+e)=0`, equivalently `He=Hx`.  Minimizing its Hamming weight gives the
first identity.  Surjectivity of `H` makes every syndrome occur among the
roots, giving the second.

If a minimum-weight correction uses two coordinates carrying the same column
type, delete both coordinates.  Their syndromes cancel over `F_2` and the
weight drops by two.  Hence a shortest correction uses each column type at
most once.  The profile depends only on `S_H`; conversely, its level-one set
is exactly `S_H`.  This proves (CR.7).

Split a correction for `[H_1 H_2]` as `(e_1,e_2)` and put `u=H_1e_1`.
The second block must contribute `s+u`, and minimization within the two blocks
gives (CR.2).  The support interpretation gives (CR.8).  Associativity is
associativity of concatenation (or of group convolution).  Inequality (CR.9)
follows because every candidate sum changes by at most the displayed
right-hand side, and taking a minimum cannot increase a uniform error.

It remains to prove operational minimality using only code environments, not
arbitrary formal terminal potentials.  For each nonzero `s in G`, let `E_s`
have one column of every type in

```math
(G\setminus\{0\})\setminus\{s\}.                    \tag{CR.11}
```

Thus `E_s` has `2^w-2` columns; its size is part of the declared unrestricted
query class and does not depend on the unknown fragment `H`.

If `s in S_H`, the composite contains every nonzero column type, so its
covering radius is one.  If `s notin S_H`, its support omits only `s`.  Every
other nonzero syndrome has distance one, while `s` has distance two: choose
`u notin {0,s}` and write `s=u+(s+u)`.  Both summands are nonzero and different
from `s`; such a `u` exists because `w>=2`.  Therefore

```math
\mathcal R_H(E_s)=
\begin{cases}
1,&s\in S_H,\\
2,&s\notin S_H.
\end{cases}                                         \tag{CR.12}
```

The complete future response determines `S_H`, hence `lambda_H`.  Conversely,
(CR.2) and (CR.6) compute the response to every `E` from `lambda_H` and
`lambda_E`.  Thus the response function and `lambda_H` determine one another,
which is exact quotient minimality. `square`

### What has synchronized

The full root identity `x` and all correction multiplicities disappear.
Translation symmetry first quotients roots to their syndrome.  Binary
cancellation then synchronizes the coset-leader profile to the Boolean support
of available generator types.  This is a finite deterministic instance in
which a rooted extremal state becomes a strict function of a smaller global
algebraic parameter.

It is important that this uses the binary Hamming metric.  With coordinate
weights, nonlinear symbol costs, labeled punctures, or fields that inspect
particular roots, multiplicities and labels can return.

## 3. Sharp information growth

### Theorem CR.2 (worst-case exact response complexity is `Theta(2^w)` bits)

Fix `w>=2`.  Recording `S_H` takes `2^w-1` bits and answers every response in
(CR.10).  Conversely, there is a fixed fragment length for which any
deterministic summary answering every response with uniform additive error
strictly below `1/2` requires at least

```math
2^w-1-w-\log_2(2^w-w)                               \tag{CR.13}
```

bits.  Hence the exact worst-case response complexity is `Theta(2^w)` bits.

#### Proof

Fix a basis `B` of `G` and let

```math
N=|G\setminus(\{0\}\cup B)|=2^w-1-w,
\qquad t=\lfloor N/2\rfloor.                         \tag{CR.14}
```

For every `t`-subset `U` of the remaining nonzero vectors, form the fragment
with one column of every type in `B union U`.  All fragments have the same
length `w+t` and full rank.  Equations (CR.7) and (CR.12) show that their
future-response functions are distinct, with sup distance one.  If two shared
one summary under error below `1/2`, the triangle inequality would put their
response distance below one, a contradiction.  Therefore at least
`binom(N,t)` messages are required.  The largest binomial coefficient obeys

```math
\binom N{\lfloor N/2\rfloor}\ge {2^N\over N+1},      \tag{CR.15}
```

which gives (CR.13).  The support bit vector is the matching upper bound up to
`O(w)` bits. `square`

The packed fragments have the common length `w+t=Theta(2^w)`.  The theorem
uses unnormalized response error below half the integer lattice spacing; it
does not give the same rate for `poly(w)`-length fragments or for distortion
growing with `w`.

This lower bound is about message bits, unlike a tropical factorization rank,
which counts real-valued terms and by itself gives no bit bound.  It is also
scale-specific: a one-bit support change causes an integer response change of
one.  If all radii are divided by a growing normalization, an approximate
rate theorem must be reproved at that normalized distortion.

Under the uniform prior on the `N` optional support bits, the complete special-
environment response vector is an isometric copy of the latent Hamming cube.
Thus the independently audited posterior-width theorem, if promoted in its
present normalization, immediately gives a Shannon version.  With average
squared error `Delta` over the `N` special environments,

```math
I(S_H;Z)\ge N\,[1-g(\min\{4\Delta,1\})].             \tag{CR.16}
```

Equation (CR.16) is a corollary of that general theorem, not needed for
Theorem CR.2.

## 4. The state is strictly smaller than the code

Take `w=2`, identify the nonzero syndromes with `1,2,3`, and use the two
length-five column lists

```math
H_A=(1,1,1,2,3),
\qquad
H_B=(1,1,2,2,3).                                    \tag{CR.17}
```

Both supports contain all three nonzero syndromes, so

```math
\lambda_{H_A}=\lambda_{H_B}=(0,1,1,1).              \tag{CR.18}
```

They are response-equivalent after every future fragment.  Nevertheless,
their code weight enumerators are

```math
W_{C_{H_A}}(z)=1+3z^2+3z^3+z^5,
```

```math
W_{C_{H_B}}(z)=1+2z^2+4z^3+z^4.                    \tag{CR.19}
```

The first follows by choosing the first three bits freely: their parity is
repeated in the last two positions.  For the second, the first pair and
second pair have the common parity carried by the fifth bit.  Different
weight enumerators rule out coordinate isometry.  Thus `lambda_H` does not
secretly reconstruct the code.

The full labeled root-distance map would reconstruct the code as its zero
set.  The syndrome profile forgets both the root-to-syndrome map and codeword
geometry not seen by the declared future covering-radius experiment.

## 5. Exact falsifier for the outer spectrum

For a full-rank binary linear code,

```math
O_{C_H}(z)
=|C_H|\sum_{s\in G}z^{\lambda_H(s)},                \tag{CR.20}
```

because every syndrome fibre has size `|C_H|`.  Thus the outer spectrum keeps
only the histogram of the syndrome profile.

Let `w=3` and take the length-four fragments

```math
H_A=(1,2,3,4),
\qquad
H_B=(1,2,4,7).                                      \tag{CR.21}
```

Indexing syndromes by `0,...,7`, their profiles are

```math
\lambda_{H_A}=(0,1,1,1,1,2,2,2),
```

```math
\lambda_{H_B}=(0,1,1,2,1,2,2,1).                  \tag{CR.22}
```

Both have histogram `(1,4,3)`, and both codes have size two, so

```math
O_{C_{H_A}}(z)=O_{C_{H_B}}(z)=2+8z+6z^2.            \tag{CR.23}
```

Append the same three-column fragment

```math
E=(1,3,5,6).                                        \tag{CR.24}
```

This appended fragment is itself full row rank.

The first composite has column support `{1,2,3,4,5,6}`; syndrome seven needs
two columns, so its radius is two.  The second has every nonzero column type,
so its radius is one:

```math
\rho(\ker[H_A\ E])=2,
\qquad
\rho(\ker[H_B\ E])=1.                             \tag{CR.25}
```

Consequently the outer spectrum is sufficient for its declared pressure and
Cartesian-product experiment, but not for parity-check fibre composition.
The minimal repair here is the group-labeled profile `lambda_H`; no labeled
table on all `2^n` roots is required.

## 6. Relation to trellises and the exact tropical-rank theorem

The closest modern theorem is K. Sheshadri,
[*Trellis State Complexity as an Exact Tropical Factorization Rank*](https://arxiv.org/abs/2607.23471)
(arXiv:2607.23471v1, 26 July 2026).  It was already reconstructed in the main
ledger, Section 10.64.4, so it is a validation rather than a new project
retrieval.

For a fixed linear code `C subset F_2^m` and coordinate cut `L sqcup R`, put

```math
W(x_L,x_R)=d((x_L,x_R),C),
```

```math
s=\dim C-\dim C_L-\dim C_R.                         \tag{CR.26}
```

The paper proves exactly

```math
\operatorname{rank}_{\min,+}W
=\operatorname{rank}_{\rm trop}W
=\operatorname{rank}_{\rm Kap}W
=2^s.                                               \tag{CR.27}
```

The upper bound is the classical door identity through the `2^s` minimal-
trellis states.  The lower bound chooses one lifted codeword per state class,
obtaining a `2^s` square submatrix with zero diagonal and positive
off-diagonal entries.  A min-plus rank-one term cannot be tight on two
diagonal cells, and the identity is the unique zero-cost tropical
permutation.  This verifies the exact factor count without a complexity
assumption.

Three distinctions are essential.

1. Equation (CR.27) counts min-plus scalar terms.  It is **not** a generic
   bit, mutual-information, or metric-entropy lower bound; one real-valued
   factor can carry arbitrary description information.
2. It concerns the complete conditional root-distance table at one fixed
   coordinate cut.  It does not assert sufficiency for named punctures,
   arbitrary bridges, or every future operation on the code.
3. In maximum-energy convention, set `K=-W`.  The min-plus factorization
   becomes a max-plus factorization of `K`; covering radius still applies an
   outer maximum to the distance profile.  Mixing these two extrema without
   the sign change reverses the claim.

The present syndrome theorem adds a different statement.  It varies the
parity-check fragment, fixes the boundary group, proves an associative
repeated-composition algebra, identifies its exact operational quotient for
all appended-fragment radius queries, and gives a sharp bit lower bound over
that response class.  It neither implies nor improves (CR.27).  Conversely,
the door identity suggests (CR.2), but the rank theorem alone does not prove
support synchronization, future-environment minimality, or (CR.13).

For a fixed code or restricted fragment family the classical minimal trellis
may quotient the interface below `G`.  Theorem CR.2 is a worst-case statement
for the universal family over a fixed labeled syndrome apparatus, not a claim
that every individual code needs all `2^w` profile entries.

## 7. Reproducibility

The script
[`verify_phase2_code_syndrome_profiles.py`](../experiments/verify_phase2_code_syndrome_profiles.py)
checks:

- all `16` ordered compositions of the four spanning supports in `F_2^2`;
- all `8464` ordered compositions of the `92` spanning supports in `F_2^3`;
- all special-environment exposing identities in those dimensions;
- the two different code weight enumerators in (CR.19); and
- the outer-spectrum collision and radii in (CR.21)--(CR.25).

The deterministic output is
[`phase2_code_syndrome_profiles_results.json`](../experiments/phase2_code_syndrome_profiles_results.json).

The computation is verification of the finite examples, not evidence for the
general theorems, whose proofs are above.

## 8. Scope and falsifiers

The state should be rejected or enlarged when:

1. a future query names a coordinate or root rather than only the syndrome
   interface;
2. finite-temperature counts matter--then the syndrome weight enumerators,
   not only their lowest exponents, compose by ordinary convolution;
3. coordinate weights distinguish duplicate column types;
4. the two blocks use incompatible or unknown identifications of `G`;
5. the environment couples through a nonlinear observable not factored by
   syndrome addition; or
6. `w=Theta(n)` and exact `Theta(2^w)` state is no longer sub-landscape at the
   desired scale.

Theorem CR.2 is also an exact-scale statement.  It does not preclude a much
smaller state at additive error `epsilon w`, nor does the tropical rank theorem
preclude approximate min-plus factorization.  Those are distinct open
questions.

## 9. Director checkpoint

### Did this produce a theorem that was previously hard to formulate?

Yes, with a restrained novelty claim.  The code algebra itself is classical
syndrome decoding.  The response framework forced three statements that are
not merely vocabulary:

1. the coset-leader profile is the coarsest exact quotient for a declared
   family of **future code compositions**, not just a decoding lookup table;
2. its exact response complexity has the sharp growth law `Theta(2^w)` bits;
3. the root-averaged outer pressure is falsified by an explicit fibre
   composition, while the needed repair remains a strict quotient of the
   labeled landscape.

This is a genuine positive known-model validation: a non-mean-field landscape
class has a composable state strictly smaller than its full labeled
root-distance landscape and code.  It is not smaller than the complete
declared future-response vector, which Theorem CR.1 proves equivalent to the
profile.  It also explains when the same state ceases to be useful: growing
syndrome width forces exponential exact information.

### Strongest next theorem

The exact theory is now complete enough that the next question should be
approximate, not another exact state variant:

> Determine the response-metric entropy of binary syndrome profiles under
> uniform additive distortion `epsilon w` for all appended-fragment covering-
> radius queries.  Either construct a quotient of size
> `exp(o(2^w))` for fixed `epsilon>0`, or prove a packing of realizable
> profiles with pairwise future-response distance `Omega(w)`.

That theorem would decide whether algebraic closure remains informative when
the interface width grows, and would connect this coding model to the
framework's extremal rate--distortion question at the correct normalized
scale.
