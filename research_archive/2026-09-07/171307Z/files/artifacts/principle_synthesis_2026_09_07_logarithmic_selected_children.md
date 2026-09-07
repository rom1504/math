# Actual selected optimal children in every bounded-cap signing

2026-09-07. Elementary positive selection theorem. This complements, rather
than contradicts, the logarithmic random-restriction obstruction. It is not
a comparable-scale hierarchy or a convergence recurrence.

## Theorem

Fix `C<infinity` and `0<delta<1/2`. For all sufficiently large `n`, EVERY
hollow complete signing `A` with `Q(A)<=C n^(3/2)` contains EVERY prescribed
labelled signing `T` of order

```math
k\le(1/2-\delta)\log_2 n
```

as an exact induced signing, with no edge edits and no switching necessary.
In particular, one can select an actual order-`k` exact minimizer inside
every order-`n` exact minimizer, and indeed inside every bounded-cap input.
The target may depend on `n`; the statement is uniform over all targets.
It does not claim that the selected target is computationally easy to find
without already knowing its signs.

## 1. A large spectrally regular core

The archived Grothendieck/diagonal-majorant core theorem gives a principal
set of size `m>=n/2` for which

```math
\|A_R\|_{\rm op}\le K_C\sqrt m
```

with `K_C` depending only on `C`. For example, deleting coordinates whose
majorant diagonal exceeds `8K_G C sqrt(n)` gives
`K_C=8 sqrt(2) K_G C`. No refill is used: this remains a genuine induced
subsigning. This fixed-fraction core theorem is already proved in
`fresh_range_and_spectral_regularization_2026_09_05.md` and reconstructed in
`decisive_independent_entropy_compatible_spectral_regularization_2026_09_07.md`.

## 2. Direct greedy embedding lemma

Let `B` be a hollow order-`m` signing with `||B||op<=K sqrt(m)`. Fix
`epsilon in (0,1)`, and put `a=(1-epsilon)/2`. For any candidate set `S`,

```math
\sum_i\left|\sum_{j\in S}B_{ij}\right|^2\le K^2m|S|.
```

Thus at most `K^2m/(epsilon^2 |S|)` vertices have absolute signed sum into
`S` greater than `epsilon |S|`. A vertex outside this bad set has at least
`a|S|-1` neighbors of EITHER prescribed sign in `S`, after deleting itself
if necessary. The minus one accounts for the hollow diagonal and exclusion
of the newly chosen vertex.

Maintain one candidate set for every unembedded target vertex. After `t`
choices, each set has size at least

```math
m a^t-\frac{1-a^t}{1-a}\ge m a^t-2.
```

If `m a^k>=4`, all relevant candidate sizes are at least `m a^k/2`. The
union of the bad sets for at most `k` future candidates has size at most
`2kK^2/(epsilon^2 a^k)`. Therefore the current candidate set contains a
vertex good for ALL future sets whenever

```math
m a^{2k}>4kK^2/\epsilon^2.                              (1)
```

Choose it and intersect every future candidate with the neighbor set of
the sign prescribed by `T`. Previously selected vertices stay excluded.
Induction embeds every edge of `T` exactly. This proof uses neither a
random-target assumption nor a sequential conditional-independence claim.

## 3. Substitution of the logarithmic scale

Take `epsilon=1-2^(-delta)`, so `a=2^(-(1+delta))`. For the stated `k`,

```math
m a^{2k}\ge\tfrac12 n^{\delta+2\delta^2}.
```

This eventually dominates the right side of (1), which is `O_C,delta(log n)`.
Also `m a^k -> infinity`. The greedy lemma on the unchanged regular core
proves the theorem.

## 4. Selection and random restriction are different quantifiers

The previously proved logarithmic local law shows that, for growing
`k<=(1-delta)log_2 n`, a UNIFORM induced child (after harmless switching)
has iid-like distribution in total variation. Its cap is typically at
least `(2/3)sqrt(2/pi) k^(3/2)`, strictly above the known all-order upper
for `M_k`. The present theorem simultaneously supplies exact optimal
children at the smaller logarithmic constant `1/2-delta`.

There is no contradiction: their selection probability can be extremely
small, and total-variation convergence does not erase events whose iid
probability also tends to zero. In particular a random-restriction no-go
must not be promoted to nonexistence of selected optimal children.

These children are only logarithmic in the parent order. They do not form
the comparable-size selected hierarchy needed for the current summable
cross-order target. This is a proved boundary between random-law and
selection statements, not an assertion that a small local state controls
the parent cap.
