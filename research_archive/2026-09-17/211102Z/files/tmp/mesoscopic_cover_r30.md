# Wave 30 root route: mesoscopic collision equivalence

Status: verified finite set-system theorem and asymptotic reduction.  It does
not prove the required collision probability.

Let `Omega=binom([n],m)`, let `mathcal A` be any finite family of eligible
row-good block cosets, and let `H_a subset Omega` be the hit set of candidate
`a`.  Write `M=|mathcal A|` and

```math
delta_* = min_{w in Delta(Omega)} max_{a in mathcal A} w(H_a),
```

the dual value in (10.873).  For an integer `r>=1`, define the common-coset
event

```math
J_r = {(S_1,...,S_r): there is a in mathcal A with S_j in H_a for all j}.
```

For every selector law `w`, the union bound and any maximizing candidate give

```math
max_a w(H_a)^r
 <= w^r(J_r)
 <= sum_a w(H_a)^r
 <= M max_a w(H_a)^r.
```

Taking roots and then the infimum over `w` proves the exact sandwich

```math
M^(-1/r) inf_w w^r(J_r)^(1/r)
 <= delta_*
 <= inf_w w^r(J_r)^(1/r).
```

Thus the fractional cover problem can be replaced, at the correct exponent,
by a mesoscopic common-coset collision problem.

For at most `k` blocks, the number of labeled vertex-to-block maps is at most
`k^n`, and each partition has at most `2^n` diagonal cosets.  Hence the number
of distinct eligible candidates obeys

```math
M <= (2k)^n.
```

Put `L=n^(3/4-c)` and choose

```math
r = ceil(n log(2k)/L).
```

Then `log(M)/r<=L`; for the entropy-matched
`k=Theta(L/log n)`, this is

```math
r=Theta(n^(1/4+c) log n),
```

the same scale as the balanced block size.  The following is therefore a
sufficient successor to (10.875): for every selector law `w`,

```math
w^r(J_r) >= exp(-K r L).
```

Indeed the sandwich yields `delta_*>=exp(-(K+1)L)` and hence
`tau_*<=exp((K+1)L)`.  Conversely, `delta_*>=exp(-K L)` implies
`w^r(J_r)>=exp(-K r L)` for every `w`.  Up to the harmless candidate-count
term `exp(-L)`, the two statements are exponent-equivalent.

This reduction makes the missing higher-order content precise.  It is not
necessary to glue all `binom(n,m)` selectors into one coset; it is enough that
an adversarial iid batch of only
`Theta(n^(1/4+c) log n)` selectors has common-coset probability
`exp(-O(rL))`.  On the other hand, pairwise or fixed-order agreement is not
enough, and Wave 29's exact A9 triangle shows that even `r=3` can fail for
exact-child witnesses at a tight finite cap.

In particular, a deterministic **mesoscopic batch lemma** would already
close the route: for every ordered `r`-tuple of selectors, choose one aligned
witness per selector so their balanced signature count is at most `k` and
some balanced refinement of the signature partition has its entire coset
under the row cap.  Then every tuple lies in `J_r`, so the preceding argument
gives `delta_*>=M^(-1/r)>=exp(-L)`.  This is much weaker than selecting one
coherent family for the whole slice, but it still requires higher-order and
full-refinement control absent from pairwise agreement.

The equally sharp falsifier is a sequence of selector laws for which

```math
w^r(J_r)=exp(-omega(rL)).
```

The upper half of the sandwich then gives `delta_*=exp(-omega(L))`.  No such
minimizer-specific law, and no proof of the sufficient batch collision, is
known.

The exact `A_9` threshold-zero `3+3+3` systems from Wave 29 provide a useful
finite check of the hierarchy.  At both `(m,C)=(5,80)` and `(6,80)`, every
unordered selector pair and every unordered selector triple (repetitions
allowed) has a common row-good hit coset, although no coset hits the whole
slice and the game values are respectively `839/995` and `21/22`.  Thus
fixed-order common cover can coexist with a nontrivial fractional game; the
mesoscopic order in the theorem cannot simply be replaced by pairwise or
triple agreement.  These facts are exactly rechecked by
`tmp/common_coset_minimax_r29.py`.
