# Audit: mergeable Boolean-port reservoir

**Object audited.** `boolean_port_mergeable_reservoir.md`, Theorem MR.1.

**Verdict.** PASS for a fixed finite labelled universe and a composition DAG
declared independently of the public randomness.  It is not a uniform theorem
over adaptively selected subsets.

## Algebra and sampling law

For one replica, a uniform total order on `U` has a unique minimum on every
nonempty subset.  Symmetry makes that minimum uniform on the subset.
Independence of the `k` total orders makes the `k` winners independent, so
the response payload is exactly a with-replacement empirical sample.  No
independence across different aggregate nodes is asserted or needed.

For each order, the reservoir component obeys

```math
\min(E\cup F)=\min(\min(E),\min(F)).
```

Coordinatewise minimum is associative, commutative and idempotent.  The rank
and row type stored in MR.1 suffice to perform it; the full constituent set is
not reconstructed.  The complete carrier also stores `|E|`.  Under disjoint
composition this count adds, so the unnormalized response
`L_E=p|E|R_E` remains recoverable.  The full count-plus-reservoir carrier is
associative and commutative but, correctly, is not idempotent as a multiset
state.

## Probability constants

RC.1 gives expected supremum error at most `4/sqrt(k)` for each fixed node.
Changing one of its `k` iid winners changes every normalized response by at
most `1/k`, hence changes their supremum error by at most `1/k`.  McDiarmid
therefore has denominator

```math
\sum_{ell=1}^k(1/k)^2=1/k
```

and upper tail `exp(-2kt^2)`.  A union bound over `T` prescribed nodes is
valid despite cross-node dependence.  Setting `t=sqrt(log(T/delta)/(2k))`
gives failure probability `delta`.  The two sufficient conditions

```math
k>=64/eta^2,
\qquad k>=2\log(T/\delta)/eta^2
```

make the expectation and deviation terms at most `eta/2` each, verifying
MR.8.

## Scope and information

The result removes **depth accumulation** for a declared union tree.  If only
the root matters, its guarantee is independent of both depth and node count.
Uniform accuracy at all declared nodes introduces `log T`, not a sum of local
errors.  An aggregate chosen after inspecting the sketch is not covered;
certifying all possible subsets by substituting `T=2^|U|` would destroy the
small-state conclusion.

Exact priorities do not require unbounded real precision in the finite model:
one may sample a random permutation and store its integer rank.  With public
randomness the stored merge metadata costs `k ceil(log_2|U|)` bits in addition
to the `k(p-1)` response-payload bits and one
`ceil(log_2(|U|+1))`-bit cardinality.  This rank cost, and not response error,
grows with the labelled universe.

The reservoir component is a congruence for set union; the full state with its
additive count is a congruence for disjoint unions of labelled row occurrences.
The theorem makes no claim for overlapping multiset union, weighted mixtures,
max-plus switching, or quadratic interaction.
