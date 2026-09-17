# Independent audit of same-window charging memo

## Verdict

The algebraic identities (R12.1)--(R12.12), the factor `4` in the
residual harvest, the finite max-flow/min-cut formula, and the displayed
order-five centered demand are correct.  The result remains a conditional
finite reduction and not a descent theorem.  Several scope qualifications
are needed before ledger integration.

## Reconstruction

On each temporal path,

```math
\sum d=q_n-Q(A_L)=q_n-q_m-\varepsilon_L,
\qquad
\sum\lambda=q_n-(m/n)^{3/2}q_n.
```

Subtracting and averaging proves (R12.1).  Equation (10.329) is
`2a-c=E(d+g)`, which proves (R12.2); positive/negative parts then prove
(R12.3).  Parent coarsening is indeed the least positive mass by Jensen.

For `r_b=[c_b-partial_b]_+`, path-cover domination gives

```math
\sum_b\omega_b(a_b/c_b)r_b
\le\int\sum_{b\in\pi}\omega_b r_b\,d\nu(\pi)
\le4\sup_\pi\sum_{b\in\pi}\omega_b r_b.
```

The verified local inequality `c_b<=partial_b+Q(D_b)` gives
`r_b<=Q(D_b)`.  A chain uses one orientation per edge, its peeled siblings
are disjoint, and
`Q(D_b)<=|D_b|(|D_b|-1)<=(s-1)|D_b|`; hence (R12.5)--(R12.6) have exactly
the stated factor `4`, with no missing orientation factor.

The LP (R12.7) is a standard capacitated bipartite flow with one additional
universal source of capacity `C`.  Its minimum unserved mass is exactly

```math
\max_{J\subseteq I}[z(J)-C-\kappa\widehat a(N(J))]_+.
```

Summing the primal constraints and using (R12.6) proves (R12.11) with
`4 kappa(s-1)n`.  If `kappa` is uniformly bounded and
`s(n)<=n^{1/2-delta}` uniformly over all allowed comparisons, this residual
part contributes `O(n^{-delta})` after normalization and is summable down a
canonical dyadic chain.  The Hall deficiency must separately obey the full
tail condition (10.500).

For the displayed `A_5`, exact enumeration confirms `P=N=Q=8`.  With
`p=(1,1,1,1,1)` and `n=(1,1,1,1,-1)`, the split is `{4}|{0,1,2,3}`;
the capacities are `(8,8)` toward the singleton and `(0,0)` toward the
four-block, while the corresponding decrements are `8` and `0`.  Every
residual is zero.  Deleting vertex four gives

```math
d=0,\quad g=8,\quad a=4,\quad c=0,\quad\varepsilon_L=0,
```

and therefore (R12.12).  So the full Hall cut fails for this prescribed
endpoint pair/tree and its Section 10.62 allocation.

## Required scope corrections

1. The symbol `r_v` is used first for current order and then for negative
   credit.  Rename one of them.

2. The Hall theorem is exact for the stated **relaxed** LP, but its common
   credit sink permits every negative temporal increment and terminal excess
   to pay every demand, regardless of time or geometry.  Thus (R12.13) is the
   weakest condition only after unrestricted global cancellation by `C`.
   A genuinely local/causal transport theorem must make negative credits
   separate atoms with their own compatibility edges.  Parent-level demand
   atoms likewise pre-cancel incompatible deletion outcomes; transition
   atoms are required if outcome geometry matters.

3. The full-set cut `J=I` is sufficient for the scalar comparison, but it is
   essentially that scalar comparison restated after subtracting the small
   residual capacity.  It does not construct a transport.  All subset cuts
   construct endpoint transport only within the relaxed universal-credit
   model above.

4. The order-five obstruction is **not** optimized over endpoint ties.  Of
   the 25 positive/negative endpoint pairs, exactly five singleton/four-block
   pairs have zero residual everywhere.  The other twenty `2+3` pairs have
   positive residual (total residual capacity four at the root).  Therefore
   the example falsifies an error-free theorem for every prescribed legal
   tree, and the particular same-cut pointwise proposal, but it does not show
   that every endpoint choice/allocation for `A_5` has zero total residual.

5. `s=o(sqrt(n))` and `E_kappa=o(n^{3/2})` suffice for one comparison only.
   Bridging an arbitrary liminf requires uniformity over every allowed
   `(n,m)` and (10.500), as the memo later correctly notes.  The dyadic
   example applies only to the residual term unless the same tail estimate is
   imposed on `E_kappa`.

6. `check_residual_small_r12.py` produces the claimed integer data, but it
   has no assertions and does not check the temporal values, allocation,
   Hall LP, or tail statement.  Describe it as an enumeration/printout unless
   assertions for those claims are added.  Its `A_8` enumeration is not used
   in the stated order-five proof.

