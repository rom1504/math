# Orientation-balanced near-minimizers exist at all orders

2026-09-06. Director derivation; **independently audited** in
`transfer_seed_orientation_repair_independent_audit_2026_09_06.md`. This is a
one-sided-orientation repair theorem, not convergence and not isotropy.
It answers a limitation of the newly proved ALL-near-minimizer counterexample:
failure for every near-minimizer need not preclude selecting a balanced one.

Write `P(A)=max H_A`, `R(A)=max(-H_A)`, and `Q(A)=max(P,R)` in the
original single-copy normalization. All matrices below are genuine hollow
sign matrices, including edges not used by the intended gadget.

## 1. A finite random-remainder bound

On `N` vertices let `E` be any prescribed set of `e` edges. Assign independent
unbiased signs there and zero elsewhere, producing the weighted remainder `V`.
For any fixed Boolean spin, its energy is a sum of `e` independent signs.
The elementary exponential-moment bound therefore gives

```math
\Pr\{Q(V)>t\}\le 2^{N+1}\exp[-t^2/(2e)].
```

For `e>0`, some deterministic choice thus satisfies

```math
 Q(V)\le L(N,e):=\sqrt{2e(N+2)\log 2}.
```

For `e=0` set `L=0`. This is an existence argument with a checkable full
Boolean cap, not a polynomial-time deterministic construction claim.
For every fixed internal matrix `W`, separately in both orientations,
`|P(W+V)-P(W)|<=Q(V)` and `|R(W+V)-R(W)|<=Q(V)`.

## 2. Repair a given bounded-cap seed at a nearby order

Reverse all signs if necessary so `P(A)=Q(A)>=R(A)`. Set
`Delta=P(A)-R(A)>=0`. Choose the smallest nonnegative integer `k`
such that `binom(k,2)>=Delta`, taking `k=0` when `Delta=0`.
When `Delta>0`, `k<=sqrt(2Delta)+2` and
`0<=binom(k,2)-Delta<k`.

Adjoin a negative clique on `k` new vertices. Before filling the bridge,
the block-diagonal weighted signing `W=A direct_sum (-J_k+I_k)` has

```math
 P(W)=P(A)+\lfloor k/2\rfloor,\qquad
 R(W)=R(A)+\binom k2.
```

Here the maximum negative-clique energy is `floor(k/2)`, including odd `k`;
the minimum is `-binom(k,2)`. Hence both oriented extrema of `W` lie between
`Q(A)` and `Q(A)+k`. Fill all `nk` bridge edges using Section 1. The actual
full signing `B` at `N=n+k` satisfies

```math
 Q(A)\le Q(B)\le Q(A)+k+L(N,nk),
 \quad |P(B)-R(B)|\le k+2L(N,nk).                    (1)
```

The lower cap bound follows from exact principal restriction, not the
two-sided perturbation estimate. For `Q(A)<=C n^(3/2)`,
`k=O_C(n^(3/4))`, `L=O_C(n^(11/8))`. Thus the cap loss is power-saving
and the change of order is sublinear. Every original near-liminf sequence
can be replaced by a nearby sequence with normalized orientation imbalance
`O_C(n^(-1/8))` and the same limiting normalized cap.

## 3. All-order near-minimal balanced realization

Assume only the already proved all-order bound `M_j<=C j^(3/2)` for large `j`.
Choose a constant `K>sqrt(2C)+3`. For each sufficiently large target order `N`,
put `r=ceil(K N^(3/4))`, `m=N-r`, and take an EXACT order-`m` minimizer `A`.
The integer `k` from Section 2 is at most `r` for large `N`.

On the `r` additional vertices prescribe a negative clique on `k` of them;
all other edges not in the old `m`-block are left to the random remainder.
There are at most `Nr` such random edges. Prior to filling, the remaining
`r-k` vertices are isolated weighted vertices and change neither oriented cap.
Section 1 yields a genuine target-order signing `B_N` with

```math
 M_N\le Q(B_N)\le M_m+k+L(N,Nr)
                    \le M_N+O_C(N^{11/8}),
 \qquad |P(B_N)-R(B_N)|=O_C(N^{11/8}).                (2)
```

The last comparison uses only monotonicity `M_m<=M_N`. The lower inequality
is the definition of `M_N`; no assertion that the old minimizer remains
minimizing after insertion is made. Consequently an orientation-balanced
near-minimizing signing exists at EVERY sufficiently large order, with
normalized cap error and orientation imbalance `O_C(N^(-1/8))`.

## 4. Exact scope and the remaining original-problem gap

This is not a cross-order estimate with a comparable split: the constructed
matrix uses an already available minimizer at order `N-O(N^(3/4))`.
It gives no upper construction from one fixed seed to unbounded target orders,
and no convergence of `M_N/N^(3/2)`. Its actual removed obligation is only
the EXISTENCE of a balanced near-minimizer sequence, if a future argument
requires that property alone at the displayed tolerance.

Balanced oriented extrema do not give a common isotropic spin law, radial
subgradient, low cap response complexity, or bounded normalized operator norm.
The planted clique can itself have norm of order `N^(3/4)`. In particular,
this result cannot be substituted into a theorem requiring an `O(sqrt(N))`
operator bound or `o(sqrt(N))` near-ground slack.

The finite one-edge minimality condition and the separate conjecture that
EXACT minimizers have orientation gap at most two remain distinct questions.
This proof neither assumes nor proves that conjecture.
