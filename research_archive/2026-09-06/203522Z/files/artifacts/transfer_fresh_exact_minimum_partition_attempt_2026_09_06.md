# Exact-minimum partition selection: the exchange calculation and its gap

Date: 2026-09-06. Fresh convergence track. This is a bounded attempt at
the surviving ACTUAL-minimum recurrence, after the pointwise obstructions
were banked. The exchange identities below are proved; the desired
partition theorem is not.

## 1. The exact selection target

Write `b(A)=Q(A)^(2/3)` and `b_n=M_n^(2/3)`. At each prescribed pair
of comparable orders `m,n`, `1/2<=m/n<=2`, the director's sufficient
target is the existence of at least ONE exact minimizing parent `A`
of order `N=m+n` and an `m`-subset `T` with

```math
b(A_T)+b(A_{T^c})\le b_N+C\sqrt N.                       (P)
```

The parent may depend on `(m,n)`. Since actual child caps dominate their
order minima, (P) implies the reverse-Fekete inequality in
`transfer_fresh_reverse_fekete_2026_09_06.md`, and hence convergence.
No claim over all near-minimizers, all exact minimizers, or random
partitions is needed.

Choose `(A,T)` lexicographically: first minimize `Q(A)` over all
order-`N` signings, then minimize `b(A_T)+b(A_(T^c))` over all such
parents and all prescribed-size partitions. This minimum exists because
the search space is finite.

## 2. Exact two-vertex exchange

This identity was independently derived on the seed track and checked
here. Put `U=T^c`, choose `v in T,w in U`, and let
`C=A_(T\{v})`, `D=A_(U\{w})`. For a core `C` and signed row `a`, set

```math
E_C(a)=\max_x\{|H_C(x)|+|a\cdot x|\}.
```

This is the EXACT cap of the one-vertex extension, by optimizing the
new Boolean coordinate. Minimality of the partition gives

```math
E_C(A_{v,T\setminus\{v\}})^{2/3}
+E_D(A_{w,U\setminus\{w\}})^{2/3}
\le E_C(A_{w,T\setminus\{v\}})^{2/3}
+E_D(A_{v,U\setminus\{w\}})^{2/3}.                       (1)
```

The edge `vw` appears in none of these four extensions. The sizes are
unchanged. This inequality needs only an optimal partition of a fixed
parent, not global minimizing status of the parent. Averaging donor
rows does not currently help: `E_C` is convex in its row, and Jensen
does not give the upper exchange estimate required for (P).

## 3. A genuinely exact-parent consequence: exact ground blockers

The following improvement DOES use the lexicographic exact-parent
choice. Suppose an edge `e={i,j}` internal to `T` can be flipped so
that `Q((A_T)^e)<Q(A_T)`. Then

```math
Q(A^e)=M_N+2,                                            (2)
```

and there is an oriented EXACT parent ground `(sigma,z)` satisfying

```math
\sigma H_A(z)=M_N,\qquad
\sigma A_{ij}z_i z_j=-1.                                (3)
```

Meanwhile EVERY oriented exact child ground `(tau,x)` satisfies

```math
\tau A_{ij}x_i x_j=+1.                                  (4)
```

Proof. A single edge flip changes every Boolean energy, and therefore
the cap, by at most two. All possible cap values of a fixed order have
the parity of its number of edges. Global minimality gives
`Q(A^e)>=M_N`. Equality would preserve the first lexicographic objective
and strictly lower the second, which is impossible. Hence (2).

An oriented witness for (2) obeys
`sigma H_A(z)-2 sigma A_ij z_i z_j=M_N+2`. Since
`sigma H_A(z)<=M_N`, both asserted equalities in (3) are forced.
For an exact child ground, a negative value in (4) would instead
increase its oriented child energy by two, contradicting the assumed
strict child-cap improvement. This proves (4).

For a batch `F` of `k` INTERNAL edge flips that strictly improves the
secondary partition score, the same lex argument gives a parent witness
with

```math
\sigma H_{A^F}(z)\ge M_N+2,
\qquad M_N-\sigma H_A(z)\le2k-2.                         (5)
```

This is sharper by two than the ordinary global-minimum replacement
witness, which only has to show that the changed signing still has cap
at least `M_N`. Equation (5) permits changes in both internal blocks,
provided the secondary score improves and bridge edges are unaltered.

## 4. Why this has not proved selection

The blocker in (3) can depend on the edge. Equations (3)--(4) do not
supply one probability law on parent grounds that simultaneously
separates all child-improving directions. A child with cap above its
standalone order minimum may also be a strict local minimum under
small edge edits, so there need not be an improving edge at all.

This latter caveat is substantive, not an artificial finite exception.
The already archived Paley radial comparison plus parity in
`resumed_convergence_paley_radial_hamming_law_2026_09_06.md` gives
coefficient-Hamming local optimality through order `sqrt N` changes,
while the saturated Paley family remains extensively above the known
strict-subhalf minimum benchmark. The bordered version is also used in
`cross_order_exposed_shell_sparse_repair_no_go.md`. A fresh derivation
of the odd-flip special case was screened against that archive and was
not banked as a new theorem.

Replacing a high-cap child outright by a standalone minimizer does
produce an improving secondary score, but typically changes order
`m^2` coefficients. Then (5) gives a slack allowance of order `N^2`,
larger than the entire relevant `N^(3/2)` cap scale. The parent states
active after replacement need not be active before it. Merely saying
that the parent is globally optimal does not close this gap.

Equivalently, an internal block of an exact minimizing parent minimizes
a CONSTRAINED COMPLETION objective with the exterior held fixed; it is
not thereby a minimizer of its standalone Boolean cap. No monotonicity
of that completion objective in the scalar standalone cap has been
proved.

## 5. Finite information and current proof obligation

The separate screen
`transfer_seed_existential_partition_screen_2026_09_06.md` shows that
the current archived representatives pass (P) with `C=1` at every
prescribed comparable split tested. The order-12 representative has
no four/eight split with BOTH exact children, despite passing (P).
Thus exact heredity of every chosen parent is a stronger assertion than
needed even in these finite data.

The newly replayed external order-15 and order-16 matrices have favorable
balanced partitions too; their global minimum status was NOT imported
as a proved fact for this test. Those computations were already done
on the seed track and were not duplicated here.

The precise outstanding ingredient is therefore one of the following,
not an established conclusion: a simultaneous use of the exact blockers
in (3)--(5); a profile-sensitive cap-preserving large-block replacement;
or another mechanism producing one partition satisfying (P). The
minimum-value reverse-Fekete inequality remains unproved.
