# Exact minimizer optimized-bridge response audit

**Status.** Frozen exact finite experiment.  Query orders one and two are
proved by exhaustive bridge enumeration.  Query order three is a
solver-certified exact computation with independently re-evaluated witnesses
but no standalone CP-SAT proof object.  This note gives finite
hypothesis-ranking evidence, not an asymptotic theorem.

The frozen program and output hashes are

```text
program 620eb3941aad86c70079a52fc77d3af3bf037aabd4a05cb04c593199219cddc1
result  24047fb3563b51c983eea3e66907393319dae463fdbaaa63f570452a607e84a5
```

## 1. Correct gauge quotient

For a hollow signing `A` of order `n` and a hollow sign query `C` of order
`k`, define

```math
F_C(A)=\min_{B\in\{\pm1\}^{n\times k}}
Q\!\begin{pmatrix}A&B\\B^\mathsf T&C\end{pmatrix}.                 \tag{OB.1}
```

If `S,T` are diagonal sign matrices, then

```math
\begin{pmatrix}S&0\\0&T\end{pmatrix}
\begin{pmatrix}A&B\\B^\mathsf T&C\end{pmatrix}
\begin{pmatrix}S&0\\0&T\end{pmatrix}
=
\begin{pmatrix}SAS&SBT\\TB^\mathsf TS&TCT\end{pmatrix}.          \tag{OB.2}
```

Because `B -> SBT` bijects the complete bridge fibre, `F_C(A)` depends only
on the child switching classes.  The program verifies (OB.2) as an exact
integer matrix identity for every projective pair of child switches and every
saved optimal bridge.  It also checks permutation and global-sign transport.

There is one switching/permutation query class at `k=1,2`.  At `k=3` there
are two, classified by the triangle product
`tau(C)=c_12 c_13 c_23`.  The source minimizer classification also identifies
`A` with `-A`, while

```math
F_C(-A)=F_{-C}(A).                                                   \tag{OB.3}
```

Consequently the correct class invariant sorts the two order-three
coordinates.  It is not legitimate either to count byte-distinct matrices as
classes or to compare the two oriented coordinates without this swap.

## 2. Exact protocol

The authoritative exhaustive class representatives are the committed
`m3_minimizer_orbits.json` through `m8_minimizer_orbits.json` files.  Their
hashes and the hashes of the exact `M_3,...,M_11` sources are embedded in the
result.

- At `k<=2`, all `2^(nk)` bridges are evaluated in vectorized integer chunks
  against every projective parent spin.  This is complete enumeration.
- At `k=3`, all parent-spin inequalities are put into a deterministic,
  one-worker CP-SAT model.  Feasibility starts at certified `M_(n+3)` and
  advances by the forced parity step two.  Every infeasible target and all
  solver statistics are saved; every feasible bridge is independently
  evaluated over all parent spins.
- The CP model removes only proved symmetry: signed-permutation
  automorphisms of `C` and `B -> -B`.  It does not identify arbitrary bridge
  matrices.

The exact reproduction command and evidentiary qualifications are in
`experiments/exact_minimizer_optimized_bridge_response_protocol.md`.

## 3. Frozen response table

Write the canonical signature as

```text
(F_k=1, F_k=2, sorted(F_tau=+1,F_tau=-1)).
```

The complete finite result is:

| `n` | signature | exhaustive exact-minimizer classes |
|---:|---:|---:|
| 3 | `[4,4,5,7]` | `[0]` |
| 4 | `[4,5,9,9]` | `[0]` |
| 5 | `[5,9,10,10]` | `[0]` |
| 6 | `[9,10,12,12]` | `[0]` |
| 7 | `[10,12,13,15]` | `[0,2]` |
| 7 | `[12,12,13,15]` | `[1]` |
| 8 | `[12,15,17,17]` | `[0]` |
| 8 | `[12,15,19,19]` | `[1]` |

Thus the language is nonconstant on exact minimizer classes.  At order seven,
the one-vertex response already separates class 1 from classes 0 and 2; none
of the tested queries separates classes 0 and 2.  At order eight, both
classes have identical optimized responses to every query of order at most
two, but the order-three query separates them by two.

## 4. An exact order-eight witness

Take the positive-triangle query

```math
C=\begin{pmatrix}0&1&1\\1&0&1\\1&1&0\end{pmatrix}.                \tag{OB.4}
```

Use class 0 and class 1 from the authoritative order-eight file, with orbit
identifiers

```text
class 0  a5bbc9a3785f85e929367e670d5e7e0bf6bc46cec302dc054c01b5eee2d07fc9
class 1  75397bf3565083fefc1f5be5b402c0f0bb92c607871ebad3d2a01d53f39c7c5d
```

For class 0, the bridge

```math
B_0=\begin{pmatrix}
 1& 1&-1\\
-1& 1&-1\\
-1&-1& 1\\
 1&-1&-1\\
 1&-1& 1\\
 1& 1&-1\\
 1&-1&-1\\
 1& 1& 1
\end{pmatrix}                                                       \tag{OB.5}
```

has parent cap 17.  Since the independently certified global value is
`M_11=17`, this proves `F_C(A_0)=17`.

For class 1, deterministic CP-SAT proves that cap 17 is infeasible (125,241
branches and 14,793 conflicts in the frozen run).  The bridge

```math
B_1=\begin{pmatrix}
 1& 1& 1\\
 1& 1&-1\\
-1& 1&-1\\
 1& 1& 1\\
 1&-1& 1\\
 1& 1&-1\\
 1& 1&-1\\
 1&-1&-1
\end{pmatrix}                                                       \tag{OB.6}
```

has cap 19.  Order-eleven energies are odd, so the infeasibility at 17 and
this witness give `F_C(A_1)=19`.  The negative-triangle coordinate gives the
same pair `(17,19)`, as self-complementarity of both order-eight classes also
predicts.

This is an exact computational separation modulo the stated solver-certificate
qualification.

## 5. Archive comparison after freezing

The order-seven coordinate is not new.  `artifacts/scale_transfer_profile_no_go.md`
already proves, by exhaustive one-row insertion, that two order-seven
minimizers with the same signed magnetization-extrema profile have optimized
extension caps 12 and 10.  The present run recovers that obstruction in the
authoritative orbit quotient and shows that the third order-seven class joins
the cap-10 response cell.

The archive also already distinguishes the two order-eight classes under a
*fixed universal double*: their caps are 40 and 32 despite equal complete
energy histograms.  That is a fixed-composition separation, not (OB.1).
The present finite result is different: it allows the bridge to be freely and
separately optimized, removes child switching gauge, finds no distinction for
queries of orders one or two, and then finds an exact distinction at order
three.  No earlier archived optimized-query table of this form was found.

This survives the new gauge scope correction: optimized bridges erase a
switching label, but they do **not** erase all information distinguishing
inequivalent minimizer classes.

## 6. Research judgment

The finite evidence ranks the residual optimized bridge/holonomy fibre above
pure switching-orbit response as a potentially meaningful state.  It also
gives a useful minimal finite falsifier: any proposed gauge-invariant state
that identifies the two order-eight minimizer classes cannot predict all
optimized three-vertex continuations to additive error below one.

What it does **not** show is more important:

1. the gap two is fixed, not macroscopic;
2. only two order-eight classes are available, so there is no packing-rate
   lower bound;
3. query order three may be an isolated finite threshold;
4. evaluating (OB.1) still performs full bridge optimization and is not a
   compressed composition theorem.

Therefore this result warrants the residual fibre as a sharply defined next
finite/asymptotic question, but it does not establish scaling, convergence,
or even an asymptotic information lower bound.
