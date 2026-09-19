# Native chiral control versus its joint fractional resource

Status: exact identities and verified limitations; the proposed native bound
below is OPEN. These observations do not count as a completed construction.

## 1. An exact bilinear reformulation and its limitation

Let D=[[A,C],[C,-A]], where A is symmetric hollow and C is symmetric
(its diagonal may contain matching signs). For J subset[n], I=J^c, define

```math
 C_J=\begin{pmatrix}C[I,I]&A[I,J]\\
 A[J,I]&-C[J,J]\end{pmatrix},\qquad
 \beta(C_J)=\max_{u,v\in\{\pm1\}^n}|u^TC_Jv|.
```

With this partition order,

```math
 Q(D)=\max_J\beta(C_J).                            \tag{1}
```

Proof: quarter-rotate the coordinate pairs in J by
x_i=y'_i,y_i=-x'_i. The resulting full matrix is still chiral and has
bridge C_J. Global reversal of either spin block gives the general bridge
lower bound Q(D)>=beta(C_J). Conversely choose a maximizing parent spin
(x,y), and set J={i:x_i=-y_i}. In the rotated coordinates both halves are
the SAME spin x, so the internal energies cancel and the value is
x^T C_J x. This attains Q(D) in absolute value. Both inequalities prove(1).

In particular a parent ground-state gauge gives beta(C_J)=Q(D), attained
with equal left and right spins. After a simultaneous switch, these can
be made all-positive. The resulting bridge row sums c=C_J1 and internal
row sums a=A_J1 satisfy c_i>=|a_i|: each single-coordinate flip of either
parent half cannot increase the positive maximum. More generally,

```math
 |H_{A_J}(u)-H_{A_J}(v)|
 \le Q(D)-|u^T C_Jv|\quad\hbox{for every Boolean }u,v. \tag{2}
```

Equation(2) is equivalent to the full cap assertion, not a smaller sufficient
state. Equation(1) still maximizes over2^n cuts and full bilinear responses.
We do NOT count either as a strict reduction or an efficient algorithm.
Pair rotations need not retain the restricted relation B~A, as the exact
order4 counterexample in the adversarial follow-up shows.

## 2. What the sharp diamond inequality actually supplies

For |a_i|+|b_i|<=1 and |c_i|+|d_i|<=1 put

```math
 \Delta(A)=\max|H_A(a)-H_A(c)+b^TAd|.
```

The uniform track proves Delta(A)<=4beta(A)/3, with a14-term rational
certificate and a sharp weighted witness. The certificate was discovered
by a finite linear program, but its proof is an integer matrix identity
with no solver assumption. The director reran that identity and its
1900 lifted integer checks; the adversarial track independently checked
every scalar identity and the vertex/gauge reduction.

This resource has a physical implementation on balanced twins. It is NOT
the native optimized twisted cap F(A). Indeed the order6 child has
Delta=12 while its entire native core family has minimum18. More strongly,

```math
 \Delta(A)\le n\|A\|_{op},                         \tag{3}
```

by bounding the quadratic terms with half their squared Euclidean norms,
and b^TAd by ||A||op||b||||d||. The diamond constraints give
||a||²+||b||²<=n and ||c||²+||d||²<=n; the remaining difference is
-||A||op(||b||-||d||)²/2.

For symmetric conference seeds, (3) is <=n sqrt(n-1). Thus a universal
native rounding theorem F(A)<=Delta(A)+o(n^(3/2)) would force parent
normalized caps <=1/(2sqrt2), contradicting the archived original lower
bound (any verified lower constant greater than1/(2sqrt2) suffices).
The resource-to-native loss is genuinely leading-order, not just parity.

Completed2-clone children are themselves excluded from the near-minimizer
regime: Q(clone(A))>=4Q(A)-n, so their normalized liminf is at least
sqrt2 times the original lower constant, greater than the known upper.
An o(n^(3/2))-norm perturbation cannot close that gap.

## 3. A precise stronger native question, still unresolved

The logically different statement

```math
 F(A)\le\tfrac43\beta(A)+o(n^{3/2})                \tag{4}
```

is not falsified by (3). Its ZERO-error core version is false already
at order6: beta=12 but the minimum core cap is18>16. No O(n)-error
version follows from the diamond certificate or its sharpness example.

If(4) held uniformly on symmetric conference seeds at all sufficiently
large admissible orders, beta(A)<=n sqrt(n-1) would give parent normalized
cap at most2/(3sqrt2). Paley conference orders are ratio-dense, by the
prime number theorem in the fixed progression1 modulo4, and principal
restriction fills the intervening parent orders. Thus this would be a
genuine new all-order upper bound, not yet a convergence proof. The
Paley/PNT normalization and realization are already reconstructed in
ar_design_quasirandom_literature_toolkit.md.

Conversely the archived parent lower constant c_- implies that ANY
infinite full-sign symmetric seed family with

```math
 \limsup \frac{\beta(A_n)}{n^{3/2}}
 <\frac{3\sqrt2}{2}c_-\simeq0.919226
```

would falsify(4), independently of the twist. Finite small bilinear norms
are not such a family. Tensor amplification cannot be assumed to preserve
the normalized bilinear norm.

Current attempts remove neither obligation. Coarse block approximation
gives only O(n²/sqrt(log n)) error, larger than the required scale.
The retrieved Gale--Berlekamp literature supplies no applicable infinite
symmetric full-sign family below the displayed threshold. This is a
bounded diagnostic, not a new unconditional bound or an automatic campaign.

## 4. Source-level reconstruction of the key finite claim

The director reread the independent order10 census and the frozen parent
evaluator, rather than using audit verdicts as proof steps. The census
fixes the root row by switching, leaving36 signs; its two completion
records cover disjoint halves of all2^36 masks. Every accepted mask is
checked against all512 projective spins using independent signed half-edge
tables. Both shards see minimum13 and181440 members of each of the two
classes. The other census uses different energy arithmetic.

The parent search enumerates all rooted permutation copies and all root-fixed
switches. Signed automorphisms/antiautomorphisms reduce the search, with
orbit-size weights summing to185794560 for each order10 child. Antiautomorphisms
are legitimate because exchanging the two parent halves restores A. All
matchings are tested by exact Hamming-distance masks. Parent cap parity is
fixed, so exclusion at bestcap-2 rules out every smaller cap. The direct
profile recurrence, distance cutoffs, and representative restriction
x0=y0=+1 were checked algebraically. Each full completion checks9292808
representatives and returns44; explicit witnesses are independently evaluated.

This is an exact computational classification, not a formal proof-assistant
certificate or a newly established global order20 optimum. The recorded full
runs were not redundantly rerun during the closing source inspection.
