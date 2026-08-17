# Independent audit: mesoscopic bipartite-core third witness

**Verdict: PASS after three required repairs.**

This audit freezes
`mesoscopic_bipartite_core_third_witness.md` at SHA-256
`2e642cf5f4d713e21f6cf7628f54982638218de400d2e9c50762845b1ed2b37f`
and checks its conventions against the repaired PP.4 draft at SHA-256
`01e1080a6b5efaceb73ceabea776ee96378947d736c781a5dde60fb7ce189532`.
The main conditional conclusion is correct: two nearly antipodal positive
shell poles with a common-correct reservoir of size at least `M-2s` force a
third projective shell line, at deficit `O(n log^2 n)`, whose projective
distance from both poles is `(1/2-o(1))M`.  This remains only
`Theta(n^(-1/2))` in normalized edge distance and does not prove Cut-DH(3).

The frozen statement has two correctable logical defects and one material
scope overclaim.  They do not affect Corollary MB.3 after it is stated with
the intended choice of `r`.

## 1. PP.4 gauge and the common-correct reservoir: PASS

In the PP.4 gauge,

```math
z_0=\mathbf1,\qquad z_1=-c(v_S),\qquad P=\delta(S).
```

Both signed words equal `+1` on `P` and are opposite off `P`.  With

```math
Z=\{e\in P:a_e=+1\},\qquad B=\{e\in P:a_e=-1\},
```

the code dictionary indeed makes `Z` common-correct and `B` common-wrong.
Moreover

```math
C=\sum_{e\in P}a_e=|Z|-|B|
 ={\langle a,z_0\rangle+\langle a,z_1\rangle\over2},
```

so `C>=M-2s`, and

```math
p=|Z|={D+C\over2}\ge C\ge M-2s.
```

Here `C>=0` uses the hypothesis that the two poles are positively oriented;
it is not implied by the numerical lower bound when `2s>M`.  The draft's
terminology supplies this hypothesis, but stating it explicitly would avoid
an ambiguity.

The complemented-cut convention is also correct.  PP.4 assumes actual code
distance `E-D` and projective distance `D<E/2`; the agreement coordinates
of the two **oriented sign words** are the `D` coordinates in `P`.

## 2. Lemma MB.1: PASS

For `F subseteq Z`, every flipped coordinate has `a_e=+1`.  If the
maximizer of `A^F` has `q` negative entries on `F` and
`d=M-\langle a,z\rangle`, then exact global minimality gives

```math
M\le Q(A^F)=\langle a^F,z\rangle
=M-d-2r+4q.
```

Consequently `q>=r/2+d/4`; because `0<=d` and `q<=r`, also `d<=2r`.
All signs and factors of two agree with the exact flip certificate (2.3)--
(2.5) in `nearmin_deterministic_inequalities.md`.  This lemma is precisely a
localized use of that archived certificate, not a new optimality principle.

## 3. Hypergeometric selection in MB.2: PASS

For a fixed augmented cut with at most `theta p` negative coordinates on
`Z`, sampling an `r`-subset without replacement gives

```math
\Pr\{X\ge r/2\}
\le \exp[-2(1/2-\theta)^2r].
```

This is the standard Hoeffding comparison for hypergeometric sampling.  A
union bound over at most `2^n` augmented cuts is legitimate and, under
(MB.7), produces one `F` simultaneously excluding every concentrated cut.
The adaptively selected maximizer from MB.1 has `q>=r/2`, so it must have
more than `theta p` negative entries on all of `Z`.  Strict versus weak
inequalities and odd `r` cause no rounding problem: `X>=r/2` means
`X>=ceil(r/2)`.

This is a genuine geometric deduction beyond FB.3's cardinality statement,
although it uses the same exact-flip primitive.

## 4. Required repair 1: the selected orientation need not be positive

Theorem MB.2 currently calls the selected `z` a **positive** signed-cut word
without assuming `2r<M`.  MB.1 proves only

```math
\langle a,z\rangle=M-d\ge M-2r.
```

Thus the actual maximizing orientation can have negative original energy
when `2r>=M`.  Reorienting it may destroy the signed-core conclusion
(MB.8): if `z` is negative on `t>theta p` coordinates of `Z`, then `-z` is
negative on `p-t`, which need not exceed `theta p`.

There are two clean repairs.

1. Add `2r<M` to MB.2.  This is already true in both asymptotic uses in
   Corollary MB.3.
2. More generally, state that MB.2 first produces an *oriented response*
   `z` satisfying (MB.8).  Its projective line has a positive representative
   in the deficit-`2r` shell: if `d>=M`, then the positive representative
   `-z` has deficit `2M-d<=M<=2r`.  Projective distances are unchanged, but
   (MB.8) belongs to the original response orientation, not necessarily to
   that positive representative.

The sentence after (MB.10) notices the issue but does not repair the theorem
as written.  Corollary MB.3 is unaffected because its choices have
`r=o(M)`, hence `2r<M` eventually.

## 5. Projective-distance calculation: PASS

Let `h_j=d_E(z,z_j)`.  Since both poles are `+1` on `Z`, the selected
orientation gives `h_j>theta p`.  Coordinatewise, the sum of two energies
is at most twice their number of agreement coordinates, so

```math
2M-d-2s
\le\langle a,z\rangle+\langle a,z_j\rangle
\le2(E-h_j).
```

Therefore

```math
E-h_j\ge M-s-d/2,
```

and taking `min(h_j,E-h_j)` proves (MB.9).  Monotonicity under `d<=2r`,
`p>=M-2s`, and PP.4's `d_P(z_0,z_1)=D>=M-2s` gives (MB.10) with all
constants correct.

This complementary-distance inequality is exactly the two-word positivity
bound AO.20, with unequal deficits retained.  It certifies only the energy
scale `Theta(M)`, not edge scale `Theta(E)`.

## 6. Required repair 2: Corollary MB.3 over-quantifies `r`

Part 1 says that “any integer `r>8n log 2`” yields deficit `O(n)` and
distance `(1/4-o(1))M`.  That is false as written: one may take, for
example, `r=n^{4/3}`, which is not `O(n)`, or `r` comparable to `M`, for
which the term `M-s-r` need not be `(1-o(1))M`.

Replace this by a concrete choice such as

```math
r=\lfloor8n\log2\rfloor+1,
```

or require simultaneously

```math
8n\log2<r=O(n).
```

Then `r<=p` eventually, the exponent is strictly larger than `n log 2`,
the deficit is `O(n)`, and (MB.16) follows.  Part 2 is correct: with
`theta_n=1/2-1/log n`, its exponent is greater than `n log 2`, while
`r_n=O(n log^2 n)=o(M)` and `r_n<=p` eventually.  Hence (MB.18) follows.

## 7. Required repair 3: the proof does not use bipartite cut geometry

After PP.4 supplies two nearly antipodal positive poles, the proof uses only

```math
|Z|\ge M-2s,
```

exact finite-edge flip minimality, the fact that there are at most `2^n`
responses, and two-word positivity.  It never uses that `P=delta(S)`, that
`P` is complete bipartite, the shore size, or any cut-incidence identity.
Indeed the same theorem holds for any response family of cardinality at most
`2^n` with such a common-correct reservoir.  The generic PP.2 code is
compatible with the conclusion: its thin shell can contain three directions
at `Theta(M)` separation while its whole projective diameter remains `o(E)`.

Accordingly, the following phrases should be narrowed:

* the complete-bipartite identity does not explain the reservoir's **size**;
  two nearly antipodal positive words already do;
* PP.4's `delta(S)` identity is not a cut-specific input used by MB.2;
* the argument is a conditional, generic third-witness theorem applied to
  the PP.4 branch, not yet a theorem about mesoscopic shore structure.

The architecture-ceiling claim is justified only in the precise sense that
**MB.1 + this hypergeometric exclusion + AO.20 certify no more than
`Theta(M)` projective distance**.  Equation (MB.13) is a lower bound, not an
upper bound on the true complementary distance, so it does not prove that
every nonrecursive multi-edge argument must fail.  A further cut-specific
fact could still show `E-h_j=Theta(E)`.  The generic PP.2 countermodel does
show that no strengthening follows from the generic premises alone.

Thus replace “decisive no-go for this nonrecursive proof architecture” by a
narrower statement: the displayed one-reservoir entropy/positivity proof has
reached its certifiable ceiling, and fixed projective scale requires a new
cut-specific or recursive input.

## 8. Final classification

| Component | Verdict |
|---|---|
| PP.4/code/sign conventions | PASS |
| MB.1 exact flip algebra | PASS |
| MB.2 hypergeometric union bound | PASS |
| MB.9 projective constants | PASS |
| MB.3 asymptotics | PASS after repairing the `r` quantifier |
| positivity wording | REPAIR REQUIRED |
| novelty claim | conditional geometric consequence; not cut-specific |
| fixed-scale implication | none; separation is `(1/2-o(1))M=Theta(n^{3/2})=o(E)` |
| architecture no-go | valid only for the displayed generic proof ingredients |

No counterexample was found to the repaired theorem.  I also exhaustively
checked the positivity-form conclusion over all 72 all-one-gauged exact
order-six minimizers for every applicable large-`r` parameter on a fine
`theta` grid; this finite check is diagnostic only and is not used in the
proof.  The rigorous verdict rests on the calculations above.
