# Self-audit: approximate residual-shell law

**Scope.**  This audit checks
[`approximate_residual_shell_law.md`](approximate_residual_shell_law.md)
against Theorems 17.1l, 17.1u, and 17.2.  It is not an independent audit.

## 1. ARS.1

The proof uses the stronger and unambiguous factorization radius

```math
\inf_{a,p}\|P-a\otimes p\|_\infty,
```

not an unspecified pairwise Hilbert row diameter.  Given
`T_v(k,j)=a_k+p_j+E_kj`, a left prefix takes a maximum of quantities whose
errors all lie in `[-epsilon,epsilon]`; it cannot amplify this interval.
This proves (ARS.6).  A max changes by at most a sup-norm perturbation, so
subtracting two maxima costs the factor two in (ARS.7).  The factor is not
silently folded into the projective metric.

The conclusion is only for the normalized terminal functional `R_z`.  It
does not approximate the accumulated max-plus scalar.  Therefore it neither
assumes nor proves the semiconjugacy in Theorem 17.2.  The suffix itself is a
finite state; a smaller code is reusable only under the explicitly stated
right-congruence condition.

## 2. ARS.2 and ARS.2a

For rank-one matrices,

```math
(a\otimes p)(b\otimes q)
=a\otimes q+\max_j(p_j+b_j),
```

so the cyclic closing term is `max_j(p_last(j)+a_first(j))`.  This verifies
the orientation in (ARS.12).  The cycle-coboundary criterion is exactly the
finite reward criterion already present in Theorem 17.1l and is labelled as
such; the note does not claim it as a new general theorem.

If each block factor changes entrywise by `epsilon_v`, successive tropical
multiplication changes the product by at most `sum epsilon_v`.  Max-plus
spectral radius is one-Lipschitz in the entrywise norm, since every directed
cycle mean changes by at most that norm.  Hence (ARS.15b) is sound.  A cyclic
walk decomposes into simple cycles, proving the `Delta` bound.  A terminal
remainder of fewer than `D` original letters contributes one fixed-alphabet
absolute boundary constant, not a second extensive term.

The state `(left type,right type,scalar)` is a genuine strict quotient for
rank-one products.  Merging only right profiles need not respect its scalar
composition law because the omitted quantity is the pairing with the next
left profile.

## 3. Sharp two-letter counterexample

The two right profiles are `(0,0)` and `(0,delta)`.  Their projective distance
is half the oscillation, `delta/2`; their midpoint has radius `delta/4`, and
the triangle inequality proves optimality.

The compatibility table is

```math
delta [[0,1],[1,1]].
```

Thus a cyclic word receives one `delta` on every adjacency except `AA`.
For a mixed word, if `k` is the number of cyclic `A`-runs, then
`N_AA=N_A-k`.  The proposed tolls `(delta/4,5delta/4)` have defect
`delta(k-t/4)`, whose absolute value is at most `delta t/4` because
`0<=k<=t/2`.  The same formula with `k=0` handles the two constant words.

For the lower bound, an error rate `d` on `A`, `B`, and `AB` implies

```math
|g_A|<=d,
|g_B-delta|<=d,
|g_A+g_B-2delta|<=2d.
```

The first two give `g_A+g_B<=delta+2d`, while the third gives
`g_A+g_B>=2delta-2d`; hence `d>=delta/4`.  This proves exact optimality, not
just failure of one guessed toll.

Every matrix and product is rank one, so its projective action is constant.
The positive rate therefore survives maximal forgetting and cannot be
attributed to a failure of contraction.  Both matrices are source-plus-target
binary field tables, hence a width-one Ising transfer benchmark with zero
spin-spin bond.

## 4. Small-shell universality

For every nonempty product and each fixed initial endpoint, changing the
terminal endpoint changes only the last edge and costs at most `alpha`.
Comparing in both directions bounds each row range by `alpha`.  Choosing its
midpoint as the free left scalar gives `rad(P;0)<=alpha/2`.  Ordinary
positive scaling commutes with both max and addition, so the entire word
response is exactly scaled by `alpha`.

This proposition asserts no exponential lower bound by itself.  Its precise
claim is that every bounded all-finite weighted-automaton response algebra
can be embedded in a fixed-radius one-profile shell.  Any complexity lower
bound for a chosen base alphabet transfers at the same relative distortion.

## 5. Support-core guardrail

For (ARS.24), direct multiplication gives `pT_a=pT_b=p`.  Since each matrix
is all finite, the common finite left eigenprofile has eigenvalue zero and
every word has spectral radius zero.  The threshold relations are

```math
R_a={(1,1)},\qquad R_b={(2,1),(2,2)}.
```

The descending one-context core goes `I -> {1} -> empty`.  The note draws no
semantic conclusion from this failure and uses it only to enforce the
separation between support leakage and scalar residual error required by
Theorem 17.1u.

## 6. Evidentiary boundary

The genuinely new conclusions relative to the canonical theorem file are:

1. the bounded-delay terminal-shell theorem without multiplicative
   contraction;
2. the exact `delta/4` reset-versus-rate separation;
3. the fixed-one-profile small-shell universality obstruction;
4. the resulting trichotomy between terminal residual, scalar compatibility,
   and witness support.

The rank-one multiplication identity, cycle-coboundary test, and entrywise
perturbation bound are classical ingredients.  They are not advertised as
new in isolation.
