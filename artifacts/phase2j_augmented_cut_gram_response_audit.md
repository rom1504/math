# Augmented-cut Gram and conditional-response audit

Date: 2026-07-31. This is an agent-authored exact finite audit of candidate
composition states. It gives several rigorous falsifiers and one finite
survivor, but no uniform composition lemma.

## 1. The augmented-cut marginal already fails

For a signing `A`, distances from its edge word to the augmented cut code are
an affine rewriting of the Boolean energy multiset. The two exact order-eight
minimizer classes have the same complete energy/distance enumerator:

```text
-10:4, -8:10, -6:12, -4:16, -2:16, 0:12,
  2:16,  4:16,  6:12,  8:10, 10:4.
```

Nevertheless their fixed universal doubles have exact caps 40 and 32. Thus
the one-coset augmented-cut distance enumerator does not predict even this
specific finite composition operation. The failure is exact and occurs on
two genuine minimizer classes.

## 2. Gram-response candidates

For a bridge `C` between children `A` and `B`, uniform averaging gives

```math
\mathbb E_y(x^{\mathsf T}Cy)^2=x^{\mathsf T}CC^{\mathsf T}x,
\qquad
\mathbb E_x(x^{\mathsf T}Cy)^2=y^{\mathsf T}C^{\mathsf T}Cy.             \tag{GR1}
```

This motivates three progressively richer exact states:

1. the full matrices `(CC^T,C^TC)`;
2. their Boolean quadratic-value distributions conditioned on the internal
   child energies;
3. for every energy-shell pair `(e,f)`, the exact conditional second moment
   `sum_(x,y)(x^TCy)^2`.

Every sign bridge was exhaustively enumerated for both child orientations in
the cases `3+3`, `3+4`, `3+5`, and `4+4`: respectively 512, 4,096, 32,768,
and 65,536 bridges per orientation. Parent caps were evaluated exactly over
all projective spin pairs.

All three states fail. The smallest explicit collisions, for the positive
`3+3` case, are:

### Full Gram pair

```text
C_13 =  [ 1 -1  1]       cap 9
        [ 1 -1 -1]
        [-1 -1 -1]

C_66 =  [-1  1 -1]       cap 11
        [-1 -1 -1]
        [ 1 -1 -1]
```

These bridges have exactly the same `(CC^T,C^TC)`. Their common-state hash is

```text
7f11081ad7082d3b12c03906ea2a53d16c394a71e59a87118692118d9f0617b0
```

Thus even the **full** two-sided bridge Gram data misses composition-relevant
phase alignment with the child energy shells.

### Energy-conditioned marginal variances

Bridge codes 10 and 11 have the same two conditioned marginal variance
profiles but caps 11 and 9. The common-state hash is

```text
9e5d2ab69804093cc3e55b4befefad4d4b07822add3732112d8bed4cbcd766e1
```

### Energy-pair-conditioned second moments

Bridge codes 78 and 85 have the same conditional second moment in every
internal-energy-shell pair but caps 11 and 9. Their common-state hash is

```text
a408ea28af1b7dbca06907ab7a0324b5f31ca917e1a70b4dd49484e7cdbf9ce5
```

Independent collisions recur in every exhaustive held-out case. These are
not solver outcomes or sampled coincidences; the JSON certificate contains
both bridge matrices, state hashes, and exact caps.

## 3. A finite survivor: conditioned second and fourth moments

Adding the conditional fourth sum

```math
\sum_{H_A(x)=e,\,H_B(y)=f}(x^{\mathsf T}Cy)^4                         \tag{GR2}
```

to each shell pair produced no same-state/different-cap collision in the
preceding exhaustive cases. An additional held-out test exhausted all
`2^18=262144` bridges for `3+6` with positive child orientation. It yielded
484 distinct second-plus-fourth states and again no collision.

This is reproducible finite prediction, not a theorem. It should not be
promoted for three reasons:

1. no identity shows that two moments determine the maximum cross response;
2. a moment-to-maximum bound pays for the number of states in the shell and
   therefore has an exponential entropy loss unless new structure is proved;
3. the state is not smaller than the exact shellwise maximum envelope in
   coordinate count.

If the children have `L_A,L_B` energy levels, the second-plus-fourth state
stores two nontrivial integers in each of `L_A L_B` bins. The exact envelope

```math
R_C(e,f)=\max_{H_A(x)=e,\,H_B(y)=f}|x^{\mathsf T}Cy|                  \tag{GR3}
```

stores one. In the exhaustive cases the numbers of distinct observed
second-plus-fourth states versus envelope states were

```text
3+3: 18 vs 18
3+4: 151 vs 148
3+5: 336 vs 261
4+4: 2752 vs 2720.
```

The `3+6` second-plus-fourth test had 484 states; no asymptotic inference is
made.

## 4. Exact envelope and missing uniform lemma

Because global negation of one child preserves its internal energy and flips
the cross term, the parent cap is exactly

```math
\max_{e,f}\bigl(|e+f|+R_C(e,f)\bigr).                                \tag{GR4}
```

The program verifies (GR4) for every enumerated bridge. Hence the vector
`R_C` predicts finite composition perfectly and survives all held-out
minimizer cases by identity. But it is a repackaging of the full bridge
obligation: constructing `C` with

```math
R_C(e,f)\le T-|e+f|\quad\text{for every }(e,f)                        \tag{GR5}
```

is equivalent to satisfying every original parent spin constraint. No
obligation has been removed.

A genuine composition advance would require a uniform theorem bounding
(GR3) from a smaller state—perhaps conditioned fourth moments with no shell
cardinality loss. The exact collisions above prove that Gram data and second
moments cannot supply such a theorem alone. The current fourth-moment survival
is too finite and too large-dimensional to justify that conjecture.

## Reproduction

```bash
PYTHONPATH=computations .venv/bin/python \
  computations/phase2j_bridge_gram_response_collision.py \
  --output computations/results/phase2j_bridge_gram_response_collision.json
```

The exact result and concise log are:

- `computations/results/phase2j_bridge_gram_response_collision.json`;
- `computations/results/phase2j_bridge_gram_response_collision.log`.
