# Wave 29 root route: exact common-coset minimax

Status: the finite-game identity is verified.  The `A_9` values below have
exact rational primal and dual certificates after exhaustive enumeration;
they are finite evidence only.

Fix `(A,m,t)` and any finite family `mathcal A` of block cosets whose whole
support satisfies the chosen row-square cap.  For `a in mathcal A`, let

```math
H_a=\{S:Q(P_a^T D_a H_S D_a P_a)\ge Y_A(S)-t\}.
```

The best common law in (10.838) has value

```math
\delta_*=
\max_{\mu\in\Delta(\mathcal A)}\min_S\mu\{a:S\in H_a\}
=\min_{w\in\Delta(\binom{[n]}m)}\max_{a\in\mathcal A}w(H_a).
```

This is finite LP duality.  Equivalently, if `tau_*` is the minimum total
weight of a fractional cover `(u_a)` satisfying
`sum_{a:S in H_a}u_a>=1` for every selector, then
`delta_*=1/tau_*` (with the usual zero/infinity conventions).  Thus the
missing asymptotic theorem has a sharp adversarial formulation: every
selector law `w`, not just the uniform law, must put
`exp(-O(n^(3/4-c)))` mass in one row-good coset hit set.

Pointwise planting alone gives only the trivial exponential guarantee.  If
each selector `S` has one assigned hit coset `a_S`, the law induced by a
uniform selector has `mu(H_S)>=1/binom(n,m)` because it contains the atom
`a_S`.  Its negative log is `Theta(n)` at fixed density, rather than the
required `O(n^(3/4-c))`.  Any proof through this game must therefore establish
genuine overlap against arbitrary selector weights; self-collision is not
enough.

The exhaustive checker `common_coset_minimax_r29.py` enumerates all translates
of all `3+3+3` block subgroups of the exact minimizer `A_9`, retains only
cosets whose entire projective spin support is under the row cap, constructs
the threshold-zero incidence matrix, and solves both sides.  Because every
hit margin is separated from zero and the zero-one LP solutions rationalize
and verify exactly, it gives:

| instance | good cosets | deterministic full cosets | best one-coset coverage | uniform-law minimum | exact game value |
|---|---:|---:|---:|---:|---:|
| `A_9,m=5,C=80` | 4342 | 0 | 124/126 | 1816/4342 | 839/995 |
| `A_9,m=6,C=80` | 4342 | 0 | 83/84 | 1472/4342 | 21/22 |

Commands:

```bash
.venv/bin/python tmp/common_coset_minimax_r29.py A9 5 80 3 3 3 \
  --expect-fraction 839/995 --expect-cosets 4342 --expect-full 0
.venv/bin/python tmp/common_coset_minimax_r29.py A9 6 80 3 3 3 \
  --expect-fraction 21/22 --expect-cosets 4342 --expect-full 0
```

These examples show that fractional multi-coset mixing can remove almost all
of a deterministic compressed-alignment failure.  They neither give an
asymptotic lower bound nor overcome the arbitrary-selector dual.  The exact
next lemma is a minimizer-specific statement that for every selector law `w`
there is a balanced row-good coset with
`w(H_a)>=exp(-O(n^(3/4-c)))`, uniformly in one active ratio window.
