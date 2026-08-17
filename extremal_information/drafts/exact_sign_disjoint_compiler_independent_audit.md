# Independent audit: character-preserving disjoint compiler

**Verdict:** CORRECT AFTER ONE WORDING REPAIR.  Lemmas/Theorems CD.1--CD.4
and every displayed inequality are correct.  The supplied finite verifier
passes.  However, the prose after CD.3 and in CD.4 must say that the maximum
defect is **at least** `Omega(k^(3/2))`, not that it *is*
`Theta(k^(3/2))`; the same architecture can have a quadratic worst defect.

I ran

```text
python3 experiments/verify_exact_sign_disjoint_compiler.py
```

which reports `exact-sign disjoint compiler checks passed: 11157`.

## 1. CD.1: scale transfer

If a query `t` separates `z,z'`, applying CD.4 to each of the two responses
and using the triangle inequality loses at most `2eta`, giving CD.6.  If
`N/k` tends to infinity, an `O(k^(3/2))` unscaled signal divided by
`N^(3/2)` is `O((k/N)^(3/2))=o(1)`.  Thus both directions of the stated
scale test are correct.

The edge-ownership convention is compatible with EL.1: the old block is
child-owned, the rank-one bridge is query-owned, and the new clique is
public.  The query-only calibration cancels in contextual distances but not
in an uncalibrated parent cap, exactly as stated.

## 2. CD.2: exact character-pullback classification

Fix output coordinate `0` and let `F_a=E_(0a)` for `a!=0`.  Since
`phi_0^2=1`, uniqueness of Boolean Fourier characters gives

```math
E_(ab)=F_a\mathbin\triangle F_b.
```

The `F_a` are distinct: equality would make `E_(ab)` empty, contrary to
CD.8.  Since each `F_a` and every `F_a triangle F_b` has cardinality two,
each pair of sets intersects in exactly one element.

For `k>=5`, there are `k-1>=4` such sets.  A pairwise-intersecting family of
at least four distinct two-subsets is necessarily a star.  Indeed, after
`{u,v}` and `{u,w}` are chosen, a set avoiding `u` must be `{v,w}`; those
three triangle edges admit no fourth distinct edge meeting all three.  Hence

```math
F_a=\{u,v_a\},
```

with all `v_a` distinct.  The `k-1` leaves exhaust the other input
coordinates, so `a -> v_a`, together with `0 -> u`, is a permutation.

For signs, set `s_0=1` and `s_a=epsilon_(0a)`.  Multiplying the two base
identities gives `epsilon_(ab)=s_as_b`.  Defining

```math
g(x)=phi_0(x)/x_u
```

then gives `phi_0=gx_u` and `phi_a=gs_ax_(v_a)`, hence CD.9.  The arbitrary
global gauge is the only freedom not visible to quadratic pair characters.

The restriction `k>=5` is genuine for this proof and statement: with only
three nonbase sets, the triangle family is possible.  The supplied verifier
checks the finite star classification at `k=5`; the preceding combinatorial
argument proves all larger orders.

## 3. CD.3: constants and the global antipode

Pointwise,

```math
||R^Tx||_1=\max_y x^TRy\ge |x^TRphi(x)|,
```

because both `phi(x)` and its global antipode are admissible new-shore
states.  Thus the defect is nonnegative and the absolute value exactly
accounts for the gauge `g(x)` in CD.9.  With `C=RDP`, the intended reward is
`|x^TCx|`; right multiplication by a signed permutation preserves complete
unit signs.

For uniform `x`, each of the `k` bridge columns has a length-`k`
Rademacher sum, so the sharp `p=1` Khintchine inequality gives

```math
E||R^Tx||_1\ge k\sqrt{k/2}.
```

The exact Fourier expansion is

```math
x^TCx=tr(C)+\sum_(i<j)(C_ij+C_ji)x_ix_j.
```

Consequently

```math
E(x^TCx)^2
=(tr C)^2+\sum_(i<j)(C_ij+C_ji)^2
\le k^2+4{ k\choose2}<3k^2.
```

Cauchy--Schwarz, subtraction of expectations, and
`max Delta>=E Delta` prove CD.13 with the stated constants.  No independence
between columns is used, and no additional antipodal assumption is hidden.

The required wording repair is here.  CD.13 proves

```math
max_x Delta_(R,phi)(x)
\ge(1/\sqrt2-o(1))k^{3/2},
```

not equality or an `O(k^(3/2))` upper bound.  For example, take `R=J` and
`phi(x)=x`.  Writing `s=|sum_i x_i|` gives

```math
Delta(x)=ks-s^2,
```

whose maximum is `k^2/4+O(1)`.  Therefore the sentence “the worst locking
defect is ...” should read “the worst locking defect is **at least** ...”.
The no-`o(k^(3/2))` conclusion remains fully valid.

## 4. CD.4: replication and multiplicities

Let `n_j=|pi^(-1)(j)|`.  The proof uses exactly the declared assumptions:

```math
1\le n_j\le L,\qquad \sum_jn_j=m.
```

The lower bound does not need the first inequality, although surjectivity is
natural if every old coordinate is meant to be represented.  With

```math
C_ij=\sum_(a:pi(a)=j)R_ia s_a,
```

one has

```math
||C||_F^2\le k\sum_jn_j^2\le kLm,
\qquad |tr C|\le\sum_jn_j=m.
```

For a nonsymmetric `C`, orthogonality gives

```math
E(x^TCx)^2
=(tr C)^2+\sum_(i<j)(C_ij+C_ji)^2
\le m^2+2||C||_F^2
\le m^2+2kLm.
```

Each of the `m` bridge columns still contributes at least `sqrt(k/2)` to
the expected unrestricted roof.  This proves CD.17.

When `m=Theta(k)` and `L=O(1)`, surjectivity and the fibre bound also imply
`k<=m<=Lk`.  Hence the first term of CD.17 is `Theta(k^(3/2))` while the
square-root subtraction is only `O(k)`.  The theorem therefore proves an
`Omega(k^(3/2))` defect.  It does **not** prove a matching upper bound; the
same all-ones example can again be quadratic.  Replace “the defect is still
`Theta(k^(3/2))`” by “the defect is still
`Omega(k^(3/2))`” or “the defect is at least of order `k^(3/2)`.”

CD.4 is a separate replicated-coordinate architecture, not literally an
extension of CD.2's pair-character hypothesis: two outputs in the same
fibre have constant product rather than a two-coordinate character.  Its
scope is nevertheless clear and its proof is correct.

## 5. Scope relative to EL.1

There is no conflict with the positive coordinate-pin compiler EL.1.

- EL.1 chooses a **query-dependent** rank-one bridge and needs to expose one
  coordinate response per context.  It does not ask one fixed bridge to
  realize the entire pair-character algebra on its new shore.
- CD.3 fixes one complete bridge and assumes a state encoding whose every
  pair character pulls back exactly.  It rules out optimizing that bridge
  separately as an `o(k^(3/2))` equality lock.
- CD.3 does not cover a child-dependent active witness, correlated
  auxiliary interactions, a narrower structured query family, or joint
  cancellation before the bridge reward is paid.

Thus the hierarchy in Section 4 is logically consistent.  After replacing
the two equality/order phrases by lower-bound language, the draft is
theorem-level and accurately delimits the remaining joint same-switch
problem.
