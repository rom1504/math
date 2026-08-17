# Independent audit: mesoscopic affine scale-amplifier no-go

**Verdict: PASS AFTER TWO MATERIAL SCOPE REPAIRS.**

The contraction theorem and all displayed constants are correct.  The
task-local source report had two overstatements which matter to the research
frontier:

1. FR.5 proves a lower pairwise separation of order `M_n`, not an upper
   `O(M_n)` diameter.  Therefore the no-go applies to the branch whose source
   pairs remain `o(n^2)` apart; it does not prove that the whole newly
   constructed FR.5 family lies in that branch.
2. Escaping the no-go need not broadcast the label specifically into the
   old--old block.  It may broadcast it into any child-owned coefficient
   block before the common future is attached.  What is necessary is leading
   uniform distance between the *complete child landscapes*.

Both repairs are explicit in the canonical task-local draft
`mesoscopic_scale_amplifier_no_go.md`.  No canonical theorem, axiom, or
frontier file was changed in this audit.

## 1. Frozen inputs

```text
original task-local report
  /home/math/quadra/tmp/scale_amplifier_falsifier_report.md
  sha256 b3ddaad94c7b3f576a0d3ad89b6708fa60822af0cbc036a964796d897c8465bf

repaired theorem draft
  extremal_information/drafts/mesoscopic_scale_amplifier_no_go.md
  sha256 fcd3f98c1c4cadb585bca2f26013fbf5c3ebaa34e36a4ff4751e67a588e6a28b
```

The repository sources used to reconstruct the claim were principally:

* `nearmin_absolute_overlap_physical_compiler.md`, especially AO.5--AO.7
  and the independently audited Theorem AO.2;
* `fractional_reservoir_localized_flip.md`, especially FR.3--FR.5;
* `orientation_visibility_threshold.md`, OV.1--OV.3;
* `bounded_cap_contextual_metric_compiler.md`, BCX.1--BCX.3;
* `universal_pin_cap_barrier.md`; and
* Axiom 61, the archived public-interaction non-amplification principle.

## 2. Boolean normalization and diagonal audit: PASS

For a hollow coefficient matrix `D`, the repository convention is

```math
 H_D(x)=\sum_{i<j}D_{ij}x_ix_j.
```

There is no factor of two in this edge-vector convention.  If
`z=\sigma c(x)`, then

```math
 \langle D,z\rangle=\sigma H_D(x),
```

and maximization over all augmented cuts is exactly

```math
 \max_z\langle D,z\rangle
 =\max_x|H_D(x)|=\|D\|_{\rm B}.
```

For a block parent, each old--new edge occurs once as `x^TBy`; the new
quadratic term is `H_C(y)`.  Thus no diagonal term or hidden factor two is
missing from SA.1--SA.4.

The argument actually allows real hollow old blocks.  Exact signs enter only
when identifying the result as an obstruction for physical parent signings.

## 3. Common-future nonexpansiveness: PASS

Let

```math
 F_D(x,y)=H_D(x)+J(x,y).
```

For the same `J`, pointwise

```math
 |F_D(x,y)-F_{D'}(x,y)|
 =|H_{D-D'}(x)|\le\|D-D'\|_{\rm B}.
```

The elementary inequality

```math
 |\max |F|-\max |F'||\le\|F-F'\|_\infty
```

proves SA.5.  The analogous inequalities for a one-sided maximum and a
minimum are identical.  Crucially, this proof never compares the two
optimizers.  The optimizer may change adversarially with the child, so the
claim does not hide a common-maximizer assumption.

The continuation may contain arbitrarily many auxiliary variables and an
arbitrary joint landscape.  It may also vary with the public query label
`q`: applying the same inequality separately to the common continuation
`J_q` gives SA.6.  What is forbidden is allowing the coefficients called
the query to depend jointly on both `q` and the child label; that would no
longer be a common continuation.

This is the exact Boolean analogue of the general statement in Axiom 61.
It is not a scalar-channel bound: cancellation inside `J` occurs before the
outer maximum and absolute value.

## 4. Literal `m`-shore Lipschitz law: PASS

For

```math
 G_\sigma(x)=\max_y\sigma\{x^TBy+K(y)\},
 \qquad |B_{ij}|\le1,
```

comparison at a fixed `y` gives

```math
 |(x-x')^TBy|\le2m d_H(x,x').
```

Because `K(-y)=K(y)`, the substitution `y\mapsto-y` proves
`G_\sigma(-x)=G_\sigma(x)`, giving the projective minimum and SA.8.
For the opposite outer channel,

```math
 G_-(x')=\max_y\{x'^TBy-K(y)\},
```

so the additional cost is exactly at most `2\|K\|_infty`.  The coefficient
`2m` and the internal-landscape coefficient `2` are both correct.

The evenness of `K` is essential here.  It holds for every hollow quadratic
new-spin block, but would fail if external linear fields were added.  It is
not needed for the stronger common-future nonexpansiveness lemma.

If two augmented cuts have a common sign label,

```math
 z=\sigma c(x),\qquad z'=\sigma c(x'),
```

and the projective vertex distance is `d`, their exact edge Hamming distance
is

```math
 h=d(n-d).
```

Since `d\le n/2`, one has `d\le2h/n`, hence the `4mh/n` bound in SA.11.
Two positive deficit-`s` energies lie in an interval of length `s`, so the
old-energy contribution to their designated-witness score difference is at
most `s`, proving SA.12 and the rearranged lower bound SA.13.

There is no claim here about optimized response: a remote old optimizer can
beat both displayed anchors.  That loophole is closed only in Section 4,
where the entire child landscapes are shown close.

## 5. Orientation and projective conversion: PASS WITH EXPLICIT SCOPE

The common augmented orientation is not cosmetic.  If
`z=\sigma c(x)` and `z'=\sigma c(x')`, then

```math
 h=d(n-d)\le\lfloor n^2/4\rfloor.
```

Writing `d_P=min(h,E-h)`, the largest possible excess of `h` over `d_P` is

```math
 2\lfloor n^2/4\rfloor-E=\lfloor n/2\rfloor.
```

This proves SA.19.  Thus same-orientation projective distance `o(n^2)`
implies actual augmented-word distance `o(n^2)`.

By contrast, words with opposite augmented orientations can have actual
distance `E-O(M_n)` while their projective distance is only `O(M_n)`.  The
affine terms `p z^u` and `p z^v` can then be a leading distance apart.  Such
orientation-sensitive behavior is not ruled out; OV.1--OV.3 identify the
separate internal-cap budget relevant to exposing it.

A growing family has a growing same-orientation subfamily because there are
only two values of `sigma`.  This justifies applying the no-go to a
same-orientation hard branch, but it does not prove that branch has
`o(n^2)` diameter.

## 6. Reconstruction of AO.6 and the affine collision: PASS

AO.6 states, uniformly over every augmented cut `z`,

```math
 \left|\langle b^u,z\rangle-
 \{(1-p)\langle a,z\rangle+p\langle z^u,z\rangle\}
 \right|\le\rho.
```

Subtracting the instances for `u` and `v` cancels the entire base term.
The two approximation errors cost `2rho`, while

```math
 |\langle z^u-z^v,z\rangle|
 \le\|z^u-z^v\|_1=2h_{uv}.
```

Maximization over augmented cuts is exactly the Boolean norm, so

```math
 \|b^u-b^v\|_{\rm B}\le2p h_{uv}+2\rho.
```

No orientation, optimizer, or relaxation assumption is used in this
inequality itself.

For AO.2,

```math
 p=\alpha n^{-1/2},
 \qquad
 \rho=O(\sqrt\alpha\,n^{5/4}+n).
```

At fixed `alpha`, `h=o(n^2)` makes `2ph=o(n^(3/2))`, and the uniform error
is also `o(n^(3/2))`.  At the narrower scale `h=O(M_n)=O(n^(3/2))`, the
affine term is only `O(alpha n)` and the AO sampling error is
`O(n^(5/4))`.  Applying SA.5 gives the same ceiling after every common
future, independently of its order or cap.

The deterministic extension is also correct.  If

```math
 \|b^u-((1-p)a+pz^u)\|_B\le\rho,
 \qquad Q(b^u)=O(n^{3/2}),
```

and `z^u` is positive for `a`, evaluation at `z^u` gives

```math
 Q(b^u)\ge(1-p)\langle a,z^u\rangle+pE-\rho
 \ge pE-\rho.
```

Hence `p=O(n^(-1/2))` when `rho=o(n^(3/2))`, and subtraction yields SA.24.
This extension assumes one common `p`.  A label-dependent mixing strength
could itself carry information and is not silently included.

## 7. The fractional-reservoir collision: original overclaim repaired

FR.5 proves a family of size

```math
 L_n=\Theta(\log n/\log\log n)
```

in an `o(M_n)` shell, with minimum pairwise projective distance

```math
 (1/4-o(1))M_n.
```

Its proof gives no matching upper bound on pairwise distance.  A newly
created optimizer can be macroscopically remote from the earlier anchors;
FR.11 only lower-bounds both a distance and its complement.  Therefore the
original report's unqualified statement that “the current energy-scale
shell packing cannot be amplified” was too broad.

The exact valid statement is SA.4:

* if a same-orientation FR subfamily has diameter `o(n^2)`, its AO affine
  response code is uniformly `o(n^(3/2))` under all common futures;
* if it contains a growing fixed-fraction projective packing, audited AO.2
  already supplies the desired leading contextual separation.

This is a useful dichotomy, not a proof that the first branch occurs.

## 8. Archive collision and escape audit

* **Axiom 61 / TC-type public-interaction ceilings.** SA.1 is the exact
  same-old-domain Boolean specialization.  It adds the sharp AO.6 metric
  calculation.
* **UP.1.** Universal pinning is stronger than needed here and forces
  quadratic cap.  SA.1 permits arbitrary child-dependent optimizer
  switching, so it does not rely on UP.1's universal-witness hypothesis.
* **BCX.2--BCX.3.** These do not contradict the no-go.  Their common-query
  response gap is already `Theta(n^(3/2))`; by SA.1 the switched old blocks
  must themselves be separated at that scale.  They do not pass through the
  collapsed AO affine metric.
* **OV.1--OV.3.** These concern opposite old orientations and show why an
  unrestricted orientation no-go would be false.
* **Scalar, polarization, and same-map rounding barriers.** None is used.
  The future is optimized jointly before applying one uniform contraction.

The theorem rules out only a public continuation appended to children that
are already close in uniform Hamiltonian norm.  It does not rule out:

1. fixed-scale shell geometry;
2. a non-affine state encoding with leading child distance;
3. state broadcast into child-owned old, auxiliary, or interface blocks;
4. jointly child--query-owned coefficients; or
5. an orientation mechanism paying the OV internal-cap cost.

Accordingly the repaired `L_broadcast` is phrased for complete child
landscapes, not only for old coefficient matrices.

## 9. Final pass/fail table

| Item | Verdict | Comment |
|---|---|---|
| Boolean norm and augmented-cut normalization | PASS | no missing factor two or diagonal |
| common-future nonexpansiveness | PASS | arbitrary optimizer switching included |
| bank of common queries | PASS | coordinatewise contraction; joint ownership excluded |
| `m`-shore Lipschitz constant | PASS | `2m` per old-spin flip |
| projective reduction | PASS | requires even `K`; same orientation handled exactly |
| opposite-channel term | PASS | exact extra ceiling `2||K||_infinity` |
| AO.6 subtraction | PASS | gives `2ph+2rho` exactly |
| AO fixed-`alpha` scaling | PASS | `h=o(n^2)` implies subleading child distance |
| deterministic affine extension | PASS | positivity and common `p` are explicit |
| claim about all FR.5 pairs | FAIL AS ORIGINALLY WORDED | repaired to a conditional diameter dichotomy |
| claim that only old-block broadcast escapes | FAIL AS ORIGINALLY WORDED | repaired to complete child-owned state broadcast |
| final architecture-level no-go | PASS AFTER REPAIR | affine mesoscopic branch only |

The repaired theorem is a rigorous scalable no-go for the live AO affine
implementation.  It is not a no-go for all low-cap near-minimizer response
packings and does not by itself settle the fixed-ambient shell question.
