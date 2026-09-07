# Joint complex-cut completion: active direct-transfer arm

2026-09-07. Status: exact interfaces and rejected shortcuts only. No new
asymptotically lossless selectable completion has been proved in this note.

The actual target remains a value transfer starting from `M_n`, not merely
another bounded-cap family. Global old-edge changes and selection among
optimal children are allowed. The restricted skew lift below is a sufficient
route, not a necessary form of every possible solution.

## 1. The skew completion really needs a joint magnitude bound

For `L_-=[[A,C],[-C,-A]]`, `A` symmetric hollow signs and `C` skew hollow
signs, write `x=z+w, y=z-w`, with disjoint signed supports `z,w` partitioning
the coordinates. Then

```
H_{L_-}(x,y)=2 z^T(A-C)w,
Q(L_-)=2 max_{z,w} (|z^TAw|+|z^TCw|).
```

The last equality uses cut reversal. Thus a sufficient selectable estimate
is

```
min_{Q(A)=M_n, C skew signs}
 max_{z,w} (|z^TAw|+|z^TCw|)
 <= sqrt(2) M_n + O(n^(3/2-eta)).
```

One must not replace this by the stronger Euclidean-completion requirement

```
max_{z,w} sqrt((z^TAw)^2+(z^TCw)^2) <= M_n+o(n^(3/2))
```

without a separate proof. Cauchy--Schwarz makes the latter sufficient; it
does not make it equivalent. Nor does a spectral completion, a Gaussian
covariance choice, or cancellation on a selected extremal cut establish the
displayed joint magnitude estimate.

The finite selectable computation is preserved separately in
`principle_synthesis_2026_09_07_selectable_anti_lift.md`. The actual all-order
subhalf anti-family theorem is in
`principle_synthesis_2026_09_07_anti_invariant_weave.md`; it does not land an
arbitrary low optimal value as the first-half child value.

## 2. Exact four-state edge insertion, and the remaining first response

Use the directed-cut form `H=4 sum_{i->j} a_ij z_i w_j`, where each site has
one nonzero component, equal to a sign. For a cavity with one original edge
deleted, let

```
p = Pr(z_i != 0, w_j != 0),
q = Pr(w_i != 0, z_j != 0),
s = E[z_i w_j],       t = E[w_i z_j].
```

At inverse coupling `lambda`, the four choices of oriented signed edge
multiply the cavity partition by

```
1+(cosh(4lambda)-1)p +- sinh(4lambda)s,
1+(cosh(4lambda)-1)q +- sinh(4lambda)t.
```

If both orientation and sign are freely optimized for this cavity, the
optimal increment before applying the logarithm is exactly

```
min{(cosh(4lambda)-1)p-sinh(4lambda)|s|,
    (cosh(4lambda)-1)q-sinh(4lambda)|t|}.
```

This formula is not an edgewise characterization of the selectable problem
when `A` is constrained to stay an exact child optimizer: flipping its sign
can leave that constraint set. It is an exact unconstrained cavity interface.
Even there, anti-invariance of the complete energy range does not remove the
absolute first-response terms or force `p=q`. Averaging an optimizer orbit
may symmetrize averages, but does not justify a pointwise cavity condition.

## 3. The larger signed-quarter-turn family is already in the archive

There is another exactly balanced full-sign family,

```
L_+ = [[A,C],[C,-A]],
```

where `C` is symmetric with every entry, including its diagonal, a sign.
The signed permutation `J=[[0,I],[-I,0]]` has `J^2=-I` and
`J^T L_+ J=-L_+`. There are no missing matching edges. Its exact interface is

```
H_{L_+}(x,y)=2z^TAw+z^TCz-w^TCw.
```

Thus the second matrix contributes internal-shore energies, whereas it
contributes another cross-shore response in the skew family. For
`u=x+iy`, the same identity is
`H_{L_+}=Re[u^T(A-iC)u]/2`.

This representation is **not new**: it is explicit in
`dependent_lift_analytic_audit.md`, `chiral_scale_preserving_lift.md`, and
`decisive_audit_same_spin_unitary_phase_geometry_2026_09_07.md`. Those notes
already preserve finite dependent successes, the exact-compression centering
obstruction, a specific Clifford-family failure, and the non-lossless
continuous-phase-to-Boolean conversion. None of those narrower results
proves or refutes an unrestricted asymptotically selectable symmetric-`C`
completion. The canonical choice `C=A+I` is much narrower and has separate
bounded-operator subhalf counterexamples in
`decisive_independent_bounded_operator_h2_counterexample_2026_09_07.md`.

Balancing the energy endpoints of all principal child blocks would still
not by itself settle this interface. An internal state replacement changes
the cross response `z^TAw`; separate endpoint data omit exactly the joint
query that the completion needs. This is the same quantifier issue exposed
by the earlier profile collisions.

## 4. Literature scope check, not a new imported bridge theorem

Balla--Hambardzumyan--Tomon,
[Factorization norms and an inverse theorem for MaxCut](https://arxiv.org/abs/2506.23989)
(published in 2026), gives Boolean inverse statements for bounded
factorization norm or bounded normalized trace norm. It does not supply the
needed compactness here. The existing fourth-moment cap control and Schatten
interpolation imply, for **every** signing with `Q(A)=O(n^(3/2))`,

```
||A||_* >= ||A||_F^3 / ||A||_S4^2 = Omega(n^(3/2)).
```

Consequently its normalized trace norm, and hence its factorization norm,
grow at least as `sqrt(n)`. Passing to a Boolean `0/1` encoding changes the
trace norm by at most a rank-one `O(n)` term, so it does not repair that
scope mismatch. Its small ordinary-MaxCut-surplus hypothesis is likewise
not our two-sided `n^(3/2)` discrepancy regime.

No cited external theorem in this section is being used to assert a new
cross-order result. The direct selectable completion remains an open proof
obligation within the current campaign.
