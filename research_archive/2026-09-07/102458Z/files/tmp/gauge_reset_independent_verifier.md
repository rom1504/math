# Independent verification of Theorem 16.14

**Files checked**

* `extremal_information/theorems.md`, Theorem 16.14;
* `/home/math/quadra/tmp/gauge_reset_necessity_report.md`;
* `/home/math/quadra/tmp/gauge_reset_self_audit.md`;
* `extremal_information/experiments/verify_tropical_absorption.py`;
* the associated summary language in `README.md`, `axioms.md`, and
  `examples.md`.

**Verdict.** The mathematical claims in integrated Theorem 16.14 are
correct with their declared quantifiers.  The integration appropriately
implements the self-audit downgrade: the clamp refutes only the narrow
`entrywise kernel gauge or O(delta) full-image reset` dichotomy.  It does not
claim to escape an arbitrarily broad zero-holonomy recurrence category.

I found no theorem-level error.  Two wording edits would make the proof fully
literal:

1. in the gauge sentence, replace “after the first tangent reset” by “once
   the syndetic hypothesis supplies a completed tangent reset in the final
   length-`L` window”; one isolated old reset does not control residuals
   accumulated arbitrarily long afterward;
2. in the selector lower-bound proof, replace “its two preimages” by “the
   two distinct input coordinates selected by that output pair.”

Neither changes a bound or conclusion.

## 1. Finite-semigroup hybrid bound

Let a shortest representative of the element represented by `w` be

```math
v=a_m\cdots a_1,
\qquad m\le L.
```

Because both families are exact actions of the same semigroup,

```math
F_w=F_v,
\qquad G_w=G_v.
```

Introduce the `m+1` hybrid products obtained by replacing the factors of
`F_v` by the corresponding factors of `G_v` one at a time.  Two adjacent
hybrids apply `F_(a_i)` and `G_(a_i)` to the same intermediate point.  Their
distance there is at most `epsilon`, and every common output-side suffix is
nonexpansive.  Thus each hybrid step costs at most `epsilon`, giving

```math
sup_x d(F_wx,G_wx)\le m\epsilon\le L\epsilon.
```

No commutativity, cancellation, invertibility, or compactness is used.  The
hypothesis that the two fixed families satisfy the *same exact semigroup
relations* is essential; one-step closeness alone does not imply it.

For the clamp example, both maps represent the generator of the monoid
`{1,p}` with relation `p^2=p`.  The largest shortest normal-form length is
one, so the theorem gives exactly the observed one-step error at every
positive depth.

## 2. Clamp orientation and all constants

The repository uses the column-output convention

```math
(F_Su)_b=\max_a(u_a+S_(ab)),
\qquad u=(0,z),
\qquad z'=(F_Su)_2-(F_Su)_1.
```

For

```math
S_0=\begin{pmatrix}0&0\\-1&0\end{pmatrix}
```

the two outputs are `max(0,z-1)` and `max(0,z)`, hence

```math
z'=\max(0,z)-\max(0,z-1)=clip(z,0,1).
```

For

```math
S_\delta=\begin{pmatrix}0&\delta\\-1&0\end{pmatrix}
```

they are `max(0,z-1)` and `max(delta,z)`, hence

```math
z'=clip(z,\delta,1)
```

for `0<delta<1`.  Both clips are idempotent and nonexpansive.  In the
half-oscillation convention,

```math
d_H(z,z')={1\over2}|z-z'|,
```

so

```math
sup_z d_H(P_0z,P_\delta z)=\delta/2.
```

The projective image diameters are respectively `1/2` and
`(1-delta)/2`, and remain so under every positive power.  Thus neither is an
`O(delta)` full-image reset as `delta -> 0`.

The kernel difference is

```math
E=S_\delta-S_0=
\begin{pmatrix}0&\delta\\0&0\end{pmatrix},
```

and its alternating rectangle is

```math
E_(11)+E_(22)-E_(12)-E_(21)=-\delta.
```

It is therefore not row-plus-column separable, which is precisely the
entrywise endpoint-gauge obstruction used by the earlier theorem.

The drifting comparator is also oriented correctly:

```math
\widehat S_\delta=
\begin{pmatrix}0&0\\-1+\delta&\delta\end{pmatrix}
```

induces `clip(z+delta,0,1)`.  Starting from zero its coordinate displacement
is `min(t delta,1)`, and its Hilbert displacement is half of that.  Its
kernel difference from `S_0` is the row potential
`((0,0),(delta,delta))`; the drift is an adjacent-interface mismatch, not a
rectangle defect.

This confirms the exact scope of the counterexample.  The paired clamp
orbit is fixed after the first application, so finite-semigroup absorption
can equally be regarded as a finite certificate for zero-increment paired
recurrence.  The integrated theorem, README, axioms, and Example 30 all state
this qualification accurately.

## 3. Selector-reset upper bound

For

```math
e_t=P_te_(t-1)+\eta_t,
```

unrolling gives

```math
e_T=\sum_(s=1)^T P_TP_(T-1)\cdots P_(s+1)\eta_s.
```

Take a reset factor occupying times `a,...,b` in the final length-`L`
window.  Every disturbance with `s<a` traverses the zero projective product
`P_b...P_a` and is killed.  A disturbance inserted at time `a` is added
*after* `P_a`, so it need not be killed.  Exactly the terms from `a` through
`T` can survive, and

```math
T-a+1\le L.
```

Selector products are nonexpansive, so their total Hilbert norm is at most
`L epsilon`.  This verifies both the orientation and the constant.  The
factorial-language assumption ensures that the final window is itself an
allowed word.  Before length `L`, the direct bound is `T epsilon`.

## 4. Reset-free adversarial lower construction

For a reset-free word define

```math
A_s=P_T\cdots P_(s+1),
\qquad 1\le s\le T,
```

with `A_T` the identity.  Since `r>=2` and every nonempty suffix is a
reset-free contiguous factor, each `A_s` is a nonconstant selector.  Choose
an ordered pair of *final output coordinates* `(j_s,k_s)` which `A_s` maps
to two distinct input coordinates `(alpha_s,beta_s)`.  Among the
`r(r-1)` possible ordered output pairs, one pair `(j,k)` occurs at least

```math
ceil(T/[r(r-1)])
```

times.  (The theorem uses the weaker floor.)

For each occurrence put

```math
\eta_s(alpha_s)=+\epsilon,
\qquad
\eta_s(beta_s)=-\epsilon,
```

and set every other disturbance to zero.  Then

```math
||\eta_s||_H=\epsilon,
\qquad
(A_s\eta_s)_j-(A_s\eta_s)_k=2\epsilon.
```

Every selected term has the same sign on the same final coordinate pair, so
their sum has Hilbert norm at least one `epsilon` per occurrence.  Hence

```math
||e_T||_H
\ge ceil(T/[r(r-1)])\epsilon
\ge floor(T/[r(r-1)])\epsilon.
```

The advertised converse follows exactly: if every residual sequence obeys
`||e_T||_H <= C epsilon`, a reset-free word must satisfy

```math
floor(T/[r(r-1)])\le C,
```

and therefore `T < (C+1)r(r-1)`.

## 5. Endpoint gauges

With an endpoint term, the exact recursion is

```math
e_t=P_te_(t-1)+h_t-P_th_(t-1)+\eta_t.
```

Writing `bar e_t=e_t-h_t` gives

```math
\bar e_t=P_t\bar e_(t-1)+\eta_t,
\qquad \bar e_0=e_0-h_0.
```

Thus before a completed reset one must pay

```math
||P_T\cdots P_1(e_0-h_0)||_H
```

in addition to accumulated residuals, and converting back to `e_T` costs
`||h_T||_H`.  Once the syndetic hypothesis supplies a reset factor in the
final length-`L` window, the initial term and all earlier residuals are
killed, yielding

```math
||e_T||_H\le ||h_T||_H+L\epsilon.
```

This is the content of integrated (16.109).  The phrase “after the first
tangent reset” is harmless in the surrounding syndetic context but is not
literally sufficient without continued bounded reset gaps; the more precise
final-window wording above is preferable.

## 6. Affine-selector cycle criterion

For `A(e)=P_sigma e+b`, direct iteration gives

```math
A^k(e)=P_\sigma^ke+\sum_(t=0)^(k-1)P_\sigma^tb.
```

Coordinate `j` samples `b` along the functional-graph orbit

```math
j,\sigma(j),\sigma^2(j),...
```

Transient contributions are bounded.  On a directed cycle `C`, the linear
growth rate is its mean

```math
\beta_C=|C|^(-1)\sum_(j\in C)b_j.
```

If two cycle means differ, coordinates in their basins differ by

```math
k(\beta_C-\beta_(C'))+O(1),
```

so the Hilbert norm grows at least
`k |beta_C-beta_(C')|/2-O(1)`.  Conversely, when all means equal `beta`, the
cycle equations

```math
p_j-p_(sigma(j))=b_j-\beta
```

are consistent.  They define `p` on every cycle and then recursively on its
incoming trees.  Hence

```math
b=p-P_\sigma p+\beta 1
```

and telescoping yields

```math
A^k(e)=p+P_\sigma^k(e-p)+k\beta 1.
```

The selector term is projectively bounded, proving sufficiency.  All signs,
transport directions, and constants in (16.110)--(16.111) are correct.

I additionally exhaustively checked all selector maps on `r=2,3,4` with
translation coordinates in `{-2,-1,0,1,2}`: 163,475 affine-selector cases
agreed with the common-cycle-mean criterion and its predicted drift.

## 7. Executable verification

Running

```text
python3 extremal_information/experiments/verify_tropical_absorption.py
```

returned

```text
exact clamp/idempotence checks: 1275
reset-free selector lower-bound checks: 1811
```

The checker has the same product orientation as (16.112).  Its `r=3` word
test is sampled rather than exhaustive, so the analytic proof above—not the
sample—is the verification of the general lower bound.

## 8. Final scope judgment

Theorem 16.14 safely establishes three facts under distinct quantifiers:

1. common finite-semigroup factorization is a sufficient coherent
   bounded-normal-form mechanism;
2. endpoint gauge plus syndetic tangent reset is quantitatively complete for
   fresh adversarial residuals on a fixed factorial selector language;
3. a recurrent affine selector has bounded projective drift exactly when its
   translation is a selector-twisted coboundary.

It does **not** classify arbitrary coherent all-finite max-plus families,
switching/tied selector cells, near conjugacies, or general infinite compact
semigroups.  Most importantly, it does not overstate the idempotent clamp:
the example only falsifies the kernel-gauge/full-image-reset formulation and
is compatible with a broad zero-increment recurrent interpretation.

