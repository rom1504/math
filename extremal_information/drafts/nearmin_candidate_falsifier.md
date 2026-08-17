# Adversarial audit of P1, P2, and radial replacement

Date: 2026-08-17.

Status: independent counterexample-builder report.  The candidate conclusions
below were frozen after reading only the campaign prompt, the contextual
targets, the blind structural audit, and the deterministic inequalities.  No
unrelated theory archive was consulted.  The finite computations are exact;
no scalable counterexample to P1 or P2 is claimed.

## 1. Bottom line

I did **not** find a scalable certified near-minimizer counterexample to the
edit-thick P1 or P2 statements.  Sparse edits, planted blocks, and local
surgery cannot provide one by themselves: the very same edit set supplies the
allowed peeled core.  This is a rigorous barrier to this counterexample
method, not evidence that P1 or P2 is true.

I did find a sharp finite falsifier to the naive radial principle.  There is
an order-8 signing with

```math
Q(A)=12=M_8+2
```

whose edge-Hamming distance from **every** exact order-8 minimizer is three
and whose response distance

```math
d_\square(A,\mathcal E_8)
:=\min_{B:Q(B)=M_8}\max_x|H_A(x)-H_B(x)|
```

is six.  Thus neither

```math
d_H(A,\mathcal E_n)\le {Q(A)-M_n\over2}
```

nor

```math
d_\square(A,\mathcal E_n)\le Q(A)-M_n
```

is true even at small order.  Any uniform linear replacement theorem needs
constant at least `3/2` in Hamming distance, or at least `3` in response
distance, under these normalizations.  This does **not** disprove an
unspecified `O(Q-M_n)` theorem, and the data do not establish growth of the
constant with `n`.

The main adversarial conclusion is consequently a narrowing:

> A scalable attack on P1 or P2 must either find bad **exact** minimizers, or
> prove that near-minimizers can lie response-far from every good exact core.
> Ordinary sparse surgery, block implants, and geodesic planting do neither.

## 2. The peeling barrier for sparse-edit counterexamples

Let `s_n=n^(3/2)`, let `\mathcal E_n` be the exact minimizer set, and let
`A` differ from `A_0 in \mathcal E_n` on `t` unordered edges.  Then

```math
Q(A)-M_n\le2t,
\qquad
d_\square(A,A_0)\le2t.                            \tag{2.1}
```

Therefore any thick property of the form

```math
\exists B\in\mathcal G_n:\quad d_\square(A,B)\le h
```

is automatically satisfied with `B=A_0` and `h=2t`, provided the exact core
belongs to `\mathcal G_n`.  In particular, an `o(s_n)` clique implant,
arbitrary surgery on `o(s_n)` edges, and the geodesically planted face are
all absorbed by precisely the residual that P1 and P2 were designed to
allow.

This remains true at the level of the P1 selector quadratic.  In the notation
of P1, fix a port frame and endpoint, and write

```math
y=Za^\epsilon.
```

Because the Fourier sum evaluates the Boolean selector rowwise,
`y_i in {+-1}`.  If `B` and `B'` differ on `t` edges and use the same
normalization `r n`, then

```math
|\Delta_\tau(B,W)-\Delta_\tau(B',W)|
\le {4t\over rn}.                                 \tag{2.2}
```

Indeed each edited unordered edge contributes at most four to
`|y^T(B-B')y|`.  With `r=Theta(sqrt(n))`, an `o(s_n)` edit changes every fixed
selector defect by `o(1)`.  Conversely, creating a fixed defect from a good
core by this route costs `Omega(s_n)` edits and therefore has fixed normalized
optimality excess under the only available certificate.

Equation (2.2) explains why none of the requested sparse-edit constructions
produces the desired discontinuity.  It also shows what would be needed:
either an exact minimizer with a fixed P1 defect, or a near-minimizer whose
best core remains defective after the existential peeling.

## 3. Exact finite radial audit

Two reproducible programs were added:

```text
extremal_information/experiments/nearmin_radial_distance_small.py
extremal_information/experiments/nearmin_radial_distance_n8.py
```

The first enumerates every root-gauged signing through order 7, expands the
complete minimizer set by switching, and computes distance to it by an exact
multi-source Hamming-cube BFS.  The maximum ratios

```math
{d_H(A,\mathcal E_n)\over (Q(A)-M_n)/2}
```

seen at orders `4,5,6,7` are respectively

```text
1, 1, 3/2, 2.
```

At order 7, some cap-`M_7+2` signings require two edge edits rather than one.

The second program uses the authoritative two-orbit classification of all
4,200 root-gauged order-8 minimizers, closes it under all 128 switchings, and
compares 89 distinct saved cap-12 signings against all 537,600 labelled exact
minimizers.  The exact distance distribution is

```text
distance 1: 41 signings
distance 2: 40 signings
distance 3:  8 signings.
```

One distance-three witness has SHA-256

```text
20bb81c068b556af607a07e73d10497887d9599c099b4b90d8bc050c493d8e11
```

and upper-triangle encoding

```text
+-+++++/--++-+/+-+--/+-++/-+-/-+/-
```

The response-distance certificate is also exact.  If two signings differ on
an edge set `F`, then orthogonality of the Boolean characters `x_i x_j`
gives

```math
\mathbb E_x\left(\sum_{e\in F}c_ex_e\right)^2=|F|,
```

so their response distance is at least `2sqrt(|F|)`.  Hence response distance
below six can occur only at Hamming radius at most four.  The program checks
every exact minimizer in that radius explicitly; the minimum is six.  All
farther minimizers have integer sparse cap at least three, hence response
distance at least six.

This is a genuine finite falsifier, not a scalable one.  An `O(Q-M_n)` claim
permits constants at least those above.  Replicating the order-8 gadget in
blocks does not certify near-minimality for the complete-graph signing problem,
because no corresponding formula for `M_{8k}` or for the cross-block edges is
known.

## 4. Attack on P1: collective selector coercivity

### What survives the attack

P1's thick form survives every sparse/local construction I could certify.
The peeling barrier (2.1) removes the surgery, and (2.2) prevents an
`o(s_n)` surgery from leaving a fixed selector defect after the good core is
restored.  Tensoring or taking block sums does not help: without a theorem
relating the resulting signing to `M_N`, it ceases to be a certified
near-minimizer; with only the edit-Lipschitz certificate, it again comes with
its original exact core.

### A formulation weakness

P1 is not yet an intrinsic property of `A`: it depends on a declared port
language `mathscr W`.  Non-emptiness alone does not fix this.  A benign
language can avoid a bad frame, while an adversarial language can include
irrelevant low-quality frames.  The optional all-frame condition (2.5) is
intrinsic, but the actual SML (2.13) asks for a nonvacuous language generated
by the intended composition without defining that language independently of
the desired core.

This is not a counterexample, but it blocks a clean falsification theorem.
Before P1 is promoted, the port-generation rule should be frozen without
reference to the candidate core.  Otherwise `near-minimality => P1` can be
made weak or strong by changing the query family.

### Sharpest remaining falsifier

The exact finite falsifier is a supplied frame `W` and endpoint `epsilon`
whose joint defect exceeds the budget after optimizing over every allowed
core in the specified response ball.  No current finite dataset supplies the
last universal-over-cores certificate.  The order-8 radial witness shows that
the response ball cannot be assumed to contain an exact core at unit slope,
but its P1 defect has not been certified.

Verdict: **not falsified; sparse-edit counterexample class exhausted by an
exact peeling argument; port language must be made intrinsic.**

## 5. Attack on P2: contracting fibre plus sparse recurrent memory

P2 is still broader and therefore harder to falsify.  The same edit residual
absorbs all local surgery.  A few persistent planted bits may simply be stored
as recurrent cohomology, exactly as the definition permits.

There is also a semantic circularity risk.  For any finite declared query
language, one may always use singleton fibres indexed by its exact contextual
Myhill--Nerode quotient.  The centred contraction is then vacuous; all the
difficulty has moved into proving that the quotient and recurrent rank have
`o(n)` description.  Thus the P2 SML is not independently checkable until a
carrier map, or a restricted carrier construction class, is specified without
first enumerating the response language.  As written, a lower bound against
P2 is essentially the requested contextual information lower bound itself.

I tested the natural amplification ideas conceptually:

* many microscopic edit bits can be decoded by full Boolean pinning, but full
  pinning reconstructs the entire landscape and is explicitly not the
  intended composition language;
* repeating fixed gadgets in `n/b` blocks creates total response scale
  `n sqrt(b)`, which is `o(n^(3/2))` whenever `b=o(n)` and can therefore be
  discarded at the target accuracy;
* using blocks of linear size preserves target-scale response but permits only
  `O(1)` independent blocks;
* random sparse edits can have Hamming size much larger than their individual
  edge scale, but their same-spin response is of the same order as the best
  available cap-increase certificate, so the residual pays for them.

These calculations explain why a scalable P2 falsifier needs a genuinely
global response packing, not a direct sum of local memories.  The current
multi-selector lower bound is growing but sublinear and hence does not yet
contradict (3.4).

Verdict: **not falsified; current semantic form is too broad for a structural
near-minimizer test.  Freeze an explicit carrier class or prove a linear
physical selector packing.**

## 6. Radial replacement: exact boundary of the result

Three different statements should not be conflated.

1. **Hamming radial repair**

   ```math
   d_H(A,\mathcal E_n)\le C[Q(A)-M_n].             \tag{6.1}
   ```

   The finite audit forces `C>=3/2`; no scalable verdict is known.

2. **Response radial repair**

   ```math
   d_\square(A,\mathcal E_n)\le C[Q(A)-M_n].       \tag{6.2}
   ```

   The finite audit forces `C>=3`; no scalable verdict is known.  This is the
   version actually relevant to P1/P2.

3. **Vanishing normalized repair**

   ```math
   Q(A_n)-M_n=o(s_n)
   \Longrightarrow d_\square(A_n,\mathcal E_n)=o(s_n). \tag{6.3}
   ```

   This is strictly weaker than a dimension-free linear inequality and is the
   minimal useful radial statement for thick response cores.  None of the
   current constructions falsifies it.

Generic max-of-affine convexity does not prove any of these: a 1-Lipschitz
integer function on a cube may have broad near-minimal plateaux.  The cut-code
geometry must supply a genuine error-bound theorem.  Conversely, a large
Hamming distance alone need not refute (6.3), because character orthogonality
only converts distance `r` into response separation `Omega(sqrt(r))`.

## 7. Director-facing recommendation

The strongest falsifiable next statement is (6.3), not the coefficient-one
or Hamming versions.  It is both weaker than identifying a structured P1/P2
core and strong enough to remove every `o(s_n)` radial surgery.  A decisive
experiment should compute or bound

```math
\inf_{B\in\mathcal E_n}d_\square(A,B)
```

for certified one-step near-minimizers at the largest orders where the full
minimizer set is classified.  At present that means order 8, and the exact
answer above is six.  Larger-order saved representatives cannot certify the
infimum because the complete exact minimizer set is unknown.

For theory selection:

* P1 remains the more concrete live structure, but only after freezing an
  intrinsic port language.
* P2 should not be selected as a near-minimality lemma in its current semantic
  form; it is closer to the desired conclusion than to a checkable matrix
  property.
* Radial response recovery (6.3) is a legitimate alternative structural
  lemma.  The order-8 witness warns that any quantitative proof needs a real
  stability argument and cannot use a unit-slope descent heuristic.

Classification: **one exact finite radial falsifier; one rigorous
counterexample-method ceiling; no scalable P1/P2 counterexample.**
