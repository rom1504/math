# Archive verification: contrarian nonconvergence

## Verdict

The normalization and the conditional implications in the frozen report are
correct.  None of the three proposals proves, or nearly proves, genuine
nonconvergence.  In every case the required quantifiers are

```math
\exists A_r^-\ \forall x\quad Q(A_r^-)\le c_- (n_r^-)^{3/2},
\qquad
\forall A_r^+\ \exists x\quad Q(A_r^+)\ge c_+(n_r^+)^{3/2},
```

with `c_-<c_+` uniformly on two infinite subsequences.  The arithmetic and
tensor identities in the report instead give selected-family statements of
the form `exists A exists x`.  They do not support the all-signings high
theorem, and absence of one displayed state is not an upper bound over all
spins.

The elementary implication is sound: `0.478=239/500` and
`0.479=479/1000`, so the two clauses of either `L_FF` or `L_KP` would give

```math
\liminf_n M_n/n^{3/2}\le0.478<0.479\le
\limsup_n M_n/n^{3/2}.
```

The tensor constants likewise differ by exactly `delta_H/8` in the original
half-quadratic normalization.

## Independent Paley calculation

Let `q=5^m`, let `B(a,b)=chi_q(a-b)`, and put
`x_a=phi(Tr(ta))`, where `sum_(u in F_5) phi(u)=1`.  Then
`sum_a x_a=q/5`.  The additive Fourier support consists of the DC mode and
the four modes indexed by `k t`, `k in F_5^*`.  If `m` is even, then
`chi_q(k)=chi_5(k)^m=1`, so those four modes have one common Paley eigenvalue
of magnitude `sqrt(q)`.  Parseval therefore gives exactly

```math
|x^T Bx|=\sqrt q\left(q-{(q/5)^2\over q}\right)
={24\over25}q^{3/2}.
```

For the bordered conference matrix, choose the infinity spin to align with
the core.  In matrix normalization the border adds `2q/5`; hence, in the
problem's normalization,

```math
Q(C_q)\ge {12\over25}q^{3/2}+{q\over5}.
```

Thus the stated `0.48+o(1)` selected-state calculation is correct.  For odd
`m`, `chi_q|_(F_5^*)=chi_5`, so this particular four-mode state splits across
the two spectral halves.  That proves only the loss of this witness; it does
not bound the maximum over all Boolean vectors.

This calculation is not a new low-side theorem.  Ledger Section 10.120 and
`artifacts/constructive_family_phase2_report.md` prove the stronger result
that every square-field Paley conference `PC(r^2+1)` has an exact Boolean
eigenvector and cap `r(r^2+1)/2`, i.e. normalized cap tending to `1/2`.
Here `5^{2r}=(5^r)^2`.  The new trace-fibre identity is a correct weaker
selected-state formula, in the wrong direction for a low construction.

## Candidate 1: fixed-characteristic Paley alternation — **C**

Mechanism stripped of terminology:

* proved arithmetic: for one selected even-tower matrix there exists a spin
  of normalized energy `0.48+o(1)`;
* missing low theorem: the selected odd Paley matrix has **all** Boolean
  energies at most `0.478 n^{3/2}`;
* missing high theorem: **every** signing at the even-tower order has some
  Boolean energy at least `0.479 n^{3/2}`.

The latter two bullets are the separated-subsequence obligation itself; the
finite-field calculation supplies neither.  This is an implication-level
collision with ACTIVE_STATE item 8, ledger Sections 1.12, 10.18, 10.120, and
10.129.4, and `artifacts/prime_paley_cosine_saturation.md`,
`artifacts/paley_resonance_gadget.md`, and
`artifacts/bent_parity_nonconvergence_audit.md`.  Prime-Paley resonance does
not formally falsify the fixed-characteristic odd tower, so the proper grade
is `C`, not `D`: the specific odd-Paley upper bound remains an unproved,
falsifiable family question.  What is not new is the asserted route from a
selected algebraic parity effect to a universal order obstruction.

The archive's padding theorem does not itself kill these towers: adjacent
epochs differ by a factor asymptotic to `5`, not by `o(n)`.  Conversely,
ordinary restriction of the larger tower loses a fixed leading factor and
does not prove the proposed transfer.  Dense Paley restrictions give only
the all-order upper bound `1/2+o(1)`; because `0.479<1/2`, they do not
contradict `L_FF`.  The decisive failure is lack of any implication to the
universal high clause (and, separately, lack of the odd-family cap theorem).

Surviving exact issue: determine the Boolean cap of `PC(5^{2r+1}+1)`; this is
a legitimate selected-family problem but cannot establish nonconvergence
without an independent all-signings theorem.  Minimum revision: split off
and prove the odd-Paley upper theorem before proposing any high-side
architecture.

## Candidate 2: Arf/Krawtchouk sector cancellation — **D**

The displayed Krawtchouk identity and the sufficient condition `(KP)` are
correct.  If `(KP)` holds and `P` is nonpositive throughout the central
interval, the coset sum is positive for every translate, so every coset has
a word outside that interval.  This really would be an all-coset theorem and
does not reconstruct a complete coset histogram.

The advertised mod-eight mechanism, however, is rigorously erased before
that certificate is reached.  For even `n`,
`q_n(z)=sum_(i<j)z_i z_j` has polar form `I+J` and the stated Gauss sum is
correct.  Adding an edge word contributes the linear form

```math
z\longmapsto \sum_i \deg_a(i)z_i \pmod2.
```

Even the word consisting of one edge changes the refinement by
`z_i+z_j`, which flips the Arf sign.  Thus arbitrary cosets realize both
signs within either order sector.  Correspondingly, `(KP)` replaces every
coset phase by its absolute value.  That is exactly what makes it uniform,
and exactly what deletes the proposed Arf-sign advantage.

After removing the Arf story, the remainder is the archived radial signed
dual-certificate obligation: ledger Section 3.15 and
`artifacts/eulerian_free_energy_identity.md` show that fixed degree misses
planted extreme cosets, while a resonance-sensitive growing signed hierarchy
recovers the complete energy histogram.  The absolute-dominance variant is
an even stronger, separately-paid phase condition; no archived theorem says
it is feasible at `0.479`.  It is a generic Delsarte certificate search, not
a sector-separated mechanism.  The low Paley clause is also wholly unproved.

Decisive reason for `D`: the only claimed source of separation, the Arf sign,
is not uniform over cosets and is explicitly discarded by `(KP)`.  Surviving
issue: a phase-blind radial LP could still be studied as a generic universal
lower-bound method, but it is no longer this nonconvergence architecture.

## Candidate 3: odd/even symmetric-Hadamard tensor depth — **C**

The exact algebra is correct.  A symmetric Hadamard matrix of order `12`
cannot have a Boolean extremal eigenvector because `sqrt(12)` is irrational.
At even depth,

```math
(G\otimes G)\operatorname{vec}(G)=d\operatorname{vec}(G)
```

gives a Boolean spectral eigenvector.  In fact the report's diagonal loss can
be sharpened: a symmetric order-12 Hadamard has trace zero, so
`tr(H^{\otimes k})=0` and `x^T S_k x=x^T K_k x` exactly for every Boolean
`x`.

But the seed gap does not imply the odd-depth gap for arbitrary entangled
Boolean tensors.  Ledger Section 3.8 gives an explicit failure of quadratic
tensor submultiplicativity, and
`artifacts/switching_natural_product_obstruction.md` isolates the uncontrolled
entangled tensor channel.  Conversely,
`artifacts/fixed_fiber_boolean_channel_test.md` already proves the displayed
two-level vectorization/persistence phenomenon: a missing Boolean eigenvector
can delay a channel by one level but does not remove it.  Those theorems do
not formally disprove the special odd-depth inequality, so `D` would be too
strong.

The high clause again reverses `exists A exists x` for the selected even
tensor into `forall A exists x` over all order-`12^{2r}` signings.  No tensor
or Hadamard theorem links arbitrary signings to that selected algebra.  The
high constant is below `1/2`, so all-order Paley upper propagation does not
contradict it; it simply supplies no support for it.  Accordingly the full
proposal is `C`: the even-depth fact is archived, the odd-depth cap is a new
unproved entangled-norm question, and the universal high theorem is the
original nonconvergence obligation reinserted as a clause.

Surviving issue: the odd-depth Boolean norm of one fixed nonsquare symmetric
Hadamard is a crisp stand-alone problem.  It is not yet a useful low-side
construction because no uniform upper gap has been proved.  Minimum revision:
first prove a fixed odd-depth gap against all entangled spins; then present a
separate, genuinely universal high-side mechanism.

## Cross-audit and useful mathematics

* Restriction/padding continuity rigorously eliminates finite residue and Witt
  stories on dense order sets, but not the factor-`5` or factor-`12` towers.
* All proposed high constants are strictly below `1/2`; hence dense
  conference upper propagation does not directly falsify them.  It does rule
  out treating selected high-energy conference/Hadamard matrices as universal
  lower bounds.
* The archive's symplectic/Witt audit is stronger against residue-only
  explanations: Witt class constrains a selected Cayley matrix, while
  arbitrary edge signings have no such class and the augmented cut code has
  only bounded divisibility.
* No candidate supplies a proved new low-side construction.  Candidate 1's
  valid trace identity is high-energy and weaker than the archived exact
  square-field theorem; Candidate 2 reuses the same unproved odd-Paley upper;
  Candidate 3 gives a new test family but no all-spin upper bound.  The only
  genuinely reusable standalone item is the correct trace-fibre identity,
  useful as a finite-field spectral diagnostic rather than as evidence for
  separated values of `M_n`.

