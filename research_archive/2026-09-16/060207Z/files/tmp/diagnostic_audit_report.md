# Diagnostic 1: adversarial audit of selected priors and the common active face

Status: independent audit of the committed ledger through the consolidation
checkpoint.  This report does not alter any tracked file.

Notation used below is as in (10.792)--(10.795).  For an exact order-`n`
minimizer `A`, put `q_n=Q(A)`, let `p=m/n` and
`p_2=m(m-1)/(n(n-1))`, and for an oriented cut `d=(tau,x)` put

```math
\Delta(d)=q_n-\tau x^T A x,qquad R(d)=x^T A^2x,
```

```math
\widehat\ell(S,d)=\delta_S(d)-p_2\Delta(d)-B_{n,m},
\qquad B_{n,m}=(p^{3/2}-p_2)q_n.
```

The favorable event is `F_t={widehat ell<=t}` and
`u(d)=U_m(F_t^d)`.  The recurrence scales are

```math
H=n^{3/4-c},\qquad T_n=n^{3/2-c},\qquad
R_*=n^{9/4-c},\qquad 0<c<1/4.
```

## Executive conclusion

The selected-prior target is, at the scales used in the recurrence, exactly
the bare fixed-cut tail (10.795) written in averaged form.  It is not a weaker
lemma.  The reference-free information formulation has the same collapse to
one deterministic column.  Captured parent deficit is not needed for
(10.795), because the deficit is already part of `widehat ell`; only favorable
mass and captured row are needed.

The common active-face theorem (10.1340)--(10.1343) is a genuine theorem from
exact minimality and is not circular.  However, it controls polynomial-scale
first moments of selected edge blocks, while the desired favorable collision
has mass `exp(-Theta(H))`.  Its actual completion condition is not weaker than
the bare tail:

- if the prior is unrestricted, favorable mass plus captured row rounds back
  to one bare-tail cut with only constant losses;
- if the prior is required to be the particular active-face law and the event
  must also be an escape event, that support/intersection requirement is an
  *additional* condition not implied by a bare-tail cut.

Thus the positive active-face theorem has removed a witness-multiplicity
defect inside one proposed implementation, but it has not removed an
obligation in the original bare-tail theorem.  The route currently lacks a
mechanism with resolution at the entropy scale `exp(-Theta(H))`.

### Crisp bridge verdicts

| proposed bridge | exact relation to (10.795) | reverse implication / countermodel | leverage supplied by (10.1340) |
|:--|:--|:--|:--|
| unrestricted selected prior: saved `F_t` mass plus captured row | Equivalent up to fixed constants (Theorem A). | A bare-tail cut gives the point-mass prior exactly. | None is needed; this is a change of language. |
| saved bare favourability alone, with no active-law or escape restriction | Equivalent to the row-free statement `max_d u(d)>=exp(-O(H))`; strictly omits the row clause. | The maximum is at least every prior average, and a maximizing cut gives a point-mass prior. | No theorem proves this mass. |
| saved favourability in `F_t cap escape` under an active-face law | Implies the row-free bare fibre, but adds support and escape compatibility. | A bare-tail cut need not lie in that law or escape.  An abstract two-state model can put a perfect active/escape certificate on an unfavorable state and the bare-tail fibre on an unsupported state.  The finite `A_9` instance supplies an actual pointwise escape/favourability separation, though not an asymptotic counterexample. | No entropy-scale leverage: the certificate is stable under arbitrary `e^{-Theta(H)}` perturbations. |
| captured project row on the same saved favourable event | Together with the preceding saved mass, equivalent up to constants to the full bare-tail cut. | Conversely the bare-tail cut gives a point-mass law with the same captured row; only the extra active/escape support need not follow. | None: (10.1340) controls first edge-sign marginals, not row conditioned on a rare event. |

The abstract two-state countermodel in the third row is deliberately only a
logical countermodel to deduction from the stated active-face conclusions:
take `mu=delta_(d_0)`, give `d_0` zero deficit and signs making every named
block escape, declare its favorable fibre empty, and give an unsupported
state `d_1` a full favorable fibre and low row.  It is not asserted to be a
global exact-minimizer construction.  The point is that any reverse theorem
must use additional signing/minimizer structure absent from (10.1340)--
(10.1343); the ledger contains no such theorem.

The strictness claim for the row-free statement also has a simple abstract
column-table witness: give one cut `u(d)=1` and row `n^(5/2)`, and give every
cut of row at most `R_*` an empty fibre.  Saved favorability holds while
(10.795) fails.  The scalable hub-avoidance wall (10.1373) expresses the same
incidence phenomenon with the project exponents, but neither example is an
exact-minimizer sequence.  Thus the row-free clause is logically weaker on
general incidence data, while no verified theorem says it is easier on the
required class of exact minimizers.

## 1. Exact deterministic-extraction theorem

Fix an exact minimizer `A`, a selector slice `Omega`, a threshold `t`, and a
set `D` of oriented global cuts.  Write

```math
u(d)=U_m\{S:(S,d)\in F_t\},\qquad R(d)=x^T A^2x.
```

For a prior `nu` on `D`, put

```math
Z_\nu=\sum_d\nu(d)u(d),\qquad
\overline R_\nu=\frac{\sum_d\nu(d)u(d)R(d)}{Z_\nu}.
```

### Theorem A (one-cost equivalence, with explicit constants)

If `Z_nu>0`, then for every `a>1` there is a cut `d` such that

```math
u(d)\ge (1-a^{-1})Z_\nu,
\qquad R(d)\le a\overline R_\nu.
```

Conversely, if one cut has `u(d)>=z` and `R(d)<=r`, its point-mass prior has
`Z_nu>=z` and `Rbar_nu<=r`.

Proof.  Under the output marginal conditioned on favorability,

```math
\pi(d)=\frac{\nu(d)u(d)}{Z_\nu},
```

Markov gives `pi{R<=a Rbar}>=1-1/a`.  If every cut in that set had
`u(d)<(1-1/a)Z_nu`, its `pi`-mass would be strictly less than

```math
(1-a^{-1})\sum_{R(d)\le a\overline R_\nu}\nu(d)
\le 1-a^{-1},
```

a contradiction.  The converse is immediate.

In particular, with `a=2`,

```math
Z_\nu\ge e^{-CH},\quad \overline R_\nu\le C_RR_*
\quad\Longrightarrow\quad
u(d)\ge \tfrac12e^{-CH},\quad R(d)\le2C_RR_*.
```

The multiplicative `1/2` changes `-log u` by only `log 2`, and the row factor
two is harmless in the `O(R_*)` hypothesis.  This is exactly (10.795) at the
scales used by (10.1023)--(10.1024).

### Theorem B (two costs)

If additionally

```math
\overline\Delta_\nu
=Z_\nu^{-1}\sum_d\nu(d)u(d)\Delta(d),
```

then for any `a,b>1` with `theta=1-1/a-1/b>0`, some cut has

```math
u(d)\ge\theta Z_\nu,
\quad R(d)\le a\overline R_\nu,
\quad\Delta(d)\le b\overline\Delta_\nu.
```

This is the same union-bound/Markov proof.  It confirms (10.1353), but the
third coordinate is an optional strengthening: (10.795) itself imposes no
separate deficit bound.  Its event is

```math
\widehat\ell(S,d)=\delta_S(d)-p_2\Delta(d)-B\le t,
```

so parent deficit has already been retained inside the favorable fibre.
The row-only version of Theorem A is sufficient.

### Corollary (selected priors are not a strict reduction)

At `z=exp(-O(H))` and `r=O(R_*)`, existence of a successful unrestricted
selected prior and existence of a successful deterministic cut are
equivalent up to fixed constants.  The variational equality (10.1356) is the
penalized exact version of the same fact: the optimum is a constant output.

Likewise, for the reference-free formulation, (10.1355) gives
`K(P)>=E h_D`, where `h_d=log(1/u(d))`.  If `K(P)<=K_0` and
`E R(D)<=R_0`, Markov on `h_D` and `R(D)` gives a cut with, for example,

```math
u(d)\ge e^{-3K_0},\qquad R(d)\le3R_0.
```

Conversely, the uniform law on one cut's favorable fibre has
`K(P)=h_d` and mean row `R(d)`.  Thus the low-information collision language
also makes no asymptotic reduction.

## 2. What the active-face support restriction changes

Let `A_{eta,r}` denote laws satisfying the exact active-face certificate

```math
\Phi(\mu):=\mathbb E_\mu\Delta
+4h_r((\mathbb E_\mu s_e)_e)\le\eta.
```

There are two logically different readings of the proposed successor.

### Unrestricted-prior reading

If one uses (10.1340) only as motivation and asks for *some* prior with saved
bare mass and captured row, Theorem A shows exact equivalence with (10.795).
The words "selected prior" remove bookkeeping but no mathematical
obligation.

### Particular-active-law reading

The aggregate proposal (10.1361) asks for the particular law delivered by
the active-face theorem to satisfy

```math
(U_m\otimes\mu)(F_t\cap E)\ge e^{-C'H}
```

and to have captured row `O(R_*)` on that same event, where `E` is canonical
block escape.  This implies (10.795) by Theorem A after replacing `u(d)` by
the smaller fibre `U_m(F_t^d cap E_d)`.  The converse is unavailable:

1. a bare-tail cut need not belong to the support of the minimax law;
2. it need not escape the chosen canonical child block; and
3. it need not be near-ground.  Indeed (10.795) permits arbitrary parent
   deficit because `-p_2 Delta` is retained in `widehat ell`.

The finite `A_9` diagnostic in (10.1358)--(10.1361), where escape and
zero-slack favorability are disjoint for a parent ground, already falsifies
the pointwise implication.  It is not an asymptotic falsifier, but it confirms
that the extra intersection is real.  Therefore the actual common-law
completion implies the bare tail while its converse is unavailable; it is
potentially strictly stronger, and is not a verified weakening of it.

## 3. Entropy-scale blindness of the active-face certificate

There is a precise scale obstruction.  Since `h_r` is the support function
of the fractional block polytope, it is convex and positively homogeneous.
For any state `d`, any active law `mu_0`, and
`mu_alpha=(1-alpha)mu_0+alpha delta_d`,

```math
\Phi(\mu_\alpha)
\le(1-\alpha)\Phi(\mu_0)
+\alpha\{\Delta(d)+4h_r(s(d))\}
\le(1-\alpha)\eta+\alpha\{\Delta(d)+4r\}.
```

For `r=Theta(n^(3/2))`, every deficit is `O(n^(3/2))`.  At
`alpha=exp(-CH)`, the change is `exp(-CH) poly(n)`, negligible compared with
`eta=Theta(n^(5/4))` and with every tolerance used in the common-face
corollary.  More generally, changing a law by total variation `alpha` changes
the certificate by at most `O(alpha n^(3/2))`.

Hence the polynomial-precision moment statement (10.1340) cannot by itself
distinguish arbitrary behavior placed on exactly the rare mass
`exp(-Theta(H))` which the selected-prior lemma needs.  If one relaxes the
certificate to its asymptotic conclusions (`E Delta=o(T_n)` and
`h_r=o(r)`), an arbitrary cut can be inserted at the desired entropy scale
without violating them.  If instead one freezes the exact minimax optimizer
and its conditioned support, one has added the unsupported compatibility
condition described above.  Neither interpretation removes the tail
obligation.

This scale mismatch does not prove that no stronger exact-minimizer theorem
can use the same law.  It proves that the currently stated active-face
certificate and marginal escape conclusion have no resolution at the
required collision scale.

## 4. Audit of the two advertised missing bridges

### Saved bare favourability

Without the words "active-face" and "escape", saved average favorability
`Z>=exp(-O(H))` is exactly the row-free fixed-cut statement
`max_d u(d)>=exp(-O(H))`, since an average is at most its maximum.  This is a
proper relaxation of (10.795) only because it omits row.

The actual proposed bridge is stronger: it asks for saved favorability in the
intersection with canonical escape under a law already selected by a
different minimax problem.  A bare-tail cut does not imply that statement.
No theorem currently shows that this support/intersection constraint makes
the proof easier.

There is also a circular scalar route to avoid.  By (10.796), for every fixed
cut, and hence for every selector-independent prior,

```math
\mathbb E_{U_m\otimes\nu}\widehat\ell(S,D)
=\mathbb E_{U_m}Q(A[S])-p^{3/2}q_n.
```

The right side contains the unknown restriction excess and is independent of
the prior.  More explicitly, an upper bound

```math
\mathbb E_{U_m\otimes\nu}\widehat\ell(S,D)\le CT_n
```

would immediately imply

```math
q_m-p^{3/2}q_n
\le\mathbb E_{U_m}Q(A[S])-p^{3/2}q_n\le CT_n,
```

which is the desired recurrence before any tail extraction is used.  A
concentration argument centered at that mean has the same problem: locating
the threshold within `O(T_n)` of the mean first requires this recurrence-size
upper bound.  This is genuine circularity of the scalar-averaging proof plan,
not merely a technically difficult estimate.  A noncircular successor would
need a one-sided tail mechanism not obtained by first controlling this mean.
The active-face theorem does not change the identity.

### Captured project row on the same event

This is not even a nonvacuous stand-alone bridge when favorable mass is zero,
and it is not an independent weaker bridge once saved favorability is assumed.
Theorem A says precisely that saved mass plus captured row is a single
low-row saved column.  Conversely that column gives the point-mass prior.
Thus this second bridge completes the exact bare-tail equivalence, up to
constant factors.  Marginal row on all escape events would be weaker, but it
would not help: favorability and escape can be disjoint, and conditioning can
raise the row arbitrarily.  What recurrence needs is row captured on the
same favorable mass.

## 5. Recurrence-dependency and circularity audit

The proved chain from (10.795) is logically sound:

1. `Q(A[S])>=q_m` is the definition of the order-`m` optimum, not a
   recurrence assumption.
2. `widehat ell<=t` implies the centered fixed-cut payoff
   `X_d>=G_(n,m)-t`, where `G_(n,m)=q_m-p^(3/2)q_n`.
3. Fixed-slice Hanson--Wright (10.1022), the hypotheses
   `R(d)=O(R_*)`, `u(d)>=exp(-O(H))`, and `t=O(T_n)` bound
   `G_(n,m)=O(T_n)` in (10.1023).
4. The auxiliary inputs `q_n=O(n^(3/2))` and
   `||A||_op=O(n^(3/4))` are previously verified exact-minimizer bounds; they
   do not assume convergence.
5. Uniform fixed-density recurrence, geometric-window summability, and exact
   landing then imply convergence.

The common active-face theorem is also noncircular:

1. fractional attenuation plus simultaneous sign rounding is compared with
   the exact minimum cap `q_n`;
2. minimax gives (10.1340);
3. conditioning and the unconditional bounds `q_m=Theta(m^(3/2))`,
   `q_n=Theta(n^(3/2))` make the canonical block size `Theta(n^(3/2))` and
   yield (10.1342)--(10.1343).

Those order bounds come from the established universal lower and constructive
upper bounds, not from the desired recurrence.  I found no hidden use of the
convergence recurrence in (10.1340)--(10.1343).

The circularity is instead in possible *successor arguments*: controlling the
mean of `widehat ell` invokes exactly the unknown endpoint excess, and proving
saved favorable mass plus captured row is, by Theorem A, already proving the
bare-tail lemma.

## 6. Research judgment

Concrete result of this diagnostic: **Theorem A is an exact reverse
implication showing that the selected-prior/captured-row package is
asymptotically equivalent to the bare-tail statement.**  The common
active-face theorem remains a correct positive structural fact, but its
present conclusions are polynomial-scale marginals blind to the target rare
mass.  Requiring the rare favorable/low-row mass to lie on its canonical
escape support adds a compatibility obligation rather than deleting one.

Accordingly this route does not justify another ordinary wave unless a new
theorem is proposed that has entropy-scale resolution and is stated directly
in a quantity not already equivalent to `u(d)` and `R(d)`.  Examples of what
would count are an exact-minimizer inequality lower-bounding a favorable
low-row collision from a polynomial active-face slack, or a multiscale
transport theorem that reduces row while losing at most `O(H)` log-mass.  No
such inequality is presently in the ledger.

## Explicit limitations

1. This audit proves an equivalence of sufficient conditions and a scale
   limitation of the *stated* active-face certificate.  It does not prove
   that every possible strengthening of the active-face method must fail.
2. The two abstract countermodels separate logical interfaces only.  They are
   not exact-minimizer signings and do not constitute asymptotic route
   falsifiers.  The `A_9` separation is an exact finite mechanism falsifier,
   not a fixed-density asymptotic construction.
3. No claim here improves the rigorous interval, proves convergence, or
   proves nonconvergence.  The new concrete theorem is the explicit reverse
   extraction/equivalence and the correction that captured deficit is not a
   clause of (10.795).
4. `tmp/diagnostic_audit_extraction_check.py` passed 200,000 arbitrary finite
   table instances.  That computation is only a sanity check; the proof of
   Theorems A and B is the displayed Markov/union-bound argument.
