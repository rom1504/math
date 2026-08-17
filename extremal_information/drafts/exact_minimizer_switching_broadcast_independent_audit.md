# Independent audit: exact-minimizer switching broadcast

**Overall verdict: CORE THEOREMS PASS; CENTRAL PHYSICAL-EXPOSURE CLAIM FAILS
AS WRITTEN AND REQUIRES A SCOPE REPAIR.**

SB.1--SB.7 are correct after one minor asymptotic quantifier repair in SB.2.
The tail hypothesis SB.8 correctly invokes Theorem 21.8 to produce a
linear-rate packing of *boundary response profiles*.  It does not, by
itself, produce a bank of all-spins-free bounded-cap parent-cap gaps.  The
source draft conflates those two contextual languages in Sections 5--6.

One can fill the bridge to a complete order-`2n` exact-sign parent with cap
`O(n^(3/2))`, and the separated conditional response profiles survive that
fill.  Optimizing the boundary spins, however, can erase their projective
profile separation.  A further low-cap boundary-selector/anti-pin theorem is
needed for the claimed physical exposure.

No canonical theorem, axiom, or frontier file was edited in this audit.

## 1. Frozen source

```text
extremal_information/drafts/exact_minimizer_switching_broadcast.md
sha256 87f517bdf945f71e2c45a82bcdb75bee4c55d5498b067a86845339b0b7a5c5ea
```

The principal archive sources used for comparison were Theorems 21.8,
21.29, 21.34, and 36.3, together with
`bounded_cap_linear_response_rate.md`,
`nearmin_pinned_response_packing.md`, and
`bounded_cap_switching_response_probe.md`.

## 2. SB.1: bipartite Boolean norm

**Verdict: PASS.**

For a sign matrix `R in {+-1}^{k times ell}` and uniform Boolean `q`, the
sharp real Khintchine inequality at exponent one gives, row by row,

```math
 \mathbb E_q |(Rq)_i|
 \ge {1\over\sqrt2}
       \left(\sum_jR_{ij}^2\right)^{1/2}
 =\sqrt{\ell/2}.
```

Therefore

```math
 \mathbb E_q\|Rq\|_1\ge k\sqrt{\ell/2}.
```

For a realizing `q`, choosing `p_i=sign((Rq)_i)` makes
`p^TRq=\|Rq\|_1`.  This proves SB.1 with no polarization factor.  Applying
the same proof to `R^T` gives the advertised transposed bound.

## 3. SB.2: random projective code

**Verdict: PASS AFTER A MINOR QUANTIFIER/ROUNDING REPAIR.**

For independent uniform labels,

```math
 d_{\rm P}(u,v)={n-|u\mathbin\cdot v|\over2}.
```

Hoeffding gives

```math
 \Pr\{|u\mathbin\cdot v|>n/2\}
 \le2e^{-n/8}.
```

Sample

```math
 N=\lceil e^{cn}\rceil
```

labels.  The pair union bound tends to zero whenever `2c<1/8`; the source's
`c<1/32` is safely inside that range.  On the good event no two labels are
equal or antipodal, so projection to
`{+-1}^n/{+-1}` loses nothing, and every pair obeys

```math
 n/4\le d_{\rm P}(u,v)\le n/2.
```

The theorem should say **for every sufficiently large `n`**.  As written,
“every order-`n` signing” is formally false at `n=1`, and the random-code
argument is only asymptotic.  Using `ceil(e^(cn))` also makes the literal
cardinality `|U|>=e^(cn)` rather than only its exponential rate explicit.

These are statement-hygiene repairs, not changes to the asymptotic theorem.

## 4. SB.2: switching difference and the factor `1/4`

**Verdict: PASS.**

Let `s=u odot v` and choose the projective representative with

```math
 S=\{i:s_i=-1\},\qquad k=|S|\in[n/4,n/2],
 \qquad \ell=n-k\ge n/2.
```

Conjugation by a Boolean diagonal is an isometry for `\|.\|_B`, because it
is the change of variables `x\mapsto v odot x`.  Hence

```math
 \|A^u-A^v\|_B=\|A^s-A\|_B.
```

The coefficient difference is zero on edges internal to `S` and to its
complement, and is exactly `-2A_ij` on every crossing edge.  If `R` is that
`k by ell` sign block, then for an old spin `(p,q)`,

```math
 H_{A^s-A}(p,q)=-2p^TRq.
```

There is no matrix-versus-edge factor missing: the Hamiltonian sums each
unordered crossing edge once.  SB.1 now gives

```math
 \begin{aligned}
 \|A^u-A^v\|_B
 &\ge2k\sqrt{\ell/2}\\
 &\ge2{n\over4}\sqrt{{n/2\over2}}
 ={1\over4}n^{3/2}.
 \end{aligned}
```

Switching acts freely on complete signings modulo the global label for
`n>=2`: a nonconstant switch changes every coefficient across a nonempty
cut.  Switching also merely permutes the Boolean energy table, so every
child has exactly the same cap.  If `A` is an exact minimizer, every `A^u`
is therefore an exact minimizer.

This establishes the universal target-scale old-block packing independently
of shell geometry.

## 5. SB.6--SB.7: exact amplitude-`n` pin

**Verdict: PASS.**

If `x` differs from `r` on `d` vertices, exactly `d(n-d)` quadratic
monomials change.  Thus, for every complete signing `C`,

```math
 H_C(x)-H_C(r)\le2d(n-d).
```

The field changes by

```math
 nr^Tx-nr^Tr=-2nd.
```

The total objective change is at most

```math
 2d(n-d)-2nd=-2d^2\le0.
```

Therefore

```math
 R_C(nr)=n^2+H_C(r)
```

exactly, including ties at the declared boundary amplitude.  A Boolean-norm
witness for `A^u-A^v` can be sign-oriented by swapping `(u,v)`, so SB.7
follows with the same `n^(3/2)/4` constant.

The collection of pair witnesses is a common, predeclared query bank.  If
two switched minimizers shared one summary state, their common decoder value
at the corresponding `nr_(uv)` could not approximate both responses below
half the gap.  Hence the exponential coefficient packing indeed yields
`Omega(n)` bits of exact pinned-response entropy.

This remains a one-sided amplitude-`n` response theorem.  As the source
correctly says in Section 4, realizing the pin by free exact-sign auxiliary
spins incurs a quadratic calibration in the known universal construction.

## 6. What Theorem 21.8 actually gives from SB.8

### 6.1 Boundary-response implication

**Verdict: PASS.**

After replacing `A` by `-A` when necessary, an exact minimizer satisfies

```math
 P=\max_xH_A(x)=Q(A)=M_n.
```

SB.8 is exactly hypothesis (21.38), uniformly over this class.  Theorem
21.8 then supplies:

1. one exact sign bridge `B in {+-1}^{n times n}` with
   `\|B\|_(2 to2)=O(sqrt n)`;
2. `exp(gamma n)` switches `s`; and
3. pairwise separation

   ```math
   d_{\rm proj}(P_BH_{A^s},P_BH_{A^t})
   \ge d n^{3/2},
   ```

   where

   ```math
   (P_BH_{A^s})(y)
   =\max_x\{H_{A^s}(x)+x^TBy\}.
   ```

All children remain exact minimizers.  This is a rigorous bounded-operator
exact-sign **boundary-response** packing.

### 6.2 Exact complete parents and their cap

**Verdict: REPAIR NEEDED, THEN PASS FOR CONDITIONAL PROFILES.**

The bridge alone does not make a complete signing on `2n` vertices because
the new--new edges are absent.  Choose one common hollow sign matrix `C_n`
with

```math
 Q(C_n)=O(n^{3/2}).
```

Such a matrix exists at every sufficiently large order by the elementary
random-sign union bound.  Define

```math
 P_s=
 \begin{pmatrix}
 A^s&B\\ B^T&C_n
 \end{pmatrix}.
```

This is a complete hollow exact signing of order `2n`, and

```math
 \begin{aligned}
 Q(P_s)
 &\le Q(A^s)+
       \max_{x,y}|x^TBy|+Q(C_n)\\
 &\le M_n+n\|B\|_{2\to2}+Q(C_n)
 =O(n^{3/2}).
 \end{aligned}
```

Its conditional one-sided boundary profile is

```math
 F_s(y)=H_{C_n}(y)+(P_BH_{A^s})(y).
```

Because the first term is common, every pairwise projective response
distance from Theorem 21.8 survives exactly.  Thus SB.8 really does produce
complete bounded-cap parents with separated *conditional boundary profiles*.
The source should add this common fill and state the response language.

### 6.3 All-spins-free physical exposure

**Verdict: FAIL.**

The scalar all-spins-free response is obtained by optimizing `y` as well as
`x` (and, for the signing cap, taking the outer absolute value):

```math
 Q(P_s)=\max_{x,y}
 |H_{A^s}(x)+x^TBy+H_{C_n}(y)|.
```

Large projective distance between the functions `F_s(y)` and `F_t(y)` does
not imply that their maxima, minima, or maximum absolute values are
separated.  Each switched child may choose a different boundary optimizer,
and the common optimization can erase the profile difference.  Theorem 21.8
contains no selector forcing the query-linked `y` to remain active.

This limitation is explicit in the archived
`bounded_cap_linear_response_rate.md`: its response lower bound applies to
the declared external Boolean query `y`, and is *not automatically* a lower
bound for a context class which forbids pinning fields.  The source draft's
claim that Theorem 21.8 “already supplies the restricted compiler” for
bounded-cap all-spins-free physical exposure is therefore unsupported.

The correct implication chain is

```text
L_tail
  + Theorem 21.8
  => linear-rate exact-sign bounded-operator boundary-response packing

that packing
  + MISSING low-cap boundary selector / restricted anti-pin
  => all-spins-free bounded-cap physical contextual packing.
```

The missing selector must expose a chosen boundary coordinate without the
quadratic calibration of the exact amplitude-`n` pin.  It is precisely the
kind of restricted anti-pin which the source says is still needed before it
prematurely declares it supplied.

## 7. Is `L_tail` genuinely weaker and sufficient?

### Information content

SB.8 retains only the cardinality of one fixed-width upper level set of a
supplied exact minimizer.  It does not identify all response values, compare
orders, or ask for `M_n` numerically.  In that representation sense it is
strictly less information than the Boolean landscape and is a natural,
finite-testable structural hypothesis.

Its proof complexity is not known to be smaller.  It quantifies over every
exact minimizer and uses its top value `P`; no theorem currently derives it
from exact minimality.  Thus “strictly less information” should not be
upgraded to “demonstrably easier theorem.”

### Sufficiency classification

* For Theorem 21.8's external/boundary response packing: **sufficient**.
* For complete bounded-cap parents with separated conditional profiles:
  **sufficient after the common `C_n` fill above**.
* For the claimed all-spins-free physical response packing: **not
  sufficient**; a boundary selector remains missing.
* For convergence or cross-order transfer: **not sufficient and not
  claimed to be**.

The construction in Theorem 21.8 chooses one source maximizer `u_*` to link
switches to queries.  This is an existence witness at the supplied order,
not optimization at a new target order, so it does not create target-order
circularity.  It does mean the construction is not an optimizer-free
algorithmic compiler.

## 8. SB.9: bounded operator norm implies the tail deficit

**Verdict: PASS AFTER MAKING ONE UNIFORM LOWER-BOUND STEP EXPLICIT.**

For a hollow sign matrix,

```math
 \|A\|_F^2=n(n-1),
 \qquad \mathbb EH_A(U)=0.
```

After orienting the exact minimizer positively, use the known universal
lower bound

```math
 P=M_n\ge c_*n^{3/2}
```

for all sufficiently large `n`, with an absolute `c_*>0`.  Choose
`0<d_0<c_*`.  Membership in the SB.8 shell implies

```math
 H_A(U)>(c_*-d_0)n^{3/2}.
```

If `\|A\|_(2 to2)\le C\sqrt n`, the Rademacher Hanson--Wright inequality
gives

```math
 \Pr\{H_A(U)>(c_*-d_0)n^{3/2}\}
 \le2\exp\left[-c
 \min\left\{
 {n^3\over\|A\|_F^2},
 {n^{3/2}\over\|A\|_{2\to2}}
 \right\}\right]
 \le e^{-\kappa n}.
```

Multiplication by `2^n` proves SB.8.  The source's claim is therefore
correct, but it should mention the positive universal lower bound and choose
`d_0` below it.  The sparse-edit warning is appropriately scoped and does
not affect the exact-minimizer implication.

## 9. Archive duplication audit

* The bipartite Khintchine mechanism in SB.1 appears repeatedly in the
  archive.  Its application to the difference of two arbitrary switchings
  is elementary but I found no archived theorem stating SB.2's universal
  exponential target-scale norm packing.
* The switching action's freeness is already used in Theorem 21.34, but
  there it counts orbits of all signings rather than giving a separated
  subcode inside one orbit.
* SB.6 is the same exact amplitude-`n` pin identity used in Theorem 36.3.
  SB.2 plus SB.6 is new only as the exact-minimizer switching specialization.
* Theorem 21.29 already embeds *any* old-block Boolean norm isometrically
  into complete exact-sign coordinate contexts, but with a quadratic common
  calibration.  Combining it with SB.2 gives another high-cap exact
  contextual packing, not the desired bounded-cap result.
* SB.8 is exactly a proposed uniform exact-minimizer instance of Theorem
  21.8's hypothesis.  The tail-to-response mechanism itself is archived and
  should not be presented as a new compiler theorem.
* BCX supplies a genuinely all-spins-free bounded-cap anti-pin only for a
  special regular-Hadamard switching code.  It does not furnish the missing
  arbitrary-exact-minimizer selector.

## 10. Final result table and required source repairs

| Claim | Verdict | Required repair or scope |
|---|---|---|
| SB.1 Khintchine bound | PASS | none |
| random projective code | PASS | say “sufficiently large `n`”; use `ceil` or weaken literal cardinality |
| SB.2 Boolean-norm factor `1/4` | PASS | none |
| switching freeness/cap preservation | PASS | `n>=2`, automatic asymptotically |
| SB.6 exact pin identity | PASS | none |
| SB.7 pinned response packing | PASS | one-sided amplitude-`n` query language |
| SB.8 invokes Theorem 21.8 | PASS | only for boundary-response profiles |
| exact complete parent order/cap | PASS AFTER REPAIR | add one common bounded-cap `C_n`; state conditional profile |
| bounded-cap all-spins-free exposure | FAIL | requires an additional boundary-selector/anti-pin lemma |
| SB.9 spectral sufficient condition | PASS AFTER REPAIR | invoke uniform `M_n>=c_*n^(3/2)` and choose `d_0<c_*` |
| “old-block clause of `L_broadcast` solved” | PASS METRICALLY | exposure clause remains open in physical language |
| `L_tail` is the remaining physical SML | FAIL | `L_tail + L_selector` is the honest sufficient pair |

The source can be canonicalized after changing its headline and frontier
classification to distinguish:

1. the proved universal exact-minimizer switching norm and pinned packings;
2. the conditional bounded-cap **boundary-response** packing from `L_tail`;
3. the still-open all-spins-free physical selector.

Without that repair, the main mathematical novelty SB.2 remains valid, but
the claimed completion of physical exposure does not.

## 11. Post-repair disposition

**Final verdict on the repaired source: PASS.**

The source was rechecked after the audit repairs.  Its post-repair frozen
hash is

```text
extremal_information/drafts/exact_minimizer_switching_broadcast.md
sha256 5278e6cb96a3a554141fe52cfd31dbb1ca38cf7b1260a33a554c116bf6074e8f
```

This differs from the audited post-repair hash only by the final status line,
which now records this audit disposition; no mathematical content changed.

The repaired version now:

1. states SB.2 only for sufficiently large `n`;
2. labels SB.8 as sufficient for boundary-response exposure, not for a
   scalar all-spins-free cap gap;
3. adds a common bounded-cap new--new fill before claiming complete
   order-`2n` exact-sign parents;
4. explicitly retains the low-cap boundary selector/restricted anti-pin as
   an additional missing lemma for all-spins-free exposure; and
5. invokes a uniform lower bound `M_n>=c_*n^(3/2)` and chooses `d_0<c_*`
   in the Hanson--Wright implication from SB.9.

The stale Section 5 heading and opening sentence were also repaired after a
first post-repair read.  The source no longer claims that Theorem 21.8 alone
solves physical exposure.  All formulas SB.1--SB.9, the conditional-profile
parent cap/order accounting, and the final frontier classification now agree
with this audit.
