# Archive verification: restriction and action proposals

Date: 2026-08-16.  Scope: frozen reports
`reports/banach_proposals.md` and `reports/graph_limits_proposals.md`, checked
against the complete ledger and the cited artifacts.  Here `Q(A)` is the
one-copy cap `max_x |sum_{i<j}a_ij x_i x_j|` and
`b_n=M_n/n^(3/2)`.

## Bottom line

No proposed architecture is `A`.  The classifications are:

| Frozen item | Class | Decisive implication-level reason |
|---|---|---|
| Banach I: PRS | **C** | Its `n >> m` quantifier already contains `b_m <= liminf b_n+o(1)`; the archive proved that far-down competitive parents contain every fixed signed pattern, so this is not a new restriction mechanism. |
| Banach II: joint block gluing (4) | **C** | It is the archived `M_n^(2/3)` almost-subadditive target with the exact full parent bridge objective still present.  The joint statement is not disproved, but no obligation is removed. |
| Banach III: tensor-power regularization | **D** as a route | The imported tensor theorem controls the independent bilinear norm only to a root-exponential rate on one power subsequence.  Fixed same-spin polarization loss, `N^{o(1)}` prefactor uncertainty, and lack of all-order/minimum control rigorously prevent the advertised theorem from implying the target. |
| Graph A, clause SR alone | **B** | Fixed-`C`, `o(n^(3/2))` spectral regularization is strictly stronger than the archived two-limit tradeoff and is neither proved nor falsified. |
| Graph A, clause AR alone | **C** | This is verbatim the archived all-order bounded-operator signed-action realization obligation, only with a stronger universal quantifier. |
| Graph A, combined `L_FA` | **C** | Its implication is correct, but it leaves the old AR obligation intact and strengthens rather than reduces the old purification clause. |
| Graph B: IR | **C** | By the archive's uniform induced-pattern theorem, IR is equivalent to convergence (for a fixed competitive bound `R`), not a new exchangeability theorem. |
| Ordinary graphon / bare-action alternative | **D** | The zero graphon loses the `n^(3/2)` coefficient, and a scalable Boolean spike is action-invisible without uniform integrability. |

There is one result worth retaining independently: the graph report's
quantitative action-continuity bound is a **verified new quantitative
refinement** of an archived qualitative continuity statement.  It is not a
convergence architecture and therefore is not an `A` route.

## Cross-cutting gate audit

| Item | Full parent/coset optimization | Full coset histogram | Channels / leading loss | All orders / near-minimizers | State size / conference scope |
|---|---|---|---|---|---|
| PRS | Recovers the unknown scalar optimum comparison in its far-down regime; no maximizing spin is recovered. | No; the order-six noninjectivity witness is correct but irrelevant to the scalar equivalence. | PRS itself pays jointly and has no fixed loss; Rudelson--Vershynin leaves a fixed `Theta(m^(3/2))` term and fixed norm constants. | Uniform over every `Q(A)<=Cn^(3/2)` and every `n>=m`; this excessive uniformity is part of the collision. | No advertised finite state; not conference-only. |
| Joint gluing | **Yes:** `Q(C_D)` is the full parent maximum and bridge selection retains its state-dependent landscape. | Not formally proved necessary, but exact scalar child caps do not close; exact profile closure reads the full growing Boolean action. | The proposed lemma is joint.  Triangle/KSZ and known factorization attempts separately pay a leading rectangular channel. | All orders but only stated for exact optimal children. | Exact bridge constraints are indexed by exponentially many `(x,y)` unless a new closure theorem is supplied; not conference-only. |
| Tensor powers | Controls the wrong, independent-input optimization; it does not recover the same-spin parent minimum. | No explicit histogram. | **Yes:** real polarization leaves a fixed interval, while the root theorem leaves `N^{o(1)}` rather than `1+o(1)`. | One fixed-base subsequence `N=n^k`, not all orders or arbitrary near-minimizers. | Large tensor norm is an exponential Boolean optimization; the theorem only compresses its exponential rate.  Algebraic/power-subsequence, not specifically conference-only. |
| SR | References only `M_n`, `Q`, and operator norm; no parent maximizer/histogram. | No. | SR itself has no loss.  The archived construction has a fixed loss for every fixed operator bound. | Existential at every order, not uniform over arbitrary near-minimizers. | Matrix state; not conference-only. |
| AR / `L_FA` | Does not recover exact finite optimization, but all-order recovery of an extremal limit is exactly the missing cross-order obligation. | No: `o(n^(3/2))` edge changes alter all labeled coset distances while vanishing in the action quotient. | Joint absolute functional; no separate-channel factor. | AR explicitly requires every sufficiently large order. | The action profile is an infinite hierarchy with no proved finite/polynomial closure, but no exponential lower bound is proved.  Not conference-only. |
| IR | Far-down local universality recovers `M_n` itself from the family of restrictions. | Not one parent's full energy histogram; it does expose all `2^(n choose 2)` fixed signed patterns. | No separate payment or fixed loss in IR. | Uniform over all competitive hosts once `N>=N_0(n,R,epsilon)`. | Its far-down mechanism ranges over the full signed-pattern universe; not conference-only. |

## 1. Banach Architecture I: PRS -- class C

The frozen theorem is (`banach_proposals.md`, Section 4):

```math
Q(A[S])
<=m^{3/2}\left({Q(A)\over n^{3/2}}+\omega_C(m)\right),
\qquad \omega_C(m)\to0,                              \tag{PRS}
```

for every `C`, every `n>=m`, and every complete signing
`Q(A)<=Cn^(3/2)`, with some `|S|=m`.  Equivalently the report defines

```math
\Delta_C(m)=\sup_{n\ge m}\sup_{Q(A)\le Cn^{3/2}}
\left[\min_{|S|=m}{Q(A[S])\over m^{3/2}}
-{Q(A)\over n^{3/2}}\right]_+\to0.                  \tag{B.6}
```

The conditional implication printed in the report is correct.  Its claimed
strictness as a research reduction is not.  Choose `C` above the uniform
upper bound and take optimal ambient orders `N_j` with `b_{N_j}->ell`, where
`ell=liminf b_N`.  Since every restriction has cap at least `M_m`, already

```math
\Delta_C(m)\ge [b_m-\ell]_+.
```

Thus the required vanishing defect contains the desired scalar limsup/liminf
comparison before any restriction-selection theorem is used.

The exact prior appearance is ledger Section 10.113.2, “Far-down principal
restriction,” equations (10.1392)--(10.1394).  For every competitive parent,
the archive proves

```math
\left|\sum_{i\in U,j\in V}a_{ij}\right|
\le {3\over2}Q(A_N)                                  \tag{10.1393}
```

for disjoint shores.  Hence its signed graph is uniformly cut-quasirandom at
error `O_R(N^(-1/2))`.  The induced counting lemma then gives, for every fixed
`m` and all sufficiently large `N`,

```math
\min_{|S|=m}Q(A_N[S])=M_m.                           \tag{10.1394}
```

So at `N >> m`, PRS does not harvest inherited structure: the parent contains
every order-`m` signing, including an optimizer, and PRS reduces exactly to
the unknown scalar comparison.  At fixed density it returns to the archived
optimized principal-restriction target, ledger Section 10.23 equation
(10.60), and the hidden-versus-revealed selector problem summarized in ledger
Sections 10.72--10.113.  This is an exact implication collision, not a
vocabulary collision.

The report's strictness witness is nevertheless arithmetically correct.  I
independently enumerated its two order-six matrices: both have `Q=11` and
principal minima `(1,3,4,6,11)`, while their absolute-energy counts are
respectively `(28,20,10,2,2,2)` and `(28,16,14,4,0,2)` at energies
`(1,3,5,7,9,11)`.  It proves that these scalars do not recover the histogram;
it does not prove that PRS is weaker than convergence.  The finite anchor is
also correct:

```math
{17\over11^{3/2}}-{21\over14^{3/2}}
=0.0650802151\ldots .                                \tag{B.8}
```

The Rudelson--Vershynin specialization is normalized correctly.  Their
hollow principal-submatrix theorem gives

```math
E||A[Q]||_C
<=K\{p^2||A||_C+p^{3/2}(||A||_{Col}+||A^T||_{Col})\}, \tag{B.9}
```

and for a full sign matrix the column term is `Theta(m^(3/2))`.  Cut/bilinear
to same-spin conversion also carries fixed constants.  Therefore this import
does not prove PRS; it has exactly the fixed leading error PRS forbids.

## 2. Banach Architecture II: joint block gluing -- class C, not D

The report proposes, for optimal children and a chosen sign bridge `D`,

```math
Q(C_D)^{2/3}
<=Q(A)^{2/3}+Q(B)^{2/3}+c(n+m)^{1-\delta}.           \tag{B.4}
```

Its almost-Fekete deduction is correct: with `a_n=M_n^(2/3)`, the defect
`O(N^(1-delta))` satisfies the geometric summability criterion.

This is already the exact archived target.  Writing locally
`u_n=M_n^(2/3)`, ledger Section 10.115.2 proves that

```math
u_{m+n}<=u_m+u_n+e(m+n),\qquad
\sum_{j\ge1}{e(2^jk)\over2^jk}\to0                 \tag{10.1406--10.1407}
```

implies convergence, and notes that `e(N)=O(N^(1-delta))` suffices.  Ledger
Section 10.113.2 equation (10.1391) had already proposed exactly this
`2/3`-power architecture.

More importantly, ledger Sections 1.6 and 10.115.3 give the exact parent
identity

```math
Q(A*_DB)=\max_{x,y}
\bigl(|H_A(x)+H_B(y)|+|x^TDy|\bigr).                \tag{10.1410}
```

This is the report's equation (3) after optimizing the relative global block
sign.  Thus (B.4) retains the complete state-dependent parent bridge
optimization.  Ledger Section 10.140.1 classifies the corresponding
block-completion candidate as an “equivalent sufficient reformulation”; the
artifact `consolidation_blank_slate_diagnostic_2026_08_15.md`, equations
(1.1)--(1.4), gives the aligned derivation.

The known no-go is narrower than the joint assertion, so `D` would overstate
the archive.  A random/KSZ bridge has rectangular norm
`Theta(sqrt(nm(n+m)))`, a leading channel (ledger Sections 1.8 and 10.115.5),
and triangle inequality or vector-valued factorization pays it separately.
But no archived theorem falsifies a genuinely state-correlated bridge.
Therefore the right class is `C`: old full-parent obligation, not a proved
impossibility.

## 3. Banach Architecture III: full tensor powers -- class D as a route

The exact imported statement is

```math
\lim_{k\to\infty}
\|\phi_T^{\otimes k}:
(\ell_\infty^n)^{\otimes_\varepsilon k}
\to(\ell_1^n)^{\otimes_\pi k}\|^{1/k}
=\gamma_2^*(\phi_T).                                \tag{B.5}
```

The Banach-space identifications with `ell_infinity^(n^k)` and
`ell_1^(n^k)` are correct.  So is the diagonal observation: hollowing
`T^(tensor k)` changes the quadratic polynomial by the constant
`tr(T)^k/2`, whose magnitude is at most `N/2=o(N^(3/2))`.

The theorem nevertheless controls the independent-input norm.  The report's
own exact audit gives

```math
2Q(A)<=||A||_{infinity to 1}<=4Q(A).                \tag{B.1}
```

That fixed interval cannot identify a same-spin limiting coefficient.
Moreover a `k`th-root limit permits an `exp(o(k))=N^{o(1)}` prefactor, far
larger than `1+o(1)`, and concerns only orders `N=n^k` for one fixed base,
with no lower control on the minimum.

The closest exact archive obstruction is ledger Section 3.8:

```math
Q(A\otimes B)>=176>144=Q(A)Q(B),                    \tag{ledger 3.8}
```

which falsifies the tempting same-spin tensor submultiplicativity.  It does
not falsify the Aubrun--Mueller-Hermes bilinear theorem; it proves that the
missing same-spin transfer is substantive.  The fixed polarization boundary
is also proved in
`artifacts/symmetric_grothendieck_krivine_obstruction_audit.md`.  Hence the
import is valid, but the architecture as frozen has no implication to
convergence and is rigorously insufficient.

## 4. Graph Architecture A: spectral action recovery

### 4.1 SR alone -- class B

The report asks for fixed `C` and `epsilon_n->0` such that at every order

```math
Q(A_n)<=M_n+\epsilon_n n^{3/2},\qquad
||A_n||_{op}<=C\sqrt n.                              \tag{SR}
```

This exact statement is not in the archive.  The nearest proved theorem is
ledger Section 10.29 equation (10.75), equivalently the boxed tradeoff in
`artifacts/optimality_excess_and_regularization.md`, Section 2:

```math
Q(A'_n)<=M_n+O(K^{-1/2}n^{3/2}),\qquad
||A'_n||_{op}=O(K\sqrt n).                           \tag{archive SR tradeoff}
```

For each fixed operator bound this leaves a fixed leading loss.  Sending the
loss to zero sends the operator bound to infinity, so a diagonal choice does
not prove SR.  Conversely, the unbounded-spike example in
`artifacts/optimality_excess_and_regularization.md`, Section 1, only shows
that arbitrary asymptotically minimizing sequences need not be spectrally
regular; SR is existential and survives it.  SR is therefore a genuinely
stronger variant with one isolated open issue: `B`, not `D`.

SR does not reconstruct a coset histogram and pays no separate channels in
its statement.  A known proof attempt purifies exceptional blocks and cross
edges separately and consequently retains the fixed `O(K^(-1/2))` loss.

### 4.2 AR alone -- class C

The frozen clause is

```math
\forall T\in\mathscr L_C\ \forall\varepsilon>0\ \exists N\
\forall m\ge N\ \exists B_m\in\mathfrak A_m:
||B_m||_{op}<=C_T\sqrt m,\quad d_M(T_{B_m},T)<=\varepsilon. \tag{AR}
```

This is the exact missing theorem already isolated in ledger Sections 10.29
and 10.31.6.  `artifacts/bounded_op_signed_realization.md`, Section 4, states:

> Every extremal bounded-operator action limit arising from symmetric
> off-diagonal sign matrices is approximable, with objective convergence, by
> such sign matrices at every sufficiently large order.

The report strengthens “extremal” to every `T in L_C`; it does not remove the
obligation.  The residual-variance identity in that artifact,

```math
||R_ij||_F^2=k^2-r_ij^2=(1-o(1))k^2,                \tag{artifact (2)}
```

rigorously defeats naive rational sign blow-ups, independent residuals, and
Hadamard residuals, but does not falsify an unknown absorbing realization
theorem.  Thus AR is `C`, not `D`.

### 4.3 Combined `L_FA` -- class C; implication verified

The proof from SR+AR to convergence is correct.  Backhausz--Szegedy
compactness applies to `T_A=A/sqrt(n)` under the fixed `2 to 2` bound; the
limit is representable and self-adjoint.  Continuity gives
`Phi(T)=2 liminf b_n`, and AR gives exact sign realizers at every sufficiently
large order with objective tending to that value.

This is nevertheless the same architecture as ledger Section 10.29:
purification, compactness/continuity, and the still-missing all-order
realization.  In fact the archived two-limit purification is enough if AR is
available for every fixed bounded class: prove a limsup bound with fixed
purification error and then send that error to zero.  Replacing it by the
stronger SR does not reduce AR.  The combined route is therefore `C`.

The normalization and stability constants in the frozen report are correct:

```math
\Phi(T_A)={2Q(A)\over n^{3/2}},\qquad
||T_A||_{infinity to 1}={B(A)\over n^{3/2}}
<={4Q(A)\over n^{3/2}},                              \tag{G.1--G.2}
```

and if `A,A'` differ on `r` unordered edges,

```math
|Q(A)-Q(A')|<=2r,\qquad
d_M(T_A,T_A')<=6\sqrt{r/n^{3/2}}.                   \tag{G.8--G.9}
```

The last bound follows from the exact `infinity to 1` difference
`<=4r/n^(3/2)` and Backhausz--Szegedy Lemma 2.19.  It validates the claim that
the action quotient cannot reconstruct the exact labeled energy/coset table.

## 5. New supporting result: quantitative action continuity -- verified

The report proves, for `||S||_(2 to 2),||T||_(2 to 2)<=C`,

```math
|\Phi(S)-\Phi(T)|<=5C\sqrt\delta+\delta,
\quad
\delta=d_H^{LP}(\overline{S_1(S)},\overline{S_1(T)}), \tag{G.5}
```

and consequently

```math
|\Phi(S)-\Phi(T)|
<=5C\sqrt{2d_M(S,T)}+2d_M(S,T).                     \tag{G.6}
```

I checked the proof.  Under a Strassen coupling, truncate the output at
`R`; the two tails cost `2C^2/R`, the good event costs `(R+1)delta`, and the
exceptional event costs `2R delta`, giving

```math
2C^2/R+(3R+1)delta.
```

With `R=C/sqrt(delta)` this is exactly `5C sqrt(delta)+delta` (with `C=0`
separate).  Uniform `L^2` bounds make `xy` uniformly integrable also on the
closed profiles, so taking Hausdorff-matched suprema is legitimate.  Finally
the `k=1` term has weight `1/2` in `d_M`, hence `delta<=2d_M`.  For large
`d_M`, the displayed bound extends trivially from `Phi<=C`; only the local
modulus is used in the implication.

The archive already asserted qualitative continuity in ledger Section 10.29
and `artifacts/concentration_compactness_boolean_profiles.md`, Section 6.
The explicit modulus and direct closed-profile proof are new and correct.
They should be retained as a proved supporting theorem, but they do not solve
SR or AR.

The report's bare spike test is also correct:
for `k=floor(n^(3/4))` and a padded clique `S_n`,

```math
d_M(S_n,0)<=k/n->0,\qquad
\Phi(S_n)>=k(k-1)/n^{3/2}->1,
```

while `||S_n||_(2 to 2)~n^(1/4)`.  The stronger full-sign version is already
proved in `artifacts/action_convergence_boolean_spikes.md`, equations (4)--(7).

## 6. Graph Architecture B: IR -- class C

The frozen lemma is

```math
\forall R,\varepsilon\ \exists n_0\ \forall n\ge n_0\ \exists N_0
\ \forall N\ge N_0\ \forall A:\quad
q_N(A)\le R
\Longrightarrow
\exists |S|=n:\ {Q(A[S])\over n^{3/2}}\le q_N(A)+\varepsilon. \tag{IR}
```

This is the archive's far-down restriction candidate (ledger Section
10.113.2, equation (10.1392)) without an explicit rate.  Equations
(10.1393)--(10.1394), quoted above, imply uniformly for every such `A` and
fixed target `n` that, once `N` is large enough,

```math
\min_{|S|=n}Q(A[S])=M_n.
```

Thus IR becomes simply

```math
b_n<=q_N(A)+\varepsilon.
```

Taking `A` optimal and `N` down a liminf sequence proves convergence.  In the
other direction, if `b_n->L`, choose `n_0,N_0` so that
`b_n<=L+epsilon/3` and `b_N>=L-epsilon/3`; since every `q_N(A)>=b_N` and an
optimal target pattern occurs as an induced restriction, IR follows (after
harmlessly reallocating `epsilon`).  With any fixed `R` above the known
competitive bound, IR is therefore equivalent to convergence plus the
already verified local-universality theorem.

Diaconis--Janson exchangeability does not supply the missing growing-`n`
inequality.  Fixed restrictions see only the fair signed graphon; allowing
`n` to grow enough to control the `2^n` spin supremum is exactly where the
scalar optimum reappears.  IR is `C`, with no separate-channel or conference
escape.

## 7. Ordinary graphons and Gamma limits -- class D (correct negative audit)

For every competitive signing, the frozen report derives

```math
||W_A||_square<=B(A)/n^2<=4Q(A)/n^2=O(n^{-1/2}).     \tag{G.3}
```

Hence all bounded signed graphons converge to zero, and any bounded-weight
Gamma limit is normalized at `n^2` and loses the target coefficient.  This is
exactly archive collision-index item E15 and
`artifacts/action_convergence_boolean_spikes.md`.  Scaling to
`sqrt(n)W_A` restores the action size but destroys tightness unless a moment
bound such as SR is supplied.  Therefore the report is right not to offer an
ordinary-graphon third architecture.

## Final verifier judgment

The archive comparison leaves no executable `A` architecture.  PRS and IR
are especially deceptive: their scalar/non-histogram formulations are real,
but their far-down quantifiers encode the desired optimum comparison.  Joint
block gluing remains an open old parent-bridge problem, not a disproved one.
The tensor-power theorem is valid but rigorously too coarse.  The action
proposal contributes a sound quantitative continuity lemma and a genuinely
stronger standalone SR target, but its convergence route still terminates at
the exact archived AR obligation.
