# Independent all-order recovery proposals: rounding, designs, and conference carriers

Date: 2026-08-16.

Scope: this is a fresh, ledger-blind proposal note based only on the design/quasirandom
toolkit, the matrix-rounding toolkit, and the minimal all-order recovery note.  The
three architectures below are frozen for this round; no fourth architecture is being
held implicitly in reserve.

> **Post-freeze verifier correction.**  Section 2.2 correctly observes that a
> bounded-degree `O(n)` leave has `o(1)` normalized operator defect, but bounded
> degree is not necessary for the recovery objective or action profile.  Any
> `O(n)`-edge repair has `O(n)` objective error,
> `O(n^{-1/2})` normalized `L^infinity -> L^1` error, and only an `O(1)`
> normalized operator-norm increment by Frobenius.  It is therefore already a
> valid repair module under a common finite operator bound.  This strengthens
> the repair half and does not supply the missing universal-profile enforcement
> half.

## 0. Problem and decision

For a symmetric hollow sign matrix (A) of order (n), write

\[
 H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
 Q(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|,\qquad
 M_n=\min_A Q(A).
\]

Put (T_A=A/\sqrt n), with the uniform probability measure on the vertices, and

\[
 \Phi(T_A)=\sup_{|f|\le1}|\langle f,T_Af\rangle|
           ={2Q(A)\over n^{3/2}}.
\]

The objective is to recover a selected purified liminf cluster densely enough in
the orders to prove that (M_n/n^{3/2}) converges.

The frozen ranking is:

1. **Low-variance weighted outer recovery followed by scalar/bilinear rounding.**
   This is the best reduction.  Its rounding half is already a theorem; its weighted
   recovery half is a clean new lemma.
2. **Finite-gadget outer-profile packing followed by bounded-degree design repair.**
   The repair half is stronger than a mere edge-count estimate, but the universal
   profile property is not supplied by any design theorem in the toolkit.
3. **Conference-fibre/quasirandom lifts of a cluster approximant.**
   This is an attractive all-order carrier but the fibre modes are a serious, readily
   falsifiable obstruction.  It should not be developed before passing the proposed
   finite tests.

## 1. Proposal I (selected): low-variance weighted outer recovery

### 1.1 Exact known inputs

The proposal uses four inputs, three already proved in the recovery note and one
an elementary specialization of scalar Bernstein rounding.

**Purified cluster input.**  If

\[
 L=\liminf_n {M_n\over n^{3/2}},
\]

then, for every fixed tolerance (eta>0), there is a bounded action-limit object
(T_\eta), obtained from exact signings, such that

\[
 2L\le \Phi(T_\eta)\le 2L+\eta.                    \tag{1.1}
\]

The norm bound may depend on (eta).

**Directed one-profile continuity.**  If

\[
 \|S\|_{2\to2},\|T\|_{2\to2}\le D,
 \qquad \partial_1(S,T)\le\delta,
\]

then

\[
 \Phi(S)\le\Phi(T)+5D\sqrt\delta+\delta.          \tag{1.2}
\]

Only the outer inclusion from the one-profile of (S) toward that of (T) is
used.

**Near-order deletion.**  If (N\le m), every (N\)-vertex principal
submatrix (C[S]) of a signing (C) satisfies

\[
 Q(C[S])\le Q(C).                                  \tag{1.3}
\]

Indeed, extend a fixed spin on (S) by independent uniform signs and average.
Consequently a construction on an upward ratio-dense set of orders transfers to
all orders without an additive energy loss.

**Biased sign rounding with an explicit bilinear bound.**  Let (W=(w_{ij}))
be symmetric, hollow, and in ([-1,1]).  Independently on unordered edges choose

\[
 A_{ij}=1\quad\hbox{with probability }{1+w_{ij}\over2},
 \qquad A_{ij}=-1\quad\hbox{otherwise},
\]

and put

\[
 E=A-W,\qquad
 V(W)=\sum_{i<j}(1-w_{ij}^2),\qquad
 B(E)=\max_{x,y\in\{\pm1\}^m}|x^{\mathsf T}Ey|.
\]

Scalar Bernstein and a union bound over (4^m) sign pairs give

\[
 \Pr\{B(E)\ge t\}
 \le 2\,4^m\exp\!\left[-{t^2\over2(4V+4t/3)}\right].       \tag{1.4}
\]

For example, with

\[
 a_m=2m\log2+\log4,
 \qquad
 t_m=\sqrt{8V(W)a_m}+{8\over3}a_m,                         \tag{1.5}
\]

the standard Bernstein inversion makes the right side of (1.4) at most
(1/2).  Hence a supported exact signing exists with (B(A-W)\le t_m).
Since (q(E):=\max_x|x^{\mathsf T}Ex|\le B(E)),

\[
 Q(A)\le Q(W)+{t_m\over2}.                                \tag{1.6}
\]

Thus

\[
 V(W)=o(m^2)\quad\Longrightarrow\quad
 B(A-W)=o(m^{3/2})\quad\hbox{and}\quad
 Q(A)\le Q(W)+o(m^{3/2}).                                 \tag{1.7}
\]

This is the exact scale needed here.  It is also obtainable in the form
(O(\sqrt{mV}+m)) from the toolkit's scalar rounding theorem.  If only (r)
edges are genuinely fractional, Spencer's discrepancy theorem together with
Lovasz--Spencer--Vesztergombi shifted discrepancy gives the parallel
(O(\sqrt{rm})) conclusion.

For comparison, spectral rounding is unnecessarily strong.  The
Bandeira--van Handel independent-entry bound gives an outcome with

\[
 \|A-W\|_{\rm op}\le C(\sqrt{v_{\max}}+\sqrt{\log m}),
 \qquad
 v_{\max}=\max_i\sum_{j\ne i}(1-w_{ij}^2),                 \tag{1.8}
\]

and Wang--Lau--Zhou gives a deterministic supported outcome with leading term
(2\sqrt{v_{\max}}) plus its explicit lower-order terms.  Those theorems reach
(o(\sqrt m)) under sublinear **row** variance.  Total variance
(V=o(m^2)) does not imply (v_{\max}=o(m)), and (1.8) is not needed for the
scalar convergence argument.

### 1.2 The exact new lemma

Call a set ({\cal N}\) upward ratio-dense if its least member (s(N)\ge N)
satisfies (s(N)/N\to1).

> **Lemma L (low-variance weighted outer recovery, LV-WOR).**  There is a
> null sequence (eta_\ell\downarrow0) and, for every (ell), one purified
> cluster object (T_\ell=T_{\eta_\ell}) satisfying (1.1), an upward
> ratio-dense set ({\cal N}_\ell), and symmetric hollow matrices
> (W_m\in[-1,1]^{m\times m}), (m\in{\cal N}_\ell), with numbers
> (D_m\ge\max\{\|T_{W_m}\|_{2\to2},\|T_\ell\|_{2\to2}\}), such that
>
> \[
> \boxed{
> \begin{gathered}
> \delta_m:=\partial_1(T_{W_m},T_\ell)\longrightarrow0,
> \qquad D_m\sqrt{\delta_m}\longrightarrow0,\\[2mm]
> V(W_m):=\sum_{i<j}(1-w_{ij}^2)=o(m^2)
> \quad(m\to\infty,\ m\in{\cal N}_\ell).
> \end{gathered}}                                           \tag{L}
> \]

The order of limits is important: fix (ell), let (m\to\infty), and only
then send (ell\to\infty).  No norm bound uniform in (ell) is requested.

### 1.3 Boxed implication (L\Rightarrow) convergence

For a weighted symmetric hollow (W), the maximum of the absolute
multi-affine quadratic form over ([-1,1]^m) occurs at a cube vertex.  Hence

\[
 \Phi(T_W)={2Q(W)\over m^{3/2}}.                            \tag{1.9}
\]

At each (m\in{\cal N}_\ell), apply (1.4)--(1.6) and select an exact
signing (A_m).  Lemma L and (1.2) give

\[
 \begin{aligned}
 {Q(A_m)\over m^{3/2}}
 &\le {1\over2}\Phi(T_{W_m})+{t_m\over2m^{3/2}}\\
 &\le {1\over2}\Phi(T_\ell)
       +{5\over2}D_m\sqrt{\delta_m}+{\delta_m\over2}+o(1)\\
 &\le L+{\eta_\ell\over2}+o(1).                          \tag{1.10}
 \end{aligned}
\]

For arbitrary large (N), take (m=s_\ell(N)\) and an (N\)-vertex
principal submatrix of (A_m).  By (1.3),

\[
 {M_N\over N^{3/2}}
 \le\left({m\over N}\right)^{3/2}{Q(A_m)\over m^{3/2}}.
\]

Ratio density and (1.10) imply

\[
 \limsup_N {M_N\over N^{3/2}}\le L+{\eta_\ell\over2}.
\]

Letting (ell\to\infty) and using the definition of (L) yields

\[
 \boxed{
 \text{purification + directed continuity + Lemma L + biased rounding}
 \quad\Longrightarrow\quad
 {M_n\over n^{3/2}}\text{ converges}.}                     \tag{1.11}
\]

### 1.4 Why total fractional variance is a genuine reduction

The variance condition is not cosmetic.  It changes an exact-sign recovery
problem into an almost-integral weighted recovery problem, and the remaining
integrality step is completely discharged by (1.4).

More is true than the objective estimate.  By convexity of the cube,

\[
 \sup_{|f|\le1}\|Ef\|_1=B(E).
\]

With normalized vertex measure,

\[
 \sup_{|f|\le1}
 \|T_Af-T_Wf\|_{L^1}
 ={B(E)\over m^{3/2}}=o(1).                                \tag{1.12}
\]

Under the identity coupling, Markov's inequality therefore puts the laws of
((f,T_Af)) and ((f,T_Wf)) within
(O(\sqrt{B(E)/m^{3/2}})=o(1)) in Levy--Prokhorov distance, uniformly in
(f).  Thus the same rounding preserves the one-profile in the topology
actually used here even when it does not give spectral (o(\sqrt m)).

There are two necessary qualifications.

* (V=o(m^2)) plus **externally certified weighted recovery** is a real
  reduction.  (V=o(m^2)) alone says nothing about the target.
* If “weighted recovery” is defined only by
  (Q(W_m)\le(L+o(1))m^{3/2}), it is the scalar missing inequality in new
  notation.  The directed target condition in Lemma L is what keeps the
  state tied to independently obtained limit information.

The distinction from operator rounding is also real.  A subquadratic total
variance can be concentrated in a few rows, leaving (v_{\max}=\Theta(m)).
The matrix-norm theorems then stop at their natural \(\Theta(\sqrt m)\) scale,
whereas (1.4) still gives (o(m^{3/2})) bilinear error.

### 1.5 In what sense Lemma L is strictly weaker

Lemma L is a strict recovery quotient of full exact action recovery:

* only one selected purified liminf cluster is treated at each tolerance;
* recovery is needed only on an upward ratio-dense order set;
* the intermediate object may be weighted;
* only directed one-profile inclusion is asserted, with neither reverse
  inclusion nor joint (k\)-profiles;
* after rounding, the convergence proof uses only the scalar/bilinear defect;
  it never reconstructs the map (x\mapsto H_A(x)).

Full exact action recovery would imply Lemma L by taking (W_m=A_m), for
which (V=0).  Lemma L does not assert full action recovery.  It also does not
ask that (W_m) or the rounded (A_m) minimize the target-order objective.

This is strictness of the requested information and quantifiers, not a claim
that Lemma L is already easy.  Once combined with the proved reductions it
implies scalar convergence, so a bare existential proof may still hide the
original difficulty.  A satisfactory proof must construct (W_m) from a
finite approximation to (T_\ell), without querying (M_m) or an
order-(m) minimizer.

### 1.6 Finite and scalable falsifier

An asymptotic existential lemma with no rate cannot literally be disproved by
one finite instance.  A proposed constructor or a quantitative version of
Lemma L can, however, be attacked by the following finite protocol.

Fix a finite signing (B) far enough along the sequence defining (T_\ell),
a claimed accuracy (epsilon), and a candidate (W_m).

1. **Cheap necessary checks.**  Compute (V(W_m)/m^2) and a certified norm
   upper bound.  Solve the weighted MaxQP
   (max_{x\in\{\pm1\}^m}|x^{\mathsf T}W_mx|).  A spin exceeding the
   continuity-predicted threshold is an immediate, independently verifiable
   counterexample.  Exact branch-and-cut gives a finite certificate; an SDP
   supplies a scalable search and sometimes a certified bound.
2. **Finite profile check.**  Round inputs to a grid (G_h\subset[-1,1]).
   The proved discretization estimate

   \[
   d_{\rm LP}(\mathcal L(f,Sf),\mathcal L(g,Sg))
   \le h+\sqrt{Dh}
   \]

   fixes (h=\Theta(\epsilon^2/(1+D))), so
   (q=|G_h|=O((1+D)/\epsilon^2)).  The two sets

   \[
   \{\mathcal L(g,T_{W_m}g):g\in G_h^m\},\qquad
   \{\mathcal L(u,T_Bu):u\in G_h^{|B|}\}
   \]

   are finite.  Exhaustive directed containment, with the displayed
   discretization allowance, is therefore an exact finite falsifier.
3. **Scalable separation.**  Do not enumerate all (q^m) colorings.  Search
   for an adversarial grid coloring and a piecewise-linear
   (1)-bounded, (1)-Lipschitz test function separating its empirical law
   from the finite target net.  The labels and a piecewise-linear test can be
   encoded in a mixed-integer program; column/cut generation adds only the
   colorings actually found.  A returned coloring plus separating test is a
   short certificate whose empirical integrals are checked in (O(m^2)).
4. **Scaling test.**  Run the same constructor at geometrically increasing
   orders and decreasing (h).  Nondecaying lower bounds for either
   (V/m^2) or the directed profile defect falsify its claimed modulus.  They
   do not, by themselves, disprove the unquantified existential Lemma L.

This protocol is especially suited to detecting a hidden fibre coloring or
absorber coloring that fixed local-density checks miss.

### 1.7 Plausible theorem route

A plausible proof program is a separation-oracle polarization scheme.

1. At fixed (eta), replace the selected one-profile by its finite
   (q)-grid, bounded-Lipschitz description at accuracy (epsilon).
2. Build a weighted order-(m) lift of a sufficiently late finite cluster
   approximant, using randomized microcells or a colored-gadget packing.
3. Repeatedly expose a violating grid coloring and separating test, add its
   linear/bilinear constraint, and polarize edge variables while maintaining
   all accumulated constraints.
4. Prove an entropy or potential-drop lemma saying that the oracle process
   ends with total fractional variance (o(m^2)) and no new separating
   coloring.  Then invoke (1.4), not a matrix-norm theorem, for the final
   exact rounding.
5. Let (epsilon=\epsilon(m)\downarrow0) slowly and diagonalize, holding
   (eta) fixed during the order limit.

Spencer/LSV discrepancy shows that a final (o(m^2)) fractional support can
be rounded at the right scale.  It does **not** by itself provide the needed
potential-drop lemma: applying its worst-case (O(\sqrt{rm})) bound from
(r=\Theta(m^2)) starts at (O(m^{3/2})), not little-(o).  The new content
must come from slack or structure in the finite target constraints.

### 1.8 Circularity and scale traps

A proof of Lemma L must avoid all of the following.

1. Choosing (W_m) or (A_m) by minimizing (Q) at order (m).
2. Using the desired asymptotic value of (M_m) to certify the weighted
   outer profile.
3. Taking (W_m=A_m) from unknown target-order minimizers; (V=0) then makes
   the rounding condition vacuous and hides the whole theorem.
4. Storing all (2^m) values (H(x)), or an equivalent complete histogram,
   as the “finite” target state.
5. Proving concentration for each fixed coloring but omitting the universal
   outer quantifier over (q^m) colorings.
6. Reversing the directed inclusion.  For an upper bound, profiles of the
   recovered object must lie near the target; merely realizing target
   profiles does not exclude a new high-energy profile.
7. Confusing (V=o(m^2)) with (v_{\max}=o(m)), or quoting an
   (O(m^{3/2})) discrepancy theorem as a little-(o) result.
8. Using unbiased rounding at (W=0), where (V=\binom m2) and every theorem
   correctly stops at the natural, non-negligible scale.
9. Recovering only on the original possibly sparse liminf subsequence.
10. Sending (eta\downarrow0) before the fixed-(eta) order limit, thereby
    silently assuming a uniform operator bound.

## 2. Proposal II: outer-profile gadget designs plus bounded-degree repair

### 2.1 Architecture

At fixed accuracy, encode the finite grid approximation to the target
one-profile in a fixed library of colored/directed signed gadgets.  Pack the
gadgets into a typical dense host so that the signed packed core has the
required directed outer profile.  Leave a bounded-degree graph uncovered and
fill its edges arbitrarily at the end.  Accuracy then tends to zero slowly.

The all-order appeal is genuine: the near-packing theorem does not require
the final order to satisfy the exact design divisibility conditions, and the
leave can be filled because the desired final object is merely a sign matrix,
not an exact design.

### 2.2 What the known design theorem actually supplies

For a fixed graph (F), the Glock--Kuehn--Lo--Osthus typical-host packing
theorem gives, under its explicit fixed typicality hypotheses, an (F)-packing
whose leave (L_m) satisfies

\[
 \Delta(L_m)\le C(F,p),\qquad e(L_m)=O(m).                 \tag{2.1}
\]

Keevash and Glock--Kuehn--Lo--Osthus give zero-leave decompositions at all
sufficiently large admissible orders.  Colored and directed design machinery
can encode a fixed finite menu of local statistics.  None of these theorems
states the universal outer-profile conclusion.

The maximum-degree part of (2.1), not merely the (O(m)) edge count, is
important.  If a zero/missing core is completed on (L_m), the perturbation
(E) obeys

\[
 \|E\|_{\rm op}\le\max_i\sum_j|e_{ij}|\le C,              \tag{2.2}
\]

or \(2C\) if the operation flips already assigned signs.  Therefore
\(\|E/\sqrt m\|_{\rm op}=O(m^{-1/2})\), and the completion preserves all
fixed profiles as well as changing (Q) by only (O(m)).

A leave known only to have (O(m)) edges gives the scalar (O(m)) repair
bound, but the generic Frobenius estimate is only
(|E\|_{\rm op}=O(\sqrt m)).  Thus “(O(m))-leave” and
“bounded-degree (O(m))-leave” must not be conflated.

### 2.3 Verdict on universal-profile enforcement

Bounded-degree repair can **preserve** a universal outer profile already
proved for the packed core.  It cannot **enforce** that profile.  Fixed
template counts, typicality, bounded conflicts, high girth, and spread
inclusion probabilities do not control every one of the (q^m) grid
colorings.  An additional theorem is required:

> for every fixed finite target profile certificate, there is a distribution
> or a deterministic choice of gadget packing whose signed core satisfies all
> outer-profile inequalities simultaneously, while retaining a bounded-degree
> leave.

This “uniform outer design” statement is the entire missing bridge.  It could
plausibly be attacked by a high-entropy random packing/absorption process,
a finite separating-test net, and exponential concentration strong enough to
union bound over (q^m) colorings.  Existing spread bounds control inclusion
events, not the needed concentration, so they cannot be cited as this bridge.

The scalable falsifier from Section 1.6 applies directly.  In addition, color
vertices according to gadget roles, absorber membership, or decomposition
orbits.  Such colorings are invisible to uncolored local counts and are the
first likely witnesses of an outer-profile violation.

**Frozen status:** retain as a promising last-mile architecture only after a
uniform outer-design lemma is formulated and passes finite coloring tests.
The bounded-degree repair module is valid; it is not the enforcement module.

## 3. Proposal III: conference-fibre/quasirandom lifts

### 3.1 Architecture

Take a late finite signing (B_j) from a selected purified cluster.  To reach
a much larger order, replace each vertex by a fibre and each base edge by a
signed quasirandom block.  The block row sums should carry the coherent base
action at the (\sqrt r) scale, while the orthogonal fibre directions should
cancel.  Unequal fibres and the final (o(m)) vertices would be handled by
principal compression or balanced random padding.

The exact carrier theorem is strong.  A Paley symmetric conference matrix of
order (q+1), for (q\equiv1\pmod4), satisfies

\[
 C^2=qI,\qquad \|C\|=\sqrt q.                              \tag{3.1}
\]

For every sufficiently large integer (r), take the least prime
(p\equiv1\pmod4) with (p+1\ge r) and an (r\)-vertex principal
compression.  The prime number theorem in this fixed progression and the
variational bound for a principal compression give an exact symmetric hollow
sign matrix \(C_r\) with

\[
 \|C_r\|\le\sqrt p=(1+o(1))\sqrt r.                       \tag{3.2}
\]

Thus optimal-order sign carriers exist at every large fibre size.

### 3.2 Missing lemma and principal obstruction

The required new statement would be a **fibre-collapse lemma**: every bounded
coloring of the lifted matrix could be replaced, up to vanishing profile
error, by information carried by the base fibres.  Operator norm and cut
discrepancy do not imply this.  A normalized conference block has spectral
modes of order one, not vanishing modes, and an adversarial coloring may align
with or threshold an internal mode.  These modes can create profile elements
not present in the base target.

Principal compression also forfeits the exact conference identity and does
not automatically provide a prescribed regular row sum.  Neither Paley nor
Chung--Graham--Wilson quasirandomness supplies the required coherent
(\sqrt r) row sums together with collapse of every orthogonal coloring.

The decisive finite test is therefore simple: for small base signings and
available conference blocks, build the proposed fibre lift and run the grid
outer-profile separation oracle of Section 1.6.  Seed it with colorings given
by signs/thresholds of extreme fibre eigenvectors and with colorings that vary
inside only one fibre.  A persistent separated law or an energy ratio above
the base threshold kills that lift rule before any all-order bookkeeping is
attempted.

**Frozen status:** high risk.  Conference compression is an excellent
all-order low-norm comparison family and possibly a padding/filler component,
but no theorem in the toolkit makes it a target-dependent recovery family.
Develop this architecture only if the fibre-mode falsifier finds collapse in
a concrete, more strongly mixed block rule.

## 4. Final assessment

The matrix-rounding proposal isolates the smallest credible theorem target:
recover a selected cluster by almost-integral weighted matrices in directed
one-profile, with total fractional variance (o(n^2)).  Once that is done,
exact signing and the correct little-(o(n^{3/2})) error are already supplied
by scalar Bernstein rounding.  Requiring spectral purification of the
rounding residual would strengthen the hypothesis without helping the scalar
convergence proof.

Design absorption and conference matrices remain useful modules, but neither
currently supplies the universal outer quantifier.  A bounded-degree design
leave is action-stable repair, not profile enforcement; a conference carrier
is spectrally optimal, not profile selective.  The campaign should therefore
test concrete constructors against finite adversarial colorings while treating
Lemma L, rather than an all-purpose exact action-recovery theorem, as the
frozen primary target.
