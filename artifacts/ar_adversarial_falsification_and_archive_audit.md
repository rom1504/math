# Adversarial falsification and archive audit for all-order action recovery

Date: 2026-08-16.

Status: adversarial theorem audit.  The complete repository was available.
This note distinguishes a fixed-positive all-order obstruction from a
reformulation or an order effect that vanishes asymptotically.

## 1. Claims being audited

For a symmetric hollow signing (A) of order (n), put

```math
T_A=A/\sqrt n,
\qquad
\Phi(T_A)=\frac{2Q(A)}{n^{3/2}}.
```

The proved minimal structural hypothesis
\(\mathrm{AR}_{\min}^{\to}\) asks, for one selected purified liminf
cluster (T) at each of a null sequence of tolerances, for exact signings
on an upward ratio-dense set of orders such that

```math
\delta_m=\partial_1(T_{A_m},T)\to0,
\qquad
D_m\sqrt{\delta_m}\to0,
```

where (D_m) bounds both normalized operator norms.  Quantitative
one-sided action continuity and lossless principal deletion then prove
convergence.

The weighted objective condition (WR.9) asks instead for symmetric hollow
\([-1,1]\)-matrices (W_m) with

```math
\Phi(T_{W_m})\le\Phi(T)+o(1),
\qquad
V(W_m):=\sum_{i<j}(1-w_{ij}^2)=o(m^2).              \tag{1.1}
```

The directed weighted variants additionally require outer one-profile
recovery and the operator bounds needed for continuity.

## 2. Executive verdict

1. No scalable parity, spectral-multiplicity, conference-order, design
   divisibility, or Witt-class obstruction to selected
   \(\mathrm{AR}_{\min}^{\to}\) was found.  Every concrete arithmetic
   obstruction located changes only (O(1)) vertices or admits a
   ratio-dense carrier, hence has vanishing action cost.  There is presently
   no example of a signed action cluster (T) with a fixed positive recovery
   gap on infinitely many target orders.
2. Objective sign-near weighted recovery is **logically equivalent** to
   exact objective recovery on the same orders.  It is a valid and useful
   rounding module, but not a strict reduction of the remaining convergence
   obligation.
3. A sign-near weighted matrix exposes an exact sign skeleton on
   (\binom m2-o(m^2)) edges.  Without a separately specified low-complexity
   constructor it is not demonstrably less informative than an exact
   signing.
4. The new bilinear rounding argument is correct, including its constants
   and its directed one-profile consequence.  It does not control the
   rounded operator norm, but the proposed scalar convergence implication
   does not need that control.
5. The canonical profile-pressure and hard-shell entropy lemmas are
   stronger than all-order recovery, not reductions of it.  They collide
   with the archived signing-pressure/microcanonical program and inherit its
   endpoint and cross-block obstruction.
6. Sparse design repair is a genuine, separately valid last-mile lemma, and
   the maximum-degree hypothesis proposed for it is unnecessary: an
   (O(n))-edge leave is already negligible in directed one-profile and keeps
   the normalized operator norm bounded.  The universal-profile enforcement
   preceding repair remains the whole AR obstruction and is not supplied by
   design theory.
7. Ordinary conference-fibre, independent balanced, and separately paid
   mean-plus-residual lifts are already obstructed by the exact block
   Frobenius and microcanonical ANOVA identities in the archive.
8. Uniform mesoscopic induced sampling is decisively false: every
   operator-bounded sign parent is dense-quasirandom, so its fixed-size
   induced samples tend to iid signs and incur an explicit
   (0.063846\ldots) gap in \(\Phi\) above the extremal (1)-scale.

The campaign has therefore produced a useful sharpening of the boundary of
AR and the strict profile-information quotient in Section 2.1, but no
realization theorem for that quotient and no real epsilon falsifier.

### 2.1 The strongest set-valued weakening

The fixed selected cluster in \(\mathrm{AR}_{\min}^{\to}\) is not needed.
Fix one tolerance and one purified near-liminf sequence with a common
normalized operator bound (C_\eta), and let \(\mathcal K_\eta\) be its
compact set of action cluster points.  Quantitative continuity gives

```math
\sup_{T\in\mathcal K_\eta}\Phi(T)\le2L+\eta.             \tag{2.1}
```

There are two successive weakenings.  First, the target operator may vary
with the order: it would suffice that

```math
\inf_{T\in\mathcal K_\eta}\partial_1(T_{A_m},T)\to0.    \tag{2.2}
```

More strongly, the target may vary with the **individual source profile**.
Define the closed extremal profile envelope

```math
\mathcal E_\eta
:=\overline{\bigcup_{T\in\mathcal K_\eta}\mathcal S_1(T)}
```

in Levy--Prokhorov distance, and define its directed excess by

```math
\partial_{\mathcal E_\eta}(S)
:=\sup_{\mu\in\mathcal S_1(S)}
   \inf_{\nu\in\mathcal E_\eta}d_{\rm LP}(\mu,\nu).     \tag{2.3}
```

This condition does not require all profiles of (S) to coexist in any one
target operator.  Each source profile may be matched to a profile from a
different member of \(\mathcal K_\eta\), or to a limit of such profiles.

> **Proposition 2.1 (envelope continuity).**  Suppose
> \(\|S\|_{2\to2}\le D_S\), every member of
> \(\mathcal K_\eta\) has norm at most (C_\eta), and
> \(\partial_{\mathcal E_\eta}(S)\le\delta\le1\).  Put
> \(D=\max\{D_S,C_\eta\}\).  Then
>
> ```math
> \boxed{
> \Phi(S)
> \le\sup_{T\in\mathcal K_\eta}\Phi(T)
>    +5D\sqrt\delta+\delta.}                             \tag{2.4}
> ```

**Proof.**  Take any profile law \(\mu=\mathcal L(X,Y)\) of (S) and
choose \(\nu=\mathcal L(X',Y')\in\mathcal E_\eta\) at LP distance at
most \(\delta+o(1)\).  The output second moments are at most (D^2).
Strassen coupling makes the two pairs coordinatewise \(\delta+o(1)\)-close
outside an event of probability \(\delta+o(1)\).  On the good event use

```math
|XY-X'Y'|\le |Y-Y'|+|Y'||X-X'|,
```

and on the bad event use Cauchy--Schwarz.  This gives, with room in the
constant,

```math
|\mathbb E[XY]-\mathbb E[X'Y']|
\le5D\sqrt\delta+\delta+o(1).                            \tag{2.5}
```

If \(\nu\) is only in the closure, approximate it by profile laws from
members of \(\mathcal K_\eta\).  Their uniform output (L^2) bound makes
\(X'Y'\) uniformly integrable, so

```math
|\mathbb E[X'Y']|
\le\sup_{T\in\mathcal K_\eta}\Phi(T).
```

Take the supremum over \(\mu\). \(\square\)

It is therefore enough, on an upward ratio-dense set, to construct exact
signings (A_m) with

```math
\delta_m:=\partial_{\mathcal E_\eta}(T_{A_m})\to0,
\qquad
D_m\sqrt{\delta_m}\to0,                                 \tag{2.6}
```

where (D_m\ge\max\{\|T_{A_m}\|_{2\to2},C_\eta\}\).
Equations (2.1), (2.4), and principal deletion prove the usual limsup bound.
Call (2.6) **extremal-envelope recovery (EER)**.

The implication chain among structural targets is

```math
\text{one fixed target}
\Longrightarrow\text{one target varying with }m
\Longrightarrow\text{profilewise envelope recovery}
\Longrightarrow\text{scalar recovery}.                 \tag{2.7}
```

Neither reverse implication is known.  In particular, scalar convergence
does not force target-order low-objective signings to have profiles in the
envelope of a previously chosen liminf sequence.  Thus envelope recovery is
not equivalent to scalar (MR); it remains a genuine structural sufficient
condition.  It discards cluster identity and all information about which
profiles coexist in one operator.  It is the weakest action-profile target
identified in the audit short of retaining only the expectation bound, which
is exactly scalar recovery.  The outer universal quantifier over source
profiles nevertheless remains.

One can write a formally weaker condition by matching only an absolute
maximizing profile of each (A_m).  That is not an admissible information
reduction here: defining the selected profile requires solving the complete
Boolean maximum (or storing its ground-state layer), and its only used datum
is the maximal expectation.  Under the project's noncircularity criterion it
is scalar (MR) with an attached witness.  Envelope recovery is the weakest
version found that avoids selecting a target-order optimizer or ground state.

## 3. Objective WAR is equivalent to exact objective recovery

Fix an action object (T) and an upward ratio-dense set of orders
\(\mathcal N\).  Consider the following two statements.

* \(\mathrm{ER}_{\rm obj}(T)\): there are exact signings (A_m),
  (m\in\mathcal N), with
  \(\Phi(T_{A_m})\le\Phi(T)+o(1)\).
* \(\mathrm{WAR}_{\rm obj}(T)\): there are weighted matrices (W_m)
  satisfying (1.1).

> **Proposition 3.1.**  On the same order set,
>
> ```math
> \boxed{
> \mathrm{ER}_{\rm obj}(T)
> \quad\Longleftrightarrow\quad
> \mathrm{WAR}_{\rm obj}(T).}                         \tag{3.1}
> ```

**Proof.**  Exact recovery implies weighted recovery by taking (W_m=A_m),
for which (V(W_m)=0).  Conversely, scalar biased rounding gives a signing
\(A_m\) at the same order with

```math
Q(A_m)\le Q(W_m)+C\bigl(\sqrt{mV(W_m)}+m\bigr)
        =Q(W_m)+o(m^{3/2}).                              \tag{3.2}
```

The identity \(\Phi(T_W)=2Q(W)/m^{3/2}\) now gives exact recovery. \(\square\)

For a null sequence of purified clusters, the family version of
\(\mathrm{ER}_{\rm obj}\) is the scalar statement (MR) from
`minimal_all_order_action_recovery.md`.  Purification and deletion make (MR)
equivalent to convergence: convergence supplies target-order minimizers,
and (MR) supplies the missing limsup inequality.  Consequently the family
version of objective WAR is also equivalent to convergence.

This does not diminish the rounding theorem: it completely removes the
integrality **substep** once a suitable weighted construction is already in
hand.  It does mean that existence of the weighted sequence, without an
independent construction principle, cannot be counted as a strict
mathematical reduction.

The directed LV-WOR proposal contains additional structural information and
is not asserted to be logically equivalent to scalar convergence.  However,
its weighted relaxation by itself does not make that information smaller:
exact directed recovery witnesses it with (W=A), while its remaining
outer-profile condition is precisely the unresolved universal quantifier.

## 4. Sign-nearness retains almost the complete sign phase

Let (N=\binom n2).  For fixed (0<\gamma<1), define

```math
E_\gamma(W)=\{e:|w_e|\le1-\gamma\}.
```

Every edge in this set contributes at least (2\gamma-\gamma^2) to
\(V(W)\).  Hence

```math
|E_\gamma(W)|
\le\frac{V(W)}{2\gamma-\gamma^2}=o(n^2).                \tag{4.1}
```

At every fixed numerical precision resolving the margin (\gamma), (W)
therefore reveals the signs of (N-o(n^2)) individual edges.  Equivalently,
the map (W\mapsto\operatorname{sgn}W) retains an unrestricted
\(2^{N-o(n^2)}\)-sized family of possible sign skeletons.  A particular
algebraic formula may of course describe its skeleton succinctly, but that
must be proved from the constructor; it is not a consequence of
\(V=o(n^2)\).

The entropy boundary can be sharpened from the Pinsker estimate recorded in
the weighted-recovery note.  Let (\mu) be a law on the (N) edge signs,
let (U) be uniform product measure, and let (w_e=\mathbb E_\mu A_e).
Put (q_e=(1-|w_e|)/2).  Marginal entropy is (h(q_e)), and

```math
q_e\le\frac{1-w_e^2}{2}.
```

Since binary entropy is increasing on \([0,1/2]\), subadditivity of entropy
and Jensen's inequality give

```math
\begin{aligned}
H(\mu)
&\le\sum_e h(q_e)
 \le\sum_e h\!\left(\frac{1-w_e^2}{2}\right)\\
&\le N h\!\left(\frac{V(W)}{2N}\right),\\[1mm]
D(\mu\|U)
&\ge N\left[\log2-h\!\left(\frac{V(W)}{2N}\right)\right]. \tag{4.2}
\end{aligned}
```

Thus (V(W)=o(N)) forces

```math
\boxed{D(\mu\|U)\ge N\log2-o(N),\qquad H(\mu)=o(N).}    \tag{4.3}
```

A microcanonical law whose barycenter is sign-near pays asymptotically the
maximal possible relative entropy against iid signs, not merely a positive
quadratic-rate cost.  This does not obstruct a deterministic construction;
it does rule out any claim that a low-entropy-cost tilt of iid signs will
produce WAR through its barycenter.

## 5. Independent audit of bilinear/profile rounding

The claim in `ar_design_rounding_independent_proposals.md` is correct.  Let
\(E=A-W\) under independent biased rounding and put

```math
B(E)=\max_{x,y\in\{\pm1\}^m}|x^{\mathsf T}Ey|.
```

For fixed (x,y), expansion over independent upper-triangular edges gives

```math
x^{\mathsf T}Ey
=\sum_{i<j}E_{ij}(x_i y_j+x_jy_i).                       \tag{5.1}
```

Each summand is bounded by four and the total variance is at most (4V(W)).
Bernstein and a union bound over (4^m) pairs give

```math
\Pr\{B(E)\ge t\}
\le2\,4^m
\exp\!\left[-\frac{t^2}{2(4V(W)+4t/3)}\right].          \tag{5.2}
```

With (a_m=2m\log2+\log4) and
\(t_m=\sqrt{8V(W)a_m}+(8/3)a_m\), the right side is at most (1/2).
Thus some supported outcome has (B(E)=O(\sqrt{mV(W)}+m)).  All factors of
two are consistent with the half-quadratic definition of (Q).

Duality and extremality of the cube give the exact identity

```math
\sup_{|f|\le1}\|Ef\|_1=B(E).                            \tag{5.3}
```

On the uniform probability space,

```math
\sup_{|f|\le1}\|T_Af-T_Wf\|_{L^1}
=\frac{B(E)}{m^{3/2}}.                                  \tag{5.4}
```

The same-input coupling and Markov's inequality therefore put every
one-profile of (T_A) within
\(\sqrt{B(E)/m^{3/2}}=o(1)\) of the corresponding one-profile of (T_W).

This estimate does **not** imply
\(\|A-W\|_{op}=o(\sqrt m)\), nor a bounded normalized operator norm for
the rounded matrix.  The scalar argument remains valid because it applies
action continuity to (W_m\) and (T) before rounding, then uses (5.2)
directly on (Q).  It should not be cited as a proof of the quantitative
operator-bound clause in \(\mathrm{AR}_{\min}^{\to}\).

## 6. Canonical pressure and shell entropy add a stronger obligation

The proposed profile mismatch Hamiltonian is

```math
R_{T,D}(A)=\partial_1(T_A,T)
 +\min\{1,(\|T_A\|_{2\to2}-D)_+\},                      \tag{6.1}
```

with minimum (m_n) and normalized pressure

```math
F_n(\beta)=-\frac1{\beta N_n}
\log\left(2^{-N_n}\sum_A e^{-\beta N_nR_{T,D}(A)}\right). \tag{6.2}
```

The Laplace sandwich

```math
m_n\le F_n(\beta)\le m_n+\frac{\log2}{\beta}             \tag{6.3}
```

is correct, and convergence of (F_n(\beta)) for every fixed integer
\(\beta\), together with a subsequence on which (m_n\to0), would imply
all-order recovery by sending (\beta\to\infty) after (n\to\infty).

However pressure convergence is not a weaker statement than recovery.
Here is an exact abstract falsifier to that logical interpretation.  Let a
finite landscape have (2^{N_n}) states and minimum zero at **every** order.
At even (n), give (2^n) states energy zero and all other states energy
one.  At odd (n), give half the states energy zero and half energy one.
The zero set may be grouped into blocks of size (2^n), matching the
switching-orbit quantization for all sufficiently large (n).  Then
all-order recovery is already complete, while

```math
\lim_{\substack{n\to\infty\\n\ \mathrm{even}}}F_n(\beta)
=\min\left\{1,\frac{\log2}{\beta}\right\},
\qquad
\lim_{\substack{n\to\infty\\n\ \mathrm{odd}}}F_n(\beta)=0. \tag{6.4}
```

Thus (F_n(\beta)) can fail to converge for every fixed (\beta>0) even
when a zero-energy state exists at every order.  Likewise, at any shell
radius below one, normalized shell entropy alternates between zero and
\(\log2\).  The hard-shell entropy limit is also strictly stronger than
nonemptiness at all orders.

This is not a counterexample inside the particular matrix Hamiltonian; it is
a proof that the pressure/shell assertions add an independent entropy-density
regularity obligation.  No cited theorem supplies that obligation at the
\(A/\sqrt n\) fluctuation normalization.

The archive collision is exact rather than terminological:

* `good_signing_entropy_threshold.md` already proves canonical and refined
  microcanonical convergence criteria on edge-signing space, then isolates
  the support-edge/cavity obstruction;
* `microcanonical_disorder_counting_composition.md` proves an exact
  lower-tail product theorem, but its parent composition contracts the
  project temperature and misses the endpoint by a leading cross term; and
* the block identity in `good_signing_entropy_threshold.md` shows that
  cross edges cannot cancel internal energy after the relevant sign flip.

Replacing the objective by directed profile mismatch removes useful algebraic
structure while retaining the supremum over all colorings.  It does not evade
those endpoint issues.  The canonical and hard-shell proposals are therefore
class C: clean sufficient reformulations, but not strict reductions and not
currently backed by an interpolation or LDP at the required scale.

## 7. Design repair must be separated from profile enforcement

The design proposal contains one valid theorem-level module, in a slightly
stronger form than stated there.  Suppose a packed signed core already has
the desired outer profile and its leave has (r=O(n)) unordered edges.
Arbitrarily signing missing edges or flipping assigned edges produces a
perturbation (E) with entries of magnitude at most two and

```math
\|E\|_F\le2\sqrt{2r}=O(\sqrt n).                         \tag{7.1}
```

Thus the normalized operator perturbation is (O(1)), so a uniformly bounded
core remains uniformly bounded.  More importantly, directly for every
\(|f|\le1\),

```math
\|T_Ef\|_{L^1}
\le\frac{1}{n\sqrt n}\sum_{i,j}|e_{ij}|
\le\frac{4r}{n^{3/2}}=O(n^{-1/2}).                       \tag{7.2}
```

The same-input coupling therefore gives directed one-profile error
\(O(n^{-1/4})\) by Markov, while the edge-edit inequality changes (Q) by
only (O(n)).  Consequently every (O(n))-edge repair, not only a
bounded-degree one, is action-stable and objective-negligible.  Bounded
maximum degree improves (7.1) to (O(1)) before normalization, but that extra
strength is unnecessary for the recovery implication.  Sparse repair is
grade A as a last-mile statement.

It does not enforce the profile of the packed core.  At fixed accuracy the
profile admits a finite alphabet and a finite net of output laws, but the
outer condition still says that **every** one of (q^n) vertex colorings
must land near that target net.  Keevash-type decomposition, nibble,
absorption, and typical-host packing theorems enforce a fixed menu of local
densities; none supplies this universal local-global coloring assertion.
The proposed “uniform outer design” theorem is therefore the entire missing
AR bridge.  Without a new separation-oracle potential or a same-coloring
cancellation theorem, it is class C, not a design-theoretic reduction.

This correction matters for archive classification: typical-host packing
with an (O(n))-edge leave already solves the repair part.  It does not solve
the enforcement part.

## 8. Blow-up and conference-fibre archive comparison

For an (r\times r) sign block whose coarse row sum is
\((1+o(1))\sqrt r\), the orthogonal mean/residual decomposition gives

```math
\|R\|_F^2=r^2-(1+o(1))r=(1-o(1))r^2.                    \tag{8.1}
```

Thus almost all microscopic mass remains in fibre modes.  This is
`bounded_op_signed_realization.md`, not a heuristic.  The stronger exact
microcanonical ANOVA identity in `regular_microblock_absorption_audit.md`
shows that, for arbitrary fixed sign blocks and fixed coarse
magnetizations, the macro, one-fibre, and two-fibre terms contribute as a
sum of nonnegative squares.  A static row/column “Onsager” correction cannot
cancel the macro term; moving Frobenius mass out of the two-fibre residual
only creates nonnegative one-fibre variance.

Consequently:

* independent biased residuals are killed by the leading centered supremum
  in `random_biased_lift_no_go_phase2.md`;
* regular Hadamard/conference-type orthogonal blocks are killed by their
  surviving fibre spectral modes and the ANOVA identity;
* separately paying macro and residual channels has a fixed leading loss;
  and
* an unspecified “joint residual absorber” or “fibre-collapse lemma” is just
  directed recovery for that lift family until an independently defined
  nonlinear mechanism is supplied.

The all-order existence of Paley conference **carriers** is not enough.  It
controls operator norm but does not make their nonconstant fibre modes vanish
from the action profile.  The conference-fibre proposal is class D for the
ordinary regular/independent implementation and class C only if “absorber” is
left as an unspecified new theorem.

## 9. The independent action proposals

### 9.1 Mesoscopic induced sampling is rigorously false

The proposed lemma \(L_{\rm samp}\) is killed by a universal
quasirandomness consequence of the operator bound.

> **Theorem 9.1 (mesoscopic iid universality).**  Let (A_N) be any
> symmetric hollow sign matrices satisfying
> \(\|A_N\|_{op}\le C\sqrt N\).  For fixed (m), choose a uniform ordered
> (m)-tuple of distinct vertices and let (B_{N,m}) be the induced labeled
> signing.  Then
>
> ```math
> \boxed{
> d_{\rm TV}\bigl(\mathcal L(B_{N,m}),
>                  \mathcal L(R_m)\bigr)\longrightarrow0
> \quad(N\to\infty),}                                  \tag{9.1}
> ```
>
> where (R_m) has independent Rademacher upper-triangular entries.

**Proof.**  Form the ordinary graph with off-diagonal adjacency

```math
G_N=\frac{J-I+A_N}{2}.
```

For all vertex sets (S,T),

```math
\left|\sum_{i\in S,j\in T}
 \left((G_N)_{ij}-\frac12\right)\right|
\le\frac12\|A_N\|_{op}\sqrt{|S||T|}+O(N)
=O(N^{3/2}).                                            \tag{9.2}
```

Hence the dense cut distance from (G_N) to the constant-(1/2) graphon is
\(O(N^{-1/2})\).  The dense counting lemma (equivalently, the standard
quasirandom-graph equivalence) implies that every fixed labeled induced graph
on (m) vertices has density

```math
2^{-\binom m2}+o_N(1).                                  \tag{9.3}
```

There are finitely many such patterns for fixed (m), so (9.3) is exactly
the total-variation assertion.  Sampling an unordered subset and randomly
labeling it is equivalent. \(\square\)

This theorem says that the iterated sampling limit forgets **every** bounded
signed action object, not only conference structure.  Its zero-temperature
effect has a quantitative gap.  For iid (R_m), expose vertices in order and
set

```math
x_1=1,
\qquad
x_j=\operatorname{sign}\sum_{i<j}(R_m)_{ij}x_i.
```

The resulting energy is a sum of independent simple-random-walk absolute
values.  Its mean divided by (m^{3/2}) tends to

```math
c_{\rm iid}=\frac23\sqrt{\frac2\pi}
=0.5319230405\ldots,                                    \tag{9.4}
```

and its variance is (O(m^2)).  Chebyshev therefore gives

```math
\frac{Q(R_m)}{m^{3/2}}
\ge c_{\rm iid}-o_{\mathbb P}(1),
\qquad
\Phi(T_{R_m})
\ge\frac43\sqrt{\frac2\pi}-o_{\mathbb P}(1).          \tag{9.5}
```

Put

```math
\gamma:=\frac43\sqrt{\frac2\pi}-1
=0.0638460810\ldots.                                    \tag{9.6}
```

Now let (T) be any bounded action limit with \(\Phi(T)\le1\), as occurs
for the relevant null sequence of purified extremal clusters because
\(L\le1/2\).  If \(L_{\rm samp}\) held with sample norm bound (D), choose
fixed \(\delta>0\) so small that

```math
5\max\{D,\|T\|_{2\to2}\}\sqrt\delta+\delta<\gamma/3. \tag{9.7}
```

With probability tending to one in the iterated limit, the sample would
have norm at most (D\sqrt m) and directed profile error at most (\delta).
One-sided continuity would force its \(\Phi\) to be at most
\(1+\gamma/3\).  Equations (9.1) and (9.5) force it to be at least
\(1+2\gamma/3\) with probability tending to one.  This is a contradiction.

Thus

```math
\boxed{
L_{\rm samp}\text{ is false for every bounded }T
\text{ with }\Phi(T)\le1.}                              \tag{9.8}
```

For clusters with \(\Phi(T)\le1+\eta\), the same proof applies whenever
\(\eta<\gamma\), which is enough for the null tolerance sequence in the
convergence program.  The fourth-moment conference calculation in the
proposal is consistent with this stronger no-go but is no longer needed.
Mesoscopic uniform sampling is grade D.

The order of limits matters.  The theorem kills the exact iterated statement
\(L_{\rm samp}\) and every diagonal lying inside a regime where induced
subgraphs have already mixed to iid.  It does not kill principal sampling at
relative size (1-o(1)), which is the proved near-order theorem, nor does it
by itself rule out a carefully tuned diagonal (m=m(N)\) above the parent's
quasirandom mixing scale.  Such a diagonal result would have to use
quantitative structure beyond the operator bound and the limiting action
object; it cannot follow from the proposed uniform mesoscopic invariance.

### 9.2 Exposed-face insertion is an exact archive collision

The identities

```math
Q(A^b)=\max_x\bigl(|H_A(x)|+|b\cdot x|\bigr),
\qquad
\min_bQ(A^b)=Q(A)+\Delta(A)                              \tag{9.9}
```

are correct.  But they, the slack-weighted constraint over every spin, the
near-cap code, and the exact derivative coefficient are already developed in
`cap_discrepancy_insertion.md` and `near_cap_insertion.md`.  This is exactly
ledger Section 10.44, equation (10.254), not merely a similar identity.  In
particular,
the desired leading increment is precisely

```math
\frac32 q\sqrt n,
```

and the archive proves that exact ground-state balancing does not control
the thick cap.  `scale_transfer_profile_no_go.md` also gives finite matrices
with the same coarse insertion profile but different optimal extensions.

The proposed chain with error (r_n\sqrt n), (r_n\to0), would indeed be
sufficient: weighted Cesaro summation makes
\(\sum_{k\le n}r_k\sqrt k=o(n^{3/2})\).  The missing lemma, however, is the
archived sharp star-insertion obligation, not a consequence of action
recovery or a new imported theorem.  Its definition still minimizes a row
against the full slack landscape.  It is class C: a valid exact sufficient
route, but not a new AR reduction.

### 9.3 Microcanonical no-gap is the same stronger entropy obligation

The proposed (L_{\rm ent}\) forces nonempty directed-profile shells at
every large order by requiring a normalized shell-entropy no-gap inequality.
This is quantitatively stronger than all-order nonemptiness.  It is the same
architecture audited in Section 6 and collides with
`good_signing_entropy_threshold.md` and
`microcanonical_disorder_counting_composition.md`.  Without an independently
proved rate/interpolation formula, it is class C.

## 10. Search for a real all-order obstruction

The required falsifier would be a signed action cluster (T) and
\(\epsilon>0\) such that for infinitely many prescribed orders every exact
signing sufficiently close to (T) has
\(\Phi(T_A)\ge\Phi(T)+\epsilon\).  No such example was found.

The following apparent obstructions do not qualify.

* **Parity and fixed congruence.**  They change (O(1)) vertices.  Principal
  compression and random balanced completion have (o(n^{3/2})) objective
  cost and vanishing profile cost under a common operator bound.
* **Conference/Paley orders.**  Primes in a fixed progression give larger
  Paley carrier orders (n+o(n)).  Principal compression removes their
  arithmetic restriction.  This does not recover an arbitrary cluster, but
  it prevents conference divisibility itself from being an epsilon
  obstruction.
* **Spectral multiplicities.**  Rounding a fixed finite collection of
  multiplicity proportions changes empirical spectral/rooted laws by
  (O(1/n)).  Exact even multiplicity cannot by itself produce a fixed
  action gap.
* **Fixed design divisibility and Witt residues.**  All located constraints
  either have bounded modulus or admit all sufficiently large admissible
  orders with bounded/vanishing repair.  No invariant continuous in the
  directed one-profile topology was found that separates two infinite order
  classes by a positive amount.
* **The zero operator.**  It is indeed separated from every exact signing in
  directed one-profile distance, but it is not a signed action cluster.  The
  Paley--Zygmund row-output argument therefore shows non-universal density,
  not failure of AR for the selected class.

Projective exchangeability is a real architectural obstruction but not a
counterexample to AR: spectral tightness forces iid signs, whose greedy
energy constant is at least (0.531923\ldots>1/2).  Order-dependent
microcanonical or deterministic realizers remain outside that theorem.

The absence of an epsilon example is substantive.  Any such example for a
purified liminf cluster would amount to a structural order-class separation
at the same scale as the original convergence problem.  Finite residue
effects cannot substitute for it.

## 11. Proposal grades and stopping judgment

| Proposal | Grade | Adversarial judgment |
|---|---:|---|
| Extremal-envelope recovery (EER) on ratio-dense orders | B | Strongest proved sufficient information quotient: profiles may match different low-objective cluster phases; realization remains open. |
| Fixed-cluster directed recovery | B/C | Stronger than necessary; subsumed by extremal-envelope recovery. |
| Objective sign-near weighted recovery | C as reduction; A as rounding module | Equivalent to exact objective recovery by Proposition 3.1. |
| Directed LV-WOR | B/C | Structurally phrased, but sign-nearness retains almost a full edge phase and universal outer enforcement remains. |
| Canonical profile-pressure limit | C | Stronger entropy-regularity obligation; archived pressure collision. |
| Hard microcanonical shell limit | C | Stronger than shell nonemptiness; archived microcanonical collision. |
| (O(n))-edge design repair | A module | Valid negligible scalar and directed-profile repair; bounded degree is not required. |
| Uniform outer-profile design | C | The entire missing universal-coloring bridge; no imported theorem proves it. |
| Ordinary conference/regular fibre lift | D | Exact residual and ANOVA obstructions apply. |
| Unspecified joint residual absorber | C | Renames directed recovery until an independently checkable mechanism is stated. |
| Projective exchangeable sampler | D | Rigorous iid/greedy obstruction. |
| Mesoscopic uniform induced sampling | D | Operator-bounded parents have iid fixed-size induced limits, causing the explicit 0.063846 Phi gap. |
| Exposed-face one-vertex absorption | C | Exact sharp-insertion route already archived; full slack landscape remains. |
| Parity, multiplicity, or design-divisibility falsifier | D | Effects vanish; no fixed epsilon. |

At this checkpoint AR remains a coherent possible architecture.  EER is a
strict information quotient of fixed-cluster AR and a sharper sufficient
target, but no proposal realizes it; hence the original convergence problem
is not closed.  The strongest positive module is sign-near rounding.  The strongest
new negative theorem is the mesoscopic iid-sampling no-go (9.8), supplemented
by the essentially maximal entropy cost (4.3).  Further execution should
proceed only if an
independent proposal supplies either

1. a constructor for the outer profile whose state is provably smaller than
   an (N-o(N))-bit sign skeleton; or
2. a continuous order-class invariant yielding a real fixed-positive AR
   obstruction.

Absent one of those, further relabeling of recovery as pressure,
microcanonical entropy, weighted polarization, or design enforcement should
stop.
