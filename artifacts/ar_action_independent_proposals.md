# Independent proposals for all-order action recovery

Date: 2026-08-16.

## Verdict

There are two concrete geometric ways to try to turn one purified liminf
representative into all-order realizers, and one genuinely different
statistical-mechanics route:

1. **top down:** take mesoscopic induced submatrices of a much larger good
   representative;
2. **bottom up:** add one vertex at a time by discrepancy-coloring only the
   exposed near-maximizer face; and
3. **order by order:** prove a no-gap theorem for the microcanonical entropy
   of a directed action-profile neighborhood.

None is currently supplied by the cited literature.  Each has a finite,
testable failure mode.  In particular, ordinary graphon sampling, ordinary
Gamma recovery, and a projectively exchangeable infinite signing do not
establish any of them.

The sign-near weighted target is **not presently a construction mechanism**.
It is an excellent terminal rounding interface, but its missing input still
contains the outer universal Boolean/profile quantifier.  Existentially, its
objective form is asymptotically equivalent to exact-sign scalar recovery.

Throughout, \(T_A=A/\sqrt n\),

```math
\Phi(T_A)=\frac{2Q(A)}{n^{3/2}},
\qquad q(T):=\frac12\Phi(T).
```

For a selected purified cluster \(T\), one has
\(q(T)\le L+o(1)\), where \(L=\liminf M_n/n^{3/2}\).  Thus any mechanism
below which supplies exact signs \(B_m\) with a common normalized operator
bound and \(\partial_1(T_{B_m},T)\to0\) gives

```math
\limsup_m\frac{M_m}{m^{3/2}}
\le \limsup_m\frac{Q(B_m)}{m^{3/2}}
\le q(T)
```

by the one-sided action-continuity estimate.  Ratio-dense orders suffice by
lossless principal deletion; all three proposals aim at every large order.

## 1. Mesoscopic induced-submatrix invariance

This is the most literal rescaled-graph-limit mechanism.  It uses no new edge
optimization: obtain an order-\(m\) signing by uniformly sampling \(m\)
vertices from a much larger already-good signing, and normalize the induced
matrix by \(\sqrt m\), not by the parent scale.

### Exact missing theorem \(L_{\rm samp}\)

Let \(A_{N_j}\) be one bounded purified sequence with
\(T_{A_{N_j}}\to T\).  If \(S_{j,m}\) is a uniform \(m\)-subset of
\([N_j]\), require a finite \(D\) such that, for every \(\epsilon>0\),

```math
\lim_{m\to\infty}\ \limsup_{j\to\infty}
\Pr\!\left(
 \partial_1\!\left(T_{A_{N_j}[S_{j,m}]},T\right)>\epsilon
 \text{or}
 \|A_{N_j}[S_{j,m}]\|_{op}>D\sqrt m
\right)=0.                                                \tag{S}
```

The iterated limit is intentional: the sample fraction may tend to zero.
A diagonal choice of \(j=j(m)\), followed by selection of one supported
sample, gives exact hollow signings at every order with directed profile
error tending to zero.  Hence (S) implies convergence by the common
implication above.

### Why this is a strict reduction

The data are one selected liminf sequence and a uniform vertex sample.  No
target-order minimizer is selected, and no target-order edge is optimized.
Only the one-sided one-profile quotient must be sampling-stable; labels,
reverse inclusion, and all joint profiles are discarded.  This is therefore
strictly less information than the maps
\(x\mapsto H_A(x)\) for all target-order signings.  It may nonetheless be a
hard theorem.

### Finite/structural falsifiers

For the exact directed statement, a direct finite test is

```math
F_{j,m}:=\mathbb E_S\Phi\!\left(T_{A_{N_j}[S]}\right).
```

If, with uniformly bounded normalized operator norms,
\(\liminf_m\liminf_jF_{j,m}>\Phi(T)+c\), then (S) is false by one-sided
continuity.

There is also a sharp warning against strengthening (S) casually to full
operator/action sampling.  If a symmetric conference matrix satisfies
\(C_N^2=(N-1)I\), and \(B=C_N[S]\) for a uniform \(m\)-subset, direct closed
walk counting gives

```math
\mathbb E\operatorname{tr}B^4
=2m(m-1)^2-m(m-1)
-\frac{m(m-1)(m-2)(m-3)}{N-3}.                            \tag{S.1}
```

Thus \(N^{-3}\operatorname{tr}C_N^4\to1\), whereas in the mesoscopic
regime \(m/N\to0\),
\(m^{-3}\mathbb E\operatorname{tr}B^4\to2\).  Global algebraic
cancellation need not survive rescaled induced sampling.  The selected
extremal cluster would have to possess extra sampling self-similarity; action
convergence alone does not provide it.

### Circularity and source boundaries

- Choosing the subset by minimizing its Boolean objective is target-order
  optimization in disguise; (S) must follow from concentration under the
  uniform sample.
- Bounded-graphon sampling cannot be inserted here.  The step kernel is
  \(\sqrt N A_N\); the bound in Theorem 4.7 of
  [Borgs--Chayes--Lovasz--Sos--Vesztergombi](https://arxiv.org/abs/math/0702004)
  grows like \(\sqrt{N/\log m}\).  The fixed-\(L^p\) results of
  [Fekete--Kunszenti-Kovacs, Theorems 2.1, 2.3 and Corollary 2.5](https://arxiv.org/abs/2203.07581)
  have the same normalization mismatch.
- The normalization error \(A[S]/\sqrt N\) versus \(A[S]/\sqrt m\) is the
  whole problem, not a harmless rescaling.
- Section 11 of
  [Backhausz--Szegedy](https://arxiv.org/abs/1811.00626) gives deterministic
  full-sign representatives at each order (Lemma 11.2), but Proposition
  11.1 still passes to subsequences and explicitly does not prove an
  all-natural-orders limit.

## 2. Exposed-face one-vertex absorption

This mechanism reduces the next-order problem to a discrepancy statement
about near-maximizing spins of one already-good matrix.  It has an exact
finite identity, so failure cannot be hidden in limit terminology.

For an order-\(n\) signing \(A\), put

```math
s_A(x):=Q(A)-|H_A(x)|,
\qquad
\Delta(A):=\min_{b\in\{\pm1\}^n}
 \max_{x\in\{\pm1\}^n}\bigl(|b\cdot x|-s_A(x)\bigr).      \tag{E.1}
```

If \(A^b\) is obtained by adjoining one vertex with incident signs \(b\),
then

```math
Q(A^b)=\max_x\bigl(|H_A(x)|+|b\cdot x|\bigr),
\qquad
\min_bQ(A^b)=Q(A)+\Delta(A).                              \tag{E.2}
```

The first equality is just
\(\max_{y=\pm1}|H_A(x)+y\,b\cdot x|=|H_A(x)|+|b\cdot x|\).

### Exact missing theorem \(L_{\rm face}\)

For every selected purified cluster \(T\), every \(\epsilon>0\), and
\(q=q(T)\), some sufficiently large seed from its realizing sequence starts
a nested extension chain \(A_n\) for which

```math
\Delta(A_n)
\le \frac32(q+\epsilon)\sqrt n+r_n\sqrt n,
\qquad r_n\longrightarrow0.                              \tag{E}
```

By (E.2) and
\(\sum_{k<n}\sqrt k=(2/3)n^{3/2}+o(n^{3/2})\),

```math
Q(A_n)\le(q+\epsilon)n^{3/2}+o(n^{3/2}).
```

Starting the chain arbitrarily far out removes the seed error.  Sending the
purification tolerance and then \(\epsilon\) to zero proves convergence.
No all-order action convergence is needed after the seed.

### Why this is a strict reduction

At step \(n\), only \(n\) new edge signs are chosen.  More importantly,
states far below the maximum receive the slack \(s_A(x)\) and impose no
effective constraint.  The relevant object is the exposed near-maximizer
face of the one-profile energy functional, not the full Boolean landscape
and not a fresh minimization over \(\binom{n+1}{2}\) edges.  In action
language, (E) asks for a discrepancy bound on the exposed face

```math
\left\{\mu\in\mathcal S_1(T):
 \left|\int uv\,d\mu\right|=\Phi(T)\right\},
```

which is a much smaller target whenever that face has controlled metric
complexity.

### Finite/structural falsifier

For \(t\ge0\), define

```math
\mathcal E_t(A)=\{x:s_A(x)\le t\sqrt n\},
\qquad
\kappa_t(A)=\frac1{\sqrt n}
 \min_b\max_{x\in\mathcal E_t(A)}|b\cdot x|.
```

Then, exactly,

```math
\frac{\Delta(A)}{\sqrt n}\ge\kappa_t(A)-t.               \tag{E.3}
```

Consequently a sequence of reachable matrices for which
\(\kappa_t(A)>\frac32(q+\epsilon)+t+c\) falsifies (E).  This is a finite
covering/discrepancy computation on the near-ground-state code.  A highly
covering exposed face is the structural obstruction.

### Circularity and Gamma hazards

- Taking \(A_{n+1}\) to be an order-\(n+1\) minimizer proves nothing; the
  new row must be obtained from an independent discrepancy theorem for the
  exposed face.
- Action convergence of \(A_{N_j}\) alone does not give convergence of
  near-maximizer sets.  One needs an upper, equicoercive max-Gamma statement
  (hypoconvergence of the exposed faces), including uniform integrability of
  the output-energy moment.
- Theorems 4.3--4.4 of
  [Braides--Cermelli--Dovetta](https://www.numdam.org/item/10.1051/cocv/2019029.pdf)
  recover a prescribed vertex state after a cut-convergent kernel sequence
  is supplied.  They do not control **every** near-maximizer, choose the new
  edge row, or recover the kernels.  Remark 4.2 explicitly rules out weak
  \(L^1\) convergence as a substitute.
- Replacing \(\Delta(A)\) by its definition and then brute-force minimizing
  it merely repackages the universal quantifier.  A proof must bound it from
  geometry (entropy, Gaussian width, or vector discrepancy) of
  \(\mathcal E_t(A)\).

## 3. Nonprojective microcanonical no-gap entropy

This is the only proposal here that does not couple different orders by
submatrix or extension.  It seeks an order-dependent, globally conditioned
law at each size.

For \(D<\infty\), let

```math
\mathcal M_m^\epsilon(T,D)=
\left\{A:\ A\text{ is a symmetric hollow order-}m\text{ signing},
\ \|T_A\|_{2\to2}\le D,
\ \partial_1(T_A,T)<\epsilon\right\},
```

and set \(h_m^\epsilon=m^{-2}\log|\mathcal M_m^\epsilon|\), with
\(\log0=-\infty\).

### Exact missing theorem \(L_{\rm ent}\)

For every bounded subsequential sign limit \(T\), there is a finite \(D\)
such that for every \(\epsilon>0\) there is
\(0<\epsilon'<\epsilon\) with

```math
\liminf_{m\to\infty}h_m^\epsilon(T,D)
\ge
\limsup_{j\to\infty}h_{N_j}^{\epsilon'}(T,D)-o_\epsilon(1),              \tag{M}
```

where \(N_j\) is any bounded realizing subsequence and
\(o_\epsilon(1)\to0\) as \(\epsilon\downarrow0\).  The right side is
finite (indeed nonnegative once a microstate exists).  Since an empty set
has entropy \(-\infty\), (M) forces \(\mathcal M_m^\epsilon(T,D)\ne\varnothing\)
for every sufficiently large \(m\).  A diagonal choice in \(\epsilon\)
then gives all-order directed recovery and hence convergence.

The intended proof mechanism is an independently defined pressure/rate
functional plus interpolation or approximate amalgamation showing that its
finite-rate domain has no order gaps.  Merely declaring the left side to be
the rate function is not a mechanism.

### Why this is a strict reduction

Only one selected cluster and its one-sided one-profile neighborhood enter;
there is no objective threshold and no target-order optimizer.  The finite
state at fixed accuracy is the profile quotient, not the labelled table
\(x\mapsto H_A(x)\).  The entropy conclusion is quantitatively stronger
than bare existence, but the *information being transported* is strictly
coarser than Boolean optimization.  This distinction is essential: without
an independent pressure/interpolation formula, (M) would just rename AR.

### Finite/structural falsifier

For fixed \((T,D,\epsilon)\), enumerate or certify
\(\mathcal M_m^\epsilon\).  Infinitely many empty sizes together with a
nonempty subsequence falsify (M).  Structurally, any macroscopic divisibility
invariant that is open in the directed profile topology creates such a gap.
The failure of universal finite approximability in
[Kun--Thom, Theorem 1.3](https://arxiv.org/abs/1901.03963) warns that no
version of (M) can hold for arbitrary graphops/PMP actions; it must use the
special exact-sign cluster and the selected profile.

### Circularity and exchangeability hazards

- The uniform law on \(\mathcal M_m^\epsilon\) cannot be used to prove that
  \(\mathcal M_m^\epsilon\) is nonempty.
- Conditioning on \(Q(A)\le(q+o(1))m^{3/2}\), using an order-\(m\)
  minimizer, or defining the rate through the desired all-order liminf is
  circular.
- The measures must be genuinely order-dependent.  A single projectively
  consistent exchangeable signing with tight normalized operator norm is
  forced by the Aldous--Hoover representation to be iid signs; see
  [Diaconis--Janson, Sections 5--6](https://arxiv.org/abs/0712.2749).
  The greedy exposure bound then gives normalized objective at least
  \((2/3)\sqrt{2/\pi}=0.5319\ldots\), above the available \(1/2\) upper
  scale.  Exchangeability cannot supply the no-gap theorem for free.
- A low-cost tilt around iid signs is also implausible: a globally
  sign-near barycenter pays order-\(m^2\) relative entropy.  A valid
  microcanonical construction must allow quadratic conditioning cost and
  strong global edge correlations.

## The sign-near weighted target: rounding interface, not recovery

For \(W\in[-1,1]^{m\times m}\), put
\(V(W)=\sum_{i<j}(1-w_{ij}^2)\).  Biased independent rounding gives an exact
signing \(A\) with

```math
Q(A)\le Q(W)+C\bigl(\sqrt{mV(W)}+m\bigr).
```

Thus \(V(W)=o(m^2)\) is exactly enough for scalar recovery.  After deleting
\(o(m)\) high-variance rows, the spectral estimate based on
[Bandeira--van Handel, Corollary 3.6](https://arxiv.org/abs/1408.6185)
also gives \(\|A-W\|_{op}=o(\sqrt m)\), hence profile recovery.

What is not constructed is \(W\).  The proposed weighted input already asks

```math
\Phi(T_W)\le\Phi(T)+o(1),
\qquad V(W)=o(m^2),                                       \tag{W}
```

on dense orders.  The first line is still a maximum over all Boolean states
(separate affinity moves extrema to cube vertices), while the second says
that all but \(o(m^2)\) total fractional variance is already resolved toward
signs.  Exact sign realizers are a special case \(V=0\), and the rounding
theorem gives the converse at the objective scale.  Hence existence of (W)
is not a strict existential reduction of scalar recovery; it moves the
universal quantifier to a nearly integral weighted matrix.

No cited Gamma or sampling result supplies (W):

- [Braides--Cermelli--Dovetta, Theorems 4.3--4.4](https://www.numdam.org/item/10.1051/cocv/2019029.pdf)
  and [Zhang--Scott--Du--Porter, Theorem 6.1](https://arxiv.org/abs/2408.00422)
  recover vertex states for a kernel sequence already assumed to converge;
- [Le--Jegelka, Theorems 2 and 4](https://arxiv.org/abs/2306.04495) obtain
  every-order arbitrary weighted discretizations only for regular
  \(L^2\)-operators, with no near-sign saturation; and
- for the present scaling the step-kernel \(L^p\) norms grow as \(\sqrt m\),
  so bounded or fixed-\(L^p\) graphon sampling treats the fluctuation as an
  error even though it is the target.

A real weighted mechanism would need a new **saturated Gamma-limsup theorem**:
from independently defined enriched action-limit data (including the
microscopic variance/noise component), construct every-order \(W_m\) with
directed profile convergence and \(V(W_m)=o(m^2)\), without minimizing the
target-order Boolean functional.  No primary source in the toolkit states
such a theorem.  Until that enrichment and construction are given, the
weighted route should be regarded only as the last rounding step.

## Ranking

The induced-submatrix theorem (S) is the cleanest unbounded-limit question
and has a canonical construction if true, but conference-type fourth-moment
loss shows that it needs genuine mesoscopic self-similarity.  The extension
theorem (E) is the sharpest finite reduction: its sole obstruction is the
weighted discrepancy of the exposed near-ground-state code.  The entropy
theorem (M) is the broadest and naturally avoids projective exchangeability,
but it qualifies as a mechanism only if its rate/pressure formula is derived
independently of the desired microstate counts.  These are genuinely
different pressure points; none is a consequence of existing action
compactness, graphon sampling, or Gamma recovery.
