# Growing-state and realization audit

Date: 2026-07-31. This is an agent-authored research report. It does not
modify user directives or the strategic files.

## 1. Starting point and classification

The exact order-ten collision and balanced twin blowups in
`phase2b_phi6_collision_report.md` rigorously falsify universal control of
cap by the fixed six-vertex switching profile: equal profiles can have a
`Theta(N^2)` cap gap. More strongly, the certified Sylvester--Hadamard lifts
in that report preserve the common profile, keep both caps at
`O(N^(3/2))`, and have a `Theta(N^(3/2))` cap separation. Thus fixed `phi_6`
does not control cap even throughout the correct-scale class. This still
does not exhibit a gap among unknown near-minimizers. Any near-optimal-only
repair needs an independently checkable structure; defining the relevant
fiber by cap simply restates the problem.

This note asks whether one of three larger interfaces is genuinely more
tractable:

1. a growing Boolean/action profile on a purified bounded-operator class;
2. a growing induced-subgraph or moment state;
3. an augmented cut-code response state for amalgamation.

The conclusions below are **proved no-leverage results for standard uses of
these interfaces**, not a proof that every nonlinear compression is
impossible.

## 2. Bounded operator norm does not make a direct Boolean net small

Let `A` be a symmetric zero-diagonal signing of order `n`, put

```math
H_A(x)={1\over2}x^{\mathsf T}Ax,
```

and suppose

```math
\|A\|_{op}\le C\sqrt n.                                \tag{G1}
```

If `x,y` differ in `d` coordinates, symmetry gives

```math
\begin{aligned}
|H_A(x)-H_A(y)|
 &= {1\over2}|(x-y)^{\mathsf T}A(x+y)|\\
 &\le 2C\sqrt{n d(n-d)}.
\end{aligned}                                           \tag{G2}
```

Thus a Hamming net of radius `delta n` controls every energy from its sampled
values only to accuracy

```math
2C n^{3/2}\sqrt{\delta(1-\delta)}.                      \tag{G3}
```

Conversely, any radius-`delta n` covering set `S` of the Boolean cube obeys
the elementary sphere-covering bound

```math
|S|\sum_{j\le\delta n}{n\choose j}\ge 2^n.             \tag{G4}
```

For `0<delta<1/2`, this implies

```math
\log |S|\ge n(\log2-h(\delta))-O(\log n),              \tag{G5}
```

where `h` is binary entropy in natural units. Therefore a Lipschitz/net
certificate with error `epsilon_n n^(3/2)`, where `epsilon_n -> 0`, requires
`delta_n -> 0` and

```math
|S|\ge \exp((\log2-o(1))n).                             \tag{G6}
```

This has essentially the full exponential rate of the Boolean state space.
The result is deliberately scoped: it rules out tractability obtained merely
by operator purification followed by a Hamming net. It does not rule out an
algebraic description that optimizes over exponentially many spins
implicitly.

## 3. Generic tail control needs linear moment order

For uniform Rademacher `X`, a standard quadratic-chaos moment estimate under
(G1) is

```math
\|H_A(X)-\mathbb EH_A(X)\|_q
 \le K\bigl(n\sqrt q+C\sqrt n\,q\bigr)                 \tag{G7}
```

for `q>=2`, with an absolute `K` after harmless convention changes. For
`q<=n`, the right side is `O_C(n sqrt(q))`.

This estimate cannot certify a Boolean maximum at the `n^(3/2)` scale from
moments of order `q=o(n)`. Indeed, if any spin exceeds a threshold, then the
uniform-spin tail has mass at least `2^(-n)` (or `2^(-(n-1))` after
projectivizing). Markov at `t=c n^(3/2)` and (G7) can supply at best an
exponent of order

```math
q\log\sqrt{n/q}.
```

If `q=o(n)`, write `n/q=r -> infinity`; then

```math
q\log(n/q)=n{\log r\over r}=o(n),                      \tag{G8}
```

so this cannot beat the necessary `exp(-Theta(n))` threshold. Linear moment
order `q=Theta(n)` is necessary for this generic moment/entropy mechanism.

The switching profile through order `k` determines every energy moment up to
order `k`: in the expansion of `H_A^q`, every surviving monomial is an
even-degree multigraph with `q` edge occurrences and hence at most `q`
nonisolated vertices. Consequently, using local profiles to implement the
generic moment route requires

```math
k=Theta(n),                                             \tag{G9}
```

not a bounded or slowly growing local state. A special rigidity theorem for
competitive signings could bypass (G8); without one, “take a larger local
profile” has not reduced the entropy obligation.

## 4. Size of an induced switching-profile state

Root gauge identifies labeled switching classes on `k` vertices with the
`2^(binom(k-1,2))` signings on the remaining pairs. Quotienting by vertex
permutations gives at least

```math
{2^{\binom{k-1}{2}}\over k!}                            \tag{G10}
```

unlabeled classes and at most `2^(binom(k-1,2))`. Hence a complete induced
switching-profile vector through order `k` has

```math
2^{Theta(k^2)}                                          \tag{G11}
```

possible class coordinates at its top level. At the linear scale forced by
(G9), explicitly carrying this vector is far larger than the original edge
description. Carrying only the scalar moments avoids (G11), but `Theta(n)`
moments have `Theta(n^2 log n)` elementary bit size and, more importantly,
still face the sharp constants and truncated-moment ambiguity in the
extreme `2^(-n)` tail.

Thus neither the full local state nor its generic moment projection is a
bounded-complexity landing interface. A viable growing state would need
additional algebraic closure—such as an association scheme or code whose
few parameters implicitly determine the critical tail—and an all-order
realization theorem for that closed family.

## 5. What action compactness would have to prove

The purification theorem in
`concentration_compactness_boolean_profiles.md` is a verified reduction:
for each fixed accuracy, minimization may be restricted at arbitrarily small
leading cost to signings satisfying a bound of the form (G1). On that class,
action convergence is compact enough and the Boolean objective is
continuous.

The exact missing statement is nevertheless global in order:

> **All-order absorbing realization.** Every bounded-operator action limit
> arising from symmetric off-diagonal sign matrices has, at every sufficiently
> large prescribed order, a sign-matrix realization converging to the same
> action object and Boolean objective.

It would prove convergence: purify a liminf sequence, extract an action
limit, realize it along all orders, and let the purification loss tend to
zero. This is a real sufficient theorem and does not define the target-order
matrix using `M_m`. But no current construction makes it simpler than the
original scale-transfer problem.

In particular:

- proportional principal restriction rescales the normalized operator by
  `alpha^(-1/2)` and preserves the object only for `alpha -> 1`;
- a `k`-fold sign-block lift whose macro row sum is the required `sqrt(k)`
  forces `(1-o(1))k^2` Frobenius energy per block into microscopic modes;
- independent, Hadamard, and biased random residuals contribute a new
  leading-order Boolean component rather than a summable defect.

These are exact obstructions to the known sampling and lift proofs, not to
the realization theorem itself. Moreover, (G6) shows that the theorem cannot
be obtained just by adding a subexponential Boolean test net to an action
profile.

## 6. Cut-code amalgamation still needs an exponential response or a new invariant

The augmented cut-code identity and independent amalgamation inequality
proved in `second_phase_independent_abstraction.md` are

```math
M_{m+n}\le M_m+M_n+B_{m,n},                             \tag{G12}
```

where `B_(m,n)` is the minimum rectangular sign discrepancy. Bowlin's exact
theorem makes `B_(n,n)=Omega(n^(3/2))`; the independent bridge therefore
pays a leading defect.

For a bridge with `l` rows, the exact optimized response as a function of
all left spins has rank `2^(l-1)` in the column-type counts. Hence no smaller
**exact linear** state preserves the complete bridge response. A useful
approximate state must instead prove that only a compressed subset of
responses matters at `o(n^(3/2))` accuracy. Equations (G6)--(G9) explain why
generic nets and moments do not provide that compression.

One can record the child's scalar energy histogram, or even the number of
spins in each energy shell, but composition needs the joint alignment of
those spins with bridge responses. Supplying that full joint response is
equivalent to full bridge optimization. Therefore cut-code amalgamation is
not presently a distinct tractable state; it becomes one only after a new
proof-relevant invariant compresses this joint alignment.

## 7. Full-response approximation has quadratic information cost

There is a separate exact rigidity result for any state intended to retain
the **entire** Boolean energy landscape. Let `A,B` be order-`n` signings and
let `t` be the number of unordered edges on which they differ. Put `D=A-B`
and

```math
Q(D)=\max_{z\in\{+1,-1\}^n}|z^{\mathsf T}Dz|.
```

If `d_i` is the number of differing edges incident with vertex `i`, the
sharp elementary Khintchine inequality gives, for uniform Boolean `Y`,

```math
\mathbb E\|DY\|_1
 \ge\frac1{\sqrt2}\sum_i\left(\sum_jD_{ij}^2\right)^{1/2}
 =\sqrt2\sum_i\sqrt{d_i}
 \ge {2\sqrt2\,t\over\sqrt n}.                       \tag{G13}
```

Choose `y` attaining this expectation and
`x_i=sign((Dy)_i)`. For
`u=(x+y)/2` and `v=(x-y)/2`, symmetry gives

```math
x^{\mathsf T}Dy=u^{\mathsf T}Du-v^{\mathsf T}Dv.     \tag{G14}
```

Each of `u,v` is a partial sign vector with entries in `{0,+1,-1}`. Randomly
filling its zero coordinates shows that the absolute value of its quadratic
form is at most `Q(D)`. Equations (G13)--(G14) therefore prove

```math
\max_z|H_A(z)-H_B(z)|={Q(D)\over2}
 \ge {t\over\sqrt{2n}}.                               \tag{G15}
```

Suppose a finite state `s(A)` and decoder `F_s(z)` approximate every response
uniformly:

```math
\max_z|H_A(z)-F_{s(A)}(z)|\le\epsilon_n n^{3/2},
\qquad \epsilon_n\longrightarrow0.                   \tag{G16}
```

Two matrices in one state fiber then differ on at most
`2 sqrt(2) epsilon_n n^2=o(n^2)` edges by (G15). A Hamming-ball volume bound
in the `binom(n,2)` edge cube shows that every fiber has size
`2^(o(n^2))`. Since there are `2^(binom(n,2))` signings, the number of states
must be

```math
2^{\binom n2-o(n^2)},                                 \tag{G17}
```

requiring `binom(n,2)-o(n^2)` bits.

This proves that a uniformly accurate full-response state is essentially a
copy of the signing itself. Its scope is important: cap is only one scalar
functional of the landscape, so (G17) does not rule out a cap-specific
algebraic invariant. It does rule out presenting a compressed version of all
child or bridge responses as the missing bounded state.

## 8. Research judgment and exact surviving target

The following proposed states are now inactive as standalone mechanisms:

- fixed `phi_6`, by the exact scalable collision;
- `phi_k` with `k=o(n)` used only through generic moments, by (G8)--(G9);
- a subexponential Hamming test net on the purified class, by (G3)--(G6);
- the complete cut-code bridge response, because its exact rank is
  exponential and the standalone rectangular term is leading.

The only defensible version of a growing critical-scale state is an
**algebraically generated state with implicit exponential coverage**. A
concrete theorem would have to provide all three items below:

1. a family definable without `M_n`, with `O(n^(2-epsilon))` or otherwise
   demonstrably tractable parameters;
2. a uniform theorem that those parameters control the Boolean extreme tail
   to `o(n^(3/2))`, not merely finitely many moments;
3. a realization/composition operation at every large order whose accumulated
   `b_n=M_n^(2/3)` defect is geometrically summable.

Association schemes, two-graphs, and code families are plausible sources of
item 1, but the current conference and restriction examples do not supply
items 2 or 3. The all-order absorbing-realization theorem supplies exactly
the implication to convergence, but in its present unstructured form there
is no evidence that it is easier than the original cross-order landing
problem.

The next falsifiable mathematical target should therefore not be another
larger profile computation. It should be one explicit algebraically closed
family together with an order-filling operation, followed first by a proof
or counterexample to item 2 at the `o(n^(3/2))` scale. Without such an
operation, “growing state” is vocabulary for retaining more of the original
optimization rather than a reduction.
