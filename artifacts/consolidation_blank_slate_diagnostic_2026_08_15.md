# Consolidation blank-slate diagnostic

Date: 2026-08-15.

Status: **three independently generated formulations evaluated; none supplies
a strict or defensible new reduction**.

This audit followed the README stopping rule after two consecutive
substantive checkpoints without primary progress.  Each candidate was stated
from the original problem before its author consulted the project vocabulary
or ledger.  Only afterward was it compared with the accumulated evidence.

Throughout,

```math
M_n=min_{a_{ij}\in\{\pm1\}}
\max_{x\in\{\pm1\}^n}
\left|\sum_{i<j}a_{ij}x_ix_j\right|.
```

## 1. Candidate A: uniform block completion

For signings `A` and `B` of orders `p` and `q`, and a rectangular signing
`R`, let `A *_R B` retain `A,B` on its diagonal blocks and use `R` across
them.  Write

```math
D(A)=\max_x|H_A(x)|.
```

The independently proposed theorem was: there exist universal
`C<infinity` and `delta>0` such that, for every `A,B`, some `R` satisfies

```math
D(A*_RB)
\le\sqrt{p+q\over p}\,D(A)
+\sqrt{p+q\over q}\,D(B)
+C(p+q)^{3/2-\delta}.                             \tag{1.1}
```

### Exact convergence mapping

Put `u_n=M_n/sqrt(n)`.  Applying (1.1) to minimizers gives

```math
u_{p+q}\le u_p+u_q+C(p+q)^{1-\delta}.              \tag{1.2}
```

Balanced binary merging accumulates per-vertex error

```math
\sum_{j\ge0}O((2^jk)^{-\delta})=O(k^{-\delta}),
```

so the standard almost-subadditive argument makes
`u_n/n=M_n/n^(3/2)` converge.

The candidate proof mechanism was entropy-stratified partial coloring and
chaining of the rank-one bridge constraints.  A concrete falsifier is a
comparable-block family with a fixed leading gap:

```math
\min_RD(A*_RB)
-\sqrt{p+q\over p}D(A)
-\sqrt{p+q\over q}D(B)
\ge\eta(p+q)^{3/2}.                               \tag{1.3}
```

### Ledger comparison

This is not a new mechanism.  The exact identity

```math
D(A*_RB)=\max_{x,y}
\left(|H_A(x)+H_B(y)|+|x^TRy|\right)               \tag{1.4}
```

follows by independently reversing one block of spins.  It is the existing
fixed-child, state-dependent bridge objective.  At equal blocks, (1.1) asks
for `2 sqrt(2) M_n` at leading order, exactly the current
`M_n^(2/3)` composition scale.  The proposed chaining metric and its random
rectangular norm floor have already been analyzed.  Exact small bridge tests
also fail for some optimal representatives; this does not asymptotically
disprove (1.1), but weighs against its stronger `for every A,B` quantifier.

**Judgment:** reject.  The candidate hides the full bridge selection problem
and strengthens its quantifier.  It removes no obligation.

## 2. Candidate B: bounded dual-cycle determinacy

Let

```math
\mathcal C_n^+
=\{(\sigma x_ix_j)_{i<j}:\sigma,x_i\in\{\pm1\}\}.
```

With `E_n=binom(n,2)`, the exact coding map is

```math
M_n=E_n-2\rho(\mathcal C_n^+),
```

and the dual is the even-cardinality Eulerian sector

```math
D_n=(\mathcal C_n^+)^\perp
=\{F:\partial F=\varnothing,\ |F|=0\pmod2\}.
```

For a signing `a`, define

```math
P_n(a,\beta)
={1\over n}\log\left(2^{-n}\sum_x
\cosh{\beta H_a(x)\over\sqrt n}\right),
```

```math
W_a(z)=\sum_{F\in D_n}a_Fz^{|F|}.
```

The high-temperature expansion is exact:

```math
P_n(a,\beta)
={E_n\over n}\log\cosh{\beta\over\sqrt n}
+{1\over n}\log W_a\left(\tanh{\beta\over\sqrt n}\right).   \tag{2.1}
```

Let `T_{n,L}` replace the last logarithm by its Taylor polynomial through
degree `L`.  The proposed theorem was

```math
\lim_{n\to\infty}\min_aT_{n,L}(a,\beta)
\quad\hbox{exists for every fixed }L,\beta,          \tag{2.2}
```

and

```math
\lim_{L\to\infty}\limsup_{n\to\infty}
|\min_aP_n(a,\beta)-\min_aT_{n,L}(a,\beta)|=0.       \tag{2.3}
```

### Exact convergence mapping

The uniform sandwich

```math
\beta{M_n\over n^{3/2}}-\log2-o(1)
\le\min_aP_n(a,\beta)
\le\beta{M_n\over n^{3/2}}                          \tag{2.4}
```

shows that convergence of the minimum pressure for every fixed `beta`,
followed by `beta -> infinity`, forces convergence of the ground-state
sequence.

A concrete falsifier is a fixed `beta,epsilon>0` such that every fixed `L`
misses the true minimum pressure by at least `epsilon` on infinitely many
orders.

### Ledger comparison

The code identity and (2.1) are already proved in the project.  More
importantly, fixed cycle/replica depth has the exact advertised falsifier:

- deleting one vertex opens exponentially many even boundary sectors;
- fixed conference diagrams are universal while a planted Boolean
  eigenvector changes pressure by order `n`;
- fixed replicas do not determine the least output likelihood; and
- correct-scale cut-code certificates require degree `Theta(n)`.

These facts do not logically exclude an unexplained lucky selection principle
for the minimum in (2.3), but they invalidate the claimed bounded-complexity
mechanism.  Letting `L` grow far enough to see the planted resonance restores
the complete signed Eulerian/coset histogram.

**Judgment:** reject.  Fixed `L` is support-blind; the surviving growing-degree
version has not been shown weaker than the original tail.

## 3. Candidate C: symplectic/Witt nonconvergence

The provisional subsequences were

```math
n_k^-=4^k,qquad n_k^+=2\cdot4^k.                  \tag{3.1}
```

Their ratio is two, so they pass the necessary multiplicative-separation
test.  On `V=F_2^k x F_2^k`, let

```math
B((p,q),(p',q'))=p\cdot q'+q\cdot p'
```

and

```math
H_B(u,v)=(-1)^{B(u,v)},\qquad A=H_B-I.
```

Then `H_B^2=nI`.  For suitable quadratic bent phases `x`,

```math
H_Bx=-\sqrt n\,x,
```

so the cap is exactly

```math
\operatorname{cap}(A)
={n\over2}(\sqrt n+1)
=\left({1\over2}+{1\over2\sqrt n}\right)n^{3/2}.   \tag{3.2}
```

Thus the candidate good family reaches only the already known constant
`1/2`.  To separate the proposed second subsequence it would need a fixed
`epsilon>0` and the all-signings theorem

```math
M_{2\cdot4^k}
\ge\left({1\over2}+\epsilon\right)(2\cdot4^k)^{3/2}.\tag{3.3}
```

### Falsification

For every large `n`, choose a symmetric Paley conference order `N>=n` with
`N/n -> 1` and take an `n`-vertex principal submatrix.  Its operator norm is
at most `sqrt(N-1)`, hence

```math
M_n\le {n\over2}\sqrt{N-1}
=\left({1\over2}+o(1)\right)n^{3/2}.              \tag{3.4}
```

Equation (3.4) applies to every subsequence and directly contradicts (3.3).
The augmented cut code also has only bounded divisibility, so Witt parity of
a selected Cayley family imposes no restriction on arbitrary edge signings.
Exact small values show that the symplectic family is not even landing-optimal
at orders four and sixteen.

**Judgment:** close this mechanism.  Genuine nonconvergence could still occur
between two constants strictly below `1/2`, but would require both a strict
low-subsequence construction and an all-signings lower theorem on
multiplicatively separated epochs.  No arithmetic or design mechanism in the
archive provides either half.

## 4. Global research judgment

The three candidates were mathematically distinct before ledger comparison:

1. deterministic cross-order composition;
2. finite-temperature dual-cycle compactification; and
3. arithmetic nonconvergence.

After comparison, the first is an equivalent bridge obligation, the second is
a bounded-degree route already defeated by support-sensitive resonances, and
the third is rigorously falsified in its required constant range.  None
provides a strict reduction, a new scalable construction, or a theorem target
with evidence of greater tractability.

The recent moving-representation feedback has already been applied as far as
the current evidence justifies: it yielded an exact operator same-switch
inequality and strong partial-matching no-go theorems.  The suggested
executor/verifier/research-director architecture was also used in that audit.
Changing the orchestration again cannot replace a missing mathematical
mechanism.

**Decision:** pause autonomous route generation and seek external mathematical
review.  A justified restart needs at least one genuinely new input of one of
the following forms:

- a graph-orbit moving representation with coefficient one and an algebraic
  rooted mass theorem;
- a diffuse bridge law with `O(n)` entropy and a proved linear pressure gain;
- a strict reduction of the existential bridge problem to a state of
  demonstrably smaller information content; or
- a scalable construction/theorem producing two normalized constants below
  `1/2` on multiplicatively separated epochs.

Without such an input, another autonomous wave would predictably generate an
equivalent sufficient condition or a route-specific falsifier and would
violate the project's stopping discipline.
