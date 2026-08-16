# All-order action recovery: architecture audit

Date: 2026-08-16.

Status: first checkpoint draft.  Claims are labelled **proved**, **published**,
**obstructed**, or **open**.  The long project archive was consulted only
after the independent theorem packets and proposals were frozen.

## 1. Target and scale

For an order-\(n\) hollow signing \(A\), let

```math
T_A=A/\sqrt n,
\qquad
\Phi(T_A)=\frac{2Q(A)}{n^{3/2}}.
```

The exact minimal implication and all error budgets are proved in
`minimal_all_order_action_recovery.md`.  Its structural hypothesis
\(\mathrm{AR}_{\min}^{\to}\) asks only for one selected purified liminf
cluster at each member of a null tolerance sequence, directed one-profile
recovery, and an upward ratio-dense set of target orders.  If recovered
operator bounds are \(D_m\) and directed errors are \(\delta_m\), the
required scale is

```math
D_m\sqrt{\delta_m}\longrightarrow0.
```

The even weaker objective-only statement is equivalent to convergence and
is not counted as a reduction.

## 2. Architecture A: sampling plus correction

### Published input

Bounded graphon sampling realizes a fixed bounded kernel at every prescribed
large order in cut distance.  Unbounded graphon sampling likewise assumes a
fixed uniformly controlled \(L^p\) kernel.  Action convergence gives compactness
of bounded operator families, but Backhausz--Szegedy Section 11 still obtains
convergence of normalized iid sign matrices only after selecting a subsequence;
convergence along all natural orders is explicitly left open.

### Scale audit

The integral kernel corresponding to \(A/\sqrt n\) on the uniform probability
space is \(\sqrt n A\), whose \(L^p\) norm grows as \(\sqrt n\).  Ordinary
dense sampling therefore controls the wrong normalization.  Centering produces
the zero graphon and deletes the entire fluctuation object that carries
\(Q(A)/n^{3/2}\).

There is also a proved projective obstruction.  A jointly exchangeable
infinite sign array whose order-\(n\) restrictions have tight
\(\|A_n\|_{op}/\sqrt n\) must have iid Rademacher edges, and almost surely

```math
\liminf_n \frac{Q(A_n)}{n^{3/2}}
\ge \frac23\sqrt{\frac2\pi}=0.531923\ldots>\frac12.
```

See `exchangeable_recovery_obstruction.md`.

### Verdict

**Obstructed** for ordinary graphon sampling and for a single projectively
consistent exchangeable array.  **Open** for nonprojective, order-dependent,
globally conditioned laws.  Such a law must be specified by information other
than the unknown target-order optimum.

## 3. Architecture B: blow-up plus pseudorandom residual

### Proved obstruction

To preserve a base macro coefficient in a \(k\)-fold sign blow-up, every
off-diagonal \(k\)-by-\(k\) block needs row sum
\((1+o(1))\sqrt k\).  Removing that mean leaves exact Frobenius mass

```math
\|R\|_F^2=k^2-(1+o(1))k=(1-o(1))k^2.
```

Thus the microscopic residual cannot be norm-small.  Independent residuals
have a leading Boolean supremum; Hadamard residuals create a leading tensor
channel; constant blocks amplify the macro action by \(\sqrt k\).  These are
the exact obstructions in `bounded_op_signed_realization.md` and
`random_biased_lift_no_go_phase2.md`.

### Minimal theorem that would rescue the route

For a finite approximant \(A_k\) of the selected cluster and every sufficiently
large fibre size \(r\), construct sign microblocks with the forced Frobenius
mass but prove directly that the complete lifted operator has no outer
one-profile beyond the target, with error \(o(D^{-2})\).  The cancellation must
be taken before any norm or scalar-channel bound.

This statement is not yet a reduction: without an independently defined
algebraic microblock class it simply restates directed recovery.

### Verdict

**Obstructed** for independent, Hadamard, constant, or separately paid
mean-plus-residual lifts.  **Open but presently class C** for an unspecified
absorbing residual.

## 4. Architecture C: weighted realization plus global sign rounding

A weighted matrix may realize a regular operator at every order.  The final
rounding needs either

```math
\|A_n-W_n\|_{op}=o(\sqrt n)
```

or an equally strong joint same-spin/profile estimate.  Merely bounding the
final \(\|A_n\|_{op}\) by \(O(\sqrt n)\) says nothing about the rounding defect.

If a positive fraction of entries of \(W_n\) remain a fixed distance from
\(\{\pm1\}\), every sign rounding has

```math
\|A_n-W_n\|_F=\Omega(n),
\qquad
\|A_n-W_n\|_{op}=\Omega(\sqrt n).
```

Consequently a uniform little-\(o\) operator-rounding theorem is false on
genuinely fractional inputs.  Modern matrix-discrepancy theorems reach a
constant times the natural square-root scale; that still permits a fixed
leading change in \(Q/n^{3/2}\).

There is, however, a positive sign-near regime.  If the total fractional
variance satisfies

```math
V(W_n):=\sum_{i<j}(1-w_{ij}^2)=o(n^2),
```

independent biased rounding, scalar Bernstein, and a union bound over all
spins directly give an exact signing with Boolean error

```math
O\bigl(\sqrt{nV(W_n)}+n\bigr)=o(n^{3/2}).
```

If preservation of the directed profile is wanted, delete the \(o(n)\) rows
of excessive variance.  On the retained matrix, the Bandeira--van Handel
nonhomogeneous random-matrix bound gives some exact symmetric sign rounding
with

```math
\|A_n-W_n\|_{op}=o(\sqrt n).
```

The proof and the resulting weighted recovery statement are in
`sign_near_weighted_recovery.md`.

### Verdict

**Solved** for globally sign-near weighted inputs with total fractional
variance \(o(n^2)\).  **Obstructed** for arbitrary fractional matrices.  This removes the
last exact-sign rounding obligation if an all-order weighted realization can
be made sign-near; naive blow-ups have \(v=\Theta(n)\) and remain outside the
theorem.

## 5. Architecture D: approximate design and absorption

Under a fixed \(2\to2\) bound \(C\), rounding a test function to mesh \(h\)
changes its one-profile law by at most

```math
h+\sqrt{Ch}.
```

Thus the abstract profile has a finite net at every fixed accuracy, with
complexity depending on \((C,\epsilon)\) but not on \(n\).  This is a genuine
finite-state observation.  It does not turn AR into finitely many local design
constraints: the outer profile condition must hold for every coloring of
\(n\) vertices by the fixed alphabet, still exponentially many colorings.

Keevash-type decomposition, Kuperberg--Lovett--Peled balancing, nibble, and
absorption theorems impose fixed local templates/statistics subject to
divisibility and extendability.  No cited theorem controls this universal
local-global coloring set at fluctuation scale.

### Verdict

**Open**, but no strict reduction yet.  The exact missing design theorem is a
separation-oracle or absorption principle that enforces the entire outer
finite-alphabet profile without enumerating it.  Calling the fixed-accuracy
state finite does not supply that theorem.

## 6. Architecture E: near-order transfer

### Proved theorem

Principal deletion is lossless:

```math
Q(A[S])\le Q(A).
```

Randomly signing the \(r=nh+\binom h2\) new edges gives a deterministic
order-\(N=n+h\) completion with

```math
Q(B)\le Q(A)+\sqrt{2r(N+1)\log2}.
```

Therefore \(h=o(n)\) costs \(o(n^{3/2})\).  In particular, if a family of
good orders is upward ratio-dense, its normalized upper bound transfers to
all orders and convergence follows when that bound is the liminf.

### Boundary

This solves the order-mismatch step once good orders have multiplicative gaps
\(1+o(1)\).  It does not make an arbitrary liminf subsequence ratio-dense.
Arbitrary padding can cost \(nh\), so balancing is essential.

### Verdict

**Proved strict weakening of every-order recovery to upward ratio-dense
recovery.**  The remaining density theorem is open.

## 7. Gamma-convergence translation

On the disjoint union of exact signing spaces, set

```math
F_n(A)=\frac{Q(A)}{n^{3/2}}=\frac12\Phi(T_A)
```

and use bounded-operator action convergence as topology.  The proved
continuity theorem supplies both lower and upper continuity along any already
convergent uniformly bounded signing sequence.  The missing Gamma-limsup is
not recovery of a vertex spin for a supplied kernel; it is recovery of the
kernel itself by exact signs at sufficiently dense orders.

Published graphon Gamma-convergence theorems assume the finite graph kernels
already converge and then recover vertex labels or phase fields.  Their
quantifiers point in the opposite direction from directed outer-profile
recovery.  Thus Gamma terminology identifies the missing axiom but does not
currently prove it.

### Verdict

**Equivalent framework**, not a reduction, absent a new exact-sign density
theorem.

## 8. Near-minimizer rigidity

The repository proves tolerance-dependent spectral purification: for every
fixed objective tolerance there are competitive signings with a fixed
normalized operator bound.  This is enough for the minimal implication; a
single bound uniform as tolerance tends to zero is unnecessary.

No verified additional property \(P\) currently has all three required
features:

1. all near-minimizers satisfy \(P\);
2. \(P\) is action-closed; and
3. every object with \(P\) has upward-ratio-dense exact-sign realizations.

Conference identities are sufficient for selected \(1/2\)-scale families but
are not known to be forced on arbitrary near-minimizers.

## 9. Current director judgment

The campaign has produced three theorem-level advances about the AR
architecture:

1. full every-order/full-profile recovery is reduced exactly to selected,
   directed one-profile recovery on an upward ratio-dense order set, without
   fixed-\(C\) regularization; and
2. projectively consistent exchangeable recovery is rigorously impossible at
   the extremal scale; and
3. sign-near weighted recovery with total fractional variance \(o(n^2)\)
   rounds to exact sign recovery with no leading loss.

Neither result proves \(\mathrm{AR}_{\min}^{\to}\), nor has a counterexample
to that selected extremal statement been found.  The remaining possible
architectures are sign-near weighted realization, nonprojective
microcanonical realization, and a genuinely joint profile-level absorption
theorem.  They proceed to independent specialist formulation before archive
classification; ordinary sampling, independent residuals, and unrestricted
norm-only rounding do not.
