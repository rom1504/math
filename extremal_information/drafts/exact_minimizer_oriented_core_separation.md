# Tail failure forces an orientation-separated principal core

**Status.** Task-local theorem draft awaiting independent audit.  This note
combines the principal-core inverse theorem with near-order principal
heredity.  It does not exclude the resulting core.

## 1. One-sided caps and a partition identity

For a hollow real symmetric matrix `A`, write

```math
H_A(x)={1\over2}x^TAx,
\qquad
P(A)=\max_xH_A(x),
\qquad
N(A)=-\min_xH_A(x),
\qquad
Q(A)=\max\{P(A),N(A)\}.                         \tag{OC.1}
```

### Lemma OC.1 (oriented principal superadditivity; archived)

For every partition `[n]=T\mathbin\dot\cup R`,

```math
\boxed{
P(A)\ge P(A[T])+P(A[R]),
\qquad
N(A)\ge N(A[T])+N(A[R]).}                       \tag{OC.2}
```

#### Proof

Choose Boolean maximizers `u` and `v` for the two positive principal caps.
The two full spins `(u,v)` and `(-u,v)` have the same internal energies and
opposite cross energies.  One of those cross energies is nonnegative, which
proves the first inequality.  Apply the first inequality to `-A` to obtain
the second. `square`

This is stronger than monotonicity of the absolute principal cap, but it is
only a lower bound: it does not control how the cross block changes the full
optimizer.  The identity is archived in `nested_restriction_paving.md`,
equation (1), `concentration_compactness_boolean_profiles.md`, equation (8), and ledger
Section 10.13; no novelty is claimed for OC.1 itself.

## 2. Exact minimizers separate the orientations

Let `M_n` denote the minimum absolute cap over order-`n` signings.  Globally
negate an exact minimizer if necessary so that

```math
P(A_n)=Q(A_n)=M_n.                               \tag{OC.3}
```

### Theorem OC.2 (a global-scale small core forces the opposite shore)

Let `A_n` satisfy (OC.3), let `T_n subseteq[n]` have `k_n=|T_n|=o(n)`, and
put `R_n=[n]\setminus T_n`, `m_n=n-k_n`.  Suppose that for a fixed `t>0`,

```math
P(A_n[T_n])\ge(t-o(1))n^{3/2}.                  \tag{OC.4}
```

Then

```math
\begin{aligned}
Q(A_n[R_n])&=M_n-o(n^{3/2}),\\
N(A_n[R_n])&=M_n-o(n^{3/2}),\\
P(A_n[R_n])&\le M_n-P(A_n[T_n]),\\
N(A_n[T_n])&=o(n^{3/2}),\\
N(A_n)&=M_n-o(n^{3/2}).                          \tag{OC.5}
\end{aligned}
```

In particular, the small core is positive-dominant at the global energy
scale, while its near-minimal complement is negative-dominant at that same
scale.  The full exact minimizer is asymptotically saturated in both
orientations.

#### Proof

The random-bridge near-order inequality in Theorem 36.15 gives

```math
M_n-M_{m_n}=o(n^{3/2}).                          \tag{OC.6}
```

Because `A_n[R_n]` is an order-`m_n` signing and a principal restriction of
`A_n`,

```math
M_{m_n}\le Q(A_n[R_n])\le M_n.                  \tag{OC.7}
```

This proves the first line of (OC.5).  By positive superadditivity,

```math
P(A_n[R_n])\le M_n-P(A_n[T_n])
             \le M_n-(t-o(1))n^{3/2}.           \tag{OC.8}
```

Equations (OC.7)--(OC.8) force the negative orientation to realize the
absolute cap of the complement for all large `n`, and hence

```math
N(A_n[R_n])=Q(A_n[R_n])=M_n-o(n^{3/2}).          \tag{OC.9}
```

Negative superadditivity and `N(A_n)<=Q(A_n)=M_n` now give

```math
0\le N(A_n[T_n])
\le M_n-N(A_n[R_n])=o(n^{3/2}).                 \tag{OC.10}
```

The same two inequalities imply the last line of (OC.5). `square`

The exact-minimizer assumption is used twice: it fixes the cap in (OC.3)
and makes every sublinear complement near-minimal through (OC.6)--(OC.7).
For a generic bounded-cap signing there is no reason for the two principal
orientations to separate.

## 3. The intrinsic signed-cut form

Choose a positive maximizer `u` of the core and gauge it to the all-positive
word, so `b_(ij)=a_(ij)u_i u_j` on `T_n`.  For `S subseteq T_n`, put

```math
w(S)=\sum_{i\in S,\ j\in T_n\setminus S}b_{ij}.
```

The exact flip identity is

```math
H_b(\mathbf1^{S})=P(A_n[T_n])-2w(S).             \tag{OC.11}
```

### Corollary OC.3 (cut-positive, almost one-sided core)

Under OC.2,

```math
0\le w(S)\le {P(A_n[T_n])+N(A_n[T_n])\over2}
             ={P(A_n[T_n])\over2}+o(n^{3/2})    \tag{OC.12}
```

for every `S`.  Moreover

```math
\max_Sw(S)-\mathbb E_Sw(S)
 ={N(A_n[T_n])\over2}=o(n^{3/2}),               \tag{OC.13}
```

where `S` is a uniform random subset.  Thus every signed cut is
nonnegative, and a uniform cut is within `o(n^(3/2))` in expectation of
the maximum signed cut.  Equivalently, all signed row sums are nonnegative
and sum to `2P(A_n[T_n])`.

#### Proof

The lower bound in (OC.12) is maximality of the gauged all-positive word.
The upper bound follows by minimizing (OC.11).  Every edge crosses a
uniform random cut with probability `1/2`, so
`Ew(S)=P(A_n[T_n])/2`; the maximum is
`(P+N)/2`.  This proves (OC.13).  Singleton cuts are the signed row sums,
and their sum counts every edge twice. `square`

If `k_n=(\sqrt{2t}+o(1))n^{3/4}`, then the trivial edge-count lower bound
is asymptotically tight.  In that special case, after the above switching,
all but `o(k_n^2)` core edges are positive.  More generally OC.2 does not
force clique structure when `k_n\gg n^{3/4}`: the average signed row bias is
only `2P/k_n`, a vanishing fraction of `k_n` in that regime.

## 4. Consequence for the rare-tail route

Apply Corollary PC.3 from
`rare_upper_tail_principal_core_dichotomy.md` to a sequence of globally
oriented exact minimizers.  If the upper tail at any fixed threshold `t`
has vanishing exponential rate, then its core satisfies OC.2--OC.3.
Therefore the remaining structural target can be stated more sharply than
plain no-core:

> **Exact-minimizer core-orientation lemma (`L_balance`, sequential form).**
> There are fixed `0<t<c_-` and `c>0` such that, for every sequence of
> globally oriented exact minimizers and every `T_n=o(n)` with
> `P(A_n[T_n])\ge(t-o(1))n^(3/2)`, one has
> `\liminf N(A_n[T_n])/n^(3/2)\ge c`.

Any sequential formulation excluding `P(A[T])=Theta(n^(3/2))` together
with `N(A[T])=o(n^(3/2))` suffices.  PC.3 and OC.2 show

```text
L_balance
  => positive upper-tail entropy deficit
  => matched-roof code (Theorem 21.8)
  => scalar physical contextual packing (BR.2--BR.3).
```

This is a genuine narrowing of the *shape of a counterexample*: failure is
not an arbitrary spectral spike or arbitrary principal core, but a
cut-positive, strongly one-sided core glued to an oppositely oriented
near-minimal complement.  On exact minimizers OC.2 makes this balance
statement essentially equivalent, up to fixed slack, to excluding the
principal-core obstruction; it is not a new weaker/easier missing arrow.
It gives no recurrence or convergence theorem.

## 5. Archive comparison and classification

- Principal-cap monotonicity, one-sided superadditivity, and near-order
  random padding are archived.  The new synthesis is the **two-orientation**
  consequence after a fixed-level core is supplied.
- The selected-prior/common-active-face route concerned optimizer transfer;
  OC.2 uses no selected prior and identifies no parent optimizer.
- The two-cap and orientation-ceiling artifacts anticipated the importance
  of sign, but do not state the forced positive-core/negative-complement
  decomposition (OC.5).
- A switched positive clique realizes the one-sided core geometry in
  isolation, so no theorem about arbitrary sign matrices can exclude it.
  Exact minimality and the coupling to the opposite near-minimal shore are
  indispensable.

Classification:

```text
PROVES A STRUCTURAL REFINEMENT:
  tail failure -> arbitrary spectral spike -> arbitrary zero-density core
  is refined to
  tail failure -> orientation-separated cut-positive core/complement pair.

DOES NOT PROVE:
  L_balance, a scalar packing unconditionally, a cross-order recurrence,
  or convergence.
```
