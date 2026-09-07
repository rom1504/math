# An exact delocalized relaxation limit, and its additional width obligation

2026-09-07. The relaxed limit below is proved by direct sums. It is NOT the
original limit, and the strongest suggested sign-recovery hypothesis would
also prove the currently unproved original width-to-cap equivalence. This
scope issue prevents treating the relaxation as an automatically easier
replacement for the selected-original-child problem.

## 1. Fixed-cutoff thermodynamic limit

For `delta>0`, let `R_n(delta)` consist of real hollow symmetric matrices
`B` with

```math
\sum_j B_{ij}^2=1\quad\hbox{for every }i,
\qquad \max_{ij}|B_{ij}|\le\delta.
```

The class is compact and nonempty for every
`n>=n_0(delta):=ceil(delta^(-2))+1`: any complete signing divided by
`sqrt(n-1)` is feasible. Define

```math
d_n(\delta)=\min_{B\in R_n(\delta)}Q(B).
```

Direct sums preserve the row norm and the entry cutoff. The triangle
inequality therefore gives `d_(n+m)<=d_n+d_m` whenever both terms are
defined. The eventual-domain version of Fekete's argument is elementary:
for any fixed admissible `k`, write `n=qk+r` with
`n_0<=r<n_0+k`. Compose `q` minimizing `k`-blocks and one feasible remainder.
The finitely many remainder costs are bounded independently of `n`.
Consequently

```math
c(\delta):=\lim_{n\to\infty}\frac{d_n(\delta)}n
          =\inf_{n\ge n_0(\delta)}\frac{d_n(\delta)}n.     (1)
```

Since the classes decrease when `delta` decreases, `c(delta)` increases.
It is bounded above by the original normalized liminf: for a full signing
`A`, `B=A/sqrt(n-1)` gives

```math
\frac{Q(B)}n=\frac{Q(A)}{n^{3/2}}\sqrt{\frac n{n-1}}.
```

Thus the canonical number

```math
c_{\rm deloc}:=\lim_{\delta\downarrow0}c(\delta)
\le\liminf_n\frac{M_n}{n^{3/2}}                         (2)
```

exists. This proof uses no claimed sign flatification or matrix LDP.

## 2. The relaxation identifies cap and width exactly at the limiting level

Write `P(B)=max H_B`, `R(B)=-min H_B`, and `W(B)=(P+R)/2`. For every
feasible `B`, the direct sum `B⊕(-B)` is feasible at order `2n` and has

```math
P(B\oplus(-B))=R(B\oplus(-B))=P(B)+R(B)=2W(B).
```

Since `W<=Q`, (1) implies both inequalities in the identity

```math
c(\delta)=\inf_{n\ge n_0(\delta)}\;
              \inf_{B\in R_n(\delta)}\frac{W(B)}n.        (3)
```

One inequality uses `W(B)<=Q(B)`; the other applies (1) at order `2n` to
`B⊕(-B)`. In particular near-optimal values of the relaxed thermodynamic
problem have exactly balanced block representatives at even orders. The
missing cross-block entries are zeros, which are legal HERE. They are
precisely the leading-order filling issue in the original signing problem.

Let `W_n^sign` be the original minimum of half-width over complete signings.
Equation (3) also gives

```math
c_{\rm deloc}\le\liminf_n W_n^{\rm sign}/n^{3/2}.
```

This extra comparison is why full recovery of the relaxation is a stronger
assertion than existence of the original cap limit alone.

## 3. All-order delocalized representatives do exist in the relaxation

Choose `delta_j ->0`. For each `j`, a finite block in `R_(k_j)(delta_j)`
approximates the infimum in (1), or use its balanced double from (3).
Repeat it on all sufficiently large orders and use an admissible remainder
of size in `[n_0(delta_j),n_0(delta_j)+2k_j)`. Choose the order thresholds
so large that this remainder costs at most `1/j` per vertex. Let `j=j(n)`
increase sufficiently slowly. This constructs actual real matrices `B_n`
at EVERY sufficiently large order with

```math
\max_{ij}|(B_n)_{ij}|\to0,
\qquad Q(B_n)/n\to c_{\rm deloc}.
```

The operator norm at each fixed accuracy stage is bounded by the norm of
the selected finite block and of its finitely many remainders, independently
of the large replicated order. No one operator cutoff uniform in `j` is
claimed. Nor does the real matrix `sqrt(n-1) B_n` have uniformly bounded
entries as `n` grows: that would be a different, much narrower class.

The order of limits matters. If the entry cutoff is not sent to zero, a
four-vertex CHSH block already gives a spurious low asymptotic constant:
`B=[[0,H_2],[H_2^T,0]]/sqrt2` has row norm one and
`Q(B)/4=1/(2sqrt2)`. Repeating this fixed block retains its nonvanishing
normalized entry size `1/sqrt2`. It does not give a delocalized counterexample.

## 4. Exact conditional implication and why it is stronger than needed

Suppose one proved a uniform recovery theorem: for every fixed small
`delta`, every `B in R_n(delta)` of cap `O(n)` admits an actual full signing
`A_n` with

```math
\frac{Q(A_n)}{n^{3/2}}\le\frac{Q(B)}n+
            \omega(\delta)+o_n(1),\qquad
\omega(\delta)\to0,                                    (4)
```

where the error is uniform over the stated input class. Applying (4) to
the fixed-cutoff minima, first taking `n -> infinity` and then
`delta ->0`, would give `limsup M_n/n^(3/2)<=c_deloc`. Together with (2),
this proves the original limit without needing a summable per-level defect.

But (3) then also forces

```math
\lim_n\frac{W_n^{\rm sign}}{n^{3/2}}
=\lim_n\frac{M_n}{n^{3/2}}=c_{\rm deloc}.
```

So (4) contains BOTH global flat-sign recovery AND original width-to-cap
equivalence. Neither follows from the new same-order balancing operation,
which starts from the larger endpoint `Q`, not from `W`. Neither follows
from a lower bound for weighted profiles. No recovery theorem is proved
here, and the continuous minimizing-profile stationarity equations alone
have not supplied one.

A narrower original-seed replica recovery could avoid this additional
width obligation, but would still need an actual construction with joint
control of all new and rewritten old edges. Merely noting the exact relaxed
Fekete limit does not provide that construction.
