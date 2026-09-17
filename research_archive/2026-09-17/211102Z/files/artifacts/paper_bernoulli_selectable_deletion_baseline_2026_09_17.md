# Selectable exact-ground cheapness from the one-vertex identity

2026-09-17. Elementary variational baseline, not claimed as original
progress. It is a consequence of the exact insertion identity already
proved in [the cap-discrepancy archive](cap_discrepancy_insertion.md).
The subsequence and near-minimizer qualifiers are indispensable.

Write M_n for the minimum full-sign cap and
c_inf=liminf M_n/n^(3/2). Then there are orders n_j tending to infinity,
actual full signings B_j of order n_j, and physical columns h_j such that

```
Q(B_j)/n_j^(3/2)->c_inf,
0<=Q(B_j)-M_(n_j)<= (3c_inf/2+o(1))sqrt(n_j),
|h_j dot x| <= (3c_inf/2+o(1))sqrt(n_j)
                         +Q(B_j)-|H_(B_j)(x)|
                     for EVERY Boolean x.                    (1)
```

In particular the ENTIRE exact absolute ground code of B_j has
deterministic response at most (3c_inf/2+o(1))sqrt(n_j). Classical
Hadamard orders give c_inf<=1/2, so the coefficient is at most 3/4,
strictly smaller than sqrt(2/pi). This does not need the stronger
reported numerical bound or existence of a normalized cap limit.

## Proof with explicit limit order

Principal restriction and averaging give M_n<=M_(n+1). Fix delta>0
and then choose arbitrarily large N such that

```
M_N/N^(3/2)<=c_inf+delta^2,
M_k/k^(3/2)>=c_inf-delta^2  for EVERY k>=floor((1-delta)N).
```

This is possible by the definition of liminf, taking N beyond the
eventual lower-envelope threshold. Put m=floor((1-delta)N). Summing
the nonnegative increments on [m,N) shows that some n in that interval
obeys

```
M_(n+1)-M_n
 <=[(c_inf+delta^2)N^(3/2)
           -(c_inf-delta^2)m^(3/2)]/(N-m)
 <=[3c_inf/2+O(delta)+o_N(1)]sqrt(n).               (2)
```

Here the last inequality also uses n>=m and the Taylor expansion of
(1-delta)^(3/2). Moreover monotonicity gives

```
c_inf-delta^2 <=M_n/n^(3/2)
 <=(c_inf+delta^2)(N/m)^(3/2)=c_inf+O(delta)+o_N(1).
```

Take delta_j down to zero and then choose N_j sufficiently large at
each fixed delta_j; this diagonal order of choices produces (2) with
o(1) and n_j tending to infinity at the liminf value.

Take an EXACT minimizing signing of order n+1, delete any one vertex,
and call the remaining signing B and its incident column h. The exact
absolute identity is

```
M_(n+1)=max_x[|H_B(x)|+|h dot x|].                  (3)
```

Since M_n<=Q(B)<=M_(n+1), equations (2)--(3) give all of (1).
No spectral assumption, Gaussian replacement, or selection among
query laws occurs.

## Collision and precise scope

The mathematical mechanism is the archived insertion identity and
telescoping near a liminf order. This is recorded as a baseline, not
as a new route to convergence. The archive already explains why an
unsummable derivative-scale local error cannot imply convergence.

- B_j is a principal core of an exact optimizer at order n_j+1; it
  need not be an exact optimizer at order n_j.
- The conclusion holds on a selected subsequence, not at all orders.
- The cheap law here is a SINGLE physical column. Global sign
  symmetrization centers it but leaves covariance h_j h_j^T, not I.
- The all-query deficit allowance in (1) has coefficient ONE before
  repetition. Repeating the same column q times multiplies its
  deficit allowance by q and does not prove a useful macroscopic
  extension bound.
- No isotropic reusable law, entropy profile at fixed normalized
  deficit, or cross-order recurrence is supplied.

Thus selectable exact-ground cheapness is already a variational fact;
the harder reusable-law and escape questions remain substantive.
