# Correlated skew orientations: exact normalized law and bounded test

2026-09-07. **Exact identities and finite diagnostics only. No asymptotic
selectable-child transfer theorem is claimed.** This records a correlated
alternative after the actual independent-frame obstruction, preserving its
promise and its unpaid step without replacing the latter by new notation.

## 1. The original-value interface

For an actual child `A`, choosing the opposite diagonal child gives

```
L=[[A,C],[-C,-A]],       C skew hollow full signs,
Q(L)=2 max_(u,v disjoint signed supports)
                  (|u^T A v|+|u^T C v|).                 (1)
```

The supports partition the coordinates. Writing `A-C` as twice a signed
oriented adjacency matrix gives the equivalent weighted tournament
directed-cut norm. Both shores carry arbitrary signs, so ordinary
unweighted cut-flow discrepancy or an orientation balancing only vertex
degrees does not control (1). The earlier selectable anti-lift artifact
already preserves this exact interface and its finite optimizer search.

An asymptotically useful orientation theorem must pay the joint expression
in (1) for a selected exact child, at most `sqrt(2)M_n+o(n^(3/2))`.
The correlated subhalf anti-family constructed earlier does not impose
that child's optimality. Independent skew orientations alone have the
scoped lower bound

```
E Q(L)>=Q(A)+sqrt(2/pi)n^(3/2)-o(n^(3/2)),               (2)
```

using the same shared-edge field identity as the constant-two-frame
theorem, with the favorable polarity chosen in the diagonal children.
This can reject an iid orientation law near sufficiently low hypothetical
seed values; it does not reject correlated orientations.

## 2. A legal correlated law with automatic edge variance

Let `G` be skew with independent standard Gaussian upper-triangular
entries and put

```
K=A G+G A,       C_ij=sign(K_ij),       C_ji=-C_ij.       (3)
```

The raw matrix `K` is skew. For EVERY hollow full signing `A`,

```
Var(K_ij)=2(n-2)   for i!=j.                            (4)
```

Thus for `n>=3` the edge normalization is the SAME scalar
`sqrt(2(n-2))` regardless of global operator spikes. It is not a covariance
rescaled by `||A||op`, and it produces actual full skew signs almost surely.

For any real disjoint-support vectors `u,v`, the exact Gaussian response
variance is

```
Var(u^T K v)
 =||v||^2||Au||^2+||u||^2||Av||^2
        -2(u^T A v)^2+2(u^T A u)(v^T A v).              (5)
```

To derive it, write the response as
`(Au)^T G v+u^T G(Av)` and use
`E[(a^T G b)(c^T G d)]=(a.c)(b.d)-(a.d)(b.c)`.
The extra general-vector term is `-2(u.v)(Au.Av)`, which vanishes for
disjoint supports. Inserting distinct coordinate vectors gives (4),
since the two squared column norms are `n-1` and `A_ij^2=1`.

The negative cross-square in (5) is the relevant favorable feature.
For comparison, a skew COMMUTATOR driven by a symmetric GOE matrix has
the opposite signs on the two last terms and encourages the high child
cross-response in this variance test.

In the special flat case `A^2=(n-1)I`, write
`p=||u||^2`, `q=||v||^2`, `a=u^TAv`,
`b=H_A(u)+H_A(v)`, and `d=H_A(u)-H_A(v)`. Then

```
Var(u^T K v)/[2(n-2)]
 =[(n-1)pq-a^2+b^2-d^2]/(n-2),
|a|+|b|<=Q(A).                                         (6)
```

The last inequality follows from the two legal spins `u+v,u-v`.
Thus this is genuine joint geometry, rather than independent estimates
of a child cap and bridge norm.

## 3. Why this does not close the old Gaussian-parent obligation

At a flat child, the law in (3) is closely related to the earlier
opposite-child Gaussian covariance designs: it emphasizes the appropriate
same-polarity spectral pairs. That geometry is not being claimed as a new
solution. Equation (4) extends the legal normalization without a global
operator denominator, but it does not bound the whole covariance operator,
the high-order sign chaos, or the actual maximum in (1).

In particular Gaussian sign rounding does not preserve (5) verbatim.
Even after a valid whole-configuration Gaussian-sign comparison, the
favorable correlated Gaussian parent maximum still has to be bounded.
Crude union bounds using only (6) remain far above the required
`sqrt(2)M_n` joint threshold. Repeating covariance stationarity, or naming
that remaining maximum as a new kernel, would not be additional progress.

## 4. Bounded actual-cap experiment

The script
`computations/principle_synthesis_2026_09_07_skew_gaussian_law_check.py`
uses the archived actual optimal seed at each of `n=6,7,8`. For each seed
it generates 500 bridges under (3), iid skew signs, and the symmetric-GOE
commutator control. It exhausts all joint disjoint-cut tests in (1), so
each sampled parent cap is an exact integer, not a heuristic ground state.

| n | child cap | iid mean | anticommutator mean | commutator mean | smallest anticommutator cap |
|---|---:|---:|---:|---:|---:|
| 6 | 5 | 23.792 | 21.456 | 24.736 | 20 |
| 7 | 9 | 31.520 | 29.488 | 34.136 | 24 |
| 8 | 10 | 38.768 | 36.240 | 41.848 | 32 |

The law moves the finite responses in the favorable direction, but the
order-eight sample never reaches the unfilled target `Q(L)<=28`.
At order seven the cap-24 witnesses meet its finite unfilled doubling
target. These facts do not establish a scalable success or failure.

For every seed the script also checks (5) on 30 disjoint-support pairs
by summing the exact squared response over a complete Gaussian basis,
and checks (4) on every edge. All 90 response identities and all edge
variances PASS. The full cap histograms are preserved in
`computations/results/principle_synthesis_2026_09_07_skew_gaussian_law_check.json`.

The outstanding task is still an actual seed-sensitive joint orientation
estimate, allowing a selected optimizer and correlated or rare choices.
No new authority restriction, preferred covariance requirement, or
universal information-compression premise has been introduced.
