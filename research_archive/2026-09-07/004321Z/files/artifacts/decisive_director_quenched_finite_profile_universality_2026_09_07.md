# Quenched finite-profile universality for cap-bounded outer signings

Date: 2026-09-07. Status: proved finite-alphabet limit, independent audit
requested. This is a scoped obstruction to seed-sensitive bounded-profile
pressure, NOT a limit theorem for the original M_n.

## Statement

Use the finite classes, row color counts and proportions of
`decisive_director_positive_transport_limit_2026_09_06.md`. Replace its kernel
by two fixed strictly positive symmetric matrices K_+,K_-. At order n let
S_n be ANY hollow signing with Q(S_n)=o(n^2). Define the QUENCHED quantity

```math
Z_n(S_n)=\mathbb E\prod_{i<j}K_{S_{ij}}(X_{ij},X_{ji}).
```

The expectation is over independent uniform row arrangements, NOT over S.
Then its normalized logarithm converges to a value independent of S_n:

```math
\lim n^{-2}\log Z_n(S_n)
=\frac14\max_\gamma\sum_{a,b}\pi_a\pi_b
 \sum_{e=\pm1}[H(\gamma_{ab}^e)+\langle\gamma_{ab}^e,\log K_e\rangle]
-\sum_a\pi_a H(\nu_a).                                  \tag{1}
```

Here gamma_ab^e are probability matrices, gamma_ba^e=(gamma_ab^e)^T, with

```math
\frac12\sum_b\pi_b\sum_{e=\pm1}\sum_d\gamma_{ab}^e(c,d)=\nu_a(c).
```

The assignment of vertices to the finitely many classes can depend on S_n.
The result applies, in particular, to ACTUAL optimizing outer signings:
the elementary random-sign bound already gives M_n=O(n^(3/2)).

## Proof and where cap control enters

The elementary polarization bound gives
`||S_n||_(infinity->1)<=4Q(S_n)=o(n^2)`; it applies to all bounded vectors.
Hence every pair of vertex classes has asymptotically half positive and
half negative edges, with error o(n^2), uniformly in class assignment.
Diagonal slots contribute only O(n).

The entropy upper bound in the companion proof now averages edge pair laws
separately by ordered class pair AND edge sign. Compactness gives (1)'s
upper bound and its averaged constraints.

For the lower bound choose any feasible gamma and independently draw each
edge color pair from gamma_ab^(Sij). At vertex i of class a, expected color
counts differ from their balanced-class prediction by terms bounded in sum
over all i and colors by a fixed multiple of

```math
\sum_b\sum_i\left|\sum_{j\in V_b}S_{ij}\right|
\le r\|S_n\|_{\infty\to1}=o(n^2).
```

Errors in limiting proportions/profile counts add o(n^2). Independent-edge
concentration adds O(n^(3/2)sqrt(log n)) total count error with probability
tending to one. Recolor endpoints row by row, as in the companion theorem.
Only o(n^2) incidences change; fixed positivity of both kernels bounds the
log-weight loss, and the repair preimage multiplicity is exp(o(n^2)).
The usual independent-edge information and weight concentration complete
the lower bound. No same-temperature Gibbs projectivity is assumed.

## Exact centrally symmetric homogeneous simplification

Suppose r=1, the color alphabet has an involution c->bar c, nu is invariant,
and

```math
K_-(c,d)=K_+(c,\bar d),\qquad
K_+(\bar c,\bar d)=K_+(c,d).
```

Then (1) equals

```math
\frac12\max_{\pi=\pi^T,\,\pi1=\nu}
 [\langle\pi,\log K_+\rangle-D(\pi\Vert\nu\otimes\nu)].   \tag{2}
```

To verify the potentially delicate marginal step, view a feasible pair as
the joint law of (e,X,Y) with e fair and only the UNCONDITIONAL X,Y marginals
fixed to nu. Average it under (e,X,Y)->(-e,bar X,Y) and
(e,X,Y)->(-e,X,bar Y), and their products. Costs and constraints are preserved;
entropy cannot decrease. After this averaging each conditional X,Y marginal
given e is nu. Each sign sector is therefore bounded by the same self-transport
optimum, after reflecting Y for the negative sector. Conversely choose an
optimal positive-sector coupling, symmetrize under swap and simultaneous
reflection, and take its Y-reflection in the negative sector. This is feasible
and attains (2). Fixed points of the involution cause no difficulty.

For real bounded colors, K_e(x,y)=exp[-t(x-e y)^2] satisfies these hypotheses.
Thus every cap-bounded outer signing gives exactly the same leading pressure
on such homogeneous symmetric profiles, even though no outer signs are averaged.

## Crucial limitation: energetic spikes

This theorem holds for a fixed finite bounded alphabet, then the order limit.
It does not authorize truncating arbitrarily large row spectra with small
TWO-SIDED pressure loss. O(n) coordinates of magnitude sqrt(n) can change
log weight by Theta(n^2). Deleting those interactions is an upper relaxation,
not an asymptotically exact approximation. Seed information invisible in (1)
can reside in precisely that sparse energetic compatibility. The theorem
therefore narrows a particular construction mechanism; it does not say that
all optimizing signings have the same extremal value.
