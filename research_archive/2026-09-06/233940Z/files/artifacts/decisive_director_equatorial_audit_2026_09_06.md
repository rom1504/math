# Independent audit: equatorial switching

2026-09-06. The director independently reconstructed
`decisive_audit_equatorial_switching_2026_09_06.md`.
PASS for all real hollow symmetric matrices and its signing corollary.

The only Eulerian four-edge multiplicity types in the fourth-moment
expansion are one repeated edge, two doubled edges, and a simple four-cycle.
The Gaussian/Boolean difference is nonnegative for the first two types
(coefficients are squares); the cycle contributions agree, even when their
coefficient product is negative. Hence comparison with the Gaussian
quadratic form is legitimate for signed real coefficients. For
H_G=G^T A G/2, its fourth cumulant is 3 tr(A^4) and variance sigma^2,
giving E H_G^4<=15 sigma^4. No positivity assumption on A is used.

Interpolation gives E|H|>=sigma/sqrt(15); mean zero and Cauchy--Schwarz
give both strict energy signs probability at least 1/60. Along the path
from independent uniform X to independent uniform Y, the i-th field has
independent coordinates and coefficient-square sum at most v. A union
bound, not an independence assertion about different path fields, makes
the probability of a field above sqrt(2v log(7200n)) at most 1/3600.
Opposite endpoint signs have probability at least 1/1800. The difference
of these probabilities is positive, and a crossing of a jump of size at
most 2L has an endpoint of magnitude at most L.

For signs, gauge switching keeps Q exactly. Majority-edge flips bring
the gauged total S to 0 (even edge count) or +/-1 (odd edge count).
There are at most |S|/2 flips, each costing at most two in Q. Therefore

```math
0\le M_n^{\rm balanced}-M_n
\le\sqrt{2(n-1)\log(7200n)}.
```

The normalized error is O(sqrt(log n)/n), uniformly at every order.
This removes a balancing constraint, but neither asserts exact balanced
gauges nor relates optimal values at different orders. The original
convergence problem remains unchanged by this relaxation.
