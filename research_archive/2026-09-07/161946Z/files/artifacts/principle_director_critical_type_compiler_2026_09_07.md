# Critical binary type compilation — live proof attempt

2026-09-07. Status: EXPLORATORY; no unproved counting bound is a theorem.

The successful upper family averages a dense positive-kernel network subject
to row types. The new finite-type theorem shows that bounded kernels and
finite local state erase all low-cap seed distinctions at the leading m^2
pressure scale, even without independent output column phases. In the binary
Gaussian kernel, the canonical seed contrast is bounded by
O(t B^2 beta(S-T)); hence retaining a leading contrast requires t of order
sqrt(m) for low-cap seeds. The elementary independent repair bound then
becomes leading. This is a concrete critical-scale problem, not a new name
for the original maximum.

## Hard binary endpoint

Let m be odd, d=m-1 even, e=m(m-1)/2. On each edge i<j choose a sign a_e,
and give its two directed incidences the signs u_ij=a_e and u_ji=S_ij a_e.
Ask that every row contain d/2 plus signs and d/2 minus signs.

Necessary condition: multiplying all row products gives

    product_edges S_ij = (-1)^e.

This is also sufficient. Choose one Euler tour of the complete graph.
Pair successive halfedges at each visit. Require paired halfedges at a
vertex to have opposite signs, and require halfedge products on edges to
be S_ij. Propagate a sign around the tour. Its consistency is precisely
product_edges(-S_ij)=1. Each vertex then has balanced signs. This is an
actual construction, though not yet a sufficiently numerous family.

If parity fails, change one S edge for purposes of this incidence
assignment: all but that one edge can still obey the desired relation.
At a soft Gaussian kernel this costs O(t), not O(m^2).

## Counting identity to investigate

Pair all incident halfedges at each vertex arbitrarily. The resulting
transition system decomposes all edges into circuits. Let h(C) be the
product of -S on circuit C. A transition system admits a consistent
balanced incidence assignment iff every h(C)=1; in that case it has
2^(number of circuits) assignments. Every balanced assignment is compatible
with exactly ((d/2)!)^m transition systems. Therefore

    count(S) ((d/2)!)^m
       = sum_transition_systems product_circuits (1+h(C)).

The target is a uniform count lower bound of the form

    count(S) >= 2^e exp[-O(m log^c m)]

on the parity-compatible class, or the corresponding version with a
sublinear number of mismatching edges. This would pay critical-temperature
exact-type compilation far more efficiently than modifying O(m^1.5) random
incidences. It is not currently proved here. A finite zero count beyond
the stated parity condition would falsify the elementary existence proof.

Potential primary import: Schrijver (1983), Bounds on the number of Eulerian
orientations, DOI10.1007/BF02579193, and Borbenyi--Csikvari arXiv:1905.06215.
Their ordinary-orientation count is not silently assumed for signed
incidence constraints. The applicability of a signed extension is open.
