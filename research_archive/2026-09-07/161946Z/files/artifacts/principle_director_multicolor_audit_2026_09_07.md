# Independent audit: magnitude-first hard compiler

2026-09-07. Root reconstruction of
`principle_invent_2026_09_07_multicolor_hard_compiler.md`: PASS for its fixed
alphabet, fixed positive masses, and stated integrality hypotheses.

The proof first repairs undirected MAGNITUDE colors, not signed incidences.
Its simultaneous common-neighbor event concerns each ordered pair of colors
and each pair of distinct endpoints. Chernoff and a polynomial union bound
give the claimed linear candidate counts. A two-edge color exchange fixes
opposite degree discrepancies and leaves its mediator degree unchanged.
A three-edge alternating path fixes two same-sign units and leaves both
internal degrees unchanged; coincident endpoints use a triangle and change
that endpoint degree by two. The even sum of discrepancies justifies this
last operation. Previously completed colors are never touched.

Mediator loads can be capped at m^(3/4)/2. At most
O(m^(3/4)sqrt(log m)) vertices then become ineligible, while every initial
candidate family has linear size. Endpoint loads across the fixed number
of colors are O(sqrt(m log m)), since each corresponding discrepancy
decreases monotonically. Hence total incident changes at each vertex stay
below m^(3/4) eventually. The common-neighbor tests survive all repairs.
The same-color tests imply each final positive-density color graph is
connected, not merely regular. Its even target degree is forced by the
symmetric signed type. All these conditions are needed by the signed
Eulerian count.

The repair's probability-ratio cost is exponential in its edit count, and
its Hamming preimages cost exp[O(m^(3/2)log^(3/2)m)]. Exact degree colorings
all have the same independent-color probability exp[-eH(p)], so the entropy
count direction is correct. Signs within each nonzero color graph then
have the required exact local balance. One seed-parity violation per color
costs exactly 4t a_r^2; magnitudes are unchanged. Different colorings are
recovered from the output magnitudes, so their counts do not overlap.

For the upper comparison, a coupling with mismatch probability delta has
H(Y|X)<=h(delta)+delta log(L-1), and transport energy at most
-t gap^2 delta. This proves the displayed canonical upper bound. The
positive-kernel upper identity pays only the original type denominator;
the large condition number is NOT inserted into the old repair error.
Thus t tending to infinity with t=o(m^2) has the claimed uniform leading
pressure. This does not extend to arbitrary growing alphabets or moving
small type masses, and does not identify the exact parent Boolean cap.

The saved finite checker at m=257,513 is useful path bookkeeping, not a
proof of the asymptotic entropy estimate or of an original minimax claim.
