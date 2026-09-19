# Validation notes

The first finite replay stopped at an ancillary exact hard-ferromagnet
identity. The test had copied the special left-child formula
`beta+n/2+max_x[(B1).x-(sum x)^2/2]`, valid only for the all-negative
left child, while its test A was random. It was corrected to the general
identity `beta+max_x[H_A(x)+(B1).x]` for a hard-ferromagnetic right child.
This was a test setup error; the bounded-cap proof and its pointwise
inequalities use arbitrary A and were unchanged.

The second replay completed all finite checks, then stopped because a
symbolic test used structural equality rather than simplification of the
difference. Replacing that comparison by `simplify(lhs-rhs)==0` fixes
the representation-sensitive assertion. No mathematical formula changed.
