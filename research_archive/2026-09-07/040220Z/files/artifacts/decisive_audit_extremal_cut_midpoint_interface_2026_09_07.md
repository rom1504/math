# Exact extremal-cut geometry and what midpoint recovery must change

For a hollow symmetric signing A let P=max H_A, R=max(-H_A),
w=(P+R)/2 and I=(P-R)/2. Replacing A by -A if necessary gives I>=0.
Choose exact top and bottom spins, and switch the top spin to all ones.
Partition the vertices by the relative sign of the bottom spin. Write

    A=[[B,C],[C^T,D]],     h(u,v)=H_B(u)+H_D(v).

Then, exactly,

    h(1,1)=I,       1^T C 1=w,       beta(C)=w,              (1)

and for every Boolean u,v,

    |u^T C v|+|h(u,v)-I|<=w.                               (2)

Indeed the chosen endpoint energies are I+w and I-w. For an arbitrary
u,v, flipping one whole shore gives the two energies h(u,v) +/- u^T C v.
Both belong to [I-w,I+w]. This is equivalent to (2), which also proves
beta(C)<=w; the chosen endpoint pair gives equality in (1).

There is no optimizing-matrix hypothesis in this identity. It applies
in particular to a global width minimizer, but global minimality has not
yet supplied a new filling inequality beyond it.

Erasing the two internal blocks produces the weighted bipartite matrix

    A_cross=[[0,C],[C^T,0]],      Q(A_cross)=w.               (3)

Thus an apparent free midpoint repair has simply left all internal
edges unsigned. Filling those entries with signs while keeping absolute
cap w+o(n^(3/2)) is a substantive missing theorem. Equation (2) certifies
the EXISTING internal energies centered at I, not a centered-at-zero
filling. A random independent filling has leading-order fluctuations;
it cannot be omitted from the proof.

The other exact formal repair is

    A_center=A-(2I/n) Id,       Q(A_center)=w,               (4)

because its quadratic energy is H_A-I on Boolean spins. But (4) uses a
forbidden diagonal. When I has order n^(3/2), each new diagonal entry has
order sqrt(n), and the total constant correction is of leading order.
Allowing only sign diagonals could shift energy by at most n/2.

Finally, adding only o(n) auxiliary vertices while preserving the old
principal signing cannot fix a leading midpoint. If A is a principal
submatrix of a full signing A' of order n+o(n), conditional expectation
over new independent fair spins gives

    Q(A')>=Q(A)=w+I.                                       (5)

After normalization by the new order, a fixed positive I/n^(3/2)
survives. Any successful such repair must modify old coefficients,
not merely append a gadget.

These observations isolate the width-to-cap obligation but do not
settle it. An actual non-cut sparse-sign surgery under a low-rank
high-state hypothesis is proved separately in
`decisive_audit_low_rank_midpoint_sign_surgery_2026_09_07.md`.
