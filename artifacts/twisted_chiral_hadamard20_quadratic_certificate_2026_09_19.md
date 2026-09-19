# Compact quadratic-form certificate for the order20 witness

This is an alternate proof of the upper bound in the canonical
[Hadamard20 note](twisted_chiral_hadamard20_2026_09_19.md), together with
a standalone exact witness check. It does not assert global M20=40.

For a symmetric order20 Hadamard H and a sign vector z, put
r_i=z_i(Hz)_i and S=sum_i r_i. Orthogonality means any two rows differ
in10 positions, so all r_i share one residue0 or2 modulo4. Also
sum r_i^2=400. If the residue is2, summing
(r_i-2)(r_i-6)>=0 gives S<=80. If it is0, summing
(r_i-4)(r_i-8)>=0 gives S<=86+2/3; since S is divisible by4, S<=84.

To exclude84, set w=(Hz-4z)/4. This is an integer vector satisfying

    ||w||_2^2=3,  z^T w=1,  Hw=z-4w.

Thus w has exactly three nonzero entries, each +/-1, and
w^T H w=1-4*3=-11. But every entry of H has absolute value1, so
|w^T H w|<=||w||_1^2=9, a contradiction. Applying the argument to -H
proves |z^T H z|<=80.

For the retained conference child A10, use d=all+, B=A, and

    D=[[A,A+I],[A+I,-A]],    H=D+diag(-I10,+I10).

Then H^2=20I and traceH=0. Hence the hollow parent's energy equals
z^T H z/2 and has cap at most40. The explicit spin

    (1,-1,-1,-1,1,1,1,1,1,1, -1,-1,-1,1,1,1,1,1,1,1)

attains40. The switched H row sums are

    (6,6,6,6,2,2,6,2,2,6, 2,2,2,2,6,6,2,6,6,2),

namely ten2s and ten6s. Thus this exact cap certificate needs no parent
spin enumeration. The independent exhaustive checks remain separately
retained as corroboration, not a dependency of this algebraic argument.

Run the standard-library-only verifier:

    python computations/twisted_chiral_hadamard20_certificate_2026_09_19.py

It embeds the10-by10 conference child, matching, and spin, and checks
A^2=9I, H^2=20I, symmetry/full signs, trace, row sums and attained energy.
Its output is retained in
`computations/results/twisted_chiral_hadamard20_certificate_2026_09_19.json`.
