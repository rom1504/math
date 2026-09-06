# Fixed-Hadamard regularization: primary scope checkpoint

Date: 2026-09-05. No general fixed-Hadamard tensor regularization theorem is
imported here. The odd-Walsh case is now settled separately in
`fresh_schmidt_odd_walsh_regularization_2026_09_05.md`.

## 1. New exact finite nonregularizable classes

H. Kharaghani, B. Tayfeh-Rezaie, and V. Zaitsev,
*On Regular Quaternary Hadamard Matrices*, J. Combinatorial Designs (2026),
[DOI](https://doi.org/10.1002/jcd.70006),
[author PDF](https://www.cs.uleth.ca/~hadi/research/regular-quat-2025-f-1.pdf),
Section 4: real Hadamard matrices of order 36 occur in eight inequivalent
classes with maximal excess 204 and eight with maximal excess 208. Since
`36^(3/2)=216`, none is equivalent to a regular matrix. Their algorithm
exhausts all column-sign assignments and takes the sum of absolute row sums;
this is exactly the bilinear norm `beta(H)`.

The displayed examples arise from general quaternary matrices and are **not
symmetric**. The paper asserts no symmetric representative and no preservation
of the gap under tensor amplification. Its
[author-hosted data](https://www.cs.uleth.ca/~hadi/research/28%20matrices%20and%20excess.txt)
contains 28 explicit real matrices with the published excess labels.
An independent integer check of all 28 displayed arrays during the
final scope audit found none equal to its transpose. This is a statement
about the displayed representatives, not about every equivalent matrix.

## 2. Bounded exact symmetry search

`computations/fresh_hadamard36_symmetrization_probe.py` downloads that data,
checks the Hadamard equations in integers, and asks whether a signed column
permutation makes a chosen matrix symmetric. This searches all independent
row/column equivalences relevant to symmetry: if `P H Q` is symmetric, then
its congruence by `P^T` is `H Q P`, so a column-only signed permutation suffices.

After normalizing the first row of `H` to all plus, write the selected column
permutation as `p`. Global sign lets the first column sign be one. Symmetry in
row zero forces the remaining column signs to be `d_j=H[j,p_0]`. The solver
therefore imposes exactly

`H[i,p_j] d_j = H[j,p_i] d_i` for every `i<j`,

plus the permutation constraint. This does not require regularity or a
constant diagonal. Every positive result is independently checked in integers.

OR-Tools CP-SAT 9.15.6755 returned `INFEASIBLE` for dataset indices zero
(published excess 204; 17.50 seconds) and one (208; 29.65 seconds), using
two workers and a 120-second cap. These are solver outcomes rather than
portable formal unsatisfiability certificates. No conclusion about all
16 equivalence classes is drawn. No symmetric cosquare counterexample was
obtained from these two displayed examples.

This finite source/probe is not used in the subsequent universal
actual-sign weighted theorem. The separate exact symmetric cosquare
counterexamples were constructed and certified internally and have
their own asymptotic scope audit; they are not inferred from these
order-36 source matrices or the CP-SAT outcomes.

## 3. XOR norm conventions and a recent near miss

The actual regularization is

`R(B)=0.5 sup_s beta(H_s tensor B)/s^(3/2)`.

It is not the ordinary quantum XOR bias `gamma_2^*(B)/2`: the already audited
five-cycle Seidel seed has `gamma_2^*(A)=5 sqrt(5)` but
`R(A)<=T(A)=12 sqrt(5)/5 < gamma_2^*(A)/2`. A possible norm identity must
instead concern the operator Hilbert-factorization ideal
`Gamma_2(B:ell_infinity -> ell_1)=2T(B)`, or a narrower seed class.

Cleve–Slofstra–Unger–Upadhyay,
[*Perfect Parallel Repetition Theorem for Quantum XOR Proof Systems*](https://cs.uwaterloo.ca/~cleve/pubs/2008PerfectParallelRep.pdf),
Theorem 1.3, proves quantum **parity**-game bias multiplicativity. Combined
with Grothendieck, it gives a self-tensor kth-root classical limit, not the
unrooted fixed-catalyst supremum above.

Andris Ambainis,
[*Optimal bounds on the classical value of the repeated CHSH game*](https://arxiv.org/pdf/2608.16439)
(August 2026), Theorem 2, proves strictness of a multiplicative relaxation
below quantum value whenever the original game has a classical–quantum gap.
Its definition on pages 1–3 is **AND repetition**: players return answer
tuples and must win every coordinate. This is not the matrix tensor/parity
objective defining `R`; its asymptotic strict gap cannot be imported as a
Hadamard excess gap.
