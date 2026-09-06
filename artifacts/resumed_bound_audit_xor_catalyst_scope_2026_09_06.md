# Focused primary-source audit: the stabilized XOR-game interpretation

Date: 2026-09-06. Status: exact normalization and a concrete separation
of the quantum target from the required cut-majorant target. No primary
theorem found that realizes the latter with the prescribed catalysts.

## 1. Exact game normalization

For a real matrix C write

    beta_count(C)=max_(x,y Boolean) |x^T C y|.

For a nonzero finite symmetric seed B, let Z_B=sum_(i,j)|B_ij|. Its
normalized XOR game has classical bias beta_count(B)/Z_B. A regular
Hadamard H_m has beta_count(H_m)=m^(3/2), by the all-one witness and
the spectral upper bound. Its uniform-question game therefore has
bias 1/sqrt(m).

The exact bilinear-lift identity for the prescribed regularization is

    2R(B)=sup_(allowed m) beta_count(H_m tensor B)/m^(3/2).

Equivalently, R(B) is Z_B/2 times the supremum of the classical bias of
the XOR product game, divided by the Hadamard game's own bias. It is
indeed a relative classical-bias stabilization.

The prefactor is 1/(2m^(3/2)) for the ORDINARY counting-space
infinity-to-one matrix norm. A prefactor 1/(2sqrt(m)) is equivalent
only when an outer probability normalization 1/m is already built
into that norm. This distinction is an entire factor m, not notation
that can be left implicit.

## 2. The closest available primary theorem has a different target

Cleve, Slofstra, Unger, and Upadhyay prove multiplicativity of the
quantum bias under XOR product in Theorem 1.3 of
[Perfect parallel repetition theorem for quantum XOR proof systems](https://cs.uwaterloo.ca/~cleve/pubs/2008PerfectParallelRep.pdf).
The statement is about optimal vector/entangled strategies, not
classical strategies with prescribed Hadamard catalysts. Their
parallel conjunction theorem is a separate assertion and must not
be confused with the XOR product used here.

Together with real Grothendieck, this yields the elementary consequence

    beta_q(B)^r/K_G <= beta_count(B^(tensor r)) <= beta_q(B)^r,

and hence the r-th-root classical asymptotic tensor-power value is
beta_q(B). But the seed itself is being tensor-powered in that limit;
it is not fixed while only the allowed outer grows. Neither that
corollary nor bias multiplicativity is the missing R=T theorem.

## 3. A three-point example strictly separates the targets

Take B=J3-2I3. Its standard vector/quantum bilinear value is EXACTLY 6.
For the lower bound, choose three unit vectors u_i with sum zero and
put v_i=-u_i. The score is

    sum_(i,j) (1-2delta_ij)<u_i,-u_j>=6.

For the upper bound, ||B||op=2 and the two vector matrices each have
squared Frobenius norm 3, so their bilinear pairing is at most 6.

In contrast, the already independently replayed three-cell certificate
gives

    2R(B)=2T(B)=17/3 < 6.

The T upper bound follows from the feasible absolute PSD majorant
|B|=2I-J/3: its Boolean cube maximum is 17/3. The matching lower value
uses the cut covariance B^2/3 and equals Tr|B|^3/3=17/3.

Thus this is an exact counterexample to a claim that the prescribed
Hadamard stabilization universally attains the STANDARD quantum XOR
value. There is no contradiction with quantum multiplicativity or
Grothendieck; those results never assert this stronger stabilization.

Geometrically, the required T dual optimizes a common covariance in
the cut polytope. Quantum vector strategies allow general correlation
matrices. The triangle's zero-sum vector Gram has off-diagonal entries
-1/2 and is not a cut covariance, exposing the distinction explicitly.

## 4. What remains open after the bounded search

The focused primary-source search did not locate a theorem for this
fixed-seed, prescribed regular-Hadamard classical stabilization that
identifies its value with the cut-PSD-majorant T. General quantum
activation, quantum embezzlement, and tensor-power regularization do
not supply classical Boolean signed modules with the required shared
label geometry. Nor is a Grothendieck gap automatically preserved
after tensoring: that would require its own uniform certificate.

The positive binary-phase and controlled-translation modules are
therefore concrete extra realizations, but they still need a genuine
lossless composition or matching theorem. The quantum analogy by
itself neither gives that theorem nor refutes it.
