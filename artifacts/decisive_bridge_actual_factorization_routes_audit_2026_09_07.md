# Actual-signing factorization routes: checked obstructions and archive collisions

Date: 2026-09-07. These are failed routes, not a new limit theorem.

## Fractional Schur/purity route: archive collision

For a fractional B with entries in [-1,1], a vanishing quadratic-cap-norm
rounding error would require average sign defect sum(1-|b_ij|)/n² ->0.
Bipartition/polarization plus Khintchine gives the necessity; independent
biased rounding and Bernstein give sufficiency. The exact rounding bound
sqrt(2VL)+4L/3 with L=(n+2)log2 is ALREADY in
second_phase_independent_abstraction.md, section3, and
ar_matrix_rounding_literature_toolkit.md, section4. Their slab-body bridge
formulation already states the meaningful unsolved near-sign saturation
condition. This attempt is not new and was stopped after checking the archive.

## Full-response Cayley dilation: two exact restrictions

Encoding the whole finite seed response into a Cayley kernel on F_2^r does
not evade the tensor/flat-spin floor. For a symmetric signing
K(g,h)=f(g-h), every Fourier character is itself a Boolean spin vector.
Thus the Boolean quadratic maximum attains N times the largest absolute
convolution eigenvalue. Parseval for a +/-1 function f forces that eigenvalue
at least sqrt(N). Removing the diagonal costs only O(N), so the normalized
cap is asymptotically at least1/2. This holds regardless of how completely
f encodes the original seed.

Nonabelian groups avoid the first argument, because their higher-dimensional
matrix coefficients are not Boolean characters. But EVERY group-developed
Hadamard matrix still has a constant row sum r. Orthogonality applied to the
all-ones vector forces r²=N. Hence the all-ones spin already gives absolute
quadratic energy N^(3/2), before the factor1/2 and diagonal correction.
Nonabelian flat Fourier/Hadamard dilation therefore also has cap floor1/2.

An attempted matrix-unit realization through an extraspecial2-group has a
different obstruction. On its order-four elements g, inversion sends g to
zg, where z is the central involution. A function odd under z would isolate
the desired high-dimensional central-character representation, but symmetry
f(g^-1)=f(g) then forces f(g)=0 there. It is not a full +/-1 signing.
Complex phases or zero fills change the problem and cannot be silently used.

Balanced NON-Hadamard nonabelian Cayley kernels are not excluded by these
arguments. However no actual seed-controlled Boolean response inequality or
all-order realization was obtained for them. Merely replacing spectral norm
by a representation label does not supply the missing transfer theorem.
