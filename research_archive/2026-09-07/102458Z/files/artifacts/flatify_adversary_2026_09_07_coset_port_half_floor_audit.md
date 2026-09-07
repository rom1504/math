# Independent audit of the coset-port half floor after fibre replacement

2026-09-07. Section 6 of
`flatify_construct_2026_09_07_balanced_spread_audit.md` passes independent
reconstruction. The exact finite Walsh witness in
`decisive_independent_prescribed_hadamard_half_floor_2026_09_07.md`
was reread directly, rather than assuming its scope generalized.

For n=4^s and dyadic m>=2sqrt(n), the full coset-port matrix is exactly
F tensor W_m up to coordinate permutation, with seed F symmetric full
sign, including its diagonal. W_m is switching-equivalent to
(-1)^d H4^(tensor d), where m=2^d. The global sign is harmless for q.
The archive's same-spin witness appears already at outer order 4n,
and all further H4 factors preserve its normalized energy by regularity.
Thus its finite bound applies at outer order m^2, with no cofinal-depth
assumption, precisely on the claimed power-four seed orders.

The deletion of within-fibre blocks cannot simply be called norm-small.
Instead use the exact automorphism group. The kernel W_m is the sign of
the nondegenerate alternating form on w=(i,b) in F2^(2d). A symplectic
linear map preserves every kernel entry, so simultaneous permutation of
the physical w coordinate in every seed component preserves the full
matrix and the full witness energy. The group acts transitively on
nonzero vectors by extension to symplectic bases. For w!=w', their
transformed difference therefore lies in the d-dimensional fibre kernel
with probability (m-1)/(m^2-1)=1/(m+1).

The exceptional terms w=w' always lie within a fibre, including seed
off-diagonal terms and genuine matrix diagonal terms. Their sum is
(1/2) sum_w x_w^T F x_w, so its absolute value is at most m^2 q(F).
Choose the sign sigma of a full cap witness so sigma times its energy
is positive; no claim that reversing the witness changes quadratic sign
is needed. Averaging sigma times the cross-only energy gives

    (m/(m+1))*(q(F tensor W_m)-sigma*E_equal_w)
      >= (m/(m+1))*(n^(3/2)m^3/2-m^2 q(F)).

At least one literal Boolean witness attains that average lower bound.
This proves the cross-only cap inequality for the actual matrix, without
changing its blocks or substituting an arbitrary orthogonal catalyst.

Good internal completions have total cap at most C*m*(nm)^(3/2), with
universal C; subtracting this bound at the selected witness costs only
C/sqrt(m) after normalization by N^(3/2), N=nm^2. For F obtained from
an actual minimizer by any sign diagonal, q(F)<=M_n+n/2=O(n^(3/2)).
The remaining seed correction is O(1/m), uniformly for growing seeds.
Hence the completed architecture has normalized cap at least 1/2-o(1)
for growing n=4^s and m>=2sqrt(n), or for fixed such n and m tending to
infinity. There is no uniform-in-seed stabilization assumption hidden here.

The conclusion is restricted to this coset-port tensor architecture and
good internal fibre replacement. It does not establish the same floor
for odd-logarithm seed orders, arbitrary port codebooks, or unrestricted
global changes of the old internal and cross edges.
