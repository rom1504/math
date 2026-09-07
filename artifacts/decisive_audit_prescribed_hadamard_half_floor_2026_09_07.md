# Independent audit: prescribed Hadamard stabilization half floor

Date: 2026-09-07. PASS for
`decisive_independent_prescribed_hadamard_half_floor_2026_09_07.md`.
This is an actual Boolean witness obstruction to a specified tensor
architecture, not a claim about the original optimal sequence.

1. With `H4=J4-2I4` and `z=(x,x,y,y)`, the same-block contributions
   cancel exactly, including their diagonals:
   `z^T(H4 tensor C)z=8x^TCy`. Thus q is at least four times beta.
   Together with `q<=beta/2`, and the normalized outer size `4s`, this
   proves `R=.5 sup_s beta(H_s tensor B)/s^(3/2)`.

2. Every prescribed outer is regular with positive row sum sqrt(s).
   Tensoring with its all-one vector makes normalized q and beta
   nondecreasing in both generator exponents. Their shifted-quadrant
   supremum is unchanged. This proves exactly
   `R(H_t tensor B)=t^(3/2)R(B)`; a one-sided tensor inequality alone
   would not have sufficed here. Finiteness follows from the spectral
   bound stated in the source.

3. If s>=n, `M=H_s^T E` is literally a sign matrix and
   `H_s M B^T=s E B^T`. Its output entrywise norm is
   `s sum|B_ij|`, so the bilinear witness is exact at the allowed basis.
   This does not substitute an arbitrary polar matrix or an unavailable
   Hadamard into the construction.

4. The two log generators are irrationally related. For each epsilon,
   finitely many multiples of log144 modulo log4 form an epsilon-net.
   At sufficiently large target log order, the remaining nonnegative
   log4 exponent supplies an allowed order between the target and its
   epsilon-enlargement. Thus `s_+(m)/m->1`. Applying the preceding witness
   to `H_t tensor B`, not just to B, gives
   `R(B)>=sqrt(t/s_+(nt)) sum|B_ij|/2`; the limit is exactly the claimed
   `sum|B_ij|/(2sqrt(n))`.

5. The prescribed H144 is indeed symmetric regular Hadamard: direct
   orthogonality gives `K^2=144I` and
   `K vec(F)=12 vec(F)`, so the displayed diagonal conjugation gives
   positive row sum12. No existence hypothesis beyond an order12
   Hadamard is added.

6. The Walsh-only restriction is correctly separated: arbitrary seed
   orders cannot use the two-generator density step. At n=4^r the
   allowed outer s=n suffices, and the H4 conversion lands at outer4n.
   If n>=4, a balanced sign diagonal D has trace zero. For
   `B=A+D`, q(B)=Q(A) exactly, and
   `tr(H_(4n) tensor B)=tr(H_(4n)) tr(B)=0`. Hollowing therefore does not
   change any Boolean quadratic energy. The witness is `4n^3` at order
   `N=4n^2`, exactly `N^(3/2)/2`, with no diagonal error.

The final fixed-seed cofinal lower bound and the adaptive choice of outer
depth for growing seeds follow from the exact regular monotonicity. They
do not assert a uniform depth bound in arbitrary seed order. In particular,
the argument excludes strict-subhalf cap transfer by this full prescribed
Hadamard stabilization, even for actual optimal hollow seeds completed
with signs, without proving the stronger real-seed identity R=T.
