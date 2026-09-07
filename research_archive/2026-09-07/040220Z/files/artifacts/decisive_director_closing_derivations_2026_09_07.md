# Director's final-hour derivations and falsified shortcuts

This is continuing research, not a claim that the original limit is proved.
It preserves unsuccessful arguments as well as checked deductions.

## 1. Width lower bound: direct reconstruction

For a hollow symmetric matrix define P=max H, R=max(-H), W=(P+R)/2.
Every cube mean has energy in [-R,P]. For the literal marked feasible
means mu_+=F+H sign(BF) and mu_-=-F+H sign(BF), their energy difference
is 2 sqrt(n-1) sum_i H_i |(BF)_i|. Thus the old certificate bounds W,
not just Q. Independent rounding and principal restriction preserve the
whole interval. The simultaneous majorant can use beta(A)<=4W directly:
x^T A y=2[H((x+y)/2)-H((x-y)/2)]. The fixed-operator theorem is universal,
so no original-minimizer hypothesis is silently transferred to width
minimizers. The fresh lower audit gives the unchanged c_*>.4333221116640807.

Consequently, for any sequence with Q(A_n)/n^(3/2)<=u+o(1),

    min(P(A_n),R(A_n))/n^(3/2) >= 2c_*-u-o(1),
    |P(A_n)-R(A_n)|/(2n^(3/2)) <= u-c_*+o(1).

At the new all-order cap upper this gives lower .3721290983 for BOTH
polarities of any such near-minimizing sequence. It does not make the
midpoint o(n^(3/2)), and so it does not identify minimum width with M_n.

## 2. Bipartite extension: exact normalization checked independently

For A=[[0,C],[C^T,0]], C a full m by m signing, B=A/sqrt(m), each row
has squared norm one, Qcap(A)=||C||_(infinity->1), and beta(A)=2Qcap(A).
Keeping the support masks is essential. Every surviving connected diagram
has uniquely forced shore coloring, giving m choices per free label.
The rooted isomorphism normalization and doubled-tree leaf sums are
unchanged; nonbipartite diagrams vanish. This is not an arbitrary weighted
extension of sign parity.

For the first nonlinear old tree the director independently recomputed
the Gram matrix by matching the spin set {j,k,l}. At root i in the first
shore its coefficient is sqrt(2) C_ij C_kj C_lj/m^(3/2), k<l, k,l!=i.
Hence its variance is (m-1)(m-2)/m^2. Distinct same-shore roots i,r have
covariance (m-2)(m-3)Q_ir/m^2; cross-shore covariance is exactly zero
because the two monomials have different shore counts. There is no hidden
factor two in the limiting Gaussian law.

The director read the complete bridge audit and checked the equal-shore
spectral deletion and final conversion: 2m vertices, normalization sqrt(m),
and the feasible energy difference yield ||C||_(infinity->1)>=
(2c_*-o(1))m^(3/2). This is an asymptotic theorem, not a contradiction to
small-order values. The transfer and bridge proof reconstructions agree.
The director replayed the transfer checker: 64 matrices, 4384 exact Gram
entries and 64 exact cap normalizations, all PASS. The second independent
checker is separately replayed; finite checks supplement the support proof.

## 3. Bipartite feedback centering and its nonclosure

For the old odd-degree rooted trees, each whole tree has an odd number of
vertices on each shore. Removing an external L root leaves an even number
of L input spins and an odd number of R input spins. Therefore an odd
response F_L has global-input parity (+,-), while F_R has (-,+). An even
mask H has (+,+). With the ODD convention sign(0)=0, T=H sign(BF) has the
reversed characters. Consequently

    E H_A(F)=E H_A(T)=0,
    E H_A(F+T)=sqrt(m) sum_i E H_i |(BF)_i|.

Both equalities are exact before any matrix limit. A hard +1 tie convention
would invalidate the parity statement; odd soft signs or symmetric tie
randomization also work. This explains an extra bipartite symmetry but
does not automatically permit another feedback pass: F+T has mixed
characters, so the centered subclass is not closed under that update.

## 4. Failed Pythagorean block shortcut

The proposed universal inequality

    W(A)^2 >= (W(A_S)+W(A_T))^2 + ||A_(S,T)||_(infinity->1)^2

is FALSE even at an actual order-four optimum. Take the two internal edge
signs +1 and -1, and all four cross signs +1. Writing u=x1+x2, v=x3+x4,
the energy is x1*x2-x3*x4+u*v. If both pairs agree the internal energies
cancel and the energy is +/-4; if exactly one agrees the energy is +/-2;
if neither agrees it is zero. Thus W=Q=4=M_4, whereas the proposed right
side is 2^2+4^2=20. This is a finite falsifier only, not a scalable
obstruction to a near-minimizer recurrence.

## 5. Literature and remaining limit obligation

A refreshed search of Boolean Bohnenblust--Hille/Sidon results found no
verified all-order flat-coefficient theorem. Projection constants have
an exact Hermite/CLT limit; their comparison with Sidon constants has a
fixed degree-dependent factor. The latter allows arbitrary coefficients,
not just full equal-magnitude coefficients. Monotonicity of a dimension-free
BH constant therefore cannot silently establish convergence of M_n.
Primary source inspected in the search:
https://ri.conicet.gov.ar/bitstream/handle/11336/257790/CONICET_Digital_Nro.ad85c82b-29b4-42ef-b8cf-cad2db078c4c_B.pdf?sequence=2
No theorem from this source is a dependency of a new bound here.

The strongest original-limit target still requires a quantitative actual
comparison, not a new name for a state: either vanishing-loss all-order
transfer or a summable original-value recurrence. The current exact width
interpolation leaves optimizer flip-cost imbalance and the width/absolute
gap. The favorable flatification route leaves fixed-oversaturation recovery.
None is proved by the preceding scoped extensions.
