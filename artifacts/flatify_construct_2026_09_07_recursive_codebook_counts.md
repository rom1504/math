# Recursive Hadamard codebook suppression: proved counts beyond fixed Walsh frames

Date: 2026-09-07. The exact small flat-count identity has been exhaustively
verified. The asymptotic statements below are proofs, not extrapolations.

## 1. Subexponential annealed fully-flat count

Let k=2^d with d odd and D=2k. Form the full Hadamard

    F_D=[F1 F2; F1 -F2],

where F1,F2 are Walsh Hadamards of order k with independent uniform signed
ROW permutations. Additional column gauges do not affect the counts.
Let P_k be the number of 1-plateaued Boolean functions at Walsh order k,
and B(F_D) the number of Boolean x with |F_D^T x|=sqrt(D) coordinatewise.
Then the exact identity is

    E B(F_D) = P_k^2 / binom(k,k/2).                     (1)

Indeed, for x=(xL,xR), put z+=(xL+xR)/sqrt(2), z-=(xL-xR)/sqrt(2).
If the parent output is fully flat then each child output is Boolean,
so both child input energies must equal k. Thus both supports have size
k/2. There are binom(k,k/2)2^k parent Boolean inputs with this agreement
pattern. A fixed child ternary word, after a uniform signed permutation,
is uniform among binom(k,k/2)2^(k/2) such words. Exactly P_k of them have
Boolean normalized Walsh output, by inversion and the definition of
1-plateaued. Independence of the two child gauges gives (1).

Potapov, https://arxiv.org/html/2303.16547, Corollary1, supplies the exact
degree bound deg f<=(d+1)/2 for 1-plateaued functions. Hence

 log2 P_k <= sum_(j<= (d+1)/2) binom(d,j)
           = k/2+binom(d,(d+1)/2).

Using the middle binomial estimate in (1),

    log E B(F_D) <= O(D/sqrt(log D)) = o(D).             (2)

This is a genuinely stronger count for the randomized actual Hadamard
parent than merely inserting a fixed Walsh flat-row count. No bound for
approximately flat spectra has been imported.

### Exact replay

`computations/flatify_construct_2026_09_07_recursive_flat_count.py` checks
all relative permutations at child orders2 and8. At k2, P_k=4 and every
parent has8 fully-flat inputs. At k8, P_k=112; the40320 permutations give

    parent flat count 896:1344 permutations;
                      384:9408;
                      128:18816;
                        0:10752.

The mean is896/5=112^2/binom(8,4). All65536 parent Boolean inputs were
directly checked for one representative of every count. Thus some literal
order16 frames in this family have no fully-flat inputs at all.

## 2. Actual consequence for the rank-two all-fully-flat sector

Use independent copies of these order-D frames in m=D/2 fibres, with
fresh independent output column signs. N=mD. For any fixed fully-flat
row tuple, every cross edge contributes +/-sqrt(2) to H_cross/sqrt(N),
independently across edges. The expected number of such tuples is
(E B(F_D))^m=exp(o(N)), by (2). Hoeffding and a union bound give, for
any fixed epsilon>0,

 P(exists fully-flat row tuple with |H_cross|>epsilon N^(3/2))
 <= 2 exp[o(N)-epsilon^2 N],                             (3)

with the harmless finite factor N/(N-D) available if desired. Thus some
actual realization has restricted cap o(N^(3/2)) on this ENTIRE sector.
The explicit rate from (2) is O(N^(3/2)/(log N)^(1/4)). Completing fibres
adds O(N^(5/4)). The proof averages frames and gauges jointly; it does not
assume one frame simultaneously realizes every annealed count.

This improves the deterministic1/(2sqrt(2)) bound on fully-flat rows,
but does not control non-flat row tuples or prove a new global M_n bound.

## 3. General prescribed-profile overlap formula

For an arbitrary child magnitude histogram rho of length k, define
T_k(J,rho) to count words u in {0,+/-1}^k with support J for which
sqrt(2/k) H_k^T u has magnitude histogram rho. Thus the child input is
sqrt(2)u. If cbar_D(rho) counts parent Boolean rows of prescribed histogram
rho in expectation, then EXACTLY

 cbar_D(rho)=sum_(J=0)^k 1/binom(k,J)
     sum_(rho+,rho-): (rho++rho-)/2=rho
           T_k(J,rho+) T_k(k-J,rho-).                   (4)

The factor1/binom(k,J) is a real overlap suppression, not an omitted
normalization: it comes from the parent input count divided by the two
child signed-support orbit sizes. Child energies automatically force
m2(rho+)=2J/k and m2(rho-)=2(k-J)/k.

The number of possible child histograms is exp(O(k^(2/3))), by the same
integer-square partition argument as before. Consequently (4) implies

 log cbar_D(rho)
 <= max_[J,rho+,rho-] {log T_k(J,rho+)+log T_k(k-J,rho-)
                       -log binom(k,J)}+o(k).           (5)

Equation(1) is the case where both child outputs are Boolean and only
J=k/2 is possible. Formula(5) is useful whenever independent bounds for
the two sparse transformed codebooks are available. It does not replace
those bounds by the unrestricted number of input words.

## 4. A further uniform suppression for the integral-spectrum sector

There is also a nontrivial bound covering multiple output magnitudes.
Let C(F_D) count Boolean parent rows whose NORMALIZED output coordinates
are all integers. This sector includes dyadic integer amplitudes, not
just one-level spectra. Then

    E C(F_D) <= 2^((c_int+o(1))D),
    c_int=(1+p0)/2=.5550139...,
    h2(p0)=1/2, 0<p0<1/2.                              (6)

Here is a direct proof. Put a=(d-1)/2 and M=2^a. Let L_J count ternary
words u of support J with H_k u divisible by M coordinatewise. Parent
integrality is exactly this condition in both children. Therefore

    E C(F_D)=sum_J L_J L_(k-J)/binom(k,J).               (7)

Elementary modular lemma: if an integer word u on F_2^d has Walsh transform
divisible by2^a, then u mod2 has algebraic degree at most d-a. To see this,
the sum on any coordinate face of dimension r>=d-a+1 is divisible by2;
the parity of that sum is the corresponding Boolean Mobius coefficient.

Thus the support indicator s=u mod2 belongs to the Reed--Muller space of
degree<=d-a. For a fixed support S, write u=s-2n. The difference of two
possible negative-part indicators n has Walsh transform divisible by
2^(a-1), so its mod2 reduction has degree<=d-a+1. Consequently

 L_J <= min{binom(k,J),2^B0} * 2^min(J,B1),
 B0=sum_(j<=d-a)binom(d,j), B1=sum_(j<=d-a+1)binom(d,j).

Both B0 and B1 equal k/2+o(k). Inserting this in (7), writing p=J/k and
using symmetry p<=1/2, yields exponent per k bounded by

    2 min{h2(p),1/2}+1/2+p-h2(p).

Its maximum is1+p0: below p0 the expression is1/2+p+h2(p), which increases;
above p0 it is3/2+p-h2(p), a convex function whose larger endpoint value
is1+p0. Divide by D=2k to obtain (6). Endpoints J=0,k cause no exception;
the finite binomial entropy approximation is uniform.

The bound (6) is genuine entropy suppression of a multiple-magnitude
sector, but it does NOT yet annihilate every dangerous resonant profile
in the joint transport certificate. In particular, low-entropy sparse
profiles can also have small transport cost. Arbitrary fine-lattice or
Gaussian-like output profiles remain outside this integrality theorem.
