# A positive support-adapted operation for the entire paired-noise sector

Scope update: the fixed-center energy-window application is now covered
without the paired constraint by
`flatify_adversary_2026_09_07_dephased_hadamard_energy_window_audit.md`.
The exact restricted-operator construction below is still correct.

2026-09-07. The independent track observed that active/inactive column
pairing makes the bridge contractive on paired-noise subspaces. This note
extends that observation to an exact center-zero operation controlling
ALL common noise levels in that sector, with no Gaussian exception.

Let k=2^d>=4, m=k/2, n=mk, and take arbitrary actual children with chosen
positive ground centers x0,y0. In each fibre use the Walsh frame row-signed
by its center. Fix a nonzero translation a of the physical Walsh domain.
The states x=x0*eta with eta(u)=eta(u+a) form a Boolean subspace sector;
their normalized features are supported on the m active frequencies
j dot a=0. The other m frequencies are identically zero.

Pair every active column with an inactive column, putting the active
column in bit zero of its rank-two pair. Label the resulting m pairs by
opposite fibres. On the left assign the center's unique active frequency
zero to port j=i. On the right assign its center frequency zero to port
i=j+1 modulo m. Complete all other labels bijectively. These assignments
are possible for every m>=2.

Build the usual literal rank-two sign bridge C from these frames. On
active feature coordinates, K has only its active-active coefficient
1/sqrt(2) in each reciprocal pair. Thus

    P_L K P_R=(1/sqrt(2)) times an orthogonal permutation

between the two n/2-dimensional active spaces. The center supports never
meet in a reciprocal pair, so x0^T C y0=0 EXACTLY. No random signs or
probabilistic center selection are needed. All bridge entries remain
signs, and the full operator is still orthogonal after normalization.

`computations/flatify_adversary_2026_09_07_support_adapted_bridge_check.py`
and its result JSON preserve literal n=8,32 matrices, verify the center
zero and restricted orthogonality identities exactly, and exhaust all
104 equal-absolute-overlap paired-sector pairs at n=8 against the angular
bound below. These finite tests supplement, not replace, the general proof.

## Exact common-overlap bound

For paired-sector states with equal absolute overlaps alpha with their
centers, flip either full state if necessary to make both overlaps
nonnegative; internal energies and the absolute bridge energy are
unchanged. The restricted orthogonal map sends the two normalized
centers to orthogonal vectors. Spherical triangle inequality therefore
gives

    |x^T C y|/n^(3/2) <= b_pair(alpha),

    b_pair(alpha)=1/sqrt(2),                 0<=alpha<=1/sqrt(2),
                  sqrt(2)*alpha*sqrt(1-alpha^2), 1/sqrt(2)<=alpha<=1.

Indeed each candidate is at angle arccos(alpha) from its respective
center. When twice that angle is below pi/2, the two candidate directions
remain separated by at least pi/2-2arccos(alpha); otherwise use the norm
bound. Applying the same argument to the negative of one candidate
handles the absolute value. This is NOT the minimum of the two displayed
expressions on the whole interval; that false shorthand would fail near
alpha=0.

The exact maximum of b_pair(alpha)/[2(sqrt(2)-alpha^2)] is

    c_pair=sqrt(2+sqrt(2))/4 = .4619397662... .

For alpha^2<=1/2 the ratio is increasing. On the other branch its square
is v(1-v)/[2(sqrt(2)-v)^2], v=alpha^2, whose unique maximum occurs at
v=sqrt(2)/(2sqrt(2)-1) and has value (2+sqrt(2))/16.

Consequently for any c>c_pair, in ALL child-energy windows satisfying
|H_A(x)+H_D(y)|<=2c alpha^2 n^(3/2)+o(n^(3/2)), this one deterministic
actual sign bridge satisfies a strict uniform parent saving below
2sqrt(2)c n^(3/2), simultaneously for every common overlap alpha and every
paired-sector pair. In particular c=47/100 works. No exceptional
Gaussian-profile set remains inside this constrained sector.

Independent paired noise about actual optimizer centers gives the stated
child-energy windows with probability 1-o(1), by the compressed-field
variance proof in `flatify_adversary_2026_09_07_paired_noise_exceptional_entropy.md`
with every fibre paired. That typicality is not needed for the deterministic
bridge bound itself, which covers the full constrained sector.

## Scope

This is a genuine favorable change of the dense bridge construction, not
an inference that exceptional profiles can be ignored. It controls a
proper Boolean subspace sector, not arbitrary mixtures of constrained
and unconstrained fibres, unequal overlaps, or all parent spins. It does
not establish recurrence or exclude any broader operation which also
changes internal edges. The compatible-order statement above is exact;
padding/restriction can give an inherited all-order sector with vanishing
errors, but should not be confused with the unpadded paired constraint.
