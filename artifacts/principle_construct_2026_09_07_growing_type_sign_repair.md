# Quantitative sign repair for growing, rare-atom row profiles

2026-09-07. **Proved derivation; audit requested.** Fixed finite-type seed-pressure universality is already proved in `principle_invent_2026_09_07_finite_type_seed_universality.md`; the homogeneous symmetric hard compiler is in `principle_invent_2026_09_07_multicolor_hard_compiler.md`. No novelty is claimed for those finite settings. The purpose here is an explicit error bound uniform in rare atom masses, growing alphabet sizes, and EVERY macro signing. It applies to the actual raw-kernel permutation partition in the new fixed-seed mixed-orbit weave, not to its entire spin optimization.

## 1. Profiles and exact quantities

There are m vertices and D=m(m-1) directed ports. At row i prescribe a real multiset of size m-1, with at most d>=1 positive magnitudes a and |a|<=V. The all-zero case is immediate. The profiles may differ between rows. Write n_(i,a)^+, n_(i,a)^- for the two signed counts, n_(i,a) for their sum, and ignore zero coordinates when assigning signs. Define

    R = sum_(i,a) sqrt(n_(i,a)),
    A_asym = sum_(i,a) |n_(i,a)^+ - n_(i,a)^-|,
    B = ceil(R+A_asym).

In particular R<=m sqrt(d(m-1)). No lower bound on any nonzero atom mass is assumed. Let each row be an independent uniformly permuted copy of its prescribed signed multiset. For any symmetric macro signing S define

    Z_S(t)=E product_(i<j) exp[-t(u_ij-S_ij u_ji)^2].

Suppose B<=D/2. Then uniformly over ALL S and T,

    |log Z_S(t)-log Z_T(t)|
       <= E_rep :=4tV^2 B+D h(B/D)+log2.                  (1)

The result holds for every finite t>=0; its constant is explicit. If B>D/2, the proof still gives the same statement with D log2 in place of D h(B/D), but no small-error conclusion is implied.

## 2. Condition on all magnitudes

Generate a row by first uniformly assigning its magnitude multiset to ports, then independently assigning fair signs at its nonzero ports, conditioned on obtaining the prescribed signed counts in each magnitude class. This gives exactly the uniform signed-multiset permutation law.

The prior probability of that sign-count event is

    p0=product_(i,a) [2^(-n_(i,a)) binom(n_(i,a),n_(i,a)^+)],

independent of the magnitude arrangements and of S. If signs are left unconditioned, the edgewise partition is instead

    Z_fold(t)=E_magnitudes product_(i<j)
          exp[-t(a_ij^2+a_ji^2)] cosh(2t a_ij a_ji),

which is exactly independent of S.

For each fixed magnitude arrangement, tilt the independent fair endpoint signs by the edge kernel. The resulting measure factors over UNDIRECTED edges. Each endpoint sign is still fair, because simultaneous reversal of both signs of an edge preserves its weight. Signs on distinct edges are independent, so within any fixed row/magnitude class its positive count is binomial (n_(i,a),1/2), even though opposite endpoints of one edge are correlated.

The expected number of port-sign changes needed to repair every row to its target counts is at most

    (1/2)sum_(i,a) [sqrt(n_(i,a))
                            +|n_(i,a)^+-n_(i,a)^-|]
       =(R+A_asym)/2.

Hence with probability at least one half, at most B sign changes suffice. Fix a deterministic rowwise repair rule.

Each port-sign flip changes the log edge weight by at most 4tV^2. A repaired sign array has at most sum_(j<=B) binom(D,j) preimages within B flips. For B<=D/2 this is at most exp[D h(B/D)]. Comparing the Gibbs masses of an input and its repaired output therefore gives, uniformly in the magnitude arrangement and S,

    Pr_(edgewise tilted signs){all target counts}
       >=exp(-E_rep).                                    (2)

The reverse bound is one. Consequently

    exp(-E_rep) Z_fold(t)/p0 <=Z_S(t)<=Z_fold(t)/p0.        (3)

Equation (1) follows immediately by comparing the common interval (3). This is an ACTUAL partition comparison, not just equality of a variational upper bound.

## 3. Growing types and almost-symmetric profiles

Let a_m=A_asym/D. Then B/D<=O(sqrt(d/m)+a_m+m^-2). Thus the leading m^2 pressure is independent of S whenever

    d=o(m), a_m=o(1),
    t V^2 [sqrt(d/m)+a_m] ->0.                            (4)

This includes exactly symmetric signed profiles (a_m=0), arbitrary rare atom counts, heterogeneous profiles, and the loss of one diagonal coordinate per row (which contributes only A_asym<=m). It also gives a quantitative almost-symmetric extension. Constants do not deteriorate when individual positive frequencies tend to zero.

For completeness, the prior conditioning term is also subleading under the first two assumptions of (4). The elementary binomial type bound gives

    -log p0 <=sum_(i,a)log(n_(i,a)+1)+ (log2) A_asym
             <=dm log(1+(m-1)/d)+(log2) A_asym.           (5)

The inequality for the relative entropy uses h((1+s)/2)>=(1-|s|)log2. Consequently log Z_S and log Z_fold themselves differ by o(m^2) under (4), not merely the comparison between two seeds.

For bounded Hadamard spectra whose magnitudes lie on a mesh of order m^(-1/2), d=O(V sqrt(m)). At fixed V,t and A_asym=O(m^(7/4)) (in particular for exactly symmetric profiles or one-coordinate deletions), the error from (1) is O(m^(7/4)log m), with no fixed-alphabet assumption. General vanishing a_m still gives an o(m^2) error without this particular rate. This statement concerns only those row-profile sectors; it does not say every Boolean row query has such a profile.

## 4. Uniform tail truncation and exact scope

There is a simple deterministic truncation estimate. If the squared energy of all ports is at most C m^2 and the energy in ports with magnitude greater than V is at most epsilon m^2, replacing those ports by zero changes the logarithm of every raw-kernel product, and therefore its log partition, by at most

    4t sqrt(C epsilon) m^2.                              (6)

To see this, write the edge difference map as L_S u=(u_ij-S_ij u_ji)_(i<j). Its norm is at most sqrt2 for every S. Compare ||L_Su||^2 and ||L_Su_truncated||^2 by Cauchy--Schwarz. This yields (6) uniformly in all permutations.

Thus at fixed temperature, uniformly integrable squared tails can be truncated and then treated by the growing-type or finite quantization argument. Genuine macroscopic signed asymmetry and nonuniformly integrable tails remain outside this conclusion. In particular the actual mixed-orbit weave produces row profiles which need not be symmetric; this theorem does not erase its entire seed dependence or its law-level K5 statistic.

For homogeneous symmetric types the older hard compiler has much stronger low-temperature scope, including temperatures approaching m^2. Equation (4) is not a replacement for that theorem. Its additional role is uniform control of growing/rare/heterogeneous approximately symmetric types at the finite temperatures used in the strict cap certificates.
