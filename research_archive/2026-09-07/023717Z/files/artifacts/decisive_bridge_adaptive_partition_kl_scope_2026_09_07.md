# Adaptive partitions: exact freedom and the cut-switch KL scale

This examines actual globally optimized partitions and the proposed use of cut-switch KL to maintain a favorable block comparison. It does not replace the missing asymptotic kernel theorem.

## 1. Full reoptimization leaves no extra partition freedom

Let w(u,P) be a magnitude profile determined by a partition P into fixed block sizes: magnitudes depend only on u and the blocks containing the endpoints. Write

    F_n(u,P)=min_{A_e∈{±1}} log Z(A,w(u,P)).

For any two partitions P,P' with the same block sizes, a vertex permutation π carries P to P'. Relabeling both spins and A proves

    F_n(u,P)=F_n(u,P')

exactly at every u. This holds also for the quenched global-orientation field, since that scalar field is unaffected by vertex permutations.

Consequently, choosing P(u) adaptively and fully reoptimizing the signing follows exactly the same scalar envelope as a fixed partition. At every differentiability point the full derivative of any active optimizer agrees with that envelope derivative. A change of partition costs zero when accompanied by full relabeling/reoptimization, but cannot thereby improve the integrated scalar defect.

At a homogeneous starting point, tied minimizers can indeed be placed to select a favorable one-sided derivative. This does not imply an additional choice at positive u: choosing a favorable partition relative to a **frozen** A generally moves A away from the globally minimizing class. A separate bound is needed on the resulting objective increase. The exact positive derivative in §6 of decisive_bridge_actual_cavity_block_defect_2026_09_07.md demonstrates that all active optimizers can have the same unfavorable derivative at an interior point.

This does not preclude a more general interpolation in which block sizes, magnitude distributions, or auxiliary state variables evolve. It rules out obtaining a new bound merely by relabeling the fixed two-block profile while minimizing over all signings as before.

## 2. Exact cut-switch KL identity

For the homogeneous two-sided Gibbs law

    p_A(σ,x)=Z_A^−1 exp(λσH_A(x)),   λ=β/√n,

let T_S flip the spins in S and set p_A^S(σ,x)=p_A(σ,T_Sx). Its normalizer is exactly the same. Therefore

    D(p_A || p_A^S)=2λ Σ_{e∈∂S} A_e〈σxi xj〉.          (1)

This proves every cut inequality used in the cavity degree bound and strengthens its interpretation. It requires no stationarity or optimality.

For even n and a uniformly chosen balanced cut S, the probability that a fixed edge crosses is n/[2(n−1)]. Thus

    E_S D(p_A || p_A^S)
       = [n/(n−1)] λ ∂_λ log Z_A.                    (2)

By convexity of log Z_A in λ,

    λ ∂_λ log Z_A ≥ log Z_A−log Z_A(0)
                  ≥ λ Q(A)−(n+1)log2.              (3)

These identities hold in particular at every actual globally optimized signing.

## 3. The averaged KL budget is necessarily extensive

A short elementary lower bound suffices. Split the n vertices into two sets of size m=n/2. For random signs on one side and optimized signs on the other, the expected cross bilinear sum is m E|S_m|, where S_m is the sum of m independent fair signs. Hölder and its second and fourth moments give

    E|S_m| ≥ (E S_m²)^(3/2)/(E S_m⁴)^(1/2)
            =m^(3/2)/√(3m²−2m) ≥√(m/3).

The cross bilinear sum is half the difference of two quadratic energies obtained by flipping one side, so Q(A) dominates it. Hence every signing satisfies

    Q(A) ≥ n^(3/2)/(2√6).

Combining with (2)–(3), already at β=4,

    E_balanced S D(p_A || p_A^S)
       ≥ [2/√6−log2+o(1)]n
       > (.1233+o(1))n.                              (4)

The established stronger asymptotic cap lower bound improves this coefficient, but is unnecessary for the conclusion. Thus even actual global minimizers do not have a generally sublinear *averaged cut-switch KL budget*.

This does not prove that every cut is costly, nor that pressure changes equal these KL divergences: a gauge switch preserves pressure exactly despite a positive KL divergence. It shows why a direct sum/average-of-KL budget cannot by itself justify a sublinear partition-realignment error. A successful argument must select special cuts or use cancellations beyond such a budget, and must retain its favorable cavity comparison after that selection.
