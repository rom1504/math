# Actual optimized cavity kernels: exact counterexamples and the needed scale

Fix λ=β/√n and let A globally minimize the two-sided partition Z_A=Σ_{σ,x}exp(λσH_A(x)). Delete edge ij and write r^0_ij=〈σxi xj〉_0, Kij=|r^0_ij|, Kii=0, t=tanh λ. This note concerns actual global partition minimizers, not arbitrary correlation matrices.

## 1. The precise block condition and what it would buy

For a bipartition of proportions r,1−r let z_i=√((1−r)/r) on the first block and z_i=−√(r/(1−r)) on the second. Then Σz_i=0 and ||z||²=n. Define

    D(K)=Σ_{i<j} Kij[1−1{i,j in first}/r−1{i,j in second}/(1−r)].

The exact algebraic identity is D(K)=−z^T Kz/2. Consequently the sufficient asymptotic condition is

    sup_relevant_blocks D(K)_+ = o(n^(3/2)).                    (P)

A stronger sufficient condition is λmin(PKP)≥−o(√n), where P projects onto the zero-sum subspace. Exact conditional PSD is unnecessarily strong. A fixed O(1) diagonal correction is harmless, since it contributes only O(n) to these quadratic forms.

The independent agent's exact edge-optimal identity gives an insertion cost log cosh λ+log(1−tKij). Its first-order adaptive term is −tKij; its total Taylor remainder over O(n²) edges is O_β(√n), using ΣKij=O_β(n^(3/2)). Thus, **if an actual interpolation has already been justified with this edge-cost difference**, (P) makes its unfavorable first-order block term at most t·o(n^(3/2))=o(n), and the log remainder is also sublinear. This proves the implication from the kernel condition to a sublinear interpolation defect, not the interpolation itself. In particular an edge-insertion proof must handle reoptimization of all old signs; evaluating a new-edge cost at an old global optimizer gives only a one-sided inequality.

For block ratios bounded away from 0 and 1, bounded scalar interpolation coefficients do not change this conclusion. Endpoint layers or unbounded coefficients require their own integrable control.

## 2. Stronger exact degree and cut bounds

For any fixed A and any vertex subset S, scale all couplings crossing S by a real scalar a. The log partition is convex in a, and its values at a=1 and a=−1 agree by switching every spin in S. Therefore

    Σ_{e∈∂S} A_e〈σxi xj〉 ≥0.

At the global optimizer, edge optimality says A_e r^0_e=−K_e, so

    A_e〈σxi xj〉 = t−(1−t²)K_e/(1−tK_e).

Combining gives the exact family of cut inequalities

    Σ_{e∈∂S} K_e/(1−tK_e) ≤ |∂S|t/(1−t²).

In particular max_i Σ_j Kij≤(n−1)t/(1−t²)=O_β(√n). This improves the total ℓ1 control to a uniform degree bound. Gershgorin then gives λmin(K)≥−O_β(√n), precisely the critical scale, not the o(√n) bound required by (P). The cut inequalities by themselves are not the desired within-versus-cross comparison.

## 3. Exact n=4 failure of uncorrected conditional PSD

For n=4, after switching to make the first row positive, the six nonconstant choices of the remaining three signs are exactly the global partition minimizers for every λ>0. Their absolute energy histogram on the eight spin configurations with x_0=1 is

    |H|=0,2,4 with counts 2,4,2.

The two excluded choices have counts 4,3,1 at energies 0,2,6. Their cosh-sum excess is

    cosh(6λ)−2cosh(4λ)−cosh(2λ)+2
       =4[cosh(2λ)−1]²[cosh(2λ)+1]>0.

For the minimizing upper-triangle signs (1,1,1,−1,−1,1), K is a perfect-matching adjacency matrix multiplied by

    a=2t³/(1+t⁴)>0.

It has a zero-sum eigenvector of eigenvalue −a. Thus exact conditional PSD already fails at genuine global minimizers. Here I+K is PSD, so this does not challenge a bounded diagonal correction or (P).

## 4. Exact n=7 failure even after adding the identity

At n=7 and λ=log 2, all 2^15=32768 switching classes were enumerated with exact integer partition weights. There are exactly 840 global minima. In lexicographic upper-triangle edge order, one witness is

    A=(1,1,1,1,1,1,−1,−1,−1,1,1,−1,1,−1,−1,1,1,1,1,−1,−1).

For the zero-sum block vector z=(4,−3,−3,4,−3,4,−3), the cavity correlations are rational and the exact result is

    z^T(I+K)z = −1140197007043868844/228780342883779149
              < −4.9838.

This is a counterexample to diagonal-repaired conditional PSD at an actual homogeneous global optimizer. It is not a counterexample to the asymptotic o(√n) spectral relaxation.

The complete, independently rerunnable integer enumeration and rational quadratic-form check is

    computations/decisive_bridge_actual_cavity_exact_n7_2026_09_07.py.

The partition comparison uses Σ_x [2^(21+|H|)+2^(21−|H|)] over x_0=1; all integers fit in int64. Edge-deleted energies are even and bounded by 20, allowing exact rational cavity numerators and denominators. Switching exhausts all full sign matrices and preserves the cavity absolute kernel up to relabeling, so the reduced enumeration is exhaustive.

## 5. Finite exploratory interpolation tests

The script decisive_bridge_actual_cavity_psd_checks_2026_09_07.py exhaustively minimizes the numerical partition for n=3,...,7, β∈{.2,1,2,4,8}, and weights

    w_internal,i=√[1+u(1/r_i−1)],  w_cross=√(1−u),
    u∈{0,.25,.5,.75,1}.

For homogeneous weights it checks every minimizing switching class and every two-block split. For inhomogeneous weights it checks the designated split, as well as the whole zero-sum spectral subspace. These are floating exploratory tests, unlike the exact witness above. Results are in decisive_bridge_actual_cavity_psd_results_2026_09_07.jsonl.

At homogeneous n=7, β=4, the minimum zero-sum eigenvalue is approximately −2.347 and the worst normalized block defect D/n^(3/2) is approximately .4124. Small orders therefore provide no numerical evidence for a simple uniformly bounded diagonal repair, but they cannot establish or refute (P). A valid next theorem must use asymptotic global optimality beyond the exact degree/cut inequalities above.

## 6. Exact positive block derivative after arbitrary orientation synchronization

The independent agent regularizes with a quenched branch field hgσ, g standard Gaussian, minimizing E_g log Z_A(g) with A chosen before g. The following counterexample is valid for **every h≥0**, including arbitrarily strong synchronization.

Take n=6, two blocks of size3, and the variance-redistribution path

    λ_internal(u)=(β/√6)√(1+u),
    λ_cross(u)=(β/√6)√(1−u).

Fix β=(√15/2)log2 and u=3/5. Then λ_internal=log2 and λ_cross=(log2)/2. Exhaustive integer enumeration of all 1024 switching classes proves that exactly twelve matrices minimize the product Z_+(A)Z_−(A), and every one of those twelve is balanced: Z_+=Z_−. One witness has upper-triangle signs

    (1,1,1,1,1,−1,−1,1,1,1,−1,1,1,−1,−1).

This also certifies the global minimizers of the quenched objective for every h. Indeed write b_A=(log Z_++log Z_−)/2 and d_A=(log Z_+−log Z_−)/2. Then

    E log Z_A(g)=b_A+E log[2cosh(d_A+hg)].

The second term is even convex in d_A and minimized at zero. Every matrix has b_A≥b_min, while the twelve product minimizers attain both b_min and d=0. Thus they are exactly the global quenched minimizers for every h≥0.

At each of these twelve minimizers the exact full-correlation directional sum is

    Σ_internal A_e E〈σxi xj〉 −2Σ_cross A_e E〈σxi xj〉=58/119.

It is independent of h because d=0 and E tanh(hg)=0. Along the path, λ'_internal=5log2/16 and λ'_cross=−2λ'_internal. All active derivatives agree, so the minimized envelope is differentiable, and

    d/du [6^−1 min_A E log Z_A(g)] = (145/5712)log2 >0.

The exact checker is decisive_bridge_quenched_block_exact_n6_2026_09_07.py. Energies in units (log2)/2 are odd integers from −21 to21; multiplying the branch partitions by 2^(21/2) makes them integers. The minimum scaled branch product is 42284519424. All comparisons and the final derivative use integers or rational arithmetic.

Therefore neither actual global optimality nor arbitrarily strong orientation synchronization implies exact monotonicity along the block variance path. This addresses the complete derivative, not merely a truncated first-order kernel. It still does not exclude a sublinear asymptotic defect: the example is fixed order, and no globally optimal dense sign-matrix enlargement of it has been proved.
