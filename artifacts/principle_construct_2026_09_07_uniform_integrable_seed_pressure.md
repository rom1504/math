# Fixed-temperature seed-pressure universality beyond finite types

2026-09-07. **Proved extension; audit requested.** The finite-type theorem is due to the earlier `principle_invent_2026_09_07_finite_type_seed_universality.md`. I independently reconstructed its canonical dual, positivity bound, repair likelihood identity, and smooth separated-kernel comparison. This note removes the finite-alphabet and positive-mass restrictions at fixed temperature by deterministic approximation. It does not compare actual parent maxima.

## 1. Uniform bounded-profile theorem

At each of m vertices prescribe an arbitrary real multiset of length m-1, independently uniformly permuted among its directed ports. The rows may be heterogeneous and signed-asymmetric, with arbitrarily many distinct coordinates. Suppose every coordinate has magnitude at most a fixed V. For a hollow macro signing S, let

    Z_S(t)=E product_(i<j) exp[-t(u_ij-S_ij u_ji)^2].

Fix t,V,C<infinity. Uniformly over ALL such profile arrays and all full-sign macro seeds S,T with Q(S),Q(T)<=C m^(3/2),

    |log Z_S(t)-log Z_T(t)|=o(m^2).                    (1)

More generally (1) holds along any seed sequences with beta(S)+beta(T)=o(m^2). No symmetry of the signed profiles is required.

## 2. Deterministic quantization and rare-atom pruning

If a deterministic map is applied to every value of a uniformly permuted multiset, the image is EXACTLY a uniformly permuted copy of the image multiset. Indeed every image arrangement has the same product of multinomial refinement counts. The map may differ from row to row. This observation avoids any entropy fee when values are merged.

Choose a symmetric grid in [-V,V] of spacing at most delta, with L<=2ceil(V/delta)+1 values. Move every original value to a nearest grid point. For an edge, changing one endpoint by at most delta changes its log kernel by at most 4tV delta. Summing over all D=m(m-1) directed ports gives

    |log Z_S - log Z_S^(grid)|<=4tV delta D,           (2)

uniformly in S and in the original profiles.

Next choose alpha<1/L. In each row, merge every positive-frequency grid atom of frequency below alpha into a most-common atom. Such an atom has frequency at least 1/L and is not removed. At most L alpha of that row's ports change. Every retained positive frequency is now at least alpha, while all values stay in [-V,V]. A changed port changes the log kernel by at most 8tV^2, so

    |log Z_S^(grid)-log Z_S^(pruned)|
       <=8tV^2 L alpha D.                             (3)

The exact pushforward property above applies to both transformations. Thus (2)--(3) are comparisons between normalized permutation EXPECTATIONS, not just between individual configurations or unnormalized word counts.

The pruned arrays satisfy every hypothesis of the earlier finite-type theorem, with fixed alphabet L, fixed positive mass alpha, and kernel condition number controlled by t,V. That theorem gives

    |log Z_S^(pruned)-log Z_T^(pruned)|
       <=C_(L,alpha,t,V)[beta(S)+beta(T)+m^(3/2)log m].  (4)

Divide by m^2, let m tend to infinity, and use the assumed seed cut-norm condition. Equations (2)--(4) give a limsup bounded by

    8tV delta+16tV^2 L alpha.

First choose delta arbitrarily small, then alpha arbitrarily small at that fixed grid. The constants in (4) can depend on both: those parameters are fixed before m grows. This proves (1), uniformly over all original bounded profiles. There is no unjustified growing-dimension limit.

## 3. Uniformly integrable squared tails

Assume more generally that the total port energy is at most C_row m^2 and that the profile family has uniformly integrable squared tails:

    sup_arrays m^(-2) sum_ports u_port^2 1_(|u_port|>V)
       <=omega(V),    omega(V)->0.

Truncating all values of magnitude greater than V to zero changes each log partition by at most

    4t sqrt(C_row omega(V)) m^2.                       (5)

For completeness, the edge difference operator L_S sends the two ports of an edge to u_ij-S_ij u_ji and has norm at most sqrt2. Comparing ||L_Su||^2 with its truncated version gives (5), uniformly in every row permutation and every S. The pushforward of the row permutation law remains exact after truncation.

Apply (1) to the bounded truncated arrays and then let V grow. Equation (5) proves the same o(m^2) seed-pressure comparison for this entire uniformly integrable family, including genuinely signed-asymmetric profiles and arbitrary alphabets.

## 4. Corrected escape boundary

At FIXED temperature and with uniformly integrable squared port energy, growing alphabets, rare atoms, and signed asymmetry alone do not preserve a leading distinction between different low-cap seeds. The approximation above closes those apparent escapes from the finite-type theorem. This does not assert equality with the folded-sign kernel: for asymmetric profiles the common seed-free value can differ from the independently sign-averaged value, as the finite-type theorem already explains.

The remaining scope limitations are substantive: nonuniformly integrable energy tails, growing/singular temperature or kernels, or joint row dependence outside the independent permutation-type model. The approximation constants are not uniform in those regimes.

In the actual restricted Hadamard weave, the tail issue is unavoidable rather than merely hypothetical. For any selector T of k rows, the spin word equal to a selected column has normalized spectrum with one coordinate sqrt(k); its total squared norm is m. At retention k/m tending to p>0, that one coordinate carries asymptotic fraction p of the row energy. Hence the family of ALL Boolean row queries is not uniformly square-integrable. The theorem cannot erase the complete fixed-seed compiler or its law-level cycle statistic. It identifies the coherent high-amplitude sectors as the remaining place where a leading seed-value effect could occur at finite temperature.
