# Aggregate unmatched-heavy-energy audit

2026-09-07. **PASS**, independently reconstructed from the director's
proposed strengthening of the descending matching lemma.

Let e_ib be the unmatched squared energy at vertex i in bin b. Since bins
use disjoint original copies, sum_b e_ib<=Cm. Independently assign vertex i
to bin b with probability e_ib/(Cm), or leave it unassigned. At a vertex
assigned to b, select any sign of an available unmatched stub in b.

For every realized assignment the groups U_b are disjoint. With
`D_b=Cm/a_b^2`, every unused edge inside U_b has switched sign minus one,
while the number of used edges is at most |U_b|D_b/2. Therefore its
induced signed energy is at most

    |U_b|D_b - binom(|U_b|,2).

If |U_b|>=4D_b+2, this is at most -|U_b|^2/4. Flip each such group by an
independent common sign; assign independent unbiased spins to every other
vertex. All intergroup and outside contributions vanish in expectation,
while the negative intragroup energies remain. Consequently

    sum_good |U_b|^2 <=4Q(S).

Cauchy--Schwarz gives sum_good |U_b|<=2sqrt(BQ(S)); the other groups have
total size at most sum_b(4D_b+2). Thus every assignment satisfies

    sum_b |U_b| <= 4sum_b D_b+2B+2sqrt(BQ(S)).

Taking expectation of the left side gives exactly E_bad/(Cm), without
any independence assumption about different bins. Summing the geometric
sequence of D_b yields

    E_bad/m^2 <=
      4C^2/[V^2(1-(1+delta)^(-2))]
      +2CB/m+2C sqrt(BQ(S))/m.

For fixed relative bin width, B=O(log m). Therefore the complete
heavy-tail pressure and minimum-defect reductions extend from
Q(S)=O(m^(3/2)) to `Q(S)=o(m^2/log m)`, uniformly when this last ratio
is uniformly vanishing over the allowed seed class. Their bounded-profile
dependency requires only beta(S)=o(m^2), which the stronger hypothesis
implies through polarization. This audit does not assert sharpness of
the logarithmic threshold; a separate actual counterexample is needed.

## Full sharpness theorem and explicit counterexample

The subsequently completed
`principle_director_sharp_port_universality_2026_09_07.md` was read in
full and **PASS**es, including the matching-scale counterexample.

Set B=4^b, s=4^(B+2), m=Bs. The regular symmetric Hadamard tensor of
order m has constant block diagonal polarity; a global sign makes every
diagonal group block -H_s. Since B+2 is even, H_s has diagonal +1.
Hollowing the full matrix therefore changes its operator norm by at most
one, giving Q(T)<=m^(3/2)/2+m/2.

The positive graph in each -H_s block has even degree
`(s-sqrt(s))/2`. Euler orientation produces equal in- and out-degrees;
the associated regular bipartite graph decomposes into perfect matchings.
No edge is oriented in both directions, so their directed cycle covers
are ordinary 2-factors. Thus every requested even degree
`h_g=2*4^(B-g)`, bounded by h_1=s/32, is available as a positive factor.
Placing equal positive amplitudes sqrt(m/h_g) along those factors gives
exact row energy m and zero total T defect.

Changing all intragroup edges to minus one changes each block by exactly
H_s-J_s, including zero diagonal. Its all-positive energy is
`(s^(3/2)-s^2)/2`. The original T cap is lower order relative to m^2/B,
so both the displayed perturbation upper bound and this test word give
`Q(S)=(1/2+o(1))m^2/B`. Here B is proportional to log m.

Amplitudes in different groups differ by at least a factor two, while
equal nonzero amplitudes belong to the same all-negative group. Every
edge therefore pays at least one-fifth of its endpoint squared energy.
Summation counts each directed-port square once, proving D_S>=m^2/5
for every placement. A prescribed positive-factor placement for T has
probability product_g binom(m-1,h_g)^(-s); its negative logarithm is
O(s^2 log B)=o(m^2). This proves the fixed-temperature pressure gap as
well as the minimum-defect gap.

The bounded-profile dependency indeed requires only Q=o(m^2): at fixed
type and temperature its seed error is a constant times beta(S)+beta(T),
plus O(m^(3/2)log m). Quantization/pruning uses fixed parameters before
the order limit. Therefore no stronger hidden cap premise is needed for
the sharp logarithmic theorem. T is an actual bounded-cap signing, not
claimed to be an exact minimizer; S lies at the larger threshold scale.
