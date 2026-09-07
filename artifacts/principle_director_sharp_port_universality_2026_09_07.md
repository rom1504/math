# Sharp logarithmic scale for complete port-defect universality

2026-09-07. **PROVED; aggregate matching independently reconstructed by the
construction reviewer. Full counterexample and dependency audit requested.**
This strengthens the complete seed-pressure theorem and supplies a matching
scale counterexample to weakening its hypothesis to ordinary quasirandomness.
It is NOT an impossibility theorem for globally rewritten exact minimizers.

For a hollow signing S of order m and real directed-port multisets u_i of
length m-1, put

    D_S(u,pi)=sum_(i<j)(u_(i,pi_i(j))-S_ij u_(j,pi_j(i)))^2,
    Z_S(u;t)=E_pi exp(-t D_S(u,pi)),
    d_S(u)=min_pi D_S(u,pi),

where the row permutations are independent and uniform. Assume EACH row
has squared norm at most Cm, with fixed C.

## 1. The theorem and its sharpness

If Q(S_m),Q(T_m)=o(m^2/log m), then, uniformly in these row profiles,

    |log Z_S(u;t)-log Z_T(u;t)|=o(m^2)       (every fixed t>0),
    |d_S(u)-d_T(u)|=o(m^2).

The associated uniform Laplace principle holds for every t_m tending to
infinity, as in the complete seed-pressure artifact. For each fixed error
epsilon, near-optimal placements have probability at least exp(-A_epsilon m^2).
Uniformity in seeds means Q(S),Q(T)<=q_m for a fixed envelope
q_m log m/m^2 tending to zero; no unquantified uniformity over little-oh
sequences is intended.

The scale cannot be enlarged to Q=o(m^2): there are actual S_m,T_m and
profiles of EXACT row energy m with

    Q(T_m)<=m^(3/2)/2+O(m),
    Q(S_m)=Theta(m^2/log m),
    d_T(u)=0,               d_S(u)>=m^2/5,
    log Z_T(u;t)=o(m^2),    log Z_S(u;t)<=-t m^2/5.

Thus even the leading ground-defect value and every fixed-temperature
pressure distinguish two cut-quasirandom sign sequences. This does NOT
contradict universality on the actual bounded-cap minimizing class.

## 2. An aggregate negative-energy lemma

Use the descending-bin greedy matching in
`principle_synthesis_2026_09_07_complete_seed_pressure_universality.md`.
Heavy magnitudes lie in bins [a_b,(1+delta)a_b), a_b=V(1+delta)^b;
the number B of nonempty bins is O_(C,delta,V)(log m).
The used degree at bin b is at most D_b=Cm/a_b^2.

Let e_ib be the squared mass of the unmatched stubs at vertex i in bin b.
The total over bins at any vertex is at most Cm. Independently assign
vertex i to bin b with probability e_ib/(Cm), or to no bin with the
remaining probability. A selected vertex has an actual unmatched stub;
choose one of its signs epsilon_i. The resulting sets U_b are DISJOINT.

For ANY realization of the assignment, maximality gives

    E_b=sum_(i<j in U_b) S_ij epsilon_i epsilon_j
        <= |U_b| D_b-binom(|U_b|,2).

Call a bin good if n_b=|U_b|>=4D_b+2. Then E_b<=-n_b^2/4.
Apply independent uniform global sign flips to the words on the good U_b,
and independent fair spins on every remaining coordinate. Cross-bin and
outside contributions vanish in expectation. Since every full energy
is at least -Q(S),

    sum_(good b) n_b^2 <=4Q(S).

Consequently, pointwise for EVERY assignment,

    sum_b n_b <= sum_b(4D_b+2)+2sqrt(BQ(S)).

Taking expectation on the left gives E_bad/(Cm), with
E_bad=sum_ib e_ib. Therefore

    E_bad/m^2
      <= 4C^2/[V^2(1-(1+delta)^(-2))]
          +2CB/m+2C sqrt(BQ(S))/m.                 (1)

This improves the previous bound that paid sqrt(Q) separately at every
amplitude level. The gain comes from simultaneously exposing the negative
energy of DISJOINT residual sets, not an assumed cancellation in a parent.

## 3. Consequence for pressure and minimum defect

Under Q=o(m^2/log m), the last two terms of (1) tend to zero. The first
can be made arbitrarily small by choosing V large after delta is fixed.
Every remaining step of the complete compiler is unchanged: compatible
heavy placements, category entropy, tagged-zero swaps, and a profile-chosen
low-energy annulus reduce uniformly to a fixed bounded profile.

The bounded-profile theorem requires only Q=o(m^2), not the stronger
O(m^1.5) display in its first application. Indeed its fixed finite-type
comparison error is O_(types,t)[beta(S)+beta(T)+m^1.5 log m], and
beta<=4Q. Quantization and rare-atom pruning are uniform at fixed accuracy.
Thus it applies under the present hypothesis. Sending m to infinity
after fixing the finite grid, and then its accuracy to zero, proves the
fixed-temperature statement. The deterministic heavy construction and the
finite-grid entropy bound prove the minimum-defect and uniform Laplace
statements with exactly the same order of limits. No growing grid is
inserted silently into a fixed-type theorem.

## 4. An explicit logarithmic-scale counterexample

Let B=4^b tend to infinity, s=4^(B+2), and m=Bs. Split the m vertices
into B groups of size s. Write H_4=J_4-2I_4, and H_(4^k)=H_4 tensor ...
tensor H_4 with k factors. These symmetric full sign matrices satisfy

    H_d^2=dI,   H_d 1=sqrt(d)1,   diagonal(H_(4^k))=(-1)^k.

Take the appropriate global polarity of H_B tensor H_s so that each
diagonal group block is -H_s, then set the full diagonal to zero.
Call the result T. Its operator norm is at most sqrt(m)+1, so
Q(T)<=m^(3/2)/2+m/2.

Because B+2 is even, H_s has diagonal +1. Thus the positive-edge graph
of the hollow -H_s block is regular of degree

    d_+=(s-sqrt(s))/2,

which is an even integer. Every finite even-regular graph decomposes into
2-factors: orient an Euler tour in each component, giving equal in/out
degrees; split vertices into left/right copies and repeatedly use a
perfect matching in the regular bipartite graph. Each matching is a
directed cycle cover. Opposite directions of the same edge are absent,
so each cover yields an ordinary 2-factor. Hall's theorem here follows
by counting incident edges. This supplies all even-degree factors up to d_+.

For group g=1,...,B choose

    h_g=2*4^(B-g),      a_g=sqrt(m/h_g).

Each row in that group has exactly h_g copies of a_g and all other
ports zero. Its energy is h_g a_g^2=m. Moreover h_g is even,
h_1=s/32<d_+, and successive amplitudes have ratio two.

Choose an h_g-factor of positive T edges within group g and place the
nonzero ports on its incident edges. Every such edge has equal positive
amplitudes at its two ends; all remaining edges have zeros at both ends.
Thus this ACTUAL simultaneous port placement has D_T=0.

Now define S by replacing every within-group edge by -1 and leaving all
between-group edges as in T. Its perturbation from T is block diagonal
with each block H_s-J_s. Hence

    Q(S)<=Q(T)+m(s+sqrt(s))/2.

At the all-positive word the perturbation energy is exactly
B(s^(3/2)-s^2)/2, while |H_T(1)|<=Q(T). Since s=m/B and
B=Theta(log m), these upper/lower estimates give

    Q(S)=(1/2+o(1))m^2/B=Theta(m^2/log m)=o(m^2).       (2)

## 5. Every S placement pays a macroscopic defect

All nonzero amplitudes are positive. If an edge has nonzero endpoints
of equal magnitude, its endpoints belong to the same group and S_ij=-1;
its defect is (a+b)^2>=a^2+b^2. If magnitudes differ, their ratio is
at least two, and for either edge sign

    (a-S_ij b)^2 >= (a-b)^2 >=(a^2+b^2)/5.

With only one nonzero endpoint the same bound is immediate. Summing
over unordered edges counts every directed-port square exactly once.
Therefore D_S>=m^2/5 for EVERY independent choice of row permutations.

Finally the zero-defect T placement has probability at least

    product_(g=1)^B binom(m-1,h_g)^(-s).

Its negative logarithm is at most
s sum_g h_g log[e(m-1)/h_g]=O(s^2 log B)=o(m^2), because h_g decreases
geometrically. Hence -o(m^2)<=log Z_T<=0 for fixed t, whereas
log Z_S<=-tm^2/5. The same example also has a ground-defect gap.

## 6. What principle this proves

Weak cut quasirandomness preserves bounded observations but not arbitrary
energy-bounded port observations: an ordered logarithmic family of coherent
amplitude scales can expose mesoscopic sign communities. Quantitatively,
the aggregate cut budget controls those scales through sqrt(BQ)/m.
The complete universality threshold Q=o(m^2/log m) and construction (2)
match in scale. This is a statement about a particular permutation response
algebra, not a universal information barrier, and it leaves correlated
placement and the simultaneous parent-spin optimization untouched.
