# Exact rare-packet temperature tradeoff for one-hole selectors

2026-09-07. **Analytic derivation; independently matched by invention
agent; exhaustive integer regression PASS.** This
refines the stratified marked-selector theorem at k=L-1, L>=4. It is a
statement about that proved finite-block envelope, not an actual cap
lower bound and not a convergence obstruction.

Let H be any dephased Hadamard of order L. For a nonempty proper packet
S, write s=|S| and C_j=sum_(l in S) H_jl. For 1<=c<=L define

    rho(c)=min_S [sum_j min(C_j^2/c,1)-1]/s.

Then EXACTLY

    rho(c)=min(1/c, 2/L, L/c-1).                         (1)

If DC is unclipped, every coefficient is unclipped and the expression
is L/c-1/s >= L/c-1. If DC is clipped, discard its contribution one.
The non-DC squares sum to s(L-s), and their maximum is at most
min(s,L-s)^2. Their clipped sum divided by s is at least

    (L-s)/max(c,min(s,L-s)^2) >= min(1/c,2/L).

For the second inequality, (L-s)/c>=1/c; also
(L-s)/min(s,L-s)^2>=2/L, by splitting at s=L/2.
This proves the lower bound in (1).

Three exact packet types attain its three branches. A complement of
one coordinate has s=L-1 and all non-DC squares1, giving1/c whenever
this is the minimum branch. The positive support of a non-DC Hadamard
row has s=L/2 and exactly two nonzero coefficients, both of magnitude
L/2. Since c<=L<=L^2/4, it gives2/L. A singleton has all squares1 and
gives L/c-1. Thus no support classification is needed.

Equivalently R(c)=c rho(c) is the piecewise linear function

    2c/L       for 1<=c<=L/2,
    1          for L/2<=c<=L-1,
    L-c        for L-1<=c<=L.                           (2)

## Exact leading finite-block exponent

Put p=(L-1)/L, a=1-2r, and t=[k/(4c)] log(1/r). The uniform rare-source
information theorem and the group-word entropy bound imply

    Gamma_t(a) <= -(k/L) rho(c) r log(1/r)
                                      +o(r log(1/r)).

This leading coefficient is attained in the finite-block supremum.
Choose a minimizing fixed packet S with probability theta=kr/s; otherwise
take the empty packet. Conditional on empty, choose T uniformly among
k-subsets; conditional on S choose T uniformly among its k-supersets.
This is an admissible mean-a law. Its selector-adjusted entropy equals
theta log(1/r)+O_L(r), and every output obeys the matching rare-source
asymptotic. Therefore the displayed inequality is an equality to the
stated leading order.

The resulting leading variance-response coefficient of this certificate is

    b(c)=[1-R(c)/L]/(2 sqrt(p)).                        (3)

The best coefficient sqrt(p)/2 persists over the plateau c in
[L/2,L-1]. Allowing a slightly larger coefficient expands the available
temperature interval in both directions.

## Bandwidth/coefficient tradeoff

Fix sqrt(p)/2 < b <1/2 and let delta=1-2b sqrt(p), so
1-sqrt(p)<delta<1/L. Formula (3) is at most b precisely, within the
relevant range, when

    L^2 delta/2 <= c <= L(1-delta).                     (4)

For densities r_i in [r^alpha,r], a single common t gives c_i values
whose ratio is alpha. Thus the same proof certifies the band whenever

    alpha < 2(1-delta)/(L delta).                       (5)

Use strict inequalities to leave a uniform margin for finite depth and
the rare-source remainders. Constants and the usual uniform additive
repair error are handled exactly as in the canonical band theorem.
At the best coefficient, the limiting bandwidth is2(L-1)/L. As b tends
to1/2 from below, it tends to

    2 sqrt(p)(1+sqrt(p)) <4.

This sharpness concerns the fixed-temperature finite-block certificate:
the explicit rare packet policies attain its exponent. It does not
exclude a different construction, different joint estimate, or an
actual cap improvement invisible to this envelope.

Exact regression:
`../computations/principle_construct_2026_09_07_packet_curve_check.py`
checks every nonempty proper packet at L4,8,12,16 and every half-integer
c from1 toL. All2,129,624 packet/temperature pairs satisfy the integer
lower inequality, and every tested temperature has an equality packet.
The results JSON records exact rational minima and equality counts.
