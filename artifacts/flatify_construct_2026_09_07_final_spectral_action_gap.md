# What the new fourth-moment theorem pays in actual constructions

2026-09-07, final constructive assessment. This records quantitative consequences and failures, not a new auxiliary state or a claimed recurrence.

## Verified inputs

For every actual hollow full signing A with Q(A)<=C n^(3/2), the new factorization proof gives

    Tr A^4 <= D n^3,    D=16 K_G^2 C^2.

If P_K is its spectral projector onto |lambda|>K sqrt(n), then

    rank(P_K)<=D n/K^4,
    ||P_K A||F^2<=D n^2/K^2,
    ||P_K A||*<=D n^(3/2)/K^3.                         (1)

The last inequality follows by bounding |lambda| by |lambda|^4/(K sqrt(n))^3 on the tail.

## Actual bridge consequence, with all error paid

Apply the new random-gauge projector construction to the two child tail spaces. Its uniform actual-sign bridge error is

    O(n^(3/2) (log(n)/K^4)^(1/6)+n^(41/40)),            (2)

with constants depending only on the common child cap bound. Thus a growing cutoff K>>(log n)^(1/4) gives an o(n^(3/2)) realization error, and K=n^a for fixed a>0 gives a power saving. Both designated tail spaces are removed from the realized continuous bridge target.

This is a real operation on full signs, but it does not bound the total parent energy: after removal, the two child operators still have norm at most K sqrt(n), and their actual quadratic energies are unchanged. In particular the scalar residual-operator upper bound remains obstructed by the Rademacher near-orthogonal width theorem for every sublinear rank. The bridge may have favorable joint correlation with those energies, but no such paid inequality was proved.

A covariance construction whose positive-definiteness requires dividing these residual child operators by K sqrt(n) loses its fixed-strength interaction when K grows. Formula (2) therefore does not automatically transfer the fixed-ellipticity Gaussian MGF theorem to the desired original-energy variance proxy.

## Spectral surgery via nuclear-budget rounding

Attempt instead to replace A by its spectrally truncated hollow target. The nuclear tail budget in (1), combined with generic target contraction and sign rounding, has error at most on the scale

    n sqrt(||P_K A||*) = O(n^(7/4)/K^(3/2)).            (3)

This becomes a power saving against n^(3/2) only if K grows faster than n^(1/6), up to the desired error exponent. For smaller K one can always cap variance by the total number of edges, but that only returns a leading n^(3/2) rounding cost. Moreover the cap of the truncated target itself has not been shown to be at most Q(A)+o(n^(3/2)). Equation (3) is consequently not a truncation theorem.

Hollowing the tail changes only a diagonal and must be treated explicitly in any actual application; it does not resolve either of the preceding scale or cap issues.

## Coordinate removal is already scoped by an archive counterexample

The Grothendieck-Pietsch regular-core theorem in `grothendieck_pietsch_spectral_regularization.md` deletes at most epsilon n vertices to give an operator bound O_C(epsilon^(-1) sqrt(n)). Random refill can pay O(sqrt(epsilon) n^(3/2)) in cap. This is useful approximate regularization, not a fixed-cutoff, power-saving completion theorem.

The full-sign clique tower in `spectral_peeling_counterexample.md` already shows that low cap alone does not permit o(n) coordinate deletion with a fixed O(sqrt(n)) residual operator constant. The new fourth-moment theorem does not contradict that construction: its large clique modes have small aggregate spectral mass, while every fixed threshold still encounters a positive vertex fraction in some clique class.

The tower can have arbitrarily small fixed cap excess by implanting it in a good signing. It does not establish the same obstruction with o(n^(3/2)) excess, or among exact minimizers. Conversely, no theorem charging these coordinate irregularities to optimality excess was proved here.

## Strongest remaining gap

We now have actual full-sign realizations of prescribed projected Hadamard bridges through every power-sublinear rank, a dimension-free elliptic Gaussian-sign MGF expansion, and bounded normalized fourth spectral moments for all actual low-cap inputs. What is still missing is a uniform joint inequality for the unchanged or globally modified child energies and the constructed bridge. Neither a small spectral tail in Frobenius norm nor a small realization error pays this joint energy. The original favorable-flatification recurrence remains open in this work.

## Central-shell obstruction to a fixed-pair MGF union certificate

Even the exact flat Gaussian-sign covariance improvement cannot close the parent target by the shellwise independent-pair Chernoff/union certificate alone. For every full signing, uniform Boolean x satisfies

    E H_A(x)^2=binom(n,2).

Thus, for fixed epsilon>0, the shell |H_A(x)|<=epsilon n^(3/2) contains at least (1-O(1/(epsilon^2 n)))2^n spins. Its entropy per coordinate tends to log 2, for every child, without any optimality assumption.

In the flat opposite-spectral Gaussian bridge, the exact variance proxy is 1-(8rho/pi)e_A e_D+o(1). On the two central shells it tends to 1 as epsilon tends to zero. Therefore the shellwise Chernoff/union criterion at parent coefficient L requires, in this limit,

    L^2/2 > 2 log 2,   equivalently L>2sqrt(log 2)=1.665109...

The desired comparable equal-child coefficient is L=2sqrt(2)c, at most sqrt(2) for c<=1/2, and approximately 1.396 for the active upper. Hence even the central shell fails that certificate, regardless of the favorable variance reduction near high-energy child states.

This is not an actual Gaussian-bridge cap lower bound: the exponentially many events are correlated. It identifies the precise need for a joint landscape, covering, or other dependence-sensitive bound beyond the proved fixed-pair MGF, not merely sharper high-energy shell counts.

## Subsequent positive advance and refined gap

The newly completed `flatify_independent_2026_09_07_gaussian_sign_quenched_universality.md`, independently audited in `flatify_construct_2026_09_07_quenched_gaussian_universality_audit.md`, DOES supply the dependence-sensitive quenched transfer missing from the preceding fixed-pair discussion. It compares the actual Gaussian-sign parent expected cap, with arbitrary unchanged child energies, to its covariance-matched Gaussian parent at normalized error O(n^(-1/6)sqrt(log n)) for fixed ellipticity and operator constants.

Thus the central-shell paragraph remains only a limitation of the older union certificate, not a limitation of this new theorem. The current gap is evaluating the resulting covariance-matched Gaussian parent optimization favorably. In the specific opposite-spectral law normalized by p=||A||op||D||op, `flatify_construct_2026_09_07_large_denominator_gaussian_bridge.md` proves an additional limitation: when p/n diverges, that Gaussian parent is iid-bridge equivalent in expected normalized cap, at error O(sqrt(n/p)). This law needs a genuinely leading covariance effect, or a different normalization/construction, to improve the expected-cap bound.
