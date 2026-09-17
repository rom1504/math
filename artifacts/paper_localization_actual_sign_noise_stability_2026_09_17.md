# Actual-sign noise stability from cloned preparation

2026-09-17. Localization-track theorem; director independently reconstructed
the mechanism and the discrepancy researcher read Sections 1--6 and the
supporting Gaussian-edge Section 8 in full and returned PASS. This is a random SAME-ORDER
full-sign neighborhood, not a uniform statement about all flip patterns
and not a new-vertex extension. External priority is unestablished.

## 1. Quantitative statement and order of choices

Let eta=eta_N tend to zero, with

    eta>=N^(-2/3),
    N^(-1/6) eta^(-1/2) [log(e/eta)]^(-1/3) ->0.       (1)

For example eta>>N^(-1/3) suffices. Write ell=log(e/eta). Given ANY
full hollow signing A of order N, first prepare the actual full signing W
using the [cloned-block construction](paper_director_cloned_block_regularization_2026_09_17.md).
It has

    Q(W)<=Q(A)+C sqrt(eta)N^(3/2),
    b(E_W(D))<=D/r+C B for every D>=0,
    r=floor(sqrt(N)eta^(3/4)),
    p=floor(eta^(-1/2)),  B=N/sqrt(p).                (2)

Choose a sufficiently small universal c_0>0 and set

    rho=c_0 eta^(3/2)/ell,
    T=C_0 eta N^(3/2),  C=E_W(T),                    (3)

where C_0 is a sufficiently large universal constant. Both W and C
are deterministic before any flips are drawn. Independently flip each
edge of W with probability rho, obtaining an ACTUAL full signing W^rho.
Then, for all sufficiently large N,

    Pr{every absolute maximizer of W^rho belongs to C}
       >=1-exp[-c N sqrt(eta)ell],
    log|C|<=C N eta^(1/4)ell=o(N).                    (4)

There is also the cap bound

    E|Q(W^rho)-(1-2rho)Q(W)|<=C eta N^(3/2).          (5)

With the same order of failure probability, the absolute difference
in (5) is at most C eta N^(3/2). If Q(A)<=L N^(3/2), then
Q(W)<=C_L N^(3/2), and the factor (1-2rho) can be removed at a further
cost 2rho Q(W)=o(eta N^(3/2)). Thus (5) becomes actual cap stability
around Q(W), not merely around its contracted mean landscape.

The proof below works also for any smaller positive rho for which
N^(4/3)rho^(1/3)=o(T); its explicit failure exponent is
c T^2/(rho N^2+T). Fixing (3) is convenient for a single clean statement.
The parent matrix W depends on the prescribed eta sequence. No common
preparation for every vanishing precision is asserted.

## 2. Signed anchors, mean contraction, and the escape witness

Index signed states by i=(sigma,x), sigma in{+1,-1}, and define

    v_i=(sigma x_u x_v)_(u<v),
    d_i=Q(W)-W dot v_i>=0.

Fix one old signed ground state i_0, so d_(i_0)=0. Put a_i=v_i-v_(i_0);
each coefficient of a_i is in{-2,0,2} and ||a_i||_2<=sqrt(2)N.
Write s=1-2rho>=1/2 and

    W^rho=s W+Z,
    Z_e=-2 W_e(B_e-rho),  B_e iid Bernoulli(rho).      (6)

Then E Z_e=0, Var Z_e=v=4rho(1-rho),
E|Z_e|^3<=8rho, and changing B_e changes Z_e by magnitude2.

The nonnegative escape functional is

    F(z)=max(0, max_(d_i>=T)[z dot a_i-s d_i/2]).      (7)

If any signed maximizer of W^rho has d_i>=T, comparison to i_0 gives
Z dot a_i>=s d_i. Hence F(Z)>=sT/2>=T/4. This implication includes
both polarities, all old words, and even signed states whose wrong
polarity has a large deficit. An absolute maximizer outside E_W(T)
necessarily produces such a signed state. Controlling all d_i>=T
is slightly stronger than needed and handles boundary equalities.

## 3. Gaussian expectation bound: every energy shell is paid

Let G have independent standard Gaussian edge coordinates. Section 8 of
[Gaussian-edge stability](paper_localization_gaussian_edge_stability_2026_09_17.md)
gives the anchored all-energy profile

    E max_(d_i<=D)|G dot a_i|
       <=C sqrt(N ell)(D/r+B).                       (8)

To justify the absolute value, take both signs of the scalar noise
direction in the product-index comparison; the additional O(N) term
is absorbed by sqrt(N ell)B. The supremum is sqrt(2)N-Lipschitz.
For sufficiently small c_0 and large C_0 in (3), (8) implies on every
shell 2^jT<=d_i<2^(j+1)T that

    sqrt(v) E max_(d_i<2^(j+1)T)|G dot a_i|
                   <=s 2^jT/4.                     (9)

Indeed sqrt(v)C sqrt(N ell)/r is an arbitrarily small universal
constant, and sqrt(v)C sqrt(N ell)B<=c T. Gaussian concentration,
now integrated rather than just asserted with high probability, gives

    Pr{F(sqrt(v)G)>u}
       <=sum_(j>=0) exp[-c(2^jT+u)^2/(vN^2)],
    E F(sqrt(v)G)
       <=C sqrt(v)N exp[-cT^2/(vN^2)].               (10)

The sum may extend to infinity because empty high shells contribute
nothing. The ratio T/(sqrt(v)N) diverges, so the first shell dominates
and integrating the Gaussian tails proves the second line. This step
is essential: a high-probability Gaussian containment statement alone
would not justify comparing expectations of (7).

## 4. Uniform all-offset replacement in expectation

For any finite family of affine functions h_j+c_j dot z, with
|c_je|<=2 and K members, use their soft maximum with inverse
temperature lambda. Its uniform approximation error is log(K)/lambda,
and its third derivative in any edge coordinate is bounded by
C lambda^2: it is lambda^2 times a third centered moment of a variable
in[-2,2]. Sequentially replace independent Z_e by independent
sqrt(v)G_e. Means and variances match, so the first two Taylor terms
cancel. The third absolute moments of both variables are at most C rho.
With at most N^2/2 edges and log K<=C N,

    |E max_j(h_j+c_j dot Z)
       -E max_j(h_j+c_j dot sqrt(v)G)|
       <=C[N/lambda+lambda^2 N^2rho]
       <=C N^(4/3)rho^(1/3).                        (11)

The final choice is lambda=(N rho)^(-1/3). There is no restriction
on the offsets, their signs or sizes, the spectrum of W, or the
number of near-ground states beyond the full cube count. Adding the
zero affine function includes the positive part in (7). Applying
(10)--(11), then (1)--(3), gives

    E F(Z)<=C N^(4/3)rho^(1/3)+exponentially small
             =o(T).                                (12)

At (3), the ratio of the displayed error to T is precisely of order
N^(-1/6)eta^(-1/2)ell^(-1/3). This is the only reason for (1).

## 5. Bernoulli martingale concentration upgrades Markov

Changing one B_e changes every affine term of (7) by at most4, so
F is coordinatewise4-Lipschitz. Reveal the independent B_e one by
one. Conditional on all previous coordinates, let d_e be the difference
between the two future-averaged values with B_e=1 and B_e=0.
Then |d_e|<=4 and the Doob martingale increment is exactly

    (B_e-rho)d_e.

It has absolute value at most4 and conditional variance at most16rho.
Thus the total predictable variance is at most8rho N^2. The elementary
bounded-increment exponential-moment estimate gives, for t>=0,

    Pr{F-EF>=t}
       <=exp[-t^2/(16rho N^2+(8/3)t)].               (13)

For completeness, if a centered variable U obeys |U|<=b, Taylor's
series and |E U^k|<=b^(k-2)E U^2 imply
log E e^(lambda U)<=lambda^2 E U^2/[2(1-lambda b/3)]
for 0<=lambda<3/b. Conditional iteration and optimizing lambda give
the usual denominator2(V+bt/3), which yields (13). No Gaussian
approximation is involved in this concentration step.

By (12), eventually EF<=T/8. Escape forces F>=T/4, so (13) proves

    Pr{escape}<=exp[-c T^2/(rho N^2+T)].              (14)

Under (1)--(3), rho N^2/T is of order sqrt(N eta)/ell and diverges.
The exponent in (14) is therefore at least c N sqrt(eta)ell,
proving (4). The stronger concentration step was supplied independently
by the director; using only Markov would give a weaker o(1) failure.

## 6. Cap response and exact-support scope

Apply (11) to the unrestricted signed cap, now with coefficients
v_i in{-1,1} and offsets s W dot v_i. The Gaussian stability theorem
at amplitude sqrt(v)/s gives

    0<=E Q(sW+sqrt(v)G)-sQ(W)<=C T.

Jensen gives the analogous nonnegative mean increment for Z. Equation
(11) and its o(T) error therefore yield
0<=E Q(W^rho)-sQ(W)<=C T. Alternatively, the anchored residual
max_i[Z dot a_i-sd_i] is nonnegative and has expectation at most C T.
The anchor Z dot v_(i_0) has absolute mean at most sqrt(v)N, absorbed
by T. This proves (5). The same Bernoulli martingale calculation for
Q(W^rho), with coordinate Lipschitz constant2, gives its two-sided
concentration on scale T with failure exp[-cT^2/(rho N^2+T)].

Every output coefficient of W^rho is exactly plus or minus1. This
is a genuine random full-sign perturbation theorem with a fixed finite
witness, unlike the weighted Gaussian statement used in its proof.
It is not uniform over all subsets of flipped edges; it does not
produce a low-response bridge-column law or favorable new-child value;
and its window still shrinks after normalization by N^(3/2).

The director's [direct mixed-tail chaining proof](paper_director_mixed_tail_sign_stability_2026_09_17.md)
now gives the wider domain sqrt(N)eta^(3/4)>=C log(e/eta). The
localization researcher independently read Dirksen's actual Theorem 3.5
and complete proof, then reconstructed that full transfer and returned
PASS. It pays the jump geometry separately instead of a global Gaussian
replacement error. The elementary proof above remains independent of
that improvement and uses precisely the domain (1).

Finite replay: `computations/paper_localization_2026_09_17_actual_sign_noise_stability.py`
passes 4,384 exact full-sign flip worlds, 2,508 signed escape implications,
4,372 conditional martingale nodes, 120 deterministic Gaussian half-deficit
certificates, 28 integer clone-parameter checks, and8,128 binomial-mode
conditioning checks. The output is in
`tmp/paper_portfolio_2026_09_17/localization/actual_sign_noise_stability_audit.json`.
These checks validate the finite identities, not the asymptotic probability
theorem, whose proof is above.

## 7. A single witness for an entire random actual-sign Hamming path

This is a further consequence of the director's mixed-tail theorem,
not a strengthening of the Lindeberg estimate (11). Use its wider domain

    eta=o(1),  sqrt(N)eta^(3/4)>=C log(e/eta),
    rho_max=c eta^(3/2)/ell,
    T=C eta N^(3/2),  U=N sqrt(eta)ell.

Prepare the same W and deterministic code C=E_W(T). Let d=binom(N,2),
choose a UNIFORM RANDOM PERMUTATION of the d edges, and let W_k be
obtained by flipping precisely the first k edges. Then with probability
at least1-C exp(-cU), SIMULTANEOUSLY for every integer
0<=k<=floor(rho_max d),

    every absolute maximizer of W_k belongs to C,
    |Q(W_k)-(1-2k/d)Q(W)|<=C T.                     (15)

If Q(A)<=L N^(3/2), the second line also holds with Q(W) in place
of (1-2k/d)Q(W), after changing C to C_L. In particular this is a
whole path of distinct actual signings, not only one endpoint or a
weighted interpolation. It still does not cover all subsets of edges.

Here are the details of the two uniformity steps. For every deterministic
rate 0<rho<=rho_max, the director's proof has the same complexity
slopes and intercepts bounded by those at rho_max; its gamma_1 term
does not depend on rho. Therefore its failure probability is uniformly
at most C exp(-cU). Decrease its absolute amplitude constant so that
the all-energy argument yields half-deficit rather than just optimizer
exclusion. The same mixed-tail bound on the INNER level d_i<=T gives
sup|Z dot a_i|<=C T with that failure, since its complexity intercept
is at most C T and its diameter budgets at u=cU are also at most C T.
Scalar Bernstein controls the anchor |Z dot v_(i_0)|<=C T. The
nonnegative anchored residual is at most C T: the inner level is
bounded as above, and all outside levels have negative residual by
the half-deficit event. This proves the cap envelope in (15) for each
independent Bernoulli flip law, with the same uniform failure bound.

For 1<=k<d, independent flips at rate rho=k/d, conditioned on having
exactly k flips, are a uniform k-subset. The integer k is a mode of
Binomial(d,k/d), as is checked by the adjacent probability ratios;
hence its probability is at least1/(d+1). Conditioning can therefore
increase failure by at most d+1. The first k edges of a uniform random
permutation have precisely this uniform k-subset law. A union bound
over all k<=floor(rho_max d) costs at most one further factor d+1.
Since U/log N tends to infinity throughout the stated domain, this
polynomial factor is absorbed in C exp(-cU). The k=0 state is
deterministic and already lies in its own witness. No independence
between different times along the path is required.

## 8. The same actual-sign path has a finite low-temperature witness

The path event can be selected to control more than its maxima. With
the same parameter orders and failure probability as Section 7, define

    Z_beta(W_k)=sum_x exp(beta |H_(W_k)(x)|),
    Z_beta^C(W_k)=sum_(x in C) exp(beta |H_(W_k)(x)|).

Then, SIMULTANEOUSLY for every allowed k and every beta>=C ell/r,

    0<=[Z_beta(W_k)-Z_beta^C(W_k)]/Z_beta^C(W_k)
                         <=2 exp(-beta T/8).         (16)

This is a statement about actual full signings throughout a random
path. The witnesses need not be computationally enumerable in polynomial
time, and the threshold is a low-temperature threshold at the indicated
N-dependent normalization.

For the proof, decrease the absolute flip-rate constant once more so
that the mixed-tail shell event ensures

    |Z dot a_i|<=s d_i/2 for every d_i>=T,
    s=1-2k/d>=1/2.

An outside signed state's perturbed energy is then below the perturbed
old anchor by at least s d_i/2>=d_i/4. The clone profile and VC/Sauer
bound imply

    log|E_W(D)|<=C(ell/r)D for all D>=T.

There are at most2|E_W(2^(j+1)T)| signed states with
2^jT<=d_i<2^(j+1)T, including every wrong-polarity state that is
needed in the absolute partition. For beta>=C ell/r their total
relative weight is at most exp(-beta 2^jT/8), after enlarging the
universal constant to pay the factor2. The sum over all j is bounded
by2 exp(-beta T/8); beta T diverges in the present domain. The
truncated partition contains the anchor word, so this is a relative
bound with the correct denominator. It holds on one event for every
beta in the specified interval. Section 7's conditioning and prefix
union apply unchanged to this stronger event.

The full and truncated free energies consequently differ by at most
2 beta^(-1)exp(-beta T/8). At beta of order ell/r their common
free energy differs from the exact cap by O(T), because
log|C|=O((ell/r)T). No exponentially large high-energy complement
has been omitted without a deficit payment.

The discrepancy researcher independently read Sections 7--8 in full
and returned PASS, including the uniform rate budgets, modal conditioning,
random-prefix union, scalar cap anchor, signed partition count, and
simultaneous temperature interval.

## 9. Archive comparison and the retained limitation

The closest earlier edge-noise item inspected is
`fresh_edge_noise_entropy_obstruction_2026_09_05.md`. It proves an
annealed one-state entropy obstruction to improving the half-cap family
by independent edge noise. The present theorem does not contradict or
replace it: rho tends to zero, preparation is paid, and no fixed cap
improvement is claimed. Its additional output is an all-energy,
correlation-sensitive finite witness for a random sign neighborhood.

The earlier `principle_construct_2026_09_07_universal_stability_audit.md`
proves a fixed exponential entropy deficit for low-total-violation words
of bounded-cap actual matrices, uniformly over an external field. It
does not give a subexponential code at a selectable shrinking tolerance,
nor containment under a fresh noise path. These are distinct conclusions.
This comparison establishes no external-priority claim.

Numerically in scales, the certified preparation costs sqrt(eta)N^(3/2),
whereas the witness tolerance is eta N^(3/2) and the mean contraction
is only of order eta^(3/2)N^(3/2)/ell for a bounded-cap input. Thus one
cannot turn the mean contraction into a free favorable parent slope by
ignoring preparation or the finite witness tolerance.
