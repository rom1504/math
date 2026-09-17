# All-energy Gaussian-edge stability after exact-sign star regularization

2026-09-17. Localization-track extension of the director's actual-sign
random-star regularizer. This is a finite, quantitative composition theorem,
not a claim of fixed-near-level convergence or a new asymptotic cap constant.
The base matrix is an actual hollow full signing. Its Gaussian perturbations
are real weighted matrices, and are NOT represented as new full signings.

## 1. The regularizer supplies a whole energy-width profile

For a hollow full signing W of order N write

    H_W(x)=sum_(i<j) W_ij x_i x_j,
    Q(W)=max_x |H_W(x)|,
    E_W(T)={x: Q(W)-|H_W(x)|<=T}.

Let 1<=q, m=N-q, q^3<=m. For every original full signing A, the
[same-order star construction](paper_director_sparse_random_regularization_2026_09_17.md)
can be selected so that, with universal constants C,c,

    Q(W)<=Q(A)+C N sqrt(q),                            (1)
    w(E_W(T))<=C[N/sqrt(q)+sqrt(N T) q^(-1/4)]
                         for EVERY T>=0.              (2)

Only edges incident to the q rewritten vertices change. The statement is
simultaneous in T, and includes the full absolute code and all new spins.

Here is the additional argument needed beyond the one-window theorem.
Retain its notation F_y, Delta_(y,s), and put T0=N/sqrt(q). For a fixed
secondary variance s>0, uniformly over every one of the 2^q patterns y,

    E Delta_(y,s)<=kappa m s/(2sqrt(q/m))+6sqrt(m/q),
    Pr{Delta_(y,s)>E Delta_(y,s)+4sqrt(m/q)}<=exp(-2m/q^2).

Choose T_j=2^j T0, for 0<=j<=ceil(log_2 q), and

    s_j=T_j sqrt(q)/m^2.

The old-spin near-window is b_j=(T_j+2Q(D))/sqrt(m).
The secondary-field inequality is pointwise:

    sqrt(s_j) w(C_(y,j))<=b_j+Delta_(y,s_j).

Since Q(D)<=q^(3/2), q^2<=m, and T_j>=T0, the right side divided by
sqrt(s_j) is at most C sqrt(m T_j) q^(-1/4). The Gaussian width of
the union over all patterns adds at most sqrt(2mq log 2); the q free
new coordinates add kappa q. Both are absorbed in C N/sqrt(q).
The total failure probability is at most

    [1+ceil(log_2 q)] exp[-(2-log 2)m/q^2].             (3)

This tends to zero uniformly in the allowed range. The cap event has
failure exp(-c m), so the two events intersect. Monotonicity fills the
gaps between the T_j. Below T0 use E(T) subset E(T0). Above N sqrt(q)
the trivial w<=kappa N proves (2). Thus there is no uncontrolled
all-energy region in (2). A harmless extra top dyadic point covers
rounding of q and m.

## 2. From spin width to signed quadratic width

Let g=(g_ij)_(i<j) consist of independent standard Gaussians, and write

    Z_(sigma,x)=sigma sum_(i<j) g_ij x_i x_j,
    d_(sigma,x)=Q(W)-sigma H_W(x),  sigma in {+-1}.

The signed index set with d<=T projects into E_W(T). For a fixed
polarity, if x,y differ in r coordinates then

    E[H_g(x)-H_g(y)]^2=4r(N-r)<=N ||x-y||_2^2.

Gaussian comparison therefore gives expected supremum at most
sqrt(N) w(E_W(T)) in each polarity. The two-polarity union costs
at most sqrt(2 binom(N,2) log 2). Hence (2) implies

    W_edge(T):=E max_(d<=T) Z
       <=C[N^(3/2)/sqrt(q)+N sqrt(T) q^(-1/4)].        (4)

This comparison does not identify spin width with edge width. The
factor sqrt(N), and the separate polarity cost, are essential.

For completeness, the finite Gaussian comparison used here follows by
interpolating the two Gaussian vectors in the expected soft maximum.
The derivative is beta/4 times the Gibbs average of the difference of
their pairwise increment variances. Its sign gives the comparison;
letting beta tend to infinity gives the displayed maximum inequality.
Thus matching individual variances is not required.

## 3. A finite-dimensional perturbation class, not just one noise sample

Fix an integer r>=1. Independently of W, draw r independent Gaussian
edge matrices G^1,...,G^r. For theta in R^r put

    G(theta)=sum_(ell=1)^r theta_ell G^ell.

There are universal constants C,c such that the following holds for
every radius a>0. Put

    U=C[a N^(3/2)/sqrt(q)+a^2 N^2/sqrt(q)+a N sqrt(r)]. (5)

Then

    E sup_(||theta||<=a) |Q(W+G(theta))-Q(W)|<=C U.     (6)

More importantly, for every T>=U, with probability at least

    1-C exp[-c T^2/(a^2 N^2)],                        (7)

EVERY maximizer of EVERY Q(W+G(theta)), ||theta||<=a, belongs to the
single deterministic code E_W(T). Thus on this event,

    Q(W+G(theta))=max_(x in E_W(T)) |H_(W+G(theta))(x)|
                          for every ||theta||<=a.     (8)

The code is chosen before the Gaussian matrices. Its size need not be
polynomial and this is not an efficient enumeration claim.

### Proof with a ground-state anchor and all-energy shells

Choose a signed absolute ground state i0=(sigma0,x0), and let v_i be
the edge query vector (sigma x_i x_j)_(i<j). Set a_i=v_i-v_i0;
then ||a_i||<=2sqrt(d), d=binom(N,2). For each signed state define

    Y_i=(G^1 dot a_i,...,G^r dot a_i).

For any T, compare the centered Gaussian process u dot Y_i indexed
by d_i<=T and ||u||=1. Its increment variance obeys

    ||u tensor a_i-v tensor a_j||^2
        <=2||a_i-a_j||^2+8d||u-v||^2.

Comparison with sqrt(2) g dot a_i+sqrt(8d) z dot u gives

    E max_(d_i<=T) ||Y_i||
       <=sqrt(2) W_edge(T)+sqrt(8d) E||z||
       <=A+B sqrt(T),                                (9)

where A=C N^(3/2)/sqrt(q)+2N sqrt(r) and
B=C N q^(-1/4). This supremum is 2sqrt(d)<=sqrt(2)N-Lipschitz
in all the Gaussian edge coordinates together. Consequently its
upper deviation by v has probability at most exp[-v^2/(4N^2)].

An outside state can defeat the anchored old ground state for some
||theta||<=a only if a||Y_i||>=d_i. For d_i in [2^j T,2^(j+1)T],
apply (9) to the larger level 2^(j+1)T. If

    T>=16[a A+a^2 B^2],                              (10)

then a times that expected supremum is at most 2^(j-1)T for every
j>=0. Concentration and a union over all j give (7). In particular
the proof treats ALL energies, not only the small-code window.
Empty high shells cause no issue, and extending the sum to infinity
only weakens the bound. Since A>=N, the geometric tail is bounded by
a universal multiple of its first term.

For the cap bound put

    R=sup_i [a||Y_i||-d_i]>=0.

The same shell proof, including the inner level d_i<=u, gives
Pr{R>u}<=C exp[-c u^2/(a^2 N^2)] whenever
u>=16[a A+a^2 B^2]. Integration gives ER<=C[a A+a^2 B^2].
The anchor itself contributes at most
a||(G^ell dot v_i0)_ell||, whose mean is <=a sqrt(d r).
This proves (6). For a single fixed noise direction the anchor has
mean zero, and Jensen gives the sharper sign statement

    0<=E Q(W+aG)-Q(W)
       <=C[a N^(3/2)/sqrt(q)+a^2 N^2/sqrt(q)].         (11)

All assertions extend to a=0 trivially.

## 4. A cap-stable finite witness at a nontrivial perturbation scale

Fix L<infinity. Take a=L/sqrt(N) and r<=N/q. Then U<=C_L T0,
T0=N/sqrt(q), and the deterministic witness E_W(C_L T0) has

    w<=C_L N/sqrt(q),
    log cardinality<=C_L N/sqrt(q) sqrt(log(e q)).     (12)

The entropy bound is the Gaussian-information/decoding theorem in
[the external-field artifact](paper_localization_external_field_regularization_2026_09_17.md),
applied to (2). With probability at least 1-C_L exp(-c_L N/q), this
ONE code contains every optimizer throughout the Euclidean radius
L/sqrt(N) ball in the r-dimensional random Gaussian edge subspace.
The expected uniform cap variation is O_L(N/sqrt(q)).

For q of order N^(1/3), these statements read:

    exact-sign preparation cost       O(N^(7/6)),
    random perturbation dimension     r<=c N^(2/3),
    coefficient-space radius          L N^(-1/2),
    uniform cap response / witness    O_L(N^(5/6)),
    log witness size                  O_L(N^(5/6) sqrt(log N)),
    witness failure probability       exp[-c_L N^(2/3)].

For comparison, a radius N^(-1/2) Gaussian edge perturbation has
unrestricted supremum of order N. The localized response above is
smaller by N^(1/6) at the optimized q. It is not merely a restatement
of small entropy: the old energy deficits exclude the exponentially
large complement, uniformly over a continuum of perturbations.

## 5. Scope and remaining limitation

The original full signing is kept exactly on its large principal block,
and the prepared matrix W is a genuine full signing at the SAME order.
The weighted Gaussian neighborhood is an analytic perturbation class,
not a rounding theorem and not a family of new actual-sign competitors.
The code and Gaussian subspace may depend on W in the stated order:
first prepare W for A; then form its deterministic energy code; then
draw the independent Gaussian subspace. No uniformity over every
possible original signing with one shared draw is asserted.

The witness window shrinks after N^(-3/2) normalization. This does not
give control at every fixed eta, an exact low-response center cover,
or a cap-stable macroscopic vertex extension. In particular the
[Hadamard extension obstruction](paper_discrepancy_extension_obstruction_2026_09_17.md)
is compatible with (8): small Gaussian-edge susceptibility at this
variable precision does not remove the all-sign-column mean cost.

Audit status: the discrepancy researcher independently reconstructed
Sections 1--4 and 6 in full and returned PASS, including the dyadic
profile, quadratic increment comparison, signed polarity union,
tensor-index Gaussian comparison, all-energy exclusion, and uniform-ball
constants. The Bernoulli researcher independently read and reconstructed
Sections 1--7 in full and returned PASS, including the partition tail,
signed-state factor two, and doubled-radius half-deficit certificate.
The finite replay
`computations/paper_localization_2026_09_17_gaussian_edge_stability.py`
passes 4,500 edge metric identities, 4,500 tensor comparison inequalities,
14,400 shell budgets, and 56 exact finite witness certificates with
5,600 tested perturbations. This is an algebra/certificate replay, not
an empirical justification for the Gaussian concentration theorem.

## 6. The exchangeable-information preparation: arbitrary q, logarithmic loss

The director subsequently proved the stronger actual-sign preparation
[from exchangeable incident-edge information](paper_director_exchangeable_sign_regularization_2026_09_17.md).
For every 1<=q<=N/2 it supplies a same-order full signing W with

    Q(W)<=Q(A)+C N sqrt(q),
    b(E_W(T))<=T+B for every T>=0,  B=C N/sqrt(q),      (13)

where b is Bernoulli width. The SAME conditional one-vertex increment
controls every T. The localization researcher independently read and
reconstructed the complete overlapping-row information proof: each
independent edge occurs in at most two new-vertex rows, the uniform
optimizer is permutation-equivariant, and deletion plus entropy
selection gives the increment bound. No q^3 restriction remains in (13).

This gives a complementary weighted-stability regime. Put L_q=log(e q).
There are universal c,C>0 such that, for

    a<=c/sqrt(N L_q),    r<=N L_q/q,                 (14)

one deterministic code E_W(C B) contains every absolute ground state
throughout the independent Gaussian edge ball ||theta||<=a, except
with probability at most

    C exp[-c N L_q/q].                               (15)

Its logarithmic cardinality is at most C B L_q, and the expected
uniform cap variation is at most C B. Thus q may now be any diverging
o(N), at a sqrt(log q) loss in the permitted perturbation amplitude.
This does not supersede Section 4's larger amplitude in its q^3 range.

### Proof

For any nonempty Boolean code C with Bernoulli width b, Gaussian
coordinate clipping gives

    w(C)<=C b sqrt(log(eN/b)),                       (16)

with the zero-width singleton case interpreted separately. To see this,
clip each Gaussian magnitude at R. Conditional on those magnitudes,
coordinatewise convexity and sign symmetry show that the clipped
expected supremum is at most R b. The remaining contribution is at
most N E(|G|-R)_+. Take R of order sqrt(log(eN/b)). The Gaussian
tail bound gives (16). This proof does not assume that C is symmetric.

The function u sqrt(log(eN/u)) is increasing on [0,N]. Combining
(13) with (16), and using the trivial w<=N when T+B>=N, gives

    w(E_W(T))<=C sqrt(L_q)(T+B) for EVERY T>=0.

Section 2 then gives W_edge(T)<=D(T+B), D=C sqrt(N L_q);
the polarity cost O(N) is absorbed since sqrt(N)B>=N.
The Gaussian-subspace anchor process from Section 3 has expected
level supremum at most

    A+D T,   A=C sqrt(N L_q)B+2N sqrt(r).             (17)

Choose the constant c in (14) so that aD<=1/8. If T>=8aA,
the mean on the level 2^(j+1)T is at most 2^(j-1)T. Exactly the same
all-energy shell argument proves exclusion with failure
C exp[-c T^2/(a^2N^2)] and expected uniform cap variation at most
C aA. Under (14), aA<=C B; take a sufficiently large constant
multiple of B for the witness window. This proves (15) and the cap
claim. Entropy follows directly from the director's VC/Sauer estimate
log|C|<=b(C) log(eN/b(C)), not from an unproved width equivalence.

For example q=floor(N/log N) gives

    preparation cost                 O(N^(3/2)/sqrt(log N)),
    witness window / cap response    O(sqrt(N log N)),
    log witness size                 O(sqrt(N)(log N)^(3/2)),
    random perturbation dimension    r<=c(log N)^2,
    coefficient-space radius         c/sqrt(N log N),
    witness failure                  exp[-c(log N)^2].

The first cost vanishes after N^(3/2) normalization. All-energy
exclusion, not just low entropy of a selected set, is used in the
weighted-neighborhood conclusion. As before, it does not provide an
actual-sign macroscopic extension or a fixed-eta convergence theorem.

## 7. The same finite witness captures the full low-temperature partition

This is an additional all-energy consequence, rather than an inference
from ground-state containment alone. In Section 4 put

    B=N/sqrt(q),   Lambda=sqrt(log(e q));

in Section 6 put B=C N/sqrt(q), Lambda=log(e q). In the respective
regimes, both preparations satisfy

    log|E_W(D)|<=C Lambda D for every D>=B.            (18)

For Section 4, write its spin width bound as C[B+sqrt(BD)].
For D>=B, apply the information width-to-entropy inequality and use
sqrt(BD)>=B, N/B=sqrt(q), and sqrt(BD)<=D. If the displayed width
bound exceeds N, the trivial entropy bound proves (18) instead.
For Section 6, (18) follows immediately from b(E(D))<=D+B and
the VC/Sauer estimate. Thus (18) includes large-energy shells as well.

Increasing the universal constant in the witness window T=C_L B,
the same Gaussian event used above can ensure the stronger assertion

    |G(theta) dot (v_i-v_i0)|<=d_i/2
       whenever d_i>=T, for EVERY ||theta||<=a.       (19)

Indeed use twice the noise radius in the shell certificate; for
Section 6 choose its amplitude constant correspondingly smaller.
The probability bounds and orders of all parameters are unchanged.

For the absolute partition functions define

    Z_beta(theta)=sum_x exp(beta |H_(W+G(theta))(x)|),
    Z_beta^C(theta)=sum_(x in E_W(T))
                            exp(beta |H_(W+G(theta))(x)|).

There is a universal C such that, for beta>=C Lambda, on event (19)

    0<=[Z_beta(theta)-Z_beta^C(theta)]/Z_beta^C(theta)
                <=2 exp(-beta T/4)                   (20)

SIMULTANEOUSLY for every allowed theta. To prove it, compare every
outside signed energy to the same perturbed old ground-state anchor.
By (19) its loss is at least d_i/2. A dyadic shell with
2^j T<=d_i<2^(j+1)T has at most 2|E_W(2^(j+1)T)| signed states.
Equation (18) and beta>=C Lambda make its total relative contribution
at most exp[-beta 2^j T/4]. Sum this geometric-in-the-exponent series.
The truncated partition contains the anchor term, so the relative
bound follows. Summing both polarities only enlarges the outside
absolute partition, and costs the already included factor two.

Consequently the free energies differ by at most
2 beta^(-1) exp(-beta T/4), uniformly in the entire perturbation ball.
Also beta^(-1)log Z_beta(theta) differs from its exact cap by O_L(B)
at beta of order Lambda. The Gaussian parameters and the code are
unchanged as beta ranges over this entire low-temperature interval.
No efficient enumeration of the subexponential witness is asserted.
