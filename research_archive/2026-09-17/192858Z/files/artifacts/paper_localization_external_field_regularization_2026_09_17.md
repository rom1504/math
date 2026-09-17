# External-field regularization of full Boolean near-level sets

2026-09-17. The Gaussian regularization proposal is due to the director.
The localization track independently reconstructed its constants,
absolute-polarity reduction, a Gaussian-information entropy conversion,
and the scalar-sign-field transfer below. These are finite theorems for
arbitrary offsets. The physical duplicated-column compiler is maintained
by the director; this artifact does not claim that regularization alone
proves convergence of the original minimization problem.

## 1. Gaussian near-level width: expectation and exponential probability

Let H be ANY real function on {+-1}^n. For t>0 put

    M_t(g)=max_x [H(x)+sqrt(t) g dot x],
    G(t)=E M_t(g),
    C_eta(g)={x: H(x)+sqrt(t)g dot x>=M_t(g)-eta n}.

For a nonempty finite code C write w(C)=E_z max_(x in C) z dot x,
where z is a fresh standard Gaussian vector. Put kappa=sqrt(2/pi).
Then

    E_g w(C_eta(g))/n <= sqrt(2 kappa eta)/t^(1/4).   (1)

More explicitly, for any s,a>0,

    Pr{ w(C_eta(g))/n >
       [eta+kappa s/(2sqrt(t))+a]/sqrt(s) }
                          <=exp[-a^2 n/(8t)].        (2)

In particular s=4eta sqrt(t)/kappa and a=eta give

    Pr{ w(C_eta(g))/n > 2sqrt(kappa eta)/t^(1/4) }
                          <=exp[-eta^2 n/(8t)].      (3)

These bounds are uniform in H, including its magnitude, degeneracies,
and n-dependence. The width tends to zero after normalization by n as
eta tends to zero at every fixed positive field strength.

### Proof by a second independent field

Define the nonnegative conditional increment

    Delta_s(g)=E_z max_x[H(x)+sqrt(t)g dot x+sqrt(s)z dot x]
                      -M_t(g).

Restrict the first maximum to C_eta(g), then average z. Pointwise in g,

    sqrt(s) w(C_eta(g)) <=eta n+Delta_s(g).           (4)

Gaussian addition gives E Delta_s=G(t+s)-G(t). Coupling the fields
with the same Gaussian vector also gives

    0<=G(t+s)-G(t)
       <=kappa n[sqrt(t+s)-sqrt(t)]
       <=kappa n s/(2sqrt(t)).                       (5)

The function Delta_s is 2sqrt(tn)-Lipschitz in g: each of its two
maximum terms is sqrt(tn)-Lipschitz. Gaussian Lipschitz concentration
therefore gives

    Pr{Delta_s>E Delta_s+a n}<=exp[-a^2 n/(8t)].

For completeness that concentration fact follows from the Brownian
martingale representation of a Lipschitz function f(g). Its heat
extension u(r,x)=E f(x+sqrt(1-r)g) has gradient norm at most the
Lipschitz constant L. Ito's formula writes f(B_1)-Ef(B_1) as a stochastic
integral with quadratic variation at most L^2. The exponential
martingale bound gives E exp(lambda[f-Ef])<=exp(lambda^2 L^2/2),
and Chernoff gives the stated tail. Smooth approximation handles a
nonsmooth maximum without changing L.

Combine with (4)--(5) to prove (2)--(3). Averaging (4), dividing by
sqrt(s), and choosing s=2eta sqrt(t)/kappa proves (1).

## 2. Entropy from width by the original Gaussian-observation mechanism

For EVERY nonempty Boolean code C and every sigma>0 such that
sigma w(C)/(2n)<=1/2,

    log|C|/n <= (1/2)log(1+sigma^(-2))
                   +h(sigma w(C)/(2n)).              (6)

This supplies a self-contained entropy conversion without a Sudakov
minoration or an aspect-ratio condition.

Let X be uniform on C and observe Y=X+sigma g. All codewords have
squared norm n, so nearest-code decoding is

    Xhat=argmax_(x in C) Y dot x.

Comparison with the true word gives

    ||Xhat-X||^2 <=2sigma g dot (Xhat-X).

Since g dot Xhat<=max_C g dot x and E g dot X=0,

    E d_H(Xhat,X)<=sigma w(C)/2.

Gaussian channel capacity, with tr Cov(X)<=n, gives
I(X;Y)<=n log(1+sigma^(-2))/2. Coordinatewise binary conditional
entropy and concavity of h bound H(X|Y) by the second term of (6)
times n. Adding proves (6).

In particular set alpha=w(C)/n. Always 0<=alpha<=kappa<1. If alpha>0,
choose sigma=alpha^(-1/3), obtaining

    log|C|/n <= E(alpha),
    E(alpha)=(1/2)log(1+alpha^(2/3))
                           +h(alpha^(2/3)/2).        (7)

Set E(0)=0 by continuity. Thus small normalized width forces small
normalized logarithmic cardinality, quantitatively
E(alpha)=O(alpha^(2/3)log(1/alpha)). No optimality of this envelope is
claimed. In (3), this gives entropy O_t(eta^(1/3)log(1/eta)).

### 2.1 Sharper information bound from posterior width and I--MMSE

The Gaussian-information reconstruction gives a substantially sharper
version of (6): for every sigma with sigma alpha/2<=1/2,

    log|C|/n <= alpha/sigma+h(sigma alpha/2),
    alpha=w(C)/n.                                    (7a)

Here is the complete additional argument. For u>0 observe
Y_u=sqrt(u)X+g, and write m_u(Y_u)=E[X|Y_u]. Finite support gives the
smooth derivative identity

    partial_(y_i) m_(u,i)(y)
                   =sqrt(u) Var(X_i|Y_u=y).

Condition on X and apply Gaussian integration by parts to g. Summing
over i gives

    E g dot m_u(Y_u) =sqrt(u) MMSE(u).

For every realization of (X,g), the posterior mean belongs to conv(C),
so g dot m_u(Y_u)<=max_(x in C) g dot x. Hence

    MMSE(u)<=w(C)/sqrt(u).

The I--MMSE identity reconstructed in
[the primary localization track](paper_localization_2026_09_17.md)
therefore yields

    I(X;Y_u)=(1/2)integral_0^u MMSE(v)dv
                            <=w(C)sqrt(u).           (7b)

The bound is integrable at zero; finite support also permits direct
differentiation and passage to that endpoint. For the channel X+sigma g,
u=sigma^(-2), so I<=w(C)/sigma. Combine this with exactly the same
nearest-code decoder error from Section 2 to prove (7a).

For 0<alpha<=exp(-1), put L=log(1/alpha) and sigma=sqrt(2/L).
The error alpha/sqrt(2L) is below 1/2, and h(z)<=z log(e/z) gives

    log|C|/n
       <=alpha sqrt(2L)
          +alpha[1+(1/2)log(2L)]/sqrt(2L)
       =O(alpha sqrt(log(1/alpha))).                 (7c)

This replaces the weaker envelope (7) whenever alpha is small. Its
order is sharp up to constants: a Hamming ball of radius rho n has
entropy asymptotic to n h(rho), while its Gaussian width is at most
n sqrt(8rho h(rho)). The latter follows by subtracting the fixed
center and applying the Gaussian maximum bound with variance at most
4rho n. Together with (7c), as n tends to infinity and then rho to zero,
these balls have entropy of order alpha sqrt(log(1/alpha)) n.
No literature-priority assertion is made for this entropy/width relation.

## 3. Absolute quadratic energies and both polarities

If H is even, as for a quadratic energy, then

    max_x |H(x)+sqrt(t)g dot x|
       =max_x [|H(x)|+sqrt(t)g dot x].                (8)

Indeed replace x by -x when the outer sign of the absolute value is
negative, and then maximize over the two possible signs of H.
If x belongs to the absolute eta-near-level code on the left, choose
that outer sign sigma. The word z=sigma x belongs to C_eta(g) for the
offset function |H| on the right. Consequently the full absolute
near-level code is contained in C_eta(g) union -C_eta(g).

The latter union adds at most log 2 to entropy. It adds at most
sqrt(2n log 2) to Gaussian width: the two Gaussian maxima have the
same mean w(C), are sqrt(n)-Lipschitz, and their centered exponential
moments are bounded by exp(lambda^2 n/2). The log-sum bound proves
the claimed addition. Thus (1)--(3) retain the same leading asymptotic
width bound for absolute quadratic near-level sets.

## 4. Uniformity over every endpoint pattern of a replicated sign field

Let h_(i,a) be independent fair signs, 1<=i<=n, 1<=a<=K. For every
pattern v in {+-1}^K put

    Z_(v,i)=K^(-1/2)sum_a v_a h_(i,a),
    M_v=max_x[H(x)+delta Z_v dot x],
    C_(v,eta)={x: H(x)+delta Z_v dot x>=M_v-eta n}.

Fix delta>0. Suppose K=K_n tends to infinity and K_n=o(n). Then for
every fixed eta>0, with probability tending to one EXPONENTIALLY in n,

    max_v w(C_(v,eta))/n
                  <=2sqrt(kappa eta)/sqrt(delta)+o_n(1). (9)

The union of all these codes and their antipodes consequently has
normalized entropy at most

    E(2sqrt(kappa eta)/sqrt(delta)+o_n(1))
                           +(K+1)log2/n,             (10)

when its displayed width argument is below one; for larger values the
trivial log 2 entropy bound suffices. Its normalized width is bounded
by the right side of (9) plus sqrt(2(K+1)log2/n). In particular, first
send n to infinity, then eta to zero: both width and entropy vanish.

### Proof: scalar extrema comparison, not distributional approximation

For a fixed v define Delta_(v,s) by adding an independent sqrt(s)
Gaussian field to M_v, averaging that field, and subtracting M_v.
The pointwise inequality (4) remains true. Applying
[the all-offset symmetric-frame theorem](paper_symmetric_frame_universality_2026_09_17.md)
to the independent scalar drivers h_(i,a) compares M_v with the
Gaussian-field maximum at t=delta^2. Each query coefficient has
absolute value delta/sqrt(K), hence its fourth-power sum is

    W4=delta^4 n/K.

There are 2^n queries. The comparison error is therefore at most
C delta n K^(-1/4), with a universal C. The same comparison applies
after adding the extra Gaussian field, by conditioning on that field
and retaining its arbitrary offsets. Thus, uniformly over H and v,

    E Delta_(v,s)/n
       <=kappa s/(2delta)+e_n,
    e_n<=2C delta K^(-1/4)->0.                        (11)

Changing one sign h_(i,a) changes each maximum by at most
2delta/sqrt(K), so Delta_(v,s) changes by at most 4delta/sqrt(K).
Bounded differences gives

    Pr{Delta_(v,s)>E Delta_(v,s)+a n}
                                 <=exp[-a^2 n/(8delta^2)].

Union over ALL 2^K patterns costs only the factor exp(K log2).
Set a=eta and s=4eta delta/kappa. Equations (4) and (11) now give
(9), with the explicit additional term e_n/sqrt(s), and failure
probability at most exp[K log2-eta^2 n/(8delta^2)]. For fixed eta this
tends to zero exponentially because K=o(n). The entropy and width
unions then follow from (7) and the Gaussian log-sum bound.

No near-level set is transferred by weak convergence, total variation,
or pointwise coupling. Only two scalar expected extrema are compared;
the actual sign-field near-level code is controlled by its OWN
conditional increment and concentration. This distinction is necessary.

## 5. Scope and audit

The director's intended compiler repeats each column h^a a number r
of times, so endpoint patterns of the new-spin groups produce a field
of strength r sqrt(K/n) after dividing the original energy by sqrt(n).
Section 4 applies to exactly those endpoint patterns. Controlling the
new principal block and proving that arbitrary new spins reduce to
such endpoints are separate physical steps, not assumed here.

The Gaussian theorem, information conversion, absolute-polarity step,
and sign-field comparison are self-audited here and were circulated
for independent reconstruction. No new minimizer geometry or convergence
conclusion is claimed by these analytic regularization tools alone.

The Bernoulli-track researcher independently read Sections 1--4 and
the polynomial-parameter argument and returned PASS, including the
nearest-code entropy proof, all-offset driver comparison, and parent
window bookkeeping. They separately reconstructed (7a)--(7c) and
returned PASS, including the sharp displayed asymptotic upper constant
sqrt(2) in the information/decoding envelope.

## 6. Independently checked polynomial parameters for the compiler

The director selected

    K=floor(n^(2/3)),       r=floor(n^(1/24)),
    q=Kr=Theta(n^(17/24)),
    delta=r sqrt(K/n)=Theta(n^(-1/8)),
    eta=n^(-1/4),           s=eta delta=Theta(n^(-3/8)).

Integer rounding changes only bounded factors and tends to relative
one. The finite proof of Section 4 applies to these varying parameters;
one must check its errors, not invoke the fixed-parameter limit blindly.
Its expected-increment upper bound is

    E Delta_(v,s) <= kappa n s/(2delta)
                      +O(delta n K^(-1/4))
                   =O(n^(3/4))+O(n^(17/24)).

Take the concentration deviation eta n=n^(3/4). The failure exponent
eta^2 n/(8delta^2) is Theta(n^(3/4)), while the union over endpoint
patterns contributes only K log2=O(n^(2/3)). Consequently, with failure
at most exp[-c n^(3/4)], all endpoint nearcodes obey

    w(C_(v,eta))=O(n^(3/4)/sqrt(s))=O(n^(15/16)).

Let M_joint=max_v M_v and let C_joint be the old-spin code within
eta n of M_joint for at least one endpoint pattern. Since M_joint>=M_v,
C_joint is contained in the union of the individual C_(v,eta). The
width of that union is at most

    O(n^(15/16))+sqrt(2nK log2)=O(n^(15/16)),

because n^(5/6)=o(n^(15/16)). Its antipodal union has the same order.
This per-pattern proof independently checks the director's direct
joint-maximum increment proof; it does not require comparing random
near-level sets under the two driver laws.

Formula (7) now gives the concrete bound

    log|C_joint union -C_joint|=O(n^(23/24) log n).    (12)

The sharper I--MMSE entropy conversion (7c) improves (12) to

    log|C_joint union -C_joint|
                         =O(n^(15/16)sqrt(log n)).   (13)

In physical energy units the near-level window is eta n^(3/2)=n^(5/4).
A new child with cap at most q^(3/2)=O(n^(17/16)) changes that window
by o(n^(5/4)); replacing eta by a bounded multiple leaves all exponents
unchanged. The q new coordinates themselves add at most kappa q to
Gaussian width and q log2 to entropy, also lower order than the bounds
above. These checks support the director's actual-sign compiler while
keeping its separate endpoint and cap arguments explicit.

The finite replay
`computations/paper_localization_2026_09_17_external_field.py` passed
362 exactly known width/entropy family checks, 4800 deterministic
second-field inequalities, 60 absolute-polarity reductions, 60 joint-code
subset checks, and exact rational checks of all displayed polynomial
exponents. Its output is
`tmp/paper_portfolio_2026_09_17/localization/external_field_audit.json`.
No empirical concentration test is used as a proof input.

## 7. Ordinary random-star resampling: a sharper endpoint comparison

The director subsequently observed that replication is unnecessary.
An ordinary iid-sign bridge with q new vertices regularizes the same
microscopic near-level geometry, and q may grow arbitrarily slowly.
The new ingredient is a one-dimensional transport estimate, not a
stronger high-dimensional invariance assertion.

### 7.1 Primary zero-bias lemma and its exact specialization

Goldstein's [Normal Approximation for Hierarchical Structures,
Lemma 2.1](https://dornsife.usc.edu/larry-goldstein/wp-content/uploads/sites/221/2023/06/hie.pdf)
proves W1(W,G)<=2W1(W,W*) for a centered variable W, its zero-biased
version W*, and a Gaussian G of equal variance. The proof uses the
normal Stein equation and its bound ||f''||<=2 for 1-Lipschitz tests.
The original
[Goldstein--Reinert zero-bias paper, Lemma 2.1(5)](https://dornsife.usc.edu/larry-goldstein/wp-content/uploads/sites/221/2023/06/zer.pdf)
constructs W* for independent sums by replacing one variance-weighted
summand; it also explicitly identifies the zero-bias of a fair sign
as uniform on [-1,1]. Both primary proofs were retrieved and read.

For S_q=sum_(a=1)^q epsilon_a, take a uniform random index I and an
independent V uniform on [0,1], and set U=epsilon_I V. Then U is
uniform on [-1,1] and independent of all the other summands. Thus
S_q*=S_q-epsilon_I+U has the required law, while
E|S_q*-S_q|=1/2. The same sharp constant was independently obtained
by the Bernoulli track using the Stein leave-one-out expansion.
Therefore, for EVERY integer q>=1,

    W1(S_q,sqrt(q)G)<=1.                             (14)

If J:R^n->R is separately 1-Lipschitz in every coordinate, independent
coordinate replacement consequently gives

    |E J(S_q^(vector)/sqrt(n))
             -E J(sqrt(q/n)g)|<=sqrt(n).             (15)

Maxima max_x[H(x)+z dot x] have precisely this coordinate Lipschitz
constant, regardless of their offsets. Added independent Gaussian
offsets can be conditioned on. This is why (15) applies to both
endpoints of the conditional increment, including nonsmooth maxima.

### 7.2 All sign-field patterns at once

Let B be an n by q iid-sign bridge. For y in {+-1}^q consider the
scaled old-spin field B y/sqrt(n), of variance t=q/n per coordinate.
Set

    delta=sqrt(q/n),       s=sqrt(q)/n.

Define Delta_(y,s) as before. Equations (5) and (15) yield

    E Delta_(y,s)<=kappa n s/(2delta)+2sqrt(n)
                         <=(kappa/2+2)sqrt(n).       (16)

Changing one bridge entry changes Delta by at most 4/sqrt(n), so
bounded differences gives

    Pr{Delta_(y,s)>E Delta_(y,s)+L sqrt(n)}
                             <=exp[-L^2 n/(8q)].

If 1<=q<=sqrt(n), choose an absolute constant L with L^2/8>log2+1.
Union over every 2^q new-spin pattern then has failure at most
exp[-c n/q]. On the resulting event, (4) with any old scaled near-level
window O(sqrt(n)) gives, uniformly in y,

    w(C_y)<=C sqrt(n)/sqrt(s)=C n q^(-1/4).          (17)

Their union adds at most sqrt(2nq log2) to width, which is bounded by
the same order when q<=sqrt(n). Its antipodal union and all q new
coordinates also cost lower order or the same harmless constant.
Thus physical near-level windows of order n have width
O(n q^(-1/4)) and, by (7c), entropy

    O(n q^(-1/4)sqrt(log q))                         (18)

whenever q tends to infinity. No fixed positive eta limit is claimed:
this is a microscopic physical window of order n, not n^(3/2).

For the physical compiler a new principal child of cap O(q^(3/2))
changes that old window by o(n), since q<=sqrt(n). The bridge itself
has full cap O(n sqrt(q)) with high probability, by applying bounded
differences to sum_i |sum_j B_(i,j)y_j| for each y and taking the same
2^q union. The director's principal-restriction step makes this a
same-order transformation changing only edges incident to q vertices.
The total cap increase is O(n sqrt(q)); every entry remains a sign.
For example q=ceil(log n) rewrites O(n log n) edges, increases the
cap by O(n sqrt(log n)), and gives o(n) microscopic near-code entropy.

The localization track independently audited this physical bookkeeping
and the exact q-range. Selection of the random bridge can depend on
the old child; no common draw working for all old children is asserted.

## 8. Stronger star regularity from the shifted-absolute Stein lemma

The Bernoulli track subsequently proved and froze
[a uniform shifted-absolute comparison](paper_bernoulli_shifted_absolute_2026_09_17.md):

    sup_t |E|q^(-1/2)sum epsilon_i-t|-E|G-t||<=3/q.   (19)

The localization track independently reconstructed the distributional
Stein derivative, exact trapezoid kernel, and one-lattice-atom bound,
and returned PASS. This improves generic Wasserstein replacement for
the particular two-slope test generated by a Boolean coordinate.
Conditional on the other old field coordinates,

    max(A+delta z,B-delta z)
      =(A+B)/2+delta |z-(B-A)/(2delta)|.

Thus each FIXED new-spin pattern has all-offset endpoint error at most
3delta n/q=3sqrt(n/q), and the two conditional-increment endpoints
together cost at most 6sqrt(n/q). One must retain the fixed-pattern
conditioning: the joint maximum over new patterns is not generally a
two-slope function of a single normalized row field.

Take now secondary variance s=1/n, and physical near-level tolerance
T=n/sqrt(q). Then

    E Delta_(y,s)<=[6+kappa/2]sqrt(n/q).

The bounded-difference proxy remains 4q. A deviation Lsqrt(n/q) and
the union over 2^q patterns fail with probability at most

    exp[q log2-L^2 n/(8q^2)].

For q^3<=n and a sufficiently large universal L this is bounded by
exp[-c n/q^2]. The new child cap q^(3/2) is o(T) in this range.
Dividing the old nearcode estimate by sqrt(s) yields

    w(E_parent(n/sqrt(q)))=O(n/sqrt(q)),
    log|E_parent(n/sqrt(q))|
                    =O((n/sqrt(q))sqrt(log q)).       (20)

The pattern-union addition sqrt(2nq log2) and all q new coordinates
are smaller order or harmless constants in this range. Cap cost is
still O(n sqrt(q)), and restriction followed by refill preserves the
original total order. In particular q=Theta(n^(1/3)) gives cap cost
O(n^(7/6)), physical near-level window n^(5/6), width O(n^(5/6)),
and entropy O(n^(5/6)sqrt(log n)).

The canonical actual-full-sign statement and same-order optimization
consequence are maintained by the director in
[Sparse random-star regularization](paper_director_sparse_random_regularization_2026_09_17.md).
The broader bound in Section 7 is retained as a valid independent
fallback, not presented as the final exponent.
