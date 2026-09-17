# A quadratic dual excludes every isotropic low-response center cover

2026-09-17. New deduction of the localization track, combining the
discrepancy track's random-involution Boolean eigenvectors with a
quadratically corrected absolute-response dual. The director independently
reconstructed the proof and returned PASS. External novelty is not claimed.

This concerns an explicit actual full-sign family with normalized cap
tending to 1/2. It is NOT a family of asserted asymptotic minimizers and
does not disprove the conditional center-cover extension theorem.

## 1. Exact statement, with no mode, support, or aspect-ratio restriction

Let p=2^a>=4 and n=p^2. Index coordinates by (u,v) in F_2^a x F_2^a,
and define the alternating Walsh matrix and its hollowing by

    F_((s,t),(u,v))=(-1)^(s dot v+t dot u),
    A=F-I_n.

Then A has zero diagonal and every off-diagonal entry is a sign. Put

    C_p=2+p/(p-3),       b_p=1/sqrt(6 C_p).

Let E_p be the union of the Boolean +p and -p eigenvectors of F. It
suffices below to use only the random-involution subfamily defined in
Section 2. For EVERY isotropic probability law nu on {+-1}^n,

    max_(x in E_p) E_nu |h dot x| >= b_p sqrt(n).      (1)

More strongly, there is one explicit probability law pi_p supported
on E_p for which

    E_(x~pi_p) E_(h~nu)|h dot x| >= b_p sqrt(n)        (2)

for every such nu. Thus nu may be selected after seeing A, all ground
words, a proposed center code, or an arbitrary coordinate grouping.
Its support and computational description are unrestricted.

Let B be any nonempty Boolean center family with

    sup_(f in B) E_nu|h dot f| <= mu sqrt(n).

Then

    max_(x in E_p) min_(f in B) d_H(x,f)/n
                     >= [(b_p-mu)_+]^2/4.            (3)

In particular, if mu tends to zero as p tends to infinity, this radius
has liminf at least 1/72. The center family need not be subexponential.

## 2. Boolean eigenvectors from both global eigensectors

Choose a uniformly random fixed-point-free involution pi on the p
labels. Choose independently a fair GLOBAL sign sigma, and a fair
sign epsilon_e on every two-cycle e of pi. Order each pair (v,w) once,
and set z_v=epsilon_e, z_w=sigma epsilon_e. Define

    X_(u,v)=(-1)^(u dot pi(v)) z_v.                   (4)

Then z_(pi(v))=sigma z_v, for both endpoints of every pair. Summing
first over u gives

    (F X)_(s,t)=p (-1)^(s dot pi(t)) z_(pi(t))
               =sigma p X_(s,t).                    (5)

Thus every output is Boolean and is a sigma p eigenvector. The sign
sigma is global, not independently chosen per pair. This dependence
is necessary to preserve the eigenvector identity.

The positive eigensector is the previously reconstructed involution
law; including its negative counterpart is the only change. See
[the discrepancy track, Section 17](paper_discrepancy_2026_09_17.md)
for the original construction and exact coordinate-moment audit.

## 3. A pointwise quadratic lower minorant for the ground response

Fix ANY physical sign word h, and take its Fourier transform separately
inside each v-column:

    b_v(a)=sum_u h_(u,v) (-1)^(u dot a),
    D(h)=sum_v b_v(v)^2,
    W(h)=sum_(v!=a) b_v(a)^2=p^3-D(h).                (6)

The final identity is columnwise Parseval, since h has p^2 sign
coordinates. Every |b_v(a)| is at most p.

For a fixed matching pi put

    M=sum_v b_v(pi(v))^2
     =sum_(e={v,w} in pi) alpha_e,
    alpha_{v,w}=b_v(w)^2+b_w(v)^2 <= 2p^2.           (7)

Conditional on pi and sigma, the overlap in (4) is a Rademacher sum
with one independent driver per pair:

    h dot X=sum_(e={v,w}) epsilon_e
                         [b_v(w)+sigma b_w(v)].

Let its conditional variance be S_sigma. A Rademacher sum T obeys
E T^4<=3(E T^2)^2, and Holder therefore gives E|T|>=sqrt(E T^2/3).
Since S_++S_-=2M,

    E_(sigma,epsilon)|h dot X|
       >= [sqrt(S_+)+sqrt(S_-)]/(2sqrt(3))
       >= sqrt(M/6).                                (8)

This uses only an elementary fourth-moment expansion; no sharp
Khintchine constant is imported.

An edge belongs to a uniform perfect matching with probability
1/(p-1), and two disjoint edges both belong with probability
1/[(p-1)(p-3)]. Therefore, writing m=E_pi M,

    m=W/(p-1),
    E_pi M^2 <= 2p^2 m + W^2/[(p-1)(p-3)]
             =2p^2 m+[(p-1)/(p-3)]m^2
             <=p^2 C_p m.                           (9)

For the last step use W<=p^3, so m<=p^3/(p-1). The diagonal contribution
uses alpha_e^2<=2p^2 alpha_e. Dropping the restriction to disjoint
pairs only increases the nonnegative off-diagonal sum, explaining the
second-moment upper bound even when some Fourier coefficients vanish.

Holder interpolation gives E sqrt(M)>=(E M)^(3/2)/(E M^2)^(1/2).
If m=0 the desired statement is trivial; otherwise combine this with
(8)--(9). The resulting pointwise lower minorant is

    f_p(h):=E_(X~pi_p)|h dot X|
        >= [p^3-D(h)]/[p(p-1)sqrt(6 C_p)].            (10)

The correction D(h) is genuinely quadratic in h. It is essential:
for h_(u,v)=(-1)^(u dot v), every off-diagonal Fourier coefficient is
zero, and this entire fixed-point-free law has overlap zero. Thus a
positive pointwise constant lower bound for f_p(h) would be false.

## 4. Isotropy removes the quadratic correction exactly

If E_nu hh^T=I_n, then E_nu b_v(v)^2=p for every v, so

    E_nu D(h)=p^2.                                   (11)

Averaging (10) gives

    E_nu f_p(h) >= (p^3-p^2)/[p(p-1)sqrt(6 C_p)]
                 =p/sqrt(6 C_p).

This proves (2), and hence (1). It is an explicit feasible quadratic
minorant of the absolute-overlap function, exactly the type of dual
certificate in the finite isotropic-response linear programs, but
uniform in p and proved without solving those programs.

To prove (3), suppose every x in the support has a center f with
d_H(x,f)<=rn. Isotropy and Cauchy--Schwarz give

    E_nu|h dot x|
       <=E_nu|h dot f|+sqrt(E_nu[h dot (x-f)]^2)
       <=mu sqrt(n)+||x-f|| <=(mu+2sqrt(r))sqrt(n).

Average over pi_p and compare with (2). Since C_p tends to 3, the
limiting squared-radius constant is 1/(4*6*3)=1/72.

## 5. Application to the complete near-level code of actual full signs

The matrix F satisfies F^2=p^2 I and has diagonal one. Hence

    Q(A)<=n(p+1)/2.

Every sigma=-1 word from (4) attains the negative endpoint, so equality
holds. A sigma=+1 word has absolute energy n(p-1)/2, exactly n below
the cap. It follows that the whole support of pi_p is contained in

    {x: |H_A(x)|>=Q(A)-eta n^(3/2)}

whenever eta>=1/p. For every fixed positive eta, that condition holds
eventually. Thus (3) applies to the FULL near-level set at every fixed
eta and forces a positive limiting radius for every vanishing-response
center family under every isotropic physical sign law.

In particular this defeats the entire low-Phi center family for any
isotropic tensor-mode law, after any switching, permutation, grouping,
or choice of k, because its physical law h tensor g is isotropic and
its scalar absolute response is at most sqrt(n) Phi. Leftovers can
be filled by independent signs. There is no k=o(sqrt(n)) limitation.

The previous frame-specific stress test gave the stronger radius 1/2
in its restricted growing-mode range. The present result trades that
constant for removal of ALL mode-law restrictions. Its scope remains
the explicit near-half family, whose asymptotic normalized cap is 1/2,
not actual unknown minimizers at the sharper current frontier.

The discrepancy track subsequently proved
[extension stability](paper_discrepancy_extension_obstruction_2026_09_17.md):
the isotropic radius obstruction survives every principal extension with
o(n) new vertices and o(n^(3/2)) cap increase. That proof was independently
audited here. It applies in particular to the small-cost external-field
regularizers, and separates their shrinking-window width conclusion from
the stronger low-response geometry that remains unavailable on this family.

## 6. Audit and replay

### 6.1 Small width and entropy do not imply cheap isotropic response

The discrepancy-track researcher extracted the following corollary.
The support of pi_p has at most

    2 (p-1)!! 2^(p/2)

words, so its logarithmic cardinality is O(p log p)=O(sqrt(n) log n).
The Gaussian maximum bound therefore gives its width at most
sqrt(2n log|supp(pi_p)|)=O(n^(3/4)sqrt(log n))=o(n).
Nevertheless (2) shows that EVERY isotropic physical sign law has
response at least (1/sqrt(18)-o(1))sqrt(n) on some word of this code.

Thus subexponential code entropy and sublinear Gaussian width, even
together, do not imply the low-isotropic-response premise of the positive
center construction. In particular the external-field regularization
theorem cannot be combined with that construction merely by identifying
these two notions of complexity. This is an explicit finite-support
counterexample to that proposed inference, not a failure of either
regularization or the conditional center-cover theorem.

### 6.2 Independent proof audits and integer replay

The director independently reconstructed (8)--(11) and (3), including
the disjoint-matching second moment and the global eigensign factor,
and returned PASS before the canonical proof was written. The
Bernoulli-track researcher subsequently read the entire canonical proof
and independently returned PASS on the moment, interpolation, covariance,
radius, and exact-cap steps.
The discrepancy-track researcher also independently reconstructed those
steps and returned PASS; the small-width corollary above is theirs.

The replay is
`computations/paper_localization_2026_09_17_all_law_dual.py`.
At p=4 it exhausts all 32768 projective physical sign queries against
the complete 24-atom two-sector law; at p=8 it tests 722 queries against
the complete 3360-atom law. It checks the squared dual inequality by
integer arithmetic, exact matching moments, both eigensectors, cap
deficits, isotropic cancellation, and the exceptional chirp. All pass.
Output is preserved in
`tmp/paper_portfolio_2026_09_17/localization/all_law_dual_audit.json`.

## 7. Stronger pointwise dual: no sign-law restriction at all

The discrepancy-track researcher supplied the following strengthening,
independently reconstructed and audited PASS by the localization track.
It removes isotropy from the RESPONSE obstruction (not from the constant
Hamming-radius conclusion).

Alongside the fixed-point-free two-sector law pi_p, consider the identity
involution and independent signs epsilon_v:

    X^id_(u,v)=(-1)^(u dot v) epsilon_v.

Every such word is a +p eigenvector of F and has deficit n from Q(A).
For the same D(h) as in (6), its mean absolute overlap satisfies

    f_id(h)=E|sum_v epsilon_v b_v(v)|>=sqrt(D(h)/3).

Write

    a_p=p^2/[(p-1)sqrt(6C_p)],
    d_p=p^(3/2)/sqrt(3),
    omega_p=a_p/(a_p+d_p),
    c_p=a_p d_p/(a_p+d_p)=(1/sqrt(18)+o(1))p.

Mix pi_p with the identity law, putting weight omega_p on the latter.
For z=D(h)/p^3 in [0,1], (10) and the preceding bound give

    E_mixed |h dot X|
       >=(1-omega_p)a_p(1-z)+omega_p d_p sqrt(z)
        =c_p[(1-z)+sqrt(z)]>=c_p                     (14)

for EVERY physical sign word h. The last inequality uses sqrt(z)>=z.
This explicitly repairs the exceptional chirps rather than relying on
isotropic averaging. Consequently (14) also holds after averaging h
under ANY probability law whatsoever on physical sign words.

The support still has logarithmic size O(p log p), so the small-width
corollary in Section 6.1 is stronger than first stated: this explicit
low-width code cannot be uniformly cheap under ANY sign-column law.
To convert response into Hamming distance one still needs covariance
or another response-increment estimate; no covariance-free 1/72 radius
is asserted.

### 7.1 An actual arbitrary-bridge lower bound

Let P=[A,B;B^T,D] be any full-sign principal extension, with q new
vertices and arbitrary sign bridge columns B_j. For each old word x,
orient y_j by sign(H_A(x))*sign(B_j dot x). Then

    Q(P)>=|H_A(x)|+sum_j|B_j dot x|-Q(D).

Every word in the mixed law has |H_A(x)|>=Q(A)-n. Average the last
inequality and apply (14) separately to every actual sign column:

    Q(P)>=Q(A)-n+q c_p-Q(D).                         (15)

No independence, randomness, isotropy, or mode representation of B
is needed. If Q(D)=O(q^(3/2)), this gives a positive leading linear
bridge cost for small fixed q/n. For an unrestricted high-cap new
child, its displayed Q(D) term must be retained; the result does not
claim a universal positive increment after optimizing arbitrary D.

The existing replay was extended to include the full identity-involution
laws at p=4,8, their exact eigenvector identities, the conditional fourth-
moment bound, and the pointwise mixed response for every previously tested
sign query. All passed. The mixture weights involve square roots and are
evaluated as floating diagnostics; the analytic inequality (14) is the
uniform proof.
