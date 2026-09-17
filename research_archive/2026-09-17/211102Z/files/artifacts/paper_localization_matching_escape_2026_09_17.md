# Matching bridges: an all-energy mixed-tail certificate and its remaining gap

2026-09-17, localization track. This follows the director's request to
test old-word escape, rather than only prove a cheap response on a
selected nearcode. The certificate below concerns an ACTUAL full-sign
parent and includes every old word. Its final section explains exactly
what the existing microscopic preparation estimates do not supply.
It is not a claim that all possible chaining arguments must fail.

## 1. Exact increment geometry for the physical matching law

Let n=p^2, with p even and a real Hadamard kernel U of order p. Let
P_0 be the equal positive/negative fixed-point-free matching law from
[the nonlocal response theorem](paper_localization_nonlocal_sign_response_2026_09_17.md).
For ANY real physical vector w, define

    b_v(a)=sum_u w_(u,v) U_(u,a),
    a_(sigma,{v,z})=b_v(z)+sigma b_z(v),
    B(w)=max_(sigma,e)|a_(sigma,e)|,
    M(w)=max_sigma sum_e a_(sigma,e)^2.

The actual response under a prescribed pole is

    h dot w=sum_(e in matching) epsilon_e a_(sigma,e).

Column Parseval and the elementary square inequality give

    M(w)<=2p||w||_2^2,
    B(w)<=2 max_v ||w_(*,v)||_1.                    (1)

The exact covariance is

    C_0=Cov(P_0)=[p/(p-1)]I-[1/(p-1)]K,
    K=Cov(P_id),

so E|h dot w|<=sqrt(w^T C_0 w). The quantities B and M retain more
information than the upper bounds in (1); for example they vanish on
the identity-law subspace, which P_0 annihilates.

Expanding the matching exponential moment and using
P_(p,k)>=(p/e)^k gives, for B(w)>0 and every real lambda,

 log E_(P_0) exp(lambda h dot w)
 <= e M(w)/(p B(w)^2) [cosh(lambda B(w))-1].          (2)

Indeed cosh(lambda a)-1 is at most
(a^2/B^2)[cosh(lambda B)-1], coefficientwise in its power series.
The bound is uniform in the pole, hence survives the equal mixture.
If B=0 the response is identically zero.

For a symmetric real V and an independent copy V', Jensen and
||V|-|V'||<=|V-V'| show

 E exp(lambda(|V|-E|V|))
 <=E cosh(lambda(|V|-|V'|))
 <=E cosh(lambda(V-V'))=(E exp(lambda V))^2.

Thus for q independent matching columns h_j, (2) yields

 log E exp(lambda sum_j[|h_j dot w|-E|h_j dot w|])
 <=2e q M(w)/(p B(w)^2)[cosh(lambda B(w))-1].        (3)

For positive lambda B<3, use
cosh(s)-1<=s^2/[2(1-s/3)]. The centered sum in (3) is therefore
sub-gamma with variance parameter 2e q M(w)/p and scale B(w)/3.
In particular, for every L>0, with failure probability at most exp(-L),

 sum_j |h_j dot w|
 <=q sqrt(w^T C_0 w)
        +sqrt[4e q M(w)L/p]+[B(w)/3]L.             (4)

One may instead optimize the exact cosh exponent in (3); (4) is a
convenient explicit version. No balanced-column or Gaussian-increment
assumption is made.

## 2. One finite event controls the entire actual parent

Let W be ANY hollow full signing of order n; it need not be the
Hadamard signing defining P_0. Let D be any full signing of order q,
and form the actual parent using q independent columns from P_0.
For every x,y, the cross term has its literal sign entries. Therefore

 Q(parent)<=Q(D)+max_x[|H_W(x)|+sum_j|h_j dot x|].    (5)

Partition the ENTIRE old cube into finitely many nonempty bins C_i.
In each bin choose a finite chain of maps

    pi_(i,0),...,pi_(i,k_i),    pi_(i,k_i)(x)=x.

All maps have Boolean images; they can be arbitrary deterministic
cover maps and need not preserve energy. Let F_(i,0) be the initial
image and E_(i,l) the set of pairs
(pi_(i,l)(x),pi_(i,l-1)(x)), x in C_i. Put

    Z=sum_i(k_i+1),
    L_(i,0)=log|F_(i,0)|+log Z+u,
    L_(i,l)=log|E_(i,l)|+log Z+u.

Define, for w and L,

    R_q(w,L)=sqrt[4e q M(w)L/p]+[B(w)/3]L,
    m(f)=E_(P_0)|h dot f|.

A union bound using (3)--(4), and then telescoping the absolute
responses, proves with probability at least 1-exp(-u), SIMULTANEOUSLY
for every old word x in every bin,

 sum_j |h_j dot x|
 <= q m(pi_(i,0)(x))+R_q(pi_(i,0)(x),L_(i,0))
      +sum_(l=1)^(k_i) {
          q sqrt(w_(i,l)(x)^T C_0 w_(i,l)(x))
                      +R_q(w_(i,l)(x),L_(i,l))},    (6)

where w_(i,l)=pi_(i,l)(x)-pi_(i,l-1)(x). For the initial centers the
actual mean m is retained, rather than replaced by a variance bound.
For increments we use |a|-|b|<=|a-b| and (4).

Substitute (6) into (5) and maximize over ALL x. This is a full-parent
certificate with unchanged actual child D and no discarded energy
band. In particular one can retain the exact deficit
Q(W)-|H_W(x)| against the x-dependent chain cost. Alternatively,
replace each bin's old energy and chain costs by their maxima to
obtain a more elementary finite shell bound.

The same event can be combined with deterministically scheduled
identity-repair columns. Their entire contribution is the explicit
small bipartite block bound in Section 7 of the nonlocal artifact.
It need not be hidden inside a nearcode estimate.

## 3. What Hamming distance alone does and does not control

If x,f differ in m=rho n coordinates, then ||x-f||^2=4rho n, but
B(x-f) is controlled by the MAXIMUM number of disagreements in a
single p-coordinate column, not their average. Without another
structural hypothesis, only B<=4p is available. At q=epsilon n and
a residual union of entropy n s, (4) then has normalized size

    O(epsilon sqrt(rho)+sqrt(epsilon rho s)+s).      (7)

If every column instead has at most a p disagreements, the final
term improves to O(a s). Balanced columns would give a=O(rho), but
this is NOT a consequence of the total Hamming radius.

There is an exact actual-near-ground-code stress test (the positive
ground sector has absolute deficit n). Fix one positive
involution eigenword x_0, and flip the independent pair sign on one
pair {v,z} to obtain another positive eigenword x_1. They differ in
exactly 2p physical coordinates, so their Hamming fraction tends to
zero. For w=x_1-x_0, under P_0 one has EXACTLY

 P(h dot w=+4p)=P(h dot w=-4p)=1/[4(p-1)],
 P(h dot w=0)=1-1/[2(p-1)],
 ||w||^2=8p,
 E exp(lambda h dot w)
      =1+[cosh(4p lambda)-1]/[2(p-1)].              (8)

Indeed the negative pole annihilates w. At the positive pole, only
the designated matching edge has a nonzero coefficient, of magnitude
4p. This proves that there is NO dimension-free subGaussian increment
bound with variance proxy C||x-f||^2 for this physical law, even when
BOTH queries are exact positive eigenwords. The probability in (8)
is polynomial in p, whereas such a bound at 4p would be exponentially
small in p. The large-jump term in (6) cannot simply be dropped.

This does not prove that every chaining proof fails. This particular
pair cloud has only O(p) independent directions, and refined covers
may pay it cheaply. It does prove that a Gaussian metric replacement
based only on Hamming distance would be incorrect.

## 4. Why the existing microscopic star profile does not close escape

Exchangeable-star preparation gives an actual W with a simultaneous
Bernoulli-width bound

    b_Rad(E_T(W))<=T+B_0,      B_0=O(n/sqrt(s)),

where s is the number of resampled vertices and E_T is the absolute
nearcode with UNNORMALIZED energy deficit at most T. Sauer's bound
then gives a useful subexponential entropy estimate while T+B_0=o(n),
and the trivial n log 2 bound outside that range. This is a valid
all-T statement, but not a macroscopic nearlevel entropy improvement.

For an added density q=epsilon n, the natural energy scale of the
desired scalar response term is q sqrt(n)=epsilon n^1.5. At any
fixed epsilon>0 this lies ABOVE n eventually. Consequently the
supplied width/entropy profile is already trivial on the energy
band where a leading epsilon-slope argument must exclude old-word
escape. Combining just that profile with the one-scale consequence
of (4) pays an O(sqrt(epsilon)) selection term on an exponentially
large bin, not an o(epsilon) term.

This is a limitation of the currently supplied profile and certificate,
not a theorem that the actual prepared signing necessarily has such
a large escape contribution. A successful continuation needs another
piece of information: for example a stronger macroscopic profile,
actual small jump parameters along a cover chain, or a comparison
that pays the child deficit and the random bridge jointly. Neither
uniform scalar psi_1 bounds nor a cheap mean on a microscopic code
alone supplies that missing information.

## 5. Centered bridge chaining and the stronger cloned profile

The preceding star-specific limitation is not being extrapolated to
the new cloned-block preparation. That construction supplies a stronger
all-energy profile. We can now test it using centered bridge increments,
which avoid the sum of deterministic mean increments in (6).

For q independent P_0 columns put

    Z_x=sum_(j<=q)[|h_j dot x|-E|h_j dot x|].

For any x,y and w=x-y, let D(h)=|h dot x|-|h dot y|. An independent
copy D' and Jensen give

    E exp(lambda(D-ED))<=E cosh(lambda(D-D'))
       <=E cosh(2lambda h dot w).

Indeed |D-D'|<=|h dot w|+|h' dot w| and
cosh(a+b)<=[cosh(2a)+cosh(2b)]/2 for a,b>=0. Symmetry of P_0
then identifies the last expression with its linear MGF at 2lambda.
Consequently (2) shows that Z_x-Z_y is sub-gamma with variance
parameter4e q M(w)/p and scale2B(w)/3. In particular it has mixed
tails for the TWO metrics

    d_2(x,y)=C sqrt(q)||x-y||_2,
    d_1(x,y)=C max_v||x_(*,v)-y_(*,v)||_1.           (9)

These are centered ABSOLUTE bridge increments, not merely linear
column increments. Their means m(x)=E|h dot x| remain separate.

Dirksen's primary [Theorem 3.5](https://arxiv.org/pdf/1309.3522),
whose statement and full proof were independently read by this track,
therefore bounds the centered supremum on any antipodal Boolean code C
by

    sup_(x in C)|Z_x|
       <=C[sqrt(q) w_G(C)+p log(2|C|)
                         +sqrt(qn u)+p u]           (10)

with probability at least1-C exp(-u), u>=1. To see the two complexity
terms explicitly, Gaussian majorizing measures control gamma_2 for
d_2 by sqrt(q) times the anchored Gaussian spin width. The d_1 diameter
is at most C p, and a partition made discrete once its cardinality
budget exceeds |C| gives gamma_1<=C p log(2|C|). A fixed anchor is
controlled by the same scalar mixed-tail bound; its square norm n
and jump bound O(p) give the last two terms. This imports the classical
Gaussian majorizing-measure theorem, not an unproved equality of widths.

Now take n=p^2 and prepare an arbitrary actual W by the
[cloned-block theorem](paper_director_cloned_block_regularization_2026_09_17.md)
at parameter eta. Write R=floor(sqrt(n)eta^(3/4)) for the clone size,
to distinguish it from the Hadamard kernel order p, and let
B_0=Cn eta^(1/4), ell=log(e/eta). The full physical profile gives

    w_G(E_W(D))<=C sqrt(ell)(D/R+B_0),
    log(2|E_W(D)|)<=C ell(D/R+B_0).

Thus the complexity term in (10), on EVERY energy level D, is bounded by

    C[sqrt(q ell)+p ell](D/R+B_0).                  (11)

For q=epsilon n, the generic jump slope p ell/R is of order
ell/eta^(3/4), which cannot be absorbed by the old deficit as eta->0.
Even if this entire jump term were somehow eliminated, the Gaussian
slope sqrt(q ell)/R could be made small only when

    eta^(3/4)>=C sqrt(epsilon ell).

The certified preparation scale sqrt(eta) would then be at least of
order(epsilon ell)^(1/3), much larger than the desired epsilon scale.
This is an explicit obstruction to THIS combination of the supplied
profile and generic chaining budgets. It is not a lower bound on the
true supremum, nor a theorem against all matching bridge constructions.
Refined jump covers, a more efficient preparation, or a joint mean/energy
comparison may still improve it. Formula (9) identifies exactly which
additional increment geometry would be needed.

Audit note: the Bernoulli researcher independently reconstructed
Sections 1--4 in full and returned PASS, including the rare-eigenpair
law (8), then separately reconstructed all of Section 5 and returned
PASS, including its constants, two metrics, and precise obstruction scope.
The nonlocal response replay now includes180 exact centered-increment
MGF tests, in addition to240 real-vector mixed-tail tests and the two
exact rare-jump laws. All pass.
