# Nonlocal sign laws beyond the Gaussian-angle covariance loss

2026-09-17, localization track. This is an actual full-sign Hadamard
example and a precise covariance-feasibility calculation, NOT a claim
about unknown optimizing children. External novelty is not claimed.

## 1. Exact Boolean eigenlaws with near-unit covariance response

Let p=2^a>=4, n=p^2, and index coordinates by (u,v) in F_2^a squared.
Write chi(u,v)=(-1)^(u dot v),

    F_((u,v),(s,t))=chi(u,t) chi(v,s),   A=F-I.

Thus F^2=p^2 I and A is a hollow full-sign matrix. For a uniform
fixed-point-free involution pi of the p column labels, take independent
fair signs epsilon_e on its pairs. For a prescribed GLOBAL sigma in
{+1,-1}, set z_v=epsilon_e, z_w=sigma epsilon_e on an ordered pair
e={v,w}, and

    h_(u,v)=chi(u,pi(v)) z_v.                         (1)

Then h is Boolean and Fh=sigma p h. Denote its law by P_sigma. Also
let P_id be the positive-eigenvector law

    h_(u,v)=chi(u,v) epsilon_v,

with p independent fair signs. All three laws are centered.

For distinct physical coordinates (u,v),(s,t), their P_sigma covariance
is -F_ab/(p-1) if v=t, and sigma F_ab/(p-1) if v!=t. Indeed, in the
first case pi(v) is uniform outside v and the complete character sum
vanishes; in the second, the product has nonzero sign expectation only
when {v,t} is a matched pair. Under P_id the covariance is F_ab for
v=t and zero for v!=t. Consequently the laws

    nu_- = P_-,
    nu_+ = [(p-1)/(p+1)] P_+ + [2/(p+1)] P_id

have EXACT covariance matrices

    E_(nu_-) hh^T=I-A/(p-1),
    E_(nu_+) hh^T=I+A/(p+1).                         (2)

The unequal pole mixture

    nu = [(p+1)/(2p)] nu_+ + [(p-1)/(2p)] nu_-        (3)

is exactly isotropic. The normalized covariance parameters p/(p+1)
and p/(p-1) both tend to 1, strictly beyond Gaussian sign rounding's
limiting 2/pi. This is a global dependent sign law, not another
coordinatewise rounding map.

For every Boolean query x, Jensen and (2) give the fully uniform bound

 E_nu |h dot x| <= [(p+1)/(2p)] sqrt(n+2 H_A(x)/(p+1))
                  +[(p-1)/(2p)] sqrt(n-2 H_A(x)/(p-1)).       (4)

In particular, on the positive and negative Boolean eigensectors this
is respectively at most

    sqrt(n) sqrt((p+1)/(2p)),
    sqrt(n) sqrt((p-1)/(2p)).                         (5)

Both sectors lie within n of Q(A)=n(p+1)/2. On the entire absolute
eta-near-level set, first p->infinity and then eta->0, (4) yields

    lim_(eta down 0) limsup_p sup_(|H_A(x)|>=Q(A)-eta n^1.5)
             E_nu |h dot x|/sqrt(n) <= 1/sqrt(2).     (6)

At a fixed small eta, the right side can be replaced by
[sqrt(2-2eta)+sqrt(2eta)]/2+o_p(1). This is below the required near-half
extension slope 3/4 once eta is sufficiently small. It does NOT by
itself bound the full parent: new-spin optimization, concentration,
and all lower-energy words still have to be paid.

## 2. Exact support and a concentration limitation

The support of (3) has at most

    S_p=2 (p-1)!! 2^(p/2)+2^p

words, so log S_p=O(p log p)=o(n). This does not imply a dimension-free
linear subGaussian bound. In fact some atom h_0 has probability at
least 1/S_p. If E exp(t h dot v)<=exp(K t^2 ||v||^2/2) for every v,t,
the Chernoff bound at v=h_0 gives

    1/S_p <= P(h dot h_0>=n) <= exp(-n/(2K)),
    K >= n/(2 log S_p)=Omega(p/log p).                (7)

Thus the particular exact laws above cannot be substituted into a
uniform-subGaussian independent-column theorem. Isotropy does give
variance n for every Boolean query, but that alone does not pay a
union over the full cube at positive bridge density.

## 3. The exact radial feasibility window for a general actual child

Suppose an arbitrary hollow full-sign A has two centered sign laws
with covariances I+tau A/sqrt(n) and I-tau A/sqrt(n), tau>=0. Put
P=max_x H_A(x), R=max_x -H_A(x), c=Q(A)/n^1.5, and L=||A||/sqrt(n).
Necessarily

    tau <= 1/L,
    tau <= 2 min(P,R)/[n^1.5(1-1/n)] <= 2c/(1-1/n).   (8)

The first inequality is positive semidefiniteness. The second follows
from the EXACT expectations

    E H_A(h)= +/- tau n(n-1)/(2 sqrt(n)).

More precisely, both covariance targets are in the Boolean correlation
polytope if and only if, for every real hollow symmetric B,

 tau |Tr(A B)|/(2 sqrt(n)) <= min(max_x H_B(x),max_x -H_B(x)). (9)

This is finite-dimensional separation of conv{xx^T:x Boolean}; global
sign symmetrization gives centered laws without changing covariance.
Equation (9) is an exact dual criterion, not a feasible construction
for arbitrary A. Testing B=A recovers (8).

Let f(d)=[sqrt(1+d)+sqrt(1-d)]/2. For an equal pole mixture,
Jensen alone certifies ground response at most f(2c tau) sqrt(n).
Even at the optimistic feasibility ceiling tau=2c, this is below
(3c/2) sqrt(n) only if

    c > 6/sqrt(145) = 0.498272879... .                (10)

Thus the Jensen-based escape in Section 1 does not extend by itself
to a cap constant below 0.493608094.

If one could ALSO prove a uniform conditional Gaussian absolute-
response law, the candidate response would instead be

    kappa f(2c tau),              kappa=sqrt(2/pi).

Combining this desired response with tau<=2c gives the weaker but
still genuine necessary condition

    c > sqrt(72 pi/(256+81 pi^2)) = 0.462940131... .   (11)

For a fixed c in the relevant range the required parameter is

 tau > sqrt(1-[2(3c/(2 kappa))^2-1]^2)/(2c).         (12)

For c=0.493608094 this is approximately 0.7005852, while for c=1/2
it is approximately 0.6414727. Feasibility and such a response law
are separate requirements; neither follows from merely naming an
inverse Ising model. The spectrum also has to satisfy L<1/tau.

## 4. Primary-source high-temperature input: what it does and does not give

Eldan--Koehler--Zeitouni, *A spectral condition for spectral gap:
fast mixing in high-temperature Ising models*,
[arXiv:2007.08200v2](https://arxiv.org/html/2007.08200v2), Theorem 1,
proves a Poincare inequality for

    nu_(J,h)(x) proportional to exp(x^T J x/2+h dot x),
    0<=J<I.

The proof localizes while nearly preserving one chosen observable,
decreases the interaction to rank one, and uses a supermartingale
for the single-site Dirichlet form. Rank-one influence has norm at
most ||J||. For a linear observable, its Dirichlet form is at most
||v||^2, so the theorem implies

    Cov_(J,h) X <= (1-||J||)^(-1) I

uniformly in the external field h. Integrating the log-partition
Hessian along a linear field consequently gives the same centered
linear subGaussian proxy. This would supply useful concentration
IF a desired covariance had already been realized by such a J.

A diagonal shift does not change an Ising law. Accordingly, the
usable sufficient condition for a hollow interaction J_0 is spectral
diameter(J_0)<1, NOT just ||J_0||<1. For J_0=beta A/sqrt(n), a sufficient
condition is beta(lambda_max(A)-lambda_min(A))/sqrt(n)<1. On a
two-sided flat spectrum it only covers beta<1/2 asymptotically.

The theorem does not identify Cov(X) as I+tau A/sqrt(n), does not invert
the moment map, and does not supply the conditional absolute-response
CLT needed above. No claim that its high-temperature regime realizes
the targets in (12) is made here.

There is also a precise failure of the most direct Gibbs proposal in
this primary regime. On the present Hadamard signing,
A^2=(n-1)I-2A. Under the laws proportional to
exp[+/-beta H_A(h)/sqrt(n)], conditional means and
z tanh(beta z)<=beta z^2 imply

 +/- E H_A(h)<= beta E||Ah||^2/(2sqrt(n))
             = beta[n(n-1)-4E H_A(h)]/(2sqrt(n)).

Consequently their signed radial energy coefficients
tau=+/-2sqrt(n) E H_A(h)/[n(n-1)] satisfy respectively

    tau_+<=beta/(1+2beta/sqrt(n)),
    tau_-<=beta/(1-2beta/sqrt(n)).

The negative denominator is positive in the regime considered.
Here spectraldiam(A)/sqrt(n)=2 exactly, so EKZ's sufficient regime
is beta<1/2. Thus both coefficients are at most 1/2+o(1), even below
the Gaussian-angle 2/pi scale. This is a rigorous failure of the
DIRECT proportional interaction J=+/-beta A/sqrt(n) in that regime;
it is not an obstruction to arbitrary inverse-Ising interactions.
Section 11 below strengthens this statement for the Walsh kernel's
exact radial targets by using its signed symmetry group; no analogous
symmetry is assumed for a general actual child.

## 5. Relationship to existing artifacts and status

The involution eigenvectors come from the discrepancy track and the
[all-law Hadamard obstruction](paper_localization_all_law_hadamard_obstruction_2026_09_17.md).
The exact covariance correction and its positive scalar-response use
are the present deduction. The archive already contains an explicit
Paley-conference analogue with eigenlaws of covariance I+/-C/r in
[the September 6 cavity/eigenlaw artifact](transfer_fresh_cavity_and_isotropic_ground_law_2026_09_06.md),
Sections 4--6; it also warns that these laws are not uniformly
subGaussian. Thus no claim of a new general correlation-polytope
inclusion or a new basic isotropic-ground-law construction is intended.

This theorem escapes the scalar-affine Gaussian-sign barrier proved in
[the companion response artifact](paper_localization_gaussian_sign_response_2026_09_17.md)
on a concrete actual full-sign family. It neither improves the original
asymptotic cap frontier nor shows that exact minimizers have the needed
covariance geometry. Section 7 pays the finite independent-parent and
all-state costs explicitly; it does not show that they are small enough
for an unconditional asymptotic recurrence.

## 6. A uniform matching comparison recovers the Gaussian absolute factor

This section is an elementary NEW strengthening of (6); it does not
assume that the overlap itself satisfies a Gaussian CLT.

### 6.1 A signed-perfect-matching lemma

Let p>=4 be even. Give every unordered edge of K_p a real coefficient
c_e, with max_v sum_(e incident v) c_e^2<=D, where D is fixed. Choose a
uniform perfect matching M and independent fair edge signs, and put

    T=sum_(e in M) epsilon_e c_e,
    v=E T^2=sum_e c_e^2/(p-1).

Then, uniformly in the entire coefficient array,

    E|T| <= kappa sqrt(v)+r_p(D),    r_p(D)->0.        (13)

One can take r_p(D)=O_D(1/sqrt(log p)). Here is a complete quantitative
proof. For fixed t let w_e=1-cos(t c_e)>=0, d=Dt^2/2, and W=sum_e w_e.
Then max_v sum_(e incident v)w_e<=d and W<=pd/2. Expanding over
submatchings gives

 phi_T(t)=sum_(k=0)^(p/2) (-1)^k A_k/P_(p,k),
 A_k=sum_(k-edge matchings E) product_(e in E) w_e,
 P_(p,k)=(p-1)(p-3)...(p-2k+1).                       (14)

Compare this with psi(t)=exp[-W/(p-1)]. The total weight of ordered
k-tuples containing an intersecting pair is at most

    binom(k,2) [sum_v (sum_(e incident v) w_e)^2] W^(k-2).

Since the bracket is at most pd^2,

 0<=W^k/k!-A_k <= [2k^2/p] (pd/2)^k/k!.             (15)

For every k<=p/2,

    P_(p,k)>=(p/e)^k.                               (16)

Indeed, the geometric mean of a descending prefix is at least the
geometric mean of all p/2 odd factors, and concavity of log gives
sum_(j=1)^(p/2)log(2j-1)>=integral_0^(p/2)log(2s)ds
=(p/2)(log p-1). For k<=p/4 also

 0<= (p-1)^k/P_(p,k)-1 <=(3k^2/p) exp(3k^2/p).

Combine this with (15). For k>p/4 use (16) and 1<=4k/p to bound the
tail of both series; for the exponential tail k>p/2 use 1<=2k/p.
Summing k a^k/k! and k^2 a^k/k! gives the convenient finite bound

    |phi_T(t)-psi(t)| <= [30/p](d+d^2) exp(2d).       (17)

There is no assumption that all cosines are positive. The signs in
(14) are retained, while (15)--(17) use absolute coefficient errors.
Since 1-cos z<=z^2/2,

    psi(t)>=exp(-v t^2/2).

The identity E|T|=(2/pi) integral_0^infinity[1-phi_T(t)]dt/t^2 and
|phi_T|<=1 now give, for any R>0,

 r_p(D) <= (60/(pi p)) exp(DR^2)
                    [DR/2+D^2 R^3/12]+4/(pi R).    (18)

Taking R=sqrt(log p/(2D)) proves the stated rate for D>0; D=0 is
trivial. The comparison law with characteristic function psi is a
symmetric compound-Poisson law. Its Gaussian absolute upper bound
does not require its jump part to disappear.

### 6.2 Application to exact and near Hadamard eigensectors

For a Boolean x write its column Fourier matrix as

    b_v(a)=sum_u x_(u,v) chi(u,a),    C_(v,a)=b_v(a)/p.

Every row of C has squared norm 1, and ||C||_F^2=p. If Fx=tau p x,
then C=tau C^T. Under P_sigma from Section 1,

    h dot x/p=sum_(e={v,w} in M) epsilon_e
                              [C_(v,w)+sigma C_(w,v)].      (19)

The opposite pole is identically zero. At the matching pole the
coefficients have weighted row degree at most 4, and their variance
is at most 2p/(p-1). Thus (13), followed by the weights in (3), proves

 sup_(Fx=+/-p x) E_nu|h dot x|/p
       <= kappa/sqrt(2)+r_p(4)/2+1/sqrt(p).           (20)

The identity-law contribution here is at most 1/sqrt(p), by Parseval
and its weight 1/p in (3). This improves Jensen's 1/sqrt(2) to
kappa/sqrt(2)=1/sqrt(pi), approximately 0.564189584.

For completeness the ENTIRE absolute near-level set also obeys a
uniform bound, not merely the displayed eigenword family. If
|H_A(x)|>=Q(A)-eta p^3, with 0<eta<1/4, choose the appropriate pole tau.
Then u=tau x^T F x/p^3>=1-2eta, and

    ||C-tau C^T||_F^2=2p(1-u)<=4eta p.

Put S=(C+tau C^T)/2 and call a row bad if its squared norm in
C-tau C^T exceeds 1. There are at most 4eta p bad rows. On good rows
||S_(v,*)||<=3/2, while

    sum_(v bad)||S_(v,*)||^2<=10eta p.

Delete from (19), at the matching pole, every edge touching a bad row.
The retained weighted degree is at most 9. The deleted sum has variance
at most 40eta p/(p-1). The entire matching-pole variance is at most
2p/(p-1), since S is an orthogonal Frobenius projection of C. At the
opposite pole the variance is at most 2eta p/(p-1). Applying (13) to
the retained sum and Cauchy--Schwarz to the deleted and opposite sums
therefore gives the explicit full-code estimate

 sup_(|H_A(x)|>=Q(A)-eta n^1.5) E_nu|h dot x|/sqrt(n)
   <= kappa/sqrt(2)+4sqrt(eta)+r_p(9)/2+1/sqrt(p).    (21)

This proves the Gaussian-factor response discount for actual nonlocal
sign laws at near-unit covariance, without returning to the 2/pi
Gaussian-angle loss. It remains specifically a Hadamard-family result.

The [director's Gaussian--Poisson comparison](paper_director_matching_gaussian_poisson_2026_09_17.md),
including the localization track's weighted Palm refinement in Section 7,
improves r_p(D) to O_D(p^(-1/4)). Thus the error term in (21), apart
from 4sqrt(eta), can be taken O(n^(-1/8)). The independent elementary
characteristic-function proof above is retained as a separate audit.

## 7. Uniform exponential moments and a paid finite parent certificate

There is a useful distinction between the fixed-point-free laws and
the small identity repair in (3).

For EVERY Boolean x, not just a near-eigenword, the coefficients in
(19) satisfy |c_e|<=2 and sum_e c_e^2<=2p. For any real t set
w_e=cosh(t c_e)-1>=0. The expansion (14), now with only positive terms,
(16), and the power series for cosh give

 E_(P_sigma) exp(t h dot x/p)
  <=exp[(e/p)sum_e w_e]
  <=exp[(e/2)(cosh(2t)-1)].                          (22)

This is uniform in p, sigma, and ALL Boolean x. It is not a uniform
subGaussian estimate for arbitrarily large t; it is a uniform
exponential-moment estimate at each fixed t. In particular the centered
absolute response satisfies a universal Bernstein bound: for q
independent columns with arbitrary prescribed poles,

 P(sum_j (|h_j dot x|-E|h_j dot x|)
          >C sqrt(n)[sqrt(q L)+L]) <= exp(-L),       (23)

where C is universal. To derive it directly, (22) at t=+/-1 bounds
E exp(|h dot x|/p) by a universal constant. Expanding powers gives
E|h dot x/p|^k<=C_0 k!, and centering gives the usual power-series
bound log E exp(t(|h dot x|/p-E|h dot x|/p))<=C_1 t^2 for |t|<=c_1.
Chernoff optimization proves (23). No imported concentration theorem
with unverified hypotheses is needed.

The equal mixture P_0=(P_++P_-)/2 also satisfies (22), and

    Cov(P_0)= [p/(p-1)]I-[1/(p-1)]K,
    K=Cov(P_id).

Here K is block diagonal of rank p, with its p nonzero eigenvalues
equal to p. Thus Cov(P_0)<=p/(p-1) I; its failure of isotropy is an
explicit rank-p kernel. The response conclusion (21) also holds for
P_0 with a harmless o_p(1) change and without the identity term.

An exact covariance repair in the aggregate need not impose the bad
exponential moments of the rare identity component. If q is divisible
by p, schedule q/p independent identity-law columns and q(p-1)/p
independent P_0 columns. The sum of their covariance matrices is
EXACTLY q I. For arbitrary q, balanced integer counts leave an O(p)
operator-norm covariance discrepancy, negligible relative to q when
q is a fixed positive fraction of n. This is a deterministic label
schedule, not independent draws from the mixture (3).

All r identity columns together have bridge cap exactly

    p max_(s in {+-1}^p,y in {+-1}^r) s^T E y,

where E is a p by r independent sign matrix. This follows by writing
their x-response in terms of b_v(v), whose independent feasible range
includes both endpoints +/-p in each column. A union bound over the
2^(p+r) pairs consequently gives, with probability at least 1-exp(-L),

    identity bridge cap <=p sqrt(2pr[(p+r)log 2+L]). (24)

For r=O(q/p), q=Theta(n), and L=O(p), this is O(n^1.25)=o(n^1.5),
uniformly over every old word and every new-spin choice.

For the remaining q_0 P_0 columns, partition ALL old words into J
nonempty bins C_j. Suppose |H_A(x)|<=b_j n^1.5 on C_j and
E_(P_0)|h dot x|<=mu_j sqrt(n) there, and let M_j=|C_j|. Then, for
any actual full-sign child D, the actual parent obeys simultaneously
with probability at least 1-2exp(-u),

 Q(parent) <= Q(D)+identity bridge bound (24)
   +max_j {b_j n^1.5+q_0 mu_j sqrt(n)
       +C sqrt(n)[sqrt(q_0 L_j)+L_j]},
    L_j=log M_j+log J+u.                            (25)

This pays the entire new-spin maximization and every old energy level.
One may use (21) on near-level bins and the uniform covariance bound
mu_j<=sqrt(p/(p-1)) on the rest. No small entropy of the full near-level
code is asserted. In particular (25) is a genuine usable finite
certificate, but not an unconditional asymptotic minimizing recurrence.

## 8. Verification record

The discrepancy researcher independently reconstructed Sections 1--3
and returned PASS. The director and Bernoulli researcher independently
reconstructed the full matching expansion, nearlevel row deletion,
uniform exponential moments, deterministic covariance repair, and
all-bin parent certificate, and returned PASS. The finite eta constant
4 specifically uses the exact mixture weights (p-1)/(2p), not a
replacement of those weights by 1/2 before cancellation.

`computations/paper_localization_2026_09_17_nonlocal_sign_response.py`
passes exact covariance, isotropy, and response checks on all 32,768
projective queries at p=4 and 1,550 probes at p=8. It also passes 2,080
finite matching-characteristic comparisons, 1,600 exponential-moment
checks, and 20 exact identity-column/full-bridge reductions. These
finite calculations audit identities and inequalities; they are not
used to infer the asymptotic theorem.
Output: `tmp/paper_portfolio_2026_09_17/localization/nonlocal_sign_response_audit.json`.

## 9. Portability to every real Hadamard kernel, and a sharp scope barrier

Director extension, independently reconstructed by the localization
track. The character group in Section 1 is not essential. Let U be
ANY p by p real Hadamard matrix, so its entries are signs and
UU^T=U^T U=pI. On ordered pairs define

    F_((u,v),(s,t))=U_(u,t) U_(s,v),    A=F-I.

The matrix F is symmetric, has diagonal one, and F^2=p^2 I by the
two orthogonality identities. Define the matching and identity laws by

    h_(u,v)=U_(u,pi(v)) z_v,
    h_(u,v)=U_(u,v) epsilon_v,

respectively. Direct summation gives Fh=sigma p h in the first case
and Fh=p h in the second. The covariance proof of Section 1 uses only
sum_a U_(u,a)U_(s,a)=p 1_(u=s); after excluding a=v, the same-column
off-diagonal correlation is -U_(u,v)U_(s,v)/(p-1). The cross-column
calculation is unchanged. Thus EVERY exact covariance, mixture weight,
nearlevel response, exponential-moment bound, and scheduled-repair
statement above holds without modification.

For an arbitrary Boolean query, the relevant coefficient matrix is

    C_(v,a)=(1/p)sum_u x_(u,v) U_(u,a).

Its rows have squared norm one. Transforming Fx by U in its first
coordinate transposes C, so Fx=sigma p x is equivalent to C=sigma C^T.
These are precisely the ingredients of the signed-matching response
proof. No arithmetic structure, symmetry of U, or character identity
has been used. This portability does not assert Hadamard existence at
every order or infer any prime-gap padding result.

There is also an important exact limitation of this tensor construction.
Even if U is ANY p by p sign matrix, with no orthogonality assumption,
the same formula for F is symmetric with diagonal one, hence A=F-I
is still a valid full signing on n=p^2 vertices. The Boolean query
X=U, viewed as a length-p^2 word, satisfies

    x^T F x=Tr[(U^T U)^2]>=p^3,
    H_A(x)>=(p^3-p^2)/2.                           (26)

The inequality is Cauchy--Schwarz on the p nonnegative eigenvalues of
U^T U, whose sum is p^2. Therefore the ENTIRE transpose-tensor family,
even after dropping orthogonality, obeys

    Q(A)/n^1.5 >= 1/2-1/(2p).

It cannot furnish a family with limiting cap below 1/2 by merely
optimizing the sign kernel U. This is a Boolean witness for the actual
signing, not a weakness of a spectral upper bound. The nonlocal bridge
primitive may be transferable elsewhere, but this elementary kernel
generalization does not reach the known strict-below-half frontier.

The same replay additionally checks 24 nonsymmetric switched/permuted
Hadamard kernels and 80 arbitrary sign-kernel Boolean witnesses.

## 10. The constant 1/sqrt(pi) is sharp for this specific law

There are explicit flat-Fourier Boolean eigenqueries attaining the
Gaussian-factor upper asymptotically. This is sharpness of the present
matching law, NOT a lower bound for every possible nonlocal sign law.

Take p=2^(2a)=r^2 and the Walsh kernel U. Index its rows by pairs
(a,b) in F_2^a squared and set d_(a,b)=(-1)^(a dot b). Directly summing
over a gives

    sum_(a,b) d_(a,b)(-1)^(alpha dot a+beta dot b)
                  =r(-1)^(alpha dot beta).

Thus every Walsh coefficient of d has magnitude sqrt(p). The Boolean
query matrix X=diag(d)U is a +p eigenword for the transpose kernel:

    F(X)=U X^T U=p X.

Its normalized column transform C=U diag(d)U/p is symmetric and has
every entry of magnitude 1/sqrt(p), with all diagonal entries equal
to 1/sqrt(p). Therefore under P_+ the normalized overlap is EXACTLY

    h dot x/p =_law (2/sqrt(p)) S_(p/2),

where S_j is the sum of j independent fair signs. Under P_- it is
identically zero; under P_id it is S_p/sqrt(p). Writing mu_j=E|S_j|,
the corrected isotropic law (3) consequently has the exact response

    E_nu|h dot x|/p
       =[(p-1)mu_(p/2)+mu_p]/[p sqrt(p)].            (27)

The elementary central-binomial formula for mu_j, or the already
reconstructed shifted-absolute comparison, shows that (27) tends to
kappa/sqrt(2)=1/sqrt(pi). This positive eigenword is exactly n below
the absolute cap and hence belongs to every fixed eta-nearlevel
eventually. Together with (21), this identifies the iterated
nearlevel response limit for this law along p=2^(2a):

 lim_(eta down 0) limsup_p sup_(absolute eta-nearlevel)
                     E_nu|h dot x|/sqrt(n)=1/sqrt(pi).

The same limit holds for the uncorrected equal matching mixture P_0.
Neither more accurate Gaussian approximation nor removal of the
matching error term can improve this particular scalar-response
constant; a further improvement would require changing the law.

## 11. Exact Walsh symmetry closes the radial inverse-Ising loophole

This is a further scoped consequence of the primary high-temperature
mapping in Section 4. On the WALSH transpose kernel, arbitrary finite
Ising parameters cannot realize a centered radial covariance target
in some hidden high-temperature direction different from A.

Write the coordinate set as V=F_2^(2a), with nondegenerate alternating
bilinear form [z,w], so F_(z,w)=(-1)^[z,w] and n=p^2. There are two
types of signed permutation symmetries of A:

    (P_b x)_z=(-1)^[z,b] x_(z+b),
    (P_T x)_z=x_(Tz), where [Tz,Tw]=[z,w].

Direct expansion of [z+b,w+b] proves P_b^T A P_b=A. The identity for
P_T is immediate. These symmetries are vertex-transitive and transitively
map any distinct ordered pair to (0,d), d!=0, with the accompanying
sign exactly accounting for A_(z,w). The symplectic linear group is
transitive on nonzero d: extend each nonzero vector and one vector
pairing to 1 to a symplectic basis, then map the two bases.

Consider ANY finite-parameter Ising law

    nu(x) proportional to exp(x^T J x/2+h dot x),

using the unique hollow representative for J. Suppose E X=0 and
E XX^T=I+alpha A. Its parameters are then necessarily

    h=0,                   J=gamma A.              (28)

Here is a short proof that does not assume parameter identifiability.
Two full-support Ising laws with the same first and second moments
have zero symmetrized relative entropy, by subtracting their log
densities and taking expectations. They are equal. Their parameters
are equal because the Boolean Fourier functions x_i and x_i x_j are
linearly independent modulo constants. Global reversal preserves the
prescribed moments and changes h to -h, so h=0. Every signed symmetry
P above also preserves the target moments, so uniqueness gives
P^T J P=J. Translation sends a pair (z,w) to (0,w+z), with multiplier
A_(z,w); symplectic transitivity makes J_(0,d) constant over d!=0.
Consequently J_(z,w)=gamma A_(z,w), proving (28).

Even a NONUNIFORM diagonal shift cannot improve the EKZ regime in
this symmetric situation. If 0<=J+D< I for any diagonal D, conjugate
by the vertex-transitive signed symmetry group and average. Convexity
of the two matrix inequalities gives

    0<=J+d I<I

for a scalar d. Thus spectraldiam(J)<1 is necessary for any such
diagonal-shift realization, as well as sufficient. In (28), write
gamma=+/-beta/p, beta>=0. The spectral diameter is exactly 2beta,
so the primary theorem's regime requires beta<1/2.

If the covariance is I+/-tau A/p, its radial coefficient has the
same sign as gamma, since the Gibbs expectation of H_A is strictly
increasing in gamma and vanishes at zero. The conditional-mean
inequalities from Section 4 now imply the exact bounds

    tau_+ < p/[2(p+1)],
    tau_- < p/[2(p-1)].                             (29)

In particular every centered radial target realizable by ANY finite
Ising interaction inside this spectral high-temperature class has
tau<=1/2+o(1). The explicit nonlocal targets (2) have exactly twice
these limiting finite coefficients and live outside this class.

This does not say that no Ising model realizes an interior radial
target at larger beta, nor that the EKZ sufficient condition is a
necessary condition for concentration. It rules out a specific
attempt to combine an arbitrary inverse-Ising realization of the
radial target with that primary spectral high-temperature theorem.
It remains a Walsh-family symmetry theorem, not a statement about
unknown asymmetric optimizing children.
