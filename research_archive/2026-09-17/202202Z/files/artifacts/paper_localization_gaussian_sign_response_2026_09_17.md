# Actual-child Gaussian-sign response: isotropic variance mixtures and a slope barrier

2026-09-17. Localization-track reconstruction and extension. The affine
Gaussian proposal is due to the director. The new mechanism below is
latent variance mixing: covariance gains cancel, but absolute-response
gains do not. All columns are genuine sign vectors. This is not same-map
Krivine rounding, and no spectral bound for original minimizers is assumed.

## 1. Exact rounding facts and their primary source

Let A be an n by n hollow symmetric full signing, H_A(x)=x^T A x/2,
and Q(A)=max_x |H_A(x)|. Set kappa=sqrt(2/pi). For any admissible s,

    R_s=I+s A/sqrt(n)>=0,
    h_s=sign(G_s),  G_s~N(0,R_s),
    rho_n(s)=(2/pi)arcsin(s/sqrt(n)).

Coordinate variances equal one, so signs are unambiguous almost surely,
even if the full Gaussian covariance is singular. The exact identity is

    E h_s=0,    E h_s h_s^T=I+rho_n(s) A.             (1)

The real Gaussian sign identity and its planar-angle proof were read
directly in Friedland--Lim--Zhang, *An elementary and unified proof of
Grothendieck's inequality*, Lemma 2.1, pp. 330--331:
https://ems.press/content/serial-article-files/44373?nt=1
For unit Gaussian directions making angle theta, the two disagreement
wedges have total probability theta/pi; this gives correlation
1-2theta/pi=(2/pi)arcsin(r). Since A_ij is exactly +-1, oddness of arcsin
makes (1) exact, rather than a first-order covariance approximation.

If R_s<=K I, then h_s has the uniform linear MGF bound

    E exp(v dot h_s)<=exp(K ||v||^2/2).               (2)

This follows from the Gaussian Holder inequality of Chen--Dafnis--Paouris,
Theorem 1(i), with scalar blocks and all exponents K:
https://arxiv.org/html/1306.2410v2
It bounds the left side by product_i cosh(K v_i)^(1/K); use
cosh(u)<=exp(u^2/2). Their positive-semidefinite covariance hypothesis
allows singular R_s. The main theorem and the OU proof in Section 2.2
were read and reconstructed: the interpolation derivative has the sign
of K I-R_s. No coordinate independence of h_s is asserted.

For any centered scalar X with MGF proxy sigma^2, its centered absolute
value has proxy 2sigma^2. Indeed, with an independent copy X', Jensen
and symmetry give

    E exp(lambda[|X|-E|X|])
      <=E cosh(lambda[|X|-|X'|])
      <=E cosh(lambda[X-X'])<=exp(sigma^2 lambda^2).  (3)

Thus (2) also supplies concentration of independent column responses.

## 2. Uniform absolute-response comparison, including singular endpoints

The following scalar consequence is ALREADY contained in the archived
[quenched Gaussian-sign universality theorem](flatify_independent_2026_09_07_gaussian_sign_quenched_universality.md).
Its complete conditional replacement, OU Stein-kernel, and two-star
contraction proof has been reread for this application.

If epsilon I<=R<=K I, diag R=1, and S=Cov(sign N(0,R)), then uniformly
over every Boolean x,

    |E|x dot sign(G)|-kappa sqrt(x^T S x)|
      <=C K^(2/3) epsilon^(-2/3) n^(1/3)sqrt(log(n+1)). (4)

To derive (4), apply the archived pressure theorem with dimension d=n,
two configurations {+x,-x}, and zero offsets. Dividing its pressure
error by the soft-max parameter lambda gives
C K^2 epsilon^(-2) lambda^2 n log(n+1)^(3/2); the two-state smoothing
error is at most 2log(2)/lambda. Optimize lambda. The Gaussian maximum
is exactly the Gaussian absolute value in (4). In particular, the
bound is uniform over all x; no union over 2^n spin words is required.
This is a scalar application of an existing local theorem, not a new
external central limit theorem.

A generic rank-one Gaussian-subordination CLT cannot simply be invoked:
the absolute row sum of correlations of R_s can grow as sqrt(n).
Bardet--Surgailis arXiv:1104.4732v2, Section 3, was inspected; its
rank-one summability hypothesis is not supplied here. Equation (4)
avoids that incorrect shortcut.

Here the flat full-sign structure permits the singular endpoints too.
Suppose |s|<=S0 and ||R_s||<=K0, where S0,K0 are fixed. Couple

    G_e=sqrt(1-e) G_s+sqrt(e) Z,
    R_e=(1-e)R_s+eI=I+(1-e)sA/sqrt(n),

with independent standard Z. Let D=sign(G_e)-sign(G_s) and u=sqrt(1-e).
Its covariance has the exact form

    Cov(D)=d_e I+c_e A,
    d_e=2-(4/pi)arcsin(u)=O(sqrt(e)),
    c_e=(2/pi)[arcsin(s/sqrt(n))
       +arcsin(u^2 s/sqrt(n))-2arcsin(u s/sqrt(n))].  (5)

For |s|/sqrt(n)<=1/2 and e<=1/2, the bracket in (5) is bounded in
absolute value by C e^2 |s|/sqrt(n). For example apply a second
finite-difference bound to f(v)=arcsin(exp(v)s/sqrt(n)); its second
derivative has magnitude O(|s|/sqrt(n)). Since
|s| ||A||/sqrt(n)<=||R_s-I||<=K0+1,

    ||Cov(D)||<=C_(K0) sqrt(e),
    E|x dot D|<=C_(K0) sqrt(n) e^(1/4).              (6)

The sign-covariance Gaussian comparator changes by O_(K0)(e) in
operator norm, using (1). Its absolute mean changes by at most
O_(K0)(sqrt(n e)). Combining (4)--(6) and taking e=n^(-1/8) proves

    E|h_s dot x|/sqrt(n)
      =kappa sqrt(1+2rho_n(s)H_A(x)/n)+O_(S0,K0)(n^(-1/32)), (7)

uniformly over all x and all such admissible s, including singular laws.
The slower n^(-1/32) absorbs n^(-1/12)sqrt(log n). Rates are not optimized.

## 3. Exact isotropy does not erase nonlinear absolute-response reduction

Take 0<=t<=sqrt(n)/||A||. Draw a fair pole sigma in {+-1}, then draw
h from the Gaussian-sign law R_(sigma t). Denote this sign law by nu_t.
Both covariances lie between zero and 2I. Equation (1) gives EXACTLY

    E_(nu_t) h=0,    E_(nu_t) h h^T=I.                (8)

The linear MGF proxy is at most 2, by (2), including after mixing.
For e_x=H_A(x)/n^(3/2), define

    m_t(e)=kappa/2 [sqrt(1+4te/pi)+sqrt(1-4te/pi)].    (9)

Since ||A||>=sqrt(n-1), t<=sqrt(n/(n-1)). Rayleigh also ensures
|2t e_x|<=1. Equations (1),(7) and the arcsine remainder show

    sup_x | E_(nu_t)|h dot x|/sqrt(n)-m_t(e_x) |
                                      <=C n^(-1/32). (10)

The constant is universal for n sufficiently large. In particular
m_t(e)<kappa whenever te!=0. It is even in e and decreases with |e|.
Thus if Q(A_n)/n^(3/2)->c>0 and ||A_n||<=Lsqrt(n), any fixed
0<t<=1/L yields simultaneously on the FULL absolute eta-nearcode

    sup E|h dot x|/sqrt(n)<=m_t(c-eta)+o(1),           (11)

for 0<eta<c, with n tending to infinity first. Both energy polarities
benefit under the SAME actual, exactly isotropic sign-column law.

The mechanism is concavity of standard deviation under a latent
variance mixture. Conditional variances are asymptotically
n(1+-4te/pi); averaging variances gives n, but averaging Gaussian
absolute means gives (9). Replacing the mixture by one Gaussian with
its averaged covariance would incorrectly erase the gain.

This is distinct from using an affine Gram law as a signed-energy
witness. The older
[affine-Gram exclusion](decisive_audit_nearmin_scalar_affine_gram_exclusion_2026_09_07.md)
bounds E H_A(h), explicitly not E|h dot x|. It does not imply or
contradict (10). No matching old absolute-response theorem was found
in the bounded archive search recorded with the frozen initial mechanism.

## 4. A fully paid actual-parent upper certificate

Let B have q independent columns with law nu_t. Let D be ANY actual
full-sign child on q vertices. For the genuine full-sign parent P,

    Q(P)<=max_x[|H_A(x)|+sum_(j<=q)|B_j dot x|]+Q(D). (12)

Partition ALL old words into nonempty energy bins C_j with
a_j<=|H_A(x)|/n^(3/2)<=b_j, and write M_j=|C_j|; there are J bins.
By (2),(3), every centered absolute column response has MGF proxy
4n. A union over all words and bins therefore proves, for every u>0,
with probability at least 1-exp(-u),

    Q(P)<=Q(D)+max_j {
       b_j n^(3/2)+q sqrt(n)[m_t(a_j)+C n^(-1/32)]
             +sqrt(8qn[log M_j+log J+u]) }.          (13)

This includes all new-spin choices through (12), both old energy
polarities, and every old word through the bins. The child Q(D) is
explicitly paid. Equation (13) is a valid original-parent certificate;
the response gain alone does not supply favorable bin counts or a
recurrence. In particular, a covariance calculation is not substituted
for a concentration or escape estimate.

## 5. A universal slope obstruction for the symmetric affine strategy

Set

    mu_star=kappa/2 [sqrt(1-2/pi)+sqrt(1+2/pi)]
            =sqrt((1+sqrt(1-4/pi^2))/pi)
            =0.7508551242385078... >3/4.              (14)

The strict inequality is elementary: squaring reduces it to
9pi^3(32-9pi)>1024. This polynomial is decreasing on [3,22/7], and
its value at 22/7 is 2491632/2401>1024. The standard bounds
3<pi<22/7 therefore suffice; the decimal is not a proof input.

Suppose Q(A_n)<= (1/2+o(1))n^(3/2). Choose ANY symmetric mixture of
the nu_t laws, with any distribution of admissible t (depending on A
and n). Uniformly on each absolute ground state x0, (10) and
t<=sqrt(n/(n-1)) give

    E|h dot x0| >= [mu_star-o(1)]sqrt(n).             (15)

No bounded-operator-norm hypothesis is required: normalized covariance
operators of every component are still <=2. A growing operator norm
only reduces available t and weakens the response discount.

Columns may have DIFFERENT such symmetric-mixture laws. If they are
independent, (3) gives, for every fixed delta>0,

    Pr{sum_(j<=q)|B_j dot x0|
          <q sqrt(n)[mu_star-o(1)-delta]}
                                  <=exp(-delta^2 q/8). (16)

Orient each new spin so its bridge term has the sign of H_A(x0).
For every child D, even if selected after observing B,

    Q(P)>=Q(A)+sum_j|B_j dot x0|-Q(D).               (17)

Thus the scheme has asymptotic extension slope at least mu_star, not
at most 3c/2 for c<=1/2. More precisely, let q/n->epsilon>0 and
Q(D)<=L q^(3/2)+o(n^(3/2)). With probability tending to one, its
excess above the same normalized parent cap satisfies

    [Q(P)-(1+q/n)^(3/2)Q(A)]/n^(3/2)
       >=epsilon mu_star
         -(1/2)[(1+epsilon)^(3/2)-1]
         -L epsilon^(3/2)-o(1).                      (18)

For every fixed L this is strictly positive for all sufficiently small
fixed epsilon, because mu_star>3/4. This is a scheme-specific
obstruction applicable to actual minimizing children whenever the
stated cap bound holds. It does NOT exclude other sign-column laws,
asymmetric laws favoring a uniquely dominant polarity, correlated
columns, or a low-cap parent produced by rare selection.

## 6. Sharp minimax obstruction for all scalar-affine mixtures on Hadamard children

The symmetric-mixture restriction can be removed on a concrete scalable
family. Let n=p^2, p a power of two, and take the alternating Walsh
matrix F from the
[all-law Hadamard artifact](paper_localization_all_law_hadamard_obstruction_2026_09_17.md),
so F^2=p^2 I and diag F=1. Put A=F-I. Then

    Q(A)=n(p+1)/2,
    H_A(x_+)=n(p-1)/2,   H_A(x_-)=-n(p+1)/2           (19)

for explicit Boolean eigenwords x_+,x_- supplied there. Both sectors
have deficit at most n, hence lie in every fixed positive near-level
eventually. The FULL admissible affine interval for R_s is

    -p/(p-1)<=s<=p/(p+1).                            (20)

Its covariance operators are uniformly bounded. By (7), uniformly in
s in this whole interval, the average response on these two words is

    [E|h_s dot x_+|+E|h_s dot x_-|]/(2sqrt(n))
        =kappa/2 [sqrt(1+2s/pi)+sqrt(1-2s/pi)]+o(1)
        >=mu_star-o(1).                              (21)

The same is true for EVERY probability mixture over (20), with no
symmetry, zero-covariance, or isotropy restriction. Therefore its worst
response on the absolute nearcode is at least mu_star-o(1).
Conversely the equal mixture at s=+-p/(p+1) is admissible, exactly
isotropic, and has response mu_star+o(1) uniformly on every code
whose absolute energies are (1/2-o(1))n^(3/2). Consequently

    lim_(eta down to 0) lim_(p to infinity)
      inf_(all affine Gaussian-sign mixtures)
       sup_(x in E_A(eta n^(3/2))) E|h dot x|/sqrt(n)
                                                   =mu_star. (22)

The inner liminf and limsup agree by the two bounds, so the displayed
limit is justified. For independent bridge columns with arbitrary,
possibly different affine-mixture laws, average the two total response
means in (21) and choose the better fixed sector before sampling.
That sector has old energy at least Q(A)-n. Equations (16)--(18)
then hold with an additional harmless loss n, and with denominator
8p/(p-1) in the exponent of (16) (or the uniform safe denominator 16),
since the whole interval (20) has covariance norm at most 2p/(p-1).
Thus even all asymmetric
scalar-affine mixtures fail the desired small-density slope on this
actual near-half family. These children are NOT asserted to minimize.

## 7. What is new, inherited, and still open

Inherited: exact Gaussian-angle identity; Gaussian Holder; the archived
all-offset quenched sign/Gaussian comparison. Proved here: its uniform
scalar and singular-flat endpoint consequences; exactly isotropic
two-polarity response discount; paid parent certificate; and the sharp
scalar-affine response minimax/slope obstruction above.

The isotropic mixture also gives an actual-child-designed example of
critical n-coordinate blocks with covariance I and a uniform subGaussian
constant whose absolute scalar responses differ macroscopically from
the covariance-matched Gaussian law. This agrees with, rather than
evades, the campaign's critical-block universality limitation.

Nonlinear Gaussian feedback, conditioning on the sampled energy,
non-affine covariance design, or correlated columns are not covered by
the negative theorem. No bounded spectral norm is inferred for exact
minimizers, and no new original cap bound or convergence result is claimed.

Audit status: the Bernoulli researcher independently read this entire
artifact and the complete archived quenched theorem, reconstructed the
endpoint coupling, scalar rate, isotropy, parent inequalities, and sharp
minimax, and returned PASS after the displayed finite tail-constant
correction in Section 6. The correction does not change any asymptotic
slope or minimax conclusion.
