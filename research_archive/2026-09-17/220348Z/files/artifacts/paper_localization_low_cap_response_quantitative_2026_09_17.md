# Quantitative audit: low cap forces a uniform cheap physical response law

2026-09-17. Localization-track explicit margin-to-discount extension of
the director's [low-cap uniform-response theorem](paper_director_low_cap_uniform_response_2026_09_17.md).
The central structural argument is the director's proposed combination of
the discrepancy track's thin signed-energy alternative, Grothendieck
coordinate localization, and the inherited HALF-RANGE lower bound.
This track independently reconstructed that combination and supplies the
explicit sixth-power discount below. No spectral condition is assumed.

## 1. Statement, with the half-range normalization made explicit

Write H_A(x)=sum_(i<j)a_ij x_i x_j and Q(A)=max_x|H_A(x)|. Define

    R(A)=[max_x H_A(x)-min_x H_A(x)]/2.

Thus R is HALF the energy range, and R(A)<=Q(A). Let r_*>0 be any
fixed constant such that every sufficiently large full signing B of
order m satisfies R(B)>=r_* m^(3/2). The archived paired construction
proves this for every r_*<0.4297864507376; its original full-range
conclusion is twice that number. The stronger current lower bound may
be substituted, but is not needed here.

Let K be any valid real Grothendieck constant. One may use K=2, or
the explicit K=pi/[2 asinh(1)]<1.783 proved by the odd-tensor/hyperplane
argument in the [range artifact](fresh_range_and_spectral_regularization_2026_09_05.md).
Fix C,c>0 with

    g:=2c-4K(C-r_*)>0.                              (1)

For every sufficiently large n and EVERY actual hollow full signing A
with Q(A)<=C n^(3/2), there is one centered, EXACTLY isotropic physical
sign law nu on{+-1}^n, with subGaussian proxy(3/2)I, such that

    sup_(|H_A(x)|>=c n^(3/2)) E_nu|h dot x|/sqrt(n)
                              <=kappa-delta,
    kappa=sqrt(2/pi),                              (2)

where delta>0 depends only on c,C,r_*,K. The law belongs to the compact
convex hull of equal Gaussian-sign pairs with latent spectra in[1/2,3/2].
All coefficients h_i are actual signs. The same one law protects the
ENTIRE displayed macroscopic energy code, not an average over a selected
ground law. No operator bound on A or query covariance is required.

Here is a completely explicit, conservative discount. Put

    b_0=1/2-pi/12,
    C_1=sqrt[8/(kappa^5 b_0^2)],
    C_2=[373248/(25 kappa^5)]^(1/6),
    L=C_2+(5/sqrt(2))sqrt(C_1)
                          +6K(C+r_*)C_1^(1/3),
    d_0=min{1,(8C_1)^(-2),(g/(2L))^6},
    delta=d_0/2.                                   (3)

The dimension threshold ensures the uniform scalar Gaussian-sign
comparison error is at most d_0/2 and the half-range bound applies on
every induced submatrix of order at least n/2. The formula is not an
optimized numerical response bound or a favorable parent slope.

## 2. Quantifying the thin signed-energy alternative

Fix ANY query law mu on the high-energy code, and put

    sigma(x)=sign H_A(x),  Sigma=E xx^T,
    K_mu=E sigma xx^T,  m=E sigma,
    K_0=K_mu-mI,
    Q=1_(Sigma>2),  P=I-Q,
    a=tr(Sigma Q)/n,  B=PK_0P,  b=||B||_F^2/n.

The projection Q in this section is not the cap Q(A). We use the
paired-law scalar comparison reconstructed in the
[thin-energy alternative](paper_discrepancy_signed_energy_concentration_alternative_2026_09_17.md):
for a hollow H with ||tH||op<=1/2, the equal pair sign N(0,I+-tH)
has exactly covariance I and

    E_nu|h dot x|/sqrt(n)=kappa f(v_x)+error_x,
    f(v)=[sqrt(1+v)+sqrt(1-v)]/2<=1-v^2/8,
    v_x=kappa^2 x^T asin(tH)x/n,
    sup_x|error_x|<=epsilon_n,
    epsilon_n=O(n^(-1/6)sqrt(log(en))).               (4)

Arcsine is entrywise. Its uniformity requires only the latent covariance
interval[1/2,3/2], not a bound on A. The pair also has linear MGF
at most exp(3||theta||^2/4), including after mixing such pairs.

Suppose for a contradiction that EVERY such paired law has mu-average
response at least kappa-delta. Set d=delta+epsilon_n. Applying (4)
to H=Q-diag(Q), t=1/2 gives

    E_mu v_x>=kappa^2 b_0 a,
    d>=(kappa^5 b_0^2/8)a^2,
    a<=C_1 sqrt(d).                                 (5)

The first line uses entrywise positive odd powers of the PSD projector
and rank(Q)<=an/2, exactly as in the source alternative.

Since -Sigma<=K_mu<=Sigma, we have ||B||op<=3 and b<=9. Let
H=B-diag(B), so ||H||op<=6, and choose t=sqrt(b)/36 if b>0.
This always obeys t<=1/12. The exact identity
<K_mu,H>=||B||_F^2 and the Schur-power remainder bound give

    E_mu sigma(x)v_x
       >=kappa^2[tb-216t^3]
        =(5kappa^2/216)b^(3/2).

Using f(v)<=1-v^2/8 again,

    d>=25 kappa^5 b^3/373248,
    sqrt(b)<=C_2 d^(1/6).                           (6)

The b=0 case is immediate. Importantly, the same effective error d
works in both tests; neither law is being claimed to protect every
query at this intermediate stage.

Set Y=Qx and J=E_mu sigma Y^T A Y. The exact block expansion, Frobenius
Cauchy--Schwarz, and P Sigma P<=2P give

    J/n^(3/2)>=2c-sqrt(b)-(5/sqrt(2))sqrt(a).         (7)

The cross term is at most sqrt(2a)n^(3/2), the diagonal-centering
term at most sqrt(a/2)n^(3/2), and |<A,B>|<=sqrt(b)n^(3/2).
Only ||A||F=sqrt(n(n-1)) was used.

## 3. Coordinate localization: why a thin covariance subspace is enough

Because Q is a spectral projection of Sigma,

    V=E YY^T=Q Sigma Q<=Sigma,
    0<=V_ii<=1,  sum_i V_ii=an.

For a>0 set epsilon=a^(2/3) and

    S={i:V_ii>epsilon},  T=S^c.

Then |S|<=a^(1/3)n. In L2(mu), the vectors Y_i and sigma Y_i both
have norm sqrt(V_ii). Grothendieck applied separately to the S,S;
T,T; and S,T blocks therefore yields

    J<=K beta(A_SS)
                +K epsilon beta(A_TT)
                +2K sqrt(epsilon) beta(A_ST),       (8)

where beta denotes the appropriate square or rectangular Boolean
bilinear maximum. The last term has the factor2 because A is symmetric.

Two EXACT scalar facts pay these maxima with the correct normalization:

    beta(A_SS)<=4R(A_SS),
    beta(A_TT)<=4Q(A),
    beta(A_ST)<=Q(A).                               (9)

For the first, write u=(x+y)/2 and v=(x-y)/2. They have disjoint
supports in the cube and x^T A y=2[H_A(u)-H_A(v)]. Completing partial
spins independently and fairly shows H_A(u)<=max H_A and
H_A(v)>=min H_A. Thus beta<=2 full_range=4 HALF_RANGE. The second
uses the same bound by4Q on the induced block and induced-cap
monotonicity. For the rectangular third bound, reversing every spin
of one block preserves its two internal energies and reverses the
bridge; max(|a+b|,|a-b|)>=|b|.

It follows that

    J/n^(3/2)<=4K R(A_SS)/n^(3/2)
                     +KC[4a^(2/3)+2a^(1/3)].       (10)

This step does not declare the vectors Y Boolean. It converts their
small average coordinate variances to a genuine induced-coordinate
range budget using the classical vector-to-sign inequality.

## 4. Full-sign complement range closes the contradiction

Whole-block reversal proves positive and negative extrema separately
superadditive over a vertex partition. Hence

    Q(A)>=R(A)>=R(A_SS)+R(A_TT).                    (11)

When a<=1/8, the complement has at least n/2 vertices. The inherited
uniform half-range lower bound therefore gives

    R(A_SS)/n^(3/2)<=C-r_*(1-a^(1/3))^(3/2).         (12)

The case a=0 uses S empty and J=0 and obeys the same resulting bounds.
Combining (7), (10), and (12), and using
1-(1-z)^(3/2)<=3z/2 and a^(2/3)<=a^(1/3), gives

    g<=sqrt(b)+(5/sqrt(2))sqrt(a)
                            +6K(C+r_*)a^(1/3).

If d<=1, (5)--(6) and d^(1/4)<=d^(1/6) imply

    g<=L d^(1/6).                                   (13)

For delta=d_0/2 and epsilon_n<=d_0/2, we have d<=d_0. Formula (3)
ensures a<=1/8 and Ld^(1/6)<=g/2, contradicting (13). Thus every
query law admits a physical paired column law with average response
at most kappa-delta.

Finite-dimensional minimax over the finite Boolean query set and the
compact convex hull of these paired laws now supplies ONE law satisfying
(2) uniformly. Compactness follows from continuity of Gaussian orthant
probabilities on the strictly positive latent spectral interval.
Centering, exact isotropy, and the uniform(3/2) subGaussian proxy are
preserved under convex mixtures. No limiting nonphysical law is used.

## 5. Concrete scope and remaining parent cost

For a deliberately robust example, take K=2, r_*=0.429, C=0.5 and
c=0.30. Then g=0.032>0. Therefore all sufficiently large actual full
signings with cap at most0.5 n^(3/2) admit such a law protecting EVERY
word of absolute energy at least0.30 n^(3/2). This includes a fixed
macroscopic nearlevel of current nearoptimal families. The particular
discount in (3) is extremely small; this example is not a numerical
improvement of the original cap constant.

The genuine nonground high-energy obstruction with larger total cap
does not satisfy (1), so it is consistent with this theorem. The theorem
also does not assume that optimizers have bounded covariance: it handles
every dual law, including coherent sublinear covariance spikes.

An iid bridge sampled from the resulting law still requires a full
old-word/new-spin maximum estimate. Its fixed response discount alone
does not pay the entropy or all-energy escape term, and its size is far
below the known favorable extension-slope deficit. This remaining cost
must be displayed in any attempted parent construction.

Independent audit: the discrepancy researcher read Sections 1--5 in full
and returned PASS, including C_1,C_2, the5/216 remainder coefficient,
the373248 denominator, complement normalization, sixth-power choice,
and minimax uniformity. This track also read the director's complete
canonical low-cap theorem and reconstructed its asymptotic proof.

The reproducible constants replay
`computations/paper_localization_2026_09_17_low_cap_response_constants.py`
passes1,200 algebra checks for the remainder normalization, localization
exponents, complement budget, and sixth-power margin. Its output is
`tmp/paper_portfolio_2026_09_17/localization/low_cap_response_constants.json`.
It is not a finite-dimension verification of the imported asymptotic
Gaussian-sign scalar theorem.
