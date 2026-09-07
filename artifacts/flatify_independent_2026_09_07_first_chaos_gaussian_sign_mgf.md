# A paid first-chaos MGF bound for actual Gaussian-sign bridges

Later refinement: `flatify_independent_2026_09_07_exact_gaussian_sign_mgf.md`
removes the Holder inflation when the Gaussian covariance has a paid
positive spectral gap. `flatify_independent_2026_09_07_gaussian_sign_quenched_universality.md`
then transfers the full quenched pressure, including unchanged children,
to the exact-covariance Gaussian process. The theorem below remains
valid even for singular bounded-spectrum covariances and is preserved
as the earlier fully paid bound.

2026-09-07. A full exponential bound, not a variance-to-tail inference.
It retains an energy-dependent favorable term but pays an explicit
higher-chaos inflation. It does not establish the original composition.

## 1. General finite theorem

Let G be a centered Gaussian vector in R^d with correlation matrix R,
diag R=1 and ||R||op<=K, where K>=1 is fixed. Put

    a=sqrt(2/pi),       r(z)=sign(z)-a z.

For any real source vector h with ||h||infinity sufficiently small
(depending only on K),

 log E exp<h,sign G>
 <=(a^2/2)h^T R h+[K(1-a^2)/2]||h||^2
                       +C_K ||h||infinity ||h||^2.       (1)

The constant is independent of dimension, correlation coherence, and
the signs or geometry of h. Singular R is permitted.

### Exact Cameron--Martin step

Set m=aRh. Gaussian exponential tilting, including on the support of a
singular Gaussian, gives exactly

 E exp<h,sign G>
 =exp[(a^2/2)h^T R h] E exp[sum_i h_i r(G_i+m_i)].       (2)

Since ||m||^2<=a^2 K^2||h||^2, individual large tilted means ("hubs")
are allowed; no infinity bound on m is assumed.

### The imported Gaussian Holder inequality

For nonnegative measurable f_i, the matrix condition R<=KI implies

    E product_i f_i(G_i)
       <=product_i [E f_i(Z)^K]^(1/K),                 (3)

where Z is standard normal. This is Theorem 1(i) of Chen, Dafnis and
Paouris, *Improved Holder and reverse Holder inequalities for Gaussian
random vectors*, arXiv:1306.2410v2, with scalar blocks and every p_i=K.
The hypotheses were checked directly: the block marginal variances are
one and their block-diagonal P is exactly KI. The theorem allows the
nonnegative measurable functions used here; Gaussian plus bounded
linear growth makes all their required moments finite.

Primary source: https://arxiv.org/pdf/1306.2410, pp. 1--2, Theorem 1.

Apply (3) to f_i(z)=exp[h_i r(z+m_i)]. It remains to expand ONLY
one-dimensional Gaussian expectations, not the correlated sign law.

### Uniform scalar expansion, including arbitrarily large shifts

Let f(m)=2Phi(m)-1, so f'(m)=a exp(-m^2/2). For standard normal Z,

    mu(m)=E r(Z+m)=f(m)-a m,
    v(m)=Var r(Z+m)=1-f(m)^2+a^2-4a phi(m).

Elementary Taylor bounds at zero and boundedness away from zero give

    |mu(m)|<=C m^2,
    |v(m)-(1-a^2)|<=C m^2                            (4)

for EVERY real m. For example |f(m)-a m| is bounded by a constant times
min(|m|^3,|m|), which is at most a constant times m^2.

The centered random variable

    r(Z+m)-mu(m)=sign(Z+m)-f(m)-aZ

is bounded in absolute value by 2+a|Z|, uniformly in m. Thus its
exponential third moments are uniformly bounded on every fixed compact
interval of the scalar source. Taylor's formula for the log MGF gives

 log E exp[u r(Z+m)]
    <=u mu(m)+(u^2/2)v(m)+C|u|^3                    (5)

uniformly in m and |u|<=u_0. One direct justification differentiates
the log MGF three times: under every such tilt, its third centered
moment is bounded by the common Gaussian exponential-moment envelope;
the denominator is bounded below by Jensen after centering.

Use u=K h_i in (5), divide by K, and sum. The mean error in (4) is at
most C||h||infinity||m||^2. The variance error is at most
CK||h||infinity^2||m||^2. The third-order error is at most
CK^2||h||infinity||h||^2. Substituting the norm bound on m proves (1).
This proof explicitly pays the hub regime rather than silently using a
small-coordinate Cameron--Martin shift.

## 2. Opposite-spectral bridge: the general actual law

For any hollow symmetric sign A of order n, let s=sqrt(n-1) and put

    K_A=|A|+diag(s-diag|A|).

The diagonal correction is nonnegative because
|A|_ii<=sqrt((A^2)_ii)=s. Therefore K_A>=|A|>=+/-A and diag K_A=s.
For a second child D define K_D likewise. The matrix

    R=(K_A tensor K_D-A tensor D)/(n-1)               (6)

is a correlation matrix: its numerator equals

  1/2[(K_A+A) tensor (K_D-D)+(K_A-A) tensor(K_D+D)],

a sum of PSD matrices, and its diagonal is n-1. Sample a Gaussian
matrix G with vec covariance R and use the ACTUAL full-sign bridge
C_ij=sign G_ij. This constant-diagonal repair is the director's useful
normalization correction; it preserves unweighted child energies.

For any Boolean x,y, v=x tensor y has

    v^T R v=[(x^T K_A x)(y^T K_D y)
                        -(x^T A x)(y^T D y)]/(n-1).   (7)

If ||R||op is bounded by a fixed K, (1) with h=t v/sqrt(n) gives

 log E exp[t x^T C y/sqrt(n)]
 <=(nt^2/2) { a^2 (v^T R v)/n^2+K(1-a^2) }
                                                    +O_K,t(sqrt(n)). (8)

This is a rigorous energy-sensitive fixed-pair exponential estimate.
However, an arbitrary low-cap child only has ||A||op=O(n^(3/4)); that
does NOT guarantee bounded K in (6). Moreover, the positive K_A,K_D
responses in (7) still require payment. Neither condition is inferred
from child optimality here.

## 3. Exact flat specialization and its explicit inflation

If A^2=D^2=(n-1)I, the covariance simplifies to

    R_rho=I-rho A tensor D/(n-1),       0<=rho<=1,

with ||R_rho||op=1+rho. Every row and every column of G has independent
standard Gaussian coordinates. For e_A=H_A(x)/n^(3/2) and similarly
e_D, define

    V_rho(e_A,e_D)
       =1+rho(1-2/pi)-(8rho/pi)[n/(n-1)]e_A e_D.

Then (8) becomes

 log E exp[t x^T C y/sqrt(n)]
       <=nt^2 V_rho(e_A,e_D)/2+O_t(sqrt(n)).           (9)

For fixed b>0, Chernoff therefore proves the actual two-sided tail

 P(|x^T C y|>=b n^(3/2))
       <=2exp[-n b^2/(2V_rho(e_A,e_D))+O_b(sqrt(n))]. (10)

The denominator is bounded below by a positive absolute constant, so
the optimizing t is bounded for bounded b. Same-polarity child energy
reduces this paid tail proxy. The exact covariance alone would instead
have proxy 1-(8rho/pi)e_Ae_D+o(1). The additional
rho(1-2/pi) in (9) is the explicit price of Holder on the nonlinear
Hermite residual; it has NOT been removed.

At rho=1, even optimizing diagonal Holder exponents cannot reduce
their average below two: R_1^2=2R_1 and diag R_1=1; for any diagonal
P>=R_1, Tr P=Tr(P R_1)>=Tr(R_1^2)=2n^2. Thus the uniform-source
residual inflation is intrinsic to this diagonal Holder certificate,
not an accidentally loose choice of the scalar K.

The exact sign covariance in this flat case is also available:

    Cov(vec C)=I-kappa A tensor D,
    kappa=(2/pi)arcsin(rho/(n-1)).

All nonzero off-diagonal Gaussian correlations have the same absolute
value, so there is no need for an arcsine remainder estimate. This
identity by itself would not prove (10).

## 4. Direct finite parent implication, with the missing entropy explicit

For fixed actual children, the parent cap equals

    max_(x,y) [|H_A(x)+H_D(y)|+|x^T C y|].

Partition child energies into finitely many intervals before n grows.
If a bin pair has spin counts at most exp(n s_A),exp(n s_D), child
absolute-sum ceiling q, and a uniform proxy ceiling V from (8) or (9),
then its contribution is controlled at parent threshold L n^(3/2)
whenever

    (L-q)^2/(2V)>s_A+s_D,            L>q.              (11)

Strict margins pay the O(sqrt(n)) term and the finite union of bins.
This is a genuine sufficient original-parent inequality for this
actual sign ensemble. No theorem in this note supplies (11) for every
bin of arbitrary optimizing children. In particular, a variance gain
must not be substituted for the full exponential proof above, and a
bounded total child cap does not provide the necessary shell counts.
