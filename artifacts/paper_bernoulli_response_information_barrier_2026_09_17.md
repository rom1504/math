# Cheap physical response against an isotropic dual requires extensive information

2026-09-17. Scoped physical-law barrier. The concentration ingredient is
classical; the purpose here is to keep the query-law and physical-law
quantifiers separate and test the proposed low-information rewrites.

## 1. A finite theorem

Let mu be ANY probability law on Boolean queries x in{+-1}^n, and assume

    ||E_mu xx^T||op<=L.

There is no independence or energy hypothesis on mu. Let pi be uniform
on the physical cube and let mu_n=E_pi|sum_i h_i|. For ANY physical law
nu on that cube, put

    R(mu,nu)=E_(x~mu,h~nu)|h.x|,
    Delta=(mu_n-R(mu,nu))_+.

Then, with natural-log relative entropy,

    D(nu||pi)>=Delta^2/(64L)-2.                     (1)

In particular, if mu is isotropic, reducing its average normalized
response by a fixed delta>0 below kappa=sqrt(2/pi) requires

    D(nu||pi)>=(delta^2/64-o(1))n.                 (2)

The physical law need not be centered or isotropic for this obstruction.
Conversely, the query law's second-moment bound is essential to the
argument. A small-rank covered code need not admit an isotropic query
law; there is no contradiction with a positive response construction
under such a geometric premise.

## 2. Primary concentration mechanism and exact cube scaling

We use Talagrand's convex-distance theorem from the primary paper
[Concentration of measure and isoperimetric inequalities in product spaces](https://arxiv.org/pdf/math/9406212),
Section4.1, Theorem4.1.1. That entire section's theorem and induction
proof were read for this application. The theorem says, for a product
law and a set A,

    E exp[d_T(h,A)^2/4]<=1/pi(A).                  (3)

Here d_T is the Euclidean distance from zero to the convex hull of
coordinate mismatch vectors. Its proof uses fiber/projection induction:
for a fiber A(omega), projection B and mixing parameter lambda,

    d_T((x,omega),A)^2
      <=(1-lambda)^2+lambda d_T(x,A(omega))^2
                         +(1-lambda)d_T(x,B)^2.

Hölder and the induction hypothesis reduce the estimate to

    inf_(0<=lambda<=1) r^(-lambda)exp[(1-lambda)^2/4]<=2-r,

and integrating the fibers finishes by z(2-z)<=1. Thus the imported
principle is a product-space convex-distance estimate, not a claim
that arbitrary Euclidean-Lipschitz functions concentrate on the cube.

For the sign cube specifically,

    dist_2(h,conv A)=2 d_T(h,A).                   (4)

Indeed a convex combination of mismatch vectors maps to a convex
combination of words by y_i=h_i(1-2v_i). Let f be convex and
sqrt(L)-Lipschitz, and let m be a median of f under pi. For
A_t={f<=m-t}, convexity puts conv A_t inside that same sublevel set.
Every h with f(h)>=m therefore satisfies
d_T(h,A_t)>=t/(2sqrt(L)). Applying(3) on a set of probability at least
one half gives

    pi(f<=m-t)<=2 exp[-t^2/(16L)].                 (5)

Taking A={f<=m} gives the analogous upper tail. In particular

    pi(|f-m|>=t)<=4 exp[-t^2/(16L)].               (6)

This proves BOTH tails using convex sublevel sets. No concavity of f
or unjustified lower-tail bounded-differences estimate is needed.

## 3. Application to the response function and entropy duality

The function

    f(z)=E_mu|x.z|

is convex, and Cauchy--Schwarz gives
|f(z)-f(z')|<=sqrt(L)||z-z'||. Its cube mean is EXACTLY mu_n, because
each fixed Boolean query only switches the independent uniform signs.

For completeness, uniform-cube Efron--Stein improves the needed
mean/median comparison without integrating large tail constants.
Choose a subgradient g(h) with ||g(h)||<=sqrt(L). Convexity implies

    sum_i[(f(h)-f(h^(i)))_+]^2<=4L.

Uniform resampling flips each coordinate with probability1/2; the
resampled-pair symmetry in Efron--Stein therefore gives Var_pi f<=2L.
Consequently |E_pi f-m|<=2sqrt(L), since each median side has mass at
least one half. Integrating(6) yields

    E_pi exp[(f-m)^2/(32L)]<=5,
    E_pi exp[(f-E_pi f)^2/(64L)]<=5 exp(1/8).       (7)

The last step uses (u+v)^2<=2u^2+2v^2. Thus for all s>=0, Young's
inequality sy<=y^2/(64L)+16Ls^2 gives

    log E_pi exp[s(E_pi f-f)]
       <=16Ls^2+log5+1/8.                         (8)

Relative entropy's elementary variational inequality gives

    D(nu||pi)>=s(E_pi f-E_nu f)-log E_pi exp[s(E_pi f-f)].

Optimize at s=Delta/(32L). Since log5+1/8<2, this proves(1).
This argument is finite; it does not assume a limit for mu or nu.

## 4. A simpler L2 consequence and the polynomial-density scope

If w=dnu/dpi, the same variance calculation directly gives

    |R(mu,nu)-mu_n|<=sqrt(2L)||w-1||_2.            (9)

Thus bounded-Renyi2-cost rewrites cannot produce a leading normalized
discount against an isotropic query dual. The KL theorem is stronger:
even subexponential L2 costs may have small relative entropy, and(1)
still forbids a fixed discount whenever D=o(n) and L=O(1).

As a concrete corollary, every nonnegative Boolean polynomial density
w of fixed degree D, normalized by E_pi w=1, has ||w||_2<=3^D. Indeed
Boolean hypercontractivity gives ||w||_4<=3^(D/2)||w||_2; interpolation
||w||_2<=||w||_1^(1/3)||w||_4^(2/3) gives the claim. Therefore ALL
such densities, not only polynomials in the signing energy, have

    R(mu,nu)=mu_n+O_D(sqrt(L)).                    (10)

This is a different statement from the earlier fixed-energy-filter
theorem, which controls every query without a query-covariance
hypothesis. Here the law may be an arbitrary fixed-degree physical
density, but the conclusion is an average against the stated dual.

### Growing degree: a linear degree cost for a fixed discount

The director observed that combining hypercontractivity with the KL
bound is stronger than the preceding L2 comparison. Let0<=D<=n, and
let w be ANY nonnegative normalized Walsh-polynomial density of degree
at most D. Its coefficients, degree, and query law may all depend on n.
The same finite hypercontractive calculation gives

    D(nu||pi)=E_nu log w
       <=log E_nu w=log E_pi w^2<=2D log3.         (11)

The entropy inequality is Jensen under nu; zero values of w cause no
problem because nu assigns them zero mass. Combining(1) and(11),

    (mu_n-R(mu,nu))_+
       <=8 sqrt[L(2D log3+2)].                    (12)

Thus a deficit at least delta sqrt(n) from the exact independent-sign
mean requires

    D>=delta^2 n/(128L log3)-1/log3.               (13)

For any sequence with bounded L and D=o(n), no fixed normalized
response discount is possible. More generally the conclusion holds
whenever L(D+1)=o(n). This is a finite statement about actual positive
physical densities, not just fixed powers or polynomial witnesses of
an energy. It is one-sided: sublinear degree rules out cheap average
response, not necessarily an expensive response or all-query closeness.

The degree here is the Walsh degree of the NONNEGATIVE normalized
density itself. It is not the degree of a test polynomial later squared,
a formal Gaussian expression, or the original quadratic signing.
The query-covariance assumption is indispensable. These qualifications
separate this linear-degree obstruction from the earlier scalar
Chebyshev moment-recovery threshold, despite the similar linear scale.

## 5. Exact relevance and limitations

The quartic feature laws have density L2 squared below3, so(9) already
applies. Their special polynomial structure gives the still sharper
exact bound recorded in the
[quartic-law artifact](paper_localization_quartic_isotropic_laws_2026_09_17.md):
against an isotropic query law their average discount is only O(1/n)
after normalization, uniformly over rank and frame.

Hot/cold feature tilts with information O(r)=o(n), whenever that
information estimate and the physical law are established, also cannot
uniformly protect a code carrying an isotropic query dual by a fixed
discount. This does NOT prevent their useful action on a low-rank
center code. Such a code necessarily has a large covariance direction
if most mass lies near that low-rank space.

The conclusion does not prohibit successful nonlocal laws. For example,
uniform laws on small Hadamard dictionaries can have information of
order n and are outside the low-information regime. Nor is isotropic
query mass established for arbitrary exact minimizing ground codes.
The finite order12 and14 examples have such dual laws; they motivate
the test but do not prove an all-order optimizer geometry theorem.

The general lower bound is an application of a classical concentration
theorem and entropy duality. No independent novelty claim is made for
this information-concentration principle.

Both the director and discrepancy track independently read and
reconstructed the complete proof: PASS, including the convex lower
tail, uniform-cube variance factor, and entropy constants. The
growing-degree corollary uses the same finite hypercontractive bound;
it does not extrapolate an unspecified fixed-degree error constant.
