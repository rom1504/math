# Smooth finite energy witnesses cannot lower the leading sign response

2026-09-17. A scoped negative extension of the reconstructed joint
replacement mechanism in the
[fixed-power theorem](paper_bernoulli_fixed_power_tilts_2026_09_17.md).
This concerns central quadratic-energy witnesses of ACTUAL full signs,
not arbitrary features, sharp rare conditioning, or all Gibbs laws.

## 1. Uniform finite-witness theorem

Fix an integer r, constants C,M,L<infinity and m>0. For each n let
A_1,...,A_r be actual symmetric hollow full signings satisfying

    Q(A_j)<=C n^(3/2),
    Y(h)=(H_(A_1)(h)/n,...,H_(A_r)(h)/n).

The signings need not be independent, spectrally flat, or minimizing.
Let w_n:R^r->[0,M] be ANY L-Lipschitz function, possibly chosen using
the signings and n, such that Z_n=E_Uniform w_n(Y(h))>=m. The genuine
physical sign law

    dnu/dUniform(h)=w_n(Y(h))/Z_n

satisfies, uniformly over EVERY Boolean query x,

    |E_nu |h.x|/sqrt(n)-kappa|
       <=C_(r,C)(M+L) m^(-1) n^(-1/4).             (1)

Every such law is centered because its density is globally even.
Exact isotropy is not required: (1) therefore also applies to any
exactly isotropic law that happens to belong to this class. The
theorem excludes a fixed leading cheap-response gain from any
uniformly bounded smooth function of a fixed list of normalized
central energies. Mixtures and interactions among those witnesses
are allowed through the arbitrary multivariate function w_n.

The quantitative estimate also allows varying M_n,L_n,m_n whenever
(M_n+L_n)/(m_n n^(1/4))->0. It does NOT cover indicators with shrinking
transition regions unless their actual Lipschitz constants obey that
condition, and does not discard the normalization cost of rare events.

## 2. Joint replacement, including a merely Lipschitz likelihood

Fix x and put Z(U)=x.U/sqrt(n), where the independent coordinates
of U may each be a fair sign or a standard Gaussian. At coordinate i,

    Z=z+a U_i,     Y=y+b U_i,
    |a|=n^(-1/2),  b_j=n^(-1)sum_(k!=i)(A_j)_ik U_k.

The remaining variables z,y,b are independent of U_i. For every fixed
k, direct moments of independent centered bounded/Gaussian linear
sums give

    ||z||_k<=C_k,      || ||b||_2 ||_k<=C_(r,k)n^(-1/2). (2)

No bound on y is needed below because the likelihood is bounded.
Convolve w with a fixed nonnegative smooth compactly supported
mollifier of radius delta. The resulting w_delta obeys

    ||w_delta||_infinity<=M,
    ||w_delta-w||_infinity<=C_r L delta,
    ||D^j w_delta||_infinity<=C_(r,j)L delta^(1-j), j>=1.

Put phi_delta(z)=sqrt(z^2+delta^2), with 0<delta<=1. Then
|phi_delta-|z||<=delta and its derivatives of orders1,...,4 have
size at most C_j delta^(1-j). Consider the one-coordinate test

    F(t)=phi_delta(z+a t) w_delta(y+b t).

Its fourth derivative is a sum of terms with j likelihood derivatives,
0<=j<=4, carrying a^(4-j) and j copies of b. Terms1<=j<=3 have
derivative loss at most delta^(-2). The j=0 term has loss M delta^(-3).
The j=4 term has loss L delta^(-3) times
|z|+|a t|+1. Thus (2), Holder, and the fixed moments of the removed
sign/Gaussian driver imply

    E[|U_i|^4 sup_(|t|<=|U_i|)|F''''(t)|]
       <=C_r(M+L) delta^(-3)n^(-2).                (3)

The first three driver moments agree. Taylor expansion at zero,
followed by n replacements and the likelihood/absolute smoothing
errors, proves

    |E_Rad |Z| w(Y)-E_Gauss |Z| w(Y)|
       <=C_r(M+L)[delta+n^(-1)delta^(-3)].          (4)

The same bound applies without |Z|; that case is easier and uses
only the fourth derivative of w_delta. Choose delta=n^(-1/4).
This is a fully paid replacement of the literal quadratic witnesses;
no Boolean multilinearization of a composite density is evaluated at
Gaussians, and no independence among the witnesses is asserted.

## 3. The actual Boolean cap separates the Gaussian linear query

The elementary polarization argument from the fixed-power theorem
gives for each actual signing and every Boolean x

    ||A_j x||_1<=4Q(A_j),
    ||A_j x||_2^2<=4(n-1)Q(A_j)<=4C n^(5/2).       (5)

Let e=x/sqrt(n), and decompose a standard Gaussian G=Z e+V, where
Z is standard normal independent of V in e-perpendicular. Simultaneously
for all j,

    Y_j(G)=D_j+Z L_j+a_j Z^2,
    D_j=V^T A_j V/(2n),
    L_j=e^T A_j V/n,
    a_j=H_(A_j)(x)/n^2.

The vectors D,L need not be independent, but both are independent of Z.
By (5),

    ||a||_2<=sqrt(r) C n^(-1/2),
    E||L||_2^2<=4rC n^(-1/2).

The Lipschitz property therefore gives

    |E|Z|w(Y(G))-kappa E w(Y(G))|
      <=L E[(|Z|+kappa)||ZL+aZ^2||_2]
      <=C_(r,C)L n^(-1/4).                         (6)

Combining (4), its unweighted version, and (6), then dividing by the
ACTUAL sign normalization Z_n>=m, proves (1). All estimates are uniform
in x, the signings, and the chosen likelihood under their stated bounds.

## 4. What this distinguishes, and what it does not

The exact non-Gaussian laws in
[the quartic construction](paper_localization_quartic_isotropic_laws_2026_09_17.md)
illustrate the distinction. A quartic in a fixed-rank LINEAR feature
can be order one along the normalized query direction and yields a
fixed physical response discount with exact isotropy and bounded
information cost. A smooth fixed list of central ACTUAL quadratic
energies is asymptotically independent of every normalized Boolean
linear query in the specific weighted-absolute sense (1).

This does not exclude growing witness dimension, nonsmooth or rare
energy conditioning, likelihood parameters outside the quantitative
range, a thermodynamic Gibbs tilt exp(beta H/sqrt(n)), or nonlocal
laws not represented by the stated likelihood. It is a mechanism
falsifier, not a statement that marked energies contain no useful
information. In particular the quartic marked-energy law changes
energy variance by a fixed amount even while its responses stay at
kappa; its separate exact algebra does not need this bounded-likelihood
theorem.

Any vanishing-mass convex correction of one of these laws still has
response at least (1-o(1))kappa sqrt(n) at every query, simply by
nonnegativity of the added law's response. Thus a vanishing-mixture
isotropy repair cannot turn this particular central-energy mechanism
into a fixed cheap-response improvement.

Original ingredients are the elementary mixed-product moment bounds,
fourth-order moment-matching replacement, and Boolean polarization.
The fixed-power proof was read in full; the multivariate bounded-
Lipschitz extension and its exact quantifiers are proved here. External
novelty is not asserted.

The Bernoulli researcher independently read and reconstructed Sections1--4
in full, including the multivariate mollification, Gaussian separation,
normalization and vanishing-mixture scope: PASS.
