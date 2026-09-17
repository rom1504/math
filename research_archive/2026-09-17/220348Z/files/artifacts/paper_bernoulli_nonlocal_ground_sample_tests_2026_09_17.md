# Physical nonlocal ground-sample tests: exact covariance is not enough

2026-09-17. Explicit mechanism tests on complete ground codes of stored
actual finite minimizers. Global-minimality labels are imported from
archive provenance; caps, physical laws, responses and class certificates
are replayed exactly. No asymptotic exact-minimizer obstruction is claimed.

## 1. Two sample operations preserve isotropy but lose the cheap overlap law

Let mu be a centered isotropic law on Boolean words. Distinct physical
coordinates are then pairwise independent fair signs. For independent
samples X^(1),...,X^(r), each of the following outputs is a physical,
centered, exactly isotropic sign vector:

- Coordinatewise product of the samples.
- Coordinatewise majority of three samples.

Indeed the two coordinate input strings are independent, and each
output coordinate is fair. This argument does not assert independence
of all output coordinates or a sharp subGaussian proxy.

Apply these operations to the uniform, independently globally centered
complete-ground law at orders12 and14. The maximum response over ALL
ground queries is exactly as follows.

| Physical law | Order12 | Order14 |
|---|---:|---:|
| One ground sample, optimal | 9/5 | 98/39 |
| Product of two | 74/25 | 1507/507 |
| Product of three | 1179/500 | 1382773/474552 |
| Majority of three | 4887/2000 | 5446063/1898208 |

Every displayed sample operation strictly worsens the optimal ground
response. At order12 the useful target is9/4: the original law beats
it, but all three displayed operations exceed it. Products of four
also fail, with exact value14749/5000 at order12.

At order14 there is a direct reason any operation generating off-code
columns must lose in the uniform dual average: the156 ground columns
are exactly the minimizers of that dual response among ALL physical
columns. The majority operation retains only mass1409/12168 on them.
At order12 majority retains mass7/40; the product of two retains none.

These are precise finite falsifiers of response preservation under
natural physical Gaussianizing/collision operations. They do not rule
out other nonlocal functions of ground samples, and in particular they
do not contradict the cheap original ground-sampling law itself.

## 2. A genuinely different physical matching class

Let n=2m. Choose ANY perfect matching of physical coordinates and signs
s_e in {+-1} on its edges. Draw one global fair pole sigma and independent
fair pair signs epsilon_e. For an oriented matched edge e=(i,j), set

    h_i=epsilon_e, h_j=sigma s_e epsilon_e.          (1)

The law is physical and centered. Matched correlations cancel after
averaging sigma; unmatched correlations vanish by independent pair
signs. Thus Ehh^T=I exactly. Conditional on sigma, the linear MGF is
bounded by exp(||theta||^2), so every such law is2-subGaussian.
Arbitrary mixtures of matching/sign choices retain these properties.

The matching and its signs may be selected using arbitrary ground
samples or signed-covariance data BEFORE drawing the fresh fair pole
and independent pair signs. Such adaptive selection is included in
the mixture class. Pole-dependent selection or dependence of the pair
signs is not automatically included.

For a fixed Boolean query x, let N_x be the number of matching edges
with s_e x_i x_j=+1 and let a_k=E|sum_(j<=k) epsilon_j|, a_0=0.
At one pole only these N_x pairs contribute, each with coefficient2;
at the other pole only the complementary m-N_x pairs contribute.
Consequently the response has the exact finite formula

    E|h.x|=a_(N_x)+a_(m-N_x).                      (2)

This class is not confined to the positive latent spectral band of
the low-cap Gaussian-pair theorem. Its best pointwise asymptotic floor
is kappa/sqrt(2), achieved when N_x=0 or m. Thus the narrow-band floor
above3/4 does not rule it out before the energy geometry is examined.

## 3. The complete matching-mixture game at order12 is exactly9/4

For the actual stored order12 ground code, enumerate every perfect
matching (11!!=10,395) and every edge-sign vector modulo simultaneous
reversal (2^5=32). This covers ALL332,640 laws of the form(1), because
global sign reversal of the edge signs only exchanges the two poles.

Let mu be uniform on the20 projective ground queries. Integer
evaluation of (2) gives the exact inequality

    E_(x~mu) E_matching |h.x| >=9/4                (3)

for EVERY signed matching. It therefore holds for every mixture of
such laws, including choices selected adaptively from ground samples.
There are105 minimizing signed matchings. Their uniform mixture has
response exactly9/4 at EVERY one of the20 ground queries. Hence

    min_(matching mixtures) max_(ground x) E|h.x|=9/4. (4)

One minimizing matching has pairs
(0,1),(2,4),(3,7),(5,11),(6,8),(9,10) and edge signs
(+,+,-,-,-,+). Its ground-response histogram has8 values15/8 and
12 values5/2, whose average is9/4. A single matching has best possible
maximum5/2; mixing is essential to attain(4).

By contrast the unrestricted ISOTROPIC physical game is exactly9/5,
attained by uniform ground sampling. Thus physical signs, exact
isotropy, and access to arbitrary signed perfect matchings still miss
25 percent of the optimal response on this actual example. The
matching-mixture optimum equals the proposed first-order target
3Q/(2n)=9/4, so it cannot produce a STRICT favorable slope here.

This is not an asymptotic failure theorem for the class. It is a
complete finite mechanism discriminator stronger than observing that
one particular matching chosen from K performs poorly. It also does
not rule out larger blocks or genuinely nonmatching ground laws.

The same exact enumeration closes matching-mixture games at order6
and the two stored order8 classes: respectively7/4,7/4,27/16.
Their corresponding unrestricted isotropic optima are5/3,7/4,3/2.

## 4. Replay and verification boundary

`computations/paper_bernoulli_2026_09_17_ground_sample_mechanisms.py`
computes exact XOR convolutions and directly enumerates independently
centered ternary majority outputs. It verifies all degree-one/two
Fourier coefficients giving exact isotropy and every ground-query
response. Independent centering of the three input samples is retained;
it is not replaced by merely centering the final output.

`computations/paper_bernoulli_2026_09_17_signed_matching_ground_game.py`
enumerates all matchings/signs with integers, verifies the dual minimum
in(3), and sums the response vectors of all105 minimizers to certify
the matching upper bound in(4). It uses no LP solver.

Outputs are `ground_sample_mechanisms.json` and
`signed_matching_ground_game.json` under
`tmp/paper_portfolio_2026_09_17/bernoulli/`.
These finite exhaustive computations are the certificates for the
displayed finite class optima; the general law and response formula
are proved above. No finite result is promoted to a statement about
arbitrarily large exact minimizers.

The discrepancy track independently reconstructed formula(2), the
matching recursion, pole-exchange reduction, complete dual minimum,
and uniform105-law primal, and replayed the script: PASS.

## 5. An asymptotic actual-ground obstruction for coordinate matchings

The preceding finite test has a stronger principled analogue, but on
an actual near-half family rather than unknown minimizing signings.
Let p=2^a>=8, n=p^2, and use the alternating Walsh matrix

    F_((u,v),(s,t))=(-1)^(u.t+v.s), A=F-I.

Then F^2=p^2 I and Q(A)=n(p+1)/2. Let G(A) be its COMPLETE absolute
ground code, and let L_match(n) be the convex mixture class of all physical
signed-coordinate-matching laws(1). The following sharp leading-order
statement holds:

    inf_(nu in L_match(n)) max_(x in G(A)) E_nu|h.x|/sqrt(n)
        =kappa+O(1/p).                              (5)

Proof. Sample a uniform fixed-point-free involution pi of the p column
labels and independent fair signs on its pairs, with the NEGATIVE
global relation z_(pi(v))=-z_v. Set

    X_(u,v)=(-1)^(u.pi(v)) z_v.

Exact summation gives FX=-pX, so this query law mu_- is supported on
EXACT absolute ground words, not merely the opposite eigensector at
deficit n. It is the negative-pole law already constructed in
[the nonlocal sign theorem](paper_localization_nonlocal_sign_response_2026_09_17.md).

The four-distinct-coordinate absolute-moment bound in
[the discrepancy track, Section17](paper_discrepancy_2026_09_17.md)
is unchanged at this negative pole. In the multiplicity cases4 and
2+2 the cycle signs cancel entirely. In3+1 and2+1+1 the surviving
pair contributes only the global sign -1. In1+1+1+1 two internal
pairs contribute its square. Thus the same bound

    c_p=(2p-1)/[(p-1)(p-3)]

holds except for four coordinates in one physical column whose row
labels xor to zero; those exceptional moments have magnitude one.
The same exceptional-edge-pair counting therefore gives, uniformly
for EVERY signed matching of m physical-coordinate pairs,

    E_(mu_-) S_M^2<=c_p m^2+(1-c_p)mp/2,
    S_M=sum_(ij in M) s_ij X_i X_j.                 (6)

For a full matching m=n/2, put v=S_M/m. Equation(6) yields

    E v^2<=d_p:=c_p+(1-c_p)/p=O(1/p).              (7)

By(2) and the uniform scalar absolute-moment estimate
a_k=kappa sqrt(k)+O(1), including k=0,

    E_(matching law)|h.X|/sqrt(n)
       =kappa f(v)+O(1/p),
    f(v)=(sqrt(1+v)+sqrt(1-v))/2.

The error is uniform even when one branch has zero active pairs.
For instance the already proved
[shifted-absolute estimate](paper_bernoulli_shifted_absolute_2026_09_17.md)
gives a safe absolute error at most6/p after summing the two branches.
Since f(v)>=1-v^2/2 on[-1,1], averaging and using(7) gives the explicit
uniform lower bound

    E_(mu_-) E_(matching law)|h.X|/sqrt(n)
       >=kappa-(kappa/2)d_p-6/p.                   (8)

It remains valid for EVERY mixture of coordinate matchings, even when
the mixture is selected after examining A and its whole ground code.
This dual law proves the lower half of(5). For the upper half, fix
one matching and randomize all its edge labels independently and
fairly in(1). The resulting physical law is exactly iid fair signs,
whose response is a_n/sqrt(n)=kappa+o(1/p), uniformly in x.

This is a class obstruction on actual complete ground codes, not an
arbitrary-code surrogate. It does NOT assert that A is asymptotically
minimizing: its normalized cap tends to1/2. Nor does it exclude all
physical matching-based laws. The successful Hadamard-column
involution law matches LATENT column labels and creates a nonlocal
physical dependence; it is not a matching of the physical coordinates.
That law has response at most kappa/sqrt(2)+o(1) on this ground sector.
Thus the distinction is quantitatively leading-order, not vocabulary.

The discrepancy track independently reconstructed the negative-pole
fourth-moment extension, all-matching estimate, uniform scalar-error
step, and iid upper law: PASS. The separate reproduction
`computations/paper_bernoulli_2026_09_17_coordinate_matching_hadamard.py`
enumerates all1,680 negative-pole query atoms at p=8, verifies their
exact ground energy and covariance(pI-F)/(p-1), checks10,560 fourth
moments and128 signed full coordinate matchings. Its output is
`tmp/paper_portfolio_2026_09_17/bernoulli/coordinate_matching_hadamard.json`.
Those tests check the normalization, not an exhaustive test of all
coordinate matchings at order64. The uniform statement is proved by(6).

## 6. Why coordinatewise sample Gaussianization loses marked covariance

There is also an exact classical Fourier obstruction to repairing the
failed majority/product operations while preserving their conditional
covariance. Suppose X is a globally sign-symmetric physical law and
rho=E X_i X_j satisfies0<|rho|<1. Apply the SAME odd Boolean map
F:{+-1}^r->{+-1} to r independent samples at each coordinate. The
output correlation is

    g_F(rho)=sum_(d odd) W_d(F) rho^d,
    W_d(F)=sum_(|S|=d) Fhat(S)^2,
    sum_d W_d(F)=1.                                (9)

This is the usual noise-stability Fourier formula, stated in the
primary [Mossel--O'Donnell--Oleszkiewicz paper](https://www.cs.cmu.edu/~odonnell/papers/invariance.pdf),
Section1.1. Here it applies because the two coordinate sample strings
consist of independent fair correlated pairs; no independence of all
physical coordinates is assumed. Orthogonality directly proves(9).

For0<|rho|<1, equality g_F(rho)=rho forces W_1(F)=1. A Boolean
linear function must be a signed dictator: expanding F^2=1 forces
every pair of distinct linear coefficients to have product zero.
Thus the only exact covariance-preserving common odd recombination
is copying a sample, possibly with global sign reversal.

The conclusion extends to coordinate-dependent odd maps F_i. For an
edge with nonzero nonsaturated correlation rho_ij,

    E F_i(U)F_j(V)
      =sum_(d odd) rho_ij^d <Fhat_i^(d),Fhat_j^(d)>.

Cauchy--Schwarz bounds its magnitude by |rho_ij|. Preserving the
original signed value forces both maps to be the SAME signed dictator.
Therefore on each connected component of the graph
0<|rho_ij|<1, every coordinate must copy the same input sample with
the same sign. Independent random mixtures of processing rules cannot
evade this equality condition, because each rule already has the
corresponding one-sided correlation bound.

For odd majority on r inputs,

    W_1=r[binom(r-1,(r-1)/2)/2^(r-1)]^2 ->2/pi.

Thus for small conditional correlations, majority shrinks their
leading value by precisely the familiar Gaussian-angle factor. Odd
products replace rho by rho^r. Equal mixtures of opposite conditional
sector covariances can remain EXACTLY isotropic throughout, concealing
this loss of signed sector information. The finite tests in Section1
show why isotropy alone does not certify preservation of a useful
absolute-response profile.

This is a scoped application of classical Fourier/maximal-correlation
algebra, not a new theorem forbidding arbitrary nonlocal algorithms.
It does not cover coordinate rules that examine other coordinates,
global energy conditioning, or the successful latent-index Hadamard
matching law. Those operations are genuinely outside the restriction
used to derive(9).
