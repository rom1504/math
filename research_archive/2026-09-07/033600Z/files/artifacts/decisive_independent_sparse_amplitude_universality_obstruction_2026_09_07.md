# Adaptive sign minimization is not governed by variance-only universality

Status: proved, with a precise triangular-array scope. This does NOT settle
the fixed Gaussian-magnitude model and does NOT obstruct a one-sided
flatification inequality in the favorable direction.

For a real symmetric hollow array g let

    T_n(g)= n^(-3/2) min_(A_e in {+1,-1})
                       max_x |sum_e g_e A_e x_i x_j|.

Only the magnitudes of g matter. In particular, for iid fair sign g,
T_n(g)=M_n/n^(3/2) deterministically.

## 1. A signed high-girth Gaussian-wave lower bound

Let G be a d-regular simple graph on n vertices, d>=2, with girth
greater than 2r+3, where r>=2. For EVERY signing A of its edges,

    max_x sum_(ij in E(G)) A_ij x_i x_j
      >= (2/pi) n sqrt(d-1) cos(pi/(r+1)).             (1)

Thus (1) also bounds the absolute cap from below.

Here is an explicit proof. Let s_0=1 and
s_k=d(d-1)^(k-1) for k>=1. Put a_0=0 and

    a_k=sqrt(2/(r+1)) sin(k pi/(r+1)), 1<=k<=r.

Then sum a_k^2=1 and
sum_(k=1)^(r-1) a_k a_(k+1)=cos(pi/(r+1)). For independent standard
Gaussians (g_v), set

    Z_i=sum_(k=1)^r a_k/sqrt(s_k)
                    sum_(dist(i,v)=k) A(path(i,v)) g_v.

All indicated shortest paths are unique. Var(Z_i)=1. For adjacent i,j,
the two branches of their joint tree neighborhood give EXACTLY

    Cov(Z_i,Z_j)
      = A_ij [2sqrt(d-1)/d] cos(pi/(r+1)).             (2)

Indeed, on either side of the edge the number of vertices at distances
k,k+1 from its endpoints is (d-1)^k, for k>=1. Dividing by
sqrt(s_k s_(k+1)) gives sqrt(d-1)/d. The products of the two path
signs equal A_ij, since their symmetric difference is that edge. There
is no k=0 contribution because a_0=0.

Round x_i=sign(Z_i). The Gaussian sign identity gives

    E[A_ij x_i x_j]=(2/pi) arcsin(rho),
    rho=[2sqrt(d-1)/d] cos(pi/(r+1))>=0.

Since arcsin(rho)>=rho, summing nd/2 edges proves (1). The randomness
is used only to exhibit a spin configuration for each signing; every
constant and every structural condition is independent of that signing.

## 2. IID prescribed magnitudes: a genuine Lindeberg-scale counterexample

Let d_n=(log n)^3, p_n=d_n/n (for sufficiently large n), and let

    g_ij = epsilon_ij xi_ij / sqrt(p_n),

where the epsilon are independent fair signs and the xi are independent
Bernoulli(p_n), all mutually independent. Then E g_ij=0,
E g_ij^2=1, and

    max_ij |g_ij|/sqrt(n) <= 1/sqrt(d_n) -> 0.

Nevertheless,

    for every eta>0,
    P(T_n(g)>=2/pi-eta) -> 1.                         (3)

Proof: fix r>=2 before taking n to infinity. The support graph is
G(n,p_n). With probability tending to one, all degrees belong to
[(1-epsilon_n)d_n,(1+epsilon_n)d_n], where
epsilon_n=(log n)^(-1/2). This follows directly from Chernoff and the
union bound, since epsilon_n^2 d_n=(log n)^2.

Call a vertex bad if its radius-(r+1) ball contains a cycle. For fixed r,
the expected number of bad vertices is bounded by a polynomial in d_n.
One elementary count enumerates a rooted breadth-first-search collision:
its connected witness has at most 2r+3 edges and equally many vertices,
so its expected number of embeddings is at most a constant depending
on r times sum_(ell<=2r+3) d_n^ell. Consequently the number of bad
vertices is at most n^(1/4) with probability tending to one. Together
with the degree bound, only o(nd_n) edges have a bad endpoint.

For a good vertex i define instead

    Y_i=sum_(k=1)^r a_k d_n^(-k/2)
                     sum_(dist(i,v)=k) A(path(i,v)) g_v.

On its tree ball, the degree bounds imply, uniformly in i and in ALL
edge signings,

    Var(Y_i)=1+o(1).

For an edge with both endpoints good, the common tree-ball calculation
from (2) gives, uniformly,

    Corr(Y_i,Y_j)
      = A_ij (2/sqrt(d_n))
                     [cos(pi/(r+1))+o(1)].            (4)

The error here is relative at the 1/sqrt(d_n) scale: all finite-depth
branch counts are (1+O_r(epsilon_n+1/d_n)) times d_n^k.

At every bad vertex use a fresh independent Gaussian, independent of
all the g_v and of the other fresh Gaussians. Its rounded spin has
zero correlation with every other rounded spin, so edges with a bad
endpoint have expected signed energy ZERO, not an uncontrolled negative
contribution. On good-good edges use sign(Y_i). Equations (4), the
Gaussian sign identity, and |E_good|=(1+o(1))nd_n/2 imply

    max_x sum_(ij in E(G)) A_ij x_i x_j
      >= (2/pi-o(1)) n sqrt(d_n) cos(pi/(r+1)).        (5)

The support-graph event suffices for (5) simultaneously for every
signing. Multiplication by p_n^(-1/2), division by n^(3/2), and then
letting the fixed integer r tend to infinity prove (3).

## 3. What universality assertion is falsified

The iid fair-sign amplitude model has T_n=M_n/n^(3/2), with the campaign's
rigorous all-order limsup below .494515125. Equation (3), whose constant
2/pi is approximately .636619772, proves that adaptive optimal signing
cannot have a universal asymptotic value depending only on entry variance,
even for independent centered variance-one triangular arrays whose largest
normalized interaction tends uniformly to zero.

This example also meets the third-moment scale used in the ordinary
finite-temperature SK Lindeberg replacement:

    E|g_ij|^3/sqrt(n) = 1/sqrt(d_n) -> 0.

Indeed the usual coordinate third-derivative bound for normalized
fixed-sign spin pressure gives replacement error
O_beta((E|g_ij|^3+E|G|^3)/sqrt(n))=o(1). That smooth fixed-sign-pressure
argument is legitimate here. It cannot be carried through an adaptive
minimum over all coefficient signs without a new estimate.

The conclusion is NOT that the fixed standard-Gaussian amplitude model
has value at least 2/pi. Its support is dense, and the locally tree-like
argument does not apply. Nor does this refute the potentially useful
ONE-SIDED construction bound

    M_n/n^(3/2) <= T_n(g)+o(1).

Our sparse model has the larger value, which is compatible with that
inequality. It refutes equality/two-sided replacement, not favorable
flatification from an already good weighted model.

Finally, the absolute value (or an auxiliary global orientation spin)
must be retained in the original minimax model. The literal flat
one-sided minimum min_A max_x sum A_ij x_i x_j is at most n/2 by taking
all signs negative, and therefore has normalized limit zero. The sparse
lower bound above holds even for its positive maximum; omitting the
absolute value changes the comparison substantially.

## 4. The natural Gordon block comparison fails its covariance hypothesis

This is a separate exact obstruction, not a proof that every possible
Gaussian comparison fails. It already occurs for the POSITIVE-maximum
subproblem (before adding the second orientation). Take an even N and balanced shores, and write
z_i=+1 on one shore and -1 on the other. Compare the Gaussian fields

    X_(A,x)=N^(-1/2) sum_e g_e A_e x_i x_j,
    Y_(A,x)=sqrt(2/N) sum_(e within shores) h_e A_e x_i x_j + gamma/sqrt(2),

where the g,h,gamma are independent standard Gaussians. The common scalar
gamma makes all point variances equal and has zero effect on the expected
min-max value. If R_e=A_e B_e and u_i=x_i y_i, then

    Cov(Y_(A,x),Y_(B,y))-Cov(X_(A,x),X_(B,y))
      = 1/2 + (1/N) sum_(i<j) z_i z_j R_ij u_i u_j.   (6)

For the SAME outer signing A=B, this equals

    (z dot u)^2/(2N) >= 0.                           (7)

This is the usual favorable squared-overlap term. But Gordon's min-max
comparison requires the OPPOSITE covariance ordering when A and B are
different outer indices. Let B differ from A on just one edge and choose
u=z. Then (6) equals

    N/2 - 2/N > 0.

It has the same sign as (7) and macroscopic size. Thus the required
cross-outer-index inequality is false. The example survives switching
and global-orientation quotienting: at N>=4 a single-edge reversal is
neither a cut nor the complement of a cut, so it does not identify the
two outer switching classes. The case N=4 can be checked directly and
the same statement for N>=5 follows from cut cardinalities.

For the absolute problem one must in addition include s=+1,-1 in the
maximizing index. Using Y_(A,s,x)=s(H_block(x)+gamma/sqrt(2)) changes
the expected block absolute cap by at most E|gamma|/sqrt(2), but the
same-outer-index covariance difference is now
s t (z dot u)^2/(2N). It already has both signs when the two global
orientations differ. The one-sided calculation above is therefore not
being used to assert a valid two-sided within-row covariance ordering.

For source verification, Theorem 1 of Akhtiamov et al.,
https://arxiv.org/html/2402.07356v3 , explicitly lists equal variances,
one covariance ordering for a fixed minimizing index, and the reverse
ordering for distinct minimizing indices. Its two-sided CGMT statements
add convexity of the index sets and convex-concavity of the payoff.
The minimizing signing set here is not convex; replacing it by its cube
convex hull includes A=0 and changes the absolute-cap optimum to zero.

The recent finite-third-moment min-max comparison of Chen, Chen, and Wei,
https://arxiv.org/html/2411.08303v1 , does not remove the size issue:
its smoothing lemma has errors log(I)/(beta delta) and log(J)/beta,
where here log(I)=Theta(N^2), while log(J)=Theta(N). Its main theorem
also retains explicit moment/error terms. Neither a dimension-free
adaptive universality statement nor an all-order limit is imported.
