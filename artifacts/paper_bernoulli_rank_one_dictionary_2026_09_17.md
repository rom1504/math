# Sharp asymptotic mean response of the whole Boolean rank-one dictionary

2026-09-17. Bernoulli-track derivation. This is a finite-code response
theorem, a concrete larger conditional extension class, and a precise
limitation of low entropy alone. It does not assert that a minimizing
quadratic signing has this near-level geometry. External novelty is
not established.

## 1. A support-free asymptotic minimax theorem

Let n=kp and

```
F_{k,p}={c tensor v: c in {-1,1}^k, v in {-1,1}^p}.
```

This code has cardinality 2^(k+p-1). For any probability law nu on
physical full-sign columns z in {-1,1}^n define

```
R(nu)=max_{x in F_{k,p}} E_nu |z dot x|/sqrt(n).
```

No isotropy, symmetry, or support restriction is imposed. Write
`a_d=E|epsilon_1+...+epsilon_d|` and kappa=sqrt(2/pi). Then

```
2/pi-sqrt(2log(2k)/k)-sqrt(2log(2p)/p)
 <= inf_nu R(nu) <= a_k a_p/sqrt(kp).                 (1)
```

The upper bound is attained as a response value by the law
z=h tensor g with h and g independent uniform sign vectors. This law
is centered and has covariance I_n, and EVERY x in the dictionary has
the same mean response a_k a_p. Therefore, if k,p tend to infinity,

```
inf_nu R(nu) -> 2/pi.                                (2)
```

In particular exp(o(n)) cardinality does not imply vanishing mean
absolute response under a suitably selected sign law, even when all
sign laws are allowed. This is not merely a failure of one isotropic
rounding algorithm.

### Proof of the lower bound

Identify a deterministic sign column z with its k by p matrix Z.
Average the response over independent uniform c and v. The scalar
Stein bound, reconstructed in Section 28 of
[the Bernoulli record](paper_bernoulli_2026_09_17.md), gives for arbitrary
real coefficients b_i

```
d_W(sum_i b_i epsilon_i, sqrt(sum_i b_i^2) G)
 <= sum_i |b_i|^3/sum_i b_i^2 <=max_i|b_i|.           (3)
```

The zero-variance case costs zero. Conditional on v, apply (3) to the
coefficients Zv, replacing c by a standard Gaussian vector g. Since
each coordinate of Zv is a p-sign sum, the exponential-moment maximum
bound yields

```
|E|c^T Z v|-E|g^T Z v||
 <= E max_i |(Zv)_i| <=sqrt(2p log(2k)).              (4)
```

Now condition on g and replace v by a standard Gaussian vector h.
Each coordinate of Z^Tg is a centered Gaussian of variance k, so

```
|E|g^T Z v|-E|g^T Z h|| <=sqrt(2k log(2p)).          (5)
```

No independence among the row maxima or column maxima is needed.
For singular values s_i of Z,

```
E|g^T Z h|=kappa E sqrt(sum_i s_i^2 h_i^2)
         >=kappa^2 sqrt(sum_i s_i^2)= (2/pi)sqrt(kp).
```

Indeed the expected square root is a concave function of the
nonnegative squared singular values with fixed sum; its minimum over
that simplex is attained at a vertex, where it equals kappa times the
square root of the sum. Combine (4)--(5), then average over any nu.
The maximum over dictionary words is at least their uniform average.

For the upper bound, `(h tensor g) dot (c tensor v)=(h dot c)(g dot v)`.
The factors are independent sign sums, proving their common absolute
mean a_k a_p. Covariance is I_k tensor I_p. Finally a_d/sqrt(d) tends
to kappa by the ordinary bounded-variable central limit theorem and
uniform integrability, or by the scalar Wasserstein estimate (3).

## 2. One actual empirical frame protects the entire dictionary

Let q grow with k, suppose k log(4k)=o(q), and take q independent
uniform sign labels h_j in {-1,1}^k. There is one deterministic list
such that

```
sum_j h_j h_j^T <=q(1+o(1))I_k,
sup_c q^(-1)sum_j |h_j dot c| <=a_k+O(k/sqrt(q)).     (6)
```

The first event is the previously reconstructed rank-one matrix
Chernoff estimate. For each c the second sum has mean q a_k; changing
one of its qk independent sign bits changes it by at most two.
Its centered MGF proxy is qk. A union over the 2^k choices of c gives
an O(k sqrt(q)) deviation with fixed positive failure probability.
Choose the constants so each event fails with probability at most 1/4.

Build the physical bridge column j as h_j tensor g_j, where the g_j
are independent fair-sign vectors of length p. For every rank-one old
word c tensor v its mean column response is a_p |h_j dot c|. Thus if
k,p tend to infinity and k=o(q), the center contribution is uniformly

```
sum_j E|column_j dot (c tensor v)|
 <= [2/pi+o(1)]q sqrt(n).                            (7)
```

The same empirical list has aggregate subGaussian bridge proxy
q(1+o(1))||w||_2^2 for ANY residual vector w and ANY fixed new-spin
assignment. The inner g_j signs remain independent. Both assertions
are uniform in the whole dictionary, whose logarithmic cardinality is
o(n) when k,p tend to infinity.

Consequently the already proved anchored Hamming-cover shell compiler
applies with center coefficient 2/pi rather than zero or
kappa/sqrt(k). In its notation

```
K_H=limsup_{eta down to 0} r(eta) h(r(eta))/eta,
```

and a full near-level cover by this one dictionary would suffice for
an extension slope strictly larger than

```
2/pi+2K_H.                                          (8)
```

This is a reference to the existing complete shell proof, not a new
independent all-energy argument. Take k,p tending to infinity with
k log k=o(n) before the small extension-ratio limit, and retain the
same dictionary across its finitely many hierarchy levels. Leftover
coordinates for n not divisible by k may be included with all their
sign patterns and independent physical bridge signs: their entropy
is o(n) and their mean response is at most q sqrt(ell)=o(n^(3/2))
when ell<k=o(n). Thus arithmetic divisibility is not essential.

At the currently reported lower constant c_*=0.4333221116640807,
2/pi<3c_*/2. Therefore a zero-K_H whole-nearcode rank-one cover would
contradict asymptotic minimality by the existing extension argument.
This use of c_* is conditional on the project's established lower
bound. The elementary kappa/4 lower bound alone is not strong enough.

The dictionary is substantially larger than the k specified Hadamard
modes: it contains ALL 2^(k-1) sign directions c. Conversely, the
previous arbitrary-rank-one Hadamard ground-code obstruction shows
that this premise is not automatic even on the known near-half family.
Neither that family nor the code-only minimax theorem settles its
applicability to actual asymptotically minimizing signings.

## 3. Exact Gaussian child-value control in strong independent fields

This adjacent theorem identifies what a Gaussian branch-value proof
must pay. Let D be a real symmetric hollow q by q matrix and s a
polarity. Let

```
Psi_{D,s}(sigma)=E max_y [sH_D(y)+sum_i sigma_i G_i y_i],
L=lambda_max(sD)>=0,  r_i^2=sum_j D_ij^2.
```

For strictly positive sigma_i,

```
kappa sum_i sigma_i
 <= Psi_{D,s}(sigma)
 <= kappa sum_i sigma_i+kappa sum_i (L^2+r_i^2)/sigma_i.  (9)
```

For arbitrary nonnegative sigma_i and any a>0,

```
Psi_{D,s}(sigma)
 <=kappa sum_i sigma_i+kappa sum_i(a-sigma_i)_+
     +kappa sum_i (L^2+r_i^2)/max(sigma_i,a).          (10)
```

Proof: put z_i=sign(G_i), independent of the half-normal magnitudes.
Expand the quadratic at z. For any y, let S be its flipped coordinates.
The spectral upper bound gives

```
sH_D(y)-sH_D(z)<=2sum_{i in S}[L-s z_i(Dz)_i].
```

The linear fields lose 2sum_{i in S}sigma_i|G_i|. Optimizing over S
therefore adds at most
`2sum_i[L-s z_i(Dz)_i-sigma_i|G_i|]_+` above the z energy. For A>=0,
the half-normal density is at most kappa, so

```
E(A-sigma|G|)_+<=kappa A^2/(2sigma).
```

Moreover E_z[L-s z_i(Dz)_i]_+^2<=L^2+r_i^2; the row response has
mean zero and second moment r_i^2. The baseline expected child energy
is zero. This proves (9), including its lower bound by using z itself.
Increasing an independent centered Gaussian coordinate's variance can
only increase the expected convex maximum. Replace sigma_i by
max(sigma_i,a) and apply (9) to obtain (10).

In particular full sign D has r_i^2=q-1. Children with
||D||_op<=sqrt(2q)+1 exist at every q: restrict the next larger
Sylvester Hadamard matrix to q coordinates and remove its diagonal.
For flat fields sigma_i=sqrt(n), (9) gives

```
Psi_{D,s}=kappa q sqrt(n)+O(q^2/sqrt(n)).             (11)
```

The leading flat-field response kappa q sqrt(n) is too large for the
required derivative 3c/2 throughout the reported c<=0.493608094.
Small old-code entropy and accurate Gaussian replacement cannot
remove that VALUE obstruction. No assertion that all minimizing
near-level words have flat variance profiles is made.

For coherent rank-one profiles sigma_i/sqrt(n) distributed as |G|,
floor at a=sqrt(q), with epsilon=q/n<=1. The elementary bounds

```
E(t-|G|)_+<=kappa t^2/2,
E[1/max(|G|,t)]<=kappa[2+log(1/t)]  (0<t<=1)
```

inserted in (10) instead give

```
Psi_{D,s}/n^(3/2)
 <=(2/pi)epsilon+O(epsilon^2 log(e/epsilon)).         (12)
```

The empirical uniform-sign mode construction realizes the needed
profile simultaneously in c: first fix epsilon>0, then let k,p,q
grow as above. The normalized functions |u|, (sqrt(epsilon)-|u|)_+,
and 1/max(|u|,sqrt(epsilon)) have bounded Lipschitz constants for
fixed epsilon. Scalar Stein, bounded differences, and the same 2^k
union make their empirical averages converge uniformly to the stated
Gaussian averages. The small-epsilon limit is taken only afterward.

Formula (12) controls the exact child optimization, not merely its
separately paid cap. Nevertheless the direct physical center estimate
(7) already proves the corresponding anchored-cover first-order
coefficient. Neither result supplies the missing minimizing-family
cover or excludes arbitrary outside-window old states.

## 4. Weighted independent-block hierarchy

The director proposed the following strengthening. Partition the N
physical coordinates into rectangular blocks of dimensions k_a by p_a,
a=1,...,r, with n_a=k_a p_a and sum_a n_a=N. Let C be the code whose
restriction to each block is an arbitrary Boolean rank-one array,
independently of all other blocks. Define

```
lambda_a=n_a/N,
eta=sum_a n_a(k_a+p_a)/N^2,
beta(lambda)=kappa E sqrt(sum_a lambda_a G_a^2).       (13)
```

Then the support-free minimax response obeys the finite bounds

```
beta(lambda)-8 eta^(1/4)
 <=inf_nu max_{x in C} E_nu|z dot x|/sqrt(N)
 <=beta(lambda)+8 eta^(1/4).                          (14)
```

The infimum again ranges over ALL physical sign-column laws. The
upper construction is independent uniform rank-one sign arrays in
each block; it is centered, isotropic, and has a response independent
of the selected codeword. Thus (14) is asymptotically sharp whenever
eta tends to zero, including growing numbers of blocks. No minimum
block-dimension assumption is required for this formulation.

### Uniform fourth-order replacement

For an arbitrary physical sign array Z, write its uniform-code query as

```
Y=N^(-1/2)sum_a c_a^T Z_a v_a.
```

Replace all independent sign entries of c and v successively by
independent standard Gaussians. Conditional on the other variables,
Y is affine in the driver currently replaced, and its coefficient A
has fourth moment at most 3p_a^2/N^2 for a c-driver, or at most
3k_a^2/N^2 for a v-driver. The opposite-side variables are independent
centered signs/Gaussians with fourth moments at most three, so this
bound holds at every replacement stage. Summing gives

```
sum_drivers E A^4 <=3 eta.
```

Smooth absolute value by f_rho(t)=sqrt(t^2+rho^2). It differs from
|t| by at most rho, and direct differentiation gives

```
f_rho''''(t)=3rho^2(4t^2-rho^2)/(t^2+rho^2)^(7/2),
||f_rho''''||_infinity<=12/rho^3.
```

The first three sign/Gaussian moments agree. Taylor's fourth-order
remainder, using E epsilon^4+E G^4=4, bounds each replacement by
`||f''''|| E A^4/6`. The total error is at most 6eta/rho^3.
Returning to absolute value adds at most 2rho. Choose rho=eta^(1/4)
to obtain the stated 8eta^(1/4) error, uniformly in Z and the number
and aspect ratios of its blocks.

### Gaussian minimization under the block constraints

After replacement, condition on all Gaussian right-side vectors.
The normalized mean absolute query is

```
kappa E sqrt(sum_{a,l} s_{a,l}^2 G_{a,l}^2/N),
```

where the squared singular values in block a sum to n_a. The
expectation is concave in all squared singular values. Minimize over
the product of their fixed-sum simplices: a minimizing vertex has
one nonzero squared singular value n_a in EACH block. Its value is
exactly beta(lambda). Average the original response over the code,
then over any physical column law, to get the lower bound in (14).

For the upper bound, the independent block-rank-one column law has
the same response at every codeword, namely the absolute mean of
`N^(-1/2)sum_a S_(k_a) S_(p_a)`, with all factors independent.
Apply the same replacement bound to rank-one sign blocks. Their
Gaussian value is beta(lambda), giving the upper side of (14).

The code's exact cardinality is

```
|C|=2^(sum_a(k_a+p_a-1)).                            (15)
```

In particular subexponential code size implies eta=o(1):
`k_a+p_a<=2(k_a+p_a-1)` implies sum_a(k_a+p_a)=o(N), and
eta<=max_a(k_a+p_a)/N<=sum_a(k_a+p_a)/N. (When r is comparable
to N, this entropy premise fails; (14) itself still applies.)

### Equal blocks and the sharp structural transition

For r equal-area blocks, beta(lambda)=beta_r, where

```
beta_r=kappa E chi_r/sqrt(r)
      =2 Gamma((r+1)/2)/(sqrt(pi r) Gamma(r/2)).       (16)
```

In particular

```
beta_1=2/pi,  beta_2=1/sqrt(2),
beta_3=4/(pi sqrt(3)),  beta_4=3/4,
beta_r increases to kappa as r tends to infinity.
```

Monotonicity follows from concavity of square root and the
leave-one-out representation of an average of r independent squares
as the average of its r leave-one-out (r-1)-averages. The limit follows
from the law of large numbers and uniform integrability. More
generally beta(lambda) is symmetric and concave in lambda; concentrating
all variance in one block minimizes it, while max_a lambda_a tending
to zero forces beta(lambda) to tend to kappa, by the variance bound
`Var(sum_a lambda_a G_a^2)=2sum_a lambda_a^2`.

Thus the first-order mean-response compiler cannot make four or more
equal independent blocks cheaper than 3/4. This exceeds 3c/2
throughout the reported c<0.493608094. One block is cheap enough
relative to the reported lower bound; two and three blocks have
thresholds c>sqrt(2)/3 and c>8/(3pi sqrt(3)), respectively.
These are sharp limitations of this uniform mean-response obligation,
not impossibility theorems for cancellation-sensitive joint mechanisms.

### Joint block response has an actual simultaneous empirical realization

Set K=sum_a k_a and P=sum_a p_a. Sample for each column j and block a
an independent uniform label h_(j,a) in {-1,1}^(k_a). Conditional on
these labels, the physical column in block a is h_(j,a) tensor g_(j,a),
with all inner signs independent. Put

```
delta=max_a sqrt(3k_a log(4r k_a)/q),
q>=3 max_a k_a log(4r k_a).
```

There is one deterministic list of all labels such that every block
Gram is at most q(1+delta)I and, simultaneously for every codeword x,

```
q^(-1)sum_j E_inner |column_j dot x|/sqrt(N)
 <=m_0+sqrt(2[Klog2+log4]/q),                        (17)
```

where m_0 is the common normalized response of the independent uniform
block-rank-one column law and |m_0-beta(lambda)|<=8eta^(1/4).

For the Gram claim apply matrix Chernoff in block a with failure
budget 1/(4r), then union over a. For (17), the inner expectation
depends on c_a but not v_a, because the latter can be absorbed into
the independent inner signs. Changing one label bit in block a
changes this normalized expectation by at most
`2a_(p_a)/sqrt(N)<=2sqrt(p_a/N)`. The sum over q columns therefore
has centered MGF proxy q, since sum_a k_a p_a=N. Union over the 2^K
choices of c and use the displayed deviation. Its failure is at most
1/4; both events thus intersect with probability at least 1/2.

For any residual vector w and any new-spin assignment, the independent
inner bridge signs have aggregate proxy at most
q(1+delta)||w||_2^2. Independence between the blocks makes their
covariance contributions block-diagonal; no cross-channel response
has been paid separately. The random sum of absolute column responses
at a center has this same proxy, by bounded differences. Union over
the full center dictionary costs only log|C|=(K+P-r)log2.

Consequently, whenever q is a fixed positive fraction of N,
K+P=o(N), delta=o(1), the entire dictionary has center coefficient
beta(lambda)+o(1) and sharp aggregate residual proxy. The existing
anchored-cover compiler then achieves the slope beta(lambda)+2K_H.
Fixed r and balanced growing block dimensions satisfy these conditions;
many growing-r balanced configurations do as well. This is a genuinely
joint response: separate payments would instead produce
`(2/pi)sum_a sqrt(lambda_a)`, which is not the value in (13).

## 5. Actual quadratic-energy scope: a whole dictionary is an overcover

The discrepancy track identified and the Bernoulli track independently
checked an important additional limitation. Under the uniform law on
the block-rank-one dictionary, form the quadratic feature vector
z(x)=(x_i x_j)_(i<j). Then

```
||E z(x)z(x)^T||_op <=max_a{k_a,p_a,2}.              (18)
```

Each edge feature is a character of the underlying independent c and
v signs. Within a block, a same-row pair has character v_j v_l and
multiplicity k_a; a same-column pair has multiplicity p_a; a pair
using distinct rows and columns has multiplicity two (the opposite
rectangle diagonals). Cross-block edge characters are unique. Distinct
characters are orthogonal, so the covariance is a direct sum of
all-ones matrices whose sizes are these multiplicities. This proves
(18), including blocks with dimension one.

If the ENTIRE dictionary were contained in a full signing's absolute
near-level code at deficit T<Q(A), averaging H_A(x)^2 would give

```
(Q(A)-T)^2 <=binom(N,2) max_a{k_a,p_a,2}.             (19)
```

For a subexponential block dictionary, the right side is o(N^3).
Thus the whole dictionary cannot itself be a shrinking-relative-window
near-ground code for a positive normalized cap. Its mean-response
barrier concerns the obligation to protect every center in a chosen
OVERCOVER; some of those centers necessarily have low original energy.
The sufficient Hamming-cover premise remains meaningful, because its
centers need not all be near-ground words. However restricting to just
the used/near-energy centers can invalidate the whole-dictionary lower
barrier and may permit cheaper response. No actual-minimizer obstruction
is inferred merely from beta_4=3/4.

## 6. Audit boundary

All proofs above are self-contained modulo the already reconstructed
scalar Stein, matrix Chernoff, and anchored-shell lemmas. The Gaussian
value theorem concerns arbitrary deterministic children with the
specified spectral quantity; it does not import the Parisi formula
for random SK interactions into a deterministic extremal problem.
No claim of an exact finite k,p minimax value is made beyond (1) and
(14). An initial cross-aspect audit concern about (4)--(5) was retracted:
after normalization their errors are sqrt(log k/k) and sqrt(log p/p),
not cross-denominator expressions. The two-block replacement proof
does not need an aspect-ratio restriction. Formula (14) independently
handles arbitrary block shapes and growing r with its explicit eta.

The discrepancy track independently reconstructed Sections 1, 3, and
4 and returned PASS, including all variance normalizations, the
spectral positive-part proof, fourth-moment replacement, weighted
product-simplex minimization, and growing-block scope. Its small
entropy-to-eta clarification has been included above.

Reproducibility: `computations/paper_bernoulli_2026_09_17_rank_one_dictionary.py`
enumerated 4,716 physical row/column-switching classes across nine
finite dictionaries (including heterogeneous multi-block examples).
All exact minima in these small cases were attained by block-rank-one
columns; this is diagnostic evidence only, not a finite theorem.
The same replay passed 960 pointwise spectral maximum certificates
and 480 exact-sign-law row-moment checks. It passed Python compilation.
The output is preserved at
`tmp/paper_portfolio_2026_09_17/bernoulli/rank_one_dictionary_audit.json`.
