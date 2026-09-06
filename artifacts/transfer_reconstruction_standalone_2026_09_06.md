# Standalone reconstruction of the strict all-order upper bound

Date: 2026-09-06. This proof reconstructs the mathematical steps; prior audit
verdicts are not premises. It gives the same upper bound as the previous
construction, with a shorter analytic closure and an elementary order supply.
The original convergence question is not resolved.

## 1. Statement and conventions

For a hollow symmetric sign matrix `A` of order `n`, write

```math
Q(A)=\frac12\max_{x\in\{-1,1\}^n}|x^TAx|,\qquad
M_n=\min_A Q(A).
```

Set

```math
p=\frac{31}{32},\quad t=4,\quad
a=\frac{91470529542342299}{20460000000000000000}.
```

Then

```math
\limsup_{n\to\infty}\frac{M_n}{n^{3/2}}
\le\frac12-\frac{a}{8\sqrt p}<0.499432220485404.
\tag{1}
```

All logarithms are natural. All information is relative entropy. For a
finite-second-moment probability law `nu` on Euclidean space define

```math
F_t(\nu)=\inf_{\gamma\in\Pi(\nu,\nu)}
\{D(\gamma\Vert\nu\otimes\nu)+t\mathbb E_\gamma|X-Y|^2\},
\qquad \Phi_t=-F_t/2.
\tag{2}
```

Its centered one-dimensional Gaussian value is `g_t(v)=Phi_t(N(0,v))`:

```math
g_t(v)=-tv(1-\rho)+\frac14\log(1-\rho^2),\qquad
2tv=\frac\rho{1-\rho^2},\quad g_t(0)=0.
\tag{3}
```

To verify (3), a self-coupling with covariance `rho v` has information at
least `-log(1-rho^2)/2`, by comparison to the bivariate Gaussian of the
same covariance. Gaussian couplings attain equality; minimize the resulting
one-variable function. Put `K_t=-g_t`. Direct differentiation gives

```math
g_t'(v)=-t(1-\rho),\qquad
vg_t'(v)=-\frac\rho{2(1+\rho)}.
\tag{4}
```

Thus `g_t` is decreasing, convex and `t`-Lipschitz, and is concave as a
function of `log v`. Equivalently `K_t` is increasing and concave, vanishes
at zero, and is `O_t(log(1+v))` at infinity.

For a symmetric real source `nu`, define the Bellman operator by

```math
(\mathcal Bf)(\nu)=\sup_\pi\left\{
\frac{f(\nu_+)+f(\nu_-)}2
-\frac12D(\pi\Vert\nu\otimes\nu)\right\}.
\tag{5}
```

Here the average absolute marginal of `(A,B)~pi` is the absolute law of
`nu`, and `nu_+,nu_-` are the symmetrized laws of `(A+B)/sqrt(2)` and
`(A-B)/sqrt(2)`. Averaging `pi` under simultaneous reversal and input
interchange preserves each child absolute law and decreases its entropy
cost. Hence the supremum is unchanged if `pi` has both signed marginals
equal to `nu`, is invariant under simultaneous reversal, and is invariant
under interchange. Its cost is then `I(A;B)`. Independent reversal of
only one input is not used.

## 2. Self-transport facts, including the unbounded equality case

We first establish the facts needed for the analytic stopping argument.

For a finite empirical law `nu_m`, the method of permutation contingency
tables gives the Gaussian permanent formula

```math
\frac1m\log\frac{\operatorname{per}
 [e^{-t|x_i-x_j|^2}]_{i,j=1}^m}{m!}\longrightarrow-F_t(\nu)
\tag{6}
```

when its fixed finite alphabet and frequencies tend to those of `nu`.
Indeed a table `gamma` has probability
`exp(-m D(gamma||nu_m tensor nu_m)+O(log m))`; there are polynomially many
tables, and feasible real tables admit integral approximation with bounded
count errors on a fixed alphabet.

For positive-semidefinite Gram matrices `C,D` of order `m`,

```math
\frac{\operatorname{per}(C\circ D)}{m!}
\ge\frac{\operatorname{per}C}{m!}
     \frac{\operatorname{per}D}{m!}.
\tag{7}
```

To prove this, write the two Gram matrices as Gram matrices of vectors,
take their tensor products over the coordinate index, and compare the
projection averaging simultaneous permutations with the smaller projection
averaging both factors separately. The difference of the projections is
positive semidefinite. Taking (6) for a finite joint source gives

```math
F_t(\operatorname{law}(X,Y))\le F_t(\operatorname{law}X)
                                      +F_t(\operatorname{law}Y).
\tag{8}
```

Equality holds for independent coordinates, by the product coupling and
the information inequality in the opposite direction. These statements
extend from finite laws in `W_2`: push an old self-coupling through source
transport kernels at its two endpoints. Data processing decreases entropy;
the quadratic cost changes by at most
`4t W_2(mu,nu)(sqrt(m_2(mu))+sqrt(m_2(nu)))`. Reverse the comparison.

We also need strict source concavity of `F_t`, not merely its finite-state
version. Here is a direct proof with no compact-space transport theorem.
Write `k(x,y)=exp(-t|x-y|^2)`. For a finite probability measure `alpha_n`,
the minimum in (2) exists and its matrix entries are positive: mixing with
the product measure has entropy derivative `-infinity` at a zero entry.
Lagrange multipliers and uniqueness then give a positive measure `mu_n`
such that

```math
\alpha_n=(k\mu_n)\mu_n,\qquad \|\mu_n\|_k^2=1.
\tag{9}
```

Choose finite approximations `alpha_n` converging in `W_2` to an arbitrary
finite-second-moment `alpha`. Every unit cube `Q` satisfies
`mu_n(Q)<=exp(td/2)`, since `k>=exp(-td)` on `Q times Q`. A subsequence
therefore converges vaguely to a positive Radon measure `mu`. The uniform
cube bound makes the Gaussian tail summable, uniformly on compact sets;
hence `k mu_n -> k mu` locally uniformly. Passing (9) against compactly
supported continuous functions and then using Tonelli gives

```math
\alpha=(k\mu)\mu,\qquad \|\mu\|_k^2=1.
\tag{10}
```

Kernel Cauchy--Schwarz gives `0<k mu<=1`. A fixed ball with positive
`mu` mass gives, for `f=-log(k mu)`,

```math
0\le f(x)\le t(|x|+R)^2-\log\mu(B_R).
\tag{11}
```

Thus `mu=e^f alpha`, `f` is `alpha`-integrable, and
`gamma(dx,dy)=k(x,y)mu(dx)mu(dy)` is a probability self-coupling. For every
competitor of finite entropy, the exact identity

```math
D(\theta\Vert\alpha\otimes\alpha)+t\mathbb E_\theta|X-Y|^2
=2\int f\,d\alpha+D(\theta\Vert\gamma)
\tag{12}
```

proves optimality and `F_t(alpha)=2 integral f d alpha`.

Gaussian energy is strictly positive on a nonzero difference `sigma` of
such scaling measures. One verification writes its energy as a positive
constant times the squared `L^2` norm of convolution with
`exp(-2t|.|^2)`. The cube bound makes the convolution well-defined and the
identity follows by polarization from the positive-measure identity.
If this convolution vanishes, its Fourier transform is a nonvanishing
Gaussian times the tempered distribution `sigma`; local division by this
smooth nonzero multiplier gives `sigma=0`. Equivalently use the entire
Laplace transform of the finite signed measure `exp(-t|x|^2)sigma`.

Let `alpha_0!=alpha_1` have scaling measures `mu_0,mu_1`. These measures
must differ by (10). For their mixtures `alpha_q,mu_q`, the log-sum
inequality and strict Gaussian-energy convexity give

```math
I(\alpha_q\Vert\mu_q)+\tfrac12\|\mu_q\|_k^2
<(1-q)[I(\alpha_0\Vert\mu_0)+\tfrac12]
 +q[I(\alpha_1\Vert\mu_1)+\tfrac12].
\tag{13}
```

Here `I(alpha||mu)=integral log(dalpha/dmu)dalpha` without a mass
correction. It is finite even if `mu` has infinite total mass, by (11);
the mixture inherits a quadratic bound on its logarithmic density. For
any such auxiliary measure, entropy comparison with the finite measure
`k mu tensor mu`, followed by `log z<=z-1`, gives

```math
\Phi_t(\alpha)\le I(\alpha\Vert\mu)
                         +\tfrac12\|\mu\|_k^2-\tfrac12.
\tag{14}
```

Equality holds for the scaling measures in (10), by (12). Consequently
`Phi_t` is strictly convex and `F_t` strictly concave on all
finite-second-moment source laws.

In (8), the product law is a maximizer among joint laws with fixed
marginals. Strict concavity therefore proves that equality in (8) holds
only for the product law.

For a safe Bellman pair with finite `I=I(A;B)`, every self-coupling of
`(A,B)` obeys

```math
I(AB;A'B')-I(A;A')-I(B;B')
=I(AA';BB')-2I\ge-I.
```

After adding costs, `F_t(pi)>=2F_t(nu)-I`. Combine this with (8) for the
orthogonal outputs to obtain the useful quantitative drift inequality

```math
\Phi_t(\nu)-\left\{\frac{\Phi_t(\nu_+)+\Phi_t(\nu_-)}2
-\frac I2\right\}\ge\frac I4.
\tag{15}
```

In particular `B Phi_t<=Phi_t`. Appending arbitrarily many iid pairing
levels, using the finite-variance `W_2` central limit theorem, proves

```math
g_t(m_2(\nu))\le(\mathcal B^r\Phi_t)(\nu)\le\Phi_t(\nu)\le0
\tag{16}
```

for symmetric sources. Translation invariance gives
`g_t(Var nu)<=Phi_t(nu)` without symmetry.

Finally, `Phi_t` is weakly continuous on each second-moment ball. For a
binary source mixture, concavity and disclosure of its mixture label give

```math
(1-q)F_t(\mu)+qF_t(\eta)
\le F_t((1-q)\mu+q\eta)
\le h(q)+(1-q)F_t(\mu)+qF_t(\eta).
\tag{17}
```

For a source of second moment at most `C`, remove its tail `|X|>R`, of
mass `q<=C/R^2`. Equations (16)--(17) bound the absolute change in `F_t`
by `h(q)+2q[K_t(C/q)+K_t(C/(1-q))]`, which tends uniformly to zero.
On a compact interval, weak convergence is `W_2` convergence. Truncation
at continuity radii completes the proof.

It follows that `B Phi_t` is weakly upper semicontinuous on moment balls:
almost-maximizing pair laws have weakly convergent subsequences, their
children have moment at most `2C`, and mutual information is lower
semicontinuous. Thus a maximizing pair exists and
`Gamma=Phi_t-B Phi_t` is lower semicontinuous. If `Gamma(nu)=0`, (15)
forces the maximizing inputs iid. Their outputs have a common symmetric
law, and zero drift forces equality in (8), hence independent outputs.
Writing `phi` for the real characteristic function of `nu`, independence
gives `phi(a+b)phi(a-b)=phi(a)^2phi(b)^2`; setting `a=b` and iterating
`phi(z)=phi(z/2^j)^(4^j)` identifies the centered Gaussian by its
finite-variance expansion at zero. Conversely the iid Gaussian policy
has zero drift. Therefore

```math
\Gamma(\nu)=0\quad\Longleftrightarrow\quad
\nu\text{ is a centered Gaussian, including }\delta_0.
\tag{18}
```

## 3. The conditional-variance envelope and a direct stopping closure

Define, now on every finite-second-moment real source,

```math
T_t(\nu)=\sup_L\{\mathbb E_Lg_t(\operatorname{Var}(X\mid L))
-I(X;L)\},\qquad J_t=-T_t.
\tag{19}
```

The supremum permits arbitrary standard-Borel labels; infinite information
has value `-infinity`. On a finite source it is unchanged by restricting
to finite labels: the posterior-simplex reward is continuous, and its
concave hull has a finite supporting mixture. The trivial label and
conditional self-coupling construction give

```math
g_t(\operatorname{Var}\nu)\le T_t(\nu)\le\Phi_t(\nu).
\tag{20}
```

For the upper comparison, a mixture of conditional self-couplings has
information at most `I(X;L)+E F_t(nu_L)` after including cost. Thus
`Phi_t(nu)>=E Phi_t(nu_L)-I(X;L)/2`, which implies (20).

The following continuity fact removes the need for an auxiliary smaller
envelope or a Gaussian-boundary fixed-point construction. Since `K_t` is
increasing and concave, conditional variance decomposition implies

```math
K_t(\operatorname{Var}(X\mid L))
\ge\mathbb E[K_t(\operatorname{Var}(X\mid L,R))\mid L].
\tag{21}
```

For a mixture label `R`, append it conditionally independently of `L`
given `X`, so `I(X;L)>=I(X;L|R)`. Taking infima in
`J_t=inf_L[I+E K_t(Var)]`, and then disclosing `R` for the opposite
bound, proves the binary mixture sandwich (17) with `F_t` replaced by
`J_t`. In particular `J_t` is source-concave.

To compare two sources in `W_2`, retain an old label through an optimal
source coupling. Information decreases by data processing, while

```math
\mathbb E|\operatorname{Var}(X\mid L)-\operatorname{Var}(Y\mid L)|
\le W_2(\mu,\nu)(\sqrt{m_2(\mu)}+\sqrt{m_2(\nu)}).
\tag{22}
```

This is Cauchy--Schwarz applied to the two centered conditional random
variables. The `t`-Lipschitz property of `K_t` proves `W_2` continuity of
`J_t`. Since `0<=J_t(nu)<=K_t(m_2(nu))`, the same tail argument as in
(17), now with only one copy of `K_t`, proves weak continuity of `T_t`
on each second-moment ball.

Next `B T_t<=T_t`. Choose a safe pair `(A,B)` and finite-information child
channels `M|U,N|V`, independently conditional on `(U,V)`, and set
`L=(M,N)`. Use labels `(L,B)` for `A` and `L` for `B`. The information
identity is

```math
I(A;L,B)+I(B;L)=I(A,B;L)+I(A;B)
\le I(U;M)+I(V;N)+I(A;B).
\tag{23}
```

Conditional on `L`, let `Sigma` be the covariance of `(A,B)`, and put
`b=Sigma_22`, `c=det(Sigma)/b`. Conditional linear regression gives
`E_{B|L}Var(A|B,L)<=c`. Each of `b,c` lies between the two eigenvalues
of `Sigma`, and their product is the eigenvalue product. Concavity of
`g_t` in log variance therefore gives
`g_t(c)+g_t(b)>=sum_i g_t(lambda_i)`. Convexity and diagonal majorization
give `sum_i g_t(lambda_i)>=g_t(Var(U|L))+g_t(Var(V|L))`. Jensen and
monotonicity thus give

```math
\mathbb E_{B\mid L}g_t(\operatorname{Var}(A\mid B,L))
+g_t(\operatorname{Var}(B\mid L))
\ge g_t(\operatorname{Var}(U\mid L))
   +g_t(\operatorname{Var}(V\mid L)).
\tag{24}
```

If `b=0`, use `c=Sigma_11`; if `b>0` and the determinant is zero, use
`c=0`. These also follow by continuity. Refining labels increases the
expected `g_t` reward, by variance decomposition and convexity. Average
(24), coarsen the two child rewards to `M` and `N`, subtract (23), and
divide by two. Each parent marginal is `nu`; taking the child and pair
suprema proves the supersolution assertion. All integrals are finite by
the second-moment bounds. In particular the argument is valid both for
finite reachable sources and for the unbounded extension (19).

Here is the direct quantitative closure. Put

```math
f_r=\mathcal B^r\Phi_t,\quad
D=\Phi_t-T_t,\quad E_0=m_2(\nu_0),\quad
B_0=\Phi_t(\nu_0)-g_t(E_0).
```

By continuity, compactness of moment balls, and (18)--(20),

```math
\kappa(\epsilon,C)=\inf\{\Gamma(\nu):m_2(\nu)\le C,
D(\nu)\ge\epsilon\}>0
\tag{25}
```

when the set is nonempty. Along a uniform branch of a depth-`r` policy,
the second moments form a nonnegative martingale. For a policy within
`zeta` of `f_r(nu_0)`, (16) and telescoping its actual drift show
`E sum_d Gamma(nu_d)<=B_0+zeta`. Stop upon `D<=epsilon` or second moment
exceeding `C`, otherwise at depth `r`. The surviving probability is at
most `(B_0+zeta)/(r kappa)`. Any remaining continuation is bounded above
by `Phi_t`. At the three stop types its excess above `T_t` is bounded
respectively by `epsilon`, `K_t(second moment)`, and `K_t(C)`.
Bounded optional stopping preserves expected second moment `E_0`.
Iterating `B T_t<=T_t` only to the bounded stopping time and taking
`zeta` to zero therefore proves

```math
f_r(\nu_0)-T_t(\nu_0)
\le\epsilon+E_0\sup_{s>C}\frac{K_t(s)}s
 +\frac{K_t(C)B_0}{r\kappa(\epsilon,C)}.
\tag{26}
```

Let first `C` grow, then `epsilon` decrease, then `r` grow. Thus
`lim_r f_r(nu_0)<=T_t(nu_0)`. No claim that `f_r>=T_t` was used; the
drift budget is `Phi_t-G_t`, not `Phi_t-T_t`. This is precisely the
inequality needed below.

## 4. Exact root certificate

For `nu_p=(1-p)delta_0+(p/2)(delta_{-1/sqrt(p)}+delta_{1/sqrt(p)})`, a
posterior is described, after reflecting and retaining its reflection
label, by `z,s in [0,1]`. Its variance and entropy are

```math
v(z,s)=\frac{z-z^2s^2}{p},\qquad
H(z,s)=h(z)+z h((1+s)/2).
```

Reflecting each posterior with equal probability automatically enforces
the two signed barycenter constraints; the sole remaining constraint is
`E z=p`. Hence the exact offset is

```math
p\log2+T_t(\nu_p)+t(1-\sqrt p)
=-h(p)+t(1-\sqrt p)
 +\sup_{\mathbb Ez=p}\mathbb E[H(z,s)+g_t(v(z,s))].
\tag{27}
```

The complete integer/rational verifier is
`computations/continued_feedback_conditional_variance_exact_certificate_2026_09_06.py`.
It was rerun for this reconstruction, without changing the canonical
result file, using

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python computations/continued_feedback_conditional_variance_exact_certificate_2026_09_06.py --output /home/math/quadra/tmp/transfer_reconstruction_exact_certificate_2026_09_06.json
```

For completeness, the finite computation and all interpolation payments
are as follows. Let `N=2500`, `S=10^12`. At `z=i/N,s=j/N`, use outward
upper entropies and the upper Gaussian value at variance

```math
\frac1{1000}\left\lfloor
\frac{32000(iN^3-i^2j^2)}{31N^4}\right\rfloor.
```

This rounds variance downward, so is an upper bound since `g_t` decreases.
Take the largest value for each of the `N+1` choices of `i`, and its exact
least concave majorant at mean `p`. The decreasing-slope stack computes
that majorant using rational cross-products. Its offset upper value is

```math
-\frac{13307450618189}{1056000000000000}.
```

Round arbitrary `z` randomly to adjacent grid points preserving its mean;
round `s` to the nearest grid point. Since
`|partial_z v|<=1/p`, `|partial_s v|<=2/p`, binary entropy has modulus
`h(delta)`, and `g_t` is `t`-Lipschitz, the entire continuous-to-grid
payment is at most

```math
\omega=h(1/N)+(\log2)/N+h(1/(4N))+2t/(pN)
\le\frac{630156538579809}{77500000000000000}.
```

The sum is exactly `-a`. For the transcendental enclosures, square roots
are bracketed by integer square root at scale `S`; logarithms use binary
range reduction followed by
`log y=2 sum_{j=0}^{24}w^(2j+1)/(2j+1)+R`, `w=(y-1)/(y+1)`, with
`0<=R<=2w^51/[51(1-w^2)]` for `1<=y<=2`, then outward rounding.
In (3) the upper `rho` endpoint is used in its increasing linear term,
the lower endpoint in its decreasing logarithmic term. Every grid
operation uses signed 64-bit integers; the asserted bounds
`32000 N^4<2^63`, `N S<2^63` cover the largest products. There are no
floating transcendental inputs to the certificate.

Consequently (26) implies that for every `0<a'<a` there is a finite
integer `r` such that

```math
p\log2+(\mathcal B^r\Phi_4)(\nu_p)+4(1-\sqrt p)\le-a'.
\tag{28}
```

The depth is fixed before any matrix order tends to infinity. It need
not be effective for the claimed existence theorem.

## 5. A uniform orbital estimate and the exact finite-depth type recursion

Let `G_m` be the signed-permutation group and write

```math
P_t(v)=\mathbb E_{g\in G_m}e^{-t\|v-gv\|^2},\qquad L_t(v)=\sqrt{P_t(v)}.
```

The signed kernel is
`K_t(a,b)=(e^{-t(a-b)^2}+e^{-t(a+b)^2})/2`, so
`P_t(v)=per K_t[|v|]/m!`. For a fixed finite symmetric type its logarithmic
limit is `-F_t`; optimize the relative sign conditionally on the two
magnitudes in (6). The common global sign is fair and causes no extra
entropy cost.

For every orthogonal `U` and `||v||^2<=Cm`,

```math
\mathbb E_gL_t(Ugv)\le\exp(O_{t,C}(\sqrt m))L_t(v),
\tag{29}
```

uniformly in `U,v`. To prove this, use the symmetric-Fock feature

```math
\phi(v)=e^{-t\|v\|^2}\bigoplus_{d\ge0}
\frac{(2t)^{d/2}}{\sqrt{d!}}v^{\otimes d}.
```

If `Q_G` projects to the invariant subspace, `P_t(v)=||Q_G phi(v)||^2`.
The covariance `E_g |phi(gv)><phi(gv)|` has operator norm `P_t(v)`:
its finite orbit Gram matrix has nonnegative entries and constant row
sum `P_t(v)`. Degrees at most `D` have invariant rank at most
`sum_{j<=D/2}p(j)<=exp(pi sqrt(D/3))`; invariant polynomials are symmetric
polynomials in the squared coordinates. The last bound follows from
the partition generating function and
`log product_{j>=1}(1-e^{-sj})^{-1}<=pi^2/(6s)`.
The omitted squared feature norm is the Poisson tail with mean
`2t||v||^2`. Choose `D=ceil(e^2(2tC+1)m)`; its Chernoff bound is at most
`exp(-(4tC+1)m)<=P_t(v)`. Projection, the covariance norm, and
Cauchy--Schwarz now give (29).

Invariant-subspace inclusion for a split list gives
`L_t(v_+,v_-)<=L_t(v_+)L_t(v_-)`. Permanent terms fixing a removed
coordinate and `K_t(a,a)>=1/2` give

```math
\max_i L_t(v\setminus v_i)\le\sqrt{2m}\,L_t(v).
\tag{30}
```

Build a random normalized Hadamard basis recursively, at each internal
node of size `s=2q`, by

```math
U_s=\operatorname{diag}(U_q^{(1)},U_q^{(2)})
\frac1{\sqrt2}\begin{pmatrix}I&I\\I&-I\end{pmatrix}g_s.
\tag{31}
```

Every node has a fresh independent uniform signed permutation. Every
terminal node is an arbitrary deterministic normalized Hadamard followed
on its input by its own independent signed permutation. Child matrices
are mutually independent and independent of all parent permutations.
Equation (31) is orthogonal and all its entries have magnitude `s^-1/2`.

For a signed input with symmetrized finite type `nu`, the number of words
with its absolute type is `s! 2^(number nonzero)/product c_a!`. A pair
table has `(s/2)!/product c_ab!` ordered-pair arrangements. Stirling's
formula therefore gives probability
`exp(-s D(pi||nu tensor nu)/2+O(log s))` even when the pair's individual
signed marginals are unequal. The average absolute marginal is exactly
the constraint in (5).

Fix the depth `r`. Apply (30), the split inequality, the pair-table
probability, and (29) at the terminal nodes. Summing over all intermediate
types gives

```math
\limsup_{m\to\infty}\frac1m\log\mathbb E L_t(U_m v_m)
\le(\mathcal B^r\Phi_t)(\nu).
\tag{32}
```

Here the input alphabet is fixed finite (or uniformly convergent to one),
its type tends to `nu`, and its energy is `O(m)`. All intermediate
alphabets are finite at this fixed depth. Their number of types and
factorial errors are polynomial in `m`; there are `2^r` orbital losses,
each uniform in its terminal Hadamard. The total logarithmic error is
`O_{r,t,C}(sqrt(m))+O_r(log m)=o(m)`. Compactness of the finite table
polytopes, entropy continuity at zero, and continuous kernel entries
justify taking limiting root frequencies and amplitudes. No estimate
uniform in growing depth is asserted.

## 6. Exact all-spin weave counting

Choose an independent order-`m` Hadamard `H_i` in each fibre, and a
symmetric full sign matrix `S`. Define the symmetric full sign matrix

```math
W_{(i,a),(j,b)}=S_{ij}H_i(a,j)H_j(b,i).
```

Retain `k` rows in each fibre, with selectors `T_i`. Its principal matrix
`W_T` has order `N=mk`. For block spins `x_i`, put
`h_i=H_i[T_i,:]^T x_i`, `E=x^TW_Tx`. Row orthogonality gives
`sum_i ||h_i||^2=m^2k`, and symmetry of `S` gives the exact identity

```math
D_\sigma(x):=\sum_{i,j}(h_i(j)-\sigma S_{ij}h_j(i))^2
=2(m^2k-\sigma E),\qquad\sigma\in\{-1,1\}.
\tag{33}
```

Independently permute output columns in every fibre and choose all
off-diagonal signs of `S` independently fair. Conditional on the absolute
diagonal coordinate selected from each row spectrum, discard the
nonnegative diagonal defect terms. For each undirected edge, averaging
its sign in `exp(-tD_sigma/(2k))` gives exactly
`K_t(|h_i(j)|/sqrt(k),|h_j(i)|/sqrt(k))`: the edge occurs twice in (33).

The remaining entries of each row are a uniform permutation of its
multiset. The edge kernel is positive semidefinite (expand its `cosh`
factor). Split a Gram factor into each endpoint tensor and apply graph
Cauchy--Schwarz: a tensor contraction on a graph is bounded in absolute
value by the product of the Euclidean norms of its vertex tensors. Each
edge variable belongs to precisely two tensors, so this is the product
Hölder/Finner inequality with exponent two. For a row multiset `b` of
length `m-1`, that squared norm, after normalizing the permutation tensor,
is `per K_t[b]/(m-1)!`. Repeated coordinate values give the same formula,
because each distinct arrangement has the same number of permutations.

Writing

```math
Z_i=\sum_{x_i\in\{-1,1\}^k}
\max_j L_t\left((H_i[T_i,:]^Tx_i/\sqrt k)\setminus j\right),
```

Markov's inequality and the union over every spin and both signs give

```math
\mathbb P\{\exists x,\sigma:D_\sigma(x)\le2\gamma m^2k\}
\le2e^{t\gamma m^2}\prod_{i=1}^m Z_i.
\tag{34}
```

Use `H_i=sqrt(m)(U_m^{(i)})^T` from independent copies of (31), in
addition to the fresh output permutations and signs already used.
For fixed selectors of size `k`, input signed-permutation invariance,
(30), and (32) yield, when `k/m -> p`,

```math
\limsup\frac1m\log\mathbb EZ_i
\le p\log2+(\mathcal B^r\Phi_t)(\nu_p).
\tag{35}
```

In fact the normalized vector entering (32) is
`v=(1_T x)/sqrt(k/m)` and has second moment exactly one. The independent
fibre choices in (34) produce `(E Z_i)^m`, not `E[Z_i^m]`.

Fix `0<eta<a'/4` and take
`gamma=1-sqrt(k/m)+eta`. Equations (28), (34), and (35) give a failure
probability `exp((-a'+4 eta+o(1))m^2)`, hence below one for all sufficiently
large admissible `m`. An actual full sign matrix therefore exists with

```math
\max_x|x^TW_Tx|\le m^2k(\sqrt{k/m}-\eta).
```

Deleting its diagonal costs at most `N/2` in `Q`, so its hollow signing
satisfies

```math
\frac{Q(A)}{N^{3/2}}
\le\frac12-\frac\eta{2\sqrt{k/m}}+\frac1{2\sqrt N}.
\tag{36}
```

This checks both factors of two: the undirected-edge exponential weight
in (34) and the conversion `Q=|x^TAx|/2` in (36).

## 7. Elementary all-order realization and the ordered limits

Only two fixed terminal Hadamards are needed: the order-two matrix and
an order-twelve matrix. For the latter, on `F_11` let `chi` be the
quadratic character with `chi(0)=0` and `C_xy=chi(x-y)`. Direct character
summation gives `C^T=-C`, `C1=0`, `CC^T=11I-J`; alternatively these three
finite identities are checked with integer arithmetic. Thus

```math
R=\begin{pmatrix}0&1^T\\-1&C\end{pmatrix},\qquad
H_{12}=I+R
```

has sign entries and `H_12 H_12^T=12I`. Their Kronecker products provide
all terminal orders `s=2^a12^b`, `a,b>=0`.

These orders are multiplicatively asymptotically dense. Indeed
`alpha=log(12)/log(2)` is irrational, since an equality `12^u=2^v` would
force `3^u=1`. For every `delta>0`, finitely many nonnegative multiples
`b alpha` form a `delta`-net modulo one. One elementary proof uses
pigeonhole to find `0<||q alpha||<delta`, then its consecutive multiples
to cover the circle. For every sufficiently large target logarithm `y`,
choose one of those finitely many `b` and a nonnegative integer `a` with
`y<=a log2+b log12<=y+delta log2`. Hence the least available terminal
order above any real `x` has ratio to `x` tending to one. No prime-number
theorem or short-interval result is needed.

For the fixed depth `r`, use `m=2^r s`, `k=floor(pm)`. The available
orders `N=m floor(pm)` are likewise multiplicatively asymptotically
dense, and (36) holds at all sufficiently large such orders. If `A_I`
is any principal `n`-vertex restriction of a hollow signing `A`, then
`Q(A_I)<=Q(A)`: extend its fixed spin by independent mean-zero outside
spins and average the full Hamiltonian. Taking the least constructed
`N>=n` therefore extends (36) to every `n`, since `N/n ->1`.

For each fixed `a'<a` and `eta<a'/4`, first choose the finite depth from
(28), and then take the all-order limit. Finally let `eta` increase to
`a'/4` and `a'` increase to `a`. This gives exactly (1). The upper
square-root endpoint `246062746063/250000000000` for `sqrt(31/32)` gives the
rational outward cap

```math
\frac{80459630021641337701}{161102201102367360000}
<\frac{499432220485404}{10^{15}}.
```

The construction controls every spin and every sufficiently large order,
but does not compare the limsup to the original liminf. Its seed remains
a Hadamard basis, not an arbitrary minimizing signing. The strict upper
theorem survives reconstruction; the direct-`T` stopping argument and
two-generator order supply reduce its analytic and arithmetic dependencies.

## 8. Independent exact regression checks

Run

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python computations/transfer_reconstruction_exact_checks_2026_09_06.py
```

The checker verifies the order-twelve Hadamard in integers, all 240 signed
ternary words and their 50 pair tables, the safe symmetrizations and both
child absolute laws, 24 permanent-projection cases, a nontrivial four-vertex
PSD-kernel contraction with repeated coordinates, the weave/defect identity
on all 4096 spins of a restricted example, and 138 independent rational
logarithm/square-root interval comparisons. Its decimal-cap comparison is
also exact. Adding `--full-certificate` reruns the complete root certificate.
These finite checks guard normalizations and arithmetic; the uniform and
asymptotic assertions are proved above, not inferred from these tests.
