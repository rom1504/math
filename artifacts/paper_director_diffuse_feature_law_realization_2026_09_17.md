# Physical realization of non-Gaussian laws on diffuse finite-rank features

2026-09-17. **Proved; independently reconstructed by all three paper tracks.** This combines finite-feature
replacement, bounded-likelihood change of measure, and an exact sign-
covariance repair. It is not an assumption that arbitrary near-minimizer
ground sets have finite-rank structure. External novelty is unestablished;
the finite-dimensional tilt and Lindeberg ingredients are classical.

## 1. Precise realization theorem

Fix r and an even nonnegative function w:R^r->R whose derivatives
through order five, including w itself, are bounded and continuous.
Smooth compact support is sufficient but not necessary. For G standard
Gaussian in R^r assume

    E w(G)=1,       E w(G) GG^T=I_r.                        (1)

Let U=U_n be an n-by-r matrix with U^T U=I_r, and let u_i be its row
vectors. Suppose d_n=max_i ||u_i||_2 ->0. There are centered physical
sign laws nu_n on {+-1}^n, with EXACT covariance I_n and a common finite
subGaussian proxy K_w I_n, such that uniformly over ALL Boolean x,

    E_(nu_n)|h.x|/sqrt(n)
      = E |a_x.Y+sqrt(1-||a_x||^2)Z|+O_(r,w)(sqrt(d_n+n^(-1/2))),
    a_x=U^T x/sqrt(n),                                     (2)

where Y has density w relative to standard Gaussian measure, and Z is
an independent standard scalar Gaussian. The law can have positive mass
on every physical sign vector. Neither U's entries nor its columns need
be signs; no Hadamard, block, lattice or permutation symmetry is assumed.

The rank and w are fixed BEFORE the order limit. Constants can be very
large. This statement does not assert a uniform growing-rank theorem.
The whole Boolean query set in (2) does not introduce a union entropy:
the replacement bound is deterministic and uniform in its coefficients.

## 2. The initial cube tilt and its covariance

Let epsilon be uniform independent signs, X=U^T epsilon, and

    Z_n=E w(X),    d(tilde_nu)/d(uniform)=w(X)/Z_n.

Lindeberg replacement, matching moments through order three and using
bounded fourth derivatives of w, gives

    |Z_n-1|<=C_(r,w) sum_i ||u_i||^4<=C_(r,w)r d_n^2.       (3)

Thus this is a genuine nonnegative probability law for large n.
It is centered by global sign symmetry. Its covariance diagonal is
exactly one; write S_n=Cov_(tilde_nu)(h)=I+E_n.

We prove ||E_n||op<=C_(r,w)d_n. For i!=j, condition on all other signs
and put X_ij=sum_(k!=i,j)u_k epsilon_k. Exact double differencing gives

    E epsilon_i epsilon_j w(X)
      =(1/4) integral_[-1,1]^2
        u_i^T E Hess w(X_ij+s u_i+t u_j) u_j ds dt.          (4)

Smooth Lindeberg replacement for Hess w costs at most
C_w sum_k||u_k||^3<=C_w r d_n; this uses bounded fifth derivatives.
Removing the deterministic shifts and restoring the two omitted Gaussian
coordinates costs at most C_w(||u_i||+||u_j||). The full Gaussian
sum sum_k u_k g_k is standard in R^r. Gaussian integration by parts
and (1) give

    E Hess w(G)=E[(GG^T-I)w(G)]=0.

Consequently the absolute value of (4) is at most
C_(r,w)d_n ||u_i||||u_j|. After division by Z_n, for any real v,

    |v^T E_n v|
      <=C_(r,w)d_n [sum_i |v_i| ||u_i||]^2
      <=C_(r,w)r d_n ||v||^2.                              (5)

The zero covariance diagonal is retained exactly in this estimate.
It is not enough to show entrywise covariance convergence: (5) supplies
the operator control needed for the physical repair.
In fact the same entrywise majorant proves the STRONGER estimate
||E_n||F<=C_(r,w)d_n, since sum_i||u_i||^2=r. This permits the
entropy-optimal repair in Section7.

## 3. Exact covariance repair by an actual sign law

Set e_n=||E_n||F. If e_n>0, choose

    theta_n=4e_n/(1+4e_n),
    D_n=-(1-theta_n)E_n/theta_n,
    R_n=I+sin((pi/2)D_n),                                  (6)

where sine is entrywise and D_n has zero diagonal. Then
||D_n||F=1/4 and ||D_n||op<=1/4.
For any matrix D, D^(entrywise k) is a compression of D^(tensor k).
Taking absolute values in the sine series therefore proves

    ||R_n-I||op<=sinh(pi/8)<1/2.

So R_n is a correlation matrix with spectra in [1/2,3/2]. For
h_repair=sign N(0,R_n), the exact arcsine identity yields

    Cov(h_repair)=I+D_n.                                   (7)

Indeed every entry of (pi/2)D_n lies in [-pi/8,pi/8], on which
arcsin(sin t)=t. Mixing tilde_nu with this repair law at weight theta_n
gives covariance exactly I. If e_n=0 no repair is needed.

This is an exact sign realization, not independent rounding of a weighted
matrix. The correction mass is O_(r,w)(d_n). Its normalized scalar
response is at most sqrt(5/4) by (7), so its cost in (2) is O(d_n),
not O(sqrt(n)d_n). The uncorrected law also has bounded covariance.
An additional mixture of weight n^(-1) with independent signs supplies
full physical support if desired, preserving exact isotropy and the
displayed error rate.

The tilted law has density at most M=max(1,2||w||infinity) for large n.
Symmetry and the Rademacher even-moment bound imply for every v

    E_(tilde_nu) exp(v.h)
      <=1+M[exp(||v||^2/2)-1]
      <=exp(M||v||^2/2).

Gaussian product Holder bounds the repair law by the same inequality
with M replaced by3/2. Mixing preserves the common bound. Thus one may
take K_w=max(M,3/2); no covariance-only subGaussian inference is used.

## 4. Uniform physical-query replacement, including singular joint limits

For Boolean x, consider the (r+1)-dimensional vector

    (U^T epsilon, x.epsilon/sqrt(n))
      =sum_i (u_i,x_i/sqrt(n)) epsilon_i.

Its covariance is exactly the covariance of
(G,a_x.G+sqrt(1-||a_x||^2)Z). The latter may be singular when
||a_x||=1; the proof does not invert this covariance.

Let h=d_n+n^(-1/2), and smooth the absolute value by
phi_tau(v)=sqrt(v^2+tau^2), with tau=sqrt(h). Its uniform approximation
error is at most tau. For F(z,v)=w(z)phi_tau(v), fourth derivatives are
bounded by C_w tau^(-3)(1+|v|), when tau<=1. In Lindeberg replacement,
the scalar v always has uniformly bounded moments, including intermediate
mixed Gaussian/sign inputs. Taylor remainders and matching the first
three moments therefore give a uniform bound

    C_(r,w) tau^(-3) sum_i ||(u_i,x_i/sqrt(n))||^4
      <= C_(r,w) tau^(-3) h^2 = C_(r,w)sqrt(h).             (8)

Terms with the unbounded replacement coordinate are paid by its fixed
Gaussian moments; the extra fifth-power coefficient sum is smaller.
The total squared coefficient sum is r+1. Undo smoothing and divide
by Z_n from (3). This proves (2) for the tilt; Section3 repairs it with
a smaller error. No independent replacement of the physical query
directions, no selection of one favorable spin, and no growth of r with n
has been smuggled into (8).

## 5. Genuine non-Gaussian variance allocation

The target law in (1) can have uniformly small absolute projections
while retaining covariance I. Choose a smooth even compactly supported
random vector V with Cov(V)=I_r. Fix p in (0,1) and small a>0. Let

    Y=aV with probability1-p, and Y=bV with probability p,
    b=sqrt([1-(1-p)a^2]/p).

This has smooth compact density, zero mean and covariance I. Its density
relative to Gaussian measure is an admissible w. For every real t,

    E|t.Y|<=(a+sqrt(p))||t||.                              (9)

Hence the realized physical law satisfies simultaneously

    E|h.x|/sqrt(n)
      <=(a+sqrt(p))||a_x||+kappa sqrt(1-||a_x||^2)+o(1).     (10)

In particular all Boolean queries within relative squared distance eta
of the diffuse feature space have response at most
a+sqrt(p)+kappa sqrt(eta)+o(1). This can beat the sharp bounded-band
Gaussian-mixture floor, without abandoning physical signs or exact
isotropy. The price is a potentially large fixed subGaussian constant
and restriction to finitely many diffuse features.

## 6. What this does and does not transfer

This realizes an entire non-Gaussian finite-dimensional response profile,
not just a desired covariance. It separates two previously conflated
requirements: exact physical sign realization is possible here; discovering
a feature space capturing the actual high-energy queries is a separate
geometric theorem. The construction generalizes coordinate-block slice
laws to arbitrary diffuse real feature spaces.

No such fixed-rank or uniformly manageable growing-rank feature space is
proved for original near-minimizers. The actual low-cap theorem earlier
in this campaign controls high-energy query laws without that hypothesis,
but its covariance projector need not have fixed rank or delocalized rows.
These facts cannot currently be composed into an original recurrence.

This theorem also does not approximate a Gaussian target with covariance
I by a Gaussian mixture and claim a free gain: the target Y is deliberately
non-Gaussian. Fixed-rank Gibbs conditioning/replacement are classical
mechanisms; the exact isotropic physical repair and all-query statement
are the explicitly proved combination under consideration.

## 7. Exact isotropy has asymptotically no additional information cost

The same construction satisfies the sharp entropy identity

    D(nu_n || uniform signs) -> D(law(Y) || N(0,I_r)).       (11)

Moreover ANY sequence of sign laws whose U_n projections converge weakly
to Y has liminf relative entropy at least the right side, even without
an isotropy requirement. Thus enforcing physical signs and EXACT
isotropy does not add an asymptotic cost in this fixed-rank regime.

For the initial tilt the exact divergence is

    D(tilde_nu||uniform)=E[w(X)log w(X)]/Z_n-log Z_n.

The bounded continuous function w log w, with value zero at zeros of w,
has Gaussian CLT convergence, proving the desired limit for the tilt.
The Frobenius choice in Section3 is important: it gives

    ||R_n-I||F<=pi/8,   ||R_n-I||op<1/2.

For Gaussian laws, the exact entropy formula and Tr(R_n-I)=0 give

    D(N(0,R_n)||N(0,I_n))
      =(1/2)[Tr R_n-n-log det R_n]
      <=||R_n-I||F^2<=pi^2/64.

Data processing under the coordinate sign map bounds the physical
repair entropy by the same constant. Convexity of relative entropy and
theta_n=O(d_n) therefore preserve the upper limit in (11). A vanishing
independent-cube mixture for full support has zero added entropy.

For the lower limit, use the finite Gibbs variational inequality with
bounded continuous test functions phi(U_n^T h). The uniform reference
projections converge to G, while the proposed laws' projections converge
to Y. Thus

    liminf D(nu_n||uniform)>=E phi(Y)-log E exp(phi(G)).

Taking the supremum over bounded continuous phi gives D(Y||G). In this
specific density setting it suffices to use
phi_M=log(max(w,exp(-M))) and send M to infinity AFTER n. These tests
are bounded and continuous, including outside the support of w.

## 8. Sharp leading information price for uniformly cheap projections

Let I_r(epsilon) be the infimum of D(Y||G_r) among centered laws with
Cov(Y)=I_r and

    sup_(||t||=1) E|t.Y|<=epsilon,   0<epsilon<=1.

There is a UNIVERSAL constant C such that

    r log(1/epsilon)+r[log(pi/2)/2-1/2]
      <= I_r(epsilon) <= r log(1/epsilon)+Cr.                (12)

The upper bound may be attained up to this additive Cr by smooth compact
even densities, so Section7 physically realizes these costs on every
fixed-rank diffuse interface. This is a constructive application of
classical entropy inequalities, not a claim that maximum-entropy duality
itself is new.

For the lower bound, each coordinate has mean absolute value at most
epsilon. The entropy-maximizing scalar law under that constraint is the
two-sided exponential density; equivalently its nonnegative relative
entropy gives h(Y_i)<=1+log(2epsilon). Subadditivity and Cov(Y)=I give

    D(Y||G_r)=r log(2pi e)/2-h(Y)
      >=r log(1/epsilon)+r[log(pi/2)/2-1/2].

Singular or infinite-divergence laws cannot defeat the lower bound.
This also follows directly by testing the product two-sided exponential
law in the KL variational formula, avoiding differential entropy when
individual entropies are not finite.

For the upper bound take the product of r copies of one fixed smooth
compact even variance-one scalar law V_1 with D(V_1||G_1)=D0<infinity.
In Section5 set a=epsilon/2 and p=epsilon^2/4. Then b>1 and (9) meets
the required response bound. By convexity and the exact scaling formula,

    D(Y||G_r)
      <=rD0-r[(1-p)log a+p log b]
      <=r log(1/epsilon)+r(D0+log2).

The averaged quadratic terms cancel because (1-p)a^2+pb^2=1.
The small active component supplies variance, not a negligible uncharged
tail. The state remains isotropic while its mean absolute projections
are small; its information cost grows logarithmically with their inverse.

## 9. Sharp rank limitation for bounded-covariance query laws

Independent scope audit by the discrepancy track. The response profile
in (2) has the pointwise lower bound

    F_Y(a)=E|a.Y+sqrt(1-||a||^2)Z|
      >=kappa sqrt(1-||a||^2)>=kappa(1-||a||^2).             (13)

Indeed, conditional on Y, the absolute first moment of a shifted centered
Gaussian is minimized at shift zero. This argument does not require Y
to be isotropic or even centered; those hypotheses serve other parts of
the realization theorem.

Let mu_n be ANY probability law on the Boolean query set, and suppose
Sigma_mu=E_mu xx^T has operator norm at most L. If the feature space has
rank r, then

    E_mu ||a_x||^2=Tr(U^T Sigma_mu U)/n<=Lr/n.

Averaging (2) and (13) therefore gives

    E_mu E_(nu_n)|h.x|/sqrt(n)
      >=kappa(1-Lr/n)-O_(r,w)(sqrt(d_n+n^(-1/2))).          (14)

Thus a fixed-rank realization cannot give a nonvanishing uniform response
discount on any code admitting a uniformly bounded-covariance query law.
More generally, the same conclusion holds for ranks r=o(n) whenever a
separately justified all-query realization error tends to zero uniformly;
the present theorem alone does not justify such growing ranks.

This is a concrete obstruction to inserting the construction directly
into the low-covariance branch of the signed-energy argument. There the
compressed signed covariance B obeys ||B||op<=3. If the branch supplies
||B||F^2>=b0 n for fixed b0>0, then

    rank(B)>=||B||F^2/||B||op^2>=b0 n/9.

Its relevant signed-energy component can therefore have macroscopic
rank. The exact finite-rank physical realization and optimal information
cost proved above do not establish a cheap law for this bulk component.
Neither exact isotropy repair nor a larger fixed subGaussian constant
removes the rank restriction in (14).

## 10. A codewise covariance-capture consequence

For a nonempty Boolean code C_n and fixed positive integer r, define

    s_r(C_n)=inf_(mu on C_n) [sum_(j=1)^r lambda_j(E_mu xx^T)]/n,

with missing eigenvalues counted as zero. For every fixed tolerance
theta>0 there is a finite K=K(r,theta) such that, for all sufficiently
large n, ONE centered exactly isotropic K-subGaussian physical sign law
of full support satisfies

    sup_(x in C_n) E|h.x|/sqrt(n)
       <=kappa sqrt(1-s_r(C_n))+theta.                    (15)

The threshold in n and K are uniform over the code. This is a finite-rank
feature-capture theorem, not an assertion that minimizing ground codes
have s_r bounded away from zero for fixed r.

Proof. Fix a query probability mu and its covariance Sigma. Among its
top r eigenvectors retain only eigenvalues at least eta n, for a fixed
eta>0. Their span captures mass at least s_r(C_n)-r eta. It is uniformly
diffuse: if u is a retained eigenvector, Sigma_ii=1 implies
u_i^2<=1/(eta n), and the retained frame's row norms are at most
sqrt(r/(eta n)). Apply Sections1--5 to a fixed target law with every
absolute projection at most epsilon. Formula (10) and Jensen give

    E_mu E_nu |h.x|/sqrt(n)
      <=epsilon+kappa sqrt(1-s_r(C_n)+r eta)+o(1).

All constants are uniform in mu: maximize over the finitely many retained
ranks 0,...,r, and use the displayed leverage bound. Choose epsilon,eta
so the excess is less than theta/2. The finite-cube set of centered,
exactly isotropic K-subGaussian laws is compact and convex. Finite minimax
therefore exchanges its infimum and the supremum over query probabilities,
giving one law for ALL queries. Add uniform mass n^(-1) if necessary;
the strict tolerance absorbs its vanishing cost. This proof was
independently reconstructed by the discrepancy track.

## 11. Reproducible finite diagnostic, not an asymptotic certificate

`computations/paper_director_2026_09_17_diffuse_feature_realization.py`
uses a single flat feature and a bounded smooth Gaussian-density weight.
The target mixes N(0,1/16) and its degree128 even polynomial tilt with
weights113/128 and15/128. Its variance is exactly one and its absolute
mean is approximately .5081999212494555. Binomial summation through
n=16384 checks normalization, Frobenius error and exact covariance repair;
the repaired response upper bound there is approximately .50831641.
These floating diagnostics are not proof of a rate or a statement about
original minimizing matrices. Their parameters and outputs are preserved.
