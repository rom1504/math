# Exact-covariance Gaussian-sign MGF: independent reconstruction

2026-09-07. Root's smoothing/Stein proposal independently passes. This removes the Holder inflation under a FIXED positive lower spectral bound; it does not establish the original child composition.

## 1. A dimension-free smooth Gaussian lemma

Let g be standard Gaussian in any finite dimension and U be smooth with

    sup ||grad U||<=L,  sup ||Hess U||op<=b<1.

Then

    |log E exp U - E U - Var(U)/2|
       <= b L²/[2(1-b)].                              (1)

Let P_s be the standard Ornstein--Uhlenbeck semigroup and define

    tau(g)=grad U(g) dot integral_0^infinity e^(-s) P_s grad U(g) ds.

Gaussian integration by parts and the semigroup resolvent identity give E tau=Var U and

    E[(U-EU) phi(U)]=E[tau phi'(U)].

Differentiating the displayed tau, using grad P_s=e^(-s)P_s grad, gives the uniform bound

    Lip(tau)<=bL integral(e^(-s)+e^(-2s))ds=(3/2)bL.

For t in [0,1], the probability density proportional to exp(tU(g)-||g||²/2) has potential Hessian at least (1-b)I. It consequently has Poincare constant at most (1-b)^(-1). For completeness, this follows directly from its Langevin semigroup: the derivative flow contracts by exp(-(1-b)s), so the gradient contraction and variance dissipation identity give Var(f)<=E||grad f||²/(1-b). Smoothness and the bounded derivatives here justify the semigroup argument; standard smooth truncation also gives the identity for the Lipschitz functions used below.

Writing E_t for this tilt, Cauchy--Schwarz and that Poincare estimate yield

    |d E_t tau/dt|=|Cov_t(tau,U)|
       <=(3/2)bL²/(1-b).

The Stein identity gives (log E exp(tU))'=EU+t E_t tau. Integrating the preceding derivative bound first from 0 to t, then integrating t from 0 to 1, proves (1). No independence of U's summands is required.

## 2. Smoothing the actual sign law

Suppose R is a correlation matrix with epsilon I<=R<=KI for fixed positive epsilon,K. Write

    G=W+sqrt(epsilon)z, W=(R-epsilon I)^(1/2)g,

with independent standard Gaussian g,z. Put f(w)=2Phi(w/sqrt(epsilon))-1. Conditional on W, the sign coordinates are independent with means f(W_i). Their scalar log-MGFs have the uniform expansion

    log E_z exp(sum h_i sign G_i)
      =F(g)+V(g)+E_h(g),
    F=sum h_i f(W_i),
    V=(1/2)sum h_i²[1-f(W_i)²],
    |E_h(g)|<=C sum |h_i|³<=C ||h||infinity ||h||².

Uniformity includes conditional means arbitrarily close to +/-1: the scalar third derivative is a third centered moment of a variable bounded by one, hence universally bounded.

For U=F+V, bounded first and second derivatives of f and ||R-epsilon I||op<=K imply, when ||h||infinity<=1,

    ||grad U||infinity,Euclidean <=C_(epsilon,K)||h||,
    ||Hess U||infinity,op <=C_(epsilon,K)||h||infinity.

Thus (1) applies for sufficiently small source sup norm and has remainder O_(epsilon,K)(||h||infinity ||h||²). The pointwise conditional remainder changes the log expectation by at most the same order.

Crucially F is odd in g and V is even. Therefore EF=0 and Cov(F,V)=0 EXACTLY. Also the ordinary Gaussian Poincare inequality gives

    Var V<=C_(epsilon,K)||h||infinity² ||h||².

The conditional-variance identity supplies

    Var(sum h_i sign G_i)=Var F+2 EV.

Combining these facts proves the two-sided approximation

    log E exp(sum h_i sign G_i)
      =(1/2)h^T Cov(sign G)h
         +O_(epsilon,K)(||h||infinity ||h||²).          (2)

No small-coordinate bound on the conditional or tilted means is used. Applying the smooth lemma directly to F+V avoids an extra tilted perturbation argument for V.

## 3. Opposite-spectral bridge scope

For flat children A²=D²=(n-1)I and fixed 0<=rho<1,

    R=I-rho A tensor D/(n-1)

has lower spectral bound 1-rho and upper bound 1+rho. Its exact sign covariance is I-kappa A tensor D, where kappa=(2/pi)arcsin(rho/(n-1)). With h=t(x tensor y)/sqrt(n), formula (2) has error O_(rho,t)(sqrt(n)) and quadratic coefficient

    (nt²/2)[1-4n kappa e_A e_D]
      =(nt²/2)[1-(8rho/pi)e_Ae_D+O(1/n)].

This is the actual full exponential estimate with exact covariance, not a variance heuristic. The fixed-parameter argument permits first taking n to infinity at fixed rho<1, then letting rho increase to one. The following explicit bookkeeping also permits a particular polynomial approach; it does not assert the singular law itself satisfies a uniform error bound.

## 4. Paid polynomial approach to the singular covariance

For 0<epsilon<=1, the scalar smoothing function obeys

    sup |f'|<=C epsilon^(-1/2),
    sup |f''|<=C epsilon^(-1).

Consequently, for ||h||infinity<=1,

    L<=C sqrt(K/epsilon)||h||,
    b<=C(K/epsilon)||h||infinity.

If ||h||infinity<=c epsilon/K, the smooth lemma's remainder is therefore at most

    C K² epsilon^(-2)||h||infinity ||h||².

The conditional scalar remainder has an absolute constant. Also Var V<=C K epsilon^(-1)||h||infinity²||h||², which is absorbed by the preceding bound. Thus this displayed dependence on epsilon is valid for the full exact-covariance approximation.

For the flat bridge, take epsilon_n=n^(-1/6) and rho_n=1-epsilon_n. At every bounded fixed Chernoff source t, the remainder is O_t(n^(5/6)). The difference between its exact covariance coefficient and the limiting proxy 1-(8/pi)e_A e_D is O(epsilon_n+1/n), uniformly in Boolean x,y by the flat spectral bound. This costs another O_t(n^(5/6)) in the log MGF. The Hessian condition is satisfied since b=O_t(n^(-1/3)). Hence this actual nonsingular sign ensemble has

    log E exp[t x^T C y/sqrt(n)]
      =(nt²/2)[1-(8/pi)e_A e_D]+O_t(n^(5/6)),

uniformly in the Boolean pair. This supplies an explicit power-saving MGF approximation to the singular proxy; it does not eliminate the shell-entropy or actual-child spectral hypotheses.

The result still requires payment of every relevant child-energy shell. Arbitrary actual minimizers need not satisfy the flat spectral hypothesis, and their repaired opposite-spectral covariance need not have bounded operator norm. These remain separate construction obligations.

### A simpler spectrally normalized law valid for ANY children

For arbitrary hollow sign A,D, including actual minimizers, let p=||A||op||D||op and use instead

    R_rho=I-rho A tensor D/p,  0<=rho<1.

This is always a correlation matrix with spectrum in [1-rho,1+rho]. Its sign covariance is exactly

    I-kappa A tensor D,  kappa=(2/pi)arcsin(rho/p),

because all nonzero off-diagonal entries of the tensor are +/-1. Therefore the same theorem, without any additional spectral premise, gives

    log E exp[t x^T C y/sqrt(n)]
      =(nt²/2)[1-4n kappa e_Ae_D]+O_(rho,t)(sqrt(n)).

The polynomial approach rho_n=1-n^(-1/6) remains valid with O_t(n^(5/6)) error. The price is in the reward, not an unpaid covariance response: its leading energy coefficient is 8rho n/(pi p). Since p>=n-1, flat children maximize this normalization; the current general cap bound only gives p=O(n^(3/2)), so this reward need not stay bounded away from zero. The variance proxy is uniformly positive: |(x^T A x)(y^T D y)|<=p n² and hence 4n kappa |e_Ae_D|<=p kappa<=rho for p>=1. For growing n it is in fact at most 2rho/pi+o(1).

This is an actual energy-sensitive sign bridge law for every pair of children. It neither proves a favorable parent cap nor solves the central-entropy problem below.

## 5. The direct scalar energy-bin union has a central-entropy barrier

For every hollow full-sign child, uniform Boolean spins satisfy E H_A²=binom(n,2). Consequently, for each fixed delta>0, at least a fraction 1-O(1/(delta²n)) of all spins have |H_A|<=delta n^(3/2). Both children's central energy window therefore has joint entropy 2 log 2 per n, regardless of optimality.

The exact covariance proxy in this window tends to one as delta tends to zero. A direct scalar Chernoff union at parent coefficient L can only certify the central window if, in that limit,

    L²/2>2 log 2,

requiring L>=sqrt(4 log 2)=1.6651... . The desired balanced coefficient 2sqrt(2) times the current child upper coefficient is below 1.397. Thus the straightforward energy-bin MGF/union certificate cannot close even after removing Holder inflation. This is not an impossibility result for the correlated sign law: a structural bound on its central maximum could replace the union bound. The exact exponential theorem and its genuine energy-sensitive gain remain valid.
