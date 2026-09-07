# Independent audit: uniformly elliptic Gaussian-sign small-source MGF

2026-09-07. Root construction, independently reconstructed by construct agent. **PASS**, with fixed ellipticity constants. This is an MGF theorem, not by itself a favorable-flatification theorem.

## Statement

Let G be centered Gaussian in any dimension d, with covariance Sigma having diagonal 1 and epsilon I <= Sigma <= K I, where epsilon>0 and K are fixed. Put X=sign(G). Uniformly over real h with a=||h||_infinity sufficiently small,

    log E exp(h.X) = (1/2) h^T Cov(X) h
                      + O_{epsilon,K}(a ||h||_2^2).

Constants do not depend on d, h, or Sigma. In particular d=n^2 and h coordinates of size O(n^(-1/2)) give an O(sqrt(n)) remainder. This does not permit epsilon to approach zero without tracking its constants.

## Smooth Gaussian lemma with explicit error

For standard Gaussian z, let U be smooth with ||grad U||<=L and ||Hess U||op<=H<1 globally. Then

    |log E exp(U) - EU - Var(U)/2|
       <= H L^2 / [2(1-H)].

Let P_s be the Ornstein-Uhlenbeck semigroup and define

    tau(z)=integral_0^infinity e^(-s)
           <grad U(z), P_s grad U(z)> ds.

Gaussian integration by parts and the OU resolvent imply

    E[(U-EU) f(U)] = E[f'(U) tau],   E tau=Var(U).

Moreover |tau|<=L^2. Differentiating the two gradient factors and using grad P_s=e^(-s) P_s grad gives

    Lip(tau) <= HL integral e^(-s)(1+e^(-s)) ds = 3HL/2.

The tilted law dmu_t proportional to exp(tU-|z|^2/2) dz, 0<=t<=1, has potential Hessian at least (1-tH)I. The Brascamp-Lieb variance inequality and Cauchy-Schwarz therefore give

    |Cov_mu_t(tau,U)| <= Lip(tau) L/(1-tH).

For psi(t)=log E exp(tU), the Stein identity gives

    psi'(t)=EU+t E_mu_t tau.

Subtract EU+t E tau, integrate first in the tilt parameter and then in t, and use integral_0^1 t^2 dt=1/3. The displayed error follows. All exponential integrability and differentiation are legitimate because U is Lipschitz; its growth is at most linear. In the application below U is even bounded for each fixed finite d.

Primary research source checked directly: Carlen, Cordero-Erausquin and Lieb, *Asymmetric covariance estimates of Brascamp-Lieb type and related inequalities for log-concave measures*, Ann. IHP Probab. Statist. 49 (2013), equations (1.3)-(1.4), pp. 1-2: https://www.numdam.org/item/10.1214/11-AIHP462.pdf . Their variance bound is Var_mu f <= E <grad f,(Hess potential)^(-1) grad f>. Only this standard bound, not their stronger asymmetric extension, is needed here.

## Conditional-noise reduction

Write G=W+sqrt(epsilon) z', with independent standard Gaussian z' and W of covariance Sigma-epsilon I. This latter covariance may be singular. Write W=Bz with z standard and ||B||op<=sqrt(K); no inverse of B is used.

Given W, the signs are independent with means

    m_i=erf(W_i/sqrt(2 epsilon)).

Uniformly in m in [-1,1] and small real u,

    log(cosh u+m sinh u)
      = u m + (u^2/2)(1-m^2) + O(|u|^3).

Thus define F=sum h_i m_i and V=(1/2)sum h_i^2(1-m_i^2). The conditional log MGF differs from U=F+V by at most C sum |h_i|^3 <= C a ||h||^2, deterministically. This bound passes through the outside expectation and logarithm without a tail event.

Bounded first and second derivatives of erf at this fixed epsilon, and ||B||op<=sqrt(K), give

    ||grad F|| <= C ||h||,   ||Hess F||op <= C a,
    ||grad V|| <= C a ||h||, ||Hess V||op <= C a^2.

Hence the smooth lemma yields

    log E exp(F+V)=E V+(1/2) Var(F+V)+O(a ||h||^2).

F is odd and V is even under z->-z, so EF=0 and Cov(F,V)=0 exactly. Gaussian Poincare gives Var(V)<=C a^2 ||h||^2. Finally conditional variance gives

    Var(h.X)=Var(F)+E sum h_i^2(1-m_i^2)=Var(F)+2EV.

Combining these identities proves the claimed formula.

## Scope checks

* No independence between original coordinates is assumed; their covariance enters through B and the exact final variance.
* No uniform-tail estimate or cumulant summability is being assumed.
* All sources may have arbitrary signs.
* The constant K must be uniform. A diverging covariance operator norm is not covered by the fixed-constant statement.
* Hard singular endpoints (epsilon=0) and optimization over noise removal remain separate problems. The theorem supplies an actual correlated-sign MGF estimate, but no child/bridge payment is asserted here.
