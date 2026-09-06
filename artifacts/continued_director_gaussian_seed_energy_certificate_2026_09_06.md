# An explicit Gaussian-seeded feedback energy certificate

Date: 2026-09-06. Status: director derivation, submitted for independent audit.
This deliberately changes the RANDOMIZED ALGORITHM: its initial seed is
Gaussian. It does not replace the law of the Boolean-seeded coherent return.
No statement about minimizing spectra or the limiting optimum is assumed.

## 1. Setup and quantitative certificate

Let B=A/sqrt(n-1) be a hollow symmetric signing with fixed operator cap L.
Let xi be an independent standard Gaussian vector. Define

```math
Q=B^2,\quad G=B\xi,\quad V=Q\xi,\quad
b=\mathbb E[Nf(N)],\quad r(t)=f(t)-bt,\quad Z=Br(G),
```

where f is fixed bounded odd and Gaussian-a.e.-continuous. Assume H is
fixed bounded even C2 with its first two derivatives bounded, and psi is
fixed odd C2 with psi,psi',psi'' bounded. Put

```math
R=\sum_{p\ge3\ {\rm odd}} f_p^2 Q^{\circ p},\quad T=BRB,
\quad \sigma_i^2=T_{ii},\qquad C_i=H(G_i)\psi(bV_i+Z_i).
```

For a scalar Gaussian eta independent of the coherent fields, define

```math
c_i^0=H(G_i)\mathbb E_\eta\psi(bV_i+\sigma_i\eta),
\quad a_i=\mathbb E[H(G_i)\psi'(bV_i+\sigma_i\eta)],
\quad u_i=\mathbb E[H'(G_i)\psi(bV_i+\sigma_i\eta)].
```

These are expectations against an EXACT Gaussian law with
Cov(G_i,V_i)=(B^3)_ii and Var(V_i)=(B^4)_ii. No local covariance inverse
is introduced, including when this two-dimensional Gaussian is singular.
Let D_a,D_u be the corresponding diagonal matrices and set

```math
M=D_uB+bD_aQ,\qquad K_{ij}=\mathbb E[c_i^0 r(G_j)].
```

Then the claimed actual expected feedback half-energy is

```math
\boxed{
\frac{\mathbb E C^TBC}{2n}
=\frac{\operatorname{Tr}(BMM^T)}{2n}
 +\frac{\langle BD_aB,K\rangle_F}{n}
 +\frac{\operatorname{Tr}(BD_aTD_a)}{2n}+o(1).
} \tag{1}
```

Here <X,Y>_F=sum_ij X_ij Y_ij. Every K_ij is a fixed-dimensional Gaussian
integral: the variables (G_i,V_i,G_j) have covariance

```math
\begin{pmatrix}
1&(B^3)_{ii}&Q_{ij}\\
(B^3)_{ii}&(B^4)_{ii}&(B^3)_{ij}\\
Q_{ij}&(B^3)_{ij}&1
\end{pmatrix}, \tag{2}
```

and the extra eta is independent. Formula (2) is positive semidefinite
because its variables are literal linear functions of xi. The residual
cross is retained; it is not inferred from one-root Gaussianity.

If H>=0, |f|+H<=1, |psi|<=1 and y psi(y)>=0, define

```math
j_n=\frac1n\sum_i\mathbb E
 H(G_i)(bV_i+\sigma_i\eta)\psi(bV_i+\sigma_i\eta).
```

Both +/-f(G)+C lie in the cube. With e_n the three terms on the right of
(1), the endpoint identity gives the actual-signing lower certificate

```math
\Lambda(B)\ge j_n+
 \left|\frac{b^2\operatorname{Tr}(B^3)}{2n}+e_n\right|-o(1).
\tag{3}
```

This is a quantitative energy certificate, not an assertion that its value
exceeds the established universal bound. All functions, their smoothness
bounds, and L are fixed before n tends to infinity.

## 2. Proof of the reduction

The retained-coherent-return proof in
`continued_feedback_coherent_return_energy_projection_2026_09_06.md`
applies with Gaussian instead of Boolean input. Coherent polynomials now
have ordinary Gaussian chaos decompositions rather than Boolean Walsh
decompositions. Their root kernels are finite products of the bounded-op
rows B,Q; their contractions give bounded row scalars. The global root
maps in every positive degree are bounded by Schur products of their Gram
matrices. No Lindeberg replacement of V is needed or performed.

The nonlinear r(G) branch has only odd degrees at least three. Its proper
flattenings have norm O(n^-1/2). Full contractions with coherent rows are
small locally by the flat outer B factor and column Cauchy--Schwarz.
In the cross-root energy, exactly the bare nonlinear/coherent stars survive;
all terms with an additional bridge have two bounded covariance factors
and vanish against the flat B entries. Thus the same proof gives

```math
\frac{\mathbb E C^TBC}{2n}
=\frac{\mathbb E(c^0)^TBc^0}{2n}
 +\frac{\mathbb E(c^0)^TBD_aZ}{n}
 +\frac{\operatorname{Tr}(BD_aTD_a)}{2n}+o(1). \tag{4}
```

The local Gaussian-noise comparison also identifies the cross gain j_n.
Finite polynomials precede the order limit; bounded approximation follows
as in the audited source, using the smooth H,psi hypotheses stated here.

Gaussian integration by parts gives the first xi-chaos row of c_i^0:

```math
\mathbb E[\xi_k c_i^0]=u_i B_{ik}+b a_i Q_{ik}=M_{ik}. \tag{5}
```

This uses an ambient derivative, not inversion of the local covariance.
For a fixed polynomial coherent response, every higher local chaos has
degree at least three by oddness. Its cross-root covariance consists of
products of at least three entries from Q,B^3,B^4, times bounded row
coefficients. All these covariance matrices have bounded operator norm and
Frobenius norm O(sqrt(n)). Two factors in Frobenius Cauchy--Schwarz, remaining
ones entrywise bounded, and max_(i!=j)|B_ij|=O(n^-1/2) show their normalized
energy is O(n^-1/2). B_ii=0 removes diagonal terms exactly. Therefore

```math
\frac{\mathbb E(c^0)^TBc^0}{2n}
=\frac{\operatorname{Tr}(BMM^T)}{2n}+o(1). \tag{6}
```

Bounded smooth passage follows by Gaussian L2 approximation and contractivity
of projection onto the first global chaos. Equivalently use (5), bounded
derivatives, and dominated convergence on compact regions before removing
tails. The possible singularity of the local Gaussian is harmless.
Finally Z=Br(G), so the middle term in (4) is EXACTLY
<BD_aB,K>_F/n. Equations (4)--(6) prove (1). The scalar old self-energy is
b²Tr(B³)/(2n)+o(1); the exact two-endpoint identity proves (3).

## 3. What this removes, and what it does not

For this intentionally Gaussian-seeded algorithm, every retained energy
term is determined by matrix powers through B^4, the entrywise Hermite
kernel R, and Gaussian integrals of dimension at most four. No Boolean
energy table or arbitrary bridge maximization is used to evaluate it.
This is substantially less than a full landscape, but it is only one
certificate family, not a sufficient statistic for the optimum.

In particular <BD_aB,K> need not be negligible merely because r has no
first Gaussian coefficient: Q or B^3 can have order-one off-diagonal
entries on actual bounded-op signings. Neither its sign nor cancellation
with the other two terms is assumed. A universal quantitative bound on
this explicit matrix certificate is an open next question, not a proof of
convergence. General near-minimizers need not have a fixed operator cap;
using a universal estimate in the original problem would still require
the separately audited principal spectral-deletion step and limit order.

## 4. A completely explicit sine test, without high-dimensional integration

Choose fixed 0<kappa<1, h=1-kappa, f(g)=kappa sin(g), H=h, and
psi(y)=sin(t y), t>0. The endpoint feasibility holds although psi need not
be sign-preserving. The average cross gain below is positive, which suffices:
the exact endpoint identity and Jensen give |E cross|+|E common|.
All matrix functions named sinh here are ENTRYWISE, not spectral functions.

Set b=kappa exp(-1/2) and

```math
R=b^2(\sinh^{\circ}(Q)-Q),\quad T=BRB,\quad
U=T+b^2Q^2=b^2B\sinh^{\circ}(Q)B,\quad w_i=U_{ii},
\quad a_i=h t\exp(-t^2 w_i/2).
```

Here w_i is the total Gaussianized return variance. The first Gaussian
coefficient of the sine is b, and its odd Hermite coefficients are
kappa exp(-1/2)(-1)^((p-1)/2)/sqrt(p!). This proves the displayed R
identity with absolute convergence on |Q_ij|<=1. Formula (2) and the exact
Gaussian sine product identity give

```math
K_{ij}=\frac{b a_i}{t}
 \left[\sinh(t b(B^3)_{ij})-t b(B^3)_{ij}\right]. \tag{7}
```

Since u_i=0, the first and third terms of (1) combine. For the ACTUAL
C=h sin(t B[kappa sin(B xi)]), the claim becomes

```math
\frac{\mathbb E C^TBC}{2n}
=\frac{\operatorname{Tr}(BD_aUD_a)}{2n}
 +\frac{b}{tn}\sum_{ij}(BD_aB)_{ij}a_i
  [\sinh(tb(B^3)_{ij})-tb(B^3)_{ij}]+o(1). \tag{8}
```

The cross gain is exactly identified asymptotically by

```math
j_n=\frac1n\sum_i a_i w_i>0. \tag{9}
```

Consequently (3) holds for this sine choice despite the absence of global
sign preservation. Equations (7)--(9) use only matrix products and scalar
elementwise functions. They are an explicit, finite-dimensional actual
feedback energy test, not an implicit Boolean optimization. They need the
same audited asymptotic projection as (1), and currently supply no improved
universal lower coefficient.
