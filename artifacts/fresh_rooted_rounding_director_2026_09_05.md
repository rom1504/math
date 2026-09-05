# Rooted rounding: independent director derivation

Status: **independent derivation, subsequently audited**. The canonical
rooted Gaussian proof and two independent audits establish the universal
lower bound 0.3396496212118657. This note records the director's separate
derivation. The original convergence problem is not solved.

## 1. Normalization and the concrete target

Let A be a hollow symmetric signing, m=n-1, B=A/sqrt(m), and C=B^2.
Write q(B)=max_x |x^T Bx|/2. We only need consider q(B)<=Kn, for a fixed
K, because larger caps already beat the current lower bound.

There is an elementary improved spectral bootstrap. If Av=lambda v,
choose independent signs with mean v/||v||_infinity. Then

```math
\|A\|_{\infty\to1}\ge
|\lambda|\frac{\|v\|_1}{\|v\|_\infty}\ge\lambda^2,
```

where the second inequality follows from the eigenvector equation and
|a_ij|<=1. Ordinary polarization, used only for this bootstrap, gives
||A||_(infinity->1)<=4q(A). Consequently

```math
\|B\|_{op}^2=O_K(\sqrt n),\qquad
\operatorname{Tr}B^4=O_K(n^{3/2}),\qquad C_{ii}=1.
```

Take independent Rademachers S and G=BS. For an odd smooth response psi,
define

```math
f(g)=\frac{\psi(g+t)+\psi(g-t)}2,\quad
h(g)=\frac{\psi(g+t)-\psi(g-t)}2,
\quad u_i=f(G_i),\quad v_i=S_i h(G_i).
```

Thus the two smoothed initial means are mu_sigma=sigma u+v. The proposed
lower-response target is

```math
\frac1n\mathbb E\sum_i\max(|(Bu)_i|,|(Bv)_i|)
\ge 2c_*+\delta
```

for sufficiently small fixed Gaussian dither, with one fixed delta>0,
where c_*=0.336493364431... is the current undoubled constant. The paired
initial-energy limit tends to c_*. A partial best-response update then
turns a fixed response excess into a fixed, possibly tiny, energy gain.

## 2. Exact Gaussian identities for the rooted channel

For this section only, Z is standard Gaussian, G=BZ,
V_j=Z_j h(G_j), Y=BV, b=E h(G_0). Hollowness makes Z_j independent G_j.
For smooth even h, two Gaussian integrations by parts give exactly

```math
\operatorname{Cov}(V_j,V_k)=
\begin{cases}
\mathbb Eh(G_0)^2,&j=k,\\
m^{-1}\mathbb E[h'(G_j)h'(G_k)],&j\ne k.
\end{cases}
```

Let K_h be the Gaussian covariance kernel of h'. It is PSD. Its Hermite
expansion is a positive combination of odd Schur powers of C, and hence
||K_h||op<=E h'(G_0)^2 ||C||op (Schur multiplication by a correlation
matrix is a unital positive contraction in operator norm). It follows that

```math
\mathbb EY_i^2=\mathbb Eh(G_0)^2+o(1),\qquad
\mathbb E[Y_iG_i]=b,
```

uniformly in i for fixed h. Thus R_i=Y_i-bG_i has variance tending to
E h^2-b^2. The same calculation for h minus a fixed polynomial is the
prospective L2 tail-transport estimate; it does not pay ||B||op times an
uncontrolled approximation error.

Let D_i denote Gaussian directional differentiation along row B_i. Then

```math
D_i(Y_i-bG_i)
=\frac1m\sum_{j\ne i}h(G_j)-b
 +\sum_j B_{ij}C_{ij}Z_jh'(G_j).
```

The first term has variance O_h(Tr(C^2)/m^2)=o(1), because h is even.
For the second, the covariance matrix of Z_j h'(G_j) has diagonal
E h'^2 and off-diagonal m^-1 E h''(G_j)h''(G_k). Bounding the latter
kernel entrywise by E h''^2 gives a variance bound

```math
O_h\left(\sum_j B_{ij}^2 C_{ij}^2\right)
=O_h((B^4)_{ii}/m)=o(1).
```

Conditional Gaussian Poincare therefore makes R_i L2-close to a function
of the projection of Z orthogonal to B_i. That function is independent
of G_i. This is an independence conclusion, not a claim that R_i is
Gaussian.

## 3. A cubic witness prevents residual mass from escaping

Put F_i=(B W)_i, W_j=Z_j H_2(G_j), H_2(g)=g^2-1. Since Z_j is
independent G_j, F_i is pure third Gaussian chaos. The preceding covariance
calculation gives

```math
\mathbb EF_i^2=2+o(1),\qquad
\mathbb E[R_iF_i]=d+o(1),\qquad d=\mathbb E[h(G_0)H_2(G_0)].
```

More precisely the latter covariance is
d(1-2/m+2(B^4)ii/m). Hypercontractivity gives
E F_i^4<=3^6(E F_i^2)^2. For thresholds near t_*, d tends to
-2t_* phi(t_*)=-0.4763..., and remains bounded away from zero.

Together with E R_i^2<=1+o(1), this forces a fixed positive mass of
|R_i|>=r_0 for some fixed r_0>0. Indeed

```math
|\mathbb E R_iF_i|
\le r_0(\mathbb EF_i^2)^{1/2}
 +(\mathbb ER_i^2)^{1/2}
   (\mathbb EF_i^4\Pr(|R_i|>r_0))^{1/4}.
```

This avoids assuming fourth-moment boundedness of the full nonlinear
residual. A fixed-degree witness supplies the needed uniform integrability
information operationally.

## 4. Why this would force a strict response gain

Leaving S_i unreplaced, Taylor expansion in its direct influence on G_j
should give

```math
(Bu)_i=a S_i+W_i+o_{L^2}(1),\qquad a=\mathbb Ef'(G_0),
```

with (W_i,(Bv)_i) asymptotically independent of S_i. The coefficient is
the average m^-1 sum_j f'(G_j), which concentrates by the row-Gram bound.
The direct S_i influence on Bv is m^-1 sum_j S_j h'(G_j), which has
vanishing variance. These leave-one-out steps require a separate
Rademacher audit.

Conditional on all other variables, averaging S_i gives

```math
\frac12\max(|a+W|,|Y|)+\frac12\max(|-a+W|,|Y|)
\ge\max(a,|Y|).
```

If the rooted Gaussian comparison in Section 2 transfers to Rademachers,
the right side asymptotically has Y=bG+R, with R independent G. Define
j(r)=E_G max(a,|bG+r|). It is even and strictly convex. The cubic witness
forces E j(R)>=j(0)+delta for a uniform delta>0. At zero dither and
t=t_*, a=t_* b and j(0)=2ab=2c_*.

## 5. Proof obligations being checked before promotion

1. Justify the Rademacher-to-Gaussian comparison for the rooted channel,
   jointly with G_i, while retaining S_i as a Rademacher.
2. Prove L2 polynomial tail transport in the Rademacher model by the
   direct-coupling Taylor expansion, not by a growing operator norm.
3. Check all smooth approximation and first-n-then-dither limits uniformly
   over the original low-cap sign matrices.
4. Audit the paired actual-Boolean best-response update and its relation
   to the smoothed-mean response lower bound. Jensen in the independent
   initial dithers should only improve the latter, but the normalization
   must be explicit.

For each fixed polynomial, its multilinear rooted chaos coefficients have
maximum coordinate influence O_D(1/n), because all B entries have magnitude
1/sqrt(m). Repeated-index terms have averaged L2 error O_D(1/n) before
transport; the bootstrap ||B||op^2=O(sqrt n) still makes their transported
normalized error vanish. This is the proposed finite-degree invariance
mechanism. It is being independently checked rather than assumed.

No new universal constant should be recorded until these items are closed.
