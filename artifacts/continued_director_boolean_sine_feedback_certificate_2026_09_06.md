# An explicit actual Boolean-seeded feedback certificate

Date: 2026-09-06. Status: algebra independently checked; proof dependencies
are the audited coherent-return energy projection and higher-Walsh energy
reduction. This note does not Gaussianize the coherent return QS.

## 1. Statement

Let B=A/sqrt(n-1) be an arbitrary hollow symmetric signing with a fixed
operator cap L, and let S have independent unbiased sign coordinates.
Fix kappa,h>0 with kappa+h<=1, and t>0. Set

```math
G=BS,\quad f(g)=\kappa\sin g,\quad b=\kappa e^{-1/2},
\quad r(g)=f(g)-bg,\quad Z=Br(G),\quad V=B^2S,
\quad C=h\sin(tBf(G)).
```

These are literal, actual vectors, not an auxiliary Gaussian simulation.
Write Q=B² and

```math
R=b^2(\sinh^{\circ}Q-Q),\quad T=BRB,\quad
\sigma_i^2=T_{ii},\quad d_i=h e^{-t^2\sigma_i^2/2},
\quad \theta_{ik}=tbQ_{ik}.
```

The symbol sinh^circ denotes ENTRYWISE sinh. Define finite products

```math
P_i=\prod_k\cos\theta_{ik},\qquad
L_{ik}=\sin\theta_{ik}\prod_{\ell\ne k}\cos\theta_{i\ell},
\quad M_{ik}=d_iL_{ik},\quad a_i=t d_iP_i,
```

and the n-by-n matrix

```math
K_{ij}=\frac{\kappa d_i}{2}
 \left[\prod_k\cos(\theta_{ik}-B_{jk})
       -\prod_k\cos(\theta_{ik}+B_{jk})\right]
       -b(MB)_{ij}. \tag{1}
```

Then the actual feedback energy has the asymptotic formula

```math
\boxed{
\frac{\mathbb E C^TBC}{2n}
=e_n+o(1),\qquad
e_n=\frac{\operatorname{Tr}(BMM^T)}{2n}
 +\frac{\langle BD_aB,K\rangle_F}{n}
 +\frac{\operatorname{Tr}(BD_aTD_a)}{2n}.
} \tag{2}
```

The actual cross gain satisfies

```math
\frac{\mathbb E f(G)^TBC}{n}=j_n+o(1),\qquad
j_n=\frac1n\sum_i d_i
 \left[b\sum_k Q_{ik}L_{ik}+t\sigma_i^2P_i\right]. \tag{3}
```

In particular the original-signing objective obeys

```math
\Lambda(B)\ge |j_n|+
 \left|\frac{b^2\operatorname{Tr}(B^3)}{2n}+e_n\right|-o(1).
\tag{4}
```

If 0<tb<pi/2, then j_n>0, so its absolute value can be removed.
All parameters and L are fixed before the matrix order tends to infinity.
No universal improvement over .4333221116640807 is claimed.

## 2. Proof and exact dependencies

The residual sine's Gaussian Hermite coefficients give
R=sum_(p>=3 odd)f_p²Q^(circ p)=b²(sinh^circQ-Q). The general actual
retained-return theorem gives the same energy as c0+D_a Z, where

```math
c_i^0=h\mathbb E_\eta\sin(t[bV_i+\sigma_i\eta])
      =d_i\sin(tbV_i),\quad
a_i=h t\mathbb E_{S,\eta}\cos(t[bV_i+\sigma_i\eta]).
```

Gaussian eta is a comparison for the ZERO-FIRST noise only; V=QS remains
Boolean-linear. Its exact characteristic function is the finite product
P_i. Differentiating one seed factor gives

```math
\mathbb E[S_k c_i^0]=d_i\sin(tbQ_{ik})
                    \prod_{\ell\ne k}\cos(tbQ_{i\ell})=M_{ik}.
```

The coherent-Walsh energy theorem proves that all higher Boolean Walsh
degrees of c0 contribute o(1) to its normalized flat-B energy; hence its
energy is Tr(BMM^T)/(2n)+o(1). This does NOT identify its full covariance
or erase its higher-degree variance.

The retained cross is E[(c0)^T B D_aB r(G)]/n. The elementary identity
sin u sin v=[cos(u-v)-cos(u+v)]/2 and independent sign inputs compute
E[c_i^0 kappa sin(G_j)] as the two products in (1). Also
E[c_i^0 G_j]=(MB)_ij exactly. Thus K_ij=E[c_i^0 r(G_j)] and (2) follows.

The tested one-root comparison gives (3) from
E[h(bV_i+sigma_i eta)sin(t[bV_i+sigma_i eta])]. Equivalently it is
-h times the t derivative of exp(-t²sigma_i²/2)P_i. Here sigma_i is
independent of the chosen outer frequency t. Since |Q_ik|<=1, the
condition 0<tb<pi/2 makes every cosine positive, and every summand
Q_ik sin(tbQ_ik) nonnegative. The diagonal Q_ii=1 makes j_n positive.

Both f(G)+C and -f(G)+C lie in the cube because kappa+h<=1. The exact
two-endpoint identity, expectation, and Jensen yield |E cross|+|E common|.
The known old sine self-energy is b²Tr(B³)/(2n)+o(1), proving (4).
Sine is not globally sign-preserving; the explicit sign of j_n replaces
that optional sufficient hypothesis. No false pointwise sign is used.

## 3. Computational and conceptual scope

All terms in (1)--(4) use matrix multiplication and O(n³) scalar products
and elementary function evaluations. Products excluding one index can be
computed by prefix/suffix products, without division by possibly zero
cosines. Thus no 2^n seed sum, full Boolean optimization, or fitted latent
distribution is required to evaluate the certificate. Arithmetic accuracy
must be controlled separately if a finite numerical value is certified;
the o(1) error in the asymptotic theorem is not a sampling error.

This removes a concrete actual-feedback evaluation obligation. It does not
yet prove that these explicit matrix terms force a stronger universal cap,
give a subsequent mixed-query law, or establish a cross-order theorem.
The higher-dimensional colored-seed theorem remains separate from an
unproved extension to arbitrary old nonlinear tree fields.

Proof dependencies:

- `continued_feedback_coherent_return_energy_projection_2026_09_06.md`;
- `continued_audit_coherent_return_energy_projection_2026_09_06.md`;
- `continued_audit_coherent_walsh_energy_reduction_2026_09_06.md`.

The sine-product algebra was independently reconstructed by the audit agent.
