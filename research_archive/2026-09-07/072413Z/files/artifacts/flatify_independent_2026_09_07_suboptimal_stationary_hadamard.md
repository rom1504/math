# Explicit suboptimal signings with stronger-than-required global stationarity

Date: 2026-09-07. Status: exact construction and obstruction for the
stationarity-only mechanism. These matrices are NOT claimed to be global
minimizers; the current all-order upper theorem proves they are eventually
suboptimal. This is not nonconvergence of M_n.

## Exact family

Let n=4^d, d>=1, and index vertices by F_2^(2d). Define

```math
f(z)=(-1)^(z_1 z_2+z_3 z_4+...+z_(2d-1) z_(2d)),
H_xy=f(x+y),   A=H-I.
```

Then A is a hollow full signing. For a character chi_t(x)=(-1)^(t dot x),
the product of d elementary two-bit Walsh transforms gives

```math
H chi_t = sqrt(n) f(t) chi_t.
```

Consequently H^2=nI, H has diagonal 1, and A has eigenvalues
sqrt(n)-1 and -sqrt(n)-1. Both eigenspaces contain Boolean characters.
The spectral upper bound is attained by a negative character, proving

```math
Q(A)=n(sqrt(n)+1)/2,
q=Q(A)/n^(3/2)=1/2+1/(2sqrt(n)).                    (1)
```

Thus the asymptotic cap is 1/2, strictly worse than the proved original
all-order limsup bound <.493608094.

## A single global near-ground law with optimal-scale covariance balance

Choose t uniformly among the n characters, put x=chi_t and sigma=f(t).
Since sum_t f(t)=sqrt(n), E sigma=1/sqrt(n). The signed normalized energy is

```math
sigma H_A(x)/n^(3/2)=1/2-sigma/(2sqrt(n)),
E[sigma H_A(x)/n^(3/2)]=1/2-1/(2n).
```

Every selected signed state is within 1/sqrt(n) of the absolute ground cap.
The mean slack is exactly 1/(2sqrt(n))+1/(2n). Fourier diagonalization gives

```math
C=E[sigma xx^T/n]=H/n^(3/2).
```

All n singular values of C are exactly 1/n. Therefore

```math
||C||op=1/n,  ||C||F=1/sqrt(n),  ||C||_(r)=r/n.
```

This ONE law simultaneously controls every feature subspace. Its covariance
is much smaller than the error bounds currently deduced from genuine global
optimality by sparse finite-rank/nuclear-budget surgery. Thus those covariance
conditions by themselves do not even distinguish asymptotic cap 1/2 from
the strictly better original minimizers.

## Macroscopic nuclear budget is necessary for a leading direct descent

Let G be ANY symmetric matrix with ||G||op<=1 and ||G||*<=r, with no rank
restriction. Evaluate the same signed character law on
A-theta sqrt(n) offdiag(G), theta>=0. Since H/sqrt(n) has operator norm 1,

```math
E[sigma H_(A-theta sqrt(n) offdiag G)(x)]/n^(3/2)
 = 1/2-1/(2n)
   -theta Tr[(H/sqrt(n))G]/(2n)
   +theta Tr(G)/(2n sqrt(n))
 >= 1/2-1/(2n)-theta r(1+1/sqrt(n))/(2n).
```

As a maximum dominates this expectation,

```math
Q(A-theta sqrt(n) offdiag G)/n^(3/2)
 >= q-[1/(2sqrt(n))+1/(2n)
       +theta r(1+1/sqrt(n))/(2n)].                  (2)
```

For fixed theta and r=o(n), no such perturbation can lower the cap by a
fixed positive multiple of n^(3/2). This is substantially stronger than an
o(sqrt(n))-rank obstruction: it covers the entire convex full-rank
nuclear-budget ball of radius o(n).

The scale is sharp. G=H/sqrt(n) has operator norm 1 and nuclear norm n;
its off-diagonal part is A/sqrt(n), so the perturbed matrix is exactly
(1-theta)A. For 0<=theta<=1 its cap decreases by theta Q(A).
This last real-coefficient descent is not itself an actual-sign rounding
theorem: the known bias and variance payments remain essential.

## Scope

There is also a direct ACTUAL edit-density consequence, observed by the root.
Let C be any full signing obtained from A by flipping exactly E unordered
edges. Since each off-diagonal H_ij equals A_ij, each flip decreases
Tr(C H) by 4. The same character law therefore yields

```math
Q(C)>=Tr(C H)/(2sqrt(n))
     =(n^(3/2)-sqrt(n))/2-2E/sqrt(n).                (3)
```

In particular, Q(C)<=(1/2-delta)n^(3/2), delta>0, requires

```math
E >= delta n^2/2-n/4.                               (4)
```

Thus a fixed improvement below 1/2 needs a positive density of actual edge
changes, not merely a large formal perturbation rank. Equation (4) concerns
distance from this explicit Hadamard family; it is not a lower bound on the
distance between arbitrary optimal signings at different orders.

This example does not refute a theorem that uses additional global
optimality information. It does show that even dramatically sharpening the
present covariance/nuclear stationarity errors cannot, on its own, supply
the desired cross-order mechanism. A decisive use of global optimality must
obtain information absent from this explicit suboptimal family, or use an
operation with genuinely macroscopic nuclear budget.
