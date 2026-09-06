# Tested coherent return on the first marked history

Date: 2026-09-06. Director derivation submitted for independent audit.
This is a weak, bounded-test regression, not a finite-order conditional-law
or conditional-L² statement. Its gain consequence recovers a known scoped
response lower bound; no new universal coefficient is claimed.

## 1. Exact fields and statement

Let B=A/sqrt(n-1) be an actual symmetric hollow signing, with fixed
operator bound L. Put Q=B², use independent signs S, and define

```math
G=BS,\qquad D_i=S_i h_2(G_i),\qquad Y=BD,\qquad V=QD,
\qquad \gamma_i=(B^3)_{ii}.
```

For every fixed bounded continuous test phi on {−1,1} times R²,

```math
\frac1n\sum_i\left|\mathbb E\left[
 ((QS)_i-S_i-\gamma_i G_i)\,\phi(S_i,G_i,Y_i)\right]\right|\to0, \tag{1}
```

and

```math
\frac1n\sum_i\left|\mathbb E\left[
 (V_i-D_i-\gamma_iY_i)\,\phi(S_i,G_i,Y_i)\right]\right|\to0. \tag{2}
```

The same assertions hold for any fixed polynomial test, and for finite
catalogs of tests with uniformly bounded coefficients. Along a subsequence
where the local parameters converge, they identify the conditional mean in
the limiting local experiment. They do not allow tests changing with n at
unbounded resolution, nor identify the remaining conditional distribution.

## 2. Source-degree proof

The established marked local law makes (S_i,G_i,Y_i) asymptotically a
Boolean sign and two independent standard Gaussians. Its finite Hermite
replacement keeps every source label distinct. In particular

```math
S_i^\delta h_p(G_i)h_q(Y_i),\qquad \delta\in\{0,1\},
```

is, up to a fixed-degree local L² error tending to zero, a homogeneous
source-Walsh tensor of degree delta+p+3q. The own label i can be separated
from Y by its O(1/n) influence; G_i never uses S_i. Partial contractions
between old branches are small by the old proper-cut estimates. These are
the same finite Hermite replacements used in the audited local-noise theorem.

Since QS has degree one, its only nonzero leading matches are S_i and G_i.
Their exact inner products are Q_ii=1 and (QB)_ii=gamma_i. All other
source degrees are orthogonal. This proves (1) on Hermite tests.

The field V has degree three exactly. The only degree-three old tests are
Y_i, S_i h2(G_i)=D_i, and h3(G_i). Write m=n-1. The exact covariance
Cov(D)=(1-3/m)I+2Q/m gives

```math
\mathbb E[V_iY_i]=\gamma_i+O_L(1/n),\qquad
\mathbb E[V_iD_i]=1+O_L(1/n).
```

The third candidate is negligible. If R3_i is the top Boolean component
of h3(G_i), exact enumeration of the common marked triple gives

```math
\mathbb E[R3_iD_j]
=\sqrt3 B_{ij}\left(Q_{ij}^2-\frac{n-2}{(n-1)^2}\right).
```

Therefore, using max|B_ij|=(n-1)^(-1/2), the bounded Q row ℓ² norm,
and |Q_ij|<=1,

```math
|\mathbb E[V_i h_3(G_i)]|
\le C_L n^{-1/2}
 \left(\sum_j|Q_{ij}|^3+n^{-1}\sum_j|Q_{ij}|\right)
 +O_L(n^{-1/2})=O_L(n^{-1/2}).
```

The small lower Boolean component of h3(G_i) is orthogonal to V in fact;
the displayed harmless error also covers the local replacement convention.
This proves (2) on all Hermite tests. The comparison expressions S+gamma G
and D+gamma Y have exactly the matching local moments just listed.

Fixed polynomial moments are uniformly integrable. Gaussian L² approximation
in the old (G,Y) marginal, bounded old-test truncation, and Cauchy–Schwarz with
the uniform second moments of QS,QD extend the identities to the stated
bounded tests. Finite approximation is chosen before matrix order tends to
infinity. The conclusion is averaged as stated; no hidden finite-order
conditional-expectation convergence is used.

## 3. What this licenses for an actual feedback gain

Take fixed bounded odd f(G,Y), bounded nonnegative even H(G,Y), and
|f|+H<=1, with the regularity required by the local-noise theorem. Let
b0=E[G f(G,Y)], b1=E[Y f(G,Y)], and let tau²>0 be the squared norm of
the remaining odd Gaussian Hermite part. The actual field is

```math
Bf(G,Y)=b_0QS+b_1QD+Z.
```

The audited local-noise theorem separates Z from the literal coherent
return, with variance sigma_i²=T_ii. Combining its bounded-tested version
with (1)–(2), uniform integrability, and conditional Jensen in any limiting
local experiment gives

```math
\frac1n\sum_i\mathbb E[H(G_i,Y_i)|Bf(G,Y)_i|]
\ge\frac1n\sum_i\mathbb E H(G,Y)
\left|S[b_0+b_1h_2(G)]+\gamma_i[b_0G+b_1Y]+\sigma_iN\right|-o(1). \tag{3}
```

The expectations on the right use independent S,G,Y,N. Equivalently one
can prove (3) by fixed smooth subgradients and the tested identities, then
remove their smoothing; it is not necessary to assert finite-n conditional
Jensen against an unproved exact regression.

For each fixed G,Y, the law of S[b0+b1h2(G)]+sigma_i N is symmetric.
Adding the deterministic shift gamma_i[b0G+b1Y] cannot decrease its expected
absolute value. Finally t maps to E H E_N|b0+b1h2(G)+tN| is convex and
nondecreasing on t>=0. The odd transported covariance theorem gives
avg sigma_i>=tau. Consequently

```math
\liminf\frac1n\sum_i\mathbb E[H(G_i,Y_i)|Bf(G,Y)_i|]
\ge\mathbb E H(G,Y)\,\mathbb E_N|b_0+b_1h_2(G)+\tau N|. \tag{4}
```

The exact feasible endpoint identity turns (4) into an original-signing
lower certificate. This recovers the corresponding two-coordinate creation
certificate, rather than improving the banked universal .4333221116640807.
The information added here is the explicit gamma-dependent coherent drift
in (1)–(3), valid beyond involutions. Cross-root feedback energy is still
a separate obligation, and (4) is not an all-order convergence mechanism.
