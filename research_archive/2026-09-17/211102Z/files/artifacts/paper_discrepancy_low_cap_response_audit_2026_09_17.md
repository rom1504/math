# Independent audit: uniform cheap sign response on complete low-cap high-energy codes

2026-09-17. Independent reconstruction of the director's
[low-cap uniform-response closure](paper_director_low_cap_uniform_response_2026_09_17.md).
This audit supplies an explicit rational constant certificate and checks
the exact isotropy, subGaussian proxy, minimax quantifiers, full physical
support, and both energy polarities. It does not claim a parent-extension
slope below 3c/2 or a new bound for the original optimization problem.

## 1. A fully numerical corollary

For all sufficiently large n, ANY hollow symmetric full signing A with

```
Q(A)<=n^(3/2)/2
```

admits a centered physical sign law nu with

```
E_nu hh^T=I,
E_nu exp(t.h)<=exp(3||t||_2^2/4) for EVERY t in R^n,
nu(h)>0 for EVERY h in {+-1}^n,
sup_(|H_A(x)| >= (3/10)n^(3/2)) E_nu |h.x|/sqrt(n)
                  <=sqrt(2/pi)-2^(-67).             (1)
```

One law protects the ENTIRE displayed code. It is not a different law
chosen after the query, a restriction to one energy polarity, or a
Gaussian substitute for physical columns. The tiny explicit gap is a
proof-of-mechanism constant, not a quantitatively sufficient recurrence
constant. The order cutoff is supplied by uniform analytic asymptotics;
this arithmetic certificate does not compute that cutoff.

## 2. Half-range input and the finite geometric inequalities

For a principal signing B, let W(B)=(max H_B-min H_B)/2. The archived
[minimum-width theorem](decisive_audit_certified_minimum_width_lower_2026_09_07.md)
proves the universal lower bound

```
liminf_m min_B W(B)/m^(3/2) >= .4333221116640807.
```

I read its complete proof and the complete linked
[fresh lower-chain reconstruction](decisive_audit_fresh_full_lower_chain_2026_09_07.md).
The reason the same lower applies to half-range is substantive: the
marked feasible mean vectors give an ENERGY DIFFERENCE bounded by the
whole Boolean energy range. Principal half-range monotonicity is exact;
the retained-operator cutoff is fixed before the matrix-order limit.
The present audit imports that analytic/numerical lower theorem; it does
not independently rerun the large interval certificate. For (1), only
the weaker eventual bound W(B)>=(43/100)m^(3/2) is used.

Three finite geometric facts are important.

First, for any hollow B and Boolean u,v,

```
u^T Bv=2[H_B((u+v)/2)-H_B((u-v)/2)],
beta(B):=max_(u,v) |u^T Bv| <=4W(B).                (2)
```

Both cube energies lie in the Boolean energy interval by independent
rounding. The factor four is not replaced by an absolute-cap bound.

Second, if S and T partition the coordinates, reversing all spins in
one block reverses only the bridge. This proves separately that each
positive or negative endpoint of A dominates the sum of the corresponding
block endpoints, and hence

```
Q(A)>=W(A)>=W(A_SS)+W(A_TT).                         (3)
```

Third, the rectangle bilinear norm obeys the sharper inequality

```
beta(A_ST)<=Q(A),                                   (4)
```

because for fixed block words the two full energies are H_S+H_T plus
or minus their bridge value. Their maximum absolute value dominates
the absolute bridge value. This rectangle bound improves constants,
but the director's looser factor-four version also closes the theorem.

## 3. Localize an arbitrary query law, not a selected convenient one

Fix ANY law mu supported on the full code in (1). Use the notation of
the [signed-energy alternative](paper_discrepancy_signed_energy_concentration_alternative_2026_09_17.md):

```
Sigma=E_mu xx^T, K=E_mu sign(H_A(x))xx^T,
K0=K-(E_mu sign H_A)I,
Q=1_(Sigma>2), P=I-Q, R=Q Sigma Q,
M=tr R, B=PK0P, b=||B||_F^2/n,
J=E_mu sign(H_A(x))(Qx)^T A(Qx).
```

Because Q is a spectral projector of Sigma, 0<=R<=Sigma and R_ii<=1.
If M<alpha*n, choose S={i:R_ii>epsilon} and T=S^c. Then
|S|<alpha*n/epsilon. The vectors Y_i=(Qx)_i in L2(mu), and the marked
vectors sign(H_A(x))Y_i, have norms at most one, and at most
sqrt(epsilon) when i belongs to T. Bilinear Grothendieck, (2), and (4)
therefore give

```
|J_SS|<=4K_G W(A_SS),
|J_ST|<=K_G Q(A)sqrt(epsilon),
|J_TT|<=4K_G Q(A)epsilon.                            (5)
```

Here J=J_SS+2J_ST+J_TT. There is no extra factor two for the signed
mark: it is included in the second vector family in the SAME bilinear
Grothendieck application.

The Frobenius identity from the signed-energy alternative gives

```
J/n^(3/2)>=2c0-sqrt(b)-(5/sqrt(2))sqrt(alpha),
c0=3/10.                                            (6)
```

Combining (3), (5), (6), Q(A)/n^(3/2)<=C, and the eventual half-range
lower ell yields the necessary inequality

```
C >= ell(1-alpha/epsilon)^(3/2)
     +[2c0-sqrt(b)-(5/sqrt(2))sqrt(alpha)]/(4K_G)
     -C[sqrt(epsilon)/2+epsilon].                    (7)
```

Thus sufficiently small covariance mass AND sufficiently small residual
signed Frobenius energy are incompatible with the actual cap budget.
This is what closes the earlier thin-subspace alternative.

## 4. Exact arithmetic and a uniform response gap

Use the following safe constants:

```
C=1/2, c0=3/10, ell=43/100, K_G<=2,
epsilon=2^(-14), alpha=2^(-23), b0=2^(-16).
```

The Grothendieck bound K_G<=2 follows already from the archived explicit
Krivine construction K_G<=pi/(2 asinh(1)): the elementary integral
bound asinh(1)>=5/6 and pi<22/7 give K_G<66/35<2.

If M<alpha*n and b<b0, equation (7), with
(1-z)^(3/2)>=1-3z/2, forces

```
Q(A)/n^(3/2)
 >= .43+.3/4
    -129/102400-1/2048-5/32768-65/32768
 =205257/409600
 =1/2+457/409600 >1/2,                              (8)
```

a contradiction. The four subtractions pay respectively coordinate
deletion, low signed Frobenius energy, marked cross/centering, and
the Grothendieck cross/tail blocks. Therefore every query law has
M>=alpha*n or b>=b0.

In the high-mass branch, the paired covariances
I plus or minus (Q-diag Q)/2 give

```
E_mu v_x >= kappa^2(1/2-pi/12)alpha >= (5/42)alpha.
```

In the residual-Frobenius branch use H=B-diag B and t=2^(-13). Then
||H||_op<=6, 6t<=1/2, and 216t^3<=tb0/2. The exact identity
<K,H>=||B||_F^2 and Schur arcsine remainder give

```
E_mu sign(H_A(x))v_x >= kappa^2 t b0/2 >=2^(-31).
```

The high-mass bound is larger than 2^(-31). The uniform bounded-spectrum
Gaussian-sign comparison says the paired physical response is
kappa*f(v_x)+O(n^(-1/6)sqrt(log(en))), where
f(v)=(sqrt(1+v)+sqrt(1-v))/2<=1-v^2/8. Since kappa>=1/2, either
branch therefore discounts the average response by at least 2^(-66)
before that uniform error. Eventually the error is at most 2^(-67),
leaving the gap in (1).

The complete exact arithmetic is replayed by

```
.venv/bin/python computations/paper_discrepancy_2026_09_17_low_cap_response_constant.py
```

It passes using rational arithmetic only. It does not use floating
agreement as a numerical certificate or assert a practical size cutoff.

## 5. Minimax, full support, and the precise subGaussian proxy

For fixed n, let L be the convex hull of all equal paired Gaussian-sign
laws whose two covariance matrices have unit diagonal, sum 2I, and
spectra in [1/2,3/2]. This is a compact convex subset of the finite
physical sign probability simplex. The uniform comparison error above
is uniform over this entire class. For EVERY code-supported mu, Sections
3--4 provide a member of L with mu-average response at most
kappa-2^(-67). Finite-dimensional minimax therefore provides ONE member
of L meeting that bound for EVERY code word. No dependence on the
query is left after minimax.

Each Gaussian sign law is centered. For a pair R_+,R_-, the arcsine
identity cancels every off-diagonal correlation exactly, because
R_-ij=-R_+ij for i!=j. Hence every paired law, and every mixture,
has E hh^T=I exactly.

For Gaussian covariance R<=L I with unit diagonal, the Gaussian
product Holder inequality gives

```
E exp(sum_i t_i sign G_i)
 <= product_i cosh(L t_i)^(1/L)
 <= exp((L/2)sum_i t_i^2).
```

This is the covariance comparison in the primary
[Chen--Dafnis--Paouris theorem](https://arxiv.org/abs/1306.2410),
applied to coordinate exponential functions. We use L=3/2. Convex
mixing preserves that exact MGF upper bound; there is no extra mixture
entropy cost. The proxy is (3/2)I, NOT I. Exact isotropy alone must not
be substituted for a sharp subGaussian proxy.

Finally every Gaussian covariance here is positive definite. Every
orthant therefore has positive probability. All individual physical
sign words have positive probability under every pair and every
mixture. The law in (1) has full physical support; no sign patterns
were discarded to obtain isotropy or the response bound.

## 6. Value and original-problem scope

The source closure holds more generally when
C<c_half+c0/(2K_G), with a gap depending on the strict margin.
At the reported upper C=.493608094 and certified c_half=.4333221116640807,
the safe Krivine value gives the threshold c0>.214885040817.
The numerical corollary (1) deliberately uses simpler rational constants.

This removes the covariance/operator assumption from a genuine uniform
high-energy response statement on low-cap signings. It does not alone
control the maximum fluctuation of many independent bridge columns over
an exponentially large code. Nor is its certified response close to
the 3c/2 parent slope needed by the original recurrence. Both the global
extension payment and the quantitative response target remain separate
obligations. No convergence or endpoint improvement follows from (1).
