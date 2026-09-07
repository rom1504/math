# Quenched orientation synchronization and adaptive cavity expansion

Status: Sections 1--6 independently audited PASS, with the zero-weight
boundary correction incorporated; see
`decisive_audit_quenched_orientation_cavity_2026_09_07.md`. This is
an exact-optimizer interpolation lemma, not an original thermodynamic-limit
theorem. Two remaining comparison obligations are stated in Section 6.

## 1. Regularization before optimizing the signing

Let real nonnegative edge magnitudes lambda_e be prescribed, let A_e be sign
choices, and put

```math
Z_A(g)=sum_(s=+-1,x in {+-1}^n)
 exp[s sum_e lambda_e A_e x_i x_j + h g s],
Phi(A)=E_g log Z_A(g),                  g~N(0,1).
```

The signing is chosen BEFORE the independent g; optimize Phi, not the
pointwise random objective. Let h=h_n with h_n->infinity and h_n=o(n).
For every A, Jensen and the Lipschitz bound in the added field give

```math
log Z_A(0) <= Phi(A) <= log Z_A(0)+h sqrt(2/pi).       (1)
```

Thus the minimum normalized pressures differ by at most
sqrt(2/pi) h_n/n, uniformly in every magnitude profile. At uniform
lambda_e=beta/sqrt n, the unregularized minimum pressure differs from
beta M_n/n^(3/2) by a quantity between 0 and (1+1/n)log2. Consequently
limits at each fixed beta, followed by beta->infinity, would still settle
the original minimum-value convergence problem.

## 2. Uniform synchronization of the global energy orientation

Fix ANY deterministic signing and magnitudes, including any cavity system.
Write Z_+,Z_- for its branch partition functions before adding hg. The
conditional branch field is

```math
v=(log Z_+-log Z_-)/2+h g.
```

For two independent replicas under the SAME g, the probability of opposite
orientations is (1/2)sech^2(v). Since the standard normal density is at most
1/sqrt(2pi),

```math
E_g <1_{s_1 != s_2}> <= 1/[h sqrt(2pi)].             (2)
```

This is uniform in the branch free-energy gap, including gaps of order n.
For every replica observable |F|<=1,

```math
|E_g <s_1 s_2 F>-E_g <F>| <= 2/[h sqrt(2pi)].        (3)
```

In particular signed overlap squares have the usual sign up to a uniformly
vanishing error. This does NOT make the signed one-replica correlation matrix
E[s xx^T] positive semidefinite.

## 3. Expected edge optimality and the retained quadratic term

Choose a GLOBAL minimizer A of Phi. Delete a fixed edge e, preserving all
other signs and the field hg. Let r_e(g)=<s x_i x_j>_0 be the cavity
correlation, let t_e=tanh(lambda_e), and abbreviate t=t_e, y(g)=A_e r_e(g).
The exact insertion identity is

```math
log Z_A(g)-log Z_0(g)=log cosh(lambda_e)+log(1+t y(g)). (4)
```

Comparing with the opposite edge sign yields

```math
E_g artanh(t y)<=0.                                  (5)
```

For 0<t<1, the series of artanh gives

```math
E_g y <= kappa_t:=t^2/[3(1-t^2)].                    (6)
```

Write m_e=E_g r_e(g), v_e=E_g r_e(g)^2. The nonnegative sign-selection
defect d_e=A_e m_e+|m_e| satisfies

```math
0<=d_e<=2 kappa_t.                                   (7)
```

For 0<t<=1/2, Taylor's theorem in (4), with
|log(1+z)-z+z^2/2|<=|z|^3/[3(1-|z|)], gives

```math
E_g log(1+t y)
=-t |m_e|-(t^2/2)v_e+epsilon_e,
|epsilon_e|<=2t^3.                                  (8)
```

The leading absolute value is |E_g r_e|, NOT E_g|r_e|. The quadratic term
need not be negligible: averaging the optimality comparison in g does not
justify the pointwise sign condition A_e r_e(g)<=0. At t=O(n^(-1/2)), summing
the remainder over all edges costs only O(sqrt n).

At t=0 the edge sign is unused, so (5) gives NO sign-selection information:
(6), (7), and the signed-correlation expansion below are only asserted for
positive edge magnitudes. Formula (8) itself is trivially zero at t=0. One
may choose an unused sign opposite m_e, but this choice is not forced by
global optimality. This boundary distinction was identified in the independent
`transfer_seeds` audit.

## 4. Differential form for globally minimized magnitude profiles

Let R_e(g)=<s x_i x_j> be the FULL correlation. Exactly,

```math
A_e R_e(g)=(t+y)/(1+t y)
=t+y-t y^2+t^2 (y^3-y)/(1+t y).                       (9)
```

Equations (7) and (9) imply, for 0<t<=1/2,

```math
A_e E_g R_e(g)=t-|m_e|-t v_e+rho_e,
|rho_e|<=3t^2.                                      (10)
```

For completeness, the remainder bound uses
|y^3-y|<=1, 1/(1-t)<=2, and
2kappa_t<=8t^2/9. No assumption on the number of Gibbs states is used.

Take a C^1 path lambda_e(u), leaving h fixed, and let
f_n(u)=n^(-1)min_A Phi_u(A). This finite minimum is locally Lipschitz.
At every differentiability point it has the gradient of an active global
minimizer (and the same derivative if several active branches meet there).
At an interior zero of a nonnegative C^1 magnitude path, lambda'_e=0, so that
edge contributes nothing. At path endpoints use a one-sided positive-magnitude
limit or omit zero-weight edges; do not apply (10) to an arbitrary unused sign.
Thus

```math
f_n'(u)
= (1/n)sum_e lambda'_e t_e
 -(1/n)sum_e lambda'_e |m_e|
 -(1/n)sum_e lambda'_e t_e v_e + R_n(u),
|R_n(u)| <= (3/n)sum_e |lambda'_e| t_e^2.             (11)
```

If lambda_e=beta w_e(u)/sqrt n, |w_e|<=W, and
sum_e|w'_e(u)|<=binom(n,2)V(u), with V integrable, then the integrated
normalized error in (11) is O(beta^3 W^2 n^(-1/2) integral V).
This also permits integrable square-root endpoint singularities by truncation.

The cavity and full measures differ by a tilt exp(lambda_e A_e s x_i x_j).
For t_e<=1/2 their bounded-observable expectations differ by O(t_e), uniformly
in g. Consequently v_e may be replaced in (11) by E_g R_e(g)^2 at the same
O(n^(-1/2)) integrated scale for these paths. These FULL squares are genuine
replica observables in a common Gibbs measure, so (3) applies directly.

## 5. A useful simultaneous first- and second-moment budget

At equal lambda_e=lambda>0, Phi_A(lambda) is convex and even in lambda:
replace s by -s and g by -g for evenness. Hence its radial derivative is
nonnegative. This gives sum_e A_e E_g R_e>=0, with E=binom(n,2).
Using the exact identity

```math
A_e R_e=t+(1-t^2)y/(1+t y)
```

and (7), one obtains

```math
sum_e |m_e| + t sum_e E_g[r_e^2/(1+t A_e r_e)]
 <= E t/(1-t^2)+2E kappa_t.                         (12)
```

In particular

```math
sum_e |m_e| + [t/(1+t)]sum_e v_e
 <= E t/(1-t^2)+2E kappa_t.                         (13)
```

Thus sum|m_e|=O_beta(n^(3/2)) at uniform beta/sqrt n. The formula retains a
potentially leading second-moment cost instead of losing it in a triangle
bound. For heterogeneous magnitudes, simultaneous radial scaling gives the
weighted version

```math
sum_e lambda_e(1-t_e^2)|m_e|
+sum_e lambda_e t_e(1-t_e^2) E_g[r_e^2/(1+t_e A_e r_e)]
 <= sum_e lambda_e t_e +(2/3)sum_e lambda_e t_e^2.     (14)
```

All statements concern actual globally minimized finite systems, not a
postulated limiting Gibbs law.
In (14), zero-weight edges contribute zero and are omitted from the
sign-selection argument.

### Cutwise strengthening

The same budgets hold with EVERY edge sum restricted to an arbitrary vertex
cut. For fixed g, flip all x_i on one side of the cut. This leaves entropy
and hg s unchanged, and the relative entropy between the Gibbs law and its
flipped law is exactly

```math
2 sum_(e in cut) lambda_e A_e R_e(g) >=0.
```

Thus the positivity input for (12)--(14) holds cutwise, before averaging g.
Repeating their algebra proves their cutwise versions. At equal positive
couplings, taking a one-vertex cut gives

```math
sum_(j!=i) |m_ij|+[t/(1+t)]sum_(j!=i) v_ij
 <= (n-1)t/(1-t^2)+2(n-1)kappa_t.                   (15)
```

This controls degrees at O_beta(sqrt n), still the critical rather than a
vanishing conditional-spectral scale. This subsection uses the same exact
cut-positivity argument independently reconstructed by the bridge agent for
the unregularized cavity kernel.

## 6. Missing obligations before a convergence claim

The first cavity term in (11) is the directional sum of the absolute-MEAN
kernel |E_g r_e|. Neither its global l1 budget nor generic correlation-matrix
positivity proves the required cross-versus-within inequality. The retained
quadratic term has the usual replica structure after synchronization and
could conceivably compensate bad directions, but that compensation has NOT
been proved.

There is also an independent endpoint issue. At a block-diagonal system the
two-sided partition function is

```math
Z_B,+ Z_C,+ + Z_B,- Z_C,-,
```

not the product (Z_B,++Z_B,-)(Z_C,++Z_C,-). The energy orientation is common
to the children, and choosing opposite child signing orientations may exploit
cancellation of their pressure asymmetries. The random global field removes
replica sign mismatches; it does not establish child factorization or a sharp
optimized-child comparison. Both this endpoint issue and the directional
cavity-kernel inequality remain necessary in an original-limit proof.
