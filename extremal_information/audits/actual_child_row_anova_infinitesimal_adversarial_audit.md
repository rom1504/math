# Adversarial audit: infinitesimal actual-child row ANOVA

**Object audited:**
[`drafts/actual_child_row_anova_infinitesimal.md`](../drafts/actual_child_row_anova_infinitesimal.md)

**Verdict:** **PASS after restricting the global-minimizer assertion to
`lambda>0`.**  The Taylor coefficients in RA.5--RA.8 and every factor in the
actual-sector overlap formula RA.17--RA.20 are correct.  All remainders are
finite-system remainders; there is no justified uniformity in `N`, and the
draft explicitly says so.

## 1. Product variational expansion

Write

```math
{dp_i\over dU_i}=1+\lambda h_i+\lambda^2k_i+O(\lambda^3),
\qquad Eh_i=Ek_i=0.
```

For the objective multiplied by `lambda`, direct expansion gives

```math
D(p||U)+\lambda E_pL
=\lambda\mu
+\lambda^2\sum_i\left({1\over2}\|h_i\|_2^2
                         +\langle L_i,h_i\rangle\right)
+\lambda^3 R+O(\lambda^4).
```

The quadratic minimizer is `h_i=-L_i`.  At that value, all `k_i` terms
cancel, and

```math
R={1\over6}\sum_iE L_i^3
  +\sum_{i<j}\langle L_{ij},L_iL_j\rangle.
```

This verifies both the sign and the two coefficients in RA.6 and RA.13.
Combining with the ordinary cumulant expansion

```math
V_\lambda=\mu-{\lambda\over2}\operatorname{Var}L
              +{\lambda^2\over6}\kappa_3(L)+O(\lambda^3)
```

gives all three lines of RA.5, including the extra factor `lambda` in
`I_lambda^leftarrow=lambda(V_row-V)`.

## 2. Global uniqueness and the one genuine repair

For `lambda>0`, comparison with the uniform product gives

```math
D(p_\lambda^*||U)
\le\lambda(\mu-E_{p_\lambda^*}L)=O_L(\lambda).
```

Hence every global minimizer approaches the interior point `U`.  Entropy has
positive-definite Hessian on the mean-zero factor tangent space there, and
the analytic implicit-function theorem supplies a unique nearby critical
branch.  This proves global uniqueness and RA.7 for all sufficiently small
positive `lambda`.

The original wording “small `|lambda|`” was false for the stated minimization
problem: when `lambda<0`, the coefficient `1/lambda` of entropy changes sign,
and a global minimizer need not approach `U`.  I patched the theorem and
proof to say `lambda downarrow 0` and `lambda>0`.  The analytic critical
branch itself may extend through zero, but global variational minimality is
only being claimed on the positive side.

## 3. Forward total-correlation tangent

The escort score at zero is `-(L-mu)`.  Its `i`th row marginal has score
`-L_i`, so the product of its row marginals has score `-L_add`.  For two
finite-space laws with densities `1+lambda a+O(lambda^2)` and
`1+lambda b+O(lambda^2)`,

```math
D(P_\lambda||Q_\lambda)
={\lambda^2\over2}\|a-b\|_2^2+O(\lambda^3).
```

Orthogonality of row ANOVA therefore gives exactly
`TC(q_lambda)=lambda^2 sigma_cross^2/2+O_L(lambda^3)`.  RA.8 has no missing
factor, and it does not imply a fixed-`lambda` comparison of the two KL
directions.

## 4. Actual-sector overlap mapping

At zero bridge, the fixed-orientation augmented child law is proportional to

```math
e^{t\tau(H_A(X)+\epsilon H_D(Y))}.
```

The involution `X->-X` makes `Q=tau XY^T` central.  Thus

```math
L_u(B)={u^2\over2}E\langle Q,B\rangle^2+O(u^4).
```

After removing diagonal constants, the coefficient of the Walsh monomial
`B_ij B_kl` for distinct edges is

```math
u^2 E[X_iX_kY_jY_l].
```

Different rows (`i<k`, with `j,l` arbitrary) are exactly the two-row ANOVA
monomials; same row and `j<l` are exactly the singleton-row monomials.
Walsh orthonormality then proves RA.17 and RA.18 with coefficient one.
Taking `j=l` gives `n` identical copies of `E X_iX_k`, proving RA.19.

Conditional on `tau=s`, the two child laws factor as
`e^{tsH_A(X)}` and `e^{t epsilon sH_D(Y)}`.  Averaging their correlation
matrices with the true sector weights `pi_s` gives RA.20.  There is no
orientation factor or stray `tau`, since `tau^2=1` in the two-edge product.

## 5. Uniformity and research scope

The proof is exact at each fixed finite system.  Its `O_L(lambda^3)` and
`O(u^4)` constants can grow arbitrarily with `m,n`; bounded coordinate
oscillation alone does not control the cumulants needed at fixed `lambda`,
and `u=t=beta/sqrt(N)` multiplies `Theta(N^2)` coordinates.  Therefore no
interchange of `lambda->0`, `u->0`, and `N->infinity` is currently legal.

RA.17 is a genuine actual-child, low-order nonradial statistic, but it does
not yet narrow the fixed-tilt SML.  Such a reset would require a uniform
cluster/cumulant theorem or a direct physical-scale bound on
`sigma_cross^2(L_t)`.  The draft makes this limitation explicit.

## 6. Audit of the physical mixed-response identity (RA.3)

**Verdict: PASS.**  Let `R_i` denote conditional expectation over bridge row
`i` and write the independent-replacement difference as `D_i=I-R_i` on the
enlarged probability space containing the original and replacement row.
For a row-ANOVA component `L_S`,

```math
\mathbb E(D_iD_kL_S)^2
=4\|L_S\|_2^2\mathbf 1_{\{i,k\}\subseteq S}.
```

Each replacement difference contributes a factor two to the squared norm,
so the factor `1/4` in (RA.22a) is correct.  Summing over `i<k` gives
`sum_S binom(|S|,2)||L_S||_2^2`; hence both constants in (RA.23a) are exact.

For the interpolation
`B_(s,v)=B^(ik)+s(R_i-R_i')+v(R_k-R_k')`, differentiation of the parent
log-partition gives

```math
\partial_sL=t\,\mathbb E_{\nu_{s,v}}Z_i,
\qquad
\partial_v\partial_sL
=t^2\operatorname{Cov}_{\nu_{s,v}}(Z_i,Z_k).
```

The orientation variable is already included in `Z_i,Z_k`; there is no
additional sign or factor two.  The two-dimensional fundamental theorem of
calculus has the same corner signs as `square_(ik)L`, so (RA.24a) is also
correct.  The draft correctly warns that the covariance is a parent
interpolation quantity and that `J_2` overcounts ANOVA order `s` by
`binom(s,2)`.

## 7. Audit of the rare-well ceiling (RA.4)

**Verdict: PASS, with the stated asymptotic interpretation.**  A bit flip
changes `S` by two and therefore changes (RA.27a) by at most
`2 gamma/sqrt(N)`.  Since `|F_N|<=eta N` and

```math
\Pr\{|S|>aN^{3/2}\}\le 2e^{-a^2N/2},
```

one has `Var(F_N)<=eta^2N^2 2e^{-a^2N/2}=e^{-Omega(N)}` (after absorbing the
polynomial prefactor).  For fixed `c_0`, the binomial moderate-deviation
estimate at threshold `c_0N^(3/2)` is

```math
\log\Pr\{|S|\ge c_0N^{3/2}\}
=-(c_0^2/2+o(1))N;
```

this is within the moderate-deviation regime because the threshold divided
by the number `N^2` of bits tends to zero.  Therefore

```math
\mathbb EF_N+{1\over\lambda}\log\mathbb Ee^{-\lambda F_N}
\ge\left(\eta-{c_0^2\over2\lambda}-o(1)\right)N,
```

and the condition `lambda eta>c_0^2/2` is exactly the positive-gain
condition.  Such fixed positive constants exist (take `a` small and
`gamma` large relative to `eta`).  This is a method-specific ceiling only:
the construction is not asserted to be a child-induced log-partition, as
the draft explicitly records.
