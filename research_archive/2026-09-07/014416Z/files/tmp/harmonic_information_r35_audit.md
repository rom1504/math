# Independent audit of `harmonic_information_r35.md`

## Verdict

**Accepted, with two wording/numerical-scope corrections and one uniformity
qualification.**  I found no sign, factor-of-two, chart-count, or KL-chain
rule error in (R35.1)--(R35.19).  The checker passes, and I independently
rederived the endpoint, lifted-selector, orientation, migration, and
backtracking identities below.

This audit does not upgrade any finite decimal to an asymptotic theorem.

## 1. Endpoint erasure and selector cancellation

For any coordinate `j`, the relative-entropy chain rule gives

```math
D(\mu\Vert\nu)
=D(\mu_{-j}\Vert\nu_{-j})
+\mathbb E_{\mu_{-j}}D(\mu_j(\cdot\mid D_{-j})
                         \Vert\nu_j(\cdot\mid D_{-j})).
```

Thus (R35.2) and, after summing the `n` chart coordinates, (R35.3) are exact.
The context weight is the endpoint marginal `mu_(1,-j)`, equivalently the
edge mass `M_e(1)`, not the base context mass.

For `P(S,d)=q(S)mu_S(d)`, subtracting the marginal chain rule from the full
chain rule in either order gives

```math
J_j=C_j+I(S;D_j\mid D_{-j}).
```

So the sign in (R35.5) is correct: `C_j=J_j-I_cond`.  Also

```math
J_j=\sum_Sq(S)
[D(\mu_S\Vert\nu)-D((\mu_S)_{-j}\Vert\nu_{-j})].
```

The mixture identity used in (R35.7),
`sum_S q(S)D(mu_S||nu)=D(mu||nu)+I(S;D)`, is exact.

## 2. Chart coordinates and the inclusion count

The chart has `n`, not `n+1`, coordinates:

- one global-orientation bit, indexed `0` in the checker; and
- projective flips of original vertices `1,...,n-1`; original vertex `0` is
  the gauge-fixed vertex and has no projective-flip coordinate.

For a fixed selector `S`, `h_S` depends only on the oriented restriction to
`S`.  Therefore a projective bit `i>=1` has `J_(S,i)=0` when original vertex
`i` is omitted.  The possibly active coordinates number

```math
1+|S\cap\{1,\ldots,n-1\}|=m+1-\mathbf1_{\{0\in S\}}.
```

Thus (R35.6)--(R35.7) have the correct count.  The use of `0` for both the
orientation coordinate and the gauge-fixed original vertex is potentially
confusing but not mathematically wrong; a final ledger version should call
the latter `v_*` or explicitly repeat this distinction.

Conditioning further on `Z=1_{i in S}` verifies (R35.9): in the `Z=0`
mixture, every likelihood is constant across the `i` edge, so its bit
conditional is exactly the base conditional.  No independence assumption is
being made.

## 3. Orientation quotient and exact certificate

Applying the same chain rule to the orientation coordinate gives (R35.10),
including the sign `D(full)=D(quotient)+C_0`.

I reran the exact `Fraction` calculation.  At
`A_6,m=3,beta=(log 2)/2`, the states
`d=-xx^T`, `x=(1,-1,-1,-1,-1,-1)`, and `-d` have energies `10,-10` and

```math
U(d)=847888/3590575,\qquad
U(-d)=111448/963325,
```

with ratio `1165846/571171`.  Since `f` is a common normalization times
`U`, this is exactly the endpoint/base conditional-odds multiplier on that
orientation edge.  Both base masses are positive, so its conditional KL and
heat-bath energy are strictly positive.  The certificate proves that the
orientation term is not **algebraically** zero.

Wording correction: the status sentence “its edge cannot be deleted” should
read “it cannot be deleted without retaining or separately bounding `C_0`.”
The finite certificate does not rule out an asymptotic estimate making `C_0`
negligible on a particular family.  The body of the memo already uses the
proper scoped formulation.

## 4. Migration derivatives and signs

For a context edge `e`,

```math
M_e(t)=M_e(0)e^{A_e(t)-\psi(t)},\qquad
\ell_e(t)=A_e(t)-\psi(t),
```

so `M'_e=M_e(A'_e-psi')`.  Differentiating
`D^-_j=E_M ell` makes the direct derivative of `ell` average to zero and
leaves

```math
(D^-_j)'=\operatorname{Cov}_M(A'_e,A_e).
```

Since `k_e=tA'_e-A_e`,

```math
\operatorname{Cov}_M(A'_e,k_e)
=t\operatorname{Var}_M(A'_e)-(D^-_j)'.
```

Combining this with

```math
(D(\mu_t\Vert\nu))'=t\operatorname{Var}_{\mu_t}(g),
\quad
\operatorname{Var}_{\mu_t}(g)
=\mathcal E_{t,j}(g)+\operatorname{Var}_M(A'_e)
```

gives exactly (R35.14), with a plus covariance and one factor `t` on the
energy.  Thus (R35.12)--(R35.15) have the correct signs.

For backtracking,

```math
\operatorname{Cov}_M(A'_e,k_e)=\sum_eM'_e k_e,
\qquad 0\le k_e(t)\le K_e.
```

Dropping the favorable terms with `M'_e>0` gives (R35.18), and summing gives
(R35.19).  There is no missing edge mass or factor `t`.

Uniformity qualification: the constant called `K` in (R35.17) and (R35.19)
must be an `n`-independent constant uniform over the relevant exact
minimizers, ratio window, and prescribed temperature.  Otherwise those
statements do not preserve the exponent in (10.983).  Renaming it
`K_bt` would also avoid confusion with the edge costs `K_e`.

## 5. Numerical-scope correction

The memo says that sampled `A_4/A_6/A_8/A_9` paths had nondecreasing marginal
KL on the tested grids.  As currently written, the checker's migration loop
samples only `A_6,A_8,A_9`; `A_4` appears only in the endpoint-information
loop.  Either remove `A_4` from that sentence or add it to the migration
loop.  I independently ran `migration_audit(A4,0.5,3,order=96)` and obtained
a minimum sampled marginal-KL derivative about
`2.6518e-6 > 0`, so adding that case supports the stated numerical claim.

More precisely, the checker samples the analytic derivative at
Gauss--Legendre nodes; it does not prove monotonicity on the whole interval.
“All sampled marginal-KL derivatives were nonnegative” is the exact claim.

## 6. Asymptotic scope

- (R35.16) proves only the endpoint-cost and adverse-migration parts of
  (10.983).  Uniform restoring bounds and adjacent-selector Hellinger control
  remain separate, as the memo correctly states.
- (R35.17) or (R35.19), with a uniform constant, makes an independently
  proved endpoint erasure bound sufficient for migration; neither bounds the
  endpoint erasure sum itself.
- Monotonicity alone yields
  `sum_j D^-_j(1)=nD(mu||nu)-sum_j C_j`, so the memo is correct that it does
  not remove the linear loss without an additional endpoint theorem.
- The finite values in (R35.8) and (R35.20) are diagnostics only.  In
  particular, `sum weighted energy > sum C_j` for the displayed `A_9` case
  falsifies coefficient-one domination but not domination by a uniform
  larger constant.

Subject to the three scoped edits above, the memo is ready to import into the
ledger as verified mathematics plus explicitly numerical finite evidence.
