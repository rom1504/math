# Independent audit of the double-partial common-block no-go

Status: **passed**.  This note verifies
[`generic_double_partial_common_block_no_go.md`](generic_double_partial_common_block_no_go.md)
and compares it with the longer
[`generic_double_partial_product_group_row_lifetime_no_go.md`](generic_double_partial_product_group_row_lifetime_no_go.md).
The construction is correctly scoped as generic rather than an actual-child
counterexample.

## 1. Support and one-row law

The factor supports have sizes `2^(ell+1)` and `2^(s+1)`.  The outer-product
map has precisely the simultaneous-sign kernel `(X,Y)~(-X,-Y)`, so its
uniform image has size `2^(ell+s+1)`.  This proves DP.2--DP.3 and the stated
projective factor atom sizes.

For any single row, multiplication by its fair row sign transforms the
right factor into

```math
(a\mathbf1_r,c_1,\ldots,c_s),
```

where all displayed signs are independent and fair.  Averaging its binary
channel therefore leaves the last `s` outputs fair and gives exactly

```math
p_i(b)={\cosh(uS_D(b))\over(\cosh u)^r}.             \tag{VDP.1}
```

Thus the canonical factor `nu tensor U_s`, and hence R and P in DP.7--DP.8,
have the correct normalization.

## 2. Conditional log-cosh identity

Fix the exterior of `C times D` and condition the latent sum on
`a=sigma tau`.  The core contribution is `aT`; summing the strictly positive
exterior weights in the two sectors gives

```math
W_+e^{uT}+W_-e^{-uT}
=2\sqrt{W_+W_-}\cosh\!\left(uT+{1\over2}\log{W_+\over W_-}\right). \tag{VDP.2}
```

All other factors can be absorbed into an exterior-only additive term in
`log p`.  This proves DP.10 exactly, without an approximation or a missing
factor two.

Under R and P, the core is independent of its exterior: every row factor is
itself a product between its first `r` and last `s` coordinates, and rows are
independent.  Their exterior laws are identical because R and P differ only
on the first `r` coordinates of the first `k` rows.  Consequently the
exterior terms in DP.10 cancel in expectation.  Subtracting the common value
at `T=0` and using the one-Lipschitz property of log-cosh proves DP.11.

Chebyshev gives fair probability at least `1/2` to
`|S_D|<=sqrt(2r)`.  Since `r<=N`,

```math
Z_{r,N}\ge{1\over2}\cosh(\sqrt2\,\beta)^{-\lambda}=z_*>0. \tag{VDP.3}
```

The tilted density is at most `1/Z`, so `E_nu S_D^2<=r/z_*`.  Independence
and centering then give

```math
E_R|T|\le\sqrt{kr/z_*},\qquad E_P|T|\le\sqrt{kr},
```

which is exactly DP.13 after multiplication by `u=beta/sqrt(N)`.

## 3. KL comparison and its sign

For every law A,

```math
D(A\Vert q)=D(A\Vert U)+\lambda E_A\log p
             +\log E_Up^{-\lambda}.                 \tag{VDP.4}
```

Since `I^leftarrow<=D(P||q)`, subtraction has the direction

```math
J-I^\leftarrow\ge D(R\Vert q)-D(P\Vert q).          \tag{VDP.5}
```

The common normalizer cancels.  R and P agree on every nuisance row, while
each of their `k` differing rows contributes `D(nu||U_r)=d_N`; hence DP.16
follows with the displayed **plus** sign in front of
`lambda(E_R-E_P)log p`.  DP.13 can reduce this by only `O(sqrt N)`.

The triangular-array limit `uS_D -> N(0,beta^2 rho)` is valid, and both the
tilt and the tilt times `log cosh` are bounded.  The limiting KL `d_0` is
strictly positive because the Gaussian tilt is nonconstant.  Thus
`k d_N=kappa d_0N+o(N)`, proving DP.17 with no hidden linear cancellation.

## 4. Posterior invariance

The outer-product support is a multiplicative subgroup of the bridge cube.
Entrywise multiplication by a support word preserves U, fixes the mixture
likelihood, and acts transitively on the uniform latent components.  Every
bridge tilt whose density is a function of p is invariant under this action.
Changing variables in the averaged Bayes posterior therefore makes every
latent multiplier equal, and normalization proves `bar mu=mu`.  DP.18 and
the vanishing of every deterministic quotient-retuning term are correct.

## 5. Comparison with the longer alternate proof

The alternate file studies the identical double-partial prior.  Its
core-fibre lower bound and conditional-channel area estimate yield an
extensive gap only when

```math
\kappa d_0>
 {\lambda+\lambda^2\over2}\beta^2\delta_{\rm out}.
```

The shorter proof instead chooses P to agree exactly with the canonical
factors on all nuisance rows and uses DP.10 to compare only the common core.
It proves the gap for every fixed positive block fraction and all
`beta,lambda>0`.  It therefore strictly subsumes the alternate theorem's
conclusion for the same construction.

The longer file contains a potentially reusable core-fibre/area technique,
but it adds no surviving mathematical claim here and is substantially less
sharp.  Recommendation: omit it from the theorem ledger and active context;
retain it only as an archived alternative proof if that technique is wanted
for a setting where the exact log-cosh representation is unavailable.

## 6. Verdict

All requested identities, independence statements, moment estimates, KL
signs, asymptotic rates, and posterior-symmetry claims pass.  No correction
to the stronger proof is required.
