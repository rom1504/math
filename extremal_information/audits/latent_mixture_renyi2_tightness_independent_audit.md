# Independent audit: latent-mixture Renyi-2 tightness

**Frozen source:**
`extremal_information/drafts/latent_mixture_renyi2_tightness.md`

**SHA-256:**
`9ed3128f305284201fc3f158b684ea16f8ab5edc155ce5f4fd6acdbcfdb85df6`

**Verdict:** **PASS.**  The master tail envelope, compactness and
moment/Orlicz corollaries, diagonal-window construction, optimized
truncation bound, and logical no-prescribed-rate claim are correct.  No
source repair is required.

## 1. Range and master truncation

On a row cube of size `N=2^r`, any law has density `g=N nu(x)` and

```math
1\le K_2(\nu)=N\sum_x\nu(x)^2\le N.
```

Thus `0<=Y_r<=r log 2` as stated.  On `{Y_r<=t}`, the component has
`K_2<=e^t` and its expected shortfall is at most `epsilon_r(e^t)`; on the
complement the pointwise bound is `h_beta`.  Integrating the exact latent
mixture identity gives RT.4 with no missing latent-mass factor:

```math
E\mathcal S_r
\le(1-\alpha)\epsilon_r(e^t)+h_\beta\alpha
\le\epsilon_r(e^t)+h_\beta\alpha.
```

This holds for fixed or order-dependent `t`; only the first term becomes
unknown when `t` grows.

## 2. Tightness and Orlicz criteria

For every fixed `T`, fixed-`K` convergence at `K=e^T` gives

```math
\limsup_rE\mathcal S_r
\le h_\beta\limsup_r\Pr(Y_r>T).
```

Sending `T` to infinity proves RT.1.  Each criterion in RT.2 implies the
same tightness by Markov applied respectively to `Y`, `Y^p`, or
`Psi(Y)`.  The assumptions that `Psi` is nonnegative, increasing, and
diverges are exactly sufficient.

The information comparison RT.7 is also correct:

```math
D(q_r\|\bar\nu_r^{\otimes r})
\le rI(Z;R_1)
\le r\int D(\nu_z\|U_r)d\pi(z)
\le rE_\pi Y_r.
```

The middle inequality is the KL chain identity with the nonnegative term
`D(bar nu||U)` removed, and the last is `D<=D_2`.  Bounded mean `Y` therefore
allows, rather than forces, linear total correlation; the two-half example
shows that order is attainable.

## 3. Diagonal sequence

For each integer `j`, pointwise convergence at the fixed number `e^j`
permits an order `R_j` after which

```math
\epsilon_r(e^j)\le1/j.
```

The `R_j` can be enlarged to be strictly increasing.  Hence

```math
b_r=\max\{j:R_j\le r\}
```

is finite, nondecreasing, and tends to infinity.  At each order,
`r>=R_(b_r)`, so

```math
\epsilon_r(e^{b_r})\le1/b_r\longrightarrow0.
```

No monotonicity interchange or uniform rate is hidden here.  Applying RT.4
proves RT.9.  Markov gives

```math
\Pr(Y_r>b_r)\le {EY_r\over b_r}=o(1)
```

under RT.10.

Although `b_r` was not explicitly constrained by `r log2`, this creates no
problem: RT.4 is valid for all `t>=0`; thresholds above the maximum possible
`Y_r` simply have zero tail.  One could also enlarge `R_j` so that
`j<=R_j log2` if a within-range window were desired.

## 4. Tail-envelope optimization

If `Pr(Y_r>t)<=T_r(t)`, RT.4 bounds the same left side for every
`t in [0,r log2]`; taking the infimum therefore proves RT.12.  Restricting
to this interval loses nothing: beyond `r log2` the tail is already zero
while `epsilon_r(e^t)` is nondecreasing.

Thus a chosen `t_r` yields a quantitative theorem exactly when both its
component-error and latent-tail terms vanish.  The source correctly refrains
from inferring such a choice from qualitative fixed-`K` convergence alone.

## 5. No universal prescribed rate

For any prescribed divergent `a_r`, the abstract profile

```math
\widetilde\epsilon_r(K)
=h_\beta\mathbf1_{\{\log K\ge a_r\}}
```

is nonnegative, bounded by `h_beta`, and nondecreasing in `K`.  For every
fixed `K`, it is eventually zero, but at `K=e^(a_r)` it is always
`h_beta`.  Truncating `a_r` at `r log2` preserves divergence and keeps the
test point inside the realizable density range.

This proves exactly the stated logical limitation: fixed-`K` convergence
alone cannot imply any externally named growing regime.  The profile is not
claimed to arise from conference pressure, so the argument does not
overstate the obstruction.

## 6. Scope

The theorem is a qualitative compactness result for the random component
complexity `Y_r=D_2(nu_(Z,r)||U_r)`.  It neither asserts necessity of
tightness nor a quantitative component theorem.  Those limitations are
explicit and consistent with every displayed conclusion.
