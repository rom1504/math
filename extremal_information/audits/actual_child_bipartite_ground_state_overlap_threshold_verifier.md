# Independent verification of the bipartite overlap threshold

Scope: adversarial verification of
[`actual_child_bipartite_ground_state_overlap_threshold.md`](actual_child_bipartite_ground_state_overlap_threshold.md),
with particular attention to the soft maximum, Lindeberg replacement,
Sudakov--Fernique comparison, bounded-difference constant, and asymptotic
threshold.

**Disposition: PASS.**  The theorem and all displayed constants are valid.
No correction to the source theorem is required.  The global-sign duplicate
in the `(x,y)` parametrization is harmless; the check below records both
equivalent ways to handle it.

## 1. Soft maximum and third derivative

The source defines `mathcal Q_(m,n)` as the set of **distinct** matrices
`xy^T`, so its cardinality is exactly `2^(m+n-1)`.  Therefore

```math
 X(z)\le {1\over\eta}\log\sum_{Q\in\mathcal Q_{m,n}}
 e^{\eta\langle z,Q\rangle}
 \le X(z)+{(m+n-1)\log2\over\eta}.                \tag{VBG.1}
```

There is only one soft-max entropy charge, incurred when the smoothed
Gaussian maximum is replaced by its hard maximum in (BG.10).

For one coordinate, differentiation gives

```math
 \partial_e^3F_\eta
 =\eta^2E_\eta[(Q_e-E_\eta Q_e)^3].                \tag{VBG.2}
```

If `u=E_eta Q_e`, direct expansion for a sign variable gives

```math
 E_\eta[(Q_e-u)^3]=-2u(1-u^2),
 \qquad
 \sup_{|u|\le1}2|u|(1-u^2)={4\over3\sqrt3}<1.    \tag{VBG.3}
```

Thus the source's simpler uniform bound
`|partial_e^3 F_eta|<=eta^2` is correct.

## 2. Lindeberg constant and error scale

Condition on every coordinate except the one being replaced.  Taylor's
formula at zero through order two and matching of the first two moments of
a fair sign and a standard Gaussian give

```math
 \left|E f(B_e)-E f(G_e)\right|
 \le {\|f'''\|_\infty\over6}
       \{E|B_e|^3+E|G_e|^3\}.                     \tag{VBG.4}
```

Since

```math
 E|B_e|^3=1,
 \qquad
 E|G_e|^3=2\sqrt{2/\pi},                           \tag{VBG.5}
```

the per-coordinate constant is exactly

```math
 c_{\rm Lin}={1+2\sqrt{2/\pi}\over6}.             \tag{VBG.6}
```

Summing `mn` replacements gives `c_Lin mn eta^2`.  With
`eta=N^(-1/3)`, the two finite errors are

```math
 (N-1)(\log2)N^{1/3},
 \qquad
 c_{\rm Lin}mnN^{-2/3}.                            \tag{VBG.7}
```

Both are `O(N^(4/3))=o(N^(3/2))` on comparable splits.  This verifies
(BG.3), (BG.9), and (BG.10).

## 3. Gaussian increment comparison

For

```math
 Z_{x,y}=x^{\mathsf T}Gy,
 \qquad
 Y_{x,y}=\sqrt n\,g^{\mathsf T}x+\sqrt m\,h^{\mathsf T}y,
```

put `u=x^T x'` and `v=y^T y'`.  Independent direct calculation gives

```math
 \begin{aligned}
 E(Z_{x,y}-Z_{x',y'})^2&=2(mn-uv),\\
 E(Y_{x,y}-Y_{x',y'})^2&=2(2mn-nu-mv).
 \end{aligned}                                    \tag{VBG.8}
```

The second increment exceeds the first by

```math
 2(m-u)(n-v)\ge0.                                  \tag{VBG.9}
```

Sudakov--Fernique therefore has the direction used in the source:
`E max Z<=E max Y`.  One may index this comparison by all `2^(m+n)` pairs
`(x,y)`; the two representatives `(x,y),(-x,-y)` give the same `Z`, so
duplicating them does not change `max Z`.  Alternatively, choose one
representative per rank-one matrix and upper-bound its `Y` maximum by the
maximum over all pairs.  In either convention,

```math
 E\max_{x,y}Y_{x,y}
 =\sqrt{2/\pi}\{m\sqrt n+n\sqrt m\}.              \tag{VBG.10}
```

This verifies (BG.11)--(BG.13) without a hidden factor of `sqrt 2` or a
global-sign entropy error.

## 4. Concentration and entropy transport

Flipping one bridge bit changes every affine score by `+-2`, hence changes
their maximum by at most `2`.  The bounded-differences MGF bound is therefore

```math
 \log E_Ue^{\theta(X-E_UX)}
 \le {\theta^2\over8}\sum_{e=1}^{mn}2^2
 ={mn\theta^2\over2}.                              \tag{VBG.11}
```

Donsker--Varadhan transport and optimization over `theta` give exactly

```math
 E_qX\le E_UX+\sqrt{2mnD(q\Vert U)}.               \tag{VBG.12}
```

For the actual tilt, the conditional half-log-odds of one bridge bit is at
most `|s|t`; hence the chain-rule bound

```math
 D(q_s\Vert U)\le mn\,\kappa(|s|t)                 \tag{VBG.13}
```

turns (VBG.12) into the final term
`mn sqrt(2 kappa(delta t))` in (BG.17).  All factors of two agree.

## 5. Threshold algebra

Let `gamma_N=mn/N^2`.  The normalized geometric term is

```math
 {\sqrt{2/\pi}\{m\sqrt n+n\sqrt m\}\over\rho mn}
 ={\sqrt{2/\pi}\over\rho}\left({1\over\sqrt m}
                                  +{1\over\sqrt n}\right).
                                                               \tag{VBG.14}
```

Since `rho~beta/sqrt N`, its limit coefficient is

```math
 {1\over\beta}\sqrt{2\over\pi}
 \left\{{1\over\gamma_N}+{2\over\sqrt{\gamma_N}}\right\}^{1/2}.
                                                               \tag{VBG.15}
```

The expression in braces decreases with `gamma_N`, so
`gamma_N>=gamma_0` gives the source's

```math
 \beta_{\rm BG}(\gamma_0)
 =\sqrt{2\over\pi}
   \left\{{1\over\gamma_0}+{2\over\sqrt{\gamma_0}}\right\}^{1/2}.
                                                               \tag{VBG.16}
```

Also,

```math
 {R_N\over\rho mn}=O_{\beta,\gamma_0}(N^{-1/6}),
 \qquad
 {\sqrt{2\kappa(\delta t)}\over\rho}\longrightarrow\delta,
                                                               \tag{VBG.17}
```

because `kappa(a)=a^2/2+O(a^4)`.  Finally
`A_rho->1` and `C_(rho,delta)->1+delta`.  Thus the asymptotic lower bound is

```math
 {1-\beta_{\rm BG}(\gamma_0)/\beta-\delta\over1+\delta},
                                                               \tag{VBG.18}
```

positive under the stated choice
`delta<1-beta_BG(gamma_0)/beta`.  At `gamma_0=1/4`,

```math
 \beta_{\rm BG}(1/4)
 =\sqrt{2/\pi}\sqrt8={4\over\sqrt\pi}.            \tag{VBG.19}
```

The integration factor `delta/lambda` in (BG.27) is consequently correct.

## 6. Final disposition

The proof establishes a uniform obstruction for every positive-support
actual child prior because it uses only containment of that prior's support
in the full rank-one orbit.  It does not assert a sharp bipartite spin-glass
constant, a directional reverse-product theorem, or a recurrence.  Those
scope limitations in the source are accurate.

No mathematical or normalization correction was found.
