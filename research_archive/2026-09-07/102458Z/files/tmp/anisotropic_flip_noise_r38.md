# Wave 38 main continuation: an anisotropic product-noise compression criterion

## Outcome

Uniformly noising a complement-flipped ground cannot reach a power-saving
row cutoff: its required high-cap regime lies beyond the universal spectral
maximum of every structured flip.  Allowing coordinate-dependent means gives
an exact sufficient certificate in terms of expected flipped energy, expected
parent row square, and product relative entropy, but box and spectral
calculations show that this constant-probability certificate is also
impossible at the project entropy scale.  Independent product tilting reaches
the `n^(3/4)` entropy barrier but cannot supply a fixed power saving through
this mean-and-variance mechanism.

The heterogeneous formulas remain useful diagnostics: they identify exactly
what fails and show that merely choosing a direction with smaller parent row
image does not repair the entropy cost.  A successor must leave the
independent product family or exploit a rare event without paying for it by a
constant-probability product tilt.

The expectation and variance identities are verified exactly over rational
arithmetic by `tmp/anisotropic_flip_noise_r38_check.py`.

## 1. Exact heterogeneous product-noise identities

Fix a selector `S`, write `B=B_S=A^(F_S)`, and choose any oriented base state
`(sigma,x)`.  Gauge the flipped signing by

```math
C=\sigma D_xBD_x.
```

Let `z_i` be independent signs with means `rho_i in (-1,1)`, and use the
noisy parent state `x'=x circ z` with the same orientation.  Put

```math
v_i=1-\rho_i^2,
\qquad
\rho=(\rho_i)_{i=1}^n.
```

Then the expected complement-flip energy is

```math
\boxed{
\mathbb E[\sigma(x')^TBx']=\rho^TC\rho.
}
\tag{R38.A1}
```

Writing `z=rho+eta`, the centered linear and quadratic terms are
orthogonal.  Independence and the zero diagonal give the exact variance

```math
\boxed{
\operatorname{Var}[\sigma(x')^TBx']
=4\sum_i v_i(C\rho)_i^2
+4\sum_{i<j}v_iv_j
=4\sum_i v_i(C\rho)_i^2
+2\left\{\left(\sum_i v_i\right)^2-\sum_i v_i^2\right\}.
}
\tag{R38.A2}
```

The parent row square has the separate exact mean

```math
\boxed{
\mathbb E R_2(x')
=\lVert AD_x\rho\rVert_2^2
+(n-1)\sum_i(1-\rho_i^2).
}
\tag{R38.A3}
```

Indeed, the second term is the sum of the independent column variances and
every column of `A` has squared norm `n-1`.  Unlike the uniform formula,
(R38.A3) exposes the direction of `rho` relative to the parent matrix.

The relative entropy of this product law from the uniform spin cube is

```math
\boxed{
D(\rho)=\sum_i d(\rho_i),
\qquad
d(r)=\frac{1+r}{2}\log(1+r)
     +\frac{1-r}{2}\log(1-r).
}
\tag{R38.A4}
```

For example, `d(r)=r^2/2+O(r^4)` and `d(r)<=r^2` on `|r|<=1/2`.

## 2. A per-selector certificate implies fractional compression

Fix a row cutoff `R` and an entropy budget `L`.  Call `S`
**anisotropically regularizable** if there are `(sigma,x)` and
`rho in (-1,1)^n` such that

```math
\boxed{
\begin{aligned}
\rho^TC\rho&\ge4q_n,\\
\operatorname{Var}[\sigma(x')^TB_Sx']&\le q_n^2,\\
\mathbb E R_2(x')&\le R/8,\\
D(\rho)&\le L.
\end{aligned}
}
\tag{R38.A5}
```

Chebyshev makes the probability of flipped energy below `q_n` at most
`1/9`, while Markov makes the probability of row square above `R` at most
`1/8`.  Thus the product law assigns probability at least `1/2` to the set

```math
\{d:R_2(d)\le R,\ S\in I_d\}.
\tag{R38.A6}
```

This probability converts to counting mass without requiring an equal-bias
Hamming sphere.  If `P` is any law on a finite set, `U` is uniform, and
`P(G)>=alpha`, data processing of relative entropy gives

```math
D(P\Vert U)
\ge \alpha\log\frac{\alpha}{U(G)}
 +(1-\alpha)\log\frac{1-\alpha}{1-U(G)}.
```

Consequently

```math
\boxed{
U(G)\ge
\exp\left\{-\frac{D(P\Vert U)+h(\alpha)}{\alpha}\right\}.
}
\tag{R38.A7}
```

At `alpha=1/2`, (R38.A4)--(R38.A7) show that every regularizable selector is
incident to at least

```math
2^n\exp\{-2L-O(1)\}
```

redundant oriented spin states which are row-good.  Summing these incidences
against any feasible dual of (10.1046), with the same harmless two-to-one
projective accounting as in the uniform-noise proof, yields

```math
\boxed{
\sum_{S\text{ anisotropically regularizable}}y_S
\le\exp\{2L+O(1)\}.
}
\tag{R38.A8}
```

Conditionally, (R38.A5) with `L=O(R/q_n)=O(n^(3/4-c))` for every selector in
the remaining dual support would prove the desired fractional compression.
More generally, it would remove any subfamily for which the certificate can
be constructed.  In fact its displayed energy constant is impossible even
before taking limits.  With `u=D_x rho`, box multilinearity gives

```math
|\rho^TC\rho|=|u^TB_Su|
\le |u^TAu|+2|(P_Su)^TA(P_Su)|
\le q_n+2Q(A[S])\le3q_n.
```

Thus no selector satisfies the literal `4q_n` line.  The next paragraph
shows that replacing `4` by any fixed positive constant does not repair the
project entropy scale.

The variance clause is automatic at the project entropy scale.  Equations
(R38.A2), `||B_S||_op<=3||A||_op`, and
`||A||_op^2<=2q_n` give

```math
\operatorname{Var}[\sigma(x')^TB_Sx']
\le72q_n\lVert\rho\rVert_2^2+2n^2.
\tag{R38.A9}
```

Coordinatewise Pinsker gives `d(r)>=r^2/2`, so
`||rho||_2^2<=2D(rho)`.  Thus `D(rho)=O(n^(3/4-c))` makes the variance
`o(q_n^2)`.  After enlarging fixed constants, the second line of (R38.A5)
need not be imposed separately.

The same inequalities create the fatal energy wall.  Since
`||C||_op=||B_S||_op<=3sqrt(2q_n)`,

```math
\boxed{
|\rho^TC\rho|
\le\lVert C\rVert_{op}\lVert\rho\rVert_2^2
\le6\sqrt{2q_n}\,D(\rho).
}
\tag{R38.A10}
```

Reaching even mean energy `q_n` therefore requires

```math
D(\rho)\ge\frac{\sqrt{q_n}}{6\sqrt2}
=\Omega(n^{3/4}).
\tag{R38.A11}
```

This contradicts `D(rho)=O(n^(3/4-c))` for every fixed `c>0`.  Moreover, at
the project budget the mean is `o(q_n)` while (R38.A9) is `o(q_n^2)`, so
Chebyshev upper-bounds the incidence probability by `o(1)`, rather than the
constant required by (R38.A6).  Thus (R38.A5) is a precise sufficient
criterion but an empty one at the desired exponent.  This does not exclude a
different argument which extracts enough uniform incidences from a rare
product-law tail without first making that tail have constant probability.

## 3. The uniform specialization reaches the same barrier

For uniform means `rho_i=sqrt(a)`, based at a positive ground of `B_S` with
cap `M_S`, (R38.A1)--(R38.A4) reduce to the Wave 38 sphere calculation.
Keeping mean energy at least `4q_n` requires

```math
a\ge4q_n/M_S.
```

The available worst-case parent ground-row bound is `R_2(x)<=2nq_n`, so
the row estimate can be certified only when

```math
M_S\gtrsim\frac{nq_n^2}{R}.
\tag{R38.A12}
```

More simply, the triangle inequality gives, for every selector,

```math
\boxed{
M_S=Q(-A+2P_SAP_S)
\le q_n+2Q(A[S])\le3q_n.
}
\tag{R38.A13}
```

Combining (R38.A12)--(R38.A13) requires

```math
R\gtrsim nq_n=\Theta(n^{5/2}).
\tag{R38.A14}
```

Under the explicit hypotheses of the uniform-sphere theorem, `R<=8nq_n`
makes its required threshold at least `16q_n`, while (R38.A13) caps every
structured flip at `3q_n`.  Hence its high-cap subfamily is literally empty,
not merely asymptotically empty.  At the restriction row scale the mismatch
is still larger.

Uniform noise displays the same barrier geometrically.  Anisotropy does
replace the crude cost `a R_2(x)` by the directional term
`||AD_x rho||^2`, but (R38.A10)--(R38.A11) shows that the flipped-energy
entropy cost already fails before that possible row advantage matters.  The
formal quadratic feasibility problem would ask for a bounded vector `rho`
satisfying

```math
\rho^T(\sigma D_xB_SD_x)\rho\ge4q_n,
\qquad
\lVert AD_x\rho\rVert_2^2
+(n-1)(n-\lVert\rho\rVert_2^2)\le R/8,
```

with product entropy `D(rho)=O(R/q_n)`, but (R38.A11) rules this out
universally at a fixed power saving.  This falsifies the displayed
constant-probability independent-product implementation, not the underlying
fractional-cover target, a correlated/conditional noise construction, or a
direct lower bound on a rare product-law tail.
