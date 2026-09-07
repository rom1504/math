# Independent audit: anisotropic flip noise

## Verdict

The identities (R38.A1)--(R38.A4), the entropy-to-cardinality implication
(R38.A7), the dual double count (R38.A8), and the variance estimate (R38.A9)
are correct.  Both exact checkers pass.

However, the proposed anisotropic certificate (R38.A5) is asymptotically
empty at the required entropy budget.  The same spectral ceiling that kills
uniform noise also kills every heterogeneous product law once its relative
entropy is kept at the project exponent.  Therefore the memo's claim that
anisotropy makes the argument nonvacuous must be replaced by an obstruction.

## 1. Exact identities

For `z=rho+eta`, with independent centered `eta_i` and
`v_i=E eta_i^2=1-rho_i^2`, symmetry and the zero diagonal of `C` give

```math
\mathbb E z^TCz=\rho^TC\rho.
```

Moreover

```math
z^TCz-\mathbb Ez^TCz
=2\eta^TC\rho+2\sum_{i<j}C_{ij}\eta_i\eta_j.
```

The two terms are orthogonal.  Distinct quadratic monomials are orthogonal,
and `C_(ij)^2=1`, so

```math
\operatorname{Var}(z^TCz)
=4\sum_i v_i(C\rho)_i^2+4\sum_{i<j}v_iv_j.
```

This proves (R38.A1)--(R38.A2).  For the parent row square,

```math
\mathbb E\lVert AD_xz\rVert_2^2
=\lVert AD_x\rho\rVert_2^2
+\sum_i v_i\lVert Ae_i\rVert_2^2
=\lVert AD_x\rho\rVert_2^2+(n-1)\sum_i v_i,
```

which is (R38.A3).  Product additivity gives (R38.A4).

The exact rational checker
`tmp/anisotropic_flip_noise_r38_check.py` passes on all three test cases.
The corrected uniform checker `tmp/dual_flip_compression_r38_check.py` also
passes.

## 2. Entropy cardinality and dual count

If `u=U(G)` and `P(G)>=alpha`, binary data processing gives

```math
D(P\Vert U)
\ge-h(\alpha)-\alpha\log u-(1-\alpha)\log(1-u)
\ge-h(\alpha)-\alpha\log u.
```

This rearranges to (R38.A7).  At `alpha=1/2` and `D(P||U)<=L`, the number
of good full-spin representatives at the selected orientation is at least

```math
K\ge 2^n e^{-2L}/4.
```

The good set is invariant under global spin reversal, so it represents at
least `K/2` actual projective states at that orientation.  Summing a feasible
dual over all `2^n` actual oriented states gives

```math
\frac K2\sum_Sy_S
\le\sum_d\sum_{S\in I_d}y_S
\le2^n.
```

Consequently `sum_S y_S<=8e^(2L)`, which is exactly (R38.A8) up to its
stated `O(1)` constant.

## 3. Variance bound

Since `C` is orthogonally similar to `B_S`,

```math
\lVert C\rVert_{op}
=\lVert B_S\rVert_{op}
\le3\lVert A\rVert_{op}
\le3\sqrt{2q_n}.
```

Thus the first term in (R38.A2) is at most

```math
4\lVert C\rho\rVert_2^2
\le72q_n\lVert\rho\rVert_2^2,
```

and the second is below `2n^2`.  This proves (R38.A9).  The displayed
assumption `max_i |rho_i|<=1/2` is not needed for this variance inequality.

## 4. Decisive entropy--spectral obstruction

For the one-bit entropy in (R38.A4),

```math
d''(r)=\frac1{1-r^2}\ge1,
\qquad d(0)=d'(0)=0.
```

Therefore

```math
\boxed{D(\rho)=\sum_i d(\rho_i)\ge\frac12\lVert\rho\rVert_2^2.}
\tag{Audit.1}
```

Combining (Audit.1) with the same spectral estimate gives

```math
\rho^TC\rho
\le\lVert C\rVert_{op}\lVert\rho\rVert_2^2
\le6\sqrt{2q_n}\,D(\rho).
\tag{Audit.2}
```

Hence the first clause of (R38.A5), `rho^TC rho>=4q_n`, necessarily implies

```math
\boxed{
D(\rho)\ge\frac{2}{3\sqrt2}\sqrt{q_n}
=\Omega(n^{3/4}).
}
\tag{Audit.3}
```

At the project cutoff

```math
R=O(n^{9/4-c}),\qquad q_n=\Theta(n^{3/2}),
```

the requested entropy budget is only

```math
L=O(R/q_n)=O(n^{3/4-c})=o(\sqrt{q_n}).
```

For every fixed `c>0`, (Audit.3) contradicts this budget for all sufficiently
large `n`.  This is independent of the row and variance clauses: no selector
can be anisotropically regularizable in the stated sense.

The same conclusion holds if the constant mean `4q_n` is replaced by any
fixed constant greater than one times `q_n`.  Escaping it would require a
non-product/counting argument not certified through an `O(R/q_n)` product
relative entropy, or a structural improvement to the spectral comparison
itself.

## 5. Uniform spectral comparison

The corrected uniform memo has the right scope.  Its conditional theorem
requires

```math
t\ge K\frac{nq_n^2}{R},
```

whereas every structured complement flip satisfies

```math
M_S\le n\lVert B_S\rVert_{op}
\le3n\sqrt{2q_n}.
```

Thus

```math
\frac{t}{\max_SM_S}
\ge c_1\frac{q_n^{3/2}}R
=\Omega(n^c)
```

when `R=O(n^(9/4-c))`.  The high-cap class is empty eventually.  Equivalently,
compatibility of the two inequalities would require
`R=Omega(q_n^(3/2))=Omega(n^(9/4))`, with no power saving.  The memo now
correctly labels the finite conditional theorem asymptotically vacuous and
does not claim that it narrows the remaining selector family.
