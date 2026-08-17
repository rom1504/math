# Fixed power-exponential tangent types close only in the Gaussian case

Status: rigorous scoped theorem draft.  The rate-function calculation is
elementary and the probabilistic rigidity step is a direct consequence of
the central limit theorem.  The result falsifies the proposed fixed
power-law tangent carrier; it does not exclude every possible finite
non-Gaussian semigroup.

## 1. Leading roofs close for every power

Fix `p>1` and `a,b>0`.  Write

```math
I_(p,a)(x)=a|x|^p.
```

### Theorem NG.1 (power roofs form an infimal-convolution semigroup)

For every real `z`,

```math
inf_x {a|x|^p+b|z-x|^p}=c|z|^p,                       \tag{NG.1}
```

where

```math
c=\left(a^(-1/(p-1))+b^(-1/(p-1))\right)^(-(p-1)).    \tag{NG.2}
```

Thus the leading exponential roof has a two-parameter finite algebra for
every fixed `p`.

### Proof

By homogeneity it suffices to take `z=1`.  Strict convexity puts the unique
minimizer in `(0,1)`, where

```math
a x^(p-1)=b(1-x)^(p-1).
```

Solving for the two pieces and substituting gives (NG.2).  The cases
`z<=0` follow by symmetry and homogeneity. `square`

## 2. The tangent exponent already changes with the query

For integer `n`, define the full-lattice arrays

```math
A_n^a(k)=exp\{-a|k|^p/n^(p-1)\},\qquad k in Z.         \tag{NG.3}
```

At the central output,

```math
(A_n^a*A_n^b)(0)
=sum_(k in Z)exp\{-(a+b)|k|^p/n^(p-1)\}
\sim n^(1-1/p)\int_R e^{-(a+b)|u|^p}du.               \tag{NG.4}
```

At a fixed nonzero macroscopic output `floor(nz)`, the unique saddle in
(NG.1) has a positive finite Hessian, and ordinary one-dimensional Laplace
asymptotics instead give a factor `Theta(n^(1/2))`.  Unless `p=2`, one
power-law roof therefore has at least two tangent-mass exponents:

```math
1-1/p\quad(z=0),\qquad 1/2\quad(z!=0).                 \tag{NG.5}
```

Equation (NG.4) follows directly by the Riemann sum with mesh
`n^{-(1-1/p)}`.  The off-centre assertion follows from strict convexity and
the nonzero optimal split.  This is a finite stratification, but the next
theorem shows that its exact central tangent shapes still do not close.

## 3. Gaussian rigidity of the power-exponential tangent family

Let

```math
g_(p,a)(x)=Z_(p,a)^(-1)e^(-a|x|^p)                    \tag{NG.6}
```

be the normalized continuous tangent density.

### Theorem NG.2 (closure occurs only at the Gaussian power)

If, for some `p>1` and `a,b>0`,

```math
g_(p,a)*g_(p,a)=g_(p,b),                              \tag{NG.7}
```

then `p=2` and `b=a/2`.  Conversely this identity holds for `p=2`.

Consequently, for `p!=2`, the finite tuple consisting of the leading
power-roof coefficient, tangent exponent, and scalar amplitude is not an
exact reusable convolution state.  Already the self-composition produces a
new tangent shape.

### Proof

Let `X,X_1,X_2` be iid with density `g_(p,a)`.  Every member of (NG.6) is a
scale copy of every other member with the same `p`, so (NG.7) says

```math
X_1+X_2 =_d sX
```

for some `s>0`.  These variables have finite nonzero variance, hence
`s=sqrt(2)`.  Iteration gives

```math
2^(-r/2)(X_1+...+X_(2^r)) =_d X                      \tag{NG.8}
```

for every `r`.  The central limit theorem says the left side converges in
distribution to the centered Gaussian with the variance of `X`.  Therefore
`X` itself is Gaussian.  Comparing the logarithm of its positive density
with `-a|x|^p-\log Z_(p,a)` forces `p=2`.  Variance, or direct Gaussian
convolution, then gives `b=a/2`.  The converse is the usual Gaussian
convolution identity.  This rigidity step is the classical Gaussian
stability consequence of the central limit theorem. `square`

The continuous tangent law is genuinely inherited by the discrete arrays.
Writing `s_n=n^(1-1/p)`, a dominated Riemann-sum argument gives, locally
uniformly for real `u`,

```math
{1\over s_n}(A_n^a*A_n^a)(floor(s_n u))
\longrightarrow
\int_R e^(-a|v|^p)e^(-a|u-v|^p)dv.                  \tag{NG.10}
```

Indeed put `k=s_nv`; the mesh is `1/s_n`, and on every compact `u`-set the
summands have an integrable uniform power-exponential envelope.  Thus the
failure of (NG.7) is also a failure of the corresponding discrete tangent
carrier, not merely an analogy between continuous densities.

## 4. What the falsifier teaches

The leading rate function in Theorem NG.1 discards the normalized central
tangent density.  That density is precisely what the next convolution uses.
For `p!=2`, repeated composition generates the convolution powers

```math
g_(p,a_1)*...*g_(p,a_m),                              \tag{NG.9}
```

not a fixed power-exponential family.  One may retain the whole density as
an infinite-dimensional carrier, or retain the entire list of factors, but
neither is the sought strict finite state.

The classical central limit theorem does give an **approximate typical-scale
collapse** after many comparable factors.  It does not control the
large-deviation tails or all exposed macroscopic queries, where Theorem NG.1
and query-dependent saddle geometry remain active.  Hence typical tangent
mixing is not automatically an extremal response quotient.

The project-level conclusion is more than the earlier isolated quartic
exponent calculation:
it identifies Gaussian stability as the exact reason Theorem 32.2 closes and
rules out the entire fixed-`p` generalized-Gaussian family as the next finite
extension.  It does not prove that no other finite stratified non-Gaussian
family exists.
