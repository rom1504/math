# Wave 51C: spectral-excess extraction and the exact cap obstruction

## Status

No asymptotic overlap theorem or convergence proof is obtained.  The main
advance is a verified sharpening of the coarea target: triple retention above
`lambda_1` is much stronger than necessary.  A project-row class already
supplies a useful center when its aggregate retention exceeds `lambda_2` by
the required saved amount.  This changes the repair problem materially.

The derivations below also show that, for this sharper certificate, adding
`K_0` can never improve on the better of the pure kernel and the tautological
cubic-ratio bound.  Finally, a scalable abstract singleton-family model shows
that box positivity and cap nesting alone do not force spectral excess.  That
model is **not** asserted to arise from an exact minimizing signing.

All finite arithmetic is reproduced by
`tmp/coarea_repair_r51_check.py`; its captured output is
`tmp/coarea_repair_r51.out`.

## 1. Exact spectral-excess extraction

Fix a selector-independent center class `C`, a common-core kernel `K`, and
write

```math
D_C=\mathbb E[\mathbf1_Ca_z^2]>0,\qquad
P=\frac{\mathbb E[\mathbf1_Ca_z^2p_z]}{D_C},
\qquad M=\max_{z\in C}a_z,
```

where `p_z=<f_z,Kf_z>/a_z`.  Put `delta=1-lambda_1`,
`gamma=lambda_1-lambda_2`, and `kappa=gamma/delta`.  Since
`B_z=a_z(1-p_z)`, the aggregate boundary ratio is exactly

```math
\mathcal R_C
=\frac{\mathbb E[\mathbf1_Ca_zB_z]}{\delta D_C}
=\frac{1-P}{\delta}.                                      \tag{R51C.1}
```

Averaging the verified pointwise inequality (10.1220) under the law
proportional to `1_C a_z^2`, whose mean degree is at most `M`, gives

```math
\mathcal R_C\ge1+\kappa-(1+\kappa n)M.
```

Substitution of (R51C.1) and exact cancellation prove the **Verified
spectral-excess bound**

```math
\boxed{
M\ge
\frac{[P-\lambda_2]_+}
{(1-\lambda_1)+n(\lambda_1-\lambda_2)}.}                  \tag{R51C.2}
```

This does not appear elsewhere in the ledger.  The older one-step collision
criterion (10.1176) uses the `lambda_1` threshold; (R51C.2) gains the interval
between `lambda_2` and `lambda_1` by using the sharp Boolean level-one bound
(10.1219).  Equation (10.1221) is a special stronger-margin consequence, not
the optimized statement (R51C.2).

For the positive semidefinite common-core kernels the denominator in
(R51C.2) is at most `n+1`.  Therefore, at a project-row cap, the exact
convergence-scale sufficient target is only

```math
P_\ell(R)-\lambda_2(\ell)
\ge \exp\{-O(n^{3/4-c})\}.                               \tag{R51C.3}
```

Together with nonzero box mass, (R51C.2) gives a center with the row bound
and selector degree `exp{-O(n^(3/4-c))}`, hence the bare restriction tail.
No constant `kappa`, non-strict coarea inequality, or slice FKN theorem is
needed for this implication.

The Wave 50 cap-ten example illustrates the correction sharply.  Its degree
is `M=1/42`.  Although all five positive-core scales have `P<lambda_1`,
scale one has `P-lambda_2=1/42` and (R51C.2) gives `M>=1/70`; scale two has
`P-lambda_2=1/75` and gives `M>=14/2475`.  Thus that example is a wall only
for the overstrong non-strict coarea/FKN package, not for the actual
spectral-excess target.

## 2. `K_0` admixture is redundant for (R51C.2)

Let `r=P_0=sum a_z^3/sum a_z^2` and mix a positive-core kernel with
independent resampling:

```math
K_\theta=(1-\theta)K_0+\theta K.
```

Its retention and two leading nonconstant eigenvalues are

```math
P_\theta=(1-\theta)r+\theta P,\qquad
\lambda_{j,\theta}=\theta\lambda_j\quad(j=1,2).
```

Applying (R51C.2), the resulting degree lower bound is

```math
L(\theta)=
\frac{(1-\theta)r+\theta(P-\lambda_2)}
{(1-\theta)+\theta\{(1-\lambda_1)+n(\lambda_1-\lambda_2)\}}. \tag{R51C.4}
```

The right side is a denominator-weighted convex combination of the endpoint
bounds

```math
L(0)=r,\qquad
L(1)=\frac{P-\lambda_2}
{(1-\lambda_1)+n(\lambda_1-\lambda_2)}.
```

Consequently the **Verified exact optimization** is

```math
\boxed{\sup_{0\le\theta\le1}L(\theta)=\max\{L(0),L(1)\}.} \tag{R51C.5}
```

The `K_0` endpoint is merely `M>=r`, already true by definition.  Hence
`K_0` mixing supplies no new spectral-excess mechanism; the live question is
whether a pure positive-core kernel has saved excess above `lambda_2`.

For comparison, if one insists on the stronger condition `P_theta >=
lambda_(1,theta)` and a uniform harmonic ratio, put
`d=lambda_1-P>0`.  Retention forces
`theta<=theta*=r/(r+d)`, and monotonicity of the harmonic ratio gives

```math
\kappa_{\max}
=\frac{r(\lambda_1-\lambda_2)}
{d+r(1-\lambda_1)}.                                    \tag{R51C.6}
```

Thus a target `kappa_0` is attainable exactly when

```math
d\le r\left\{\frac{\lambda_1-\lambda_2}{\kappa_0}
-(1-\lambda_1)\right\}.
```

This recovers the Wave 50 degree-circularity warning and identifies its exact
missing relative-deficit inequality, but (R51C.2)--(R51C.5) make that whole
stronger repair package unnecessary for convergence.

## 3. What cap inflation does, and what it does not do

For a signing `A`, a full cut `z`, and a coordinate set `U`, let `z^U` flip
the spins in `U`.  If `r=|U|`, then

```math
Az^U=Az-2A[:,U]z_U,\qquad
\lVert A[:,U]z_U\rVert_2\le r\sqrt n,
```

and hence the **Verified row-neighborhood bound**

```math
\boxed{R_2(z^U)\le(\sqrt{R_2(z)}+2r\sqrt n)^2.}          \tag{R51C.7}
```

If `z[S]` is a child ground and `U subset S^c`, then `z^U[S]=z[S]`, so the
same selector remains favorable.  Thus controlled cap inflation includes a
sublinear-radius ball of outside completions around a box witness.  This
increases the number of centers but need not create a second favorable
selector for any one center, which is what spectral excess needs.

The exact cap bookkeeping makes the missing input explicit.  For active row
levels `R_j`, define

```math
E_\ell(R_j)=\sum_{R_2(z)\le R_j}a_z^2
             \{p_\ell(z)-\lambda_2(\ell)\},\qquad
D(R_j)=\sum_{R_2(z)\le R_j}a_z^2.
```

Then (R51C.3) is `E_ell(R_j)/D(R_j)>=exp{-O(H)}`.  Passing to the next cap
only adds the signed shell sum

```math
E_\ell(R_j)-E_\ell(R_{j-1})
=\sum_{R_2(z)=R_j}a_z^2\{p_\ell(z)-\lambda_2(\ell)\}.    \tag{R51C.8}
```

Neither `D(R)>0` nor its monotonicity controls the sign in (R51C.8).  A valid
telescoping proof must establish a signing-specific positive shell budget or
directly force the normalized excess at some prefix.

At `ell=1`, write `r_{z,i}=|{S in F_z:i in S}|`.  Exact common-core counting
and `sum_i r_(z,i)=m r_z` give

```math
p_1(z)=a_z\,
\frac{n\sum_i r_{z,i}^2}{m^2r_z^2},\qquad
1\le\frac{p_1(z)}{a_z}\le\frac nm=\frac1p
```

with the last display interpreted as `n/m=1/p`.  Therefore

```math
r_C\le P_1(R)\le r_C/p.                                 \tag{R51C.9}
```

Since `lambda_2(1)=0`, scale one is quantitatively equivalent, up to a
density constant, to the already tautological cubic ratio.  A genuinely new
overlap mechanism must use `ell>=2` and produce excess over its positive
`lambda_2`.

## 4. A precise abstract obstruction (not an exact-minimizer example)

Take `m=3n/5`, fix one selector `S_*`, and let every center in a low-row
cluster have the singleton favorable family `F_z={S_*}`.  Put any additional
centers needed to cover the other selectors above the allowed cap.  The
low-row cluster may contain the whole outside-flip ball forced by (R51C.7);
all of those centers still have the same singleton family.  Hence every
allowed cap has

```math
a_z=\binom nm^{-1},\qquad
p_\ell(z)=h_\ell:=\binom{n-\ell}{m-\ell}^{-1}.           \tag{R51C.10}
```

For fixed density and all sufficiently large `n`, exact factorial
cancellation gives

```math
\frac{h_1}{\lambda_1}
=\frac{m}{\binom{n-2}{n-m-1}},
```

so `h_1=exp{-Theta(n)}`.  For `ell>=2`, similarly

```math
\frac{h_\ell}{\lambda_2(\ell)}
=\frac{m(m-1)}
{\ell(\ell-1)\binom{n-\ell-2}{n-m-2}}\le1,              \tag{R51C.11}
```

where the final inequality holds uniformly for large `n` at this fixed
density.  One elementary uniform proof splits at
`ell=O(sqrt(n))`: above that point the binomial is at least `n-m-1`, which
already makes the denominator dominate `m(m-1)`; below it, both arguments of
the binomial remain linear in `n`, so its standard entropy lower bound is
exponential.  Thus scale one
has only exponentially small positive excess, and every scale `ell>=2` has
nonpositive excess.  The model defeats (R51C.3) at every scale despite box
positivity, cap nesting, selector coverage, and the forced outside-flip
cluster.

This is only a falsification of incidence-only or cap-monotonicity proofs.
Real ground families of exact minimizing signings obey additional coupled
energy constraints, and no scalable exact-minimizer obstruction is known.
The sharply reduced **Open target** is therefore:

> From exact-minimizer structure plus a project-scale box witness, prove that
> some controlled cap and some `ell>=2` have
> `P_ell(R)-lambda_2(ell)>=exp{-O(n^(3/4-c))}`.

Equivalently, prove a positive signed-shell budget in (R51C.8).  Any proof
using only additional low-row completions of the same child ground cannot
work; it must force overlap among distinct favorable selectors for a common
center.
