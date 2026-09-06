# A joint nonlinear response direction and an explicit finite-degree bound

Date: 2026-09-05. This combines the all-odd weighted theorem into one
normalized direction. The channels are combined before either clipping
or estimating the Boolean energy.

## 1. The finite joint projection theorem

Let `F` be a fixed bounded odd response of a finite old Gaussian tree
family, with edge coordinate `G0`. Fix a finite set `D` of odd integers
at least three and write

```
f_r=E[F h_r(G0)],   s_D^2=sum_(r in D) f_r^2 > 0,
h(g)=sum_(r in D) (f_r/s_D) h_r(g),
R_h=sum_(r in D) (f_r^2/s_D^2) Q^(circ r),
Z=B h(BS),   v_i=(B R_h B)_ii.
```

For a fixed bounded local response `M` and bounded deterministic weights
`d_i`, linearity of the all-odd weighted theorem gives

```
n^-1 E[d circ M(X) circ Z]^T B F(X)
 = E M * s_D * n^-1 sum_i d_i v_i + o(1).             (1)
```

Indeed every term has coefficient `f_r^2/s_D`, exactly `s_D` times its
weight in `R_h`. This identity does not estimate each channel separately.
The matrix is an actual hollow signing, with `||B||op<=L` fixed, and the
finite family and degrees are fixed before its order grows.

The finite vector of unmarked chaos components has asymptotically joint
Gaussian marginals: each component has vanishing contractions, distinct
degrees are orthogonal, and all have vanishing covariance with the old
family. Thus `(X_i,Z_i)` has limiting law
`(gamma,sqrt(v_i)N0)`, with independent `N0`. This follows from the same
multivariate fixed-chaos theorem used for each single degree, not from
an assertion that uncorrelated random variables are independent.

For every odd `r`, Schur Jensen gives
`tr Q sqrt(Q^(circ r))>=n`. Operator concavity of the square root and
the convex weights defining `R_h` therefore give

```
tr Q sqrt(R_h)>=n,  sum_i sqrt(v_i)>=n.                (2)
```

Consequently `d_i=1/sqrt(max(v_i,eta))` gives local variance at most one
and `n^-1 sum_i d_i v_i>=1-sqrt(eta)`. The useful normalized response
in (1) is the full nonlinear conditional Hermite mass `s_D`, rather than
the largest individual Hermite coefficient. The previously proved local
slack, clipping, and exact paired-means argument then applies unchanged.

## 2. An explicit degree choice for scalar-mask responses

Consider `F=sign(W)1{|V|>alpha}` for a centered jointly Gaussian pair
with `Var(V)=1` and `Var(W)>0`, and an additional standard Gaussian `G0`.
No nondegeneracy of the entire triple is assumed. Put

```
a=E G0 F,  s^2=sum_(r>=3, r odd) |E F h_r(G0)|^2.
```

If `|a|>=a0>0`, boundedness of `f(g)=E[F|G0=g]` and Hermite completeness
imply the explicit positive bound

```
s^2=E(f(G0)-aG0)^2
   >= u(a0):=E(a0 |G0|-1)_+^2 >0.                   (3)
```

There is also a dimension-free Hermite tail estimate

```
sum_(r>D) |E F h_r(G0)|^2 <= 5/sqrt(D),  D>=2.       (4)
```

Here is a direct proof with deliberately loose constants. Couple two
copies of the entire Gaussian input by correlation `q` in `[0,1)`, and
write `theta=arccos(q)`. A threshold indicator of a standard Gaussian
changes with probability at most `theta/pi`. To see this without a
stability theorem, represent the two coordinates as projections of a
rotationally invariant planar Gaussian onto axes separated by `theta`.
Conditioning on its radius makes the event an angular arc (or an empty
or full set); rotating an arc changes at most `2theta` angular measure.

The sign of `W` changes with probability `theta/pi`. The central interval
indicator in `V` has two thresholds, so changes with probability at most
`2theta/pi`. If `F,F'` are the two responses, it follows pointwise that

```
(F-F')^2 <= 4*1{sign(W) differs}+1{central mask differs},
E(F-F')^2 <= 6theta/pi.
```

The Gaussian Hermite expansion, including all input coordinates, gives
`2 sum_k(1-q^k)||F_k||_2^2=E(F-F')^2`. The coefficients involving only
`G0` are a subset of its orthogonal coefficients. With `q=1-1/D`,
`1-q^(D+1)>=1-e^-1` and `theta<=2/sqrt(D)`. Thus the left side of (4)
is at most `6/[pi(1-e^-1)sqrt(D)]<5/sqrt(D)`.

Choosing any integer `D >= (10/u(a0))^2` therefore guarantees

```
s_D^2=sum_(3<=r<=D, r odd) |E F h_r(G0)|^2 >= u(a0)/2. (5)
```

In the near-optimal scalar hierarchy the previously proved uniform
linear coefficient is at least `a0=1/200`. Formula (5) gives an explicit,
although extremely large, common degree bound in place of the compactness
selection argument. The degree is fixed independently of matrix order.
This is not a claim of computational practicality or a sharp degree rate.

Subsequent improvement: `fresh_general_mask_nonlinear_direction_2026_09_06.md`
and its independent monomial audit replace this very large degree bound
by degree 417 for every marked mask with certificate at least .43.
The threshold-noise estimate above remains valid but is not needed for
that stronger finite-selection conclusion.

## 3. Scope

The analytic ingredients in (3)--(4) are elementary conditional projection
and Gaussian noise estimates. Their role here is quantitative: they give
a single finite jointly evaluated response with a uniform positive energy
direction, which survives variance normalization and arbitrary fixed
operator caps. They do not show that the enlarged response family attains
the true optimum, nor provide all-order upper recovery or convergence.

Dependencies: `fresh_all_odd_weighted_projection_2026_09_05.md`,
`fresh_limit_odd_schur_standard_deviation_2026_09_05.md`, and
`fresh_uniform_scalar_hierarchy_escape_2026_09_05.md`.

Independent audit: the literature agent checked the complete finite mixture,
the Gaussian threshold-rotation estimate, the Hermite-tail constants, and
the ordering of fixed degree and matrix size. The director reconstructed
the argument. No extra numerical increment is claimed from the enormous
explicit degree bound.
