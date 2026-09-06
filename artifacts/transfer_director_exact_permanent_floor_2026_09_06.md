# An exact finite permanent floor for the full-spin weave certificate

Date: 2026-09-06. **Proved certificate obstruction**, independently checked
by the reconstruction track. This is NOT a lower bound on actual signings
and does not change the rigorous original interval. Unlike a floor for an
approximate Bellman supersolution, it applies to the exact row quantities
before any recursive or Gaussian-boundary upper approximation.

## 1. Finite orbital inequality

For `t>0` and `v in R^m`, put

```math
K_{ij}=\tfrac12\{e^{-t(v_i-v_j)^2}+e^{-t(v_i+v_j)^2}\},\qquad
P_t(v)=\frac{\operatorname{per}K}{m!},\qquad L_t(v)=\sqrt{P_t(v)}.
```

Let `nu_v=(1/(2m)) sum_i(delta_(v_i)+delta_(-v_i))`, and use the
self-transport and Gaussian potential from the standalone proof:

```math
F_t(\nu)=\inf_{\gamma\in\Pi(\nu,\nu)}
\{D(\gamma\Vert\nu^{\otimes2})+t\mathbb E(X-Y)^2\},\quad
\Phi_t=-F_t/2,\quad g_t(e)=\Phi_t(N(0,e)).
```

Then the following inequalities hold at EVERY finite order, without
bounded-coordinate or fixed-alphabet assumptions:

```math
\boxed{\quad L_t(v)\ge e^{m\Phi_t(\nu_v)}
                  \ge e^{m g_t(\|v\|^2/m)}.\quad}                 (1)
```

To prove the first inequality, first replace coordinates by their absolute
values. Folding the two signs in the transport problem gives exactly

```math
F_t(\nu_v)=\min_{\gamma\mathbf1=\gamma^T\mathbf1=\mathbf1/m}
             \sum_{ij}\gamma_{ij}\log\frac{m^2\gamma_{ij}}{K_{ij}}.
```

Indeed, for each pair of magnitudes, minimize the conditional relative-sign
cost against a fair sign; its partition function is `K_ij`. A common fair
sign enforces both signed marginals. Duplicate coordinates and zero atoms
may be kept as distinct index labels: uniform conditional index allocation
has no extra information cost, and data processing proves the reverse
comparison.

The kernel is strictly positive. The finite minimizer is positive (mixing
with the product law has entropy derivative minus infinity at a zero),
unique, and symmetric. Its multipliers therefore give

```math
\gamma_{ij}=m^{-2}K_{ij}e^{f_i+f_j},\quad
D_{ij}=m\gamma_{ij},\quad
F_t(\nu_v)=\frac2m\sum_i f_i.
```

Here `D` is entrywise positive and doubly stochastic, with row and column
sums ONE. The van der Waerden permanent inequality gives

```math
\frac{\operatorname{per}K}{m!}
=\frac{m^m}{m!}\operatorname{per}D\ e^{-mF_t(\nu_v)}
\ge e^{-mF_t(\nu_v)}.
```

The imported theorem and its normalization are stated in Gurvits,
*Van der Waerden/Schrijver-Valiant like Conjectures and Stable Homogeneous
Polynomials*, Corollary 2.5 and Example 2.6(1), pp. 8–9 in the PDF
[arXiv:0711.3496v2](https://arxiv.org/pdf/0711.3496).
For completeness, the polynomial to which that result applies is
`p_D(z)=prod_i sum_j D_ij z_j`: it is homogeneous of degree `m`, has
nonnegative coefficients, and is stable in the right half-plane because
each factor has positive real part there. Weighted AM–GM gives capacity
at least one; `z=1` gives equality. Its full squarefree coefficient is
`per D`, so the corollary yields precisely `m!/m^m`.

The second inequality in (1) is the already reconstructed Gaussian
extremality `Phi_t(nu)>=g_t(Var nu)` for symmetric sources; its proof uses
orthogonal tensorization, independent averaging, and the finite-variance
central limit theorem. Here `Var nu_v=||v||^2/m` exactly. No empirical
asymptotic or exchange of limits occurs in (1).

## 2. Deletion and the exact Hadamard row sum

The Gaussian Fock feature is a tensor product under a coordinate split.
Invariants of the full signed-permutation group form a subspace of the
invariants of the subgroup which preserves that split. Orthogonal
projection therefore gives

```math
L_t(v)\le L_t(v\setminus v_i)L_t(v_i)\le L_t(v\setminus v_i).       (2)
```

If `H` is any Hadamard of order `m`, `T` retains any `k` rows, and
`x in {+-1}^k`, then `v=H[T,:]^T x/sqrt(k)` has squared norm `m`.
Thus the EXACT row partition appearing in the original full-spin bound
satisfies

```math
Z_H=\sum_x\max_i L_t(v\setminus v_i)
       \ge 2^k e^{m g_t(1)}.                                  (3)
```

This is pointwise in each basis and selector. No randomness, independence,
recursion depth, Gaussian approximation, or optimality of a seed is needed.

## 3. Extension to arbitrary sign bases, with variable energy retained

For an arbitrary `k by m` sign matrix `H_T`, not necessarily orthogonal,
define the tilted quantities

```math
L_t^+(v)=e^{t\|v\|^2}L_t(v),\qquad
Z_H^+=\sum_{x\in\{\pm1\}^k} L_t^+(H_T^T x/\sqrt k).
```

For uniform independent spins, `E_x ||H_T^T x/sqrt(k)||^2=m`, because
the sum of squared entries of `H_T` is exactly `km`. The function
`e -> t e+g_t(e)` is convex. Applying (1) and then Jensen twice proves

```math
\boxed{\quad Z_H^+\ge 2^k\exp\{m[t+g_t(1)]\}.\quad}             (4)
```

Thus nonorthogonal input signings cannot evade this floor merely by
retaining the fluctuating norm in the row partition. Orthogonality is
not a hypothesis of (4).

## 4. A sharp floor for this certificate, not for the constructed matrices

Put `p=k/m`. To certify full quadratic cap at most `c (mk)^(3/2)`
(the convention `Q=|x^TAx|/2`), the original Hadamard certificate uses

```math
2\exp\{t(1-2c\sqrt p)m^2\}\prod_{i=1}^m Z_{H_i}.               (5)
```

Equation (3) bounds (5) BELOW by

```math
2\exp\{m^2[p\log2+g_t(1)+t(1-2c\sqrt p)]\}.                  (6)
```

The arbitrary-sign tilted certificate has the same lower exponent, by
(4), with the additional positive factor `(2m)^(m/2)`. Averaging correlated
or independent choices of bases cannot reduce these pointwise lower bounds.

The exact optimization is

```math
\inf_{0<p\le1,\ t>0}
\frac{p\log2+t+g_t(1)}{2t\sqrt p}
=\frac{\sqrt{15}}8.                                         (7)
```

Write `2t=rho/(1-rho^2)`, `0<rho<1`. Differentiation gives
`t g_t'(1; t)-g_t(1)=-(1/4)log(1-rho^2)`, where the derivative here
is with respect to TEMPERATURE `t`, not variance. Equivalently, differentiate
the explicit Gaussian formula: `d_t[t+g_t(1)]=rho`.
For fixed `p`, the unique minimizing temperature obeys
`1-rho^2=2^(-4p)` and the minimum is

```math
\frac{\sqrt{1-2^{-4p}}}{2\sqrt p}.
```

Its square decreases in `p`: for `a=4 log2`,
`d_p[(1-e^(-ap))/p]=[e^(-ap)(1+ap)-1]/p^2<0`.
It is therefore minimized at `p=1`, where `rho=sqrt(15)/4` and
`t=2sqrt(15)`. This proves (7).

For EVERY finite `m,k,t`, if `c<=sqrt(15)/8`, exponent (6) is
nonnegative. The right side of the full-spin certificate is at least
two and cannot prove a positive probability of success. This statement
even permits `p,t` to vary with order. Hollowing changes the normalized
objective by only `O((mk)^(-1/2))` and cannot remove a fixed asymptotic floor.

## 5. The precise next obligation

The proved all-order upper bound below `0.499433` remains valid and is
strictly above this floor. Nothing here says the actual weave has cap at
least `sqrt(15)/8`, nor that the original liminf is below it. The equality
in (7) is an optimum of a lower bound on a proof certificate, not a
construction attaining that value.

This removes the possibility of obtaining an arbitrarily seed-preserving
theorem from BETTER ROW POTENTIALS ALONE within the same exact certificate.
To transfer seeds whose normalized cap is below this floor, one would need
an exponentially significant improvement before or instead of rowwise
PSD/Finner factorization or the full-spin union bound. For example, at
fixed `p,t,c`, improving the logarithm of (5) by `Delta m^2` can only help
when it beats the nonnegative exponent in (6); a loss `o(m^2)` cannot
remove a fixed positive difference in this exponent. Such an improvement
must be proved for the actual joint spin/row law, not assumed from a
seed's small cap. The original convergence problem remains open.

`computations/transfer_director_permanent_floor_checks_2026_09_06.py`
provides exact rational folded-permanent and projection checks, with
separately labeled high-precision numerical sanity checks of (1) and (7).
