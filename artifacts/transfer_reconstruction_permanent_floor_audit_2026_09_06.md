# Independent audit of the finite permanent floor and the full-row barrier

Date: 2026-09-06. The director's proposed finite floor and its
`sqrt(15)/8` certificate obstruction are correct. No matrix-order,
recursion-depth, typical-profile, or uniform-CLT limit is needed for the
floor. This is not a lower bound on actual quadratic signings.

## 1. The finite normalization

For `v in R^m`, `t>0`, define

```math
K_{ij}=\frac12\left(e^{-t(|v_i|-|v_j|)^2}
                        +e^{-t(|v_i|+|v_j|)^2}\right),\qquad
L_t(v)^2=\frac{\operatorname{per}K}{m!}.
```

Let `nu` be the symmetrized empirical law of `v`, including its zero
mass without doubling it. Folding relative signs in its self-transport
problem gives the finite transport problem on the uniform coordinate
indices with cost `-log K_ij`. Duplicated magnitudes are harmless:
refining an optimal value coupling uniformly over its repeated indices
preserves the value, and index-to-value data processing gives the reverse
comparison. At a zero magnitude the sign cost is constant, so no phantom
sign entropy is charged.

Strictly positive finite-kernel entropy minimization supplies a symmetric
scaling

```math
\gamma_{ij}=\frac{K_{ij}e^{f_i+f_j}}{m^2},\qquad
\sum_j\gamma_{ij}=\frac1m,\qquad
F_t(\nu)=\frac2m\sum_i f_i.
```

Thus `D=m gamma`, namely `D_ij=K_ij e^(f_i+f_j)/m`, is a positive
bistochastic matrix. The permanent scaling identity is exactly

```math
\frac{\operatorname{per}K}{m!}
=\frac{m^m}{m!}\operatorname{per}D\,e^{-mF_t(\nu)}.
\tag{1}
```

The sole imported permanent theorem is `per D>=m!/m^m` for a
nonnegative bistochastic matrix. It is explicitly obtained from the stable
polynomial bound in Gurvits, Corollary 2.5 and Example 2.6(1), printed page
7. For the product polynomial `product_i sum_j D_ij z_j`, positivity of
the real parts gives stability, the coefficient of `product_j z_j` is
the permanent, and weighted AM--GM gives capacity one.
[Gurvits, Electronic Journal of Combinatorics 15 (2008), R66](https://emis.dsd.sztaki.hu/journals/EJC/Volume_15/PDF/v15i1r66.pdf).

Combining (1) with the independently reconstructed Gaussian maximum of
self-transport gives the exact pointwise floor

```math
\boxed{\quad L_t(v)\ge e^{m\Phi_t(\nu)}
\ge \exp\left[m g_t\left(\frac{\|v\|^2}{m}\right)\right].\quad}
\tag{2}
```

No `O(log m)` or `o(m)` term occurs. The direction is a lower bound on
`L`, not an upper estimate on a partition sum.

## 2. Removing the diagonal coordinate does not defeat the floor

The invariant projection for all signed permutations has range contained
in that for a coordinate split. Therefore

```math
L_t(v)\le L_t(v\setminus i)L_t(v_i)\le L_t(v\setminus i),
\tag{3}
```

because `L_t(v_i)^2=(1+e^{-4tv_i^2})/2<=1`. This is a different
inequality from the polynomial-loss upper estimate on deletion; both
directions are valid in their respective senses. Formula (3) retains the
original dimension `m` in (2), avoiding any assumption that the removed
coordinate has small energy. Empty lists have `L=1`.

For every order-`m` Hadamard fibre and every retained spin on `k` rows,
the normalized full row spectrum `v=H[T,:]^T x/sqrt(k)` has
`||v||^2=m`. Thus every row partition sum in the original full-spin
Finner/Markov certificate obeys, deterministically,

```math
Z_i=\sum_{x\in\{-1,1\}^k}\max_jL_t(v\setminus j)
\ge 2^k e^{m g_t(1)}.
\tag{4}
```

It follows pointwise for arbitrary fibre bases, not merely on their
average. Correlated choices of the bases cannot lower the product of
these certified row norms below the product of (4).

## 3. Optimization of the obstruction constant

Let `p=k/m`. A target full-matrix quadratic cap `c N^(3/2)` in the
convention `Q=|x^TAx|/2` requires normalized defect
`gamma=1-2c sqrt(p)`. The exact original first-moment right-hand side is
at least

```math
2\exp\{m^2[p\log2+g_t(1)+t(1-2c\sqrt p)]\}.
\tag{5}
```

For `2c sqrt(p)<1`, differentiate in `t`; the minimizing Gaussian
correlation is `rho=2c sqrt(p)`. The minimum exponent is exactly

```math
\inf_{t>0}\{p\log2+g_t(1)+t(1-2c\sqrt p)\}
=p\log2+\frac14\log(1-4c^2p).
\tag{6}
```

At `c_0=sqrt(15)/8`, the right side is a concave function of `p` which
equals zero at both endpoints `p=0,1`, hence is nonnegative on `[0,1]`.
It only increases when `c` decreases. Thus at every `c<=c_0`, every
`0<p<=1`, every `t>0`, and every finite order, (5) is at least two.
The full-row bound cannot prove its failure probability less than one.
This remains true if `p`, `t`, recursion depth, terminal bases, or their
joint law depend on the order.

Deleting the signing's final diagonal changes normalized `Q` by
`O(N^-1/2)`, so the same obstruction applies to an asymptotic target
strictly below `c_0`. It does not prohibit a different counting
inequality, an exact contraction keeping correlations lost by Finner,
profile-specific tilts, or a different signing construction.

## 4. Exact finite regression checks

Run

```sh
.venv/bin/python computations/transfer_reconstruction_permanent_floor_exact_checks_2026_09_06.py
```

At `t=log 2` and integral coordinates the folded kernel entries are
rational. The independent checker computes 80 such permanents exactly,
checks (2) using outward rational logarithm/Gaussian intervals, and
checks all 360 coordinate deletions in (3) exactly. Forty positive
rational bistochastic matrices independently verify the scaling identity
and permanent factor in (1). It also checks 127 interior retention
points in (6) with outward rational intervals; endpoint equalities are
symbolic. All passed. These finite checks protect arithmetic and
normalization; the all-parameter barrier follows from (1)--(6).
