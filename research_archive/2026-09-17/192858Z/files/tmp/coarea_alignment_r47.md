# Wave 47B: active-ground alignment and clipped-deficit coarea

All algebraic identities and the `A6/A8/A9` rational values below are
**Verified** by `tmp/coarea_alignment_r47_check.py`.  The random-signing
comparison is **Numerical**.  No asymptotic minimizer-specific coarea theorem
is proved.

## 1. Exact adjacent-port decomposition

Let two adjacent selectors be `U+a` and `U+b`.  For a port `v`, an orientation
`tau in {+-1}`, and a fixed center `z`, write

```math
h_v(z_U)=\sum_{i\in U}a_{vi}z_i,
\qquad
Q_v^\tau=\max_{y_U}
\left\{\tau y_U^{\mathsf T}A[U]y_U+2|h_v(y_U)|\right\},
\qquad Q_v=\max_\tau Q_v^\tau.
```

Define three nonnegative regrets

```math
\begin{aligned}
o_v^\tau&=Q_v-Q_v^\tau &&\text{(orientation gap)},\\
g_v^\tau(z_U)&=Q_v^\tau-
 \{\tau z_U^{\mathsf T}A[U]z_U+2|h_v(z_U)|\}
 &&\text{(common-core response gap)},\\
r_v^\tau(z)&=2\{|h_v(z_U)|-\tau z_vh_v(z_U)\}
 &&\text{(port-spin regret)}.
\end{aligned}
```

Direct expansion gives the **Verified exact identity**

```math
\boxed{
Q(A[U+v])-\tau z_{U+v}^{\mathsf T}A[U+v]z_{U+v}
=o_v^\tau+g_v^\tau(z_U)+r_v^\tau(z),
}
```

and hence

```math
\boxed{
D_z(U+v)=\min_{\tau=\pm1}
\{o_v^\tau+g_v^\tau(z_U)+r_v^\tau(z)\}.}
\tag{R47B.1}
```

This is the cancellation which is invisible if cap motion and center-energy
motion are bounded separately: the entire common-core energy cancels before
the three nonnegative port regrets appear.  In particular, for the minimizing
orientation at an active selector, each of the three terms is at most the
deficit.  The checker verifies (R47B.1) for every center, adjacent port, and
orientation in the three stored instances.

The identity does **not** by itself prove coarea.  On an edge from an active
port to an inactive port it controls the three regrets at the active end but
does not control growth of the orientation or optimized common-core response
gap at the other end.  The exact remaining alignment statement can be phrased
as a one-sided clipped-gradient estimate for the minimum in (R47B.1); generic
separate `L^2` bounds do not imply it.

## 2. Exact finite coarea ratios

For all three stored matrices the deficits are in `4 Z`, while the actual
caps are respectively

```text
A6: 0.940591076...,  A8: 2.276095342...,  A9: 2.462452602....
```

Thus at the actual cap, and throughout `0<H<4`,
`u_z(S)=H 1{D_z(S)=0}`.  Every coarea ratio is exactly the hard escape ratio

```math
\frac{\frac12\mathbb E|u(S)-u(T)|}
{(1-\lambda_1)\mathbb E u}
=\frac{\epsilon_\ell}{1-\lambda_1}.
\tag{R47B.2}
```

For uniform centers the useful near-core values are

| instance | core | normalized ratio |
|:--|:--|:--|
| `A8,m=6` | `ell=4` | `239/250` |
| `A8,m=6` | `ell=5` | `111/100` |
| `A9,m=7` | `ell=4` | `8413/9072` |
| `A9,m=7` | `ell=5` | `13609/13608 = 1+1/13608` |
| `A9,m=7` | `ell=6` | `433/378` |

Therefore `A8` at one replacement, and `A9` at both one and two
replacements, are exact-minimizer finite obstructions to a universal
uniform-center version of (10.1198).  The `A9,ell=5` failure is tiny but
strict.  Increasing the cap past the next deficit lattice point improves all
of these finite ratios; the checker audits every interval using caps
`1,5,9,...`.  This does not help the actual finite cap and gives no
asymptotic permission to enlarge the target tolerance.

## 3. Two canonical ground-coupled laws and the exact missing correlation

Let

```math
a_z=U_m\{S:D_z(S)=0\},\qquad
B_z=\langle 1_{D_z=0},(I-K_\ell)1_{D_z=0}\rangle.
```

The **ground-lift law** is `nu_lift(z)=a_z/E a_z`: equivalently sample
uniformly from all zero-deficit incidence pairs `(z,S)` and forget `S`.
It is canonical and selector-independent; it is not selection of the best
center after seeing the test selector.  Under this law the exact coarea
criterion is

```math
\boxed{
\frac{B_{\nu_{\rm lift}}}{J_{\nu_{\rm lift}}}
=\frac{\mathbb E[a_zB_z]}{\mathbb E[a_z^2]},
\qquad
\mathbb E[a_zB_z]
\le(1-\lambda_1)(1-\eta)\mathbb E[a_z^2].}
\tag{R47B.3}
```

Thus the genuinely missing fact is degree--escape alignment: after two
ground incidences size-bias a center, its Johnson escape must be smaller than
the spectral gap.  This is a concrete correlation inequality, not the old
separate cap/energy `L^2` estimate.

At `0<H<4`, the exact ground-lift ratios are

```text
A6 ell=2,3,4: 7/24;
A8 ell=3,4,5: 1401/1760, 1321/1540, 1209/1232;
A9 ell=4,5,6: 50383/59346, 80845/89019, 10154/9891.
```

Hence size bias repairs the `A8` one-replacement and `A9` two-replacement
failures, but `A9` one-replacement still fails.

Uniform full-parent absolute grounds give a second selector-independent law.
Its corresponding ratios are

```text
A6 ell=2,3,4: 1/6;
A8 ell=3,4,5: 73/110, 467/660, 35/44;
A9 ell=4,5,6: 4321/5184, 6919/7776, 433/432.
```

It too fails strictly on `A9` at one replacement.

## 4. Admissibility and falsification judgment

The two structured laws are selector-independent, but neither is known to be
supported on the required project-row class.  A full parent ground has local
margins `ell_i>=0`, `sum_i ell_i=q_n`, and therefore only

```math
R_2(z)=\sum_i\ell_i^2\le(n-1)q_n=O(n^{5/2}),
```

which misses `O(n^{9/4-c})`.  A random outside extension of a child ground has
the same `O(n^{5/2})` generic scale from (10.804).  In the three finite audits
every center happens to satisfy the old `2n(n-1)` cap, so these are admissible
finite diagnostics only.  Row-truncating either law preserves selector
independence, but then both the denominator and (R47B.3) are open.

The improvement from size bias is not by itself minimizer-specific.  The
identity (R47B.3) holds for every family, and the fixed-coordinate cylinder
from Wave 46 defeats it when all degrees/escapes are equal.  In an additional
300-sample random-signing diagnostic, ground-lift ratios were usually *better*
than the stored minimizers (for `n=8,m=6,ell=5`, 89.7% of random samples were
below one; for `n=9,m=7,ell=5`, all were below one).  This is numerical only,
but it rejects interpreting the finite improvement as evidence of special
minimizer geometry.

The sharp surviving theorem is therefore (R47B.3) **after project-row
truncation**, or its positive-cap analogue expressed through the port regrets
in (R47B.1), together with project-scale nonzero mass.  Exact minimizer
minimality has not supplied this correlation.  The finite data sharply
falsify a law-free or nearest-neighbor theorem, while leaving a two-or-more
replacement, row-qualified, ground-incidence alignment theorem open.

