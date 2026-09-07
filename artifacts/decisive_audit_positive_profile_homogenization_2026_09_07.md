# Positive profile homogenization with an actual source-enumerator upper bound

Date: 2026-09-07. Status: finite theorem proved below. This is a positive
global transport bound, not a local Hilbert-norm relaxation. It removes
heterogeneous profile allocation from an UPPER pressure bound. It does
not yet control the surviving source entropy using the original seed cap.

## 1. Finite theorem, with no bound on the number of profiles

Let `K(c,d)>0` be symmetric on a finite alphabet of size `D`. At every
vertex `i` of a complete graph of order `n>=2`, independently choose a
uniform arrangement of `n-1` outgoing colors with exact counts `k_ic`.
Set `nu_i(c)=k_ic/(n-1)`, `nu_bar=n^(-1)sum_i nu_i`, and

```math
 Z(\nu_1,\ldots,\nu_n)
   =\mathbb E\prod_{i<j}K(X_{ij},X_{ji}),\qquad
 {cal T}_K(\nu)=
   \max_{\gamma=\gamma^T,\,\gamma1=\nu}
       \{\langle\gamma,\log K\rangle
                         -D(\gamma\Vert\nu\otimes\nu)\}.
```

Then the following finite, uniform bound holds:

```math
 \log Z(\nu_1,\ldots,\nu_n)
 \le {n\choose2}{\cal T}_K(\bar\nu)
       +nD\log n+\frac n2.                               (1)
```

Thus the leading upper bound is the HOMOGENEOUS self-transport value
at the averaged profile, even when all `n` profiles are distinct.
No concavity of an actual source entropy or symmetry of the profile list
is assumed.

### Proof

Let `Omega` be the set of legal directed-color arrays. Its cardinality is
`prod_i multinomial(n-1;k_i)`. For any law `P` on `Omega`, denote the
pair law on edge `ij` by `gamma_ij(c,d)`. Entropy subadditivity gives

```math
 H(P)\le\sum_{i<j}H(\gamma_{ij}).                         (2)
```

Choose a UNIFORM ordered pair of distinct vertex labels `(I,J)`,
independently of an array from `P`, and retain the endpoint colors
`C=X_IJ,D'=X_JI`. Its joint law `Gamma` on `(I,C,J,D')` has both
endpoint marginals

```math
                  \alpha(i,c)=\nu_i(c)/n.                (3)
```

This uses the exact row counts, not independence under `P`. A direct
entropy calculation gives

```math
 D(\Gamma\Vert\alpha\otimes\alpha)
 =\log\frac n{n-1}+\frac2n\sum_iH(\nu_i)
       -\frac1{n(n-1)}\sum_{i\ne j}H(\gamma_{ij}).        (4)
```

The color-only marginal `bar_gamma` is symmetric and has marginal
`nu_bar`. Relative-entropy data processing under deletion of the two
vertex labels gives

```math
 D(\Gamma\Vert\alpha\otimes\alpha)
       \ge D(\bar\gamma\Vert\bar\nu\otimes\bar\nu).     (5)
```

Combine (2)--(5) and average the edge energy. With `e_n=binom(n,2)`,

```math
 H(P)+\mathbb E_P\sum_{i<j}\log K(X_{ij},X_{ji})
 \le e_n\big[{cal T}_K(\bar\nu)+\log(n/(n-1))\big]
                        +(n-1)\sum_iH(\nu_i).            (6)
```

The finite Gibbs variational principle maximizes the left side to the
logarithm of the UNNORMALIZED legal-array partition function. Subtract
`sum_i log multinomial(n-1;k_i)`. The elementary type-class bound gives

```math
 0\le(n-1)H(\nu_i)-\log\binom{n-1}{k_i}\le D\log n.    (7)
```

Finally `e_n log[n/(n-1)]<=n/2`. This proves (1). Zero-count colors
are harmless throughout: their endpoint mass and all associated terms
are zero.

## 2. Exact inequality for the director's heterogeneous limit

For a finite class distribution `(pi_a,nu_a)`, let `F` be the exact
limiting pressure in
`decisive_director_positive_transport_limit_2026_09_06.md`. Then

```math
                 F\le\frac12{\cal T}_K(\bar\nu),
                 \qquad\bar\nu=\sum_a\pi_a\nu_a.         (8)
```

There is also a direct limiting proof. For an admissible family
`gamma_ab`, put

```math
 \alpha(a,c)=\pi_a\nu_a(c),\qquad
 \Gamma(a,c,b,d)=\pi_a\pi_b\gamma_{ab}(c,d).
```

The NEIGHBOR-AVERAGED constraints in that theorem make both marginals
of `Gamma` equal to `alpha`. Its full objective, including the row
multinomial normalization, is exactly

```math
 \frac12\{\mathbb E_\Gamma\log K(C,D)
                              -D(\Gamma\Vert\alpha^2)\}.
```

Discarding the class labels by data processing gives (8). In particular,
allowing classes to allocate colors differently toward different future
neighbors does not invalidate the homogeneous UPPER bound. Equality is
not asserted; the class-pair marginal constraint can make the inequality
strict.

## 3. Concrete spin-source enumerator consequence

Use the actual source-enumerator interface in
`decisive_bridge_positive_homogeneous_pressure_2026_09_06.md`. For each
legal profile `nu` with denominator `n-1`, let `N_n(nu)>=0` be its
weighted number of actual group spin tuples (including any stipulated
independent basis randomness and omitted-self-block probability).
Assume the same enumerator is used at all groups. The exact positive
spin-summed partition under consideration is

```math
 Y_n=\sum_{\nu_1,\ldots,\nu_n}
            \left[\prod_iN_n(\nu_i)\right]
                         Z(\nu_1,\ldots,\nu_n).           (9)
```

Define `s_n(nu)=n^(-1)log N_n(nu)` on profiles with positive weight,
and let `s_n^conc` be its least concave majorant on their convex hull.
Equivalently it is the supremum of `sum_a pi_a s_n(nu_a)` over finite
mixtures with the indicated average profile. Then

```math
 \frac1{n^2}\log Y_n
 \le\sup_\nu\left\{s_n^{\rm conc}(\nu)
                      +\frac{n-1}{2n}{\cal T}_K(\nu)\right\}
              +\frac{2D\log n+1/2}{n}.                  (10)
```

Indeed there are at most `n^D` legal row profiles, and thus at most
`exp(Dn log n)` profile tuples. Bound the sum by that number times its
largest summand, apply (1), and use the defining concave-majorant
inequality to the `n` profiles in that summand. This proves (10) with
the displayed finite error. Fractional enumerators cause no change.

This is an ACTUAL upper bound for the selected positive defect partition
once (9) has been established for the literal weave. The seed reflection
`R_A` remains in the positive kernel `K`; there is no edgewise unitary
norm erasure. If a known majorant for `s_n` is already concave, it can
replace `s_n^conc` directly in (10). Merely knowing the total number of
spins only yields the seed-blind constant entropy majorant.

The unresolved original-problem bridge is now sharply located: control
this source-entropy majorant and self-transport expression from the
ACTUAL seed cap, with vanishing normalized loss, or derive a different
all-order comparison for actual minimizing signings. Neither (8) nor
(10) supplies that missing estimate or proves convergence of `M_n/n^(3/2)`.
