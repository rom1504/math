# Adversarial audit: orientation-uniform cluster promotion from row Renyi two

Object audited:
[`../drafts/actual_child_orientation_uniform_cluster_promotion.md`](../drafts/actual_child_orientation_uniform_cluster_promotion.md)

Verdict: **PASS, with a scope qualification.**  Central symmetry plus a
dimension-free row `D_2` bound really does imply the dimension-free linear
MGF estimate used by the random-row-cut/Hanson--Wright proof.  The constants
in OU.1--OU.10 are consistent with SP.2.  Thus the conditional cluster
promotion can be applied in either orientation, and a target-reaching
orientation can be analyzed directly.  This removes balanced-orientation
target relevance for that conditional promotion route; it does not prove
that the balanced orientation itself reaches the target or decide the
remaining product phase.

## 1. The Renyi-two subgaussian lemma

Put `f=dP/dU_n`, `K=exp(C/2)`, and `a=||v||_2^2`.  Central symmetry gives

```math
 E_Pe^{\langle v,R\rangle}=E_P\cosh\langle v,R\rangle.
 \tag{A.OU.1}
```

Because `E_Uf=1` and `E_Uf^2<=e^C`,

```math
 E_P\cosh\langle v,R\rangle
 \le1+K\{E_U(\cosh\langle v,R\rangle-1)^2\}^{1/2}. \tag{A.OU.2}
```

The subtraction of one is essential: applying Cauchy--Schwarz directly to
`cosh` would leave a fixed multiplicative factor at `v=0` and would not be
a centered subgaussian estimate.

For `X=<v,R>`,

```math
 \begin{aligned}
 E_U(\cosh X-1)^2
 &=\tfrac12E_U\cosh(2X)-2E_U\cosh X+\tfrac32\\
 &\le\tfrac12e^{2a}-2(1+a/2)+\tfrac32\\
 &={e^{2a}-1-2a\over2}
 \le a^2e^{2a}.
 \end{aligned}                                      \tag{A.OU.3}
```

Here `E cosh(2X)=prod_j cosh(2v_j)<=e^(2a)`, while
`E cosh X>=1+E X^2/2=1+a/2`.  The last inequality is the ordinary
exponential Taylor remainder.  Hence

```math
 E_Pe^{\langle v,R\rangle}\le1+Kae^a.
 \tag{A.OU.4}
```

For `a>=0`, `ae^a<=e^(2a)-1`; for `K>=1`, convexity of `x mapsto x^K`
gives `1+K(e^(2a)-1)<=e^(2Ka)`.  Therefore

```math
 \boxed{E_Pe^{\langle v,R\rangle}
 \le e^{2e^{C/2}\|v\|_2^2}.}                       \tag{A.OU.5}
```

In the convention `Ee^(<v,R>)<=exp(sigma^2||v||^2/2)`, this is exactly
`sigma^2=4e^(C/2)`.  No coordinate independence under `P` was used.
Central symmetry is indispensable because otherwise the linear term at the
origin need not vanish.

## 2. Application to the actual canonical row

For either orientation and either transpose direction, the erased-row
forward likelihood satisfies `z(-b)=z(b)` by the global spin flip of the
base child.  Its inverse escort is therefore centrally symmetric.  Flipping
one row bit changes `log z` by at most `2u`; the conditional escort lemma
then gives

```math
 D_2(r_{\epsilon,u}\Vert U_n)
 \le n\log\{1+\tanh^2(\lambda u)\}
 \le\lambda^2u^2n.                                  \tag{A.OU.6}
```

At `u=beta/sqrt(N)`, one may use `C=lambda^2 beta^2`, uniformly in the row
width.  Thus

```math
 \sigma_*^2=4e^{\lambda^2\beta^2/2}.                \tag{A.OU.7}
```

The random-row-cut proof in SP.2 requires only independence between row
blocks and (A.OU.5).  With a general proxy `sigma^2`, its determinant bound
gives

```math
 |\theta|\|M\|_{op}\le{1\over2\sqrt2\sigma^2},
 \qquad
 \log Ee^{\theta H}\le4\sigma^4\theta^2V.          \tag{A.OU.8}
```

Substitution of (A.OU.7) yields exactly

```math
 a_*={1\over8\sqrt2e^{\lambda^2\beta^2/2}},
 \qquad b_*=64e^{\lambda^2\beta^2}.                 \tag{A.OU.9}
```

There is no hidden dependence on the row dimension.

## 3. Physical scaling and orientation scope

For the sector--Gram quadratic chaos, `V=K_epsilon` and
`||M||op<=sqrt(2K_epsilon)`.  At
`theta=-lambda beta^2/N`, the determinant premise becomes

```math
 \lambda\beta^2\sqrt{2\kappa}\le a_*
 \quad\text{when}\quad K_\epsilon\le\kappa N^2,    \tag{A.OU.10}
```

which is OU.9.  The quadratic contribution is
`b_*lambda^2 beta^4 K_epsilon/N^2=O(1)`.  Combining it with the exact
orientationwise cluster remainder from SP.1 proves OU.10--OU.11 with no
normalization loss.

If the joint two-orientation soft bridge reaches a target, at least one
sector does so up to the harmless one-bit `O(1)` difference: the joint soft
minimum lies between the better sector value and that value plus
`(log 2)/lambda`.  OU.2 may be applied in that sector.  Therefore a separate
theorem comparing the balanced and target-reaching orientations is not
needed **for this conditional cluster-promotion strategy**.

The scope should not be enlarged further.  OU.2 still assumes its Gram
threshold and the same all-order absolute cluster tail in the target
orientation.  It neither proves that the balanced orientation is
target-reaching nor transports a reverse-product or coherent-retuning
certificate between orientations.  The valid frontier change is removal of
the target-relevance premise from the conditional promotion branch, not a
decision of `L_balanced-product-phase` or a Level-6 recurrence.
