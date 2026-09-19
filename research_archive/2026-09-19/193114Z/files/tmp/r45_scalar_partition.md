# Wave 45A: the abundance partition and central block Harnack are too strong

## Status

**Verified obstruction.**  The row-weighted incidence partition proposed in
Idea 4 always has

```math
Z_\lambda\le \exp\{-\Omega(n^{3/4})\}
\qquad(\lambda\ge0)
\tag{R45.1}
```

on exact minimizers.  It therefore cannot have the project-saving lower bound
`exp{-O(n^(3/4-c))}` for any fixed `c>0`.  This abundance statement is
strictly stronger than the existential scalar target and is retired.

**Verified central-orbit obstruction.**  For `k=floor(n/2)`, uniformly in the
base cut `d`,

```math
\mathbb E_{|U|=k}\alpha_{d^U}
\le \exp\{-\Omega(n^{3/4})\},
\qquad
\mathbb E_{|U|=k}h_{d^U}
\ge\Omega(n^{3/4}).
\tag{R45.2}
```

Consequently a cut already having project surprise
`h_d=O(n^(3/4-c))` cannot also satisfy the project geometric-mean Harnack
premise in (10.1159) on the central block orbit.  The proposed half-block
implementation is falsified.

**Open.**  Blocks `k=kappa n` sufficiently far from `n/2` retain a first
chaos depending on the flipped-signing row, and high starting slack can defeat
contraction when `(1-2kappa)^2>=1/3`.  No theorem or scalable obstruction is
obtained for those surviving regimes.  The existential scalar target
(10.1149) is not falsified.

## 1. Uniform incidence is already at the no-saving exponent

Fix an exact minimizer `A`, with `Q(A)=q=q_n`.  For each selector let

```math
B_S=-A+2P_SAP_S.
```

This is again a zero-diagonal signing.  Its norms obey

```math
\lVert B_S\rVert_F^2=n(n-1),
\qquad
\lVert B_S\rVert_{\rm op}
\le3\lVert A\rVert_{\rm op}
\le3\sqrt{2q}.
\tag{R45.3}
```

The last inequality is the exact-minimizer spectral bound already used in
the ledger.  Let `d=(sigma,x)` be uniform over the `2^n` oriented projective
cuts.  Rademacher Hanson--Wright gives, uniformly in `S`,

```math
\begin{aligned}
\Pr_d\{\langle B_S,d\rangle\ge q\}
&\le \Pr_x\{|x^{\mathsf T}B_Sx|\ge q\}\\
&\le 2\exp\left[-c\min\left\{
\frac{q^2}{n(n-1)},
\frac{q}{\lVert B_S\rVert_{\rm op}}
\right\}\right].
\end{aligned}
\tag{R45.4}
```

The verified lower bound `q=Omega(n^(3/2))` makes the two exponents
`Omega(n)` and `Omega(n^(3/4))`, respectively.  Therefore

```math
\boxed{
\mathbb E_{S,d}\mathbf1_{\{S\in\mathcal I_d\}}
=\mathbb E_d\alpha_d
\le C e^{-c n^{3/4}}.}
\tag{R45.5}
```

This is an upper bound on abundance, not on the largest column.  It is fully
compatible with a rare exceptional column of project mass.

## 2. The row-weighted partition cannot prove the scalar target

For `lambda>=0`, define exactly as in Idea 4

```math
Z_\lambda
=\mathbb E_d[\alpha_de^{-\lambda R_2(d)}].
\tag{R45.6}
```

The row factor is at most one, so (R45.5) immediately proves (R45.1).
Since

```math
n^{3/4}=\omega(n^{3/4-c}),
```

no fixed constant `C` can make
`Z_lambda>=exp{-C n^(3/4-c)}` hold for all large `n`.

The exact reason this does not hurt (10.1149) is the state-count dilution.
With

```math
F_\lambda(d)=-\log\alpha_d+\lambda R_2(d),
\qquad F_\lambda^*=\min_dF_\lambda(d),
```

one has

```math
\boxed{
2^{-n}e^{-F_\lambda^*}
\le Z_\lambda\le e^{-F_\lambda^*}.}
\tag{R45.7}
```

Thus a project scalar value only forces
`Z_lambda>=exp{-n log 2-O(n^(3/4-c))}`.  Requiring `Z_lambda` itself at
project scale asks for stretched-exponential abundance of good cuts, whereas
the leading route asks for only one.

## 3. Central macroscopic blocks have an intrinsic entropy floor

Fix any base oriented cut `d=(sigma,x)` and put `k=floor(n/2)`.  If `U` is
uniform among `k`-sets, then `x^U` is a uniform Rademacher word conditioned
to lie on one central Hamming layer about `x`.  If `V` is an iid fair vertex
flip set, then

```math
\Pr\{|V|=k\}\ge c n^{-1/2}.
```

Condition (R45.4) on this event and then average over `S`.  This proves

```math
\boxed{
\mathbb E_{|U|=k}\alpha_{d^U}
\le C\sqrt n\,e^{-c n^{3/4}}}
\tag{R45.8}
```

uniformly in the base `d`.  Jensen, with `h=+infinity` for an empty column,
gives

```math
\boxed{
\mathbb E_{|U|=k}h_{d^U}
\ge-\log\mathbb E_{|U|=k}\alpha_{d^U}
\ge c n^{3/4}-O(\log n).}
\tag{R45.9}
```

In particular, if a scalar optimizer has the separate mass bound required by
the leading route,

```math
h_{d_\lambda}=O(n^{3/4-c}),
```

then

```math
\mathbb E_{|U|=k}
[h_{d_\lambda^U}-h_{d_\lambda}]
=\Omega(n^{3/4}).
\tag{R45.10}
```

This contradicts the `G_k=O(n^(3/4-c))` premise needed in (10.1159).
Hence central macroscopic switching cannot simultaneously certify the target
mass and regularize its row.  This is a global orbit wall, independent of the
one-vertex cliff from Wave 44.

The wall also applies to every layer `|k-n/2|=O(sqrt(n))`, with a changed
constant, because each such binomial atom is `Theta(n^(-1/2))`.

## 4. What happens away from the central layer

The exact contraction is simple.  For fixed `S,d`, uniform `|U|=k` gives

```math
\boxed{
\mathbb E_{|U|=k}L_S(d^U)=\vartheta_kL_S(d),
\qquad
\vartheta_k=1-\frac{4k(n-k)}{n(n-1)}.}
\tag{R45.11}
```

This is the same two-coordinate multiplier as the row identity (10.1158).
For `k/n -> kappa`, it tends to `rho^2`, where `rho=1-2kappa`.
Since `|L_S(d)|<=Q(B_S)<=3q`, a uniform constant mean gap below `q` is
present whenever

```math
\rho^2\le\frac13-\delta.
\tag{R45.12}
```

For clarity, the corresponding fixed-layer tail can be written explicitly.
Gauge `B_S` by the base word and define

```math
R_{B_S}(d)=\lVert B_Sx\rVert_2^2,
\qquad g_S=(q-\rho^2L_S(d))_+.
```

Conditioning independent biased signs only as a proof device for the
fixed-cardinality layer, linear concentration plus Hanson--Wright gives, for
fixed `kappa in (0,1)`,

```math
\boxed{
\Pr_{|U|=k}\{L_S(d^U)\ge q\}
\le C_\kappa\sqrt n
\exp\left[-c_\kappa\min\left\{
\frac{g_S^2}{\rho^2R_{B_S}(d)+n^2},
\frac{g_S}{\lVert B_S\rVert_{\rm op}}
\right\}\right],}
\tag{R45.13}
```

with the right side interpreted as the trivial bound when `g_S=0`.
The `O(q/n)` difference between `rho^2` and the exact multiplier
`vartheta_k` is absorbed in the constants under (R45.12).

This formula explains precisely why the central proof is stronger.  At
`rho=0`, the first-chaos row term vanishes and (R45.13) has exponent
`Omega(n^(3/4))`.  Away from the center, the universal estimate

```math
R_{B_S}(d)\le n\lVert B_S\rVert_{\rm op}^2=O(n^{5/2})
```

gives only an `Omega(n^(1/2))` first exponent, which does not contradict a
project Harnack bound.  Moreover, when `rho^2>=1/3`, a high-slack value
`L_S(d)` near `3q` can put the contracted mean at or above threshold.  Thus
small-`kappa`/high-slack orbits survive this audit.

There is an exact average identity which identifies a possible next input.
Writing `p_j=(m)_j/(n)_j`, direct expansion gives

```math
\boxed{
\mathbb E_{S\sim U_m}R_{B_S}(d)
=\{1-4(p_2-p_3)\}R_2(d)
+4(p_2-p_3)n(n-1).}
\tag{R45.14}
```

So a project parent row controls the *mean* flipped-signing row.  It does not
give the exponential incidence-conditioned control required in (R45.13): a
polynomial exceptional set can dominate a stretched-exponentially small hard
column.  A surviving noncentral Harnack theorem would need precisely such a
minimizer-specific incidence-conditioned row/slack estimate.

## 5. Exact finite audit

At the project price `lambda=n^(-3/2)`, exhaustive enumeration gives:

| signing | `m` | `Z_0` | scalar optimizer `(R_2,|I_d|)` | empty central switches |
|:--|--:|:--|:--|:--|
| `A_5` | 4 | `5/16` | `(16,2)` | `0/10` |
| `A_6` | 5 | `3/16` | `(30,3)` | `20/20` |
| `A_8` | 6 | `17/448` | `(40,2)` | `16/70` |
| `A_9` | 7 | `25/1536` | `(72,2)` | `70/126` |

The `A_6` optimizer has an entirely empty central orbit.  These finite facts
are exact mechanism checks, not asymptotic evidence for the Hanson--Wright
exponent.

## 6. Frontier

The abundance partition `Z_lambda` and central half-block Harnack should not
be pursued further.  The direct complement route still asks for the
existential scalar estimate

```math
\min_d\{-\log U_m(\mathcal I_d)+n^{-3/2}R_2(d)\}
=O(n^{3/4-c}).
```

Within switching methods, the only regime not covered by this obstruction is
a noncentral macroscopic block, especially a small fixed `kappa`, coupled to
a theorem controlling high complement slack and
`R_{B_S}(d)` on the actual incidence set.  Such a theorem is open and is not
implied by the unconditional average (R45.14).

## Verification

`tmp/r45_scalar_partition_check.py` verifies the partition identities,
central-orbit row and energy means, (R45.14), and the complete finite table on
`A_5,A_6,A_8,A_9`.  It ends with

```text
PASS r45_scalar_partition_check
```
