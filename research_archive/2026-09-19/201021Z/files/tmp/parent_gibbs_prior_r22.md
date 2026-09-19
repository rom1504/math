# Wave 22 memo: parent Gibbs priors reduce overlap to a conditional free-energy gap

## Status

All channel identities, KL formulas, pressure reductions, and the finite
$`A_9`$ counts below are **Verified**.  The checker is
`tmp/parent_gibbs_prior_r22.py`.  The asymptotic large-deviation and
conditional-overlap statements singled out below are **Open targets**.
Scalar parent pressure can close the rate half of (10.694) only through a
precise near-ceiling large-deviation lemma; ordinary pressure bounds do not
supply it.  Conditional outside free energy gives a strictly more flexible
route which scalar partition functions do not determine.  Nothing here proves
convergence.

Use the ordered pairing and oriented projective state space from (10.643).
For $`n,m\ge3`$, $`|\mathcal D_n|=2^n`$, and every fixed child cut on an
$`m`$-set has exactly $`2^{n-m}`$ full extensions.  All logarithms and KL
divergences are in nats.

## A common prior converts average KL into mutual information

Fix an order-$`n`$ signing $`A`$ and define the parent Gibbs prior

```math
\nu_\beta(d)
=\frac{e^{\beta\langle A,d\rangle}}{Z_A(\beta)},
\qquad
Z_A(\beta)=\sum_{d\in\mathcal D_n}e^{\beta\langle A,d\rangle}.
\tag{R22.1}
```

For any selector law $`\pi`$ and any family of channels $`P_S`$ on full
cuts, let $`P_D=\mathbb E_{S\sim\pi}P_S`$.  The common-prior chain rule is

```math
\boxed{
\mathbb E_{S\sim\pi}D_{\rm KL}(P_S\Vert\nu_\beta)
=I_{\rm Sh}(S;D)+D_{\rm KL}(P_D\Vert\nu_\beta).
}
\tag{R22.2}
```

Consequently

```math
I_{\rm Sh}(S;D)
\le\mathbb E_SD_{\rm KL}(P_S\Vert\nu_\beta).
\tag{R22.3}
```

This is the advantage of one parent prior: no separate output-entropy bound
is needed.  The remaining question is whether many selectors can be made
close to the same $`\nu_\beta`$.

## Hard near-ground conditioning and conditional outside sums

For an $`m`$-set $`S`$ and a child cut $`y\in\mathcal D_S`$, put

```math
c_S(y)=c_A(S,y),
\qquad
\Delta_S(y)=Q(A[S])-c_S(y).
```

The conditional outside partition sum is

```math
\boxed{
K_{\beta,S}(y)
=\sum_{d:d[S]=y}
e^{\beta[\langle A,d\rangle-c_S(y)]}.
}
\tag{R22.4}
```

It contains all cross and outside-block energy, with the child contribution
removed.  The exact parent disintegration is

```math
\boxed{
Z_A(\beta)
=\sum_{y\in\mathcal D_S}e^{\beta c_S(y)}K_{\beta,S}(y),
\qquad
\nu_\beta(d\mid d[S]=y)
=\frac{e^{\beta[\langle A,d\rangle-c_S(y)]}}
       {K_{\beta,S}(y)}.
}
\tag{R22.5}
```

For tolerance $`t\ge0`$, let

```math
N_{S,t}=\{y:\Delta_S(y)\le t\},
\qquad
G_{S,t}=\{d:d[S]\in N_{S,t}\}.
```

Choose the channel which conditions the parent prior on this event:

```math
P_{S,t}^{\beta}=\nu_\beta(\,\cdot\mid G_{S,t}).
```

It has distortion at most $`t`$, and conditioning is the least-KL law
supported on the event.  Precisely,

```math
\boxed{
\begin{aligned}
D_{\rm KL}(P_{S,t}^{\beta}\Vert\nu_\beta)
&=-\log\nu_\beta(G_{S,t})\\
&=\log Z_A(\beta)
-\log\!\sum_{y\in N_{S,t}}
e^{\beta c_S(y)}K_{\beta,S}(y).
\end{aligned}
}
\tag{R22.6}
```

Define the parent--child overlap fraction

```math
\boxed{
\Omega_{\beta,S}(t)
=\frac{
\sum_{y\in N_{S,t}}e^{\beta c_S(y)}K_{\beta,S}(y)
}{Z_A(\beta)}.
}
\tag{R22.7}
```

Then (R22.2) immediately gives the rate--distortion certificate

```math
\boxed{
R_\pi(t)
\le\mathbb E_{S\sim\pi}[-\log\Omega_{\beta,S}(t)].
}
\tag{R22.8}
```

This is already the exact overlap quantity needed by (10.694), rather than
a scalar pressure surrogate.

There is also a soft optimal form.  For fixed $`S`$ and target mean
distortion $`\delta`$, the least KL relative to $`\nu_\beta`$ is

```math
\boxed{
\inf_{P:\,\mathbb E_P\Delta_S\le\delta}
D_{\rm KL}(P\Vert\nu_\beta)
=\sup_{\lambda\ge0}
\left\{
-\lambda\delta
-\log\mathbb E_{\nu_\beta}e^{-\lambda\Delta_S(D)}
\right\}.
}
\tag{R22.9}
```

The optimizer is the exponential tilt
$`dP/d\nu_\beta\propto e^{-\lambda\Delta_S}`$.  Thus hard conditioning is
not hiding a better generic variational channel; (R22.9) is the exact
parent-prior rate--distortion curve.

## Child Gibbs law with the optimal parent outside extension

Allow the child and parent temperatures to differ.  Put

```math
\mu_{\gamma,S}(y)
=\frac{e^{\gamma c_S(y)}}{Z_S(\gamma)},
\qquad
Z_S(\gamma)=\sum_ye^{\gamma c_S(y)}.
```

Sample $`Y`$ from $`\mu_{\gamma,S}`$ and, conditional on $`Y=y`$, extend it
using the parent conditional law in (R22.5).  Call the resulting full channel
$`P_S^{\beta,\gamma}`$.  It has the ordinary child-Gibbs distortion

```math
\delta_{\gamma,S}
=Q(A[S])-\mathbb E_{\mu_{\gamma,S}}c_S(Y)
\le\frac{m\log2}{\gamma}.
\tag{R22.10}
```

Because channel and prior use the same conditional outside law, their KL is
only the KL of their child marginals.  Direct substitution gives

```math
\boxed{
\begin{aligned}
D_{\rm KL}(P_S^{\beta,\gamma}\Vert\nu_\beta)
={}&\log Z_A(\beta)-\log Z_S(\gamma)\\
&+(\gamma-\beta)\mathbb E_{\mu_{\gamma,S}}c_S(Y)
-\mathbb E_{\mu_{\gamma,S}}\log K_{\beta,S}(Y).
\end{aligned}
}
\tag{R22.11}
```

The structure is clearest after defining

```math
L_{\beta,\gamma,S}(y)
=e^{(\beta-\gamma)c_S(y)}K_{\beta,S}(y).
```

Since

```math
\frac{Z_A(\beta)}{Z_S(\gamma)}
=\mathbb E_{\mu_{\gamma,S}}L_{\beta,\gamma,S}(Y),
```

equation (R22.11) becomes the exact arithmetic--geometric mean gap

```math
\boxed{
D_{\rm KL}(P_S^{\beta,\gamma}\Vert\nu_\beta)
=\log\mathbb E_{\mu_{\gamma,S}}L_{\beta,\gamma,S}
-\mathbb E_{\mu_{\gamma,S}}\log L_{\beta,\gamma,S}.
}
\tag{R22.12}
```

At matched temperatures this is

```math
\boxed{
D_{\rm KL}(P_S^{\beta,\beta}\Vert\nu_\beta)
=\log\mathbb E_{\mu_{\beta,S}}K_{\beta,S}
-\mathbb E_{\mu_{\beta,S}}\log K_{\beta,S}.
}
\tag{R22.13}
```

This identifies exactly what scalar pressures omit.  The parent and child
partition functions determine the arithmetic mean
$`\mathbb E_{\mu}L=Z_A/Z_S`$, whereas the KL also needs the geometric mean.
For a fixed probability law and fixed arithmetic mean, a constant profile
$`L`$ gives gap zero, while putting $`L=\varepsilon`$ on a positive-mass set
and compensating elsewhere makes the gap grow like
$`\log(1/\varepsilon)`$.  This abstract profile argument is not a signing
counterexample.  It proves that scalar parent and child partition values
alone do not determine, or upper-bound at the required scale, the common-prior
rate.  A theorem must control the conditional profile $`K_{\beta,S}(y)`$.

For distortion $`t=O(n^{3/2-c})`$, (R22.10) permits

```math
\gamma\ge\frac{m\log2}{t}=O(n^{-1/2+c}).
\tag{R22.14}
```

The exact Gibbs continuation is therefore a conditional-free-energy flatness
lemma at this temperature scale:

```math
\mathbb E_{S\sim\pi}
\left[
\log\mathbb E_{\mu_{\gamma,S}}L_{\beta,\gamma,S}
-\mathbb E_{\mu_{\gamma,S}}\log L_{\beta,\gamma,S}
\right]
=O(n^{1/2-2c}).
\tag{R22.15}
```

This would give the required rate in (10.694).  It is a log-mgf bound at
parameter one for the centered conditional outside free energy, not ordinary
scalar pressure convergence.

## Temperature mixtures

A finite or polynomial grid of parent temperatures can hedge selector
dependence.  Let

```math
\overline\nu=\sum_{\ell=1}^Lw_\ell\nu_{\beta_\ell}.
```

If selector $`S`$ uses a channel constructed relative to
$`\nu_{\beta_{\ell(S)}}`$, then pointwise
$`\overline\nu\ge w_{\ell(S)}\nu_{\beta_{\ell(S)}}`$ and hence

```math
\boxed{
D_{\rm KL}(P_S\Vert\overline\nu)
\le D_{\rm KL}(P_S\Vert\nu_{\beta_{\ell(S)}})
-\log w_{\ell(S)}.
}
\tag{R22.16}
```

Uniform weights cost only $`\log L`$.  A polynomial grid therefore has
$`O(\log n)`$ overhead, smaller than $`n^{1/2-2c}`$ for $`c<1/4`$.
Temperature mixing permits an $`S`$-dependent optimizer, but it does not
remove the overlap or Jensen-gap requirement.

## What parent pressure alone would have to prove

There is one rigorous pressure-only certificate.  Fix $`y`$ on $`S`$ and
average its external energy uniformly over its $`2^{n-m}`$ extensions.  Every
cross or outside edge contains a uniform relative spin, so this average is
zero.  Jensen therefore gives

```math
\boxed{K_{\beta,S}(y)\ge2^{n-m}.}
\tag{R22.17}
```

Let

```math
\phi_A(\beta)
=\log\left(2^{-n}Z_A(\beta)\right),
\qquad
\phi_A^*(x)
=\sup_{\beta\ge0}\{\beta x-\phi_A(\beta)\}.
```

Using one exact child ground in (R22.6) gives, for every $`S`$,

```math
\boxed{
-\log\Omega_{\beta,S}(0)
\le m\log2+\phi_A(\beta)-\beta Q(A[S]).
}
\tag{R22.18}
```

Since $`Q(A[S])\ge q_m`$, optimizing one common temperature yields

```math
\boxed{
R_\pi(0)
\le m\log2-\phi_A^*(q_m)
\quad\text{for every selector law }\pi.
}
\tag{R22.19}
```

When the supremum is attained in the interior, its temperature is
characterized by

```math
\phi_A'(\beta_*)=q_m.
```

The bound has an exact ceiling.  Equation (R22.17) also implies

```math
\phi_A(\beta)\ge\beta Q(A[S])-m\log2,
```

so $`\phi_A^*(Q(A[S]))\le m\log2`$.  Thus (R22.19) is small only when the
parent energy rate function nearly saturates the full cost of fixing the
selected $`m`$-vertex oriented cut.

For a fixed-ratio power-saving comparison, parent pressure **alone** would
suffice through the following precise statement:

> **Parent-energy sharp-LDP lemma.**  For some $`0<c<1/4`$, the chosen
> target-specific order-$`n`$ minimizer and every required
> $`m=\rho n+O(1)`$ obey
>
> ```math
> \boxed{
> \phi_A^*(q_m)
> \ge m\log2-O(n^{1/2-2c}).
> }
> \tag{R22.20}
> ```

With $`\pi=U_m`$, (R22.19) gives exact-ground distortion zero,
$`R_\pi(0)=O(n^{1/2-2c})`$, and
$`D_{\rm KL}(\pi\Vert U_m)=0`$.  Hence (10.694) supplies the desired
power-saving restriction edge.  This is a concrete sufficient lemma about
the full parent energy large-deviation function, not a request that a scalar
pressure limit merely exist.

Near-ground multiplicity and outside alignment can be strictly stronger.
From (R22.17),

```math
\boxed{
-\log\Omega_{\beta,S}(t)
\le m\log2+\phi_A(\beta)
-\beta[Q(A[S])-t]-\log|N_{S,t}|.
}
\tag{R22.21}
```

The exact overlap (R22.7) can improve further when
$`K_{\beta,S}(y)`$ is large precisely on near-ground child states.  Neither
improvement is encoded by $`\phi_A`$ alone.

## The exact asymptotic overlap lemma

The least scalar-free statement which feeds the proved fixed-slice mgf is:

> **Parent-Gibbs overlap lemma.**  There exist $`0<c<1/4`$, fixed-ratio
> comparison windows, a selector law $`\pi`$, a tolerance
> $`t=O(n^{3/2-c})`$, and either one temperature or a polynomial temperature
> mixture such that
>
> ```math
> \boxed{
> \mathbb E_{S\sim\pi}[-\log\Omega_{\beta,S}(t)]
> +D_{\rm KL}(\pi\Vert U_m)
> =O(n^{1/2-2c}).
> }
> \tag{R22.22}
> ```

By (R22.8), this is exactly the missing hypothesis of (10.694), and therefore
proves the power-saving principal restriction comparison.  It can hold even
if the pressure-only sharp-LDP lemma (R22.20) fails, because it uses all
near-ground child states and their conditional outside partition weights.

The clean falsification criterion is equally exact: for some fixed ratio,
prove that every target-specific exact minimizer, every selector law with
$`D_{\rm KL}(\pi\Vert U_m)=O(n^{1/2-2c})`$, every permitted tolerance, and
every polynomial parent-temperature mixture has average negative log overlap
$`\omega(n^{1/2-2c})`$.

## Exact $`A_9`$ audit

Use the hidden-optimal deletion law

```math
\pi_*={1\over25}(4,2,4,0,4,4,2,5,0).
```

There are 512 full oriented projective cuts, with parent energy levels
$`-24,-16,-8,0,8,16,24`$ and multiplicities

```text
25, 60, 111, 120, 111, 60, 25.
```

For tolerances $`t=0,4,8,12`$, exhaustive enumeration shows that, for every
deleted coordinate, the fraction of the $`t`$-event within an energy level is
nondecreasing with the energy.  Hence
$`\nu_\beta(G_{S,t})`$ is nondecreasing in $`\beta`$ by the elementary
covariance inequality for two increasing functions, and $`\beta=\infty`$ is
the exact optimizer in this one-parameter parent family.

Among the 25 positive parent grounds, the event counts by coordinate are

| tolerance | counts on coordinates `0,...,8` | $`\mathbb E_{\pi_*}[-\log\Omega_{\infty,S}(t)]`$ |
|---:|:---|---:|
| `0` | `(4,8,4,8,4,4,8,8,8)` | `1.5830484787467...` |
| `4` | `(12,16,12,18,12,12,16,16,18)` | `0.6304036289976...` |
| `8` | `(20,22,20,24,20,20,22,22,24)` | `0.1888318865847...` |
| `12` | `(24,24,24,25,24,24,24,24,25)` | `0.0408219945203...` |

These values are exact finite sums

```math
-\frac1{25}\sum_i(25\pi_i)\log\frac{k_i(t)}{25}.
```

At tolerance zero, the optimal unrestricted shared prior from (10.676) costs
only `1.022686788870...` nats, so the parent-ground prior is useful but not
optimal.  At $`\beta=1`$, the common-prior chain rule is audited exactly
numerically as

```math
1.583853764561459
=1.330682170797721+0.253171593763739,
```

namely average conditional KL equals mutual information plus output-marginal
KL.  The finite wall therefore supports the common-prior mechanism while also
showing that the conditional overlap profile, not parent energy alone, gains
the last compression.

## Disposition

The parent Gibbs construction makes the missing theorem substantially more
precise:

1. a near-ceiling parent energy LDP (R22.20) is a purely scalar sufficient
   lemma and would already close a power-saving comparison;
2. the weaker and more natural target is the parent--child overlap estimate
   (R22.22), which includes near-ground multiplicity and conditional outside
   alignment;
3. for child Gibbs channels, the identical content is conditional outside
   free-energy flatness (R22.15), an arithmetic--geometric mean gap invisible
   to separate parent and child pressures.

Temperature mixtures cost only logarithmically and do not change these
requirements.  Further work should attack either the sharp parent energy LDP
or the selector-averaged conditional profile $`K_{\beta,S}`$; another scalar
pressure telescope without one of these quantitative statements cannot meet
$`R_\pi(\delta)+D_{\rm KL}(\pi\Vert U_m)=O(n^{1/2-2c})`$.

