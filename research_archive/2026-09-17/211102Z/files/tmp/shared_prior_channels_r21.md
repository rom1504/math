# Wave 21 memo: universal noisy-ground channels have linear capacity

## Status

The channel formulas, entropy bounds, parameter optimization, and scoped
coefficient no-go below are **proved exactly** in nats.  The accompanying
checker is `tmp/shared_prior_channels_r21.py`.  The no-go concerns the
universal capacity estimates furnished by these explicit channels; it does
not prove that their actual mutual information is linear for quadratic
minimizers.  A minimizer-specific shared-prior overlap theorem could make the
actual rate much smaller.  Nothing here proves convergence.

Throughout, $`n\ge m\ge3`$, so an oriented projective cut
$`d=\sigma xx^{\mathsf T}`$ (diagonal ignored) has one orientation bit and
$`n-1`$ projective spin bits.  Thus

```math
|\mathcal D_n|=2^n,
\qquad \log|\mathcal D_n|=n\log2.
```

All logarithms and entropies are natural.  Write

```math
h(p)=-p\log p-(1-p)\log(1-p)
```

for binary entropy in nats.

## Exact noisy-ground channel

Fix a selector law $`\pi`$.  For each $`m`$-set $`S`$, choose by a
deterministic tie rule an oriented projective ground $`d_S^*`$ of $`A[S]`$,
so

```math
c_A(S,d_S^*)=Q(A[S]).
```

Construct a full output cut as follows.

1. Independently flip every one of the $`m`$ child spins with probability
   $`p\in[0,1/2]`$.  Only the projective class of this noise matters.
2. Flip the orientation with probability $`q\in[0,1/2]`$.
3. Conditional on the resulting noisy cut on $`S`$, choose each of the
   $`n-m`$ outside relative spins uniformly and independently.

The last step is a uniform choice among the exactly $`2^{n-m}`$ full
extensions of an oriented projective child cut.  It changes neither the
restriction nor its payoff.

Put

```math
\theta=1-2p,
\qquad
\eta=1-2q,
\qquad
a=\eta\theta^2.
```

For distinct $`i,j\in S`$, independence gives
$`\mathbb E(Z_iZ_j)=\theta^2`$, while the orientation multiplier has mean
$`\eta`$.  Therefore, conditionally on every $`S`$,

```math
\boxed{
\mathbb E[c_A(S,D)\mid S]=aQ(A[S]).
}
\tag{R21.1}
```

If

```math
\overline Q_\pi=\mathbb E_{S\sim\pi}Q(A[S]),
```

then the channel distortion is exactly

```math
\boxed{
\delta(p,q)=(1-a)\overline Q_\pi.
}
\tag{R21.2}
```

There is no factor-of-two ambiguity: both $`Q`$ and $`c_A`$ use the ordered
matrix pairing.  Multiplying an oriented cut by the orientation noise changes
the sign of the entire ordered payoff.

## Projective-noise entropy and capacity

Let $`Z=(Z_1,\ldots,Z_m)`$ have independent flip probability $`p`$, and let
$`[Z]=\{Z,-Z\}`$ be its projective class.  Its exact entropy is

```math
\boxed{
H_m^{\rm proj}(p)
=-\frac12\sum_{k=0}^m\binom mk z_k\log z_k,
\quad
z_k=p^k(1-p)^{m-k}+p^{m-k}(1-p)^k.
}
\tag{R21.3}
```

The factor $`1/2`$ removes the double count of a class by $`Z`$ and $`-Z`$.
Also

```math
\boxed{
(m-1)h(p)\le H_m^{\rm proj}(p)\le mh(p).
}
\tag{R21.4}
```

The upper bound follows because $`[Z]`$ is a function of $`Z`$.  Conversely,
$`([Z],Z_1)`$ determines $`Z`$, so

```math
H_m^{\rm proj}(p)
=mh(p)-H(Z_1\mid[Z])
\ge(m-1)h(p).
```

Given $`S`$, the orientation noise, projective child noise, and outside
relative spins are independent coordinates of the full output cut.  Hence

```math
\boxed{
H(D\mid S)
=(n-m)\log2+H_m^{\rm proj}(p)+h(q).
}
\tag{R21.5}
```

Since $`H(D)\le n\log2`$,

```math
\boxed{
I_{\rm Sh}(S;D)
\le C_m(p,q)
:=m\log2-H_m^{\rm proj}(p)-h(q).
}
\tag{R21.6}
```

Thus (R21.2) and (R21.6) give the universal rate--distortion point

```math
R_\pi((1-a)\overline Q_\pi)\le C_m(p,q).
\tag{R21.7}
```

At $`p=q=0`$ this is the exact-ground channel with random outside extension:

```math
R_\pi(0)\le m\log2.
\tag{R21.8}
```

A deterministic outside extension has only
$`H(D\mid S)=H(D[S]\mid S)`$ and gives the weaker alphabet bound
$`n\log2-H(D[S]\mid S)`$.  Uniform outside extension improves the universal
certificate by exactly $`(n-m)\log2`$ without changing distortion.

The capacity estimate can be sharp in an abstract channel model.  If the
selected ground labels are uniform on $`\mathcal D_m`$, additive
orientation/projective noise makes the output uniform and equality holds in
(R21.6).  This model is not asserted to come from a quadratic minimizer.  It
shows that alphabet size and independent noise alone cannot improve the
bound; overlap among minimizers' child ground labels must do the work.

## Rigorous optimization over spin and orientation noise

For a required retained coefficient $`a\in(0,1)`$, write

```math
f(a)
=\log2-h\!\left(\frac{1-\sqrt a}{2}\right)
=\frac12\left[(1+\sqrt a)\log(1+\sqrt a)
 +(1-\sqrt a)\log(1-\sqrt a)\right].
\tag{R21.9}
```

The best finite-$`m`$ capacity certificate in the spin/orientation family is
the explicit one-dimensional optimization

```math
\boxed{
C_m^*(a)
=\inf_{\sqrt a\le\theta\le1}
\left\{
m\log2
-H_m^{\rm proj}\!\left(\frac{1-\theta}{2}\right)
-h\!\left(\frac{1-a/\theta^2}{2}\right)
\right\}.
}
\tag{R21.10}
```

Indeed $`\theta^2\ge a`$ and the orientation mean must be
$`\eta=a/\theta^2`$.  Equation (R21.4) gives the sharp leading optimization

```math
\boxed{
mf(a)-\log2
\le C_m^*(a)
\le mf(a)+\log2.
}
\tag{R21.11}
```

For the lower bound, every admissible $`\theta`$ has
$`h((1-\theta)/2)\le h((1-\sqrt a)/2)`$, and orientation noise saves at most
one bit.  For the upper bound take $`\theta=\sqrt a`$ and $`q=0`$.
Thus orientation randomization cannot improve the leading rate: independent
spin noise with

```math
p_a=\frac{1-\sqrt a}{2}
```

is asymptotically optimal in this family.

Erasure does not improve the leading curve.  Time-share a channel retaining
$`a_1\ge a`$ with a uniform full cut, using the former with probability
$`a/a_1`$.  The retained payoff is $`a`$, and mutual-information convexity
gives leading rate

```math
m\frac a{a_1}f(a_1).
```

The absolutely convergent series

```math
\boxed{
f(a)=\sum_{k\ge1}\frac{a^k}{2k(2k-1)}
}
\tag{R21.12}
```

shows that $`f(a)/a`$ is increasing.  Therefore the last display is minimized
at $`a_1=a`$.  Pure exact-ground erasure has rate $`am\log2`$, and
$`f(a)\le a\log2`$ shows directly that spin noise is at least as good.

## Child Gibbs channels

There is a separate exact channel which does not select one ground first.
On the $`2^m`$ oriented projective child cuts, define

```math
\mu_{\beta,S}(d)
=\frac{e^{\beta c_A(S,d)}}{Z_{\beta,S}},
\qquad
Z_{\beta,S}=\sum_{d\in\mathcal D_m}e^{\beta c_A(S,d)}.
```

Sample $`d`$ from this law and extend it by uniform outside relative spins.
Put

```math
\delta_{\beta,S}
=Q(A[S])-\mathbb E_{\mu_{\beta,S}}c_A(S,d),
\qquad
\mathcal Z^{\rm def}_{\beta,S}
=\sum_de^{-\beta[Q(A[S])-c_A(S,d)]}.
```

The Gibbs entropy identity is

```math
\boxed{
H(\mu_{\beta,S})
=\beta\delta_{\beta,S}
 +\log\mathcal Z^{\rm def}_{\beta,S}.
}
\tag{R21.13}
```

In particular,

```math
0\le\delta_{\beta,S}
\le\frac{H(\mu_{\beta,S})}{\beta}
\le\frac{m\log2}{\beta}.
\tag{R21.14}
```

Writing $`\delta_\beta=\mathbb E_\pi\delta_{\beta,S}`$ and
$`\overline H_\beta=\mathbb E_\pi H(\mu_{\beta,S})`$, uniform outside
extension gives the exact conditional entropy and universal rate certificate

```math
\boxed{
H(D\mid S)=(n-m)\log2+\overline H_\beta,
\qquad
I_{\rm Sh}(S;D)\le m\log2-\overline H_\beta.
}
\tag{R21.15}
```

Equivalently, by (R21.13),

```math
I_{\rm Sh}(S;D)
\le m\log2-\beta\delta_\beta
-\mathbb E_\pi\log\mathcal Z^{\rm def}_{\beta,S}.
\tag{R21.16}
```

At fixed $`\beta>0`$, (R21.14) makes the distortion $`O(n)`$, already
subleading relative to a fixed-ratio $`n^{3/2}`$ budget.  But no universal
lower bound presently forces $`\overline H_\beta`$ close enough to
$`m\log2`$.  The rate certificate can remain linear.  Deterministic outside
extension again loses $`(n-m)\log2`$ in the alphabet bound.

This is not a scalar-pressure restatement.  The exact unresolved quantity is
the entropy of the **mixture output**:

```math
\boxed{
I_{\rm Sh}(S;D)
=H(D)-(n-m)\log2-\overline H_\beta.
}
\tag{R21.17}
```

Separate child partition functions determine $`\overline H_\beta`$ but not
$`H(D)`$.  A useful theorem must prove overlap among the different child
Gibbs laws, or exhibit a common prior $`\nu`$ with small average divergence
$`\mathbb E_SD_{\rm KL}(P_{D\mid S}\Vert\nu)`$.

## Fixed-ratio coefficient audit

Let $`m=\rho n+O(1)`$ for fixed $`\rho\in[1/2,1)`$, and set

```math
b_\rho=\rho^{3/2}-\rho^2
=\rho^{3/2}(1-\sqrt\rho).
```

The leading allowance in (10.670)--(10.671) is

```math
B_{n,m}=b_\rho q_n+O(n^{3/2-c}).
```

The scalar coefficient is small:

```math
\boxed{
0<b_\rho\le\frac{27}{256}=0.10546875,
}
\tag{R21.18}
```

with equality at $`\rho=9/16`$.

There are two complementary coefficient estimates.  The first describes the
most immediate universal certificate; the second is a genuine necessary
condition for this channel.

- The crude bound $`\overline Q_\pi\le q_n`$ makes
  $`\delta\le(1-a)q_n`$; making this universal upper certificate fit even
  when the whole budget is reserved for distortion requires
  $`a\ge1-b_\rho`$.
- More robustly, $`\overline Q_\pi\ge q_m`$.  With the proved asymptotic
  bounds

  ```math
  q_m\ge(\gamma-o(1))m^{3/2},
  \qquad q_n\le(1+o(1))n^{3/2},
  \qquad \gamma=0.672986728862\ldots,
  ```

  the necessary distortion condition $`(1-a)\overline Q_\pi\le B_{n,m}`$
  gives

  ```math
  \boxed{
  a\ge a_\rho-o(1),
  \qquad
  a_\rho=1-\frac{1-\sqrt\rho}{\gamma}>0
  \quad(\rho\ge1/2).
  }
  \tag{R21.19}
  ```

Consequently (R21.11) gives

```math
C_m^*(a)\ge \rho f(a_\rho)n-O(1)=\Omega_\rho(n).
\tag{R21.20}
```

For example, at $`\rho=1/2`$, (R21.19) gives
$`a_\rho=0.5647860407\ldots`$ and
$`f(a_\rho)=0.3178609794\ldots`$ nats per selected vertex.  The cruder
$`a\ge1-b_\rho`$ gives $`f(a)=0.5704456478\ldots`$.

There is a fully proved universal fixed-cut mgf bound.  Since principal
monotonicity gives

```math
|c_A(S,d)|\le Q(A[S])\le q_n,
```

Hoeffding's lemma yields, for every selector law,

```math
\boxed{
\psi_\pi(\lambda)\le\frac{q_n^2\lambda^2}{2}.
}
\tag{R21.21}
```

Inserting the noisy-ground capacity certificate into (10.670) therefore
produces

```math
\delta+q_n\sqrt{2C_m^*(a)}.
\tag{R21.22}
```

Under (R21.19), its information term is $`\Omega(n^2)`$, while the budget is
$`O(n^{3/2})`$.  No choice of spin noise, orientation noise, or erasure makes
this **universal capacity/Hoeffding certificate** meet the fixed-ratio
coefficient budget.  The same conclusion holds in (10.671), even with the
best choice $`\pi=U_m`$ making $`D_{\rm KL}(\pi\Vert U_m)=0`$; a nonuniform
selector only introduces an additional term which the channel does not
control.

Even a hypothetical genuine proxy $`v_n=O(n^{5/2})`$ in (10.672) would give

```math
\sqrt{2v_nC_m^*(a)}=\Omega(n^{7/4}),
```

still too large.  With a much stronger $`v_n=O(n^2)`$ proxy, linear rate is
at the correct leading scale and constants would matter, but no such mgf
bound is proved.  A robust power-saving conclusion at that scale would be
obtained from actual rate $`o(n)`$.

This is deliberately a scoped no-go.  Equations (R21.6) and (R21.15) are
capacity upper bounds, not lower bounds on the actual rate.  Quadratic
minimizers could have far more shared structure than the abstract sharpness
model.  At minimum, an $`O(n^{5/2})`$ log-mgf proxy requires actual rate
$`O(\sqrt n)`$ with a constant fitting the budget left after distortion; a
genuine $`O(n^2)`$ proxy requires $`O(n)`$ with the analogous sharp constant.
The clean power-saving target at the first scale is

```math
\boxed{
I_{\rm Sh}(S;D)=o(\sqrt n)
}
```

and it is $`I_{\rm Sh}(S;D)=o(n)`$ at the second scale.  By (R21.17) and its
noisy-ground analogue, this must come from a minimizer-specific shared output
prior or overlap theorem, not from additional independent noise or separate
scalar child pressures.
