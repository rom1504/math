### Wave 22, Route 2: parent-Gibbs common priors and conditional outside free energy

**Status.** The identities and the $A_9$ audit below are verified by
`tmp/parent_gibbs_prior_r22.py`. The asymptotic overlap estimates are open.
All logarithms are natural.

For an order-$n$ signing $A$, let

```math
\nu_\beta(d)=\frac{e^{\beta\langle A,d\rangle}}{Z_A(\beta)}.
```

For any selector law $\pi$ and channels $P_S$ on full cuts, with mixture
$P_D=\mathbb E_{S\sim\pi}P_S$, the exact common-prior identity is

```math
\boxed{\mathbb E_S D_{\rm KL}(P_S\Vert\nu_\beta)
=I_{\rm Sh}(S;D)+D_{\rm KL}(P_D\Vert\nu_\beta).}
```

Thus the average KL to one parent prior upper-bounds the selector rate.

For a child cut $y$ on $S$, write $c_S(y)=c_A(S,y)$,
$\Delta_S(y)=Q(A[S])-c_S(y)$, and define the conditional outside sum

```math
K_{\beta,S}(y)=\sum_{d:d[S]=y}
 e^{\beta(\langle A,d\rangle-c_S(y))}.
```

Then

```math
Z_A(\beta)=\sum_y e^{\beta c_S(y)}K_{\beta,S}(y).
```

For $N_{S,t}=\{y:\Delta_S(y)\le t\}$, conditioning the parent Gibbs law on
$d[S]\in N_{S,t}$ gives distortion at most $t$ and the exact cost

```math
\boxed{
D_{\rm KL}(P_{S,t}^{\beta}\Vert\nu_\beta)
=-\log\Omega_{\beta,S}(t),\qquad
\Omega_{\beta,S}(t)
=\frac{\sum_{y\in N_{S,t}}e^{\beta c_S(y)}K_{\beta,S}(y)}{Z_A(\beta)}.}
```

Consequently

```math
\boxed{R_\pi(t)\le\mathbb E_{S\sim\pi}[-\log\Omega_{\beta,S}(t)].}
```

The exact missing overlap lemma is therefore: for some $0<c<1/4$, fixed
ratio $m/n$, permitted $t=O(n^{3/2-c})$, and a selector law $\pi$,

```math
\mathbb E_S[-\log\Omega_{\beta,S}(t)]
+D_{\rm KL}(\pi\Vert U_m)=O(n^{1/2-2c}).
```

A polynomial mixture of temperatures is allowed at only $O(\log n)$ KL
overhead. This statement feeds (10.694) directly.

There is an equivalent soft Gibbs formulation. Let
$\mu_{\gamma,S}(y)=e^{\gamma c_S(y)}/Z_S(\gamma)$ and extend $y$ using the
parent conditional outside law. Its distortion satisfies

```math
Q(A[S])-\mathbb E_{\mu_{\gamma,S}}c_S
\le \frac{m\log2}{\gamma}.
```

With

```math
L_{\beta,\gamma,S}(y)
=e^{(\beta-\gamma)c_S(y)}K_{\beta,S}(y),
```

the exact common-prior KL is the Jensen gap

```math
\boxed{
D_{\rm KL}(P_S^{\beta,\gamma}\Vert\nu_\beta)
=\log\mathbb E_{\mu_{\gamma,S}}L_{\beta,\gamma,S}
-\mathbb E_{\mu_{\gamma,S}}\log L_{\beta,\gamma,S}.}
```

At $\beta=\gamma$ this is
$\log\mathbb E_\mu K_{\beta,S}-\mathbb E_\mu\log K_{\beta,S}$. Parent and
child scalar partition functions determine the arithmetic mean of $L$, but
not its geometric mean. Hence scalar pressures alone do not control this KL;
one needs conditional-free-energy flatness or direct overlap control.

There is nevertheless one sharp pressure-only sufficient condition. Uniform
averaging over the $2^{n-m}$ extensions of a child cut gives

```math
K_{\beta,S}(y)\ge2^{n-m}.
```

Set
$\phi_A(\beta)=\log(2^{-n}Z_A(\beta))$ and
$\phi_A^*(x)=\sup_{\beta\ge0}(\beta x-\phi_A(\beta))$. One exact child
ground then yields

```math
\boxed{R_\pi(0)\le m\log2-\phi_A^*(q_m).}
```

Thus parent pressure alone closes the rate requirement if

```math
\boxed{\phi_A^*(q_m)\ge m\log2-O(n^{1/2-2c}).}
```

This is a near-saturation large-deviation theorem, much stronger than the
existence or regularity of a pressure limit. Near-ground multiplicity and
alignment of $K_{\beta,S}$ can make the full overlap route succeed even when
this scalar criterion fails.

**$A_9$ audit.** For the hidden-optimal deletion law
$25\pi_*=(4,2,4,0,4,4,2,5,0)$, the parent has 25 positive ground cuts. For
tolerances $t=0,4,8,12$, the event fraction is nondecreasing across every
parent energy level for every deletion, so $\beta=\infty$ exactly optimizes
this one-temperature family. The parent-ground event counts and weighted
average costs are

| $t$ | counts by deleted coordinate | average $-\log\Omega$ |
|---:|:---|---:|
| 0 | `(4,8,4,8,4,4,8,8,8)` | `1.5830484787467` |
| 4 | `(12,16,12,18,12,12,16,16,18)` | `0.6304036289976` |
| 8 | `(20,22,20,24,20,20,22,22,24)` | `0.1888318865847` |
| 12 | `(24,24,24,25,24,24,24,24,25)` | `0.0408219945203` |

At $t=0$ the unrestricted optimal shared prior costs only
`1.022686788870...` nats, so the parent-ground prior is useful but not
optimal. At $\beta=1$ the chain rule numerically audits as

```math
1.583853764561459
=1.330682170797721+0.253171593763739.
```

**Conclusion.** The strongest route exposed here is the exact
selector-averaged parent--child overlap (or equivalently the Jensen-gap
bound). The scalar parent-pressure route is valid only if one proves the
specific near-ceiling estimate above; a generic pressure telescope cannot
replace that lemma.
