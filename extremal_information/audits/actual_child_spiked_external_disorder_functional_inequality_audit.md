# Functional inequalities and the rare-bridge ceiling for spiked external disorder

Status: **rigorous dimension-free functional-inequality theorem, exact
two-replica variance reduction, and decisive falsifier for the all-tilt
superconcentration target in Theorem 37.42**.  The natural spiked row law is
not the source of a bad concentration constant: it satisfies dimension-free
Poincare and logarithmic Sobolev inequalities.  At zero external tilt, the
variance of the exact optimizing-child pressure is controlled by one
explicit edge-cavity two-replica overlap.

However, the proposed subgaussian estimate with proxy `o(N)` for *all* MGF
parameters is false.  Every actual child pair admits a rank-one bridge atom
of probability `exp(-Theta(N^2))` whose pressure is `Theta(N^(3/2))`, whereas
the mean pressure under the spiked disorder is only `O(N)`.  Any all-parameter
subgaussian proxy is therefore `Omega(N)`.  Thus a successful proof must
concentrate the channel output itself, truncate the external pressure, or use
a response-specific tilted overlap; ordinary external-disorder
superconcentration cannot close the spiked route.

Throughout, the children are the actual contracted-temperature pressure
minimizers.  No conference or surrogate child is used.

## 1. The row law has dimension-free Poincare and log-Sobolev constants

Gauge `y` to the all-one vector and write

```math
 z(b)={1\over\sqrt n}\sum_{j=1}^n b_j,
 \qquad e(b)={1+z(b)^2\over2},
 \qquad d\mu=e\,dU_n.                              \tag{FI.1}
```

For a cube function put

```math
 D_jf(b)={f(b)-f(b^{(j)})\over2},
 \qquad {\cal E}_\rho(f)=\sum_jE_\rho(D_jf)^2.    \tag{FI.2}
```

The following elementary multiplier estimate is useful in its own right.

**Lemma FI.1 (first-chaos multiplier).**  For every `g` on the fair cube,

```math
 \boxed{\|zg\|_2^2\le
 4\sum_j\|D_jg\|_2^2+2\|g\|_2^2.}                \tag{FI.3}
```

*Proof.*  Decompose `g=sum_(k=0)^n g_k` into homogeneous Walsh levels.
Multiplication by `z` is the normalized sum of the up and down incidence
operators on the Boolean lattice.  On level `k`, their squared operator
norms are respectively

```math
 {(k+1)(n-k)\over n}\le k+1,
 \qquad {k(n-k+1)\over n}\le k.
```

At each output level, use `\|u+v\|^2\le2\|u\|^2+2\|v\|^2`, then sum over
levels.  This gives

```math
 \|zg\|_2^2\le
 2\sum_k(2k+1)\|g_k\|_2^2
 =4\langle g,\Delta g\rangle+2\|g\|_2^2,
```

which is (FI.3). `square`

**Theorem FI.2 (dimension-free row inequalities).**  For every `n` and
every `y`,

```math
 \boxed{\operatorname {Var}_{\mu_y}f
 \le7\,{\cal E}_{\mu_y}(f),}                    \tag{FI.4}
```

and

```math
 \boxed{\operatorname {Ent}_{\mu_y}(f^2)
 \le41\,{\cal E}_{\mu_y}(f).}                  \tag{FI.5}
```

The same constants hold for `mu_y^(otimes m)`, with the Dirichlet form
summed over all row bits.

*Proof of (FI.4).*  In the variational definition of variance choose the
constant `E_Uf` and put `g=f-E_Uf`.  Uniform Poincare and (FI.3) give

```math
 \begin{aligned}
 \operatorname {Var}_\mu f
 &\le {1\over2}\|g\|_2^2+{1\over2}\|zg\|_2^2\\
 &\le {3\over2}\|g\|_2^2+2{\cal E}_U(g)
 \le {7\over2}{\cal E}_U(f).
 \end{aligned}                                    \tag{FI.6}
```

Since `e>=1/2`, `E_U(D_jf)^2<=2E_mu(D_jf)^2`, proving (FI.4).

*Proof of (FI.5).*  The case `n=1` is uniform, so assume `n>=2`.  A bit flip
changes `z` by `2/sqrt(n)`, while
`|(d/ds)log(1+s^2)|<=1`.  Hence

```math
 |\log e(b)-\log e(b^{(j)})|\le {2\over\sqrt n},
 \qquad {e(b^{(j)})\over e(b)}\le r_*:=e^{\sqrt2}. \tag{FI.7}
```

Let `A_jf=(f+f^(j))/2`.  The product rule

```math
 D_j(uv)=A_ju\,D_jv+D_ju\,A_jv
```

and (FI.7) imply

```math
 (A_j\sqrt e)^2\le A_je,
 \qquad
 |D_j\sqrt e|\le {1\over2\sqrt n}A_j\sqrt e.     \tag{FI.8}
```

Also, by flip invariance of `U_n`,

```math
 E_U[(A_je)(A_jf^2)]
 \le {1+r_*\over2}E_\mu f^2.                     \tag{FI.9}
```

Equations (FI.8)--(FI.9) give

```math
 \sum_jE_U[D_j(\sqrt e f)]^2
 \le2{\cal E}_\mu(f)+{1+r_*\over4}E_\mu f^2.    \tag{FI.10}
```

Since

```math
 \operatorname {Ent}_\mu(f^2)
 =\operatorname {Ent}_U(ef^2)-E_\mu[f^2\log e]
 \le\operatorname {Ent}_U((\sqrt e f)^2)
      +(\log2)E_\mu f^2,                          \tag{FI.11}
```

the fair-cube log-Sobolev inequality yields the defective estimate

```math
 \operatorname {Ent}_\mu(f^2)
 \le4{\cal E}_\mu(f)
 +\left({1+r_*\over2}+\log2\right)E_\mu f^2.     \tag{FI.12}
```

Apply (FI.12) to `f-E_muf`, then use Rothaus' lemma and (FI.4).  The
resulting constant is

```math
 4+7\left(2+{1+e^{\sqrt2}\over2}+\log2\right)<41.
```

This proves (FI.5).  Poincare and log-Sobolev tensorize over independent
rows. `square`

Thus the external row law has no growing functional-inequality cost.  Any
power saving must come from the pressure observable.

## 2. Exact pressure gradient and the two-replica proxy

Let `A,D` be the actual children of orders `m,n`, let `N=m+n`, and put
`t=beta/sqrt(N)`.  Write

```math
 L(B)=\log E_{x,z}\cosh\left(t\{H_A(x)+\epsilon H_D(z)
                                  +x^{\mathsf T}Bz\}\right).    \tag{FI.13}
```

Introduce the usual auxiliary spin `tau`, so the Gibbs weight is
proportional to

```math
 \exp\left(t\tau\{H_A(x)+\epsilon H_D(z)+x^{\mathsf T}Bz\}\right).
```

For a bridge edge `a=(i,j)`, delete that edge and let

```math
 r_a(B_{-a})=E_{\nu^0_{B,a}}[\tau x_i z_j].        \tag{FI.14}
```

**Theorem FI.3 (exact cavity-gradient variance reduction).**  If `B^a`
flips edge `a`, then

```math
 \boxed{
 D_aL(B)=B_a\operatorname {arctanh}
              \{\tanh(t)r_a(B_{-a})\}.}           \tag{FI.15}
```

Consequently

```math
 \boxed{
 \operatorname {Var}_{\mu_y^{\otimes m}}L
 \le7\sum_aE\operatorname {arctanh}^2
                   \{\tanh(t)r_a\}
 \le7t^2\sum_aE r_a^2.}                          \tag{FI.16}
```

*Proof.*  With edge `a` deleted, insertion of the sign `s` multiplies the
partition function by

```math
 \cosh t+s r_a\sinh t.
```

Taking the half log-ratio at `s=B_a` and `s=-B_a` proves (FI.15).  Since
`arctanh(r tanh t)<=r t` for `0<=r<=1`, (FI.16) follows from (FI.4).
`square`

There is an exact full-Gibbs two-replica interpretation.  Put

```math
 m_a(B)=E_{\nu_B}[\tau x_i z_j],
 \qquad c=\tanh t.
```

Insertion of edge `a` gives

```math
 m_a={r_a+B_ac\over1+B_acr_a},
 \qquad
 |r_a-m_a|\le {c\over1-c}.                        \tag{FI.17}
```

For normalized replica overlaps

```math
 R_X={1\over m}\langle X^1,X^2\rangle,
 \qquad R_Z={1\over n}\langle Z^1,Z^2\rangle,
```

one has, pointwise in `B`,

```math
 {1\over mn}\sum_a m_a(B)^2
 =E_{\nu_B^{\otimes2}}[\tau^1\tau^2R_XR_Z]\ge0. \tag{FI.18}
```

Thus (FI.16)--(FI.18) imply

```math
 \boxed{
 \operatorname {Var}_{\mu_y^{\otimes m}}L
 \le14t^2mn\,
 E_{B,\nu_B^{\otimes2}}[\tau^1\tau^2R_XR_Z]
 +{14t^2mn\,c^2\over(1-c)^2}.}                  \tag{FI.19}
```

At a balanced split the last term is `O_beta(1)`.  Moreover the normalized
cavity square mass in (FI.16) is `o(1)` if and only if the averaged overlap
in (FI.18) is `o(1)`: their difference is `O(t)` by (FI.17).  Therefore
the needed power saving in the Poincare carre-du-champ is equivalent to the
vanishing overlap target

```math
 \boxed{
 E_{B\sim\mu_y^{\otimes m},\nu_B^{\otimes2}}
       [\tau^1\tau^2R_XR_Z]=o(1).}                \tag{FI.20}
```

This is a scalar two-replica observable of the actual externally coupled
children.  It is strictly less information than the bridge pressure table.
No currently proved minimizing-child identity implies (FI.20): flip,
contraction, deletion, and sector--Gram identities concern the zero-bridge
child law, while (FI.20) concerns a positive-density external bridge
ensemble.

## 3. What log-Sobolev does and does not add

Let

```math
 \Gamma_L(B)=\sum_a(D_aL(B))^2.                   \tag{FI.21}
```

The flip ratio in (FI.7), (FI.5), and the elementary inequality

```math
 \left({e^{u/2}-e^{v/2}\over2}\right)^2
 \le{(u-v)^2\over32}(e^u+e^v)
```

give, for every real `s`,

```math
 \operatorname {Ent}(e^{sL})
 \le K_*s^2E[e^{sL}\Gamma_L],
 \qquad
 K_*={41(1+e^{\sqrt2})\over8}.                    \tag{FI.22}
```

Hence an all-tilt bound

```math
 \sup_{s\in\mathbb R}E_s\Gamma_L\le G_N,
 \qquad {dP_s\over d\mu_y^{\otimes m}}
 ={e^{sL}\over Ee^{sL}},                          \tag{FI.23}
```

would imply

```math
 \log E e^{s(L-EL)}\le K_*G_Ns^2.                \tag{FI.24}
```

By (FI.15), (FI.23) is an all-external-tilt version of the two-replica
condition (FI.20).  The next theorem shows that it cannot have the required
power saving for all `s`.

## 4. Rare rank-one bridges force a linear subgaussian proxy

**Theorem FI.4 (all-tilt superconcentration is impossible).**  Fix
`beta>0` and a balanced window
`theta N<=m,n<=(1-theta)N`.  Let `A,D` be the actual
contracted-temperature minimizing children and let `L` be (FI.13).  For
every spike direction `y`, if

```math
 \log E_{\mu_y^{\otimes m}}
   e^{s(L-EL)}\le {s^2\sigma_N^2\over2}
 \qquad(s\in\mathbb R),                           \tag{FI.25}
```

then, for all sufficiently large `N`,

```math
 \boxed{\sigma_N^2\ge c_{\beta,\theta}N.}         \tag{FI.26}
```

In particular the `sigma_N^2=o(N)` sufficient hypothesis in Theorem 37.42
is false even on the actual optimizing-child law.

*Proof.*  First the mean pressure is only linear.  Since

```math
 \|e_y\|_2^2={3\over2}-{1\over2n},
```

Cauchy--Schwarz gives, uniformly in `z in {+-1}^n`,

```math
 E_{\mu_y}e^{tB\cdot z}
 \le\sqrt{3/2}\,(\cosh2t)^{n/2}
 \le\sqrt{3/2}\,e^{nt^2}.                        \tag{FI.27}
```

After averaging the `m` rows and using
`cosh(a+b)<=2cosh(a)cosh(b)`,

```math
 \log E_{\mu_y^{\otimes m}}e^{L(B)}
 \le {m\over2}\log(3/2)+mnt^2+\log2+p_A(t)+p_D(t). \tag{FI.28}
```

Actual minimality and averaging over all internal edge signings give

```math
 p_A(t)\le {m\choose2}\log\cosh t,
 \qquad p_D(t)\le {n\choose2}\log\cosh t.        \tag{FI.29}
```

Thus Jensen and `t=beta/sqrt(N)` yield

```math
 E_{\mu_y^{\otimes m}}L\le C_{\beta}N.           \tag{FI.30}
```

Choose `x_0,z_0` with `H_A(x_0)>=0` and
`epsilon H_D(z_0)>=0`; such configurations exist because each quadratic
form has fair-cube mean zero.  Set

```math
 B_*=x_0z_0^{\mathsf T}.
```

Then

```math
 L(B_*)\ge t mn-(N+1)\log2.                       \tag{FI.31}
```

Every row atom under `mu_y` has mass at least `2^(-n-1)`, because
`e_y>=1/2`.  Consequently

```math
 \mu_y^{\otimes m}(B_*)\ge2^{-(mn+m)}.            \tag{FI.32}
```

For balanced splits, (FI.30)--(FI.31) show that
`d_N=L(B_*)-EL>=c_(beta,theta)N^(3/2)` for large `N`.  The subgaussian
Chernoff bound from (FI.25), applied to the event containing the single atom
`B_*`, gives

```math
 2^{-(mn+m)}
 \le\exp\left(-{d_N^2\over2\sigma_N^2}\right).
```

Rearrangement proves (FI.26). `square`

The obstruction also occurs directly in the tilted carre-du-champ, rather
than only through the tail consequence of a hypothetical MGF bound.

**Corollary FI.5 (the all-tilt gradient proxy is extensive).**  In the
setting of FI.4, let `B_max` maximize `L` over sign bridges.  Then

```math
 \boxed{
 \Gamma_L(B_{\max})\ge c_{\beta,\theta}N,
 \qquad
 \liminf_{s\to+\infty}E_s\Gamma_L
 \ge c_{\beta,\theta}N.}                          \tag{FI.33}
```

*Proof.*  For a fixed sign bridge `B`, extend the pressure radially:

```math
 F_B(u)=\log E_{x,z}\cosh
 \left(t\{H_A(x)+\epsilon H_D(z)+u x^{\mathsf T}Bz\}\right),
 \qquad 0\le u\le1.                              \tag{FI.34}
```

This is convex in `u`.  Equations (FI.29),
`cosh(a+b)<=2cosh(a)cosh(b)`, and the rank-one lower bound (FI.31) imply

```math
 F_{B_{\max}}(1)-F_{B_{\max}}(0)
 \ge c_{\beta,\theta}N^{3/2}.                    \tag{FI.35}
```

Consequently

```math
 t\sum_a(B_{\max})_a m_a(B_{\max})
 =F_{B_{\max}}'(1)
 \ge c_{\beta,\theta}N^{3/2}.                    \tag{FI.36}
```

Cauchy--Schwarz yields `sum_a m_a(B_max)^2>=cN^2`.  By (FI.17),
`sum_a r_a(B_max)^2>=cN^2` for large `N`.  Finally
`arctanh(c|r|)>=c|r|` and `c=tanh t` show from (FI.15) that
`Gamma_L(B_max)>=cN`.  As `s` tends to positive infinity, the finite tilted
law in (FI.23) concentrates on the maximizers of `L`.  The same argument
applies to every maximizer, proving the second assertion. `square`

The rare atom in FI.4 has probability `exp(-Theta(N^2))`; its contribution
to any fixed low-order or zero-tilt moment can vanish.  Thus FI.4 does not
rule out (FI.20), nor does it prove a linear range for the channel response
`R(v)=E[L\mid V=v]`.  It rules out using the complete external pressure MGF
as the route to that response range.  Conditional expectation through the
strict channel can erase these ultra-rare atoms even though it cannot erase
an actual coherent switching phase.

## 5. Scope and revised missing statement

The strongest unconditional concentration statement is now (FI.19): a
dimension-free Poincare theorem reduces the actual-law variance to the
zero-tilt two-replica product overlap (FI.20).  The analogous all-tilt
statement needed to derive Theorem 37.42 by ordinary Herbst is false at the
required scale by FI.4.

Accordingly the surviving spiked-response target must be response-specific:
prove directly that the conditional channel output has `o(N)` range, or
prove a truncated/moderate-tilt overlap estimate that controls its
`2^m` query extrema without charging the `exp(-Theta(N^2))` rank-one bridge
tail.  Merely asking for all-parameter external-disorder
superconcentration is now a rigorously falsified route.
