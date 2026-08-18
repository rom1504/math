# Truncating rare bridges and a moderate-tilt overlap certificate

Status: **rigorous output-specific truncation and moderate-tilt dichotomy for
the actual optimizing-child law**.  The rank-one atoms that falsify global
external-disorder superconcentration can be removed without changing any
spiked response by even a constant amount.  If the resulting spiked response
still has linear range, then one bounded, capped external tilt has a fixed
positive edge-cavity/two-replica overlap.  No `s to infinity` limit, full
bridge table, conference child, or surrogate row law appears.

This does not prove that the moderate-tilt overlap vanishes.  It replaces the
false all-tilt target by a compact-tilt scalar observable tied directly to
the channel output.

## 1. Uniform fixed-replica annealed moments

Let `A,D` be actual contracted-temperature minimizing children of orders
`m,n`, put `N=m+n` and `t=beta/sqrt(N)`, and let

```math
 L(B)=\log E_{x,z}\cosh\left(t\{H_A(x)+\epsilon H_D(z)
                                  +x^{\mathsf T}Bz\}\right).    \tag{MT.1}
```

For a spike direction `y` and word `v in {+-1}^m`, let

```math
 P_{v,y}=\bigotimes_{i=1}^m q_{v_i,y}U_n,
 \qquad
 q_{v_i,y}(b)={1\over2}
 \left(1+{v_i\langle y,b\rangle\over\sqrt n}\right)^2.        \tag{MT.2}
```

Recall that `||q_(v_i,y)||_2^2<=5/2`.

**Theorem MT.1 (fixed-replica annealed bound).**  For every positive integer
`k`, uniformly in the actual children, orientation, `y`, and `v`,

```math
 \boxed{
 \log E_{P_{v,y}}e^{kL}
 \le {m\over2}\log(5/2)+mnk^2t^2
      +k\{p_A(t)+p_D(t)+\log2\}.}                 \tag{MT.3}
```

In particular, for fixed `beta,k`, the right side is `O_(beta,k)(N)`.

*Proof.*  Represent the `k` copies of `cosh` by independent auxiliary signs
`tau_1,...,tau_k`.  Conditional on all replica spins, bridge row `i` sees

```math
 w_i=\sum_{a=1}^k\tau_ax_i^az^a\in\mathbb Z^n,
 \qquad \|w_i\|_2^2\le nk^2.                     \tag{MT.4}
```

Cauchy--Schwarz against the fair row law gives

```math
 \begin{aligned}
 E_{q_{v_i,y}U_n}e^{tB_i\cdot w_i}
 &\le\sqrt{5/2}
      \left(E_{U_n}e^{2tB_i\cdot w_i}\right)^{1/2}\\
 &\le\sqrt{5/2}\exp(t^2\|w_i\|_2^2)
 \le\sqrt{5/2}\exp(nk^2t^2).                    \tag{MT.5}
 \end{aligned}
```

The remaining `k` replica averages factor.  For one replica,

```math
 E_{\tau,x,z}e^{t\tau(H_A(x)+\epsilon H_D(z))}
 =E_{x,z}\cosh(t(H_A(x)+\epsilon H_D(z)))
 \le2e^{p_A(t)+p_D(t)}.                           \tag{MT.6}
```

Multiplying (MT.5) over the rows and (MT.6) over replicas proves (MT.3).
Actual minimality supplies

```math
 p_A(t)\le {m\choose2}\log\cosh t,
 \qquad p_D(t)\le {n\choose2}\log\cosh t,        \tag{MT.7}
```

so (MT.3) is linear in `N` for fixed `k,beta`. `square`

Theorem MT.1 is stronger than merely observing that the rank-one atom has
tiny mass: every fixed replica moment is uniformly linear on the log scale.
It does not assert such a bound when `k` grows with `N`.

## 2. A linear cap preserves every spiked response

Set

```math
 C_\beta={1\over2}\log(5/2)+{\beta^2\over2}+\log2.              \tag{MT.8}
```

The harmless `log 2` in this definition makes the following statement
valid for every `N>=1` without an asymptotic remainder.  For `delta>0`, put

```math
 T_N=(C_\beta+\delta)N,
 \qquad L_T=L\wedge T_N,                          \tag{MT.9}
```

and define

```math
 R(v)=E_{P_{v,y}}L,
 \qquad R_T(v)=E_{P_{v,y}}L_T.                   \tag{MT.10}
```

**Corollary MT.2 (uniform response truncation).**  Uniformly in `v,y` and
the actual children,

```math
 \boxed{
 0\le R(v)-R_T(v)\le e^{-\delta N}.}              \tag{MT.11}
```

Consequently

```math
 \left|\operatorname {range}R-\operatorname {range}R_T\right|
 \le2e^{-\delta N}.                              \tag{MT.12}
```

*Proof.*  The case `k=1` of (MT.3), (MT.7),
`log cosh t<=t^2/2`, and `mn<=N^2/4` give

```math
 E_{P_{v,y}}e^L\le e^{C_\beta N}.                 \tag{MT.13}
```

Since `(u-T)_+<=e^(u-T)` for `u>=0`,

```math
 E_{P_{v,y}}(L-T_N)_+
 \le e^{-T_N}E_{P_{v,y}}e^L
 \le e^{-\delta N}.
```

This is (MT.11); (MT.12) follows immediately. `square`

Thus the `exp(-Theta(N^2))` bridges responsible for the global MGF ceiling
are operationally irrelevant to the complete `2^m` spiked response, not
merely to an average query.

## 3. Linear output range forces overlap at one bounded capped tilt

Let

```math
 d\mu_y=e_y\,dU_n,
 \qquad e_y={1+z_y^2\over2}.
```

Under the exact channel coupling, `V` is uniform on `{+-1}^m`,
`B|V=v` has law `P_(v,y)`, and `B` has marginal `mu_y^(otimes m)`.  Hence

```math
 R_T(V)=E[L_T(B)\mid V].                          \tag{MT.14}
```

For a bridge edge `a`, let `r_a(B_-a)` be the edge-deleted Gibbs response
from FI.14.  Put

```math
 \Gamma_T(B)=\sum_a[D_aL_T(B)]^2.                 \tag{MT.15}
```

Clipping is one-Lipschitz, so the exact cavity-gradient identity FI.15 gives

```math
 \boxed{
 \Gamma_T(B)\le\Gamma_L(B)
 \le t^2\sum_a r_a(B_{-a})^2.}                   \tag{MT.16}
```

For real `s`, define the *capped moderate tilt*

```math
 {d\Pi_{s,T}\over d\mu_y^{\otimes m}}(B)
 ={e^{sL_T(B)}\over E_{\mu_y^{\otimes m}}e^{sL_T}}.             \tag{MT.17}
```

This is not the all-tilt target falsified by the rank-one bridge: below,
`|s|` is bounded in terms of the declared output gap and `L` is capped at a
linear level.

**Theorem MT.3 (output-to-moderate-overlap dichotomy).**  Fix `eta>0` and
suppose

```math
 \operatorname {range}_{v\in\{+-1\}^m}R(v)\ge\eta N.           \tag{MT.18}
```

For all sufficiently large `N`, there are a sign `sigma in {+-1}` and a
number

```math
 0\le s\le s_\eta:={8\log2\over\eta}             \tag{MT.19}
```

such that

```math
 \boxed{
 E_{\Pi_{\sigma s,T}}\Gamma_T
 \ge c_\eta N,
 \qquad
 c_\eta={\eta^2\over64K_*\log2},
 \qquad
 K_*={41(1+e^{\sqrt2})\over8}.}                  \tag{MT.20}
```

Consequently

```math
 \boxed{
 {1\over mn}E_{\Pi_{\sigma s,T}}\sum_a r_a^2
 \ge {c_\eta N\over t^2mn}.}                    \tag{MT.21}
```

At balanced splits, the right side is a positive constant depending only
on `eta,beta`.  Equivalently, by the cavity/full-Gibbs comparison FI.17,
the capped moderate tilt has a fixed positive two-replica product overlap:

```math
 E_{\Pi_{\sigma s,T},\nu_B^{\otimes2}}
       [\tau^1\tau^2R_XR_Z]\ge c_{\eta,\beta}-O(N^{-1/2}).       \tag{MT.22}
```

*Proof.*  By (MT.12), for large `N` the range of `R_T` is at least
`eta N/2`.  One query word therefore differs from `E_VR_T(V)` by at least
`eta N/4`, in one of the two signs `sigma`.  Since that word has mass
`2^(-m)` and `m<=N`,

```math
 \log E_V\exp\{\sigma s_\eta(R_T-E R_T)\}
 \ge s_\eta{\eta N\over4}-m\log2
 \ge N\log2.                                     \tag{MT.23}
```

Conditional Jensen in (MT.14) transfers this lower bound to the centered
MGF of `sigma L_T` under `mu_y^(otimes m)`.

The dimension-free log-Sobolev theorem FI.5 and its entropy consequence
FI.22 apply to `L_T`.  If

```math
 \psi(u)=\log E\exp\{\sigma u(L_T-EL_T)\},
```

then

```math
 {d\over du}{\psi(u)\over u}
 \le K_*E_{\Pi_{\sigma u,T}}\Gamma_T.             \tag{MT.24}
```

Integrating (MT.24) from zero to `s_eta` and using (MT.23), some
`s in [0,s_eta]` satisfies (MT.20).  Equation (MT.16) gives (MT.21).
Finally FI.17 identifies the normalized full-Gibbs square response with the
two-replica overlap in (MT.22), up to `O(t)=O(N^(-1/2))`. `square`

The same proof gives a rate-sensitive sufficient theorem.  Define

```math
 \rho_N(S)=\sup_{|s|\le S}{1\over mn}
 E_{\Pi_{s,T}}\sum_a r_a^2.                       \tag{MT.25}
```

**Corollary MT.4 (moderate-overlap response bound).**  For every
`0<s<=S`,

```math
 \boxed{
 \operatorname {range}R
 \le2e^{-\delta N}
 +2\left\{{m\log2\over s}
          +K_*t^2mn\rho_N(S)s\right\}.}          \tag{MT.26}
```

If

```math
 s_*=sqrt{{m\log2\over K_*t^2mn\rho_N(S)}}\le S,
```

then

```math
 \boxed{
 \operatorname {range}R
 \le2e^{-\delta N}
 +4\sqrt{K_*m\log2\,t^2mn\rho_N(S)}.}            \tag{MT.27}
```

In particular, at balanced splits, `rho_N(S_N)=O(N^(-alpha))` on a
window `S_N` containing `s_*=O(N^(alpha/2))` gives the power saving

```math
 \operatorname {range}R=O(N^{1-\alpha/2}).        \tag{MT.28}
```

*Proof.*  FI.22, (MT.16), and the definition of `rho_N(S)` give, for both
signs and `0<=s<=S`,

```math
 \log E\exp\{\pm s(L_T-EL_T)\}
 \le K_*t^2mn\rho_N(S)s^2.                       \tag{MT.29}
```

Conditional Jensen transfers the same upper bound to `R_T(V)-ER_T`.
Charging the mass `2^(-m)` of a single query word bounds each one-sided
extreme by

```math
 {m\log2\over s}+K_*t^2mn\rho_N(S)s.
```

Add the two sides and use (MT.12), proving (MT.26).  Optimizing at `s_*`
proves (MT.27)--(MT.28). `square`

## 4. What this changes

The false target was an `o(N)` subgaussian proxy for the untruncated
external pressure at every real MGF parameter.  MT.2 proves that the bridges
killing that target have exponentially small total effect on every output
word.  MT.3 then gives the exact response-specific replacement:

> a linear spiked response is possible only if one capped tilt in the fixed
> interval `[-s_eta,s_eta]` has positive normalized edge-cavity, equivalently
> two-replica product, overlap.

Quantitatively, MT.4 says that a power-saving overlap bound on the matching
moderate window transfers directly to a power-saving response bound.  It is
therefore already in the scale required by the existing basin recurrence;
the missing input is the actual-child overlap estimate, not another
concentration or truncation argument.

This observable is one scalar overlap curve on a compact interval for each
declared spike direction.  It is strictly less information than the full
bridge/Gibbs landscape and does not reconstruct the target-order optimizer.
What remains open is optimizer-specific: prove that these capped
moderate-tilt overlaps vanish uniformly, or exhibit an actual minimizing
child sequence for which one remains positive and then determine whether
that positive overlap creates a coherent branch-(iii) gain.  The rank-one
rare tail and the all-tilt obstruction are no longer part of that statement.
