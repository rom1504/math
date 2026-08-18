# Optimal-product score centering and the cavity-overlap ceiling

Status: **rigorous actual-law identity plus generic sharp falsifiers**.  For
the actual optimized-child escort, the likelihood score of a globally
optimal reverse-KL row product is exactly the cross-row Hoeffding component
of the bridge pressure, computed in the optimal product's own geometry.
This gives the unique exact row-additive centering and an exact interaction
interpolation from the product to the full escort.

The centering is not a lower-information observable: it requires the
unknown optimal product.  Moreover, uncentered cavity overlap has
row-additive false positives even for a centrally symmetric binary-channel
model, and a generic parity example shows that even the optimally centered
quadratic cavity mass cannot lower-bound reverse-product KL.  The parity
example is **not** asserted to be an actual-child pressure; it isolates the
additional optimizer/channel theorem that would be needed.

## 1. Exact optimal-score ANOVA identity

Let `U=\bigotimes_(i=1)^m U_i` be the fair law on bridge rows, let `L` be
the bridge pressure of two actual contracted-temperature minimizing
children, and put

```math
 {dq\over dU}={e^{-\lambda L}\over Z},
 \qquad Z=E_Ue^{-\lambda L}.                         \tag{OS.1}
```

Let `p^*=\bigotimes_i p_i^*` be any global minimizer of
`D(p\Vert q)` over row products, and write

```math
 \mathcal I^{\leftarrow}=D(p^*\Vert q),
 \qquad \mu=E_{p^*}L,
 \qquad
 \ell_i(b_i)=E_{p_{-i}^*}L(b_i,B_{-i})-\mu.          \tag{OS.2}
```

Define the row-interaction residual in the `p^*` geometry by

```math
 L_\perp=L-\mu-\sum_{i=1}^m\ell_i.                  \tag{OS.3}
```

Every optimal factor is strictly positive and its coordinate best response
is

```math
 {dp_i^*\over dU_i}(b_i)
 \propto\exp\{-\lambda E_{p_{-i}^*}L(b_i,B_{-i})\}. \tag{OS.4}
```

### Theorem OS.1 (the reverse-projection score is pure row interaction)

For every finite actual-child bridge problem,

```math
 \boxed{
 E_{p^*}[L_\perp\mid B_i]=0\quad(1\le i\le m),
 \qquad
 \log{dp^*\over dq}
 =\mathcal I^{\leftarrow}+\lambda L_\perp.}          \tag{OS.5}
```

Consequently, for every additive row statistic
`a=c+\sum_i a_i(B_i)`,

```math
 \boxed{
 E_{p^*}\left[
  \left(\log{dp^*\over dq}-\mathcal I^{\leftarrow}\right)a
 \right]=0.}                                        \tag{OS.6}
```

Equivalently, the conditional mean of the likelihood score is constant on
every row:

```math
 E_{p^*}\left[\log{dp^*\over dq}\mid B_i\right]
 =\mathcal I^{\leftarrow}.                          \tag{OS.7}
```

*Proof.*  Independence gives `E_(p^*) ell_i=0` and
`E_(p^*)[ell_k(B_k)|B_i]=0` for `k\ne i`.  Since
`E_(p^*)[L|B_i]=mu+ell_i(B_i)`, (OS.3) has zero conditional
mean on every row.  Summing the logarithms in (OS.4) and subtracting the
log density in (OS.1) shows that `log(dp^*/dq)-lambda L_perp` is a
constant.  Its `p^*` expectation is `D(p^*||q)`, proving (OS.5).
Equations (OS.6)--(OS.7) follow. `square`

This is stronger than first-order orthogonality alone.  Normalization of
`q/p^*` yields the exact nonlinear identity

```math
 \boxed{
 \mathcal I^{\leftarrow}
 =\log E_{p^*}e^{-\lambda L_\perp}.}                \tag{OS.8}
```

For `0\le s\le\lambda`, define

```math
 {d\nu_s\over dp^*}
 ={e^{-sL_\perp}\over E_{p^*}e^{-sL_\perp}}.       \tag{OS.9}
```

Then `nu_0=p^*`, `nu_lambda=q`, and twice differentiating its scalar log
partition gives

```math
 \boxed{
 \mathcal I^{\leftarrow}
 =\int_0^\lambda(\lambda-s)
       \operatorname {Var}_{\nu_s}(L_\perp)\,ds.}  \tag{OS.10}
```

The path stays inside the already proved actual-law regularity class.  If
`t=beta/sqrt(N)`, then

```math
 \log{d\nu_s\over dU}
 =\left(1-{s\over\lambda}\right)\log{dp^*\over dU}
  +{s\over\lambda}\log{dq\over dU}+\text{constant}. \tag{OS.11}
```

Both endpoint log densities have one-bit oscillation at most `2lambda t`.
Therefore every law on the path has the same bound, and after arbitrary
prefix conditioning and marginalization,

```math
 D_2(\nu_s(R_i\mid R_{<i})\Vert U_i)
 \le n\log(1+\tanh^2(\lambda t))
 \le\lambda^2t^2n.                                 \tag{OS.12}
```

Thus escaping row-filtration complexity is not hidden in the exact
optimal-product interpolation either.

There is also an optimizer-specific second-variation constraint.  For
mean-zero row functions `a_i in L^2(p_i^*)`, global minimality under the
product tilt proportional to
`p^* exp(sum_i epsilon_i a_i)` implies that the block quadratic form

```math
 \boxed{
 \sum_i\epsilon_i^2E_{p_i^*}a_i^2
 +2\lambda\sum_{i<k}\epsilon_i\epsilon_k
   E_{p^*}[L_\perp a_i(B_i)a_k(B_k)]\ge0}           \tag{OS.13}
```

is positive semidefinite.  To see this without separating spurious entropy
terms, rewrite the objective difference as

```math
D(p_\epsilon\Vert q)-D(p^*\Vert q)
=D(p_\epsilon\Vert p^*)
 +\lambda E_{p_\epsilon}L_\perp.
```

The first term has Hessian
`Var_(p^*)(sum_i epsilon_i a_i)=sum_i epsilon_i^2E a_i^2`.
Expanding the second term gives
`lambda E[L_perp(sum_i epsilon_i a_i)^2]`; its diagonal terms vanish by
(OS.5), leaving the displayed cross terms.  This controls every pairwise
tangent of the optimal score, but not its higher-row tail.

## 2. Exact cavity centering

For an edge `e=(i,j)`, let `rho=tanh t` and let `r_e(B_{-e})` be the exact
edge-deleted Gibbs response, so that

```math
 D_eL(B)=B_e\operatorname {arctanh}(\rho r_e(B_{-e})).          \tag{OS.14}
```

Only `ell_i` in (OS.3) changes when edge `e` is flipped.  Hence (OS.5)
gives the exact centered-gradient formula

```math
 \boxed{
 D_eL_\perp(B)=B_e\left\{
 \operatorname {arctanh}(\rho r_e(B_{-e}))
 -E_{p_{-i}^*}\operatorname {arctanh}(\rho r_e(B_{-e}))
 \right\}.}                                        \tag{OS.15}
```

In the expectation in (OS.15), the complete row `B_i` is held fixed; only
the other rows are averaged.  Thus the correct centering is not subtraction
of the global mean of `r_e`.  It removes precisely the part predictable
from the queried row itself.  In particular,

```math
 \sum_eE_{p^*}(D_eL_\perp)^2
 =\sum_{i,j}E_{p^*}
 \operatorname {Var}_{p_{-i}^*}
 \left(\operatorname {arctanh}(\rho r_{ij})\mid B_i\right).   \tag{OS.16}
```

Equations (OS.8)--(OS.10), rather than (OS.16), exactly determine the
directed information.  The distinction is essential.

## 3. Uncentered overlap has binary-channel false positives

The failure of raw overlap is not peculiar to an arbitrary scalar
landscape.  Fix a word `c in {+-1}^n`.  Independently for each row, draw a
fair sign `sigma_i` and transmit the planted row
`Q_i=sigma_i c` through a binary channel of amplitude `t`.  The forward
output likelihood relative to `U` is

```math
 p(B)=\prod_{i=1}^m
 {\cosh(tS_i)\over(\cosh t)^n},
 \qquad S_i=\sum_jc_jB_{ij}.                        \tag{OS.17}
```

It is centrally symmetric and is an exact finite binary-channel model.
Its negative escort `q proportional p^(-lambda)U` is a row product, so

```math
 \mathcal I^{\leftarrow}=0.                         \tag{OS.18}
```

On deleting edge `(i,j)`, however, the planted response is

```math
 r_{ij}(B_{i,-j})
 =c_j\tanh\left(t\sum_{k\ne j}c_kB_{ik}\right).    \tag{OS.19}
```

If `n/N->alpha in (0,1)` and `t=beta/sqrt(N)`, the fair overlap obeys

```math
 {1\over mn}\sum_{i,j}E_Ur_{ij}^2
 \longrightarrow E\tanh^2(\beta\sqrt\alpha Z)>0,
 \qquad Z\sim N(0,1).                              \tag{OS.20}
```

The same CLT under the fixed negative row tilts proportional to
`cosh(tS_i)^s`, `s in [-lambda,0]`, shows that their integrated overlap has
a strictly positive limit.  Thus a positive fixed negative-tilt overlap
does not even qualitatively imply irreducible row dependence.  This model
is not the optimized-child rank-one prior: its independent row signs are
precisely the row-additive feature removed by (OS.15).

## 4. Even centered quadratic overlap has a high-row-order ceiling

There is a stronger generic obstruction.  Let `d=mn`, take `m>=2`, let
`U` be fair on all `d` bridge bits, and put

```math
 \chi(B)=\prod_{i,j}B_{ij},
 \qquad t_N={\beta\over\sqrt N},
 \qquad L_N(B)=C_N+t_N\chi(B),                     \tag{OS.21}
```

where `C_N>=t_N` (the additive constant is otherwise immaterial) and
`lambda t_N<1`.  For a row product
`P=\bigotimes_iP_i`, set

```math
 a_i=E_{P_i}\prod_jB_{ij}.                          \tag{OS.22}
```

Data processing to the row parity and the binary entropy inequality give

```math
 D(P_i\Vert U_i)
 \ge {1\over2}\{(1+a_i)\log(1+a_i)
                 +(1-a_i)\log(1-a_i)\}
 \ge {a_i^2\over2}.                                \tag{OS.23}
```

Since `|prod_i a_i|<=|a_1a_2|<=(a_1^2+a_2^2)/2`, after relabeling any two
rows if needed,

```math
 D(P\Vert U)+\lambda t_N\prod_i a_i\ge0.           \tag{OS.24}
```

For `lambda t_N<1`, equality forces `P=U`.  Hence the unique optimal
reverse row product is `p^*=U`, while

```math
 {dq\over dU}={e^{-\lambda t_N\chi}\over\cosh(\lambda t_N)},
 \qquad
 \boxed{\mathcal I^{\leftarrow}
 =\log\cosh(\lambda t_N)=\Theta(N^{-1}).}           \tag{OS.25}
```

The entire nonconstant pressure is already the optimal row-interaction
residual: `L_perp=t_N chi`.  Every edge has

```math
 D_eL_\perp=t_N\chi(B)
 =B_et_N\chi(B_{-e}),                               \tag{OS.26}
```

so it has the formal cavity representation (OS.14) with
`r_e=chi(B_(-e))`.  Therefore, under **every** tilt,

```math
 \boxed{
 {1\over t_N^2mn}\sum_eE(D_eL_\perp)^2=1,}          \tag{OS.27}
```

while (OS.25) tends to zero.  Thus no universal inequality can lower-bound
`I^leftarrow` by the optimally centered quadratic edge/cavity mass at the
physical scale.  The discrepancy is exactly high row order: one amplitude
`t_N` is counted once by the response cumulant and `mn` times by the edge
Dirichlet mass.

The parity landscape has the correct one-bit oscillation `2t_N`, but it is
**not** an actual optimized-child log pressure and need not be a
fixed-amplitude rank-one-channel likelihood.  It therefore does not
falsify an optimizer-specific theorem.  It proves that such a theorem must
use more than weak-coordinate regularity, tight conditional row `D_2`, and
centered quadratic cavity overlap.  In particular it must rule out
high-row-order concentration or control the negative tail in (OS.8).

## 5. Consequence for the smallest missing lemma

The exact optimal centering exists, but it does not yet furnish the requested
low-information branch decider:

1. evaluating (OS.15) requires the complete optimal product `p^*`;
2. its quadratic mass does not control the nonlinear response (OS.8);
3. (OS.13) controls only pairwise tangent instability;
4. the directed information is carried by the negative exponential tail of
   `L_perp` along the regular path (OS.9).

A genuinely narrower next lemma would therefore have to be an
**actual-child anti-high-row-order/tail theorem**, for example a proved
decomposition of `L_perp` into bounded-row clusters whose negative
cumulants add with `o(N)` remainder, or a direct optimizer-specific bound
on the negative tail in (OS.8) from a finite child statistic.  Merely
replacing the raw overlap by the squared centered cavity response is not a
reset.

## 6. Independent audit disposition

An independent reconstruction verified OS.5, OS.8--OS.10, OS.13,
OS.15--OS.16, the pathwise conditional-`D_2` bound, and both falsifiers.
The audit caught and repaired one exposition issue in the Hessian proof:
entropy alone has extra terms along exponential factor tilts, whereas the
clean decomposition through `D(p_epsilon||p^*)` gives OS.13 directly.  No
formula, constant, or conclusion changed.
