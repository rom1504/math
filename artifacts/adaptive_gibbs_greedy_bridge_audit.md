# Adaptive Gibbs-gradient bridge audit

Status: exact deterministic finite-temperature identity, exact noisy-code
posterior representation, and a rigorous high-temperature gradient bound.
The construction is a new adaptive bridge heuristic, but no uniform
correlation theorem in the recurrence temperature window or summable
composition defect is proved.

## Exact greedy increment

Fix child blocks `A,D`, a relative child sign, and a partial bridge. With
`H(x,y)` its current parent energy, set

```math
Z=\mathbb E_{x,y}\cosh(\gamma H(x,y)),
\qquad
r_{ij}={\mathbb E[x_i y_j\sinh(\gamma H)]\over Z}.
```

The `y` average includes both relative block orientations. If an unrevealed
edge `(i,j)` is assigned sign `c`, direct expansion gives

```math
{Z_{\rm new}\over Z}
=\cosh\gamma\{1+c r_{ij}\tanh\gamma\}.              \tag{AG1}
```

Thus the locally optimal deterministic sign is `c=-sign(r_ij)`, with exact
log-partition increment

```math
J_\gamma(|r|)
=\log\cosh\gamma+\log(1-|r|\tanh\gamma).           \tag{AG2}
```

Unlike averaging both signs, whose gain is quadratic in `r`, (AG2) has a
linear correlation gain. For a reveal rule producing correlations `r_e`, the
exact recurrence-scale sufficient condition is

```math
\sum_eJ_{\gamma_N}(|r_e|)
\le \gamma_N\Delta_{m,n}+O(t_NN^{1-\delta}),        \tag{AG3}
```

where `gamma_N=t_N/sqrt(N)` and

```math
\Delta_{m,n}
=(M_m^{2/3}+M_n^{2/3})^{3/2}-M_m-M_n.
```

Together with the ordinary soft-to-ground entropy cost and
`t_N>=N^delta`, (AG3) yields a geometrically summable `b`-scale defect. Since
`-log(1-z)>=z`, the following cumulative first-moment inequality is a simpler
exact sufficient condition:

```math
\tanh\gamma_N\sum_e|r_e|
\ge mn\log\cosh\gamma_N-\gamma_N\Delta_{m,n}
   -O(t_NN^{1-\delta}).                             \tag{AG4}
```

In the range `gamma_N=o(1)`, its leading correlation requirement is

```math
\sum_e|r_e|
\ge {mn\gamma_N\over2}-\Delta_{m,n}
+\hbox{controlled lower-order terms}.              \tag{AG5}
```

This asks for correlations of order `t_N/sqrt(N)` on average, rather than the
near-unit polarization required by the sign-averaged overlap route.

## Exact noisy-code and posterior identity

There is a useful exact interpretation of every gradient in (AG1).  Introduce
an orientation spin `s in {+1,-1}` and the latent rank-one bridge word

```math
q_{ij}=s x_i y_j.
```

Let `nu_0` be the child Gibbs law

```math
\nu_0(s,x,y)
={\exp\{\gamma s(H_A(x)+H_D(y))\}
  \over
  \sum_{s',x',y'}
  \exp\{\gamma s'(H_A(x')+H_D(y'))\}}.             \tag{AG5a}
```

Put `a=tanh(gamma)`.  Define a probability law `Pi` on complete bridge sign
vectors `C` by first sampling `(s,x,y)` from `nu_0`, then independently
passing each `q_ij` through the binary channel

```math
\Pr\{C_{ij}=c\mid s,x,y\}={1+c a q_{ij}\over2}.     \tag{AG6}
```

If there are `K=mn` bridge edges, direct summation gives

```math
\Pi(C)={Z(C)\over2^K Z(0)(\cosh\gamma)^K},          \tag{AG7}
```

and, equivalently,

```math
2^K\Pi(C)
=\mathbb E_{\nu_0}\prod_e(1+aC_eq_e)
=e^{-K\log\cosh\gamma}
 \mathbb E_{\nu_0}e^{\gamma\sum_eC_eq_e}.          \tag{AG8}
```

After a prefix of bridge choices, the current lifted Gibbs law is exactly
the posterior distribution of `(s,x,y)` in this channel model.  Therefore

```math
r_e=\mathbb E[q_e\mid C_{<e}],
\qquad
\Pi(C_e=c\mid C_{<e})={1+c a r_e\over2}.            \tag{AG9}
```

The greedy sign is not merely a negative derivative: it is the less likely
next channel output under the exact posterior.

Let `C^g` be the completed greedy bridge, for any chosen edge order.  The
chain rule and (AG9) give the exact cumulative-potential identity

```math
\begin{aligned}
\mathcal S(C^g)
 &:=-\sum_e\log(1-a|r_e|)\\
 &=-\log(2^K\Pi(C^g))\\
 &=K\log\cosh\gamma
   -\log\mathbb E_{\nu_0}
      e^{\gamma\sum_eC^g_eq_e}.                    \tag{AG10}
\end{aligned}
```

Thus `S` is the surprise, relative to the uniform bridge law, of the greedy
least-likelihood path.  Combining (AG2) and (AG10), the left side of (AG3)
is exactly

```math
\boxed{
\sum_eJ_\gamma(|r_e|)
=\log\mathbb E_{\nu_0}
 e^{\gamma\sum_eC^g_eq_e}.}                        \tag{AG11}
```

Consequently AG3 is equivalent to the full Gibbs-response inequality

```math
\log\mathbb E_{\nu_0}
 e^{\gamma\sum_eC^g_esx_iy_j}
\le\gamma\Delta_{m,n}+O(t_NN^{1-\delta}).          \tag{AG12}
```

The elementary inequalities

```math
a u\le-\log(1-au)\le {a u\over1-a}
\qquad(0\le u\le1)                                 \tag{AG13}
```

show that (AG4) is precisely the small-`gamma` linearization of the surprise
requirement in (AG10), not an independent potential estimate.

This identity also fixes the information requirement.  The scalar child
free energy determines the normalizing constant in (AG5a), but the
conditional means in (AG9) require the posterior response of every latent
rank-one word to the entire prefix.  An entropy or mutual-information bound
averaged over a random channel output does not control the adaptively chosen
least-likelihood path.  A potential theorem for an objective-independent
edge order must therefore establish a least-output-likelihood bound for this
noisy rank-one Gibbs code.  Merely rewriting the right side of (AG10) is the
full Gibbs response that the route was intended to compress.

## Symmetry and the genuine high-temperature regime

At the empty bridge, every correlation vanishes exactly at every
temperature:

```math
r_{ij}=0.                                           \tag{AG14}
```

Indeed, `H_A(x)+H_D(y)` is unchanged by the global transformation
`x -> -x`, while `x_i y_j` changes sign.  Thus no positive per-edge or
positive initial-average gradient bound is possible; all useful correlation
must be created by earlier bridge choices.

There is also a scalable regime in which **all** available unrevealed
correlations stay below the leading AG4 scale, independently of the children,
signs, and reveal order.  Condition on `s` in (AG5a).  The current posterior
is a zero-field Ising model on `N=m+n` vertices with couplings of magnitude
`gamma` on the internal and already revealed edges.  Put

```math
\alpha=(N-1)\tanh|\gamma|.
```

The Dobrushin comparison bound, applied to an unrevealed edge (whose direct
coupling is zero), gives

```math
|r_{ij}|
\le {\tanh|\gamma|\,\alpha\over1-\alpha}
\qquad(\alpha<1).                                  \tag{AG15}
```

For completeness, the comparison matrix has entries at most
`tanh|gamma|` on existing edges and row sum at most `alpha`.  Its paths from
`i` to `j` start at length two because `(i,j)` is unrevealed; summing them
gives `tanh|gamma| alpha/(1-alpha)`.  The bound is uniform in `s`, so it
survives their Gibbs mixture.

In particular, if `alpha<=1/3`, every available gradient satisfies

```math
|r_{ij}|\le {1\over2}\tanh|\gamma|
=\left({1\over2}+o(1)\right)|\gamma|.              \tag{AG16}
```

This is a rigorous all-prefix low-correlation construction, but it lies in
the wrong temperature window.  For comparable blocks, the positive part of
AG4 requires

```math
mn\gamma\gg\Delta_{m,n}=\Theta(N^{3/2}),
\quad\text{hence}\quad \gamma\gg N^{-1/2}.         \tag{AG17}
```

The Dobrushin window has `gamma=O(N^{-1})`; equivalently
`t_N=gamma sqrt(N)=O(N^{-1/2})`, whereas the soft-to-ground argument requires
`t_N>=N^delta`.  Therefore (AG15)--(AG16) explain the symmetry barrier at
true high temperature but do not falsify the proposed low-temperature
composition regime.

## Potential-theorem audit

Equations (AG7)--(AG12) leave one exact possible theorem:

> For a fixed objective-independent edge order and exact-minimizer children,
> the greedy least-likelihood word of the noisy rank-one Gibbs code satisfies
> the least-output bound (AG12).

No present identity forces this.  The usual averaged information quantities
describe a bridge drawn from `Pi`, while greedy repeatedly conditions on the
less likely output.  Conversely, proving (AG12) is already an upper bound on
the final soft cap of the assembled signing, with the precise recurrence
constant.  It is not obtained from the child partition function or from
average posterior entropy alone.

The bounded stopping judgment is therefore:

- there is no objective-independent potential inequality presently reducing
  AG4 to a smaller scalar obligation;
- all-prefix small gradients are rigorously forced only in the irrelevant
  Dobrushin window; and
- beyond that window, controlling the adaptive gradients is exactly the
  least-output-likelihood problem for the complete Gibbs rank-one response.

Further work is justified only with a new noisy-code theorem exploiting a
specific exact-minimizer family.  Computing or tabulating the posterior means
in (AG9) is the stopped full-response route.

## Zero-temperature and active-face classification

For a fixed partial bridge, let

```math
\mathcal G
=\{(s,x,y):sH(x,y)=\operatorname{cap}(H)\}
```

be its oriented active face.  At fixed finite order, the formal
`gamma -> infinity` limit of (AG9) is

```math
r_e\longrightarrow
{1\over|\mathcal G|}\sum_{\omega\in\mathcal G}q_e(\omega)
\in\operatorname{conv}
 \{q_e(\omega):\omega\in\mathcal G\}.              \tag{AG18}
```

Thus the Gibbs gradient is a canonical barycenter of active cut gradients.
When `t_N -> infinity` but `gamma_N=t_N/sqrt(N)` does not itself tend to
infinity, the same statement holds with an `o(N^(3/2))` near-active window
rather than the exact finite-order face.  One should not claim exact-face
concentration from scaled low temperature alone: resolving unit energy gaps
would require a much larger inverse temperature.

This resembles the inactive common-active-face route, but the quantifiers are
different.

- A coherent-face assertion asks one fixed law, on one fixed active or
  near-active face, to control many coordinates or edge blocks
  simultaneously.
- Adaptive Gibbs greedy uses only coordinate `e` of the current barycenter,
  then changes the Hamiltonian.  The next gradient can lie in a completely
  different active-face convex hull.  AG4 is a sum
  `sum_e |r_e^(e)|` along these moving faces; it is not the `l_1` norm of one
  fixed bias vector.

Therefore path dependence genuinely avoids the *single coherent face*
obligation.  The verified common-active-face theorem does not settle AG4:
it constructs an existential balancing law with small top-coordinate bias,
whereas AG4 needs large cumulative bias for the canonical sequence of Gibbs
posteriors.  Conversely, AG4 would not produce one law satisfying the earlier
bare-favourability and project-row clauses.  The adaptive route should not be
classified as a direct reuse of that inactive lemma.

The endpoint, however, is not weaker.  Apply the Laplace principle directly
to (AG11).  If `H_0=H_A+H_D`, then

```math
\begin{aligned}
\lim_{\gamma\to\infty}{1\over\gamma}
 \log\mathbb E_{\nu_0}e^{\gamma\sum_eC_eq_e}
 &=\max_{s,x,y}s\left(H_0(x,y)+x^{\mathsf T}Cy\right)
   -\max_{s,x,y}sH_0(x,y)\\
 &=\operatorname{cap}(P_C)-\operatorname{cap}(A\oplus D).
                                                               \tag{AG19}
\end{aligned}
```

After choosing the relative child sign so that the two absolute child maxima
align, `cap(A direct-sum D)=M_m+M_n`.  Hence the zero-temperature limit of
(AG12) is exactly

```math
\operatorname{cap}(P_{C^g})
\le M_m+M_n+\Delta_{m,n},                           \tag{AG20}
```

up to the stated lower-order defect.  This is the desired full bridge cap
bound, not a relaxation of it.

The correct classification is therefore:

1. **not literally the common-active-face condition** -- moving faces remove
   the need for one coherent bias law;
2. **still a complete Gibbs-response condition** -- telescoping those moving
   biases is exactly the final soft cap, and its zero-temperature endpoint is
   full bridge optimization; and
3. **not rescued by the proved common-face law** -- that law has the wrong
   bias direction, resolution, and fixed-law quantifier for AG4.

A new theorem could still exploit path dependence, for example by showing
that anti-predicted noisy-code outputs force posterior concentration after a
controlled number of steps.  Such a theorem would be a genuinely different
ingredient.  Without it, calling the gradients "active-face biases" only
renames the full Gibbs response and does not reactivate the earlier route.

## Reproducible finite test

The program always reveals an unrevealed edge with maximum `|r_ij|`, then
uses the sign from (AG2). It enumerates all parent spin states when evaluating
the Gibbs law and verifies every final integer cap independently.

For equal exact children, minimizing over tested scaled temperatures
`t in {1,2,4,8}` and both relative child signs gives:

| child order | child caps | greedy parent cap | best rigorous parent value/bound | `b`-defect |
|---:|---:|---:|---:|---:|
| 4 | `4,4` | 10 | `M_8=10` | -0.398095 |
| 5 | `4,4` | 13 | `M_10=13` | 0.489091 |
| 6 | `5,5` | 20 | `M_12=18` | 1.520028 |
| 7 | `9,9` | 23 | `M_14=21` | -0.565918 |
| 8 | `10,10` | 32 | witness cap 30 at order 16 | 0.796191 |
| 9 | `12,12` | 39 | conference cap 33 at order 18 | 1.017349 |
| 10 | `13,13` | 48 | witness cap 42 at order 20 | 2.150159 |

The first two exact hits show that the identity is not vacuous. The held-out
orders show no improving scaling law: the absolute gap to the best saved
parent grows to six at orders 18 and 20, and the `b`-defect is not uniformly
negative or decreasing.

## Research judgment

Equations (AG1)--(AG5) are an exact cross-order mechanism with a substantially
weaker-looking correlation scale than near-polarization. What remains is a
uniform lower bound on cumulative *unrevealed* Gibbs gradients under a chosen
edge order. Computing those gradients uses the full `2^(m+n-1)` Gibbs
response, and no bounded-complexity state currently controls them. The finite
degradation means this is not yet a computational scaling law strong enough
to isolate a credible uniform lemma.

The route should be retained only if one proves a potential inequality for
the cumulative absolute gradients, or finds an objective-independent reveal
order whose normalized shortfall in (AG4) decreases with order. More greedy
finite caps alone do not qualify as primary progress.

## Reproduction

```bash
.venv/bin/python computations/phase2k_greedy_gibbs_bridge.py \
  --output computations/results/phase2k_greedy_gibbs_bridge.json
```

The order-nine and order-ten held-outs are in the correspondingly suffixed
result files.
