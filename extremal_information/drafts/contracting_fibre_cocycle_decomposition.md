# Contracting fibres and visible reward cohomology

**Status.**  Proof draft for independent audit.  The identities, norm bounds,
and finite cycle duality are exercised by
[`../experiments/verify_contracting_fibre_cocycle.py`](../experiments/verify_contracting_fibre_cocycle.py).

This note isolates a general law suggested jointly by the approximate
residual-shell and phase-refresh results.  It is not a claim that one scalar
notion of ``mixing complexity'' controls every future response.  The law is
instead a decomposition:

```math
\boxed{
\begin{array}{c}
\text{terminal centred information}\quad\longrightarrow\quad \rho^t,\\
\text{fresh centred reward error}\quad\longrightarrow\quad(1-\rho^t)/(1-\rho),\\
\text{visible scalar reward cocycle}\quad\longrightarrow\quad
\text{cycle-mean drift}.
\end{array}}
```

The first two channels live inside contracting probability fibres.  The last
channel lies in the invariant constant direction, so fibre mixing cannot
touch it.  Conversely, once scalar cycle flux is paid, arbitrary switching
of the contracting fibres creates no further extensive reward memory.  This
gives a positive compression theorem and an exact pumpable converse in one
model, rather than merely juxtaposing a terminal-response theorem with a
separate counterexample.

## 1. A finite control graph with probability fibres

Let `G=(Q,E)` be a finite directed graph.  At every control vertex `q` put a
finite probability space `(X_q,pi_q)`.  An edge `e:q->q'` carries a Markov
operator

```math
P_e:L^2(pi_(q'))->L^2(pi_q),
\qquad P_e1=1,
\qquad pi_qP_e=pi_(q').                                  \tag{CFC.1}
```

Thus `P_eh(x)` is the expected value of a target-fibre function after taking
edge `e`.  Assume a common strict centred contraction

```math
\|P_eh\|_(2,pi_q)\le rho\|h\|_(2,pi_(q'))
\quad\hbox{whenever }pi_(q')h=0,
\qquad 0\le rho<1.                                      \tag{CFC.2}
```

This allows rectangular kernels, changing fibre sizes, nonreversibility,
and arbitrary switching.  The common stationary laws, not a common matrix,
identify the constant and centred channels.

An edge has a reward residual `a_e in L^2(pi_q)`, and a terminal vertex has
a residual `u_q in L^2(pi_q)`.  These can be the differences between a raw
model and a proposed compressed model.  A reward depending on both endpoints
is included after taking its conditional expectation given the source; the
theorem concerns expected future responses, not pathwise rare deviations.
For a path
`p=e_1...e_t`, with vertices `q_0,...,q_t`, define the full expected response
residual at the initial fibre by

```math
D_p
=a_(e_1)+P_(e_1)a_(e_2)+...+
 P_(e_1)...P_(e_(t-1))a_(e_t)+P_pu_(q_t),            \tag{CFC.3}
```

where the final notation means `P_pu_(q_t)`.  Equivalently, (CFC.3) is the
expected sum of all edge residuals plus the terminal residual.

Write

```math
m_e=pi_(q)a_e,\qquad b_e=a_e-m_e1,
\qquad \bar u_q=pi_qu_q,\qquad v_q=u_q-\bar u_q1,     \tag{CFC.4}
```

and put

```math
B=max_e\|b_e\|_2,\qquad R=max_q\|v_q\|_2,
\qquad U=max_q|\bar u_q|.                             \tag{CFC.5}
```

The numbers `m_e` form a scalar one-cochain on the visible control graph.

## 2. The response decomposition

### Theorem CFC.1 (contracting-fibre/cocycle decomposition)

Under (CFC.1)--(CFC.5), every length-`t` path satisfies

```math
\left\|D_p-
 \left(\bar u_(q_t)+\sum_(s=1)^t m_(e_s)\right)1
\right\|_(2,pi_(q_0))
\le
rho^tR+B{1-rho^t\over1-rho}.                         \tag{CFC.6}
```

If every fibre atom has mass at least `mu>0`, the same right side divided by
`sqrt(mu)` bounds the pointwise centred response error.  Thus (CFC.6) is a
uniform finite-state response theorem as well as an averaged one; retaining
`L^2` makes explicit the information price of very rare hidden phases.

The coefficients and orders of the three terms are separately sharp.

For a scalar edge cochain `m`, define its cycle seminorm on the recurrent
part of `G` by

```math
chi_G(m)=max_(C)
 \left|{1\over|C|}\sum_(e in C)m_e\right|,            \tag{CFC.7}
```

where `C` ranges over directed simple cycles.  If no such cycle exists, set
`chi_G(m)=0`.  Every vertex is permitted as a start in the display below;
for a prescribed start set, replace `G` throughout by its reachable
subgraph.  Define the nonnegative supremum below to be zero when there is no
path of length `t`.  Then, provided the reachable graph contains a directed
cycle,

```math
\limsup_(t->infinity){1\over t}
 \sup_(|p|=t)\|D_p\|_(2,pi_(s(p)))
=chi_G(m).                                            \tag{CFC.8}
```

For an acyclic graph the asymptotic assertion is vacuous: all paths have
bounded length and all response error is a finite transient.

Here allowing arbitrary bounded start or terminal calibrations does not
change the limit.  More quantitatively, if `M=max_e|m_e|` and `G` is strongly
connected, then

```math
\sup_(|p|=t)\|D_p\|_2
\le t\,chi_G(m)+(|Q|-1)M+U+R+{B\over1-rho},           \tag{CFC.9}
```

and repetitions of a cycle attaining (CFC.7) give the reverse asymptotic
bound.

In particular, after retaining/calibrating the terminal means, a
depth-uniform absolute response approximation exists if and only if
`chi_G(m)=0`.  On every strongly connected component this is equivalent to a
potential

```math
m_e=psi(t(e))-psi(s(e)).                              \tag{CFC.10}
```

If `G` is strongly connected and (CFC.10) holds globally, the scalar sum
telescopes and the remaining uniform error is at most

```math
osc(psi)+U+R+{B\over1-rho}.                           \tag{CFC.11}
```

On a general graph, the same statement holds on each recurrent component;
cycle deletion adds a finite acyclic-transient constant.  Vanishing cycle
sums need not give one global potential: an acyclic diamond may have unequal
coterminal path sums, but it cannot be repeated and hence does not affect
(CFC.8).

The conclusion has a useful compression interpretation.  A simulator need
not retain the hidden fibre history.  It retains the visible control state,
the scalar toll's stationary-flow functional (equivalently its recurrent
cohomology class), and a centred residual shell at the single effective
scale `B/(1-rho)`.  At that declared approximation
scale this can be strictly less information than the complete path-response
landscape whenever the fibres or the path tree grow.  No exact
information-minimality claim is made.

#### Proof

Condition (CFC.1) sends constants to constants and centred functions to
centred functions.  Expanding (CFC.3) with (CFC.4) therefore gives the exact
identity

```math
D_p=
 \left(\bar u_(q_t)+\sum_s m_(e_s)\right)1
 +P_pv_(q_t)
 +\sum_(s=1)^tP_(e_1)...P_(e_(s-1))b_(e_s).          \tag{CFC.12}
```

The two final terms have norms at most `rho^tR` and
`sum_(s=1)^t rho^(s-1)B`, proving (CFC.6).
They are centred under `pi_(q_0)`, and hence orthogonal to the displayed
constant.  In particular, if
`s_p=bar u_(q_t)+sum_s m_(e_s)` and `Z_p=D_p-s_p1`, then

```math
\|D_p\|_2^2=|s_p|^2+\|Z_p\|_2^2.                     \tag{CFC.12a}
```

Thus the scalar cycle lower bound cannot be cancelled by a centred fibre
error.

Delete directed cycles successively from the visible control path.  The
remaining simple path has length at most `|Q|-1`; every deleted cycle has
absolute scalar sum at most its length times `chi_G(m)`.  Combining this
with (CFC.6) proves (CFC.9).  Conversely, repeat a reachable cycle attaining
`chi_G(m)`.  Its scalar part grows by exactly the cycle sum per repetition,
while (CFC.6) bounds the entire centred part independently of the number of
repetitions.  This proves (CFC.8).

Zero cycle sum is equivalent to (CFC.10) on a strongly connected directed
graph: fix a root and define `psi(q)` by the `m`-sum along any directed walk
from the root to `q`; appending a return walk and using zero cycle sums makes
the definition path-independent.  Telescoping proves (CFC.11).  The same
argument component by component leaves only a bounded acyclic transient.
`square`

### Stationary-flow dual and a cheap falsifier

Let `F(G)` be the polytope of normalized nonnegative stationary edge flows,

```math
F(G)=\{theta>=0:\sum_e theta_e=1,
\ \sum_(e:s(e)=q)theta_e=\sum_(e:t(e)=q)theta_e\}.    \tag{CFC.13}
```

Cycle decomposition of circulations gives the exact dual formula

```math
chi_G(m)=max_(theta in F(G))\left|\sum_e theta_em_e\right|.
                                                               \tag{CFC.14}
```

Both sides are defined as zero when `F(G)` is empty.

Consequently one does not need to find the worst cycle to falsify a proposed
scalar quotient.  Any stationary control law with edge flow `theta` and
nonzero observable excess `sum theta_em_e` certifies at least that much
asymptotic response rate.  Optimizing over stationary laws is exact, not
merely a lower bound.

This is the sense in which a stationary observable and a pumpable holonomy
are the same scalar obstruction.  It is also the precise limit of the
analogy: stationary averaging controls only the invariant scalar channel;
the fibre-centred channel is controlled by (CFC.2).

### Corollary CFC.1b (response entropy at the forgetting scale)

There is a quantitative description law behind CFC.1.  Assume for this
corollary that the probability fibres and kernels satisfying
(CFC.1)--(CFC.2) are fixed, that `G` is fixed and strongly connected, and
write

```math
r_G=|E|-|Q|+1.                                         \tag{CFC.14g}
```

Let `mathcal U` be a class of terminal dictionaries, metrized by maximum
`L^2` distance over vertices, and let `mathcal B` be a class of **centred**
edge-reward dictionaries, metrized by maximum `L^2` distance over edges.
Write their external covering numbers as `N_U(eta)` and `N_B(eta)`.  Let the
recurrent scalar cochains range over the radius-`L` ball of

```math
H_G=R^E/\{psi(t(e))-psi(s(e))\},
\qquad\|[m]\|_(cyc)=chi_G(m).                         \tag{CFC.14h}
```

Because `G` is strongly connected, this is an `r_G`-dimensional normed
space.  For every `epsilon>0`, the whole response class has a codebook,
composable **within this fixed certified carrier and modulo the exact
endpoint gauges**, of size at most

```math
\boxed{
N_U(epsilon)\,
N_B((1-rho)epsilon)\,
\left(1+{2L\over epsilon}\right)^(r_G)}              \tag{CFC.14i}
```

such that, for every length-`t` future control path, one code centre answers
the expected response within

```math
epsilon t+C_Gepsilon                                  \tag{CFC.14j}
```

in `L^2`, modulo the exact vertex-potential/start calibration.  Here `C_G`
depends only on the fixed graph (one may take a bounded linear section of
`H_G` and then use loop erasure), not on the fibre size, path depth, or
members of the response class.

The quotient convention is substantive.  If absolute endpoint potentials
are charged rather than treated as response gauges, their bounded class must
be covered as an additional finite-dimensional terminal dictionary; its
covering number multiplies (CFC.14i).  No continuous calibration is being
silently counted as one finite state.

Equivalently,

```math
log N_(dyn)(epsilon)
\le log N_U(epsilon)
 +log N_B((1-rho)epsilon)
 +r_Glog(1+2L/epsilon).                              \tag{CFC.14k}
```

Thus the static centred reward image is resolved at the **forgetting
scale** `(1-rho)epsilon`, while persistent composition-created information
is resolved in the cycle quotient and costs its cycle rank.  These are
additive description resources, not a product of two informal
complexities.

#### Proof

Choose the stated terminal and centred-reward nets.  If an external centre
of the latter is not centred, orthogonally project it by subtracting its
`pi_q`-mean edge by edge; this cannot increase its distance from any centred
dictionary.  The terminal error is at most `epsilon`, because every Markov
prefix is an `L^2` contraction.  By CFC.6, fresh centred reward error at most
`(1-rho)epsilon` contributes at most `epsilon` at every depth.

The quotient in (CFC.14h) has dimension `r_G`: the directed incidence map
has rank `|Q|-1`.  The standard volume argument in an `r_G`-dimensional
normed space covers its radius-`L` ball by at most
`(1+2L/epsilon)^(r_G)` balls of radius `epsilon`.  Fix a linear section
`s:H_G->R^E`.  Finite dimensionality gives a graph constant `K_G` with

```math
\|s([d])\|_infinity\le K_G\|[d]\|_(cyc).             \tag{CFC.14l}
```

Replace a cochain difference by this representative, absorbing the removed
gradient into the exact vertex/start calibration.  Loop erasure bounds its
path sum by `epsilon t+(|Q|-1)K_Gepsilon`.  Adding the two centred terms
proves (CFC.14j), for example with `C_G=2+(|Q|-1)K_G`.

The scales and the cycle exponent cannot be improved distribution-freely.
The two-state eigenmode examples in Section 4 attain the terminal
`rho^t` and reward `1/(1-rho)` factors.  Repeating an exposed cycle makes the
response-rate distance between two scalar classes exactly their cycle-norm
distance.  If code centres themselves are cochain classes, rate-error
`epsilon` requires an `epsilon`-cover of the `r_G`-dimensional cycle ball.
For a completely arbitrary encoder/decoder, two inputs assigned one code
state can only be concluded to lie within `2epsilon`, so a strict
`2epsilon`-packing gives the corresponding state lower bound.  Either form,
by the reverse volume estimate, forces
`Omega_G((L/epsilon)^(r_G))` states/centres for the full ball, with different
graph-dependent constants. `square`

This corollary does **not** discover a congruent finite state for arbitrary
hidden dynamics.  It is a response-description theorem once the common-law
contracting carrier has been certified.  Its nonclassical programmatic
content is the exact allocation of accuracy among three independently
exposed resources.

### Corollary CFC.1a (nonlinear stochastic-secant lift)

The same law applies to genuinely switching nonlinear response dynamics.
For clarity, let every hidden profile lie in `R^X` with one full-support law
`pi`, and let a visible edge `e` carry two additively homogeneous maps
`F_e,Fhat_e:R^X->R^X`.  Suppose that for every `x,y` there is a row-stochastic
secant

```math
F_e(x)-F_e(y)=P_e[x,y](x-y),                          \tag{CFC.14a}
```

such that

```math
pi P_e[x,y]=pi,
\qquad
\|P_e[x,y]h\|_(2,pi)\le rho\|h\|_(2,pi)
\quad(pi h=0),                                       \tag{CFC.14b}
```

uniformly in the edge and the secant.  Assume also that the same-input
defect

```math
eta_e(y)=F_e(y)-Fhat_e(y)                             \tag{CFC.14c}
```

has a control-visible mean and a bounded centred part,

```math
pi eta_e(y)=m_e,
\qquad\|eta_e(y)-m_e1\|_(2,pi)\le B                  \tag{CFC.14d}
```

for every `y`.  If true and approximate trajectories follow the same visible
path and their initial difference has centred norm at most `R`, then after
`t` steps their response difference, after subtracting the fully propagated
scalar component `pi(x_0-y_0)+sum_s m_(e_s)`, obeys

```math
\left\|x_t-y_t-
 \left(pi(x_0-y_0)+\sum_(s=1)^tm_(e_s)\right)1
\right\|_(2,pi)
\le rho^tR+B{1-rho^t\over1-rho}.                    \tag{CFC.14e}
```

Therefore (CFC.7)--(CFC.14) give the exact extensive error of every
repeatable visible switching word.  For paired trajectories following the
same declared visible path, optimizer switches and ties cause no extra
information growth once all the uniform common-law and visible-mean
hypotheses (CFC.14b)--(CFC.14d) are verified.

Indeed, with `z_s=x_s-y_s`, add and subtract `F_(e_s)(y_(s-1))` to obtain

```math
z_s=P_(e_s)[x_(s-1),y_(s-1)]z_(s-1)+eta_(e_s)(y_(s-1)).
                                                               \tag{CFC.14f}
```

Equations (CFC.14b)--(CFC.14d) propagate its constant part exactly and
contract its centred part, so the proof of CFC.6 applies verbatim.

For all-finite max-plus maps, stochastic secants always exist, including
across ties.  What is not automatic is (CFC.14b): selectors commonly have
centred norm one and need not share an invariant law.  Nor is (CFC.14d)
automatic: a hidden-state-dependent mean defect is precisely an additional
uncontracted scalar feature.  These two generator-level checks make the
corollary a falsifiable nonlinear lumpability criterion rather than an
assumption that the desired quotient already exists.

## 3. A one-sided semantic-response toll

The preceding theorem is an equality for expected reward recursions.  The
operator phase-refresh setting supplies only a one-sided response inequality.
There is nevertheless a query-mass-sensitive lower bound that does not use a
minimum atom mass.

### Theorem CFC.2 (variance tax for one-sided refresh)

Let `P` preserve a probability `pi` on a finite set and satisfy

```math
\|P-Pi\|_(2->2)\le rho<1.                              \tag{CFC.15}
```

Let `g,f` be real functions, `\|f-g\|_infinity<=omega`, and suppose

```math
f\le Pf+epsilon1,
\qquad epsilon>=0.                                    \tag{CFC.16}
```

Then, writing `sigma_g^2=Var_pi(g)` and `B_g=osc(g)`, one has

```math
\boxed{
epsilon\ge
{(1-rho)^2(\,sigma_g-omega\,)_+^2\over B_g+2omega}.}  \tag{CFC.17}
```

The right side is interpreted as zero when its denominator vanishes.
Thus a static recovery radius `omega`, a dynamic spectral gap `1-rho`, and
the transfer toll cannot be chosen independently.  A response phase with
macroscopic variance either survives in the state/recovery error or is paid
as fresh toll.

#### Proof

Put `h=f-Pf`.  Stationarity gives `pi h=0`, (CFC.16) gives `h<=epsilon`, and
both `f` and `Pf` lie between `min f` and `max f`, so
`h>=-B_f`, where `B_f=osc(f)`.  For a mean-zero random variable in
`[-B_f,epsilon]`,

```math
E h^2\le B_fepsilon.                                  \tag{CFC.18}
```

Indeed, when `epsilon<=B_f`, expand
`(epsilon-h)(h+B_f)>=0` and take expectations; when
`epsilon>=B_f`, use `h^2<=B_f^2<=B_fepsilon`.

For `v=f-Pi f`, (CFC.15) gives

```math
(1-rho)\|v\|_2\le\|v-Pv\|_2=\|h\|_2
\le\sqrt{B_fepsilon}.                                \tag{CFC.19}
```

Orthogonal centring is an `L^2` contraction, hence
`\|v\|_2>=sigma_g-omega`; also `B_f<=B_g+2omega`.
Substitution proves (CFC.17). `square`

The point-mass/mixing-time bound in the expander phase-refresh draft is
stronger when the nonconstant response is concentrated on an atom of mass
`Theta(1/S)`: it produces an exponential state/toll tradeoff.  CFC.2 is
different.  It is insensitive to the parametrization and remains useful
when a positive amount of stationary query mass has nonzero response
variance, even if no individual phase has a useful mass lower bound.

## 4. Why contraction cannot pay scalar holonomy

The decomposition is not an upper-bound artifact.  Each channel is attained
by a two-state or one-state system.

1. **Terminal mode.**  On two equiprobable states take the symmetric kernel
   with mean-zero eigenvalue `rho`, and take `u=(R,-R)`.  Then
   `\|P^tu\|_2=rho^tR`.
2. **Fresh centred mode.**  With the same kernel and
   `a_e=(B,-B)`, repeated use of one edge gives exactly
   `B(1-rho^t)/(1-rho)` in centred norm.
3. **Scalar cocycle.**  Use singleton fibres on a directed cycle and put a
   scalar residual of mean `c` on each edge.  Then `rho=0`, `B=R=0`, while
   the response error is exactly `t|c|`.

Cartesian products let all three obstructions coexist with independently
chosen scales (without claiming equality in the triangle inequality in
(CFC.6)).  Therefore no universal inequality can trade a smaller static
response image or stronger fibre contraction for a nonzero scalar cycle
flux.  The only
ways to remove the extensive term are to store/pay the missing visible toll,
make it a potential, or change the declared future query so that scalar
amplitude is quotiented out.

This also supplies a decisive counterexample to a tempting single-number
law such as

```math
\text{reusable complexity}
 \asymp \text{static response entropy}\times\text{forgetting time}.
```

That product can describe the centred channel, but an independent
one-dimensional scalar cocycle has zero static projective radius and zero
fibre forgetting time while retaining arbitrary positive response rate.

## 5. Benchmark consequences

### 5.1 Switching Markov rewards and hidden-state expectations

For a finite controlled Markov model whose kernels share transported fibre
laws and a centred `L^2` contraction `rho`, approximate each immediate reward
by a scalar toll plus a centred codeword of `L^2` error `B`.  CFC.1 gives one
depth-independent error shell, after the scalar terminal means are retained,

```math
B/(1-rho)+R,                                          \tag{CFC.20}
```

under every switching word.  Without that terminal calibration, add `U`.
There is no need to store the full posterior
path.  This is stronger than applying a terminal mixing estimate at the end:
it proves that *fresh reward errors at every time* remain bounded.  The only
extensive datum is the stationary-flow functional, equivalently the
cohomology class on the recurrent strongly connected components, detected
exactly by (CFC.14).  A global `C^1(G)/B^1(G)` class would overcount bounded
acyclic coterminal-path defects.

For a discounted Markov recursion
`D_s=a_s+lambda P_sD_(s+1)`, `0<=lambda<1`, the two channels have different
denominators.  A scalar fresh residual bounded by `M` costs at most
`M/(1-lambda)`, while a centred residual bounded by `B` costs at most
`B/(1-lambda rho)`.  There is no cycle drift at positive discount.  CFC.7
therefore describes a genuine transition at `lambda=1`, not merely a poor
constant in the discounted estimate.

### 5.2 Exact rank-one max-plus resets and width-one switching Ising

For exact max-plus rank-one blocks, first promote the previous/right-profile
type and its directed compatibility to the visible control graph, using
ARS.2/Theorem 17.7.  The remaining hidden fibre is then a singleton and hence
the `rho=0` endpoint.  This promotion is essential: an identical-row
max-plus secant can select a state-dependent row law and need not preserve
one common full-support `pi`.  After promotion, the scalar edge reward is the
directed compatibility

```math
m_(e,f)=max_j(p_e(j)+a_f(j)).                         \tag{CFC.21}
```

CFC.7--CFC.14 reproduce the *rate* part of the compatibility law without
retaining any microscopic path selector; they do not independently discover
the visible profile graph.  For the sharp two-letter matrices
in the approximate residual-shell draft, subtracting the optimal one-state
letter tolls `(delta/4,5delta/4)` leaves edge residuals

```math
{1\over delta}d=
\begin{pmatrix}-1/4&3/4\\-1/4&-1/4\end{pmatrix}.      \tag{CFC.22}
```

The `A`, `B`, and alternating `AB` stationary cycle flows give
`chi(d)=delta/4`.  Hence exact projective reset (`rho=0`) and an arbitrarily
small one-profile shell do not reduce the scalar rate.  This is precisely the
invariant-channel endpoint of CFC.1, not a failure of its contraction
hypothesis.

### 5.3 Expander semantic phases

Take one control vertex, no scalar reward, and let `g` be a homogeneous
same-scale semantic phase response.  If a recovered response `f` obeys the
one-sided Boolean-pullback inequality, CFC.17 gives

```math
epsilon\ge
{(1-rho)^2(\sqrt{Var_pi(g)}-omega)_+^2
 \over osc(g)+2omega}.                                \tag{CFC.23}
```

Gauge-orbit phases have zero semantic variance and escape, correctly.  A
homogeneous Walsh phase sample with positive stationary response variance
cannot escape at vanishing recovery and toll under a fixed expander gap.  In
the actual scale-varying Walsh hierarchy, applying (CFC.23) additionally
requires a uniform comparison of consecutive recovered responses; its error
is added to `epsilon`.  The time-inhomogeneous ER.2 theorem remains the
unconditional statement.  Thus (CFC.23) is a mass-sensitive,
coordinate-free complement, not a replacement.

## 6. Exact scope and novelty accounting

The ingredients have classical ancestors:

* geometric summation of centred Markov rewards and Poisson-equation bounds;
* circulation decomposition into directed cycles;
* the finite cycle/coboundary criterion already used in Theorems 17.1e,
  17.1h, and 17.1l;
* elementary `L^2` spectral-gap estimates;
* Theorem 16.18's stochastic-secant representation and suffix-gain law.

The theorem-level contribution is the **joint, independently sharp response
law** (CFC.6)--(CFC.14).  It identifies the maximal part of an additive
future-response landscape which contracting hidden fibres can erase, and
proves that what remains extensively is exactly the visible stationary-flow
functional/recurrent cohomology class.  It
also yields the variance toll (CFC.17), which complements rather than repeats
the point-mass phase-refresh lower bound.  The nonlinear lift adds a condition
not present in the generic secant-gain theorem: all realized secants must
transport one law, which canonically separates the invariant scalar channel
from the contracting centred channel.

The result is not a theorem about adversarial max-plus fibres with a merely
small projective image.  There, small residual shells can encode arbitrary
weighted automata, as ARS.4 proves.  Nor does it say that an arbitrary dense
quadratic bridge has a contracting probability-fibre presentation.  Its
falsifiable structural hypothesis is precisely (CFC.1)--(CFC.2): a common
transported family of laws whose centred modes contract under every allowed
switch.  If such laws exist, no further hidden additive memory is needed;
if they do not, this theorem supplies no compression by terminology alone.

## 7. Next discriminating question

The strongest extension is a controlled nonlinear analogue.  Suppose the
hidden update is a max-plus or Bellman map whose stochastic secants all
transport a common family `pi_q` and contract its centred `L^2` space by
`rho<1`.  Does every coherent reward residual admit the same decomposition
into a geometrically bounded centred channel and a finite visible scalar
cocycle?  Theorem 16.18 controls each realized secant path, but common
transported laws across switching ties are not automatic.  A positive
answer would move CFC.1 from linear expected rewards to genuine switching
extremal dynamics; a two-map counterexample would identify exactly where
nonlinearity creates new composition information.
