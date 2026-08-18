# Exact finite audit of actual child-induced bridge escorts

Status: **complete finite enumeration with numerical transcendental
evaluation**.  This is not an asymptotic theorem and does not use interval
arithmetic.  It uses the actual contracted-temperature pressure minimizers;
no conference, Paley, ground-state, or heuristic child is substituted.

The frozen protocol is
[`actual_child_bridge_law_protocol.md`](actual_child_bridge_law_protocol.md).
The full reproducible outputs are
[`actual_child_bridge_law_exact.json`](../../computations/results/actual_child_bridge_law_exact.json),
[`actual_child_row_product_shadow.json`](../../computations/results/actual_child_row_product_shadow.json),
[`actual_child_row_product_shadow_target_threshold_n8.json`](../../computations/results/actual_child_row_product_shadow_target_threshold_n8.json),
and the compact machine summary
[`actual_child_bridge_law_summary.json`](../../computations/results/actual_child_bridge_law_summary.json).

## 1. Law and exact finite scope

For `N=m+n`, `t=beta/sqrt(N)`, actual child minimizers `A,C`, and fixed
relative orientation `epsilon`, the enumerated pressure is

```math
f(B)=\log\left(2^{-N}\sum_{x,y}
 \cosh(t[H_A(x)+\epsilon H_C(y)+x^TBy])\right).
```

The negative escort and the positive output law are

```math
{dq_\lambda\over dU}(B)
={e^{-\lambda f(B)}\over E_Ue^{-\lambda f}},
\qquad
{dp\over dU}(B)={e^{f(B)}\over E_Ue^f}.
```

Every root-gauged child signing and every bridge was enumerated for every
split with `2<=m<=n`, `4<=N<=9`,

```text
beta   = 0.25, 0.5, 1, 2, 4,
lambda = 0.25, 0.5, 1, 2, 4,
epsilon= -1,+1.
```

There are `120` child-pair/orientation pressure cubes and `600` escort laws.
The largest is the balanced order-nine cube with `2^20` bridges.  Integer
energies and the finite enumerations are exact.  Child pressure histograms
were compared at `80` decimal digits.  Gibbs weights and information
quantities are floating evaluations.

An XOR--Walsh convolution evaluates the whole bridge cube.  Fixing the first
spin in each child gives the exact reduction

```math
\overline Z(B)
=E_{x_0=y_0=1}
  \cosh(t[H_A(x)+\epsilon H_C(y)])\cosh(tx^TBy),
```

which is a radial XOR convolution over bridge words.  Three direct spin-sum
checks per cube agree to at most `1.77e-15` in log pressure.

## 2. The finite thermal child selector is rigid

Across all `135` distinct child-order/raw-temperature queries used either as
contracted children or same-temperature targets, the minimizing set has
exactly **one** signed-permutation/global-sign class.  The largest child
enumeration contains all `32,768` root-gauged order-seven signings.  Whenever
there is a distinct competing energy histogram, the smallest 80-digit
pressure gap is `1.26818e-6`.

This is genuine finite structure of the actual thermal optimizer.  It gives
no all-order uniqueness or rigidity theorem.

## 3. Exact target reachability is already order-sensitive

The table gives, after optimizing over every split and both orientations,

```math
{\min_B f(B)-[P_m(\beta)+P_n(\beta)]\over N}.
```

| `N` | `beta=.25` | `beta=.5` | `beta=1` | `beta=2` | `beta=4` |
|---:|---:|---:|---:|---:|---:|
| 4 | +.00385642 | +.01487366 | +.05271180 | +.15889545 | +.41991476 |
| 5 | +.00298494 | +.01044068 | +.02509358 | +.00021240 | -.20475484 |
| 6 | +.00235715 | +.00690010 | +.00259527 | -.10751334 | -.53396166 |
| 7 | +.00211359 | +.00721485 | +.01594097 | -.00376594 | -.11344230 |
| 8 | +.00180533 | +.00574273 | +.00927644 | -.01989040 | -.13629062 |
| 9 | +.00158291 | +.00488847 | +.00833628 | +.00612520 | +.01136436 |

At order nine every exact bridge misses the exact child target at every
queried beta and split.  This is a finite falsifier to a zero-error
all-order bridge statement, not to the desired sublinear-error statement.
Indeed, at `beta=.25` the **unnormalized** best gap stays near `.014` over
these orders, which is compatible with `O(1)=o(N)` error.

Exactly `19` of the `120` pressure cubes have a finite inverse-disorder
temperature at which their negative-moment soft pressure equals the exact
child target.  The most informative comparable example is

```text
N=8, split 4+4, beta=4, lambda*=5.38210:
KL(q||U)/N                         = 0.725981
row total correlation / N         = 0.585143
best one-row-latent residual TC/N  = 0.304747
best-filtration max mean row D2    = 2.64269
chain-support proxy / 16 signs     = 0.756--0.764.
```

Thus one finite, genuinely target-reaching actual escort retains substantial
dependence after conditioning on any single complete row.  This does not
prove a fixed asymptotic mass of irreducible dependence.

## 4. Resource scaling at fixed `lambda=1`

For the balanced split and `epsilon=-1`, the table records the actual escort.
`s/d` is the smallest chain-support proxy among the preregistered coordinate
orders divided by `d=mn`.

| `N` | split | `beta` | `KL/N` | row `TC/N` | max conditional row `D2` | `s/d` | soft-target gap / `N` |
|---:|:---:|---:|---:|---:|---:|---:|---:|
| 6 | 3+3 | 1 | .001207 | .000966 | .020782 | .783316 | +.041418 |
| 6 | 3+3 | 2 | .030045 | .023402 | .314285 | .782597 | +.093521 |
| 6 | 3+3 | 4 | .232187 | .182470 | 1.438514 | .786276 | -.048230 |
| 7 | 3+4 | 1 | .001141 | .000863 | .048369 | .762779 | +.040171 |
| 7 | 3+4 | 2 | .023872 | .018318 | .413477 | .753860 | +.117456 |
| 7 | 3+4 | 4 | .125500 | .097394 | 1.359663 | .687798 | +.160003 |
| 8 | 4+4 | 1 | .001126 | .000932 | .071063 | .804527 | +.039845 |
| 8 | 4+4 | 2 | .023217 | .019664 | .543256 | .793865 | +.134740 |
| 8 | 4+4 | 4 | .119096 | .102772 | 1.498485 | .755836 | +.259065 |
| 9 | 4+5 | 1 | .000932 | .000814 | .059638 | .792047 | +.042152 |
| 9 | 4+5 | 2 | .020328 | .017997 | .509909 | .777137 | +.169667 |
| 9 | 4+5 | 4 | .114452 | .102644 | 1.814578 | .736222 | +.424844 |

All `600` grid escorts have positive row total correlation numerically.
At the strongest balanced point `(N,beta,lambda)=(9,4,4)`, conditioning on
the best one of the four rows still leaves `2.80438` nats of residual total
correlation.  Along its best row filtration, `90.67%` of the final-prefix
mass has conditional row `D2>2`.

The finite conditional complexities remain compatible with the proved
actual-law bound

```math
D_2(q(R_i\mid R_{<i})\Vert U_n)
\le\lambda^2\beta^2n/N.
```

The largest observed/bound ratio over every prefix and ordering is `0.5564`.
The data therefore support tight conditional row components at fixed
`beta,lambda`; they simultaneously reject exact row independence.

## 5. Sharp finite falsifier: the canonical reversed product mixture is far

At `lambda=1`, central symmetry identifies the canonical reversed-channel
bounded-component mixture with the positive output law `p`.  Hence its
distance from the actual escort `q_1` is directly enumerable.

| `N` | split | `beta` | `D(U||p)/N` | `D(q_1||p)/N` | `TV(q_1,p)` | affinity |
|---:|:---:|---:|---:|---:|---:|---:|
| 6 | 3+3 | 2 | .052973 | .148505 | .513354 | .767020 |
| 6 | 3+3 | 4 | .630535 | 1.065598 | .985961 | .085603 |
| 7 | 3+4 | 2 | .060263 | .138803 | .573154 | .727090 |
| 7 | 3+4 | 4 | .538529 | .979953 | .980510 | .077978 |
| 8 | 4+4 | 2 | .060562 | .137129 | .605090 | .695739 |
| 8 | 4+4 | 4 | .548783 | .970748 | .986221 | .053385 |
| 9 | 4+5 | 2 | .050854 | .117562 | .581045 | .707456 |
| 9 | 4+5 | 4 | .490354 | .889158 | .987994 | .051245 |

At order nine, beta four,

```math
D(U\Vert p)=4.41319,
\quad D(q_1\Vert p)=8.00242,
\quad \|q_1-p\|_{TV}=0.987994.
```

This is the sharpest finite falsifier in the audit: the natural
bounded-`D2` inverse-mixture proposal is almost disjoint from the actual
optimized-child escort.  It does **not** exclude a different latent-product
representation.

The rank-one support lower bound from Theorem IM.2 was also checked:

```math
D(U\Vert p)\ge mn\log\cosh t
-t\sqrt{2mn\log(2^{N-1})}.
```

That bound is negative throughout this beta grid (balanced asymptotic
positivity starts above about `4.71`), whereas the exact law already separates
strongly at beta two.  Thus the theorem's mechanism is valid but its support
union bound is quantitatively very loose for the actual children.

## 6. Exact row-product variational shadow

For every law through `N=8` and the balanced order-nine laws, the separate
audit minimized

```math
E_{p_1\otimes\cdots\otimes p_m}f
+\lambda^{-1}\sum_iD(p_i\Vert U_n)
```

by exact coordinate Gibbs updates.  Each run used uniform, exact escort-row
marginal, and `16` seeded random starts.  The table shows `lambda=1`,
`epsilon=-1`.  The captured fraction is

```math
{E_Uf-\text{best evaluated product objective}\over
 E_Uf-V_\lambda}.
```

| `N` | split | `beta` | captured fraction | evaluated product gain / `N` | candidate reverse-projection upper bound / `N` |
|---:|:---:|---:|---:|---:|---:|
| 5 | 2+3 | 2 | .380952 | .012694 | .020629 |
| 5 | 2+3 | 4 | .308932 | .088992 | .199070 |
| 6 | 3+3 | 2 | .183560 | .006506 | .028935 |
| 6 | 3+3 | 4 | .138654 | .038451 | .238863 |
| 7 | 3+4 | 2 | .264970 | .008160 | .022635 |
| 7 | 3+4 | 4 | .247954 | .047216 | .143207 |
| 8 | 4+4 | 2 | .183028 | .005515 | .024618 |
| 8 | 4+4 | 4 | .163892 | .030119 | .153655 |
| 9 | 4+5 | 2 | .127307 | .003317 | .022736 |
| 9 | 4+5 | 4 | .115911 | .019693 | .150207 |

For all balanced `lambda=1` rows above, every one of the `18` starts reaches
the same objective to `1e-9`, and the simultaneous best-response residual is
below `7e-8`.  Over all `500` shadow problems, `466` have a single terminal
objective to `1e-9`; four best runs reached the preregistered 100-sweep cap.
Focused 1000-sweep reruns reduced those residuals, but do not change the
nonconvex status.

A separate run evaluates the genuinely target-reaching balanced order-eight
law at its exact threshold `lambda*=5.382104`:

```text
best evaluated row-product captured fraction       = 0.763655
evaluated row-product gain / N                      = 0.338176
candidate reverse-projection upper bound / N        = 0.563306
candidate reverse-projection contribution / N/lambda= 0.104664.
```

Thus the large forward total correlation of that escort does **not** imply
that row products explain little of its pressure gain: the best evaluated
product already explains more than three quarters.  Four terminal objectives
occurred across the starts and the best-response residual is `1.66e-7`, so
this remains evidence rather than a certified optimum.

The evaluated product law is globally feasible.  Therefore its gain is a
numerically rigorous-from-the-candidate **lower bound** on the optimal
row-product gain, and its residual variational gap is an **upper bound** on
the directed reverse projection.  Coordinate descent is not a global
certificate.  The decreasing captured fraction is strong finite evidence
that the actual escort's gain is not explained by the easiest row-product
shadow, but it is not a lower bound on irreducible directed dependence.

## 7. Independent numerical sanity checks

- maximum Gibbs KL identity residual: `4.48e-13`;
- `D_2(q||U)>=D(q||U)`, `D(q||U)>=TC_rows,TC_columns`, and mean gain at least
  soft gain in every grid law;
- every conditional row `D2` obeys the proved bound above;
- at `lambda=1`, `D(q||p)>=D(U||p)` and
  `TV(q,p)>=1-exp[-D(U||p)/2]` in every cube;
- every best coordinate-descent objective trace is monotone up to maximum
  upward roundoff `1.07e-13`;
- the exact chain-support subset DP is used only through ten bridge signs.
  Four nearly uniform cases initially reversed the mathematical DP/proxy
  ordering by at most `1.82e-6`; the implementation now records and clamps
  this numerical cancellation explicitly.

## 8. Finite judgment and next discriminating experiment

The audit touches the actual child structure and supplies a real reset from
the conference surrogate.  It supports the following finite picture:

```text
fixed-beta/lambda conditional row D2 is tight,
but the canonical product mixture is far from q,
row total correlation is nonzero and can remain large after one-row latency,
local row-product shadows capture only a minority of the balanced lambda-one
gain at the largest orders tested, but can capture a majority at a finite
target-reaching disorder temperature.
```

What remains unproved is precisely directional and asymptotic: the scaling of

```math
\mathcal I_\lambda^{\leftarrow}
=\inf_{p\text{ row-product}}D(p\Vert q_{\lambda,\epsilon})
```

for actual optimized children.  Forward total correlation cannot replace
this reverse projection.

The most discriminating next exact experiment is a targeted balanced
order-ten (`5+5`) audit at `beta=2,4`, first at `lambda=1` and then at the
target-threshold lambda if it exists.  It should exploit bridge gauge/orbit
compression or streamed Walsh marginals rather than store every filtration.
The deciding measurements are the `q_1`--`p` separation and the globally
certified row-product variational value.  Merely adding another unconstrained
coordinate-descent start is not enough; the next step should pair the exact
`2^25` bridge cube with a branch-and-bound or convex lower certificate for
`V_lambda^row`.
