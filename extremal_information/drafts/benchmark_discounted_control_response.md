# Solution-hidden benchmark: discounted deterministic-control response

**Status.**  The operational state, discounted block algebra, exact response
metrics, and response-rate theorem were frozen on 2026-08-17 before the
literature comparison in Section 8.  Exact rational checks are in
[`verify_discounted_control_response.py`](../experiments/verify_discounted_control_response.py).

## 1. Model and declared continuations

Fix a finite state set `S={1,...,n}` and a discount `0<lambda<1`.  A
one-stage deterministic max-sum controller has actions `a in A_i`, immediate
reward `r(i,a)`, and deterministic successor `f(i,a)`.  Its Bellman operator
on a terminal landscape `v in R^S` is

```math
(Tv)_i=\max_{a\in A_i}\{r(i,a)+\lambda v_{f(i,a)}\}.          \tag{DC.1}
```

Equivalently, put

```math
K_{ij}=\max\{r(i,a):f(i,a)=j\},                              \tag{DC.2}
```

with `K_ij=-infinity` if no action has successor `j`.  A length-`m` block
records its best internally discounted reward `K_ij` conditional on its
initial and final states.  It acts by

```math
(U_{K,m}v)_i=\max_j\{K_{ij}+\lambda^m v_j\}.                 \tag{DC.3}
```

A **depth-`h` continuation query** consists of an arbitrary length-`h`
deterministic block and a probed initial state `i`; its response is
`(U_(K,h)v)_i`.  Hard deterministic blocks, having one allowed zero-reward
path and all other entries `-infinity`, are declared valid.  They can be
replaced by sufficiently unfavorable finite rewards on every bounded value
class.

Two landscapes are absolutely equivalent when every declared query has the
same response.  Projective equivalence is a different experiment in which
one scalar baseline is carried or calibrated separately.

## 2. Discounted max-plus block algebra

### Proposition DC.1 (exact serial composition)

If `(K,m)` is followed by `(L,ell)`, their composite is `(K odot L,m+ell)`,
where

```math
(K\mathbin\odot_m L)_{ik}
=\max_j\{K_{ij}+\lambda^mL_{jk}\}.                           \tag{DC.4}
```

Thus

```math
U_{K,m}U_{L,\ell}=U_{K\mathbin\odot_mL,m+\ell}.              \tag{DC.5}
```

#### Proof

Expanding the left side gives

```math
\max_{j,k}\{K_{ij}+\lambda^mL_{jk}
                  +\lambda^{m+\ell}v_k\},                   \tag{DC.6}
```

which is the right side.  Ordinary max-plus multiplication is recovered at
`lambda=1`; discounting makes duration part of the compositional type.
`square`

## 3. Coarsest exact state and exact response metrics

For two terminal landscapes define the depth-`h` operational distance

```math
d_h(v,w)=\sup_{(K,h),i}
 |(U_{K,h}v)_i-(U_{K,h}w)_i|.                               \tag{DC.7}
```

The supremum ranges over contexts for which both displayed responses are
finite.

### Theorem DC.2 (discounted response isometry)

For every finite `h>=0`,

```math
\boxed{d_h(v,w)=\lambda^h\|v-w\|_\infty.}                  \tag{DC.8}
```

Consequently the coarsest exact absolute response state is the full vector
`v`.  No two unequal vectors merge at any finite depth.

If a common scalar is carried separately, the projective response metric is

```math
\boxed{\bar d_h([v],[w])
=\inf_c d_h(v,w+c{\bf1})
={\lambda^h\over2}\operatorname{osc}(v-w).}                 \tag{DC.9}
```

The exact shape state is therefore `R^n/R 1`, together with the separately
retained baseline when absolute rewards are required.

#### Proof

The maximum inequality gives

```math
\|U_{K,h}v-U_{K,h}w\|_\infty
\le\lambda^h\|v-w\|_\infty.                                \tag{DC.10}
```

Choose a coordinate `j` attaining the sup norm and a zero-reward hard block
that sends the probed state deterministically to `j`.  Its response
difference is `lambda^h(v_j-w_j)`, proving equality.  Exact minimality
follows from the same coordinate probes.  Finally,

```math
\inf_c\|v-w-c{\bf1}\|_\infty
={1\over2}\bigl(\max_i(v_i-w_i)-\min_i(v_i-w_i)\bigr),       \tag{DC.11}
```

which proves (DC.9). `square`

Unlike the undiscounted max-plus response, old approximation error is not
merely nonexpanded: after `h` further stages it is attenuated by exactly
`lambda^h` in the worst case.

## 4. The bounded-reward value class is exactly a cube

Suppose every immediate reward lies in `[-R,R]` and the horizon-`H`
terminal reward is zero.  Put

```math
B_H=R{1-\lambda^H\over1-\lambda},
\qquad B_\infty={R\over1-\lambda}.                           \tag{DC.12}
```

Every horizon-`H` optimal vector belongs to `[-B_H,B_H]^n`.  Conversely,
every vector in this cube is realizable: give state `i` a single self-loop
with constant reward

```math
r_i=v_i{1-\lambda\over1-\lambda^H}.                          \tag{DC.13}
```

For the infinite horizon use `r_i=(1-lambda)v_i`.  These rewards have
absolute value at most `R`, and their returns are exactly `v_i`.  Hence the
union of all bounded-reward value vectors is not merely contained in the
cube; it equals the cube.  This makes the packing converse below intrinsic
to deterministic control rather than to an artificially enlarged table
class.

## 5. Sharp response-rate theorem

Let `M_abs(H,h,epsilon)` be the least number of summary states sufficient to
answer every depth-`h` query on every reward-`R`, horizon-`H` value vector to
absolute error at most `epsilon`.

### Theorem DC.3 (depth-resolved covering and packing)

For `B_H>0`,

```math
\left(\left\lfloor{2B_H\lambda^h\over3\epsilon}\right\rfloor+1\right)^n
\le M_{\rm abs}(H,h,\epsilon)
\le
\left(\left\lceil{2B_H\lambda^h\over\epsilon}\right\rceil+1\right)^n.
                                                                    \tag{DC.14}
```

In particular,

```math
\log_2M_{\rm abs}(H,h,\epsilon)
=\Theta\!\left(
n\log\left(1+{B_H\lambda^h\over\epsilon}\right)
\right).                                                       \tag{DC.15}
```

Here and below the `Theta` form describes the high-resolution regime, for
example `B_H lambda^h/epsilon>=6`, and in particular the regime where this
ratio diverges.  At the other endpoint,

```math
M_{\rm abs}(H,h,\epsilon)=1
\quad\Longleftrightarrow\quad B_H\lambda^h\le\epsilon.       \tag{DC.15a}
```

The upper direction uses the zero vector as the sole representative.  The
reverse direction follows because a hard coordinate query exposes the full
interval `[-B_H lambda^h,B_H lambda^h]`, which one answer can approximate
within `epsilon` only under the displayed condition.  The exact bounds
(DC.14) control the intermediate regime.

Modulo a separately carried scalar, for `n>=2`,

```math
\log_2M_{\rm proj}(H,h,\epsilon)
=\Theta\!\left(
(n-1)\log\left(1+{B_H\lambda^h\over\epsilon}\right)
\right).                                                       \tag{DC.16}
```

#### Proof

By Theorem DC.2, an input sup-norm net of radius
`rho=epsilon/lambda^h` is precisely an `epsilon` response net.  A coordinate
grid of mesh at most `rho` gives the upper bound in (DC.14).

For the converse, take the Cartesian grid in `[-B_H,B_H]^n` with mesh
`3rho`.  Distinct grid vectors have response distance at least
`3epsilon`.  If two used the same summary, one decoded answer accurate to
`epsilon` for both would force their response distance to be at most
`2epsilon`.  Thus every grid vector requires a distinct state.  Section 4
shows that every packing point is a genuine deterministic self-loop value
vector.

For the projective upper bound, anchor `v_n=0`; the other coordinates lie in
`[-2B_H,2B_H]` and can be gridded independently.  For the lower bound, fix
`v_n=0`, grid `v_1,...,v_(n-1)` in `[0,B_H]` with mesh `5rho`, and use
(DC.9).  This changes constants but gives the exponent `n-1`. `square`

### Infinite-depth consequences

Letting `H` tend to infinity in (DC.15) gives the sharp saturated rate

```math
\boxed{
\log_2M_{\rm abs}(\infty,0,\epsilon)
=\Theta\!\left(
n\log\left(1+{R\over(1-\lambda)\epsilon}\right)
\right).}                                                     \tag{DC.17}
```

The displayed entropy law is meant for, say,
`0<epsilon<B_infinity/6`; for coarser error the exact covering number may be
one.

This is finite and independent of planning depth.  In contrast, at
`lambda=1` the bounded-reward value radius grows linearly with the horizon.

The information in an already stored value vector disappears at resolution
`epsilon` after

```math
h\ge
\max\left\{0,
\left\lceil{\log(B_\infty/\epsilon)\over-\log\lambda}\right\rceil
\right\},
                                                                    \tag{DC.18}
```

because then one summary state answers every depth-`h` query within
`epsilon`.  Formula (DC.14) describes the entire decay curve, not only the
contraction endpoint.

## 6. Reusing an approximate Bellman rule forever

Let `T,T'` be `lambda`-contractive Bellman operators with fixed points
`V,V'`, and suppose

```math
\sup_v\|Tv-T'v\|_\infty\le\delta.                           \tag{DC.19}
```

Then

```math
\|V-V'\|_\infty
\le\delta+\lambda\|V-V'\|_\infty,
\qquad
\boxed{\|V-V'\|_\infty\le{\delta\over1-\lambda}.}          \tag{DC.20}
```

This denominator is sharp, not a proof artifact.  For self-loop reward
profiles `r,r'`,

```math
T_rv=r+\lambda v,\qquad
V_r={r\over1-\lambda},qquad
\|V_r-V_{r'}\|_\infty
={\|r-r'\|_\infty\over1-\lambda}.                           \tag{DC.21}
```

Consequently, for the declared `n`-coordinate self-loop reward-profile
family (or a fixed transition/action structure with one state reward per
coordinate), an infinite-use reward-rule code needs and suffices, up to
universal constants and in the nontrivial resolution regime, at

```math
\Theta\!\left(
n\log\left(1+{R\over(1-\lambda)\epsilon}\right)
\right)                                                       \tag{DC.22}
```

bits.  The upper bound quantizes each immediate reward coordinate at scale
`(1-lambda)epsilon`; the self-loop cube supplies the matching packing.
For a controller with an independently variable reward on every action, the
corresponding ambient dimension is the number of such reward parameters,
not automatically `n`.

If a fresh numerical error of size `delta` is introduced after every
Bellman update, the same calculation gives

```math
e_t\le\lambda^te_0+\delta{1-\lambda^t\over1-\lambda}.        \tag{DC.23}
```

Thus infinite repeated re-encoding requires a local error budget at most
`(1-lambda)epsilon` in the worst-case additive-error model; coherent scalar
errors attain equality in (DC.23).  A particular re-encoder may enjoy
cancellation.  This computational statement is distinct from the one-time
response code in Theorem DC.3.

## 7. Scope and caveats

1. Arbitrary coordinate-routing futures make the exact state the full value
   vector.  A fixed controller with a restricted action alphabet may admit a
   strictly smaller behavioral or bisimulation quotient.
2. Absolute responses retain the scalar baseline.  The projective exponent
   `n-1` applies only when that scalar is stored or calibrated separately.
3. The packing ranges over all bounded reward profiles.  Structured reward
   classes can have smaller metric entropy.
4. Theorem DC.3 concerns a value summary encoded once and then reused.
   Formula (DC.23) separately accounts for fresh error injected by every
   re-encoding.
5. The exact state never merges at a finite depth because `lambda^h>0`.
   One-state collapse is an `epsilon`-response statement.

## 8. Post-freeze literature comparison

The block contraction and geometric fixed-point denominator are classical.
Denardo developed monotone contraction operators as a common foundation for
discounted dynamic programming and successive approximation:
[Denardo, *Contraction mappings in the theory underlying dynamic
programming*](https://doi.org/10.1137/1009030) (1967).  Whitt subsequently
gave approximation and aggregation bounds for these monotone contraction
models: [Whitt, *Approximations of Dynamic Programs,
I*](https://doi.org/10.1287/moor.3.3.231) (1978).

Exact state aggregation in a fixed MDP is classically described through
reward- and transition-stable bisimulation partitions:
[Givan--Dean--Greig, *Equivalence notions and model minimization in Markov
decision processes*](https://doi.org/10.1016/S0004-3702(02)00376-4) (2003).
Ferns, Panangaden, and Precup replace exact bisimulation by a quantitative
fixed-point metric and bound optimal-value differences by state distance:
[*Metrics for Finite Markov Decision
Processes*](https://www.cs.mcgill.ca/~prakash/Pubs/Ferns_MetricsForMDPs.pdf)
(UAI 2004).

The independently frozen result agrees with the classical contraction
factor and `1/(1-lambda)` perturbation scale.  Its continuation experiment is
deliberately richer than fixed-MDP bisimulation: arbitrary hard routing can
probe every terminal coordinate, so no model-specific state merge survives.
The additional conclusion is the explicit, matching covering/packing rate
for this declared response experiment, including the exact `lambda^h` decay
with the age of a summary.  The checked classical sources motivate
contraction, approximation, and bisimulation metrics; no novelty claim is
made for those ingredients or for the elementary cube-entropy consequence.

## Benchmark verdict

**Pass, independently predicted.**  Operational future probes recover the
full value vector and the discounted semidirect max-plus composition before
lookup.  More importantly, contraction is converted into a sharp reusable
information law: infinite-horizon complexity saturates at the
`1/(1-lambda)` value radius, an old summary loses all `epsilon`-relevant
information after logarithmic depth, and deterministic self-loops give the
matching converse.
