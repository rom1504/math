# Adaptive deterministic bridge signs: exact gain and remaining obstruction

Status: materially stronger than random-sign interpolation. The exact greedy
rule reduces the required response from near-unit squared overlap to a small
first-moment correlation. No generic entropy or convex potential currently
forces the required cumulative correlation in the growing-temperature
regime.

## 1. Exact one-edge rule

At a bridge prefix, let

```math
r_e=\mathbb E_\mu[\tau x_i y_j],\qquad e=(i,j),
```

under the current extended Gibbs law. If the next sign is `c in {+1,-1}`,
the exact partition-function ratio is

```math
{Z_c\over Z_0}=\cosh\gamma+c r_e\sinh\gamma.         \tag{1}
```

Choose deterministically

```math
c_e=-\operatorname{sign}(r_e)                       \tag{2}
```

(either sign if `r_e=0`). The exact logarithmic increment is

```math
K_\gamma(r_e)=\log\cosh\gamma
+\log(1-|r_e|\tanh\gamma).                         \tag{3}
```

This improves strictly on the random-sign increment
`log cosh gamma+(1/2)log(1-r_e^2 tanh^2 gamma)` whenever
`0<|r_e|<1`. The gain is linear rather than quadratic in a small response.

Choose the initial orientation of the second internal block from `B,-B`.
Since the average of their *partition functions* at zero bridge is
`Z_A Z_B`, one orientation has

```math
\log Z_0\le\log Z_A+\log Z_B.                       \tag{4}
```

Therefore any reveal order followed by (2) constructs a complete signing
with

```math
\log Z_N\le\log Z_A+\log Z_B+\sum_{e\text{ bridge}}K_\gamma(r_e).  \tag{5}
```

There is no random-sign expectation or minimum/expectation interchange in
(5).

## 2. Exact calibrated sufficient inequality

Let

```math
T_{m,n}=(M_m^{2/3}+M_n^{2/3})^{3/2},\qquad
\Delta_{m,n}=T_{m,n}-M_m-M_n,                       \tag{6}
```

take `N=m+n`, `gamma_N=t_N/sqrt(N)`, and use exact minimizers for `A,B`.
The exact greedy response theorem sufficient for the desired recurrence is

```math
\boxed{
\sum_e K_{\gamma_N}(r_e)
\le\gamma_N\Delta_{m,n}+O(t_NN^{1-\delta}).}        \tag{7}
```

Together with the soft-to-ground entropy cost and `t_N>=N^delta`, (7) gives

```math
M_N\le T_{m,n}+O(N^{3/2-\delta}),\qquad
M_N^{2/3}\le M_m^{2/3}+M_n^{2/3}+O(N^{1-\delta}).   \tag{8}
```

Writing the greedy saving as

```math
D_\gamma(a)=-\log(1-a\tanh\gamma),                 \tag{9}
```

condition (7) is equivalently

```math
\sum_eD_{\gamma_N}(|r_e|)
\ge L\log\cosh\gamma_N-\gamma_N\Delta_{m,n}
    -O(t_NN^{1-\delta}),\qquad L=mn.                \tag{10}
```

Since `-log(1-z)>=z`, the simpler first-moment condition

```math
\tanh\gamma_N\sum_e|r_e|
\ge L\log\cosh\gamma_N-\gamma_N\Delta_{m,n}
    -O(t_NN^{1-\delta})                             \tag{11}
```

is sufficient. For `gamma_N=o(1)`, its leading threshold is

```math
\sum_e|r_e|
\ge {L\gamma_N\over2}-\Delta_{m,n}
    -O(N^{3/2-\delta}+L\gamma_N^3).                 \tag{12}
```

At balanced orders and growing `t_N`, this asks only

```math
{1\over L}\sum_e|r_e|\gtrsim {t_N\over2\sqrt N},  \tag{13}
```

up to the explicit child-optimum calibration. This is far weaker than the
`1-O(1/t_N)` squared polarization required by the random-sign route.

## 3. Why entropy and rank-one rate-distortion do not force (12)

At a fixed prefix, the needed potential is exactly the `l_1` norm of the
mean response matrix:

```math
\sum_{e\in U}|r_e|
=\max_{D\in\{\pm1\}^{U}}
  \mathbb E_\mu\!\left[\tau\sum_{(i,j)\in U}D_{ij}x_i y_j\right].  \tag{14}
```

Entropy and mutual information control this quantity in the wrong
direction. The uniform rank-one law has maximum entropy and every `r_e=0`.
Pinsker and correlation inequalities upper-bound individual `|r_e|` from
information; they do not give a lower bound. Low entropy is also
insufficient, since a two-cluster mixture can have small or zero mean
response.

The rank-one rate-distortion lemma proves that linear entropy forces many
`r_e` away from one. It permits all of them to be very small and hence gives
no part of (12). Thus the high-temperature nonpolarization theorem neither
proves nor falsifies the adaptive greedy target.

There is also no generic convex-gradient force. As a function of continuous
bridge couplings `J in [-1,1]^L`, `log Z(J)` is convex and its bridge gradient
vanishes at `J=0`, because the two fibers have independent global spin-flip
symmetries. Hence zero bridge is the continuous minimizer. A lower bound on
the gradients encountered while rounding it to a vertex must use special
algebra of the internal blocks; convexity alone supplies none.

## 4. The natural cubic potential and its exact limitation

There is one concrete structured signal. Let `S` be the current symmetric
interaction matrix, with candidate cross edge `e=(i,j)` still zero. Direct
expansion of

```math
r_e={\mathbb E_U[x_i y_j\sinh(\gamma H_S)]
          \over\mathbb E_U[\cosh(\gamma H_S)]}
```

shows that the linear term vanishes and the first formal numerator term is

```math
r_e=\gamma^3(S^3)_{ij}+\text{higher connected terms}. \tag{15}
```

The coefficient is exact: the surviving cubic monomials are precisely the
six orderings of length-three paths from `i` to `j`.

For conference internal blocks and a partial rectangular bridge `C`, the
full cubic cross block is

```math
A^2C+ACB+CB^2+CC^{\mathsf T}C.                    \tag{16}
```

At an unrevealed position `C_ij=0`, the two axial terms vanish, so its
cubic response is exactly `(ACB+CC^TC)_ij`. The linear part `ACB` has the
exact Frobenius potential

```math
\|ACB\|_F^2=(m-1)(n-1)\|C\|_F^2.                  \tag{17}
```

Thus after `k` signs are fixed, `ACB` has some entry of size at least
`Omega(sqrt(k))`. This is the first potential that points in the correct
direction, but it is not yet even an eligible-edge estimate: the large
entries could be concentrated on already revealed positions. A useful
deterministic lemma would have to lower-bound `ACB` on the zero support of
`C` and control its cancellation by `CC^TC`.

There is, however, an exact random-support theorem showing that neither
problem is generic at the cubic level. Give the entries of `C` the
independent law

```math
\Pr\{C_{ij}=1\}=\Pr\{C_{ij}=-1\}=\rho/2,
\qquad \Pr\{C_{ij}=0\}=1-\rho.                    \tag{18}
```

For every fixed candidate `(i,j)`, condition on `C_ij=0` and put

```math
Q_{ij}=(ACB+CC^{\mathsf T}C)_{ij}.
```

Independence and parity give the exact second moment

```math
\mathbb E[Q_{ij}^2\mid C_{ij}=0]
=(m-1)(n-1)(\rho+\rho^3).                         \tag{19}
```

Indeed, `ACB` is a signed sum of `(m-1)(n-1)` independent bridge
variables, so its variance is `(m-1)(n-1)rho`.  The cubic term is a sum
over the `(m-1)(n-1)` alternating length-three paths from `i` to `j`;
two path monomials have nonzero product expectation only when they are the
same path, giving variance `(m-1)(n-1)rho^3`.  A linear monomial times a
path monomial always leaves a variable of odd multiplicity, so the
covariance is zero.

For `rho` in any compact subinterval of `(0,1)`, degree-three
hypercontractivity for this fixed ternary product law gives a
dimension-free fourth-to-second moment bound.  Consequently

```math
\mathbb E[|Q_{ij}|\mid C_{ij}=0]
\ge c_\rho\sqrt{(m-1)(n-1)},                      \tag{20}
```

and summing over candidate positions gives

```math
\mathbb E\sum_{ij:C_{ij}=0}|Q_{ij}|=\Theta_\rho(N^3)
                                                               \tag{21}
```

at bounded aspect ratio.  Thus the cubic response has exactly the
aggregate scale suggested by the greedy target: formally,
`gamma^3 sum|Q_ij|=Theta_rho(t^3 N^(3/2))`.  This is large enough in scale
to compete with the `Theta(t N^(3/2))` threshold in (12).

It is not a theorem at the required scale. Even the exact cubic term can
cancel on a deterministic eligible support, while (18)--(21) concern an
independent partial bridge rather than the nested signs produced by the
greedy rule. More importantly, the higher connected expansion is uniform
only while `gamma||S||_op` is bounded strictly inside a high-temperature
domain. The recurrence needs `t_N=gamma sqrt(N)` to grow polynomially so
that the entropy error is summable. In that regime all path orders
contribute at leading scale and can cancel `(S^3)_ij`. Controlling those
cancellations is exactly a uniform Gibbs-response theorem, not a
consequence of (19)--(21).

## 5. Stopping judgment and precise surviving lemma

Adaptive greedy choice is a genuine improvement and supersedes the
near-polarization requirement of random bridge signs. The surviving target
is precisely (10), or the stronger but simpler cumulative first-moment bound
(11), for one reveal rule on exact minimizer or provably landing structured
blocks.

However, the available generic potentials do not force it:

- entropy and total correlation allow all absent-edge responses to be small;
- convexity starts at zero gradient;
- the random-support theorem gives the correct aggregate scale only for the
  first nonzero formal response term, outside the growing-temperature regime
  needed by (8).

Tracking the exact `r_e` at every prefix again retains the full Gibbs bridge
response. A new theorem summing the entire connected response into the lower
bound (11) would be genuinely useful and substantially weaker numerically
than full polarization, but no current compression proves it. This audit
therefore records (11) as the exact new lemma rather than claiming that the
greedy rule already gives a summable defect.
