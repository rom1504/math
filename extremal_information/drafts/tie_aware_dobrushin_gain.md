# Tie-aware tropical secants and exact response gain

**Status.** Proof source for Theorem 16.18. Exact rational checks are in
[`../experiments/verify_dobrushin_response_gain.py`](../experiments/verify_dobrushin_response_gain.py).

This note extends fixed-selector stability to switching and ties without
choosing a canonical optimizer. The tangent object is a row-stochastic
secant, and suffix-row separation gives an exact adversarial response gain.

## 1. Exact projective gain

Fix `r>=2` and put

```math
V=\mathbb R^r/\mathbb R\mathbf 1,
\qquad \|[v]\|_H={\operatorname {osc}(v)\over2}.
```

Every row-stochastic matrix `P` descends to a nonexpansive map on `V`. For
probability rows `p,q`, write

```math
\operatorname {TV}(p,q)={1\over2}\sum_k|p_k-q_k|,
\qquad \delta(P)=\max_{i,j}\operatorname {TV}(P_{i,*},P_{j,*}).
```

Consider

```math
e_t=P_te_{t-1}+\eta_t,
\qquad e_0=0,
\qquad \|\eta_t\|_H\le\epsilon_t.                 \tag{1.1}
```

Put `R_(s,T)=P_T...P_(s+1)`, with `R_(T,T)=I`.

### Theorem 1 (exact suffix-row gain)

```math
\boxed{
\sup_{\|\eta_s\|_H\le\epsilon_s}\|e_T\|_H
=\max_{i,j}\sum_{s=1}^T\epsilon_s
 \operatorname {TV}(R_(s,T)_(i,*),R_(s,T)_(j,*)) .}             \tag{1.2}
```

#### Proof

Unroll (1.1). Since

```math
\|e_T\|_H=\max_{i,j}{(e_T)_i-(e_T)_j\over2},
```

the supremum and finite maximum may be interchanged. For one ordered pair,
the residuals at distinct times optimize independently. If `d=p-q` is the
difference of two probability rows, then `sum d_k=0` and

```math
\sup_{\|v\|_H\le\epsilon}d\mathbin\cdot v
=\epsilon\|d\|_1=2\epsilon\operatorname {TV}(p,q).
```

The outer factor `1/2` leaves `epsilon TV`, proving (1.2). `square`

The optimizer may depend on the terminal row pair and every suffix. Thus
(1.2) is exact for fresh adversarial residuals, but only an upper bound for
fixed coherent nonlinear families whose secants and residuals share an
orbit.

## 2. Mixing versus finite reset

If every contiguous block of length `L` occurring in a legal trajectory has
Dobrushin coefficient at most `rho<1`, submultiplicativity and (1.2) give

```math
\|e_T\|_H\le {L\over1-\rho}\max_t\epsilon_t.       \tag{2.1}
```

This is broader than a tangent reset. Repeated use of

```math
P_\alpha=\begin{pmatrix}1-\alpha&\alpha\\
                         \alpha&1-\alpha\end{pmatrix},
\qquad0<\alpha<1,
```

has gain at most `epsilon/(1-|1-2alpha|)`, although no finite power has
identical rows unless `alpha=1/2`.

There is also a finite-semigroup converse. Let all transports belong to a
finite row-stochastic semigroup. Whenever it has a nonconsensus element, set

```math
\gamma=\min\{\operatorname {TV}(R_(i,*),R_(j,*))>0:
             R\in S,\ i\ne j\}>0.                 \tag{2.2}
```

For a factorial language, uniform fresh-residual gain is equivalent to a
bounded gap between contiguous products having identical rows. A reset-free
word of length `T` has every terminal suffix nonconsensus. Summing their
Dobrushin coefficients and pigeonholing the maximizing ordered pair gives

```math
\|e_T\|_H\ge {\gamma T\over r(r-1)}\epsilon.       \tag{2.3}
```

Conversely, an identical-row factor is a two-sided ideal under stochastic
multiplication and erases all older residuals. The selector theorem is the
case `gamma=1`. An infinite stochastic semigroup need not contain a reset:
geometric suffix-row mixing is the additional mechanism.

## 3. Max-plus switches and ties

For an all-finite max-plus map

```math
(F_Kx)_j=\max_i\{x_i+K_(ij)\},                    \tag{3.1}
```

and any `x,y`, there is a row-stochastic `P_K[x,y]` such that

```math
F_Ky-F_Kx=P_K[x,y](y-x).                           \tag{3.2}
```

Each row can be supported on indices active somewhere on the segment from
`x` to `y`: restrict one output maximum to the segment and integrate its
piecewise-constant slope. This includes optimizer switches and tie faces.

For trajectories `x_t=F_t(x_(t-1))`, `y_t=G_t(y_(t-1))`, decompose

```math
y_t-x_t=[G_t(y_(t-1))-G_t(x_(t-1))]
        +[G_t(x_(t-1))-F_t(x_(t-1))].              \tag{3.3}
```

The first bracket has form (3.2); if the same-input discrepancy in the
second is at most `epsilon_t`, Theorem 1 applies. Uniform scrambling of all
realized secants is therefore a depth-stable switching theorem.

Three limitations are essential. Merely crossing a tie gives no uniform
`rho<1`, since a secant weight may approach zero. The lower bound permits
fresh residuals and is not a coherent-kernel converse. Finally, recognizing
all dynamically realizable secants can still require paired-cell
reachability; allowing every stochastic row supported on a face is a safe
upper relaxation but may create false instability cycles.

The next finite target is thus precise: certify block scrambling for all
realizable paired-cell secants, or realize one nonmixing cycle together with
a coherent residual of nonzero transported mean.
