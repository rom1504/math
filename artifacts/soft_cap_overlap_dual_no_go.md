# Integrated overlap dual: exact compression and no-go

Status: final bounded abstraction for the soft-cap route. There is an exact
scalar compression of the reveal gradients, but its required lower bound is
a single-cluster theorem for the full Gibbs bridge response. No generic
entropy, total-correlation, or convex-dual argument supplies it.

Let a bridge have `L=mn` edges and reveal its uniformly random order. At
prefix `s`, let `U_s` be the unrevealed set and

```math
r_{s,e}=\mathbb E_{\mu_s}[\tau x_i y_j],\qquad e=(i,j)\in U_s.
```

The exact next-edge cost is

```math
I_\gamma(r)=\log\cosh\gamma
+\frac12\log(1-r^2\tanh^2\gamma).                  \tag{1}
```

## 1. Exact scalar compression

Define the unrevealed squared-overlap statistic

```math
q_s={1\over|U_s|}\sum_{e\in U_s}r_{s,e}^2.          \tag{2}
```

Conditional on the prefix, the next edge is uniform in `U_s`. Since
`u -> log(1-u tanh^2 gamma)` is concave,

```math
\mathbb E[I_\gamma(r_{s,e_{s+1}})\mid\mu_s,U_s]
\le J_\gamma(q_s),                                  \tag{3}
```

where

```math
J_\gamma(q)=\log\cosh\gamma
+\frac12\log(1-q\tanh^2\gamma).                   \tag{4}
```

Concavity a second time, over prefixes and reveal histories, gives

```math
\mathbb E\sum_{s=0}^{L-1}I_\gamma(r_{s,e_{s+1}})
\le L J_\gamma(\overline q),                        \tag{5}
```

with the single scalar

```math
\overline q={1\over L}\sum_{s=0}^{L-1}\mathbb E q_s. \tag{6}
```

Thus the `Theta(N^2)` adaptive gradients really can be compressed for the
purpose of an upper bound on total bridge cost.

## 2. Two-replica and convex-dual meanings

Let `(x,y,tau)` and `(x',y',tau')` be independent samples from `mu_s`. Then

```math
q_s=\mathbb E_{\mu_s^{\otimes2}}
\left[{1\over|U_s|}\sum_{(i,j)\in U_s}
 \tau\tau' x_ix_i'y_jy_j'\right].                  \tag{7}
```

For the complete bridge this is the product of the two fiber overlaps,
multiplied by `tau tau'`. For a partial bridge it is the corresponding
restricted Gram overlap. Equation (7) is an exact two-replica scalar, not a
heuristic replica ansatz.

There is also an exact convex dual. For real matrices `Z` supported on
`U_s`,

```math
q_s={1\over|U_s|}\max_Z
\left\{2\mathbb E_{\mu_s}\!\left[
 \tau\sum_{(i,j)\in U_s}Z_{ij}x_i y_j\right]
-\sum_{(i,j)\in U_s}Z_{ij}^2\right\}.              \tag{8}
```

This is just `||E_mu O||_F^2=max_Z(2<E_mu O,Z>-||Z||_F^2)` for the rank-one
observable matrix `O_ij=tau x_i y_j`. The unique optimizer is the complete
mean-response matrix `Z_ij=r_{s,ij}`.

## 3. Exact calibrated scalar theorem that would suffice

Retain

```math
T_{m,n}=(M_m^{2/3}+M_n^{2/3})^{3/2},\qquad
\Delta_{m,n}=T_{m,n}-M_m-M_n,                       \tag{9}
```

and let `gamma_N=t_N/sqrt(N)`. By (5), the calibrated bridge theorem from
the preceding audit follows from the one-scalar inequality

```math
LJ_{\gamma_N}(\overline q)
\le\gamma_N\Delta_{m,n}+O(t_NN^{1-\delta}).         \tag{10}
```

Solving (10) for the overlap gives the exact sufficient threshold

```math
\overline q\ge
{1-\exp\!\left(2\gamma_N\Delta_{m,n}/L
                 +O(t_NN^{1-\delta}/L)\right)
       /\cosh^2\gamma_N
 \over \tanh^2\gamma_N}.                           \tag{11}
```

When `gamma_N=o(1)`, this is

```math
1-\overline q
\le {2\Delta_{m,n}\over\gamma_NL}
+\text{allowed lower-order error}.                 \tag{12}
```

At balanced orders the right side is `Theta(1/t_N)`. Hence a polynomially
growing soft parameter requires

```math
\overline q=1-O(1/t_N)                              \tag{13}
```

with the leading constant calibrated by the two exact child optima. Equations
(5) and (11) show that this overlap theorem would genuinely imply the
summable `b_n` recurrence; no individual gradient needs to be tracked after
(11) is available.

## 4. Why standard scalar information does not prove (13)

The statistic has a direct metric meaning. If `O,O'` are two Gibbs bridge
responses, then

```math
L(1-q_s)=2\mathbb E[d_H(O,O')]                     \tag{14}
```

on the relevant unrevealed coordinates. Thus (13) says two independent
Gibbs samples differ on only `O(L/t_N)` bridge entries on average. By the
rank-one row/column-flip geometry, they lie near one common response cluster
up to the four-element kernel.

Entropy gives only the opposite useful implication: linear entropy forces
`1-q_s` bounded below. Low entropy does not force `q_s` near one. For
example, a two-point law on two far-apart rank-one matrices has entropy only
`log2` but can have a very small mean-response norm. Total correlation has
the same defect: it measures concentration dimension, not whether the
remaining clusters are mutually aligned.

The Gram formulation (7) does not fix this. A lower bound near one is a
replica-concentration or single-pure-state theorem. Knowledge of the overlap
spectrum would prove it only by already including the mass-near-one statement
that is needed.

## 5. The convex dual re-expands to full response

Although (8) is formally a minimax compression, its unrestricted optimizer
has one coordinate for every unrevealed bridge edge and is exactly the
unknown response matrix. Restricting `Z` to row sums, spectra, a bounded
number of fibers, or any other fixed-dimensional family gives a valid lower
bound on `q_s`, but there is no uniform reason for that lower bound to
approach one. A successful restriction would itself be the missing
bounded-complexity state theorem.

As `t_N` grows, the Gibbs gradient approaches the convex hull of rank-one
near-ground bridge responses. Condition (13) says that this convex hull has
Frobenius diameter `o(sqrt(L))` along almost the whole reveal path. Determining
that diameter requires knowing whether the near-maximizing cuts align under
every partial bridge. This is the soft analogue of the full adaptive
bridge/cap response, not information contained in the scalar child cap.

After comparison with the ledger, this matches the augmented cut-response
and full-bridge obstructions already found there: changing to the dual matrix
`Z` does not reduce its dimension. The integrated overlap `overline q` is a
real exact summary **if independently controlled**, but no generic convex or
information inequality controls it in the required direction.

## 6. Stopping judgment

The only surviving theorem is now precise: prove (11) from a genuinely
lower-complexity structural invariant of exact minimizers. Entropy,
mutual information, the overlap spectrum, and the unrestricted Frobenius
dual either point in the wrong direction or restate the complete Gibbs
response. No genuine compression presently survives, so this abstraction
stops here rather than relabeling (11) as progress.
