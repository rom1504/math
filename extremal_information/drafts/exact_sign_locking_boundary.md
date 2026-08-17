# Exact sparse contrast compilation and the one-layer locking obstruction

Status: rigorous task-local draft.  This separates two facts that are easy
to conflate after the short-seed Gram-broadcast theorem:

1. every same-support quadratic contrast has an exact **sparse** unit-edge
   disjoint lift;
2. the most economical complete-sign locking idea--duplicate the spins and
   use one dense bipartite sign bridge--has an unavoidable leading
   desynchronization incentive.

Neither statement rules out a multilayer or child-specific exact-sign
compiler.  Together they identify the precise scale and structural gap.

## 1. An exact sparse edge-variable compiler

Let

```math
T(x)=\sum_(i<j)t_(ij)x_ix_j,\qquad t_(ij)\in\{+-1\},
```

and introduce one auxiliary spin `y_(ij)` for each edge.  Put

```math
G_T(x,y)=\sum_(i<j)y_(ij)(x_i-t_(ij)x_j).             \tag{EL.1}
```

Every nonzero coefficient of `G_T` is exactly `+-1`, and the old and new
variables are disjoint ports.

### Proposition EL.1 (universal sparse contrast lift)

For every function `F:{+-1}^k->R`,

```math
\max_(x,y)\{F(x)+G_T(x,y)\}
={k\choose2}+\max_x\{F(x)-T(x)\}.                    \tag{EL.2}
```

Hence two orientations, applied to `(F,T)` and `(-F,-T)`, recover
`max_x|F(x)-T(x)|` exactly.  If an existing child cannot itself be
sign-switched, this is not a single-future absolute-response compiler.

#### Proof

For fixed `x`, the auxiliary spins separate and

```math
\max_(y_(ij)=+-1)y_(ij)(x_i-t_(ij)x_j)
=|x_i-t_(ij)x_j|=1-t_(ij)x_ix_j.                    \tag{EL.3}
```

Summing (EL.3), adding `F`, and maximizing `x` proves (EL.2). `square`

This is an exact unit-edge quadratic compilation, not an exact signing on a
complete graph.  It uses

```math
N=k+{k\choose2}=Theta(k^2)                           \tag{EL.4}
```

vertices and leaves all auxiliary--auxiliary and most auxiliary--old pairs
absent.  A `Theta(k^(3/2))` contrast gap therefore becomes only
`Theta(N^(3/4))`, not `Theta(N^(3/2))`, in the enlarged order.  Completing
the missing pairs by arbitrary signs introduces a new dense Hamiltonian and
is not a harmless padding operation.

## 2. A dense bipartite bridge cannot universally lock duplicated spins

The tempting linear-overhead compiler duplicates `x` as `z` and tries to
force `z=D x`, for a fixed diagonal sign gauge `D`, by a complete bipartite
sign bridge

```math
L_R(x,z)=x^TRz,\qquad R\in\{+-1\}^{k\times k}.       \tag{EL.5}
```

Let

```math
Delta_(R,D)(x)=\max_zL_R(x,z)-L_R(x,Dx)
=||R^Tx||_1-x^TRDx.                                  \tag{EL.6}
```

It is the energetic incentive to abandon the prescribed synchronized copy.

### Theorem EL.2 (one-layer complete-sign locking has leading defect)

For every `R in {+-1}^{k times k}` and every diagonal sign matrix `D`,

```math
\max_x Delta_(R,D)(x)
\ge k\,E|S_k|-k
\ge {1\over\sqrt2}k^(3/2)-k,                        \tag{EL.7}
```

where `S_k` is a sum of `k` independent Rademacher signs.  In particular no
single complete unit-sign bipartite layer makes a prescribed coordinatewise
copy even `o(k^(3/2))`-optimal uniformly over all source spins.

The exact pointwise condition is already impossible for `k>=3`: for every
column, an adversarial choice of the other `k-1` source signs can make the
prescribed target coordinate disagree with the bridge field.

#### Proof

For uniform `x`, every column dot product `(R^Tx)_j` has the distribution of
`S_k`, so

```math
E_x||R^Tx||_1=kE|S_k|.                               \tag{EL.8}
```

Meanwhile

```math
E_x x^TRDx=tr(RD),\qquad |tr(RD)|<=k.                \tag{EL.9}
```

Averaging (EL.6) and using the sharp real `p=1` Khintchine lower bound
`E|S_k|>=sqrt(k/2)` proves (EL.7).

For the pointwise claim, fix target coordinate `j`.  The desired sign is
`d_jx_j`, while its bridge field is `sum_iR_(ij)x_i`.  After fixing `x_j`,
choose every other `x_i` to make `d_jx_jR_(ij)x_i=-1`.  Then the desired-sign
field is at most `1-(k-1)<0` for `k>=3`, regardless of the diagonal sign.
`square`

### Corollary EL.3 (universal duplicate-block compilation fails at scale)

Suppose a proposed universal compiler places an arbitrary child landscape on
`x`, a query landscape on a duplicated block `z`, and relies only on one
complete unit-sign bipartite bridge to impose `z=Dx`.  Then there is a pinned
source configuration for which freeing `z` gains at least

```math
{1\over\sqrt2}k^(3/2)-k.                             \tag{EL.10}
```

Thus the bridge itself injects an error at exactly the scale the Gram
response packing is trying to preserve.  A valid compiler must use a
multilayer/code constraint, exploit a restricted source set, depend jointly
on child and query, or abandon complete unit signs.

The source pin in this corollary is a legitimate universal-context test; it
does not say that the bad `x` is an optimizer of every particular flat child.
Therefore EL.2 is a scalable no-go for the one-layer **universal** architecture,
not for all possible exact-sign realizations of the alternating-form family.

## 3. Information accounting

The query `T` in EL.1 may itself have a concise parameterization--for the
alternating-form family it costs only `h=Theta(k)` bits on top of the shared
short seed.  The obstacle is not enumerating an exponential language.  It is
realizing its action in the complete exact-sign class without either
quadratic auxiliary order or a fresh `Theta(k^(3/2))` synchronization error.

Finite exhaustive checks are in
[`../experiments/verify_exact_sign_locking_boundary.py`](../experiments/verify_exact_sign_locking_boundary.py).
