# Bounded-cap switching orbits retain large finite bridge responses

Status: one exact deterministic exposure lemma and a reproducible finite
experiment.  The scaling statement below is a conjecture, not a theorem.

## 1. Falsifiable scaling conjecture

Let

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
(P_BH_A)(y)=\max_x\{H_A(x)+x^TBy\}.
```

The experiment was fixed before its output was inspected to test:

> **Bounded-cap orbit-packing conjecture.** There are constants
> `epsilon,c,C>0`, bounded-cap hollow sign matrices `A_n` with
> `Q(A_n)<=C n^(3/2)`, and dense sign bridges `B_n`, for which the switching
> orbit
>
> ```math
> A_n^s=D_sA_nD_s,
> \qquad s\in\{-1,1\}^n/\{\pm1\},
> ```
>
> contains `exp(cn)` children whose transforms `P_{B_n}H_{A_n^s}` are
> pairwise separated by `epsilon n^(3/2)` in projective sup distance.

One seeded finite sequence cannot establish this.  The point of the probe is
to try to kill the conjecture cheaply: normalized gaps could have collapsed,
or the switching responses could have fallen into a small number of classes.
Neither happened through order 14.

The conjecture is intentionally existential.  The stronger claims “every
bounded-cap sequence works” and “a typical bridge works with a uniform
constant” have not been formulated as theorems and are not supported by the
finite data alone.

## 2. A deterministic query-linked exposure identity

Write

```math
F_A(h)=\max_u\{H_A(u)+h^Tu\},
\qquad P=\max_u H_A(u).
```

Choose a top state `u_*`, so `H_A(u_*)=P`, and define the nonnegative
**weighted top deficit**

```math
\Delta_A(h)=P+\lVert h\rVert_1-F_A(h).
```

For a convention `sign(0)=1`, it has the exact form

```math
\Delta_A(h)
=\min_u\left\{
P-H_A(u)
+2\sum_{i:u_i\ne \operatorname{sign}(h_i)}|h_i|
\right\}.                                                   \tag{BC.1}
```

### Lemma BC.1 (query-linked switching exposure)

For every query `y`, put

```math
s_y=u_*\odot\operatorname{sign}(By).
```

Then

```math
(P_BH_{A^{s_y}})(y)=P+\lVert By\rVert_1.                    \tag{BC.2}
```

For two queries `y,z`, let `R_y=P_BH_{A^{s_y}}`.  Their projective response
distance obeys

```math
d_{\rm proj}(R_y,R_z)
\ge {1\over2}\left[
\Delta_A(s_z\odot By)+\Delta_A(s_y\odot Bz)
\right].                                                    \tag{BC.3}
```

#### Proof

Changing variables `u=s\odot x` gives

```math
(P_BH_{A^s})(y)=F_A(s\odot By).                             \tag{BC.4}
```

The identity

```math
h^Tu=\lVert h\rVert_1
-2\sum_{i:u_i\ne\operatorname{sign}(h_i)}|h_i|
```

proves (BC.1).  For `s=s_y`, the field in (BC.4) is
`u_*\odot|By|`.  The same `u_*` simultaneously maximizes the quadratic and
linear terms, proving (BC.2).

At the query `y`, (BC.1)--(BC.2) give

```math
R_y(y)-R_z(y)=\Delta_A(s_z\odot By)\ge0.
```

At `z` the response difference is
`-Delta_A(s_y odot Bz)`.  The oscillation of `R_y-R_z` is therefore at least
the sum of these two quantities, proving (BC.3). `square`

The lemma isolates a concrete asymptotic obligation.  A query-linked packing
follows if random bridge fields avoid, in weighted Hamming distance, the
near-top set of `H_A`.  More precisely, (BC.1) immediately gives:

### Corollary BC.1a (near-top neighborhood criterion)

If, for a number `eta>0`,

```math
2\min_{u:P-H_A(u)<\eta}
\sum_{i:u_i\ne\operatorname{sign}(h_i)}|h_i|\ge\eta,
```

then `Delta_A(h)>=eta`.

This criterion is not a proof of the scaling conjecture: controlling all
pairwise linked fields for an exponentially large family remains open.  It
does identify the missing information more narrowly than the full response
landscape: weighted neighborhoods of the near-top set, rather than arbitrary
quadratic coefficients or every response value.

## 3. Exact finite experiment

The program
[`bounded_cap_switching_response.py`](../experiments/bounded_cap_switching_response.py)
loads saved exact minimizers, conference matrices, and one saved-cap heuristic.
It uses:

- a deterministic seeded Rademacher bridge at every tested order;
- the dense sign bridge `A+I` at every order;
- a Sylvester Hadamard bridge at powers of two;
- exact integer maximization over every projective Boolean child spin;
- the projective metric
  `osc(P_BH_{A^s}-P_BH_{A^t})/2`;
- greedy packings, explicitly labelled as certified lower bounds rather than
  maximum packings.

For the seeded Rademacher bridge, the complete switching orbit and all Boolean
queries were enumerated through order 10.  The exact results are:

| `n` | `Q(A)/n^(3/2)` | orbit size | distinct responses | minimum normalized distance | median normalized distance | packing at `0.1 n^(3/2)` |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 0.5000 | 8 | 6 | 0.0000 | 0.3750 | 6 |
| 6 | 0.3402 | 32 | 32 | 0.1361 | 0.3402 | 32 |
| 8 | 0.4419 | 128 | 128 | 0.1326 | 0.3536 | 128 |
| 10 | 0.4743 | 512 | 512 | 0.1581 | 0.4111 | 512 |

Thus, at orders 6, 8, and 10, the entire `2^(n-1)` switching orbit is an
exact `0.1 n^(3/2)` response packing for this fixed bridge.

For larger orders the program freezes a query-linked family of size
`2^floor(n/2)` and evaluates its responses on all queries in that same family.
Those restricted-query distances are rigorous lower bounds on the full
future-response metric:

| `n` | source type | family | minimum normalized distance | median normalized distance | packing at `0.1 n^(3/2)` |
|---:|---|---:|---:|---:|---:|
| 4 | exact minimizer | 4 | 0.3750 | 0.4375 | 4 |
| 6 | exact minimizer | 8 | 0.0680 | 0.2722 | 7 |
| 8 | exact minimizer | 16 | 0.0884 | 0.2652 | 14 |
| 10 | conference | 32 | 0.1265 | 0.2846 | 32 |
| 12 | saved-cap heuristic | 64 | 0.1443 | 0.3368 | 64 |
| 14 | conference | 128 | 0.1336 | 0.2864 | 128 |

The `A+I` bridge also keeps all 32, 64, and 128 displayed states separated at
the `0.1 n^(3/2)` scale for `n=10,12,14`.  The Hadamard bridge keeps all 16
displayed states at `n=8`.  These are robustness checks across three bridge
types, not independent asymptotic samples.

Every table entry is reproduced in
[`bounded_cap_switching_response.json`](../experiments/results/bounded_cap_switching_response.json).
The output records source hashes, bridge hashes, operator norms, all distance
quantiles, and packing sizes at eight fixed resolutions.

## 4. What the experiment does and does not show

The finite probe decisively falsifies the cheap hypothesis that a small
Boolean cap alone makes switching-orbit bridge responses collapse at small
orders.  Even matrices at the natural `n^(3/2)` cap scale can expose every
switching bit under a dense bridge.

It does **not** prove an `Omega(n)` asymptotic information lower bound.  The
orders are small, the matrices come from a few structured families, and the
larger-order packings deliberately test only `2^(n/2)` states.  A proof needs
a uniform lower bound on (BC.1), or another exposure certificate, for an
exponentially large code of linked fields.

There is also no contradiction with low-cap optimization.  All children are
switchings of one matrix, so they have identical isolated cap.  The bridge
reveals the hidden gauge.  Consequently any eventual compression theorem for
near-minimizers must either quotient this gauge jointly with the bridge or
prove that the relevant composition never exposes it.  “Bounded cap” by
itself is not an empirically credible compression hypothesis.

## 5. Next falsifiable theorem

The cleanest next statement is a weighted-neighborhood entropy theorem for a
conference sequence:

> For a symmetric conference signing `A_n`, a random dense sign bridge admits
> `exp(cn)` query-linked fields such that every pair satisfies the two-sided
> deficit lower bound in (BC.3) with `eta=epsilon n^(3/2)`.

The conference identity gives a sharp Hanson--Wright tail for the energy of a
random spin.  What is not yet controlled is the union over a weighted Hamming
ball around each bridge-field sign.  Proving or disproving that ball estimate
would distinguish genuine linear response information from a small-order
artifact without reopening the original convergence problem.
