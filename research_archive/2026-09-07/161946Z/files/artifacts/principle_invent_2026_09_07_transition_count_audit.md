# Independent audit: entropy-rich signed balanced incidences

Date: 2026-09-07. Status: **PASS conditional on the stated ordinary Eulerian-
orientation lower bound**. The director is retrieving that primary theorem.
The transition and local-order arguments below are reconstructed here and
do not rely on conventions in a quotation of the BEST theorem.

## 1. Statement and external input

Let m=2k+1, e=mk, and let S be any off-diagonal signing of the complete
graph K_m. A signed balanced incidence assignment chooses
`u_ij in {+-1}` such that

```math
u_{ij}u_{ji}=S_{ij},\qquad \sum_{j\ne i}u_{ij}=0\quad\text{for every i}.
```

The parity condition

```math
\prod_{i<j}S_{ij}=(-1)^e
\tag{1}
```

is necessary. Under (1), the proposed lower bound is

```math
\#\{\text{signed balanced assignments}\}
\ge\left[\frac{\binom{2k}{k}}{2^k k}\right]^m
=2^e\exp[-\tfrac32m\log k+O(m)].
\tag{2}
```

The only external input used by this audit is the ordinary Eulerian-
orientation inequality

```math
\operatorname{EO}(K_m)\ge
\left[\frac{\binom{2k}{k}}{2^k}\right]^m.
\tag{3}
```

The director attributes (3) to Schrijver and is checking its source.

## 2. Many single-circuit transition systems

A transition system pairs the 2k incident half-edges at each vertex. There
are

```math
T=((2k-1)!!)^m
```

systems. Alternating edge links and transition pairs partitions all edges
into undirected circuits. Let T_single count systems yielding one circuit.
Each single-circuit system supports exactly two ordinary Eulerian
orientations, given by its two directions.

Fix any ordinary Eulerian orientation D of K_m. It is strongly connected:
the condensation of a weakly connected balanced digraph cannot have an
edge between different strong components, since a source component would
have no incoming flow but a positive outgoing flow. Choose a directed
in-arborescence leading to a root r. Fix an outgoing edge e0 at r.

At every nonroot vertex, order its k outgoing edges with its arborescence
edge last. At the root, order its k outgoing edges with e0 first. There are
exactly `((k-1)!)^m` choices. Starting at r, repeatedly take the next unused
outgoing edge in the prescribed local order.

The walk can first become stuck only at r, by the equality of in- and
out-degrees. It cannot close before using every edge. Indeed the residual
graph after a closed trail remains balanced at every vertex. If a vertex
has a residual outgoing edge, its last-exit tree edge is also residual.
Following tree edges toward r propagates a residual incoming edge and
hence residual outgoing activity to r, contradicting that the walk stopped
there. Thus every local-order choice yields a full directed Euler circuit.

The full circuit recovers all local outgoing orders, so different choices
give different circuits. Its first edge is fixed as e0; therefore a single
transition system recovers the circuit uniquely by following that edge.
It follows that the number of D-compatible single transition systems is
at least `((k-1)!)^m`, without an unresolved root-degree factor.

Double-counting compatible orientations and single transition systems and
using (3) gives

```math
2T_{\rm single}\ge\operatorname{EO}(K_m)((k-1)!)^m,
\qquad
T_{\rm single}\ge\frac{T}{2k^m}.
\tag{4}
```

The second identity is exact algebra:
`[binom(2k,k)/2^k](k-1)! = (2k-1)!!/k`.

## 3. Signed propagation and the second double count

For a fixed single transition system, prescribe that incidence signs are
opposite across each transition pair and have product S_ij across each
edge. Propagating from one half-edge around the unique alternating circuit
returns consistently exactly when

```math
(-1)^e\prod_{i<j}S_{ij}=1.
```

Under (1), there are exactly two propagated assignments. They are balanced
because every vertex has k opposite-sign transition pairs.

Conversely, for any signed balanced assignment, the number of compatible
transition systems (not necessarily single-circuit) is exactly `(k!)^m`:
at each vertex match the k positive incidences bijectively to the k negative
incidences. Consequently

```math
\#\{\text{signed balanced assignments}\}(k!)^m
\ge 2T_{\rm single}.
```

Insert (4) and simplify to obtain (2). This part uses no property of S
beyond the one global parity condition. In particular no spectral or cap
bound is needed.

## 4. Consequence for the binary Gaussian kernel at any temperature

Let each row be a uniform permutation of k copies of +B and k copies of
-B. Its total assignment count is `binom(2k,k)^m`. For

```math
Z_S(t)=\mathbb E\prod_{i<j}
\exp[-t(u_{ij}-S_{ij}u_{ji})^2],
```

the count (2) gives the uniform lower bound

```math
\log Z_S(t)\ge-e\log2-m\log k
\tag{5}
```

for parity-compatible S, at EVERY t>=0, including the hard-constraint
limit. If parity fails, flip one edge of S to obtain a compatible S'. All
the assignments counted for S' incur exactly one mismatch for S. Hence

```math
\log Z_S(t)\ge-e\log2-m\log k-4tB^2.
\tag{6}
```

The symmetric canonical pressure is exactly

```math
e\log\left(\frac{1+e^{-4tB^2}}2\right).
```

Together with the exact-type upper identity in the symmetric-type note,
(5)--(6) show that the binary model has no order-m² seed distinction at
the proposed critical scale t~sqrt(m). More strongly, parity-compatible
seeds remain seed-universal at leading entropy even for an exact hard
constraint; for general seeds the explicit parity payment is only 4tB².

These conclusions concern a kernel partition, not the cap of every actual
weave or the convergence of the original extremal sequence.
