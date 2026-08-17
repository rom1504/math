# Connected packing of independent Walsh triangle-flux bits

Status: task-local rigorous draft.  This combines the audited scalar
triangle of Theorem 21.21 with an exact nonnegative bridge-padding lemma.
It gives `2^h` unrooted scalar-response states separated at the one-block
`n^(3/2)` scale, using only `h` queries.  Every query has the same single
connected support.  There are both maximum-degree-four and complete-graph
versions.  No rooted field, pinned spin, or coordinate-dependent query is
used.

## 1. An `h`-cube inside the off-diagonal Gram form

Fix `h>=1` and put

```math
m=4h+1,\qquad q=2^m,\qquad n=q^2.
```

Let `F_E` be the normalized Walsh involution on
`E=F_2^m direct-sum F_2^m`, and, for `a in F_2^m`, put

```math
C_a=D_{(0,a)}F_ED_{(0,a)}.                              \tag{CF.1}
```

In the `i`th four-coordinate chunk write

```math
\begin{aligned}
u_i&=e_{4i-3}+e_{4i-2},\\
v_i^0&=e_{4i-1}+e_{4i},
&w_i^0&=u_i+v_i^0,\\
v_i^1&=e_{4i-3}+e_{4i-1},
&w_i^1&=u_i+v_i^1.
\end{aligned}                                           \tag{CF.2}
```

For `sigma in {0,1}^h`, the marked tuple `a^sigma` has the three labels

```math
(a_{i1}^sigma,a_{i2}^sigma,a_{i3}^sigma)
=(u_i,v_i^{sigma_i},w_i^{sigma_i})                     \tag{CF.3}
```

in gadget `i`.  Its non-Gram data are independent of `sigma`:

```math
\begin{aligned}
a_{ir}^sigma\mathbin\cdot a_{ir}^sigma&=0,\\
\mathcal R_h
 &=\{c in F_2^{3h}:c_{i1}=c_{i2}=c_{i3}\text{ for every }i\},\\
\mathcal R_omega(a^sigma)&=\varnothing.
\end{aligned}                                           \tag{CF.4}
```

Indeed, each local triple has rank two and its sole nonzero relation is
`111`.  Different chunks are linearly disjoint and orthogonal.  The last
ambient coordinate is unused, so the characteristic vector cannot lie in
the span.  All cross-gadget Gram entries vanish.  Within gadget `i`,

```math
a_{ir}^sigma\mathbin\cdot a_{is}^sigma
=\begin{cases}
0,&r=s\text{ or }sigma_i=0,\\
1,&r\ne s\text{ and }sigma_i=1.
\end{cases}                                             \tag{CF.5}
```

Thus, after fixing every self-pairing, the complete relation kernel, the
complete characteristic-root fibre, and every cross-gadget pairing, the
remaining Gram cube has exactly one free bilinear flux bit per relation
triangle.  The three off-diagonal entries in one triangle are equal by
bilinearity and the relation, so (CF.5) represents `h` independent bits,
not `3h` independent entries.

## 2. The local scalar gap

For gadget `j`, let

```math
T_j^sigma(X)
 ={q\over2}\sum_{r=1}^3x_{jr}^TC_{a_{jr}^sigma}x_{jr}
 +q\sum_{1\le r<s\le3}x_{jr}^TF_Ex_{js}.              \tag{CF.6}
```

Theorem 21.21 gives, with

```math
M_0={9\over2},\qquad
M_1={3(1+\sqrt {17})\over4},\qquad
\delta=M_0-M_1={3(5-\sqrt {17})\over4},               \tag{CF.7}
```

the bounds

```math
\max T_j^sigma=M_0n^{3/2}\quad(sigma_j=0),
\qquad
\max T_j^sigma\le M_1n^{3/2}\quad(sigma_j=1).         \tag{CF.8}
```

The first equality has a Boolean witness `x` satisfying

```math
F_Ex=x,\qquad C_{a_{jr}^sigma}x=x\quad(r=1,2,3).       \tag{CF.9}
```

The construction in Theorem 21.21 applies after permuting its first four
coordinates into the `j`th chunk.

## 3. Nonnegative bridge padding does not erode the gap

Let `H` be any nonnegative weighted graph on all `3h` blocks, and let

```math
B_H=\sum_{e in E(H)}w_e.
```

The graph may have arbitrary size and may meet the three target vertices.
Define the padded target query

```math
Q_{i,H}^{sigma,0}(X)
=T_i^sigma(X)+q\sum_{uv in E(H)}w_{uv}x_u^TF_Ex_v.     \tag{CF.10}
```

### Lemma CF.1 (exact one-sided connective padding)

For every `i` and every `H` as above,

```math
\begin{aligned}
\max Q_{i,H}^{sigma,0}&=(M_0+B_H)n^{3/2}
 &&\text{if }sigma_i=0,\\
\max Q_{i,H}^{sigma,0}&\le(M_1+B_H)n^{3/2}
 &&\text{if }sigma_i=1.
\end{aligned}                                          \tag{CF.11}
```

#### Proof

For Boolean vectors, orthogonality of `F_E` gives

```math
|q x^TF_Ey|\le qn=n^{3/2}.                             \tag{CF.12}
```

Therefore the connector contributes at most `B_Hn^(3/2)` in either state.
This and (CF.8) prove both upper bounds in (CF.11).  In the zero-flux state,
put the common Boolean vector from (CF.9) in **every** block.  It attains
every nonnegative connector bound simultaneously, including connectors
far outside the target triangle.  It also attains all six target terms.
This proves equality. `square`

The point is that connectivity is not being purchased with vanishing edge
weights.  Unit path edges, or all unit cross-gadget edges, preserve the
full local gap `delta` exactly.  The connector mass can grow like `h` or
`h^2`; it contributes a public baseline to the favorable state and no more
than the same ceiling to the unfavorable state.

## 4. One common connected support, with every label active

To make every marked label occur in a nonzero child term, add a small total
mass of the other local triangles.  For `h>=2`, set

```math
gamma_h={1\over100(h-1)},                               \tag{CF.13}
```

and set `gamma_1=0`.  Query `i` is

```math
Q_i^sigma
=Q_{i,H}^{sigma,0}+\gamma_h\sum_{j\ne i}T_j^sigma.    \tag{CF.14}
```

This is one ordinary weighted unrooted Walsh graph.  Parallel terms in
(CF.14), if any, are combined by adding their real edge weights.

There are two useful common choices of `H`.

1. **Bounded degree.**  For `h>=2`, order the vertices
   `(1,1),(2,1),...,(h,1),(1,2),...,(h,3)` and let `H` be the unit
   Hamiltonian path in this order.  Every path edge joins different
   gadgets.  The union with the `h` local triangles is connected and has
   maximum degree at most four.  All onsite and edge weights lie in
   `(0,1]`.  For `h=1`, take `H` empty; the target triangle is connected.
2. **Dense.**  Let `H` contain every cross-gadget pair with unit weight.
   After the local triangles are added, the query graph is `K_(3h)`.
   Every cross-gadget edge has weight one, target internal edges have weight
   one, and non-target internal edges have weight `gamma_h`.

In both cases all `h` queries have exactly the same connected support,
namely all `3h` marked blocks.  For `h>=2`, every block also has a positive
onsite weight, so no label is present merely as a bridge-only port.

For every Boolean assignment, one local triangle has the uniform absolute
bound

```math
|T_j^sigma(X)|\le {9\over2}n^{3/2}.                    \tag{CF.15}
```

Consequently the perturbation in (CF.14) is at most

```math
epsilon_0n^{3/2},\qquad
epsilon_0={9\over2}(h-1)gamma_h={9\over200}.           \tag{CF.16}
```

## 5. Fixed-gap packing theorem

### Theorem CF.2 (connected `h`-bit Gram-flux packing)

Let `R(sigma)=(max Q_i^sigma)_(i=1)^h`, using either connected query family
in Section 4.  Then for all distinct `sigma,tau in {0,1}^h`,

```math
\|R(sigma)-R(tau)\|_infinity
\ge\Delta_*n^{3/2},                                    \tag{CF.17}
```

where

```math
\Delta_*
=\delta-2\epsilon_0
={3(5-\sqrt {17})\over4}-{9\over100}
=0.5676707807\ldots>0.                                 \tag{CF.18}
```

Hence any single summary and decoder answering all `h` scalar queries to
uniform error strictly below `Delta_*n^(3/2)/2` needs at least `2^h`
summary states, or `h` bits.

#### Proof

Choose `i` with `sigma_i\ne tau_i`, and orient the pair so that
`sigma_i=0,tau_i=1`.  Uniform perturbation of a function by at most
`epsilon_0n^(3/2)` perturbs its maximum by at most that amount.  Lemma CF.1
and (CF.16) therefore give

```math
\begin{aligned}
\max Q_i^sigma&\ge(M_0+B_H-\epsilon_0)n^{3/2},\\
\max Q_i^tau&\le(M_1+B_H+\epsilon_0)n^{3/2}.
\end{aligned}                                          \tag{CF.19}
```

Subtracting proves (CF.17).  If two states shared one summary codeword, the
common decoded answer to query `i` would make one of the two errors at least
half their separation. `square`

There is a slightly sharper bridge-only-port version: omit the second term
in (CF.14).  Lemma CF.1 then gives the full gap `delta` with unit-strength
connectors.  The positive `gamma_h` version is stated as the main result to
rule out the objection that non-target labels never enter their children.

## 6. Normalization and scope

The normalization is the one used in Theorems 21.19 and 21.21: `n` is the
order of **one** Walsh block and errors are measured in units of
`n^(3/2)=qn`.  The bounded-degree queries have a public connector baseline
`B_H=3h-1` for `h>=2`; the dense queries have
`B_H=9h(h-1)/2`.  These baselines cancel in the comparison (CF.19).

If instead one calls the total number of Boolean variables
`N=3hn` and normalizes by `N^(3/2)`, the displayed separation is
`Delta_*/(3h)^(3/2)` in those units.  The theorem is an `Omega(h)` response
memory lower bound at fixed one-port accuracy, not an extensive free-energy
density gap for the whole connected graph.

The construction proves scalable scalar visibility of independent
**relation-cycle fluxes**.  It does not prove that an arbitrary collection
of raw Gram entries is independently recoverable.  Nor does the argument
justify giving every non-target triangle unit weight: that is an
`Theta(h)n^(3/2)` perturbation, and the connected Boolean optimum need not
decompose into its triangle optima.  What is exact here is stronger and
narrower: arbitrary nonnegative bridge connectivity is harmless, while the
total mass of unrelated label-dependent terms is kept constant.

## 7. Verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_walsh_connected_flux_packing.py
```

The verifier checks the entire Gram/relation/root state cube through `h=5`,
the exact local sector polynomials from Theorem 21.21, both connected query
topologies, coefficient bounds, and the exact perturbation and gap
arithmetic.  The Boolean triangle bound and witness are imported from the
separately audited Theorem 21.21 verifier.
