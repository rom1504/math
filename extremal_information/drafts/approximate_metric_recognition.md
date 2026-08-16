# Approximate tropical idempotence: recognition and long-depth failure

**Status.** Theorems proved below and independently adversarially audited.
This draft tests whether the metric-shell algebra can be recognized from an
approximate idempotence law without retaining a full long-composition state.

The repair in AMR.1 is deliberately elementary rather than optimal.  The
sharp universal coefficient is
`max(1/2,(q-2)/q)`, obtained by shifting by that amount and then taking
shortest-path closure; it is proved in Theorem 16.7.

## 1. One-step recognition is dimension-free

Let `X` have `q` elements.  For a symmetric real kernel `K` with
`K(x,x)=0`, write

```math
(K\star K)(x,z)=\min_y\{K(x,y)+K(y,z)\},
```

and

```math
\delta(K)=\|K-K\star K\|_\infty.
```

Because the choices `y=x,z` give `K star K<=K`, this is a one-sided
defect.

### Theorem AMR.1 (sharp-order recognition)

The following hold.

1. `delta(K)` is the largest triangle violation:

   ```math
   \delta(K)=\max_{x,y,z}
   \{K(x,z)-K(x,y)-K(y,z)\}.                       \tag{AMR.1}
   ```

2. If `delta=delta(K)`, then

   ```math
   d_K(x,z)=
   \begin{cases}
   0,&x=z,\\
   K(x,z)+\delta,&x\ne z
   \end{cases}                                    \tag{AMR.2}
   ```

   is a metric (a pseudometric when `delta=0`), and

   ```math
   \|K-d_K\|_\infty\le\delta.                      \tag{AMR.3}
   ```

3. Conversely, if a metric `d` satisfies `||K-d||_infinity<=eta`, then

   ```math
   \delta(K)\le3\eta.                              \tag{AMR.4}
   ```

In particular symmetric diagonal-zero min-plus idempotents are exactly
finite pseudometrics, and approximate idempotence recognizes a nearby metric
with constants independent of `q`.

#### Proof

The first assertion follows by subtracting the minimum defining `K star K`.
For distinct `x,y,z`, (AMR.1) gives

```math
K(x,z)+\delta
\le K(x,y)+K(y,z)+2\delta
=d_K(x,y)+d_K(y,z).
```

Triangles with a repeated endpoint are tautological.  Taking `x=z` in
(AMR.1) and using symmetry gives `K(x,y)>=-delta/2`; hence every off-diagonal
entry of `d_K` is nonnegative (strictly positive if `delta>0`).  This proves
(AMR.2)--(AMR.3).  Finally, if `K` is within `eta` of a metric, then

```math
K(x,z)\le d(x,z)+\eta
\le d(x,y)+d(y,z)+\eta
\le K(x,y)+K(y,z)+3\eta,
```

which proves (AMR.4). `square`

If a kernel has a bijective zero-centre map `g`, relabeling row `a` by
`g(a)` reduces its recognition to this theorem.  Exact symmetry and
idempotence after that untwisting recognize a single twisted metric shell.
For a *composable fixed-metric family*, one must additionally verify that
every permitted twist is an isometry; a bijective zero pattern does not imply
this.  This is a mathematical characterization, not yet a compressed test:
reading all triangle defects can itself require the full kernel.

## 2. One-step recognition does not give long-depth stability

The positive theorem above is insufficient for composition.

### Theorem AMR.2 (linear accumulation despite metric proximity)

Fix integers `q>=3`, `1<=T<=q-1`, and reals `a>delta>0`.  On
`X_q={0,...,q-1}`, define

```math
K_\delta(i,j)=
\begin{cases}
0,&i=j,\\
a|i-j|-\delta,&i\ne j.
\end{cases}                                      \tag{AMR.5}
```

Then:

```math
\delta(K_\delta)=\delta,                           \tag{AMR.6}
```

`K_delta` is within `delta` of the path metric `a|i-j|`, but its `T`th
min-plus power is

```math
K_\delta^{\star T}(i,j)=
\begin{cases}
0,&i=j,\\
a|i-j|-\delta\min\{T,|i-j|\},&i\ne j.
\end{cases}                                      \tag{AMR.7}
```

Consequently, whenever `T<=q-1`,

```math
\|K_\delta^{\star T}-K_\delta\|_\infty
=(T-1)\delta,                                    \tag{AMR.8}
```

and already on row zero the projective shape distance is

```math
d_{\rm sh}\bigl(
[K_\delta^{\star T}(0,\cdot)],
[K_\delta(0,\cdot)]
\bigr)={(T-1)\delta\over2}.                       \tag{AMR.9}
```

Taking `delta=c/q` and `T=q-1` gives kernels whose one-step idempotence
defect and distance to an exact metric vanish, while their long-depth
response-shape drift tends to `c/2`.

#### Proof

For `i<k<j`, the direct cost exceeds the two-step cost by exactly `delta`;
all other triangle violations are no larger, proving (AMR.6).  A path using
`p` nonzero moves of integer lengths `l_1,...,l_p` has cost

```math
a\sum_r|l_r|-\delta p.
```

Its total variation is at least `|i-j|`.  If `p<=|i-j|`, this is at least
`a|i-j|-delta p`; if `p>|i-j|`, it is at least
`(a-delta)p>(a-delta)|i-j|`.  Equality is attained by a monotone path using
`min(T,|i-j|)` nonzero pieces and zero stays for the remaining factors.
This proves (AMR.7).  Equations (AMR.8)--(AMR.9) follow by comparing an
adjacent endpoint, where the error is zero, with one at distance at least
`T`, where it is `(T-1)delta`. `square`

## 3. What information composition creates

The omitted datum in (AMR.5) is a per-nonzero-transition toll.  At one step
it is only `delta`, but after `T` compositions the optimizer can collect it
once per useful segment.  The exact state must retain the available hop
budget (capped by the path distance), not merely the nearby metric.

This is a clean benchmark instance of microscopic information becoming
macroscopic under composition.  It gives two limits on the current theory.

1. Approximate tropical idempotence is a dimension-free **one-step**
   recognition criterion but not a depth-stable quotient.
2. A useful approximate shell theorem needs an additional global condition:
   for example exact projection/retraction after each step, a bounded useful
   path length, or a semilattice family whose composition is closed before
   approximation.

Checking every triangle is not itself compression, and assuming uniform
control of all powers would merely restate the desired conclusion.  The next
question is to identify a structured hypothesis that rules out transition
tolls using less than the complete kernel.
