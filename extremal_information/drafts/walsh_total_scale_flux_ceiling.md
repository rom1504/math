# A total-scale ceiling for state-local Walsh flux architectures

Status: task-local rigorous draft.  This isolates an architecture-level
obstruction to amplifying Theorem 21.22.  It applies to arbitrary public
scalar query families, arbitrary signed or dense public bridge terms, and
shared encodings of flux bits.  The essential hypothesis is that the hidden
Walsh state enters only through bounded onsite child terms.  The conclusion
does **not** cover a genuinely state-dependent cross-block interaction.

## 1. State-local response systems

Partition the Boolean variables into `k` blocks of a common order `n`, and
put `N=kn`.  A hidden state `s` supplies one Walsh label
`a_v(s)` to every active block.  The maps `s -> a_v(s)` may be arbitrary;
in particular a public query may request a fixed linear combination of a
larger source tuple.  Put

```math
H_{a_v}(x_v)={q\over2}x_v^TC_{a_v}x_v,
\qquad q=\sqrt n,                                      \tag{TC.1}
```

where every normalized Walsh involution `C_a` is symmetric and orthogonal.
Thus, for Boolean `x_v`,

```math
|H_a(x_v)|\le {1\over2}n^{3/2}.                        \tag{TC.2}
```

A **state-local public query** `theta` is any scalar landscape of the form

```math
E_s^\theta(X)
=S_\theta(X)+\sum_{v=1}^k d_{\theta v}H_{a_v}(x_v),
\qquad |d_{\theta v}|\le D.                            \tag{TC.3}
```

Here `S_theta` is completely arbitrary except that it is independent of the
marked state `a`.  In particular it may contain a complete signed weighted
graph of Walsh bridges, higher public terms, or a large public baseline;
neither positivity nor a norm bound on `S_theta` is used.

Write

```math
R_\theta(s)=\max_XE_s^\theta(X)                         \tag{TC.4}
```

and equip marked states with the scalar public-response metric

```math
d_\Theta(s,t)=\sup_{\theta\in\Theta}
 |R_\theta(s)-R_\theta(t)|.                            \tag{TC.5}
```

The same statements below hold for minima or absolute maxima.

## 2. Public bridges cannot amplify onsite state

### Theorem TC.1 (state-local total-scale diameter ceiling)

Every state-local public query family satisfies

```math
d_\Theta(s,t)
\le D\,\#\{v:a_v(s)\ne a_v(t)\}\,n^{3/2}
\le Dk n^{3/2}
={D\over\sqrt k}N^{3/2}.                               \tag{TC.6}
```

Consequently, if `k -> infinity`, the complete response image has diameter
`o(N^(3/2))`.  No subset containing even two states, and hence no
`2^{Omega(h)}` packing, can have a fixed positive separation in units of
`N^(3/2)`.

#### Proof

The public term cancels **before** optimization:

```math
E_s^\theta(X)-E_t^\theta(X)
=\sum_{v:a_v(s)\ne a_v(t)}d_{\theta v}
   \{H_{a_v(s)}(x_v)-H_{a_v(t)}(x_v)\}.               \tag{TC.7}
```

By (TC.2), one changed block contributes at most `D n^(3/2)`.
Taking the supremum over `X` and using

```math
|\max f-\max g|\le\|f-g\|_\infty                     \tag{TC.8}
```

proves (TC.6).  The same one-Lipschitz fact applies to `min f` and
`max |f|`. `square`

This proves that nonnegative connective padding in Theorem 21.22 is not the
source of the dilution.  Negative weights, a complete graph, an expander,
or an arbitrarily large public bridge baseline do not help: all such terms
are part of `S_theta` and cancel in (TC.7).  To obtain normalized separation
`epsilon` through this architecture one must have

```math
D\ge\epsilon\sqrt k.                                  \tag{TC.9}
```

Thus coefficient amplification has to diverge as `sqrt(k)`; after the
natural normalization making coefficients bounded, the separation again
vanishes.

## 3. Exact scaling for the disjoint triangle cube

In Theorem 21.22 there are `k=3h` Walsh blocks.  If `sigma_i` changes, the
label `u_i` stays fixed and only `v_i,w_i` change.  Hence for every bounded
state-local public query,

```math
d_\Theta(\sigma,\tau)
\le 2D\,d_H(\sigma,\tau)n^{3/2}.                       \tag{TC.10}
```

In particular the diameter of **every** code cut out of the flux cube is

```math
\operatorname {diam}(R(\mathcal C))
\le2Dh n^{3/2}
={2D\over3^{3/2}\sqrt h}N^{3/2}.                       \tag{TC.11}
```

Therefore no choice of a positive-rate error-correcting subcode and no
aggregate public query yields a constant total-scale gap.  For the entire
`2^h` cube, a Hamming-neighbour pair gives the still sharper packing ceiling

```math
\min_{\sigma\ne\tau}d_\Theta(\sigma,\tau)
\le2Dn^{3/2}
={2D\over(3h)^{3/2}}N^{3/2}.                           \tag{TC.12}
```

The lower bound `Delta_* n^(3/2)` of Theorem 21.22 is therefore of the
correct `h^{-3/2}` order for the full cube under unit-bounded onsite
weights.  Coding can improve the upper scale only to `h^{-1/2}`, never to a
positive constant.

## 4. Unequal disjoint cells: a metric-entropy ceiling

The preceding phenomenon is not an artifact of equal block sizes.  Suppose
a hidden bit `sigma_i` affects only a disjoint cell of `s_i` Boolean
variables and that, uniformly over public scalar queries,

```math
\|E_\sigma^\theta-E_\tau^\theta\|_\infty
\le L\sum_{i:\sigma_i\ne\tau_i}s_i^{3/2}.             \tag{TC.13}
```

Put `N=sum_i s_i`.

### Theorem TC.2 (state-local response packing has bounded entropy)

For every `epsilon>0`, a response code whose distinct pairs have distance
**strictly greater** than `epsilon N^(3/2)` has

```math
\log_2|\mathcal C|
\le \left\lceil {4L^2\over\epsilon^2}\right\rceil.   \tag{TC.14}
```

In particular unequal allocation cannot turn disjoint local bits into a
positive total-scale information rate.

#### Proof

Order the cells so that `s_1>=s_2>=...`.  Since
`i s_i<=N`,

```math
\sum_{i>r}s_i^{3/2}
\le N^{3/2}\sum_{i>r}i^{-3/2}
\le {2N^{3/2}\over\sqrt r}.                            \tag{TC.15}
```

If `|C|>2^r`, two codewords agree on the first `r` bits.  Equations
(TC.13)--(TC.15) put their response distance at most
`2LN^(3/2)/sqrt(r)`.  Taking
`r=ceil(4L^2/epsilon^2)` proves (TC.14); the case `L=0` is immediate.
`square`

The theorem is a response-rate--distortion statement, not merely a norm
estimate: at fixed macroscopic distortion, the number of distinguishable
states is bounded independently of the number of local hidden bits.

## 5. Shared Walsh labels in the unrooted equivariant language

One may try to store `h` flux bits in the quadratic Gram/relation data of a
smaller shared label set.  The following count applies only to the
coordinate-equivariant unrooted weighted Walsh-graph language of Theorem
21.18: no external pole, coordinate pin, or coordinate field is allowed.

### Corollary TC.3 (shared-label Walsh ceiling, scoped)

Suppose `r` source labels are represented in `t` actual child slots, with
`t>=r`.  Their complete unrooted orbit state `(G,R)` has at most

```math
2^{r(r+1)/2}\,2^{r^2}
=2^{(3r^2+r)/2}                                       \tag{TC.16}
```

possible values.  Therefore a family of `2^h` pairwise contextually
different states requires

```math
r\ge {\sqrt{1+24h}-1\over6}.                           \tag{TC.17}
```

For every unit-bounded state-local public query family within this language,
its total-scale response diameter is consequently at most

```math
{diam R\over(tn)^(3/2)}\le {1\over\sqrt t}\le {1\over\sqrt r}
\le
\left({6\over\sqrt{1+24h}-1}\right)^{1/2}
=O(h^{-1/4}).                                          \tag{TC.18}
```

#### Proof

There are at most `2^{r(r+1)/2}` binary symmetric Gram matrices.  Every
relation kernel is a subspace of `F_2^r`; representing a subspace by at most
`r` spanning vectors gives the crude upper bound `2^{r^2}`.  Theorem 21.18
says `(G,R)` already determines the entire unrooted weighted-graph
landscape, so distinct scalar responses cannot outnumber (TC.16).  Solving
`h<=(3r^2+r)/2` gives (TC.17), and Theorem TC.1 on the `t` actual slots gives
(TC.18). `square`

This covers sharing triangle vertices, dense Gram codes, and arbitrary
unrooted equivariant public Walsh graphs, as long as hidden-state dependence
remains in bounded onsite children and every source label occupies a child
port.  It does **not** cover coordinate-pinning baselines: those break the
ambient Witt symmetry and can distinguish labels with the same `(G,R)`.
Nor does it apply if many abstract labels are queried through only `O(1)`
derived ports.  Under the stated language and slot condition, shared encoding
can at best weaken the elementary `h^(-1/2)` diameter ceiling to
`h^(-1/4)`; it cannot yield a constant gap.

## 6. Exact boundary of the no-go theorem

The hypotheses identify what an escape must change.  At least one of the
following is necessary for positive total-scale response separation:

1. **Cross-block state dependence.**  One hidden compatibility bit must
   alter a macroscopic family of bridge coefficients, rather than only a
   bounded number of onsite children.
2. **Unbounded coefficients.**  Onsite weights must grow at least as
   `sqrt(k)`, which leaves the bounded-sign normalization relevant to the
   motivating problem.
3. **A different accuracy scale.**  The fixed one-port `n^(3/2)` scale of
   Theorem 21.22 is meaningful, but it is not the total `N^(3/2)` scale.
4. **A differently normalized local or vector output.**  A continuation may
   expose local bits one by one while normalizing only by the queried local
   order, or may retain a response vector rather than one scalar full-support
   landscape.  Those outputs lie outside TC.1's total-scale scalar metric.
   Coordinate-pinning scalar baselines remain inside TC.1 and still obey its
   diameter ceiling; they are excluded only from TC.3's orbit count because
   they break coordinate equivariance.

Thus the connected Walsh-flux packing cannot be amplified within its
current coefficient architecture.  A successful dense-interface lower
bound has to make compatibility information enter cross-block coefficients
or prove a different joint mechanism; adding more public connectors cannot
do it.

There is a quantitative version of the first escape condition.  Suppose a
more general `k`-block query is a sum of unit-weight atoms, each having
Boolean supremum at most `B n^(3/2)`.  If flipping hidden bit `i` can alter
at most `d_i` atoms in any one query, then the same pointwise argument gives

```math
d_\Theta(\sigma,\sigma\mathbin\oplus e_i)
\le2B d_i n^{3/2}.                                    \tag{TC.19}
```

Hence a full hidden cube with neighbour separation
`epsilon(kn)^(3/2)` requires

```math
d_i\ge {\epsilon\over2B}k^{3/2}                       \tag{TC.20}
```

for every bit.  A state-dependent-crossing escape must therefore broadcast
each locally flippable flux into `Omega(k^(3/2))` unit-scale interaction
atoms in at least one public query.  Merely allowing bounded-degree
state-dependent bridges still cannot work.  This incidence bound does not
rule out a deliberately nonlocal broadcast, and it should not: such a
broadcast is a genuinely different architecture whose semantic legitimacy
has to be assessed separately.

## 7. Verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_walsh_total_scale_flux_ceiling.py
```

The exact verifier checks (TC.10)--(TC.12), the tail inequality used in
(TC.15), exact Gaussian-binomial subspace counts against (TC.16), and the
quadratic inversion (TC.17).
