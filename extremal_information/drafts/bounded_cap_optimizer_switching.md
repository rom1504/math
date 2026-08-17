# Bounded-cap optimizer fibres and witness-cover complexity

**Status.** Rigorous task-local draft.  The first theorem is a scalable
counterexample to a cardinality bound on a common-optimizer fibre.  The
second gives the response-sensitive replacement: a target-scale response
packing forces an extensive approximate witness dictionary.  This is not a
restatement of UP.1.  UP.1 treats one universal pin against the entire child
class; the results below treat an arbitrary restricted family, arbitrary
many futures, and approximate optimizer covers.

Throughout,

```math
H_A(x)=\sum_{i<j}A_{ij}x_ix_j,
\qquad Q(A)=\max_x|H_A(x)|.
```

## 1. Common-optimizer fibres can remain quadratically large

Fix `u in {+-1}^k`.  Append one future spin `y`, join it to old vertex `i`
with sign `u_i`, and write

```math
P_A(x,y)=H_A(x)+y\,u\mathbin\cdot x,
\qquad
g_u(x)=\max_y y\,u\mathbin\cdot x=|u\mathbin\cdot x|.       \tag{WS.1}
```

This is a complete exact signing on `k+1` vertices.

### Theorem WS.1 (a huge bounded-cap common-witness fibre)

For every `k>=3` and every `u`, there are at least

```math
2^{{k\choose2}-k}                                  \tag{WS.2}
```

distinct complete sign children `A` such that

```math
Q(A)\le2k^{3/2},
\qquad
Q(P_A)\le2k^{3/2}+k,                               \tag{WS.3}
```

and `u` and `-u` are old-spin optimizers of

```math
H_A(x)+g_u(x).                                     \tag{WS.4}
```

Consequently bounded cap, even together with a fixed nontrivial exact
quadratic future, places no subexponential upper bound on the number of
children sharing one old optimizer.

#### Proof

Let `E=binom(k,2)`.  Vertex switching by
`s in {+-1}^k/{+-1}` acts freely on the `2^E` complete signings, so there
are exactly

```math
2^{E-k+1}                                           \tag{WS.5}
```

switching orbits.  In every orbit choose a signing `A` and a Boolean state
`z` maximizing `H_A`.  The switched signing
`A^z=D_zAD_z` has `1` as a maximizer, because

```math
H_{A^z}(x)=H_A(z\odot x)\le H_A(z)=H_{A^z}(\mathbf1).       \tag{WS.6}
```

Switch once more by `u`; this makes `u` a maximizer.  Thus every switching
orbit has at least one representative with the prescribed optimizer.

For a uniformly random signing and fixed `x`, `H_A(x)` is a sum of `E`
independent Rademacher signs.  Hoeffding and a union bound give

```math
\Pr\{Q(A)>2k^{3/2}\}
\le2^{k+1}\exp(-4k)
=2\exp\{-(4-\log2)k\}< {1\over2}.                 \tag{WS.7}
```

The event and `Q` are switching invariant, and all orbits have the same
size.  Hence more than half of the orbits are bounded-cap orbits.  Choosing
one prescribed-optimizer representative from each gives at least
`2^{E-k}` children, proving (WS.2)--the first part of (WS.3).

Both `H_A` and `g_u` are maximized at `+-u`, so (WS.4) follows.  Finally

```math
|P_A(x,y)|\le Q(A)+|u\mathbin\cdot x|
\le2k^{3/2}+k,                                     \tag{WS.8}
```

which proves the parent cap. `square`

The exponent `E-k` is not meant to be sharp.  Its role is to decisively
falsify any attempted theorem of the form “bounded parent cap forces a
small common-optimizer fibre.”  Response separation is indispensable.

## 2. One shared witness carries only constantly many macroscopic values

Let a fixed future `(B,C)` have

```math
g(x)=\max_y\{x^TBy+H_C(y)\},
\qquad
R_g(A)=\max_x\{H_A(x)+g(x)\}.                      \tag{WS.9}
```

### Lemma WS.2 (single-witness response capacity)

Suppose a family `F` has the same old optimizer `u` in (WS.9), and every
completed parent has order `N<=Lambda k` and cap at most `CN^(3/2)`.  A
subfamily whose scalar responses `R_g(A)` are pairwise separated by
`epsilon k^(3/2)` has size at most

```math
1+\left\lfloor{2C\Lambda^{3/2}\over\epsilon}\right\rfloor.
                                                               \tag{WS.10}
```

#### Proof

For every child in the common-witness fibre,

```math
R_g(A)=H_A(u)+g(u).                                 \tag{WS.11}
```

The parent cap puts every response in
`[-C Lambda^(3/2)k^(3/2),C Lambda^(3/2)k^(3/2)]`.
The packing number of this interval at the declared separation is (WS.10).
`square`

Thus Theorem WS.1 and Lemma WS.2 coexist without tension: a fibre may hold
quadratically many bits of child identity, but one fixed bounded response
coordinate reveals only `O(1)` target-scale bins.

## 3. Approximate witness dictionaries

Let `mathcal G` be any declared family of futures `g`.  Define the one-sided
response vector

```math
R_A(g)=\max_x\{H_A(x)+g(x)\}.                       \tag{WS.12}
```

A set `U subset {+-1}^k` is a **tau-witness cover** for a child family
`mathcal F` and `mathcal G` at the `k^(3/2)` scale if

```math
0\le R_A(g)-R_A^U(g)\le\tau k^{3/2},
\qquad
R_A^U(g)=\max_{u\in U}\{H_A(u)+g(u)\},             \tag{WS.13}
```

for every `A in mathcal F` and `g in mathcal G`.  This is strictly weaker
than storing an optimizer for every child--future pair: one dictionary is
shared, and only near-optimality is required.

### Theorem WS.3 (macroscopic response packing forces extensive witnesses)

Assume every child in `mathcal F` obeys

```math
Q(A)\le C_0k^{3/2}.                                 \tag{WS.14}
```

Suppose `mathcal F` is pairwise separated in declared response sup metric:

```math
\sup_{g\in\mathcal G}|R_A(g)-R_{A'}(g)|
\ge\epsilon k^{3/2}\qquad(A\ne A').               \tag{WS.15}
```

If `U` is a `tau`-witness cover with `2tau<epsilon`, then

```math
|\mathcal F|
\le
\left(1+\left\lceil{2C_0\over\epsilon-2\tau}\right\rceil
\right)^{|U|}.                                      \tag{WS.16}
```

Equivalently,

```math
|U|\ge
{\log|\mathcal F|\over
 \log(1+\lceil2C_0/(\epsilon-2\tau)\rceil)}.       \tag{WS.17}
```

In particular an `exp(alpha k)` bounded-cap response packing requires
`|U|=Omega(k)` at every fixed `epsilon,tau,C_0`.

#### Proof

Maxima are nonexpansive in the sup norm, so for every future

```math
|R_A^U(g)-R_{A'}^U(g)|
\le\max_{u\in U}|H_A(u)-H_{A'}(u)|.                \tag{WS.18}
```

Together with (WS.13), response separation implies

```math
\max_{u\in U}|H_A(u)-H_{A'}(u)|
\ge(\epsilon-2\tau)k^{3/2}.                        \tag{WS.19}
```

Thus the evaluation vectors
`(H_A(u))_(u in U)` form an `ell_infinity` packing at that separation.
Every coordinate lies in the interval
`[-C_0k^(3/2),C_0k^(3/2)]`.  Partitioning each interval into the displayed
number of half-open bins proves (WS.16), then (WS.17). `square`

The child-cap hypothesis follows automatically from a uniform parent-cap
hypothesis.  Indeed, for a quadratic parent

```math
P_A(x,y)=H_A(x)+x^TBy+H_C(y),                       \tag{WS.20}
```

uniform averaging over `y` gives

```math
H_A(x)=\mathbb E_yP_A(x,y),
```

so `Q(A)<=Q(P_A)`.  If every parent has `N<=Lambda k` and cap
`CN^(3/2)`, one may take `C_0=C Lambda^(3/2)` in WS.3.

### Corollary WS.3a (context--switching dichotomy)

Suppose the declared language has `q` futures, and for future `j` all chosen
old witnesses across the child family belong to a set `U_j`.  Exact
witnesses give `U=union_j U_j`, and hence an `exp(alpha k)` response packing
forces

```math
\sum_{j=1}^q|U_j|\ge c(\alpha,\epsilon,C_0)k.       \tag{WS.21}
```

Therefore either the compiler carries an extensive fixed dictionary of
query witnesses, or optimizer identity switches across children with
extensive total support.  If `q=O(1)`, at least one future uses `Omega(k)`
distinct old witnesses.  If every future has one common optimizer, then
necessarily `q=Omega(k)`.

This is the strongest deterministic conclusion available without charging
the query language itself.  A claim that bounded cap alone forces
child-dependent switching is false by WS.1; a compiler may instead pay an
extensive dictionary of common query pins.

## 4. Projective and absolute variants

For projective response distance

```math
d_proj(A,A')={1\over2}\operatorname{osc}_{g\in\mathcal G}
                    (R_A(g)-R_{A'}(g)),             \tag{WS.22}
```

the same proof applies after anchoring one coordinate of the evaluation
vector.  If `d_proj(A,A')>=epsilon k^(3/2)` and (WS.13) holds, then

```math
|\mathcal F|
\le
\left(1+\left\lceil{4C_0\over\epsilon-2\tau}\right\rceil
\right)^{|U|-1},                                   \tag{WS.23}
```

for `2tau<epsilon`.  The inessential factor two comes from anchoring
`H_A(u)-H_A(u_0)` in an interval of length `4C_0k^(3/2)`.

An absolute parent response can be handled by enlarging the witness set to
signed witnesses `(s,u) in {+-1} times {+-1}^k` and evaluating
`sH_A(u)`.  The proof is unchanged.  What matters is the number of reusable
signed old witnesses, not whether the outer response is one-sided or
absolute.

## 5. Consequence for flat Gram response families

The short-seed Gram family has `exp(Omega(k))` children separated at a fixed
`k^(3/2)` response scale.  Therefore any hypothetical `N=O(k)` bounded-cap
disjoint compiler which preserves that packing and whose old optimization
is uniformly approximable from one reusable witness dictionary must use
`Omega(k)` old witnesses.  This conclusion allows arbitrary auxiliary
interaction and arbitrary child-dependent tie breaking.

It does **not** rule out such a linear dictionary, nor does it prove that the
active witness must depend on the child rather than only on the query.  The
sharp remaining alternatives are now explicit:

1. `Omega(k)` common query pins;
2. `Omega(k)` total child-dependent optimizer support;
3. failure of every `o(k)` approximate witness cover, meaning the future
   genuinely uses a larger portion of the Boolean landscape.

Eliminating alternative 1 requires an additional exact-sign/cap theorem
tailored to a restricted flat child family.  UP.1 eliminates it only for a
single pin universal over **all** complete sign children.

## 6. Research judgment

The hoped-for cardinality statement is decisively false, while its
response-sensitive version is true and quantitative.  Bounded cap controls
the number of macroscopically distinguishable values per witness, not the
number of children in an optimizer fibre.  The correct conserved resource is
the size of an approximate old-witness cover, with a precise dichotomy
between common query dictionaries and child-dependent switching.

This is a genuine narrowing but not yet an impossibility theorem for the
flat exact-sign compiler: an extensive dictionary is compatible with the
linear hidden rate.  The next theorem would have to charge the realization
of `Omega(k)` common pins under bounded cap, or prove that the short-seed
family has no `o(k)` approximate witness cover even after query-dependent
gauging.
