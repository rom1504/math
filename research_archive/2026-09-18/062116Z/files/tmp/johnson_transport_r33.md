# Wave 33 root route: selector transport certificate and its uniform-law wall

## Status

The deterministic transport inequality below is verified and gives a simple
way to construct a low-row completion forest from partial favorable labels.
It is genuinely weaker than requiring every label to lie near the low-row
root: only one anchor pays a root distance, while the remaining labels pay
overlap disagreement plus Johnson movement.

However, this certificate cannot prove the adversarial-law event.  Under the
uniform selector law, the probability that even `s` points from the full
batch have a Johnson spanning tree of total length `O(k_0)` is
`exp{-Theta(ns)}=exp{-omega(rL_0)}`.  This is a scoped wall against proofs
which charge one unit whenever a coordinate enters a selector.  It does not
falsify the exact terminal-min-cut forest event, where free coordinates may
screen many such entries simultaneously.

## 1. Direct transport of partial projective words

Let `T` be a rooted tree.  Its root `0` carries a full spin
`z in {+-1}^n` and `S_0=[n]`.  Every other node `v` carries a projective word
`[y^v]` on `S_v`.  For an oriented edge `u -> v`, put

```math
I_{uv}=S_u\cap S_v,
\qquad
e_{uv}=|S_v\setminus S_u|.
```

Once the parent orientation has been chosen, choose the child orientation
which minimizes disagreement on `I_uv`.  Since the incidence graph is a
tree, these relative choices determine coherent node orientations
simultaneously.  Write

```math
\Delta_{uv}
=d_{\rm pr}(y^u|_{I_{uv}},y^v|_{I_{uv}}),
```

where the root word is `z`; empty overlaps contribute zero.  Construct full
representatives recursively: obey the oriented child word on `S_v` and copy
the parent bit on every coordinate outside `S_v`.  Along `uv`, leaving and
jointly absent coordinates cost zero, common pinned coordinates cost exactly
`Delta_uv`, and newly entering coordinates cost at most one each.  Therefore

```math
\boxed{
D_T(z;\{(S_v,[y^v])\})
\le\sum_{u\to v\in E(T)}(\Delta_{uv}+e_{uv}).
}
\tag{J33.1}
```

This is only an upper bound on the exact min-cut (10.923).  It deliberately
pays every entry separately; the min-cut may route or screen several entries
with one coordinate-edge cut.

If all non-root selectors have the same size `m`, orient one edge from the
root to an anchor `v_0` and every other edge between selector nodes.  Then
`e_(0v0)=0`, while

```math
e_{uv}=|S_v\setminus S_u|=d_J(S_u,S_v)
```

on selector edges.  Hence

```math
\boxed{
D_T\le
d_{\rm pr}(z_{S_{v_0}},y^{v_0})
+\sum_{uv\in E(T)\setminus\{0v_0\}}
\left[d_J(S_u,S_v)
+d_{\rm pr}(y^u|_{S_u\cap S_v},y^v|_{S_u\cap S_v})\right].
}
\tag{J33.2}
```

Thus a batch group with favorable labels, one low-row anchor, and right side
`O(k_0)` satisfies the Wave 31 forest hypothesis.  Adjacent selectors whose
labels agree projectively on their overlap cost only one.  This is the exact
positive content of local Johnson transport.

## 2. Counting tuples with a short Johnson tree

Let `N=binom(n,m)` and sample `s>=2` selectors independently and uniformly.
Let `MST_J(S_1,...,S_s)` be the minimum total Johnson length of a tree on the
`s` labeled sample positions; repetitions are allowed.  For a fixed rooted
labeled tree and nonnegative edge lengths `h_e` summing to at most `D`, after
choosing the root selector there are

```math
\binom m{h_e}\binom{n-m}{h_e}\le n^{2h_e}
```

choices for each child.  Cayley's formula and stars-and-bars therefore give

```math
\boxed{
\Pr\{MST_J(S_1,\ldots,S_s)\le D\}
\le
N^{1-s}s^{s-2}\binom{D+s-1}{s-1}n^{2D}.
}
\tag{J33.3}
```

The count overcounts tuples which admit several trees or length vectors, so
the inequality direction is safe.  It also covers zero-length edges and
repeated selectors.

At fixed density `m/n -> p in (0,1)`, use the project scales

```math
\Omega(n^{1/4+c_0-\eta}\log n)
\le s=\lceil r/T\rceil
\le O(n^{1/4+c_0}\log n),
\qquad
D=Ck_0=O(n^{3/4-c_0}/\log n),
```

where `eta<c_0<1/4`.  Since `D/s -> infinity`, the logarithm of the three
positive counting factors after `N^(1-s)` is

```math
O(s\log s+s\log(D/s)+D\log n)=O(s\log n+L_0)=o(ns),
```

whereas `(s-1)log N=Theta(ns)`.  Thus

```math
\boxed{
\Pr\{MST_J(S_1,\ldots,S_s)\le Ck_0\}
=e^{-\Theta(ns)}=e^{-\omega(rL_0)}.
}
\tag{J33.4}
```

The matching lower exponential order follows from the event
`S_1=...=S_s`, whose probability is exactly `N^(1-s)`.

The full mesoscopic batch has `r` positions and `J_(r,s)` may choose any
`s` of them.  A union bound costs only `log binom(r,s)=O(r)`, still
`o(ns)`.  Hence

```math
\boxed{
\Pr\{\text{some }s\text{ batch positions have Johnson MST}\le Ck_0\}
=e^{-\Theta(ns)}.
}
\tag{J33.5}
```

This is far below the required `e^{-O(rL_0)}=e^{-O(n log n)}`.
Allowing unobserved selector Steiner nodes does not change the exponent: in
any metric the terminal MST is at most twice a Steiner-tree cost (double the
Steiner tree and shortcut an Euler tour).  Thus a Johnson Steiner tree of
cost `D` implies the event in (J33.3) with `2D`.

## 3. Consequence

Equations (J33.1)--(J33.2) remain useful for selector laws already
concentrated in a Johnson cluster, and they precisely identify what local
transport buys.  But a universal proof cannot upper-bound the exact forest
cost by a statistic which always includes the selector Johnson MST: the
uniform selector law makes that sufficient event super-exponentially too
rare at the project scale.  Any successful completion theorem must exploit
free-coordinate screening in the exact terminal min-cut, or a direct
cross-selector global structure which avoids paying coordinate entry one by
one.
