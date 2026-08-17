# Unrooted Walsh graphs forget the label-space characteristic root

Status: rigorous task-local theorem with an exact verifier.  This draft
distinguishes the semantic state for **unrooted Walsh-graph landscapes** from
the finer rooted label-orbit carrier in Theorem 21.14.  It does not concern
arbitrary dense quadratic bridges.

## 1. Declared continuation class

Put

```math
V=\mathbb F_2^m,\qquad E=V\oplus V,\qquad q=2^m,
\qquad n=|E|=q^2.
```

Index the order-`n` Walsh matrix by `z,z' in E`:

```math
W_E(z,z')=(-1)^{z\cdot z'}.
```

For `a in V`, write

```math
\iota(a)=(0,a)\in E,\qquad
D_a(z)=(-1)^{\iota(a)\cdot z},\qquad
C_a=D_aW_ED_a.                                           \tag{US.1}
```

Given an ordered source tuple `a=(a_1,...,a_k)`, a declared continuation
chooses coefficient masks `c_v in F_2^k` for `1<=v<=t` and uses the
synchronously derived labels

```math
a[c_v]=\sum_i(c_v)_ia_i.                                \tag{US.2a}
```

It also chooses arbitrary real scalar onsite weights `h_v` and arbitrary
real symmetric edge weights `J_uv`.  Its whole Boolean landscape is

```math
\mathcal E_{\mathbf a,c}^{h,J}(x_1,\ldots,x_t)
 ={1\over2}\sum_v h_v x_v^TC_{a[c_v]}x_v
  +\sum_{u<v}J_{uv}x_u^TW_Ex_v.                         \tag{US.2}
```

This is the declared **unrooted weighted Walsh-graph** query class.  It
contains unweighted graphs, signed or weighted graphs, repetitions, and
linear combinations presented by the same coefficient masks on both sides.
It does not contain a fixed external pole, a coordinate-dependent linear
field, or an appended label specified in the old ambient coordinates rather
than intrinsically through the exposed tuple.

For a tuple define only

```math
G_{\mathbf a}=(a_i\cdot a_j)_{ij},\qquad
R_{\mathbf a}=\{c\in\mathbb F_2^k:\sum_i c_i a_i=0\}.   \tag{US.3}
```

## 2. Ambient-orbit collapse

### Theorem US.1 (unrooted ambient-orbit sufficiency)

Let `a=(a_1,...,a_k)` and `b=(b_1,...,b_k)` be tuples in `V`.  If

```math
G_{\mathbf a}=G_{\mathbf b},\qquad R_{\mathbf a}=R_{\mathbf b},             \tag{US.4}
```

then there is one coordinate permutation `P` of the `n` Boolean coordinates
such that

```math
PW_EP^T=W_E,\qquad PC_{a_i}P^T=C_{b_i}\quad(1\le i\le k).                  \tag{US.5}
```

Consequently, for every declared `(c,h,J)`, the blockwise coordinate
permutation `P^{\oplus t}` gives an exact isomorphism of the entire landscapes
in (US.2):

```math
\mathcal E_{\mathbf a,c}^{h,J}(X)
=\mathcal E_{\mathbf b,c}^{h,J}(P^{\oplus t}X).         \tag{US.6}
```

In particular all upper maxima, absolute maxima, minima, energy histograms,
and optimizer multiplicities agree.  Thus `(G,R)` is an exact carrier for
every declared unrooted weighted graph query.  The assertion is sufficiency,
not minimality: distinct `(G,R)` states can still have the same scalar graph
maxima.

#### Proof

Let

```math
\alpha(c)=\sum_i c_i a_i,\qquad
\beta(c)=\sum_i c_i b_i.
```

Equality of kernels makes

```math
\phi:\iota(\operatorname{im}\alpha)\longrightarrow
      \iota(\operatorname{im}\beta),\qquad
\phi(0,\alpha c)=(0,\beta c)                            \tag{US.7}
```

a well-defined linear isomorphism.  Equality of Gram matrices makes it an
isometry for the standard bilinear form on `E`.

The characteristic vector of the even-dimensional bilinear space `E` is

```math
\Omega_E=(\omega,\omega),\qquad \omega=(1,\ldots,1)\in V.                 \tag{US.8}
```

Every vector in either span in (US.7) has first component zero, whereas
`omega` is nonzero.  Hence neither span contains `Omega_E`.  Adjoin
`Omega_E` on both sides and extend `phi(Omega_E)=Omega_E`; this remains an
isometry because

```math
\Omega_E\cdot(0,a)=a\cdot a=(0,a)\cdot(0,a),            \tag{US.9}
```

and `phi` preserves self-pairings.  The characteristic-rooted Witt extension
lemma for the standard binary bilinear form now extends (US.7) to some

```math
O\in O(E),\qquad O^TO=I_E.                               \tag{US.10}
```

Only this bilinear extension lemma is being used here.  We are not applying
the Walsh-matrix corollary of Theorem 21.14 at a mismatched label dimension.

Let `P` act on functions by `(Pf)(z)=f(O^{-1}z)`.  Since `O` preserves dot
products,

```math
W_E(Oz,Oz')=W_E(z,z').                                  \tag{US.11}
```

Moreover `O iota(a_i)=iota(b_i)`, and hence
`O iota(a[c])=iota(b[c])` for every declared coefficient mask.  Thus

```math
(-1)^{\iota(b[c])\cdot Oz}=(-1)^{\iota(a[c])\cdot z}.   \tag{US.12}
```

Equations (US.11)--(US.12) give (US.5), and substituting the same permutation
in every block gives the pointwise identity (US.6). `square`

### Corollary US.2 (the root fibre has an exact scalable semantic collision)

Fix `k>=1`.  If `m>=3` is odd, set

```math
a=\omega,\qquad b=e_1;
```

if `m>=4` is even, set

```math
a=\omega,\qquad b=e_1+e_2.
```

The constant tuples `a^k` and `b^k` have identical `(G,R)`: `G` is the
constant matrix `a dot a=b dot b`, and `R` is the even-parity subspace of
`F_2^k`.  But their label-space characteristic-root fibres differ:

```math
R_\omega(a^k)=\{c:\sum_i c_i=1\},\qquad
R_\omega(b^k)=\varnothing.                              \tag{US.13}
```

Nevertheless Theorem US.1 identifies their entire landscape for every
unrooted real weighted graph, at every `k`.  This is a collision between
rooted orbit states, not merely an equality of one scalar maximum.

The collision disappears when the continuation declares an external root.
For example, the canonical pole/field used in Theorem 21.13 is tied to
`omega` in the original `(u,v)` splitting.  The larger ambient orthogonal
permutation in Theorem US.1 need not preserve that pole.  The projective
responses of the two singleton children are then separated by at least
`n^(3/2)/6`.  Thus the distinction is exact:

```math
\text{unrooted weighted Walsh graphs need no }R_\omega,
\quad
\text{rooted/external-field contexts can expose }R_\omega.                \tag{US.14}
```

## 3. What the wind tunnel shows

The accompanying exhaustive program freezes the following query ladder
before comparing states:

1. all labelled unweighted simple graphs on the existing blocks;
2. all signed edge weights in `{-1,0,1}` at `m=1`, and bridge weights
   `0,1,2` at `m=2`;
3. canonical `omega`-rooted fields (kept separate from the unrooted claim);
4. one-block extensions labelled by zero, `omega`, or a repeated label.

At `m=1`, unweighted graph maxima already separate all `8` rooted states at
`k=3` and all `16` at `k=4`.  At `m=2,k=2`, by contrast, the `10` algebraic
orbit states collapse to only `3` scalar classes even after weights `0,1,2`:

```math
(64,128,192),\qquad(64,104,160),\qquad(64,80,128).       \tag{US.15}
```

These are exact exhaustive values.  They prove neither a general lower bound
nor scalar minimality of `(G,R)`; they warn that orbit completeness can be
strictly finer than a particular scalar query family.  Theorem US.1, not the
finite data, proves the root-fibre collapse.

## 4. Theory consequence

The example isolates a query-relative symmetry enlargement.  The rooted
carrier `(G,R,R_omega)` is exact and orbit-complete for coordinate-rooted
Walsh responses, while the unrooted graph language admits the larger ambient
group `O(2m,2)` and factors through `(G,R)`.  Removing a query can therefore
delete a whole composition-created information resource without
approximating anything.

The next semantic question is not whether `(G,R)` is sufficient--it is--but
which part of `(G,R)` is exposed by scalar graph maxima of unrestricted
depth.  The `m=2` collapse in (US.15) shows that the answer can be strictly
smaller than the orbit state.
