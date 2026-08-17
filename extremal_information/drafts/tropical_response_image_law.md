# Tropical response-image dimension and reusable message complexity

Status: main-agent synthesis draft.  The tropical-polytope component is
classical; the contextual metric, rate statement, and separation of static
from reusable complexity are the project-level deductions.

## 0. General finite response-family law

Let `S` be a compact finite ordinary polyhedral complex in
`R^q/R1`, of topological dimension `d`, and let
`R_m(g)=max_j(m_j+g_j)`.  Define

```math
d_proj([m],[m'])
=inf_(c in R)sup_(g in R^q)|R_m(g)-R_(m')(g)-c|.             \tag{TR.0}
```

### Theorem TR.0 (polyhedral response rate)

Arbitrary coordinate futures make `S` the coarsest exact projective response
state,

```math
d_proj([m],[m'])=inf_c||m-m'-c1||_infinity
={1\over2}osc(m-m'),                                        \tag{TR.0a}
```

and, for fixed `S` as `epsilon downarrow0`,

```math
Cov(S,d_proj,epsilon)=Theta_S(epsilon^(-d)).                 \tag{TR.0b}
```

The same exponent governs optimal uniform response codebooks up to the usual
packing/covering factor of two.  A bounded scalar baseline contributes one
additional `O(log(B/epsilon))` bits; an unbounded baseline has no finite
absolute codebook.

#### Proof

The maximum inequality gives one side of (TR.0a), and futures strongly
pinning any chosen coordinate recover every coordinate difference.  Best
calibration is the midrange of `m-m'`.  A finite `d`-dimensional polyhedral
complex has an `O_S(epsilon^(-d))` cover, while a relatively open ball in a
maximal cell gives a matching packing. `square`

## 1. A finite-port kernel

Let `K in R^(p times q)` be finite and define its max-plus message map

```math
(U_Ku)_j=max_(i<=p)\{u_i+K_(ij)\}.                            \tag{TR.1}
```

Normalize vectors modulo the additive gauge `R1`, for example by
`max_jm_j=0`, and write

```math
\mathcal C_K=\{[U_Ku]:u\in R^p\}\subset R^q/R1.              \tag{TR.2}
```

This is the projective max-plus row space of `K`.  Let `d_K` be its ordinary
topological dimension, equivalently the maximum ordinary dimension in the
standard type decomposition.  This is independent of subdivision and is not
a claim about any of the competing tropical matrix ranks.  Necessarily
`d_K<=min(p-1,q-1)`.

### Theorem TR.1 (all-finite tropical-kernel specialization)

The set `mathcal C_K` is a compact finite polyhedral complex.  If an
arbitrary future profile `g in R^q` is allowed at the output port, then:

1. `mathcal C_K` is the coarsest exact projective response state;
2. writing `R_m(g)=max_j(m_j+g_j)`, its projective contextual metric is

   ```math
   d_ctx([m],[m'])
   :=inf_(c\in R)sup_g|R_m(g)-R_(m')(g)-c|
   ={1\over2}osc(m-m');                                      \tag{TR.3}
   ```

3. defining `N_K(epsilon)=Cov(mathcal C_K,d_ctx,epsilon)`, for fixed `K` and
   `epsilon downarrow0`,

   ```math
   N_K(epsilon)=Theta_K(epsilon^(-d_K)),
   \qquad
   log_2N_K(epsilon)=d_Klog_2(1/epsilon)+O_K(1).              \tag{TR.4}
   ```

The hidden constants depend on the scale, cell count, and conditioning of
this fixed `K`; dimension alone is not a uniform finite-error bound for a
sequence `K_n`.

#### Proof

For a selector pattern `s:[q]->[p]`, the inputs satisfying

```math
u_(s(j))+K_(s(j),j)>=u_i+K_(ij)\quad(i<=p,j<=q)              \tag{TR.5}
```

form a polyhedron, and `U_K` is affine on it.  There are finitely many
patterns.  Moreover, for every output pair,

```math
min_i(K_(ij)-K_(ik))
<= (U_Ku)_j-(U_Ku)_k
<=max_i(K_(ij)-K_(ik)).                                      \tag{TR.6}
```

The linear image of each selector polyhedron in the quotient is a polyhedron
and hence closed.  Bound (TR.6) makes it bounded, and a finite common
refinement of these images is a polyhedral complex.  This is the classical
type decomposition of a projective tropical polytope.

For a future `g`, the closed response is `R_m(g)`.  Without scalar
calibration, `sup_g|R_m(g)-R_(m')(g)|=||m-m'||_infinity` by coordinate
pinning.  Optimizing a common response offset replaces the sup norm by the
best constant approximation to `m-m'`, which is half its oscillation and
proves (TR.3).  Thus no exact quotient can identify two different projective
messages.

Theorem TR.0 now proves the response metric and (TR.4). `square`

The polyhedral-complex statement is classical tropical convexity; see
[Develin--Sturmfels](https://arxiv.org/abs/math/0308254).  Formula (TR.4) is
the response rate--distortion consequence, not a claim that tropical
dimension itself is new.

## 2. Exact tree and serial composition

At a tree node whose state alphabet has size `p`, with unary score `a_i` and
child messages `m_c=b_c1+r_c in R^p`, let its edge to the parent have kernel
`K in R^(p times q)`.  Form

```math
u_i=a_i+sum_c r_c(i),
\qquad w=U_Ku,

b_out=sum_cb_c+max_jw_j,
\qquad r_out=w-(max_jw_j)1.                                  \tag{TR.7}
```

Thus `(b_out,r_out)` is an exact congruence under arbitrary-depth max-sum
composition, including heterogeneous kernels on different edges.  If
arbitrary unary scores are available, every point of `mathcal C_K` is
realizable, so Theorem TR.1 is also a minimality theorem, not only an upper
bound.

Every max-plus map is nonexpansive in both sup norm and the projective
half-oscillation metric.  Consequently a one-time `epsilon` approximation of
a message remains an `epsilon` approximation through any exact serial
one-hole future.  Errors in several approximated child messages add before
the outgoing map.  This does **not** make an arbitrary nearest-net encoder
reusable: re-encoding after every step can add a fresh error each time.

There is an exact reusable finite **projective** carrier under an arithmetic
closure hypothesis.  If all kernel and unary entries belong to `eta Z`, then
normalized messages lie in

```math
\mathcal C_K\cap(eta Z)^q,                                   \tag{TR.8}
```

a finite set by (TR.6), and exact updates preserve the lattice; the absolute
baseline can still be unbounded.  Rounding arbitrary real data first gives
an exact carrier for the rounded kernels `K_tilde`, not for the original
ones.  If `F` microscopic factor tables are rounded once to mesh `eta`, every
complete assignment and conditional optimum changes by at most `Feta/2`.
Thus fixed total error needs `eta=O(epsilon/F)`.  This introduces no further
update-by-update drift, but a fixed mesh does not give a depth-independent
approximation as `F` grows.

This separates the resources cleanly:

```math
static lossy information = d_K log(1/epsilon)+O_K(1),

reusable dynamic information
= static image + an invariant arithmetic/semantic congruence.             \tag{TR.9}
```

## 3. Benchmarks predicted by the law

- **Viterbi / weighted automata.**  Theorem TR.0 applied to the bounded-spread
  family of unrestricted `q`-endpoint survivor vectors gives a projective
  polytope of dimension `q-1`,
  so the same response-metric proof gives its
  `Theta((B/epsilon)^(q-1))` codebook.  This is not literally a compact
  `mathcal C_K` for an all-finite kernel: the identity port uses forbidden
  transitions (`-infinity`) and its unrestricted projective image is
  unbounded.  An all-finite max-plus rank-one transition has `d_K=0` after
  that step and is exactly a memory reset.
- **Ferromagnetic Potts tree.**  With
  `K_(ij)=K 1_(i=j)`, the image is
  `{r in [-K,0]^q:max r=0}`, so `d_K=q-1`.  The clipped cavity-message law and
  its sharp rate follow immediately.
- **Parity trellis.**  After reversing max/min signs, Theorem TR.0 applies to
  the bounded unrestricted message cube realized by weighted prefixes on
  `q=2^r` syndrome states.
  When one baseline is carried separately, its projective dimension is
  `q-1`; literal absolute costs have the `q`-dimensional rate stated in the
  benchmark.  Forbidden completions use infinite semiring entries, so this
  is a response-geometry analogue rather than a literal all-finite instance
  of Theorem TR.1.
- **Grouped kernels.**  If `g` is onto and
  `K_(ij)=A_(i,g(j))+v_j` for `s` output groups, then every message has the
  form `v_j+w_(g(j))`; hence `d_K<=s-1`, with equality when the grouped image
  has interior.  The theorem discovers a strict state smaller than the raw
  `q`-entry table.

The same formula therefore distinguishes three notions that raw separator
width conflates: the number of physical labels, the dimension of the
realizable response image, and the existence of a finite reusable dynamic
carrier.

## 4. Ceiling

The theorem is finite-port.  If `q` grows with the interface, `d_K` can grow
as fast as `q-1`; parity trellises then have doubly exponential codebook size
as a function of their syndrome dimension.  Nor does tropical image
dimension control arbitrary bilinear bridges, whose query port itself has
exponential size.  The bridge-packing theorems are therefore genuine
obstructions, not missing cases of (TR.4).
