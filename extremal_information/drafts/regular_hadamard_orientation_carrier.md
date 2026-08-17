# The minimal gauge carrier for regular-Hadamard orientation systems

**Status.** Rigorous task-local theorem, gluing law, and finite falsifier.
This note treats only the restricted family in which every old block is one
of the two orientations of a single regular Hadamard signing and every
bridge is a scalar signed copy of the same Hadamard matrix.  It does not
claim a quotient for arbitrary exact-sign bridges.

The point is not the elementary Kronecker factorization by itself.  The new
content is the exact projective-switching quotient, its canonical minimal
coordinate count, and a gluing fibre law which separates two different
sources of composition-created information:

```math
\text{relative marginal orientations}
\quad\oplus\quad
\text{new cycle holonomies}.
```

## 1. The restricted multi-block family

Let `n=q^2`, and let `H` be a symmetric regular Hadamard matrix satisfying

```math
H^2=nI,
\qquad H\mathbf1=q\mathbf1,
\qquad \operatorname{tr}H=0.                         \tag{OC.1}
```

Put `A=H-diag(H)`.  Fix a finite graph `G=(V,E)`.  Give every vertex an
onsite orientation `sigma_i in {+-1}` and every edge a public bridge sign
`b_ij in {+-1}`.  The Boolean energy is

```math
\mathcal E_{G,\sigma,b}(x)
=\frac12\sum_{i\in V}\sigma_i x_i^THx_i
 +\sum_{ij\in E}b_{ij}x_i^THx_j.                    \tag{OC.2}
```

The trace assumption makes (OC.2) exactly the energy of hollow onsite
signings `sigma_i A`; no diagonal correction remains on the Boolean cube.
Define the scalar block matrix

```math
T_{ii}=\sigma_i,
\qquad
T_{ij}=T_{ji}=b_{ij}\ (ij\in E),
\qquad
T_{ij}=0\ (ij\notin E).                              \tag{OC.3}
```

Then, for the concatenated Boolean vector `X`,

```math
\mathcal E_{G,\sigma,b}(X)
=\frac12X^T(T\otimes H)X.                            \tag{OC.4}
```

This factorization is the specialization of the archived Kronecker carrier
WC.1.  We will not count it as a theorem.

## 2. Closed caps and projective switching

For `d:V->{+-1}`, write `D=diag(d_i)`.  Blockwise spin inversion gives

```math
\mathcal E_{D T D}(X)=\mathcal E_T((D\otimes I)X),    \tag{OC.5}
```

while the outer absolute value gives `Q_H(-T)=Q_H(T)`.  Thus introduce

```math
T\sim_{\rm ps}T'
\quad\Longleftrightarrow\quad
T'=\epsilon DTD
\quad\text{for some }\epsilon\in\{+-1\},\ D^2=I.    \tag{OC.6}
```

This is ordinary signed-graph switching together with one global
antipode.  The antipode is essential: it flips both all onsite orientations
and all bridge signs.

### Theorem OC.1 (exact projective-switching cap carrier)

For every regular Hadamard `H`, the complete Boolean energy multiset changes
only by a possible global sign on each projective-switching class.  In
particular,

```math
Q_H(T):={1\over2}\max_{X\in\{+-1\}^{n|V|}}
             |X^T(T\otimes H)X|                       \tag{OC.7}
```

factors exactly through `T/~ps`.

If `G` has `k` vertices, `e` edges, and `c` connected components, then the
labelled projective-switching quotient has exactly

```math
2^{e+c-1}                                               \tag{OC.8}
```

elements.  Equivalently, it has `e+c-1` independent binary coordinates.
For connected `G`, it has exactly one bit per edge, irrespective of `n`.

#### Canonical coordinates

Choose a spanning forest `F`, roots `r_1,...,r_c`, and distinguish `r_1`.
Use the global antipode to make `sigma_(r_1)=+1`; use vertex switching to
make every forest edge positive.  The remaining data are

```math
\bigl(\sigma_v\sigma_{r_1}:v\ne r_1\bigr)
\quad\text{and}\quad
\bigl(\sigma_{r_1}^{|C_e|}\prod_{f\in C_e}b_f:
             e\in E\setminus F\bigr),                \tag{OC.9}
```

where `C_e` is the fundamental cycle of the chord `e`.  There are

```math
(k-1)+(e-k+c)=e+c-1                                  \tag{OC.10}
```

such bits.  They determine the class and are independent.

#### Proof

Equation (OC.5) is a Boolean bijection.  Multiplication of the whole energy
by `-1` preserves its absolute cap, proving the first assertion.

Ordinary vertex switching has effective rank `k-c` on the edge signs: a
constant switch on a connected component fixes every edge.  The global
antipode contributes one further free action because onsite signs are
nonzero.  Hence the orbit quotient of the `k+e` coefficient bits has
dimension

```math
k+e-(k-c)-1=e+c-1.
```

For a constructive proof, first use the global antipode to normalize
`sigma_(r_1)`.  The forest recursion then uniquely normalizes all forest
edges after choosing one switch at each root.  A chord sign in this
representative is the
displayed modified cycle product: a length-`ell` cycle product acquires
`sigma_(r_1)^ell` when the antipode is normalized.  Conversely the data in
(OC.9) reconstruct the normalized representative.  This proves
(OC.8)--(OC.10). `square`

The count is minimal as a **coefficient-conjugacy carrier**: two different
coordinate words in (OC.9) are not related by any Boolean block inversion
and global output sign.  It is not asserted that the single scalar cap is
injective on all these classes.  Section 4 proves that both kinds of bit in
(OC.9) can nevertheless carry leading-order Boolean information.

## 3. Exact gluing and the compatibility fibre

Closed marginal caps discard one antipode per disconnected piece.  Those
antipodes cannot be chosen independently after the pieces interact.

Let `G_1,...,G_s` be connected, vertex-disjoint signed graphs.  Add `r`
cross edges so that the graph on the `s` pieces induced by the cross edges
is connected; necessarily `r>=s-1`.  Let `G` denote the joined graph.

### Theorem OC.2 (orientation--cycle gluing law)

Fix the projective-switching carrier of every piece `G_a`.  The fibre of
joined projective-switching carriers over those fixed marginal carriers has
exactly

```math
2^r                                                     \tag{OC.11}
```

elements.  After choosing marginal representative/gauge sections and a
cross-edge spanning tree, its `r` independent compatibility bits split as

```math
\underbrace{s-1}_{\text{relative marginal antipodes}}
\quad+\quad
\underbrace{r-s+1}_{\text{cross-edge cycle holonomies}}
=r.                                                     \tag{OC.12}
```

This is an exact gluing theorem: piece carriers, the `s-1` relative
orientations, and the `r-s+1` fundamental cross-cycle products reconstruct
the joined carrier, hence determine every Boolean cap (OC.7).  Repeated
gluing is associative because it is simply union of signed graphs followed
by the quotient (OC.6).

#### Proof

If `e_a=|E(G_a)|`, Theorem OC.1 gives `e_a` bits for the carrier of each
connected piece.  The joined graph is connected and has
`sum_a e_a+r` edges, hence its carrier has `sum_a e_a+r` bits.  Restriction
to the pieces is onto, so every fibre has `r` bits.  More explicitly, choose
one representative of the first piece.  Each of the other `s-1` marginal
classes has two antipodal representatives relative to it.  These are the
relative-orientation bits.  Component-constant vertex switches normalize
the signs on a spanning tree of `s-1` cross edges.  Each remaining cross
edge closes one independent fundamental cycle, whose product is invariant.
The dimensions of the two summands are intrinsic, while their coordinate
split depends on these choices.  This constructs all `2^r` joined classes
without repetition and proves
(OC.11)--(OC.12). `square`

The law also identifies the framing issue.  An unrooted closed carrier is
reusable across a single bridge only after supplying its relative antipode.
If several bridges join the same pieces, path transports (or, equivalently,
the new cycle products) must also be supplied.  A carrier that stores a
numeric bridge frame at every labelled port restores these data but, when
every vertex is a future port, also restores the discarded forest gauges.
Thus quotient size depends on whether contexts are gauge-covariant or use
fixed numeric port frames; (OC.11) is the intrinsic gauge-covariant law.

## 4. Both compatibility resources are Boolean-visible

The bit count above is not merely orbit bookkeeping.

### Proposition OC.3 (a relative antipode has leading response)

Join two one-vertex pieces by one positive bridge.  The coefficient matrices

```math
T_+=\begin{pmatrix}1&1\\1&1\end{pmatrix},
\qquad
T_-=\begin{pmatrix}1&1\\1&-1\end{pmatrix}             \tag{OC.13}
```

have identical closed one-piece carriers but different relative antipodes.
Their caps obey

```math
Q_H(T_+)=2qn,
\qquad
Q_H(T_-)\le\sqrt2\,qn.                                \tag{OC.14}
```

Thus the first summand in (OC.12) can change a cap by
`(2-sqrt(2))n^(3/2)`.

This is exactly the two-block holonomy calculation BH.2: the matrices in
(OC.13) have norms `2` and `sqrt(2)`, and the positive case is attained at
the regular pole.

### Proposition OC.4 (a cycle bit has leading response)

Join three positive one-vertex pieces in a triangle.  If the bridge product
is positive, switching makes all three bridges positive and

```math
T_{\rm bal}=J_3,
\qquad Q_H(T_{\rm bal})={9\over2}qn.                    \tag{OC.15}
```

If the bridge product is negative, switching gives

```math
T_{\rm unbal}=
\begin{pmatrix}1&1&1\\1&1&-1\\1&-1&1\end{pmatrix},
\qquad \|T_{\rm unbal}\|=2,                            \tag{OC.16}
```

and therefore

```math
Q_H(T_{\rm unbal})\le3qn.                              \tag{OC.17}
```

The single cycle holonomy changes the cap by at least

```math
{3\over2}qn={3\over2}n^{3/2}.                          \tag{OC.18}
```

Indeed, every separate term in the balanced triangle is maximized by the
regular pole, proving (OC.15).  The spectral bound for a Boolean vector of
squared norm `3n` proves (OC.17).  This is the second summand in (OC.12).

Consequently neither relative marginal orientations nor cycle holonomies
may be deleted from a uniform exact reusable carrier.  The result does not
say that every bit is exposed by the bare cap on every fixed graph; it says
that each resource class has a fixed-size, target-scale witness.

## 5. Spectrum is not an exact carrier

Switching class determines the spectrum of `T`, but the converse is false,
and the spectrum does not determine the Boolean cap.  Already at `n=4`, take

```math
H_4=\begin{pmatrix}
1&1&1&-1\\
1&-1&1&1\\
1&1&-1&1\\
-1&1&1&1
\end{pmatrix}.                                           \tag{OC.19}
```

It satisfies (OC.1) with `q=2`.  The two scalar block matrices

```math
T_0=\begin{pmatrix}
-1&-1&0&-1\\-1&-1&0&-1\\0&0&-1&0\\-1&-1&0&1
\end{pmatrix},
\qquad
T_1=\begin{pmatrix}
-1&-1&-1&-1\\-1&-1&0&0\\-1&0&-1&0\\-1&0&0&1
\end{pmatrix}                                            \tag{OC.20}
```

have the same characteristic polynomial

```math
\lambda^4+2\lambda^3-3\lambda^2-4\lambda,              \tag{OC.21}
```

but exhaustive Boolean evaluation gives

```math
Q_{H_4}(T_0)=32,
\qquad Q_{H_4}(T_1)=34.                                 \tag{OC.22}
```

Their zero patterns are respectively a triangle plus an isolated vertex and
a three-leaf star, so they are not projectively switching equivalent.
Equation (OC.22) is only a finite falsifier; no asymptotic assertion relies
on it.

## 6. Scope and consequence

1. **Strict state.**  On a connected `k`-block graph, the closed exact cap
   is controlled by `e` binary coordinates, independent of the Hadamard
   order `n` and exponentially smaller than the `kn`-spin landscape for
   fixed `G` (more generally when `e=o(kn)`).  This
   is much narrower than WC.1 because only the two onsite orientations of
   one fixed `H` are allowed.
2. **Dynamic repair.**  Closed absolute responses discard a global
   antipode.  Gluing `s` pieces through `r` bridges recreates exactly `s-1`
   relative antipodes and `r-s+1` cycle fluxes.  The one-edge BH.2 example
   is the smallest nontrivial fibre of the general law.
3. **Minimality scope.**  The quotient is minimal under the declared
   block-inversion/global-output conjugacies, and both kinds of compatibility
   bit are macroscopically observable in uniform families.  We do not claim
   injectivity of one scalar cap on every projective-switching class.
4. **No original-problem claim.**  Arbitrary dense bridges do not share one
   common Hadamard factor, so this theorem does not compress the original
   signing problem.  It supplies an exact model of how marginal absolute
   compression loses precisely quantifiable information under interaction.

## 7. Verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_regular_hadamard_orientation_carrier.py
```

The verifier checks projective-switching invariance, the canonical orbit
count for all graphs through four vertices, singleton-piece gluing fibres,
the orientation and triangle spectral identities, and the exact cospectral
cap collision (OC.19)--(OC.22).  The independent audit verifier additionally
checks non-singleton gluing fibres.
