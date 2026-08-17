# Linear Walsh labels: pairwise Gram data is not an extremal quotient

Status: rigorous task-local counterexample, with an exact verifier.  This
complements the constant-label commutation example in
[`walsh_family_composition_holonomy.md`](walsh_family_composition_holonomy.md).
It concerns only the explicit cap-`1/2` Walsh family.

## 1. The frozen candidate state

Put `q=2^m`, `n=q^2`, and let `F=W/q` be the normalized Walsh involution on
`F_2^m times F_2^m`.  For a linear truth-table label

```math
g_a(v)=a\cdot v
```

write

```math
\widehat C_a=D_aFD_a.                                    \tag{LG.1}
```

The corresponding hollow sign quadratic has exact cap `n^(3/2)/2`.

For a word `a=(a_1,...,a_k)`, freeze the candidate carrier

```math
\mathsf{Gram}(\mathbf a)=(a_i\cdot a_j)_{i,j\le k}.       \tag{LG.2}
```

This is a different invariant from the truth-table overlaps used by the
one-bridge packing argument.  It knows every label parity and binary
bilinear phase, and is invariant under global orthogonal changes of Walsh
coordinates.  In characteristic two it can nevertheless forget equality of
labels.  It has only `k(k+1)/2` bits, independent of `m`.

The question is whether (LG.2) determines the Boolean extremal response of a
composed linear-label system.  The answer is no, already on a three-vertex
path.

## 2. Two Gram-identical words

For every `m>=3`, use the two labels

```math
a=(1,1,1,0,\ldots,0),
\qquad
b=(0,0,1,0,\ldots,0).                                   \tag{LG.3}
```

They obey

```math
a\cdot a=b\cdot b=a\cdot b=1.                            \tag{LG.4}
```

Consequently the two path words

```math
\mathbf a^-=(a,a,a),
\qquad
\mathbf a^+=(a,b,a)                                     \tag{LG.5}
```

have exactly the same candidate state: both Gram matrices are the all-one
`3 by 3` matrix.  Notice what the Gram state has forgotten.  In the first
word all three labels are equal; in the second only the first and third are.
The difference `a+b` is a nonzero isotropic vector orthogonal to both labels,
and is therefore invisible to their Gram matrix.

This collision is **not** invisible to pairwise truth-table overlap.  For
linear labels,

```math
S(g_c,g_d)=q\,\mathbf 1_{c=d}.                            \tag{LG.5a}
```

Thus the constant word has overlap `q` on every pair, while `(a,b,a)` has
zero overlap on the two pairs crossing `a` and `b`.  The theorem below
rejects binary Gram alone; it is complementary to, rather than stronger
than, the earlier overlap obstruction.

For a label word on the path `1--2--3`, define

```math
E_{(c_1,c_2,c_3)}(x_1,x_2,x_3)
=\frac q2\sum_{i=1}^3x_i^T\widehat C_{c_i}x_i
 +q x_1^TFx_2+q x_2^TFx_3.                              \tag{LG.6}
```

Every spin has squared norm `n`.

## 3. An explicit bent witness

Write coordinates as `(u_0,...,u_{m-1},v_0,...,v_{m-1})`.  Define

```math
Q_m(u,v)
=u_0u_1+v_0v_1+u_2v_2+v_0+v_1+v_2
 +\sum_{j=3}^{m-1}u_jv_j,                                \tag{LG.7}
```

and put `x(u,v)=(-1)^(Q_m(u,v))`.

Modulation by `a` cancels the last three linear terms in (LG.7):

```math
Q_m+a\cdot v
=u_0u_1+v_0v_1+u_2v_2+\sum_{j=3}^{m-1}u_jv_j.            \tag{LG.8}
```

Each two-variable factor `(-1)^(rs)` is a `+1` eigenvector of the normalized
order-four Walsh transform.  Therefore

```math
F(D_ax)=D_ax,
\qquad
\widehat C_ax=x.                                         \tag{LG.9}
```

Modulation by `b` cancels only `v_2`:

```math
Q_m+b\cdot v
=u_0u_1+(v_0v_1+v_0+v_1)+u_2v_2
 +\sum_{j=3}^{m-1}u_jv_j.                                \tag{LG.10}
```

The factor `(-1)^(rs+r+s)` is a `-1` eigenvector of the normalized
order-four Walsh transform, while every other displayed factor has sign
`+1`.  Hence

```math
F(D_bx)=-D_bx,
\qquad
\widehat C_bx=-x.                                        \tag{LG.11}
```

The quadratic polar form in (LG.7) is nonsingular, so `y=Fx` is Boolean.
This also follows immediately by tensoring the explicit two-variable Walsh
identities just used.

Because `b\cdot b=1`, the involutions `F` and `\widehat C_b` anticommute.
Thus (LG.11) gives

```math
\widehat C_b(Fx)=-F(\widehat C_bx)=Fx.                   \tag{LG.12}
```

The spin triple `(x,Fx,x)` therefore lies in the maximizing eigenspace of
each of the three children in the word `a^+`, and it saturates both Walsh
bridges.

## 4. The leading-scale failure

### Theorem LG.1 (Gram-identical extensive separation)

For every `m>=3`, the two words in (LG.5) satisfy

```math
\mathsf{Gram}(\mathbf a^-)=\mathsf{Gram}(\mathbf a^+),    \tag{LG.13}
```

but

```math
\max E_{\mathbf a^+}=\frac72n^{3/2},                     \tag{LG.14}
```

whereas

```math
\max E_{\mathbf a^-}
\le\frac{3\sqrt3}{2}n^{3/2}.                             \tag{LG.15}
```

In particular, a state retaining all pairwise label Gram/commutation data
can miss at least

```math
\frac{7-3\sqrt3}{2}n^{3/2}
=0.901923\ldots\,n^{3/2}.                                \tag{LG.16}
```

#### Proof

The witness from Section 3 gives three child contributions `nq/2` and two
bridge contributions `nq`.  This is the sum of the separate operator-norm
upper bounds, proving (LG.14).

For the constant odd word, the normalized global block operator is

```math
\mathcal M^-=I_3\otimes\widehat C_a+A(P_3)\otimes F.     \tag{LG.17}
```

Anticommutation cancels its cross term after squaring:

```math
(\mathcal M^-)^2=(I_3+A(P_3)^2)\otimes I_n.              \tag{LG.18}
```

Since `||A(P_3)||=sqrt(2)`, one has `||\mathcal M^-||=sqrt(3)`.
The complete Boolean block vector has squared norm `3n`; multiplying its
Rayleigh bound by `q/2` proves (LG.15).  Equation (LG.4) proves (LG.13).

## 5. Characteristic-root obstruction and a quotient that survives

The counterexample says that relations matter.  Even Gram plus the full
relation kernel is not sufficient for responses in fixed Walsh coordinates:
the ambient bilinear space has a characteristic root which every orthogonal
symmetry must fix.

Let

```math
\omega=(1,\ldots,1)\in F_2^m,
```

and, for a tuple `mathbf a`, define

```math
\begin{aligned}
G(\mathbf a)&=(a_i\cdot a_j)_{ij},\\
\mathcal R(\mathbf a)&=\{c\in F_2^k:\sum_i c_i a_i=0\},\\
\mathcal R_\omega(\mathbf a)&=
 \{c\in F_2^k:\sum_i c_i a_i=\omega\}.
\end{aligned}                                             \tag{LG.20}
```

The last set is either empty or a coset of `mathcal R`.

### Theorem LG.2 (the characteristic root is response-visible)

Let `m>=3` be odd, take `a=omega` and any unit coordinate vector `b=e_j`, and
consider the singleton tuples `(a)` and `(b)`.  They have the same Gram matrix
`(1)` and the same relation kernel `{0}`, but

```math
\mathcal R_\omega((a))=\{1\},
\qquad
\mathcal R_\omega((b))=\varnothing.                      \tag{LG.21}
```

This missing bit is visible at leading scale through one Walsh bridge.  Let
`s_c(u,v)=(-1)^(u dot v+c dot v)` and `y_c=q^(-1)Ws_c`.  For the response

```math
\mathcal V_c(y)=\max_x\{H_{g_c}(x)+x^TWy\},              \tag{LG.22}
```

the matched-pole identity and the Walsh resolvent bound give

```math
\mathcal V_a(y_a)=\frac32n^{3/2},
\qquad
\mathcal V_b(y_a)\le\frac43n^{3/2}.                      \tag{LG.23}
```

Indeed, the distinct linear truth tables have correlation zero, so the
off-pole Rayleigh parameter is `rho=0`.  Reversing `a,b` gives the reverse
gap at `y_b`.  Hence their projective response distance is at least
`n^(3/2)/6`.  Thus `(G,\mathcal R)` is not a semantic quotient for rooted
Walsh futures, even though the two single children have the same unperturbed
cap.

### Theorem LG.3 (rooted relation-form orbit classification)

Suppose `m` is odd.  Two ordered label tuples `mathbf a,mathbf a'` have the
same three objects in (LG.20) if and only if there is an orthogonal linear
map

```math
O\in O(m,2),
\qquad Oa_i=a_i'\quad(1\le i\le k).                       \tag{LG.24}
```

Consequently, at fixed odd `m`, (LG.20) is an exact quotient for the Boolean
maximum of (WC.3) on **every** graph `G`: tuples with the same state give
permutation-conjugate quadratic landscapes and hence exactly equal maxima.
The same holds for response functions when their exposed Boolean coordinates
are relabelled by the same `O`.

#### Proof

Equality of relation kernels makes

```math
\phi:\operatorname{span}\{a_i\}\longrightarrow
      \operatorname{span}\{a_i'\},
\qquad \phi(a_i)=a_i'                                    \tag{LG.25}
```

a well-defined linear isomorphism.  Equality of `G` makes it an isometry.
Equality of `\mathcal R_\omega` says precisely that, if the source span contains
`omega`, then `phi(omega)=omega`; if it does not, neither does the target
span.

Here is a self-contained reduction to the symplectic Witt lemma.  Since `m`
is odd, `omega dot omega=1` and

```math
H=\omega^\perp
```

is a nondegenerate alternating space.  Every `u in F_2^m` has the unique
decomposition

```math
u=(u\cdot u)\omega+h(u),
\qquad h(u)\in H.                                        \tag{LG.26}
```

The map `h(u) -> h(phi(u))` is well defined on the image of the source span:
the only kernel issue is a multiple of `omega`, handled by the rooted-coset
condition.  It preserves the alternating form because `phi` preserves both
`u dot v` and the parities `u dot u`.  The finite-dimensional symplectic Witt
extension lemma extends it to an isometry `S` of all of `H`.  Then

```math
O(c\omega+h)=c\omega+S h                                 \tag{LG.27}
```

is orthogonal, fixes `omega`, and extends `phi`.  This proves the nontrivial
direction.  The reverse direction is immediate: every orthogonal map fixes
`omega`, because `x dot x=omega dot x` characterizes it, and it preserves all
three objects in (LG.20).

Finally, applying `O` simultaneously to the `u` and `v` coordinates is a
coordinate permutation preserving the Walsh kernel.  It sends every child
`C_{a_i}` to `C_{a_i'}` and leaves every Walsh bridge fixed, proving the
extremal claim.  It also fixes the canonical query `y_omega`, so the same
state is sufficient for graph responses with that rooted future.

The state in (LG.20) has an `O(k^2)`-bit presentation, independent of `m`, so
for `k=o(m)` it is strictly smaller than storing the `km` label bits.  This is
an actual Boolean-extremal equivalence, not merely the coefficient
factorization `R tensor K`.

It is not, however, an independently composable state.  Gluing two separately
summarized tuples requires their cross-Gram entries and new cross-relations;
neither is determined by the two isolated states.  Thus LG.3 is a globally
coordinated orbit quotient, while the missing cross-form is exactly
composition-created information.

## 6. Interpretation

This is a scalable counterexample to the binary-Gram Clifford quotient.
Binary bilinear data is not enough because the characteristic-two form can
hide label equalities and linear dependencies.  Pairwise truth-table overlap
already detects this particular equality pattern, so no conclusion against
that invariant follows from LG.1.  Composition
turns that hidden relation into compatibility of three maximizing
eigenspaces, changing the energy at leading scale.

A natural next candidate for even `m`, or for independently composed pieces,
augments the Gram state by the relation code

```math
\mathcal R(\mathbf a)
=\{c\in F_2^k:\sum_i c_i a_i=0\}                          \tag{LG.28}
```

Indeed, `\mathcal R` distinguishes the two words in (LG.5); the theorem proves
that any sufficient state must at least detect this particular hidden
relation, not that it must literally retain the entire code.  This draft does
**not** claim that `(Gram,\mathcal R)` is
sufficient.  In characteristic two, extension of a subspace isometry can
also depend on its position relative to the canonical all-one vector, and
independently summarized components still lack their cross-Gram data.

Nor does the theorem show that the `m`-bit labels themselves are
incompressible.  They already give a strict coefficient presentation of the
linear subfamily, exponentially smaller than the `q`-bit arbitrary truth
tables.  What fails is the much smaller pairwise invariant as an actual
Boolean-extremal quotient.

## 7. Exact verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_linear_label_walsh_gram.py
```

The verifier checks the Gram collision, factorized Walsh eigen-identities,
bentness, the exact good-word energy, the bad-word anticommutator and squared
operator identity, and the stated constants at `m=3,4`.  It also audits the
characteristic-root response input and exhaustively verifies the rooted orbit
classification for all ordered triples at `m=3`.  All arithmetic in the
structural checks is integral.
