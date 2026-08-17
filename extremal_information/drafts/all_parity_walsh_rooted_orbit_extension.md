# Rooted linear-label Walsh orbits in both parities

Status: rigorous task-local extension and audit.  This note audits Theorem
LG.3 of [`linear_label_walsh_gram_obstruction.md`](linear_label_walsh_gram_obstruction.md)
and proves its even-dimensional analogue.  It does not modify the canonical
theory files.

## 1. Setup

Let

```math
V=\mathbb F_2^m,qquad B(x,y)=x\cdot y,qquad
\omega=(1,\ldots,1).
```

The vector `omega` is the characteristic vector of `B`:

```math
B(x,x)=B(x,\omega)\quad(x\in V).                         \tag{AP.1}
```

Consequently every `B`-orthogonal map fixes `omega`.  For an ordered tuple
`a=(a_1,...,a_k)` put

```math
\begin{aligned}
G(\mathbf a)&=(B(a_i,a_j))_{i,j},\\
\mathcal R(\mathbf a)&=\{c\in\mathbb F_2^k:\sum_i c_i a_i=0\},\\
\mathcal R_\omega(\mathbf a)&=
 \{c\in\mathbb F_2^k:\sum_i c_i a_i=\omega\}.
\end{aligned}                                             \tag{AP.2}
```

The last set is empty or one coset of `mathcal R`.  Thus, conditional on
`(G,mathcal R)`, it records exactly the coordinate of the characteristic
root in the presented span, when that root is present.

## 2. Audit of the odd-dimensional proof

### Lemma AP.1 (the projected map in LG.3 is well defined)

Suppose `m` is odd, `U,U'<=V`, and `phi:U->U'` is an isometry.  Assume

```math
\omega\in U\Longleftrightarrow\omega\in U',
\qquad
\phi(\omega)=\omega\quad\hbox{when }\omega\in U.         \tag{AP.3}
```

Put `H=omega^perp` and

```math
h(u)=u+B(u,u)\omega.
```

Then

```math
h(u)\longmapsto h(\phi u)                                \tag{AP.4}
```

is a well-defined symplectic isometry from `h(U)` to `h(U')`.

#### Proof

Because `B(omega,omega)=1`, one has `H` nondegenerate alternating and
`ker h=<omega>`.  If `h(u)=h(v)`, then `u+v` is zero or `omega`.  The second
case can occur inside `U` only when `omega in U`, and (AP.3) then makes the
two proposed images equal.  This proves well-definedness and injectivity.

For `p=B(u,u)` and `q=B(v,v)`, direct expansion gives

```math
B(h(u),h(v))=B(u,v)+pq.                                  \tag{AP.5}
```

Both terms on the right are preserved by `phi`, so (AP.4) preserves the
alternating form. `square`

The symplectic Witt extension lemma now extends (AP.4) to `S in Sp(H)`, and

```math
O(c\omega+h)=c\omega+Sh                                  \tag{AP.6}
```

extends `phi`.  Thus the argument in LG.3 is sound.  The phrase “on the
image of the source span” should be read as “on the subspace `h(U)`”; the
kernel check above is the only edge case.  Degeneracy of `U` or `h(U)` is
harmless because the symplectic Witt lemma applies to arbitrary subspaces,
not only nondegenerate ones.

## 3. What changes in even dimension

When `m` is even, `B(omega,omega)=0`.  Hence `omega^perp` has radical
`<omega>`, so the odd proof cannot be reused: its projected ambient space is
degenerate.  The orbit theorem nevertheless remains true.  The missing
structure is an affine-symplectic shear, not another tuple invariant.

Choose a nondegenerate alternating complement `W` and a vector `e` so that

```math
V=\langle e,\omega\rangle\perp W,
\quad
B(e,e)=B(e,\omega)=1,
\quad
\omega^\perp=\langle\omega\rangle\oplus W,              \tag{AP.7}
```

where `e` and `omega` are orthogonal to `W` (but not to each other).  Such a
choice exists: lift a symplectic basis of
`omega^perp/<omega>` to obtain `W`, then adjust any vector pairing to one
with `omega` by an element of `W`.

### Lemma AP.2 (even orthogonal group in affine-symplectic coordinates)

For every `S in Sp(W)`, `t in W`, and `c in F_2`, the rules

```math
\begin{aligned}
T_{S,t,c}(\omega)&=\omega,\\
T_{S,t,c}(w)&=Sw+B(t,Sw)\omega,\\
T_{S,t,c}(e)&=e+t+c\omega
\end{aligned}                                             \tag{AP.8}
```

define an orthogonal map of `V`.  Conversely every orthogonal map has this
form, uniquely.

#### Proof

The displayed images preserve the pairings among `omega`, `e`, and `W`:
the two copies of `B(t,Sw)` cancel in `B(T(e),T(w))`, and all other checks
are immediate from alternation of `W`.  Thus `T` is orthogonal.  Conversely,
an orthogonal map fixes `omega`, induces some `S in Sp(W)` on
`omega^perp/<omega>`, and has `T(e)=e+t+c omega`.  Orthogonality to `T(e)`
forces the `omega` coefficient of `T(w)` to be `B(t,Sw)`.  This proves both
existence and uniqueness. `square`

In particular there is a bijective affine-symplectic parameterization

```math
O(m,2)\longleftrightarrow Sp(W)\times W\times\mathbb F_2. \tag{AP.9}
```

The multiplication in these coordinates contains the cocycle
`B(t,S t')`, so (AP.9) is deliberately not asserted to be the displayed
direct-product group.  It gives group orders `2` at `m=2` and `48` at
`m=4`, matching exact enumeration.

This also makes the even-dimensional one-vector orbits transparent.  For
`m>=4`, the nonzero vectors split into exactly three orbits:

```math
\{\omega\},\qquad
\{x\ne0,\omega:B(x,x)=0\},\qquad
\{x:B(x,x)=1\}.                                         \tag{AP.9a}
```

Indeed `Sp(W)` is transitive on the nonzero vectors of `W`, the shear `t`
can toggle the `omega` coordinate when the `W` component is nonzero, and
`T_{I,t,c}(e)` reaches every anisotropic vector.  Thus even before considering
tuples, the characteristic root is an orbit type that Gram parity alone
cannot name.

## 4. The all-parity orbit theorem

### Theorem AP.3 (rooted relation-form orbit classification for every `m`)

Let `m>=1`.  Two ordered tuples `a,a' in V^k` have equal triples

```math
(G(\mathbf a),\mathcal R(\mathbf a),
  \mathcal R_\omega(\mathbf a))                          \tag{AP.10}
```

if and only if an `O in O(m,2)` satisfies `Oa_i=a_i'` for every `i`.

#### Proof

Equality of relation kernels makes

```math
\phi:U=\operatorname{span}\{a_i\}\longrightarrow
U'=\operatorname{span}\{a_i'\},\qquad a_i\mapsto a_i'  \tag{AP.11}
```

a well-defined isomorphism.  Gram equality makes it an isometry.  Equality
of the rooted fibres says exactly that `omega` belongs to both spans or to
neither and, in the former case, `phi(omega)=omega`.

For odd `m`, Lemma AP.1 and symplectic Witt extension give the result.  It
remains to treat even `m`.  If `omega` is absent from `U`, enlarge both
spaces by `omega` and extend `phi` by `phi(omega)=omega`.  This remains an
isometry by (AP.1); the rooted-fibre condition ensures the same direct-sum
extension is possible on the target.  Hence assume from now on that both
spaces contain `omega`.

Use (AP.7) and put `U_0=U cap omega^perp`.  Since `omega in U`, projection
along `omega` shows

```math
U_0=\langle\omega\rangle\oplus P
```

for a subspace `P<=W`; similarly
`U_0'=<omega> direct-sum P'`.  There are a symplectic isometry
`S_0:P->P'` and a linear functional `ell:P->F_2` such that

```math
\phi(w)=S_0w+\ell(w)\omega.                              \tag{AP.12}
```

Extend `S_0` to `S in Sp(W)` by the symplectic Witt lemma.

If `U=U_0`, choose `t in W` representing the functional `ell` on `S(P)`:

```math
B(t,Sw)=\ell(w)\quad(w\in P).                            \tag{AP.13}
```

This is possible by nondegeneracy of `W`.  Then (AP.8), with arbitrary `c`,
extends `phi`.

Otherwise the parity functional `u -> B(u,omega)` has kernel `U_0`, so
choose one odd vector

```math
z=e+w_0+b\omega\in U,qquad
\phi(z)=e+w_0'+b'\omega.                                \tag{AP.14}
```

Isometry of `z` against (AP.12) gives, for every `w in P`,

```math
\ell(w)=B(w_0,w)+B(w_0',S_0w).                          \tag{AP.15}
```

Set

```math
t=w_0'+Sw_0,qquad
c=b'+b+B(t,Sw_0).                                       \tag{AP.16}
```

Then (AP.15) is exactly (AP.13), so (AP.8) agrees with `phi` on `U_0`, and
a direct substitution shows `T_{S,t,c}(z)=phi(z)`.  Thus it extends `phi`
on all of `U`.

Conversely, every orthogonal map preserves Gram and relations and fixes
`omega` by (AP.1), hence preserves the rooted fibre. `square`

In a precise orbit-theoretic sense the rooted fibre is the minimal extra
datum over `(G,mathcal R)`: Theorem AP.3 says it is sufficient, while any
invariant that classifies orthogonal orbits must distinguish two realizations
having different fibres.  Equivalently one may store the symbol “absent” or
the single class `[c] in F_2^k/mathcal R` represented by
`mathcal R_omega=c+mathcal R`; no larger rooted object is required.

### Corollary AP.4 (all-parity Walsh response quotient)

At every fixed `m`, (AP.10) is an exact quotient for the Boolean maximum of
every graph composition of the linear-label Walsh children

```math
\widehat C_a=D_aFD_a
```

with common Walsh bridges.  It also preserves rooted response functions
after the exposed coordinates are relabelled by the same orthogonal map.

#### Proof

Apply the orthogonal map simultaneously to the two Walsh coordinate
variables.  This is a permutation of the Boolean cube, preserves the Walsh
kernel, and conjugates `D_a` to `D_{Oa}`. `square`

The quotient has an `O(k^2)`-bit presentation independent of `m`: store the
Gram matrix, a basis for the relation kernel, and either “absent” or one
representative of the rooted coset.  It is therefore strictly smaller than
the `km` label presentation when `k=o(m)`.  It is a globally coordinated
quotient, not a componentwise composable one, because separate pieces do not
determine their cross-Gram entries or cross-relations.

## 5. Why the rooted fibre is genuinely necessary in even dimension

### Proposition AP.5 (scalable even-dimensional rooted collision)

Let even `m>=4`, take `a=omega` and `b=e_1+e_2`.  The singleton tuples
`(a)` and `(b)` have the same Gram matrix `(0)` and the same relation kernel
`{0}`, but only the first rooted fibre is nonempty.  They are not in the
same orthogonal orbit, and the rooted Walsh futures separate their projective
responses by at least

```math
\frac16 n^{3/2},\qquad n=2^{2m}.                         \tag{AP.17}
```

#### Proof

Both labels are nonzero isotropic and have no nonzero singleton relation.
Every orthogonal map fixes `omega`, so no such map sends `a` to `b`.

For `c in V`, put

```math
s_c(u,v)=(-1)^{u\cdot v+c\cdot v},qquad y_c=Fs_c.
```

The vectors `y_c` are Boolean.  At query `y_a`, the matched pole gives
response `3n^(3/2)/2`.  For completeness, the crossed bound is a one-line
involutory completion.  After modulating the optimization variable, divide
the crossed objective by `q=sqrt(n)` to obtain

```math
\frac12u^TFu+w^Tu,qquad \|u\|_2^2=n,qquad
w^TFw=0,                                                 \tag{AP.18}
```

where `w=s_{a+b}`; the last identity is Walsh orthogonality of two distinct
linear truth tables.  Since `F^2=I`,

```math
(2I-F)^{-1}=\frac{2I+F}{3}.
```

Completing the square on the containing Euclidean sphere gives

```math
\frac12u^TFu+w^Tu
\le n+\frac12w^T(2I-F)^{-1}w
=\frac43n.                                               \tag{AP.19}
```

Thus the crossed response is at most `4n^(3/2)/3`.  At `y_b` the roles
reverse.  The oscillation of the response difference is at least
`n^(3/2)/3`, hence its projective norm is at least (AP.17). `square`

Thus the rooted fibre is not an artifact of the odd proof.  Given
`(G,mathcal R)`, it is exactly the remaining orbit datum in every dimension,
and omitting it causes a leading-scale response error on an infinite even
subsequence.

## 6. Exact finite regression

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_all_parity_walsh_orbits.py
```

The script enumerates `O(2,2)` and `O(4,2)`, checks the affine-symplectic
parameterization, verifies that the state (AP.10) coincides exactly with
orbits for all ordered tuples through length four at `m=2` and through
length three at `m=4`, and checks the even rooted-collision/resolvent inputs
at `m=4`.  These finite checks are diagnostics; Theorem AP.3 is the general
proof.
