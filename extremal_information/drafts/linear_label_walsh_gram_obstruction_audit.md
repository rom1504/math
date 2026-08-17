# Adversarial audit: linear-label Walsh Gram obstruction

**Verdict: REPAIR.**  The algebraic counterexample in Theorem LG.1, including
its `n^(3/2)` normalization and constant gap, is correct.  The odd-dimensional
orbit classification in LG.3 is also correct.  However, the draft materially
overstates which previously considered state has been falsified: the binary
label Gram matrix is **not** stronger than pairwise truth-table-overlap data.
Indeed, the displayed collision is immediately separated by the latter.  The
result is promotable after narrowing that comparison and stating explicitly
that only the binary bilinear Gram/commutation quotient is rejected.

## 1. Independent normalization and dimension check

Let `q=2^m`, `n=q^2`, let `R` be the order-`q` Walsh matrix, and put

```math
W=R\otimes R,\qquad F=W/q.
```

Thus `F` is an `n by n` symmetric orthogonal involution on functions of
`(u,v) in F_2^m times F_2^m`.  For `a in F_2^m`, `D_a` is the `n by n`
diagonal with entry `(-1)^(a dot v)`, and

```math
\widehat C_a=D_aFD_a
```

is again a symmetric orthogonal involution.  The unnormalised child matrix is
`C_a=q\widehat C_a`; hence its quadratic contribution is
`(q/2)x^T\widehat C_ax`.  The bridge matrix is `W=qF`, hence its contribution
is `q x^TFy`.  A Boolean block has squared norm `n`, and `q n=n^(3/2)`.
These are exactly the factors in (LG.6).

The retained diagonal of `C_a` causes no error: its total trace is zero for
`m>=1`, so diagonal deletion changes the Boolean quadratic by the constant
zero.  The operator bound and a Boolean pole give child cap
`n^(3/2)/2`.

## 2. The Gram collision is exact, but its scope is narrower than claimed

With

```math
a=(1,1,1,0,...,0),\qquad b=(0,0,1,0,...,0),
```

one has over `F_2`

```math
a dot a=b dot b=a dot b=1.
```

Therefore the full ordered `3 by 3` binary Gram matrices of `(a,a,a)` and
`(a,b,a)` are both all-one matrices.  Also `a+b` is nonzero, isotropic, and
orthogonal to both `a` and `b`.  The tuples cannot be related by a global
orthogonal relabelling: their linear spans have dimensions one and two,
respectively.  Thus there is no accidental orbit equivalence hidden in the
example.

There is, however, an important scope error in Sections 1 and 6.  For linear
truth tables,

```math
S(g_c,g_d)=\sum_v(-1)^((c+d) dot v)
=q\,1_{c=d}.
```

Consequently pairwise truth-table overlaps distinguish the words:
`(a,a,a)` has overlap `q` on every off-diagonal pair, whereas `(a,b,a)` has
zero overlap on the pairs crossing `a` and `b`.  Thus binary Gram data is not
"stronger than the bias and truth-table-overlap data" used earlier.  It is a
different invariant, and in this example it has forgotten the elementary
equality/relation pattern retained by truth-table overlap.

The exact falsified state is therefore

```math
(c_1,...,c_k)\longmapsto(c_i dot c_j)_{i,j},
```

used alone as a quotient for the fixed Walsh-bridge Boolean extremum.  The
theorem does **not** falsify biases plus pairwise truth-table overlaps, Gram
plus equality data, Gram plus the relation kernel, or the full linear-label
presentation.

## 3. Bent witness and good-word saturation

Adding `a dot v=v_0+v_1+v_2` to `Q_m` leaves the disjoint quadratic pairs

```math
(u_0,u_1),\ (v_0,v_1),\ (u_2,v_2),\
(u_j,v_j)\quad(3\le j<m).
```

They partition all `2m` variables.  The sign vector of `rs` is a `+1`
eigenvector of the normalized order-four Walsh transform.  Hence

```math
F(D_ax)=D_ax,\qquad \widehat C_ax=x.
```

Adding `b dot v=v_2` instead leaves exactly one factor
`rs+r+s`, whose sign vector is a `-1` eigenvector; all other factors have
sign `+1`.  Hence

```math
F(D_bx)=-D_bx,\qquad \widehat C_bx=-x.
```

The polar form pairs all `2m` coordinates nondegenerately, so `x` is bent
and `y=Fx` is Boolean.  Equivalently this follows by tensoring the displayed
two-variable transforms.

Walsh modulation/translation gives

```math
F\widehat C_b=-\widehat C_bF
```

because `b dot b=1`.  Therefore

```math
\widehat C_b(Fx)=-F(\widehat C_bx)=Fx.
```

The triple `(x,Fx,x)` has child Rayleigh value `n` at all three sites and
bridge bilinear value `n` on both edges.  Its energy is therefore

```math
3(qn/2)+2(qn)=7n^(3/2)/2.
```

Separate operator-norm bounds give the same upper bound, so (LG.14) is exact
for every `m>=3`.

## 4. Constant odd-word bound and gap

For the word `(a,a,a)`, `a dot a=1`, so `F` and `\widehat C_a`
anticommute.  With `A=A(P_3)`,

```math
\mathcal M=I_3\otimes\widehat C_a+A\otimes F
```

satisfies

```math
\mathcal M^2=(I_3+A^2)\otimes I_n.
```

Since `||A(P_3)||=sqrt(2)`, `||\mathcal M||=sqrt(3)`.  A three-block Boolean
vector has squared norm `3n`, and the physical energy is
`(q/2)X^T\mathcal MX`; hence

```math
\max E_{(a,a,a)}\le (3\sqrt3/2)n^(3/2).
```

Subtracting from the exact good-word value gives

```math
((7-3\sqrt3)/2)n^(3/2)
=0.901923...\,n^(3/2),
```

as stated.  No absolute-value or factor-of-two loss is hidden here.

## 5. Audit of the added rooted orbit quotient LG.3

For odd `m`, the characteristic vector `omega=(1,...,1)` has
`omega dot omega=1`, and `H=omega^perp` is a nondegenerate alternating
space.  Every vector has the unique decomposition

```math
u=(u dot u)omega+h(u),\qquad h(u)\in H.
```

Equality of relation kernels makes the label map well defined; Gram equality
makes it an isometry; equality of the rooted relation coset is exactly what
ensures that membership of `omega` agrees and that `omega` is fixed when it
lies in the label span.  The induced map on the corresponding subspaces of
`H` is well defined and symplectic.  The symplectic Witt extension lemma then
extends it to `S in Sp(H)`, and `O(c omega+h)=c omega+Sh` is the required
orthogonal extension.  Conversely every orthogonal map fixes the unique
characteristic vector `omega`, so it preserves all three data.

Simultaneously permuting `(u,v)` by `(Ou,Ov)` preserves the Walsh kernel and
sends `D_a` to `D_{Oa}`.  Thus the extremal conjugacy claim follows.  An
independent exhaustive check at `m=3`, for all ordered tuples of lengths
one, two, and three, found that the `(Gram,R,R_omega)` classes agree exactly
with the orbits of all six orthogonal maps.  This computation is diagnostic;
the proof above is the general argument.

The response-function statement is correctly qualified as equivariance
under the same coordinate relabelling.  It would not give equality for
coordinate-labelled external queries left unrelabeled.

## 6. Verifier audit

Running

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_linear_label_walsh_gram.py
```

returns 181 passed checks, including 144 rooted orbit classes at `m=3`.  The implementation consistently uses
`ca=q\widehat C_a` and `w=qF`; its variable named `doubled_energy` is exactly
twice (LG.6), so the assertion `7*n*q` corresponds to energy
`7n^(3/2)/2`.  The full squared block identity is checked at `m=3`; the
factor identities from which it follows are checked also at `m=4`.

The verifier also checks the characteristic-root collision and exhaustively
compares `(Gram,R,R_omega)` classes with orthogonal-group orbits for every
ordered tuple of lengths one, two, and three at `m=3`.  This finite test is
diagnostic; LG.3 has a dimension-uniform proof.

## 7. Required repairs before promotion

1. Replace the claim that binary Gram data is stronger than pairwise
   truth-table overlaps.  Explicitly display the overlap matrices above.
2. Describe the result as complementary to, not a strengthening of, the
   constant-label overlap/commutation obstruction unless a precise stronger
   implication is supplied.
3. State in the theorem interpretation that the quotient rejected is Gram
   alone; equality/relation data already separates this example.
4. The optional small `m=3` orbit-classification test has now been added.

After repairs 1--3, the rigorous LG.1 counterexample and LG.3 orbit quotient
are suitable for promotion.
