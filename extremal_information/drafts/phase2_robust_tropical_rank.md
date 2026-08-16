# Robust extremal feature rank from tropical fooling sets

**Status:** director theorem draft; not yet promoted.  The exact code-rank
theorem used in Section 3 is imported from Sheshadri and was already
reconstructed in the main project ledger.  The robust approximation theorem
below is proved here.

## 1. Approximate min-plus rank

For a finite real matrix `M`, let `rank_min(M)` be the least `r` for which

```math
M(x,y)=\min_{1\le t\le r}\{u_t(x)+v_t(y)\}.          \tag{TR.1}
```

Define its uniform `epsilon`-approximate min-plus rank by

```math
\operatorname{rank}_{\min,+}^{(\epsilon)}(M)
=\min\{\operatorname{rank}_{\min,+}(\widetilde M):
       \|M-\widetilde M\|_\infty\le\epsilon\}.       \tag{TR.2}
```

Negating every entry gives the identical statement for max-plus response
kernels.

## 2. A robust tropical fooling-set theorem

### Theorem TR.1

Suppose `M` has `r` distinguished cells `(x_i,y_i)` and a number `G>0` such
that, for every `i!=j`,

```math
M(x_i,y_j)+M(x_j,y_i)
-M(x_i,y_i)-M(x_j,y_j)\ge G.                        \tag{TR.3}
```

Then every `Mtilde` with

```math
\|M-\widetilde M\|_\infty<G/4                      \tag{TR.4}
```

has min-plus factorization rank at least `r`.

#### Proof

In a factorization (TR.1), every term majorizes the represented matrix
entrywise, because the entry is the minimum of all terms.  At every
distinguished cell choose one term attaining the minimum.  If one term `t`
were tight at both cells `i` and `j`, separability would give

```math
\begin{aligned}
\widetilde M(x_i,y_i)+\widetilde M(x_j,y_j)
&=u_t(x_i)+v_t(y_i)+u_t(x_j)+v_t(y_j)\\
&=u_t(x_i)+v_t(y_j)+u_t(x_j)+v_t(y_i)\\
&\ge\widetilde M(x_i,y_j)+\widetilde M(x_j,y_i).
\end{aligned}                                       \tag{TR.5}
```

On the other hand, if the uniform error is `epsilon`, (TR.3) implies that
the right side of (TR.5) minus the left side is at least `G-4epsilon>0`, a
contradiction.  Hence the `r` cells require `r` distinct tight terms.
`square`

The proof is a robust version of the crossing/fooling-set lower bound.  It is
not a statement about computational hardness or the number of bits in an
arbitrary nonlinear encoding.

### Sharp scope

- The theorem controls uniform entry error.  Average error can hide a small
  distinguished block and requires an additional mass hypothesis.
- The strict threshold is necessary for this proof.  At `epsilon=G/4` the
  two sides of (TR.5) can tie.
- The hypothesis concerns a joint four-cell cancellation gap.  Bounding rows
  or columns separately is neither assumed nor used.

## 3. Linear-code conditional responses

Let `C<=F_2^m`, split the coordinates as `L disjoint_union R`, and put

```math
W(x_L,x_R)=d((x_L,x_R),C),
\qquad
s=\dim C-\dim C_L-\dim C_R.                         \tag{TR.6}
```

Sheshadri's exact theorem states

```math
\operatorname{rank}_{\min,+}(W)
=\operatorname{rank}_{\rm trop}(W)=2^s.             \tag{TR.7}
```

The proof selects one lifted codeword from every class of
`P_R(C)/C_R`.  The resulting `2^s` by `2^s` submatrix has diagonal zero and
off-diagonal entries at least one.

### Corollary TR.2 (robust trellis-state lower bound)

For every `0<=epsilon<1/2`,

```math
\boxed{
\operatorname{rank}_{\min,+}^{(\epsilon)}(W)=2^s.}
                                                               \tag{TR.8}
```

#### Proof

The transversal block satisfies (TR.3) with `G=2`, so Theorem TR.1 gives the
lower bound.  Taking `Wtilde=W` and using (TR.7) gives the upper bound.
`square`

Thus sub-half-unit accuracy cannot reduce the optimal trellis channel count:
every admissible approximant has rank at least `2^s`, and the unperturbed
table attains that rank.  The approximant itself may of course use more
terms.  This remains a lattice-scale statement: after normalizing all
distances by block length `m`, the protected error is only `1/(2m)`.  This
answers the small-uniform-error regime of the approximate-factorization
question posed in Remark 7 of the source paper.  It does not answer errors
that grow with block length, relative error, or average error.

The primary source for (TR.7) and the transversal block is Karthik Sheshadri,
[*Trellis State Complexity as an Exact Tropical Factorization
Rank*](https://arxiv.org/abs/2607.23471), arXiv:2607.23471v1 (2026), Theorem
1 and Lemmas 2--3.  The four-page proof was independently reconstructed for
this draft rather than inferred from the abstract.

## 4. Relation to extremal information

The theorem supplies a rigorous growth law for one query-generated feature
algebra:

```math
\text{linear-code state dimension }s
\quad\longmapsto\quad
\text{minimum exact or sub-half-error tropical channels }2^s.
```

This state is much smaller than the full conditional table, whose dimensions
are `2^|L|` by `2^|R|`, whenever `s<<min(|L|,|R|)`.  Conversely, families
with `s=Theta(m)` prove exponential feature-channel growth inside a restricted
linear-code class, not merely inside the universal class of arbitrary
kernels.

The conclusion is deliberately algebra-specific.  A planar signed-graph
instance can have exponential tropical rank while its ground state is
polynomial-time computable.  Tropical response rank measures the number of
separable optimization channels crossing the declared interface; posterior
width measures mutual information about a random latent instance.  Neither
quantity generally bounds the other.

### Proposition TR.3 (average error erases a full-rank resonance)

Let `D_r` be the `r` by `r` matrix with zero diagonal and every off-diagonal
entry equal to one.  Then

```math
\operatorname{rank}_{\min,+}^{(\epsilon)}(D_r)=r
\quad(0\le\epsilon<1/2),                            \tag{TR.9}
```

but the rank-one all-one matrix `J_r` has normalized Frobenius distortion

```math
{1\over r^2}\|D_r-J_r\|_F^2={1\over r}\longrightarrow0. \tag{TR.10}
```

#### Proof

The `r` diagonal cells satisfy Theorem TR.1 with `G=2`, and the standard
`r`-term representation proves (TR.9).  The matrices `D_r` and `J_r` differ
only on the `r` diagonal cells, proving (TR.10). `square`

Thus high exact tropical response rank may be carried by a zero-density set
of exposed anchors.  No average-error rank theorem can follow from exact rank
alone.  It must additionally control the query mass of the cells certifying
separation.  This is the same structural obstruction seen in the rare-fibre
synchronization example: an averaged statistic can converge while a declared
zero-temperature query keeps exposing the exceptional fibre.

## 5. Minimal next theorem

The robust theorem stops at the lattice scale, and Proposition TR.3 prevents
a distribution-free extension.  The next genuinely new target is a
**massive tropical fooling-set theorem**: impose an anti-rare-anchor condition
under a declared root measure and prove a lower bound on approximate
factorization under mean-square rather than uniform response error.  A useful
formulation must be stable under composition and nonvacuous on at least one
restricted code or CSP family.  Such a result could connect algebraic channel
growth to the posterior-width rate theorem instead of leaving the two
complexity measures in parallel.
