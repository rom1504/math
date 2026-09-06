# All global Walsh cuts of fixed centered polynomial computations

Date: 2026-09-06. Status: exact algebraic proof independently reconstructed
with the director and audit agent. This is a structural operator estimate,
not a feedback energy or Gaussian-universality theorem.

## 1. Statement and canonical application

Start with a fixed number of independent Boolean seed vectors. Form a
fixed finite computation using deterministic bounded-op matrix transports
and fixed coordinate polynomial maps. Assume every input to a matrix
transport is coordinatewise mean zero; global oddness suffices. Dimensions
are bounded multiples of `n`, and depth and polynomial degrees are fixed.

For every intermediate vector `C`, every positive Boolean Walsh degree
`q`, and EVERY partition of its `q` marked input slots with at least
one slot on the side opposite the output root, the corresponding global
kernel flattening has operator norm at most

```math
C_0(\log(n+1))^{C_1}.                                      (1)
```

The output root can be placed on either side, by transposition. Constants
may depend on the fixed computation, not on matrix order. Exact source
Walsh projections, including all internal Boolean contractions, are kept.

In particular this applies to fixed coherent polynomials of

```math
S,\quad G=BS,\quad QS,\quad D=S h_2(G),\quad Y=BD,\quad V=QD.
```

The prior gradient argument controlled only the cut `root | all inputs`.
The present statement also permits arbitrary marked inputs on the output
root's side. It does not claim a SMALL proper cut for the literal return
`QD`: those cuts may remain of order one.

## 2. Higher exact Boolean derivative matrices

For a set `U` of distinct seed labels define

```math
\Delta_U=\prod_{a\in U}\Delta_a,
\qquad \Delta_a F=\frac{F(S)-F(S^{(a)})}{2S_a}.
```

Let `J_C^(s)(S)` be the matrix with output rows `i` and columns the
increasing `s`-tuples `U`, whose entry is `Delta_U C_i(S)`.
For `s=0`, coordinate values are treated as row-diagonal factors rather
than a one-column matrix. We claim that for every fixed `s>=1` and
fixed moment order `p`,

```math
\left\|\|J_C^{(s)}\|_{op}\right\|_{L^p}
 \le C_{p,s}(\log(n+1))^{C_{p,s}}.                           (2)
```

Orders above the fixed Boolean degree vanish. Transport by `M` gives
`J_(MC)^(s)=M J_C^(s)` exactly. It remains to check a coordinate
polynomial product.

## 3. Exact cover expansion for a product

By linearity it suffices to consider `C_i=prod_{ell=1}^d W_i^ell`,
where factors can repeat. For any `T`, the exact multilinear expansion
under seed flips is

```math
W_i(S^T)=\sum_{V\subseteq T}(-2)^{|V|}\chi_V(S)\Delta_V W_i(S).
```

Also

```math
\Delta_U C_i(S)
 =2^{-|U|}\chi_U(S)\sum_{T\subseteq U}(-1)^{|T|}C_i(S^T).
```

Substitution and the alternating subset sum retain exactly the covers
`V_1 union ... union V_d=U`. Thus

```math
\Delta_U C_i
 =\sum_{\substack{V_1,\ldots,V_d\subseteq U\\\cup_\ell V_\ell=U}}
 (-1)^{|U|+\sum_\ell|V_\ell|}
 2^{\sum_\ell|V_\ell|-|U|}
 \chi_U\prod_\ell\chi_{V_\ell}\,
 \prod_\ell\Delta_{V_\ell}W_i^\ell.                         (3)
```

For fixed `s=|U|` and `d`, there are only finitely many relative cover
patterns. Empty `V_l` supply ordinary coordinate row factors. For each
nonempty pattern, take the tensor product of the corresponding derivative
matrices `J_(W^l)^(|V_l|)`. Compress its output rows to equal root indices
`i_1=...=i_d`: this is a coordinate row compression of norm one. Embed the
full distinct `U`-tuple into the tensor-product column tuples prescribed
by the `V_l`. Repeated seed labels between two factors are thus handled
by a column-diagonal isometric embedding. This embedding is injective
because the subsets COVER `U`; no differentiated seed is unused.

The accumulated character `chi_U prod_l chi_Vl` is a sign depending
only on the full `U` column. It is a right diagonal of operator norm
one. Therefore each term of (3) has operator norm bounded by

```math
C_{s,d}\prod_{\ell:V_\ell=\varnothing}\max_i|W_i^\ell|
         \prod_{\ell:V_\ell\ne\varnothing}
                    \|J_{W^\ell}^{(|V_\ell|)}\|_{op}.         (4)
```

There are no uncontrolled row/column broadcasts or diagonal-removal
approximations in this factorization.

The previously checked first-gradient/Poincare induction supplies
polylogarithmic coordinate maximal moments at every operation: polynomial
maps use fixed-degree hypercontractivity; a centered transport has row
variance bounded by its operator norm squared times the preceding
covariance operator bound. The latter follows from
`Cov(W) <= E J_W^(1)(J_W^(1))^T`. Applying Holder to (4) therefore
propagates (2) simultaneously for all finitely many derivative orders.
The initial seed vectors have first derivative matrices of norm one
and all higher derivatives zero.

## 4. Fourier orthogonality identifies every global cut

Write the exact degree-`q` component as

```math
C_{q,i}=\sum_{|T|=q}c_{i,T}\chi_T(S).
```

For `1<=s<=q`, the degree-`q−s` random component of `J_C^(s)` has
entry

```math
[J_C^{(s)}]_{q-s}(i,U)
 =\sum_{\substack{|A|=q-s\\A\cap U=\varnothing}}
                c_{i,A\cup U}\chi_A(S).
```

For any deterministic vector `v` on the `s`-tuple columns,

```math
\mathbb E\|[J_C^{(s)}]_{q-s}v\|_2^2
 =\sum_{i,A}\left|\sum_{U\cap A=\varnothing}
                    c_{i,A\cup U}v_U\right|^2.              (5)
```

The right side is precisely the squared norm of the global Walsh
flattening with row `(output i, remaining marked set A)` and column
`differentiated set U`, up to fixed factorials if symmetric ordered
tensor conventions are used. Orthogonality of random Walsh degrees
gives

```math
\mathbb E\|[J_C^{(s)}]_{q-s}v\|_2^2
 \le\mathbb E\|J_C^{(s)}v\|_2^2
 \le\mathbb E\|J_C^{(s)}\|_{op}^2\|v\|_2^2.
```

Together with (2), this proves (1) for every possible global cut.

## 5. What has and has not changed

These bounds tolerate every exact Boolean diagonal contraction generated
inside a fixed centered polynomial computation. They do not supply an
extra vanishing factor by themselves. The director's subsequent two-factor
flat-transport factorization uses the new global cuts to produce that
factor. Combined with local collision control and the separate source
split/Hall arguments, it enters the full first-marked-history energy proof
in `continued_feedback_first_marked_history_energy_projection_2026_09_06.md`.
That energy conclusion is not asserted as a consequence of global cuts
alone.

The centering restriction remains essential for the proof. Transporting
the constant vector through an actual apex signing creates a coordinate
of size `sqrt(n)`; multiplying that coordinate by its own seed then
creates a positive Walsh root map of size `sqrt(n)`. No such uncentered
transport is covered by (1).
