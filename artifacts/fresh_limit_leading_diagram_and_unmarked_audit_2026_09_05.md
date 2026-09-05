# Exact leading-diagram audit and the scope of unmarked internal vertices

Date: 2026-09-05. This is an adversarial audit of the dependency supporting
the universal lower bound. The finite diagram enumeration found no
counterexample. A separate unmarked-vertex calculation identifies both a
uniform-root obstruction and a limited positive averaged extension.

## 1. Reproducible exact diagram enumeration

The integer-only program
`computations/fresh_limit_tree_pair_partition_audit.cpp` enumerates every
admissible leading pair partition in 29 selected mixed moment and energy
cases, with at most 20 spin slots. Its exact stdout is saved at
`computations/results/fresh_limit_tree_pair_partition_audit.json`.

The ancestor-closed tree family is:

| Tree | Children of its top vertex | Nonroot vertices | Root automorphisms |
|---|---|---:|---:|
| edge | none | 1 | 1 |
| star3 | edge, edge | 3 | 2 |
| asymmetric5 | edge, star3 | 5 | 2 |
| asymmetric7 | edge, asymmetric5 | 7 | 2 |
| balanced7 | star3, star3 | 7 | 8 |

The program enumerated **208,877,244 admissible pair partitions**. Every
case passed all of the following independent checks:

- Empty parity-edge quotients are exactly whole-copy rooted-isomorphism
  pairings for moments.
- For energy, they are exactly a bridge paired to the top edge of one
  `F` tree, whole-copy pairings among the remaining `F` trees, and the
  selected tree's children paired with `H` trees with no child-child pair.
- Their counts agree with a separate weighted Wick/Hermite matching
  recursion using the rooted automorphism numbers.
- Every nonempty moment parity graph is Eulerian and has a free-free edge.
- Every empty parity quotient is a tree whose edges all occur exactly twice.

The largest tests include the fourth moment of asymmetric5 (67,003,200
partitions), the mixed product of asymmetric7, balanced7, asymmetric5, and
edge (27,518,400 partitions, with no survivor), and the energy with three
asymmetric5 inputs against edge and star3 (79,315,200 partitions).

### Independence of the checks

The matching recursion prunes only identifications forbidden by injectivity:
two vertices within one original tree copy cannot pair; in an energy
diagram, the explicit spin `Sj` cannot pair with an internal vertex of an
`H` copy rooted at `j`. It does **not** require a partner of `Sj` to be a
top vertex, and it does not enforce whole-copy pairing during enumeration.

At a completed matching, one routine builds the quotient multigraph and
reduces edge multiplicities modulo two. A different routine checks the
claimed whole-copy/selected-child structure by the explicit matching of
vertices and parent edges. A third recursion predicts the total count from
tree types and automorphisms. Thus the desired classification is not built
into the admissibility test.

Build and run without external dependencies:

```text
c++ -O3 -std=c++17 computations/fresh_limit_tree_pair_partition_audit.cpp -o /home/math/quadra/tmp/fresh_limit_tree_pair_partition_audit
/home/math/quadra/tmp/fresh_limit_tree_pair_partition_audit
```

This is a finite exact audit, not a substitute for the all-orders graph
proof. It found no hostile leading diagram invalidating that proof.

## 2. Unmarking the top vertex of the same three-star

Let `m=n-1`, `B=A/sqrt(m)`, and retain the injectivity exclusions. Remove
only the spin at the top internal vertex of star3, keeping its two leaves
marked. The resulting field is

\[
 Z_i=m^{-3/2}\sum_{\substack{j,k,l\ne i\\j,k,l\text{ distinct}}}
                  a_{ij}a_{jk}a_{jl}S_kS_l.            \tag{1}
\]

The normalization here is the unnormalized-tree convention; the marked
star3 has limiting variance two. This is the even `He2` unmarked feature,
not the odd, linear-orthogonal `He3` feature used in the separate AMP note.

Suppose the root row is all `+1`, and write

\[
 A=\begin{pmatrix}0&\mathbf1^T\\\mathbf1&C\end{pmatrix},
 \qquad C\in\{0,\pm1\}^{m\times m}\text{ hollow symmetric}.
\]

Then the exact field and variance are

\[
 Z_i=m^{-3/2}\sum_{k\ne l}(C^2)_{kl}S_kS_l,
 \qquad
 \boxed{\operatorname{Var}Z_i={2\over m^3}
           [\operatorname{Tr}C^4-m(m-1)^2].}           \tag{2}
\]

Thus the missing internal spin exposes a fourth-moment matrix invariant
that the marked-tree theorem did not need.

### Explicit flat Hadamard family: variance tends to zero

Let `n=2^r`, let `H_n` be the symmetric Sylvester Hadamard matrix with its
first row all `+1`, let `D=diag(H_n)`, and set `A=H_n-D`. This is an actual
hollow signing. The root is its first coordinate. Its normalized operator
norm tends to one, and

\[
 \|B^2-I\|_{op}\le {2+2\sqrt n\over n-1}\to0.
\]

In particular its Gram off-diagonal coherence tends to zero. A direct
entrywise sum using `sum diag(H_n)=0` gives

\[
 \sum_{k\ne l}(C^2)_{kl}^2=3(m-1)(m-2),
 \qquad
 \operatorname{Var}Z_i={6(m-1)(m-2)\over m^3}\to0.     \tag{3}
\]

For clarity, if `d_k=(H_n)kk`, then, off the removed root and off the
diagonal, `(C²)kl=-1-(H_n)kl(d_k+d_l)`. Expanding its square, using
`sum_(k≠root)d_k=-1` and the zero nonfirst Hadamard row sums, yields (3).

### Actual random-sign cores: variance tends to two

Instead choose `C` with independent symmetric signs above its diagonal.
Its fourth trace has the exact mean and variance

\[
 E\operatorname{Tr}C^4=m(m-1)(2m-3),\qquad
 \operatorname{Var}(\operatorname{Tr}C^4)=192\binom m4.
\]

Indeed its random part is eight times the sum of the distinct unoriented
four-cycle edge monomials; these monomials are orthogonal. Hence (2) tends
to two in probability. Simultaneously, there exist deterministic choices
of these signings with uniformly bounded `||B||op` and Gram coherence
`O(sqrt(log n/n))`: a fixed epsilon-net and Hoeffding bound give
`||C||op≤8sqrt(m)` with probability tending to one, while entrywise
Hoeffding bounds control the row sums and off-diagonal entries of `C²`.
Intersecting these events with the fourth-trace concentration event gives
actual deterministic signing sequences with all three properties.

The two families disprove a state-free **uniform-in-root** extension of
the marked-tree law, even under bounded operator norm and vanishing Gram
coherence. The distinguished Hadamard root is exceptional; this argument
does not disprove an averaged-root extension.

## 3. A positive diagonal-state theorem for this one quadratic field

For arbitrary `A`, put `Q=B²`, `r=B1`, and

\[
 M_i=B\operatorname{diag}(B_i)B.
\]

Let `N_i` be `M_i` with its diagonal and its `i`th row and column removed.
Then (1) is exactly `S^T N_i S`, and its variance is
`v_i=2||N_i||_F²`. If `||B||op≤L` uniformly, then

\[
 \sup_i\|N_i\|_{op}=O_L(n^{-1/2}),\qquad
 \sup_i\|N_i\|_F=O_L(1).                            \tag{4}
\]

The first bound follows from
`||B diag(B_i) B||op≤L²/sqrt(m)` and
`(M_i)kk=(r_i-Bik)/m`; removing a row/column is an orthogonal compression.
The second follows from
`||B diag(B_i)B||_F≤L||diag(B_i)B||_F=L`.

For every fixed real `t`, uniformly in `i`,

\[
 E e^{itZ_i}-e^{-t^2v_i/2}\longrightarrow0.            \tag{5}
\]

A self-contained proof replaces Rademacher coordinates by independent
Gaussians one at a time. The first three moments match. Fourth-order
Taylor errors total
`O_t(sum_k ||(N_i)_k||₂^4)=O_(t,L)(n^-1)` by (4). For the Gaussian
quadratic form, diagonalize `N_i`; its eigenvalues have maximum absolute
value `O_L(n^-1/2)`, bounded squared sum, and sum zero. Expanding its
characteristic function gives (5). This permits root-dependent or
degenerate variances, including (3).

The exact root variance has the explicit state formula

\[
 {v_i\over2}=1+[B(Q^{\circ2}-I)B]_{ii}
 -{(m-3)r_i^2+2\|r\|^2-4(Q\mathbf1)_i+3\over m^2}. \tag{6}
\]

It follows by subtracting the diagonal and the `i`th off-diagonal row
norms from `||M_i||_F²=[B Q^(circ2) B]ii`.

If additionally `eta_n=max_(j≠k)|Qjk|→0`, then

\[
 {1\over n}\sum_i|v_i-2|\to0.                        \tag{7}
\]

To see this, `K=Q^(circ2)-I` satisfies
`||K||_F²≤eta_n² Tr Q²≤eta_n² L²n`. Consequently the mean squared
diagonal of `BKB` is at most `L^6 eta_n²`. All row-sum terms in (6) have
vanishing averaged absolute value, using `||r||²≤L²n`. Combining (5)--(7)
gives the averaged-root standard Gaussian law `N(0,2)` for this field.

Thus bounded operator norm and coherence do support one explicit
unmarked quadratic extension with diagonal state, and an averaged standard
law. This does not assert a general multilevel unmarked AMP theorem or
remove the stronger all-power hypothesis of the separate imported result.
