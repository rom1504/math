# High-degree flat transport of a bounded trigonometric coherent response

2026-09-06. Verified finite-degree theorem. Independently reconstructed by
the director and audited in
`transfer_adversary_trigonometric_transport_audit_2026_09_06.md`.
This strengthens the global-cut note by directly targeting the missing
high-degree Bc0 cuts. It does not alone close every rich feedback query.

## 1. The finite theorem

Let W_i be a real Boolean polynomial of degree at most M, with any
deterministic mean. Assume every positive-degree coefficient tensor has
all GLOBAL cuts bounded by C log(n+1)^C. The root may lie on either side
of the cut, and the opposite side must contain at least one seed slot.
All degrees are fixed. Let L be a deterministic root transport with
`max_ij |L_ij|<=epsilon_n`. For fixed real t and each fixed q>M, put

```
C_i=exp(it W_i), V=L C.
```

The conclusion is that every proper seed cut of the exact
degree-q Walsh coefficient tensor of V_i is bounded by

```
C_(M,q,t) epsilon_n log(n+1)^C_(M,q,t),                 (1)
```

uniformly in i. The theorem also applies to a finite real trigonometric
polynomial of a fixed vector W whose components have degrees at most M:
apply it to each real linear combination t dot W. No exponential moment
or polynomial approximation in the value of W is assumed.

## 2. Random derivative tensors: every external-slot flattening

For a nonempty set U of s external derivative labels, Delta_U W_i is
a polynomial of degree at most M-s in the remaining random signs.
Consider ANY matrix flattening of its external labels and root, with
at least one external label opposite the root. Its operator norm has
every fixed L^p moment bounded by a fixed power of log(n+1).

Here is a direct proof from the deterministic coefficient cuts.
Fix one residual random degree h. Randomly color the original seed
indices into h classes and retain monomials having exactly one residual
label of each color. Their retention probability is h!/h^h. The
unrestricted residual polynomial equals h^h/h! times the color average
of these restricted sums. Conditional on a coloring, the h sign
families are independent. No external derivative label is reintroduced
as a random label: those entries of the coefficient tensor are zero by
exact squarefreeness.

For a rectangular Rademacher matrix series, the variance parameters are
the norms of its horizontal and vertical coefficient stacks. Its L^p
operator norm is at most a constant times
`sqrt(p+log(d_row+d_col))` times the larger stack norm. This follows by
integrating the primary rectangular tail bound in
[Tropp, User-friendly tail bounds, Theorem 1.5](https://tropp.caltech.edu/papers/Tro11-User-Friendly-preprint.pdf).
That theorem permits deterministic real rectangular matrices and
independent Rademacher signs. Both square-function orientations are
retained; one cannot keep only one of them.

Apply this bound successively to the h color classes, conditioning on
the remaining ones. At each step the new color-label axis is appended
to either the row or the column stack. At the end each of the at most
2^h stack choices is exactly a deterministic coefficient flattening
covered by the assumed all-global-cut bound. Restriction of a color
axis is a coordinate projection. Matrix dimensions are fixed powers
of n, so every logarithmic dimensional factor is O_M(log(n+1)).
Minkowski over the color average and finitely many residual degrees
proves the stated moment bound. For h=0 this is simply the given cut.

This gives MORE than the derivative root-map estimate: some external
derivative slots are explicitly allowed on the root's side. That
additional flexibility is needed in the next section.

## 3. Exact covers force two marked factors when q>M

Use the exact finite cover identity in
`transfer_seed_bounded_trigonometric_global_cuts_2026_09_06.md`:

```
Delta_U C_j=(-1)^q 2^-q chi_U exp(itW_j)
       sum_(F covers U) product_(A in F)
          [exp(it(-2)^|A| chi_A Delta_A W_j)-1].         (2)
```

Here |U|=q, the family F is a subset of the nonempty subsets of U,
and each member occurs at most once. Since W has degree at most M,
every factor with |A|>M is zero. Thus every surviving cover of q>M
has at least TWO nonempty factors.

Fix the signs S for the moment. Each factor in brackets is a rooted
tensor on its external labels A. By Section 2 every one of the
corresponding real derivative tensor's global cuts has polylogarithmic
operator moments. The entrywise exponential bound

```
||exp(itJ)-1||op <= |t| ||J||op + (t^2/2)||J||op^2
```

therefore gives the same moment bounds for every global cut of the
bracketed tensor. Characters chi_A split into signs on the two sides
of any chosen flattening, hence are row/column unitaries. The scale
(-2)^|A| is fixed.

## 4. Flat transport, overlapping covers, and expectation

Choose any proper partition U=U_left union U_right. First give each
cover factor its OWN copies of its labels. Place every copy of a label
on the same side as that label in the chosen partition. There are two
distinct factors with a nonempty left and a nonempty right label,
respectively: if only one factor met both sides, any second nonempty
factor meets one side and can be paired with the first on the other.

Apply the banked two-factor flat-transport factorization to this
disjoint-copy product, grouping the remaining factors as necessary.
The root weight is `L_ij exp(itW_j)`, whose supremum is at most epsilon_n
for every S. The operator bound is epsilon_n times the product of the
relevant global factor cut norms. One uses a genuinely nonempty seed
side in every such global factor cut; no zero-degree root vector is
assigned a bounded Euclidean norm.

Now identify copies of every repeated external label in the cover.
Since ALL copies of a label were put on the same side, these are only
row-diagonal or column-diagonal isometric compressions. They do not
create an uncharged internal trace or a cross-side partial transpose.
The remaining requirement that different external labels be distinct
is a finite sequence of coordinate projections and rectangular
pinching complements; its cost depends only on q. The character chi_U
in (2) also splits into a row and a column sign.

This bounds each proper cut of Delta_U V_i at fixed S by epsilon_n
times finitely many products of random factor norms. Hölder and
Section 2 give a polylogarithmic bound for its expectation. Finally,
the exact Boolean identity

```
E Delta_U V_i = E[V_i chi_U]
```

identifies that expectation with the degree-q Walsh coefficient tensor.
Operator norm is convex, so taking expectation preserves the bound.
This proves (1). The independent audit checks both the random-derivative
matrix-chaos estimate and the overlapping-cover factorization.

## 5. Precise prospective use and limits

For a rich finite old frame, each literal coherent primitive remains
a finite-degree polynomial with all-global coefficient cuts. The
theorem therefore gives small high-degree cuts of Bc0 for each
fixed trigonometric approximation of a bounded c0(W), without assuming
moment determinacy of W. Its averaged L2 approximation can then be
removed INSIDE fixed-degree, variance-capped contraction tests.

It does not yet establish every full-noise contraction into the mean
linear-noise return, the required bounded source approximations for all
other feedback terms, or a uniform positive high-degree source block
for a general function of W. Those assumptions must be separately
checked before claiming a rich full-return cap theorem.
