# Independent audit: high-degree transport of a bounded trig response

2026-09-06. **PASS for the finite theorem** in
`transfer_seed_trigonometric_high_degree_transport_2026_09_06.md`.
This verifies its new proper-cut conclusion, not merely the preceding
global-cut lemma. No full rich feedback or positive-gain theorem is
inferred from this step alone.

## 1. Matrix-chaos step with literal Boolean derivatives

Fix an external derivative set U and one residual Walsh degree h.
The random derivative tensor has coefficients

```
Delta_U W_i = sum_(V disjoint U, |V|=h) c_(i,U union V) chi_V
```

at that residual degree, up to fixed tensor-normalization factors.
Coloring ONLY the residual labels is sufficient. The probability that
the h distinct labels of V occupy all h classes is h!/h^h, independently
of U. Conditional on the coloring, the resulting matrix chaos has h
independent Rademacher families. Its coefficient matrices are deterministic:
entries in which a residual label equals an external label are already
zero in the exact diagonal-free coefficient tensor. External labels do
not need independent new random signs.

The primary rectangular-series ingredient was checked directly in
[Tropp's author preprint, Corollary 4.2, printed p.15](https://tropp.caltech.edu/papers/Tro11-User-Friendly-preprint.pdf).
For independent Rademacher signs and deterministic rectangular matrices,
its variance parameter includes both horizontal and vertical square
functions; its tail prefactor is the sum of the two dimensions. Tail
integration gives an Lp operator bound with a factor of order
sqrt(p+log(d_row+d_col)).

Condition on all but one residual color. Apply that bound and retain BOTH
block-stack orientations. Each orientation appends that color's label
axis to one side of the flattening. Repeat for all h colors. Every final
stack is one of the deterministic coefficient global cuts in the stated
hypothesis. The external side opposite the root was nonempty from the
start and remains so. Dimensions are fixed powers of n, and there are
only 2^h stack choices. Conditional Lp estimates, Minkowski, and finite
color averaging therefore give exactly the required polylogarithmic
moments of every random-derivative external-slot cut.

There is no dimension-free claim for these random derivative norms.
The logarithms are harmless only because the eventual root transport
supplies a vanishing entrywise factor.

## 2. The entrywise exponential and overlap factorization

For a real matrix J, the remainder of exp(itJ)-1-itJ is bounded entrywise
by t^2 J_ab^2/2. Absolute row and column sums are thus bounded by
`t^2 ||J||op^2/2`; the Schur test proves the source's quadratic matrix
bound. This applies separately to EVERY chosen flattening of the real
random derivative tensor. The seed character chi_A splits into row and
column diagonal signs. No exponential of a random operator norm is
introduced.

For derivative order q>M, every nonzero factor in the exact finite cover
formula has |A|<=M. A surviving cover of q labels therefore has at least
two factors. Fix a proper partition of the q original external labels
into left and right. Lift the cover product by giving each factor its
own copies of the labels. Put EVERY copy of a given original label on
that original label's prescribed side.

The lifted product satisfies the two-factor flat-transport theorem.
Its root weight is L_ij exp(itW_j), with absolute value at most
epsilon_n. Each factor uses a genuinely proper global cut; the standard
root-side assignment handles whole and straddling factors. Extra factors
are retained in that same product factorization.

Now restore all label identifications in the cover. Every identification
is wholly within the left side or wholly within the right side. These
are norm-one diagonal coordinate embeddings/compressions. They neither
trace over an internal surviving index nor identify a left index with
a right index. The last distinctness requirement BETWEEN different
original labels is handled by fixed-degree pinching complements and
coordinate projections.

Thus the factorization supplies epsilon_n times a finite product of
the global random-factor norms. Holder makes its expectation
`epsilon_n polylog(n)`.

## 3. Exact Walsh projection and conclusion

For every set U of q distinct labels,

```
E Delta_U (L exp(itW))_i = E[(L exp(itW))_i chi_U].
```

This identifies the EXPECTED derivative tensor with the exact
degree-q Walsh coefficient tensor. The derivative identity holds for
the entire bounded exponential on the finite cube; no truncation of
its infinite formal Taylor series in W is made. Convexity of operator
norm passes the preceding bound through expectation.

Consequently, at every fixed q>M, each proper seed cut of the exact
degree-q projection of L exp(itW) is
`O(epsilon_n polylog(n))`. Finite real trigonometric sums follow by
linearity. This verifies the proposed finite theorem.

Bounded-response approximation can be used subsequently only inside
an averaged-L2-continuous, fixed-degree, capped contraction argument.
It cannot convert arbitrary averaged L2 convergence into uniform
operator-cut convergence. The mean-linear-noise-return contractions
and existence of a positive high-degree block remain distinct
obligations in a richer feedback theorem.

No upstream file was downloaded into the workspace for this audit;
the primary theorem was read through the web tool. The independently
authored reconstruction and durable primary URL are preserved here.
