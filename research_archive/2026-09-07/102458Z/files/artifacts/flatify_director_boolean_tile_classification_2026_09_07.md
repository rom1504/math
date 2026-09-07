# Exact classification of unconstrained Boolean bilinear tiles

2026-09-07. Elementary exact theorem, independently reconstructed by the
constructive agent. This explains the rank-two operation's scope; it is
not an impossibility theorem for favorable flatification.

## Statement

Let R be a real p-by-q matrix such that

    u^T R v is in {−1,1} for EVERY u in {−1,1}^p and v in {−1,1}^q.

After deleting zero rows/columns and permuting or signing coordinates,
R is either [1] or H2/2. In particular its rank is at most two.
No assumption on the original child signings appears in this statement.

## Proof

Expand (u^T R v)^2=1 in the orthonormal Boolean characters. The constant
coefficient gives ||R||_F^2=1. Coefficients of u_i u_k and v_j v_l give
orthogonality of distinct rows and of distinct columns. For i!=k,j!=l,
the coefficient of u_i u_k v_j v_l gives

    R_ij R_kl + R_il R_kj = 0.                         (1)

If the nonzero support contains no two entries in different rows AND
columns, its bipartite support graph is a star. Row/column orthogonality
then permits only one entry, of magnitude one.

Otherwise choose a nonzero disjoint pair. Equation (1) forces the other
two corners of that rectangle to be nonzero. No complete 2-by-3 rectangle
can occur: ratios of the two row entries in its three columns would be
pairwise negatives by (1), impossible for three nonzero real numbers.
Similarly no complete 3-by-2 rectangle can occur. Any extra nonzero entry
outside the original 2-by-2 rectangle forces one of these forbidden
rectangles, again by (1). Thus the entire support is that 2-by-2 block.

Write it [[a,b],[c,d]]. Then ad=-bc, ac=-bd, ab=-cd. All four entries are
nonzero. Substituting d=-bc/a in the latter identities gives a^2=b^2=c^2,
and then d^2=a^2. Frobenius normalization gives |a|=|b|=|c|=|d|=1/2;
their sign product is negative. Row/column sign changes yield H2/2.
Both listed forms plainly satisfy the required Boolean identity.

## Meaning for construction, and limits

The previous reciprocal tile uses one coordinate per endpoint. The new
rank-two tile uses the only other possibility with unrestricted independent
Boolean port coordinates. Merely increasing the number of freely varying
binary coordinates in a bilinear tile cannot encode a general child seed.

Higher-rank tiles remain possible if physical row words form CONSTRAINED
codebooks, or if the tile is nonlinear, nonlocal, or changes with the seed.
Hadamard character cosets are one exact example of constrained row words.
Thus this theorem does not exclude the user's global-edge-changing target,
and says nothing about all possible sign-matrix constructions.

This is an explicit low-degree Boolean classification, not a claim to have
invented a new general theory. A related classical perspective is that a
degree-two Boolean function can depend on at most four variables; the
argument above independently proves the sharper bipartite coefficient
classification actually needed here, without an imported theorem.

## Adjacent literature checked, not silently applied

Ghaderpour, *Signed group orthogonal designs and their applications*,
[arXiv:1502.07668](https://arxiv.org/html/1502.07668), Theorems 3.1 and 4.2,
provides signed-group-to-real-design realization and full orthogonal designs
of types (2^j u_1,...,2^j u_k) for sufficiently large j with fixed integer
type tuple. The actual statements and their construction chain were read.
They supply orthogonality and type realization, NOT a bound preserving the
Boolean cap of arbitrary matrix substitutions. Their scalar commuting
variables must not be replaced by arbitrary noncommuting optimal children.
In particular Remark 3 explicitly warns that the circulant signed-group
construction need not work for nonabelian coefficients. No imported
flatification or seed-transfer theorem is obtained from this source.
