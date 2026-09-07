# Independent audit of actual optimizer noise response

2026-09-07. Source:
`flatify_construct_2026_09_07_optimizer_noise_response.md`. Verdict: PASS.

After ground switching, single-spin optimality gives nonnegative fields
and their sum is 2Q. Expanding independent biased signs into their means
plus centered variables leaves orthogonal linear and edge monomials,
giving the source's exact variance without Gaussian assumptions. Fixed-
size noise has edge-product mean ((n-2s)^2-n)/(n(n-1)), as stated.

The universal-vertex construction has exactly the claimed cap: triangle
inequality supplies the upper bound, and the common positive ground
center attains every new positive edge. Restriction is cap-contractive
by averaging outside spins, so the resulting matrix is within O(rn) of
the true order-n minimum without assuming convergence of M_n/n^(3/2).
Its new vertices have degree n-1. Choosing r=o(sqrt(n)) but larger than
n^epsilon proves the stated failure of an improved variance exponent
for the near-optimal class. It is not an exact-minimizer counterexample.

The identity averaging Var(X+Y) and Var(X-Y) cancels covariance exactly.
The explicit row-sum formula includes all independent within-child and
bridge quadratic monomials, with no missing cross covariance. Reversing
all spins in one child realizes the other sign of Y while preserving
both internal energies. Thus a one-orientation linear-field cancellation
cannot provide simultaneous cancellation in that fixed-child parent.

The center lower bound Q(parent)>=Q(A)+Q(D)+|z| requires the selected
centers to have positive cap-attaining energies, exactly as assumed in
the source's Section 3. When selecting actual optimal children, their
global edge polarities can arrange that setup. It should not be applied
to arbitrary fixed children whose absolute maxima have incompatible
energy polarities without checking that assumption.

Finally, these bridge-cancellation conclusions keep the internal child
edges fixed. They are not an obstruction to the broader favorable
flatification target, which permits global changes to those edges.
