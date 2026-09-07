# Independent audit: constant two-frame actual cap obstruction

2026-09-07. **PASS.** Complete read and independent reconstruction of
`principle_synthesis_2026_09_07_constant_two_frame_obstruction.md`.

## Exact conditioning

Conditional on the iid opposite child D, each unordered edge retains one
independent fair orientation variable C_ij, with
`C_ji=S_ij D_ij C_ij`. For any chosen y depending only on D, define x from
the hollow bridge fields with independent unbiased tie coins. The two row
fields at i,j share only this orientation variable. Their remaining sums
are independent sums of n-2 fair signs. If lambda_n is their central mass
(at zero for even n and at one for odd n), the conditional sign means
are lambda_n C_ij y_j and lambda_n S_ij D_ij C_ij y_i. Therefore

    E H_S(x)=lambda_n^2 H_D(y),
    E x^T C y=n E|sum_(1)^(n-1) epsilon_j|.

Both identities are exact before the n matching edges are filled. Those
edges are then paid by at most n; they must not be inserted into the
fields while retaining the same exact covariance formula.

## Polynomial Gaussian witness

For each fixed polynomial and positive regularizer, the semicircle walk
calculations in Section 3 have the required uniform integrability. In a
nonzero term of four walks from distinct i to j, a cyclic union has at
most half the total walk length in vertices. In a tree, the unique i--j
path is traversed by all four walks, giving the same bound. This proves
the O(1) summed off-diagonal fourth moment after normalization. For the
diagonal covariance, disconnected edge supports cancel; shared-edge
tree terms or cyclic terms lose at least one free label. This gives
O(1) total diagonal variance. The remaining fixed-walk moment formulas
are the usual even-edge Catalan counts.

The diagonal-normalization error is controlled by Frobenius norms:
`||B Delta||_F^2=((n-1)/n)||Delta||_F^2`. Together with the regularizer
and E tr(R^2)=O(n), this gives O(sqrt(n)) in the trace formula without
an unstated operator-norm integrability hypothesis. Cauchy--Schwarz of
the summed second and fourth moments gives an O(sqrt(n)) total cubic
arcsine error. The exact beta-integral quotient for (x+2)^r is
`2-6/(2r+3)`. Sending the order to infinity before the degree proves
the positive iid-child energy lower bound 2/pi.

The expectation lower bound is therefore
`(2/pi+sqrt(2/pi)-o(1))n^(3/2)`. Changing one directed incidence sign
changes two actual parent edges, so cap sensitivity is at most four.
Bounded differences gives the claimed high-probability normalized floor
`0.5071738708131547...` for every fixed seed in the iid law.

## Exact balance and operator coupling

The majority repair has uniform balanced-row marginal by coordinate
permutation symmetry. Different rows remain independent. Its exact
mean is `(1-p_i)eta_ij-a_i`, and its centered row has squared norm
`4r_i(1-p_i)` and covariance at most 8p_i I. The matrix mean identities
are valid on the hollow matrices, without hidden diagonal corrections.
The actual low-cap seed bound supplies ||S||op=O_K(n^(3/4)), enough to
make the S diag(a) mean term subleading.

Exposing repaired rows yields exactly the predictable star increments
recorded in Section 6. Their quadratic variation is bounded by
`(16r_max+8sum_i p_i)I`; their norm is at most 2sqrt(r_i). This covers
all interactions with previously exposed and future rows. I directly
read [Tropp, Theorem 1.2 and Corollary 1.3](https://tropp.caltech.edu/papers/Tro11-Freedmans-Inequality.pdf):
the self-adjoint and rectangular martingale hypotheses apply with these
conditional deterministic parameters. They give the stated
O(n^(1/4)log^(5/4)n) operator perturbation, hence an o(n^(3/2)) cap change.

The repaired rows define genuine orthogonal two-row physical frames
`[1;eta_i]`; no orthogonality between different eta rows is claimed.
At Hadamard orders this is exactly the ordinary-column-permutation law
for the first two rows of a dephased Hadamard. Thus the constant-multiplier
actual-parent obstruction transfers to the intended orthogonal frame law.

The theorem is uniform for each seed chosen before frame randomness.
It does not rule out specially correlated frames, favorable rare samples,
or the original globally rewritten-parent construction.
