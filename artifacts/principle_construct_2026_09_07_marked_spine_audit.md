# Independent audit: marked-column recursion and biased slices

2026-09-07. **PASS.** Complete read of
`principle_invent_2026_09_07_marked_column_bellman.md`. I additionally read
the full precision-Schur and stopped-tree dependencies during this audit;
the original recursive ensemble and Fock proof had previously been read
in full for the fixed-seed compiler audit.

## Actual ensemble, rather than just a scalar recursion

Write the signed input permutation as g with underlying permutation Pi.
If the pre-g marked row signs are r=(r1,r1), then direct coordinate
evaluation gives `g D_(row signs after g)=D_r Pi`. Input phases therefore
cancel exactly under physical dephasing. The repeated diagonal D_r
commutes through the normalized H2 mixing. This yields the displayed
marked child Ubar1 and sibling U2 D_r1.

The ordinary law is invariant under right multiplication by every fixed
diagonal sign matrix, because it has its own fresh uniform signed input
permutation. Conditional on the entire first child, U2 D_r1 consequently
has the same ordinary distribution. This conditional law is constant,
so it is genuinely independent of the dephased first child. Pi is also
independent. At a terminal, physical dephasing of fixed Vg leaves fixed
Vbar followed by an unsigned input permutation. No output phases are
silently restored at the marked coordinate.

The mixed terminal Fock estimate has unsigned input covariance and signed
output invariant projection. The unsigned orbit's positive Gram matrix
has constant row sum L_P(v)^2, its largest eigenvalue. Thus its covariance
operator norm is exactly that value. The signed-output invariant rank has
the same partition bound as in the original Fock proof. The Poisson tail
cutoff, uniform on each fixed moment ball, gives the stated subexponential
factor. This is the correct input/output orientation.

## Type probability and the mixed supersolution

For an unsigned permutation, ordered pair-table probability is

    (s/2)! product_a n_a! / [s! product_(a,b) n_ab!].

Its logarithm divided by s is exactly the entropy expression
`-D(pi||nu tensor nu)/2+O_L(log s/s)` under the average FULL-marginal
constraint. There is no unsigned source sign entropy to subtract or add.
Swap-averaging preserves the raw plus law and the symmetrized minus law,
while decreasing the charge and making both marginals equal to nu. For
centered nu, the plus and minus children remain centered and their
average second moment is the parent's moment.

The precision-Schur inequality was indeed proved for arbitrary input laws,
not just symmetric ones. It therefore applies to this raw plus child and
symmetric minus child. It bounds their envelope average minus I(A;B)/2
by the original E_t(nu), with the Gaussian reward evaluated after averaging
conditional variance. No alternative conditional reward was introduced.

## Uniform marked-spine extinction

The ordinary stopped-tree estimate, now with its temperature-alignment
hypothesis discharged by precision-Schur, gives uniform convergence on
every fixed symmetric moment ball. At marked depth ell the moment is at
most 2^ell v, and the ordinary sibling has tree weight 2^-ell. Expanding
only the first L marked nodes and telescoping the envelope bounds gives
exactly the sum of ordinary errors in (8).

Every marked Bellman value is nonpositive. Also E(mu)>=g(Var(mu)), so its
possible positive excess over E is at most K(m_2(mu)). The unresolved
marked descendant therefore costs at most `2^-L K(2^L v)`. This vanishes
uniformly for bounded v. One first fixes L, sends ordinary subtree depth
to infinity, then lets L grow. No asymmetric rigidity theorem, adaptive
CLT, or uniform growing-depth finite-type enumeration is used.

## Actual physical slice and uniformity clarification

For q selected rows and mean a, the physical centered source is
`1_T(x-a)/sqrt(p)`, p=q/m. Its raw three-atom law is precisely (10), has
mean zero, and has variance 1-a^2. Its projection on the marked all-ones
row is exactly zero. The slice contains binom(q,q(1+a)/2) physical words,
giving the entropy p h((1+a)/2), not the symmetric ternary entropy.

Uniformity over varying biases can be made explicit. At every fixed
recursion depth the number of possible numeric letters is bounded by a
constant depending only on that depth and the initial three letters.
Coincident letters may be merged. Stirling and table-count upper errors
are O_r(log m/m), uniformly in all masses, including zeros. Maximizing
over integer tables can only decrease the continuous Bellman supremum;
no lower rational approximation is required for this upper estimate.
Positions remain in a common bounded interval when p stays away from
zero. The supporting-line representation of E and the uniform W2
continuity of the rate-distortion Lagrangians give a common continuity
modulus on these bounded supports. Thus p_n->p and arbitrary varying
physical biases introduce no missing leading error. The common-depth
choice follows from the moment-ball argument, since all variances are
at most one. Summing q+1 bias types has subexponential cost.

After removing the exactly zero marked coordinate, the orbital norm grows
by at most sqrt(m): in the full signed-permutation orbital square, the
event that its distinguished zero is fixed has probability 1/m and
contains exactly the smaller orbital expectation. The canonical
sqrt(2m) allowance is therefore safe. Independent signed permutations of
the other output columns preserve the marked DC column and are legitimate
in the actual weave law.

The result is a valid biased one-row exponent, not yet a mixed-profile
cap bound or a coefficient proportional to 1-a^2. The global diagonal,
heterogeneous row profiles, common temperature, and seed contribution
still require the further calculation explicitly left open in the source.
