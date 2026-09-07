# High-entropy exceptional clouds in actual optimizer energy windows

2026-09-07. This gives many pairs outside the vanishing-W2 noisy-center
sector. It gives NO lower bound on their bridge energies, and no
obstruction to selecting a favorable dense bridge or changing internal
edges in a broader flatification operation.

Fix delta in (0,1/2), alpha=1-2delta, and epsilon in (0,1). Take Walsh
order k=2^d, m=k/2 fibres, and n=mk. Let A be ANY actual optimal child
signing at order n. Choose a cap-attaining center x0, and reverse the
global edge polarity if needed so H_A(x0)=M_n>0. This preserves exact
optimality. In each fibre use the fixed Hadamard frame diag(x0_i) H_k,
before the usual random column signs/permutations. Its center feature
word is sqrt(k) at one coordinate and zero elsewhere.

## Paired noise and its exact profile defect

In b=floor(epsilon*m) fibres, partition physical coordinates into pairs
(u,u+a) for one nonzero Walsh translation a. Give each pair one common
independent noise sign eta, with P(eta=-1)=delta. All remaining physical
coordinates get their own independent noise signs. Set x=x0*eta.
There are exactly g=n-bk/2 independent noise groups, each of size one
or two, and g/n->1-epsilon/2.

In every paired fibre, H_k eta vanishes at EVERY frequency j with
j dot a=1. The joint feature word therefore has an exact atom of mass
at least 1/2 at (0,0). This remains true after common signed column
permutations because the quotient histogram is unchanged.

The reference law retaining the actual center is
(Z,alpha Z+sG), s=sqrt(1-alpha^2). Its mass inside the quotient ball of
radius s/4 about zero is at most P(|G|<=1/4)<1/4: the center's nonzero
coefficient has magnitude sqrt(k), so cannot belong to that ball.
Here the elementary Gaussian density bound phi<1/2 suffices. At least
1/4 of the actual zero atom must move distance at least s/4 in any
coupling. Thus each paired fibre has squared quotient-W2 error at least

    s^2/64.

The averaged squared error is at least (b/m)s^2/64, converging to the
positive constant epsilon*s^2/64. The same conclusion holds at the actual
overlap parameter in a shrinking Hamming window, since it converges to
alpha. Hence this entire family lies outside any vanishing-error sector.

## Actual optimal-child energy remains in the usual noisy window

Ground-switch A by x0 and call the resulting matrix B. Its vertex fields
h_i=sum_j B_ij are nonnegative: flipping one coordinate cannot improve
the global maximum H_B(1)=M_n. Also sum h_i=2M_n. The expectation under
paired noise is exactly

    E H_A(x)=alpha^2 M_n+(1-alpha^2) sum_{paired edges {i,j}} B_ij,

which differs from alpha^2 M_n by at most n/2.

For a variance bound, compress each size-one/two noise group into one
variable. The off-diagonal compressed coefficient C_gh is the sum of
at most four signs, so sum_{g<h} C_gh^2=O(n^2). Its row sum is
sum_{i in g} h_i minus twice its possible internal edge. Consequently

    max_g |sum_h C_gh| <= 2n+2,
    sum_g |sum_h C_gh| <= 2M_n+n,
    sum_g (sum_h C_gh)^2 = O(n^(5/2)),

using the known M_n=O(n^(3/2)). Expand independent eta_g as alpha plus
a centered bounded variable. Distinct linear and quadratic monomials
are orthogonal in L2; their variances are uniformly bounded. Therefore
Var H_A(x)=O(n^(5/2)). Chebyshev, for example with window radius n^(11/8),
shows H_A(x)=alpha^2 M_n+o(n^(3/2)) with probability 1-o(1).
No generic surrogate child or unproved optimizer regularity is used.

The weighted number of physical flips has variance O(n), so its Hamming
fraction is delta+o(1) simultaneously. The same proof applies independently
to the second actual optimal child and its chosen ground center.

## Entropy, including a common exact nearby sphere

The number of negative NOISE GROUPS concentrates at delta*g. On this
typical event every group's product-law atom has probability
exp(-g h(delta)+o(n)). Intersecting with the energy and physical-Hamming
windows still has probability 1-o(1). Hence each child's good exceptional
family has at least exp[(1-epsilon/2)h(delta)n-o(n)] distinct physical
spins. Two independent children give at least

    exp[(2-epsilon)h(delta)n-o(n)]

pairs, all in the actual optimizer noisy-energy windows and all outside
the vanishing-W2 sector on both sides.

One can also select a COMMON exact Hamming distance j_n with j_n/n->delta.
The physical flip-count distribution is identical on both sides because
the group-size multiset is identical. Conditional good fractions on
each side have weighted average 1-o(1). Thus there is a common distance
in the Hamming window with both good fractions at least 1/2 and marginal
probability at least 1/[2(n+1)]. Restricting each family to this distance
loses only a polynomial factor in probability and preserves the same
exponential count. This does not assert that the prescribed distance
floor(delta*n) itself works without another local argument.

At delta=1/10 and epsilon=1/20, the pair entropy is
(39/20)h(1/10)>.633 n, while the squared profile defect on each side is
asymptotically at least .00028125. The count is therefore quantitatively
close to the full-cloud entropy 2h(1/10)n despite a fixed nonzero defect.

## Exact scope

The frames here are legitimate center-adapted Walsh frames in the actual
rank-two ensemble. This shows that arbitrary-center typicality cannot
discard its exceptions at a subexponential counting cost, even inside
actual optimal-child energy windows. The large exceptional family has
NOT been shown to carry large bridge energies. Random bridge gauges may
control it by another argument. Nor does this restrict the user's broader
flatification target, which permits changing the internal edges too.
