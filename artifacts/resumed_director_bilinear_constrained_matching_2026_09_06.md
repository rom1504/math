# A bilinear compatibility constraint does not prevent asymptotic cloud matching

Date: 2026-09-06. Status: director proof independently reconstructed by
the bound-audit researcher, including the augmenting-path and tail arguments.
This is a deterministic matching theorem and a scope check for the new
Fourier-permutation realization target. It does NOT realize such a
permutation as a vectorial dual-bent index map.

## 1. Deterministic theorem

Let N=2^r tend to infinity. Index two point clouds x_u,y_v in a fixed
Euclidean space by u,v in F_2^r. Suppose both empirical probability
measures converge in W2 to the same measure with finite second moment.
Then there are permutations pi of F_2^r satisfying

    u dot pi(u)=0 for EVERY u,
    (1/N) sum_u ||x_u-y_pi(u)||^2 ->0.                 (1)

The point-cloud values can depend arbitrarily on their labels. No
independence, random relabeling, or uniform bound on their norms is needed.

## 2. The compatibility graph has no large empty rectangle

Let E be the bipartite graph with edge u~v exactly when u dot v=0.
Every degree is at least N/2; the zero vertices have degree N. Write
H_uv=(-1)^(u dot v), so HH^T=NI and E=(J+H)/2. For subsets S,T,

    |e(S,T)-|S||T|/2| <= sqrt(N|S||T|)/2.             (2)

Consequently an empty rectangle obeys |S||T|<=N. This is an exact
spectral calculation for the ordinary Walsh matrix.

Inside any pair of vertex subsets S,T, a maximum matching leaves
unmatched subsets S0,T0 with no edge between them. Thus

    matching_size(S,T)>=min(|S|,|T|)-sqrt(N).          (3)

## 3. Every partial matching can be completed cheaply

For N>4, a matching missing k vertices on each side can be completed to
a perfect matching while deleting at most 2k of its ORIGINAL matched
edges. To prove it, take unmatched vertices u on the left and v on the
right. If either has an unmatched neighbor, augment with one edge.
Otherwise all their neighbors are matched. Let S be the left partners
of the right neighbors of u, and T the right partners of the left
neighbors of v. Both sizes are at least N/2, so (2) implies an edge
s~t between them, because N^2/4>N.

If t is matched to s, there is an alternating three-edge path from u
to v. Otherwise let r be s's matched right partner and l be t's matched
left partner. Then

    u -- r -- s -- t -- l -- v

is a five-edge augmenting path. The internal vertices are distinct in
the second case, since the two matched pairs differ. The augmentation
deletes at most two currently matched edges. Repeat k times. An original
edge deleted more than once can only improve the bound; at most 2k
distinct original edges are lost in total.

## 4. Matching the cloud geometry

Take a finite partition of a large ball into sets of diameter at most
eta whose boundaries have zero limiting measure. Match x and y labels
inside each corresponding partition cell using (3). The empirical cell
masses converge to the same limits. If there are L cells, the total
unmatched proportion is at most

    outside-ball masses + sum(cell-count discrepancies)/N
                                               + L/sqrt(N).  (4)

Choose ball radii tending to infinity and cell diameters tending to zero
slowly enough that this proportion tends to zero, by a diagonal choice.
Now complete the partial matching using Section 3. All but o(N) of the
original within-cell pairs survive. Their squared distance is at most
eta^2. The remaining o(N) pairs have total normalized squared cost o(1):
W2 convergence implies uniform integrability of the empirical squared
norms, so ANY subset of o(N) labels carries o(N) squared norm mass.
Use ||x-y||^2<=2||x||^2+2||y||^2 on the exceptional pairs. This proves (1).

The same proof works for triangular random clouds whose TWO empirical
laws both converge in W2 in probability to the SAME FIXED finite-second-
moment law: pass to deterministic realizations or an almost-sure subsequence.
Merely requiring their mutual W2 distance to tend to zero is insufficient.
For example, identical one-dimensional clouds x_u=y_u=N times the integer
label of u have zero mutual distance, while the constraint forbids the
diagonal on half the labels and forces diverging matching cost. These
clouds are not uniformly integrable. No independent Fourier rows are assumed.

## 5. Consequence for the newly proposed polar matching test

In the Fourier-cloud theorem, the active frequencies can be labeled
alpha=(1,u) and beta=(1,v). The standard bilinear restriction

    alpha dot beta=1

is exactly u dot v=0. Section 1 shows that the Gaussian-cloud polar
matching can respect this constraint with vanishing average defect.
Extend the resulting bijection arbitrarily on the inactive hyperplane;
the profile Fourier coefficients there are zero.

Thus a necessary bilinear identity for quadratic vectorial dual-bent
index maps does not, by itself, rule out profile-specific asymptotic
matching. A theorem forbidding a permutation from fixing the INACTIVE
hyperplane pointwise must not be upgraded to a no-go for the needed
profile action: that pointwise fixing is unnecessary.

This statement treats the standard compatible bilinear form only. An
arbitrary transformed form and prescribed pair of hyperplanes needs its
own degree and mixing check. More importantly, being a perfect matching
in this graph is far weaker than arising from inverse-linear spaces of
nonsingular polar matrices. The latter realization remains unproved.
