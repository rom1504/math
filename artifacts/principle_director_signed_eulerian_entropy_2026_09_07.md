# Signed Eulerian incidence entropy and a strong-coupling type compiler

2026-09-07. Complete proof attempt; independent reconstruction requested.
This is a positive counting/realization theorem for signed incidences.
Its application also rules out one proposed seed-retention mechanism.
It is NOT a theorem identifying the optimum of the original signing problem.

## 1. General signed counting theorem

Let G be a connected loopless graph with positive even degrees d_v=2k_v.
Give each edge e=vw a prescribed sign s_e. A balanced incidence assignment
assigns signs u_(v,e),u_(w,e) such that their product is s_e and each vertex
has exactly k_v positive and k_v negative incident signs. Write C_G(s) for
the number of these assignments and e=|E(G)|.

There are assignments if and only if

    product_edges s_e = (-1)^e.                              (1)

When (1) holds,

    C_G(s) >= product_v [binom(2k_v,k_v)/(2^k_v k_v)].         (2)

For G=K_m with m odd, k=(m-1)/2, this yields

    log C_G(s) >= e log2 -(3/2)m log k - O(m),                (3)

uniformly in the entire edge signing. The existence statement also follows
directly from one Euler tour; the count is the substantive additional step.

### External theorem used, with checked hypotheses

Schrijver's Eulerian-orientation lower bound states that the number EO(G)
of ordinary orientations with equal in/out degree is at least

    product_v binom(2k_v,k_v)/2^k_v.

The graph may be irregular, and all degrees must be even. The precise
statement is Theorem1.1 of Borbenyi--Csikvari,
[arXiv:1905.06215v3](https://arxiv.org/html/1905.06215), a primary paper
giving another proof of Schrijver's result. The original published paper
is [Schrijver1983](https://doi.org/10.1007/BF02579193), also available via
[CWI](https://ir.cwi.nl/pub/10053). No signed version is imported.

### Transition systems and a last-exit proof

At each vertex pair its incident halfedges. A transition system is the
collection of these pairings. It partitions all edges into closed trails;
the trails may revisit vertices. A single-trail transition system uses all
edges in one such trail. Let T_1 be the number of single-trail systems.

Every ordinary Eulerian orientation of connected G is strongly connected:
otherwise a source component of its condensation has outgoing but no
incoming edges, contrary to conservation of total in/out degrees.
Fix a root vertex r and one outgoing edge at r. Choose a directed spanning
tree toward r. At each nonroot vertex put its tree edge last in an ordering
of outgoing edges; at r put the fixed edge first. There are

    product_v (k_v-1)!

choices of the remaining orders. Follow the next unused outgoing edge,
starting at r. The last-exit tree guarantees that the walk cannot exhaust
the root's edges while unused edges remain: vertices with unused edges
would have a chain of unused last-exit edges leading back to r. Thus it
traverses every edge exactly once and determines a single-trail transition
system. With the orientation and first root edge fixed, its transition
system reconstructs the tour and all local orders; the construction is
injective.

Each single-trail transition system is compatible with exactly two ordinary
Eulerian orientations, the two directions of its trail. Consequently

    2 T_1 >= EO(G) product_v (k_v-1)!.                       (4)

For the prescribed signs, propagate a halfedge sign around any single trail,
reversing at each vertex transition and multiplying by s_e across an edge.
It closes exactly when (1) holds; in that case there are two assignments.
They are balanced at every vertex. Conversely, a balanced assignment is
compatible with precisely product_v k_v! transition systems: pair every
positive halfedge to a negative one independently at each vertex.
Counting only single-trail systems therefore gives

    C_G(s) product_v k_v! >= 2 T_1.

Together with (4) and Schrijver this proves (2). Multiplying all incidence
signs proves necessity of (1). This proof requires trails, NOT a decomposition
into simple balanced cycles; the latter is a stronger and sometimes false
property of signed graphs.

## 2. Exact-type probability at low noise

For K_m, independently at each edge draw the two incidences with probability

    P_s(a,b)=[1+rho s_e a b]/4,     a,b in {+1,-1},

where 0<=rho<1. Let T be the event that every row is balanced. If the seed
parity is compatible, all C_G(s) perfectly matched assignments contribute,
so

    P_s(T) >= C_G(s) [(1+rho)/4]^e.

With rho=tanh(2t), (3) gives the explicit uniform bound

    -log P_s(T) <= (3/2)m log k +O(m)
                              +e log(1+exp(-4t)).            (5)

If parity is incompatible, reverse one prescribed edge to make it compatible
and use the same assignments. Exactly that edge then receives the smaller
weight, with ratio (1-rho)/(1+rho)=exp(-4t). The right side of (5) acquires
only an additional 4t. This is an exact one-edge repair, not an unproved
low-cost many-row adjustment.

For t>= (1/8)log m the right side is O(m^(3/2)+m log m+4t); at the critical
t=Theta(sqrt m) it is already O(m log m), since e exp(-4t) is negligible.
Thus critical-temperature type feasibility has subleading entropy cost for
EVERY signing, without spectral assumptions. At a parity-compatible hard
endpoint rho=1, the bound remains meaningful and proves the corresponding
uniform counting statement directly.

## 3. Application to the Gaussian binary kernel

For symmetric binary row type and K_t(a,b)=exp[-t(a-b)^2], the exact canonical
edge coupling is the preceding P_s with rho=tanh(2t). Its edge variational
value is V_t=log[(1+exp(-4t))/2]. The exact partition identity is

    Z_s=exp(e V_t) P_s(T)/P_0(T),

where P_0 has independent fair incidences. The denominator has logarithm
O(m log m). Hence (5) pays the microscopic row constraints at strong coupling,
including the critical scale which the earlier independent-edit compiler
could not reach with a subleading bound.

For intermediate/bounded t, the earlier independent-edit argument gives
O((1+t)m^(3/2)log m) error, because the log ratio of largest to smallest
edge atoms is exactly4t. Combining the two bounds supplies a subleading
pressure error throughout t=o(m^2), with the parity-compatible class having
no 4t penalty. This rules out retaining a leading seed distinction merely
by increasing t to sqrt(m) in this SYMMETRIC BINARY model.

The theorem does not treat asymmetric prescribed degrees, general growing
alphabets, actual Hadamard row-spin realizability, or the full parent cap.
Those are separate obligations. The positive principle here is that a
global compatibility constraint can be compiled at subextensive entropy
cost by counting single-trail realizations; its signed cycle obstructions
collapse to one parity bit in this particular incidence model.

## Verification record

`computations/principle_director_signed_balance_2026_09_07.py` enumerates all
switching classes at m=3,5 and all edge-sign assignments, checking the exact
parity criterion. It gives positive counts2 at m=3 and between12 and24 at
m=5. These are finite checks, not inputs to the asymptotic count proof.
