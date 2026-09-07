# Independent affine-MUB partial-packing reconstruction

Status: **PASS**, with the fixed-seed quantifier retained. Independently checked
`flatify_construct_2026_09_07_affine_mub_partial_packing.md`, including its
external coloring hypotheses and arbitrary-completion cancellation. No
quantitative failure was found in the stated theorem.

## 1. Exact transform and spike upper bound

Each line matrix consists of Hadamard rows on disjoint line supports divided
by sqrt(q); hence it is orthogonal. Distinct-class rows overlap in one physical
point, so their inner product is exactly +-1/q. Thus `q a_ij U_i U_j^T`
really is a full-sign cross block. For y_i=U_i^T x_i, the hollow cross energy is
`q sum_p H_A(y(p))`, with no doubled-edge factor.

A signed Hadamard column on each physical line transforms into exactly one
spike of amplitude sqrt(q). Consequently y_i(p)=sqrt(q)z_i(p), each fibre uses
q spikes, and the energy is `q² sum_p H_A(z(p))`. Dividing pattern frequencies
by q gives a feasible solution of the stated packing LP, proving the exact
upper `sigma H_C<=q³ kappa_sigma(A)` on the pure-column sector.

This upper does not include arbitrary diagonal-fibre completions or non-spike
states; the source statement correctly does not claim either extension.

## 2. Random-support sampling and coloring

Fix k and an optimal packing solution. One maximizing signed pattern per
support is sufficient, since all patterns on that support have the same
constraint column. Supports of size zero or one have zero objective and can
be discarded. Choose fixed L>=max lambda_S and retain each e(p,S) independently
with probability lambda_S/L.

There are kq line vertices and O_k(q²) potential edges. At a fixed line,
the expected degree is `(q/L) sum_{S containing i}lambda_S<=q/L`. Bounded-sum
Chernoff estimates and a union bound over kq lines give maximum degree at
most q/L+o(q) with probability tending to one. At least one constraint is tight
at a nonzero optimum (otherwise scaling all weights improves it), so the
maximum degree is also at least q/L-o(q).

Any two vertices from different classes determine at most one point. At that
point at most 2^(k-2) supports contain both classes; vertices in the same class
have codegree zero. Distinct (p,S) cannot give the same edge when |S|>=2:
the classes determine S, and two of the lines determine p. Thus the hypergraph
has bounded rank, D tending to infinity, and maximum pair degree O_k(1)=o(D).

The primary Rutgers record of Kahn's
[Asymptotically good list-colorings](https://www.researchwithrutgers.org/en/publications/asymptotically-good-list-colorings/)
explicitly states the bounded-rank, maximum-degree, small-codegree theorem
needed here, without a minimum-degree requirement. Ordinary coloring is a
weaker consequence. Thus `(1+o(1))q/L` matchings suffice; no hidden regularity
assumption is imported into the randomly thinned general hypergraph.

For the pentagon-only uniform construction, the primary institutional record
of [Pippenger--Spencer](https://scholarship.claremont.edu/hmc_fac_pub/1041/)
requires fixed uniform rank, asymptotically equal minimum and maximum degrees,
and negligible codegree. The 3-uniform, exactly 3q-regular pentagon hypergraph
satisfies those stronger hypotheses directly.

## 3. Arbitrary internal completion: the cancellation really is exact

For a matching, use an independent fair sign on each matched component and
on each unmatched physical line. A matched component contains at most one
line from each fibre. Therefore every cross term between different components
averages to zero, including all within-fibre interactions between distinct
physical lines. Possible coincidences of selected points between components
do not change this: their independent component signs still cancel.

On a matched physical line, the remaining internal mean is the deterministic
Hadamard-column energy b_{i,L}(p), unchanged by the component's overall sign.
On an unmatched line the uniform Hadamard-column average is zero, because
`sum_p h_p h_p^T=qI` and the line restriction of the arbitrary D_i is hollow.
This proves the matching's expected signed energy formula (8) exactly.

Summing over retained hyperedges, the expected internal contribution is

```math
\sum_{i,L}{\sum_{S\ni i}\lambda_S\over L}
                 \sum_{p\in L}b_{i,L}(p)=0.
```

The cross expectation is `(q⁴/L)kappa_sigma`. The total is a sum of O_k(q²)
independent Bernoulli inclusions, each multiplied by a fixed O_k(q²) number,
so its variance is O_k(q^6), uniformly over the line Hadamards and D_i.
For example, a q^(7/2) deviation threshold is o(q⁴) and has failure probability
O_k(1/q). Intersect this event with the degree event. Dividing the resulting
sum by the number of colors gives a matching whose expected energy is
`[kappa_sigma-o(1)]q³`; some Boolean state in its law attains that expectation.
Every state in that law is still in the pure-column spike sector.

This checks the lower bound against all arbitrary completions, and the
matching upper/lower asymptotics for the zero-completion spike sector.

## 4. Pentagon exactness and normalization

For the cycle/chord pentagon, each polarity has five favorable triangles, and
every vertex belongs to three. The weight 1/3 packing has value 5. The dual
diagonal d_i=1 is feasible: supports of sizes 0,1,2,3,4,5 have absolute caps
at most 0,0,1,3,4,4. Hence both kappas equal 5.

The pentagon has cap 4. Its global minimality is self-contained: any order-five
full signing has second moment binom(5,2)=10, and every energy is even, so
its cap cannot be below 4. No external lower solver is needed.

The uniform pentagon hypergraph has 5q² edges and degree 3q. The total signed
cross weight is 15q⁴; internal terms cancel in the complete edge sum. Coloring
therefore supplies `(5-o(1))q³`. At N=5q² this is `1/sqrt(5)` after dividing
by N^(3/2). This exceeds both the naive coefficient-one seed constant
`4/5^(3/2)` and its row-normalized version `4/(5sqrt(4))=2/5`.

## 5. Independent integer replay

`tmp/flatify_adversary_2026_09_07_affine_packing_replay.py` imports no constructive
code. It enumerates all 243 ternary pentagon states, constructs the F4 net,
and verifies the line/intersection/basis identities. For 20 independently
randomized sets of line Hadamards and full arbitrary diagonal-fibre signings,
it checks both polarities, producing 40 cases at full parent order 80.

For each case it checks the complete hypergraph counts, degree and codegree,
and exact cancellation of the total local-line energy. It constructs a greedy
matching only as a finite diagnostic, builds the exact covariance of its
actual Boolean law, and verifies `sigma Tr(F Cov(x))/2` equals the matching
formula including the arbitrary D_i. All arithmetic after the finite random
choices is integer. Every case passed; parameters and outcomes are in the
adjacent JSON.

An initial Python-3.9 incompatibility (`int.bit_count`) stopped before the
first run; it was replaced by `bin(...).count('1')`, with this failure recorded
in the source. It had no mathematical consequence.

## 6. Scope and quantitative limitation

The coloring errors depend on fixed k and the fixed LP solution. There is no
uniform rate as k grows, and no such rate is asserted. The pentagon supplies
a genuine scalable obstruction to the stated affine-MUB construction class
from an actual exact optimal seed. It does not supply a bad sequence of
selected large minimizing children, nor an obstruction to non-affine MUBs,
unrestricted global flatification, convergence, or the new all-order upper.
