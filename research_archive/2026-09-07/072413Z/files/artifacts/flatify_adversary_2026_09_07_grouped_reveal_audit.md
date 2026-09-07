# Independent audit: two-orientation cap-preserving reveal

Date: 2026-09-07. The independent agent proposed the following genuinely
cap-preserving modification of the hidden-partition reveal. Algebra and
count-comparison scope check; signed endpoint drift remains open.

Let P be the set of revealed positive vertices, let H_P be the internal
P energy, and let H_R=H-H_P. Use the fixed ambient denominator and define

    G(v,P)=min_A log sum_x cosh(beta H_P/sqrt(N-1))
                            cosh(beta H_R/sqrt(N-1)).

This is the partition function with two independent orientation bits,
one multiplying internal P edges and one multiplying every other edge.
Initially H_P=0, giving original optimized absolute pressure exactly.
At the terminal split, H_P and H_R are the two disjoint child energies.
The spin sum factors, so the terminal optimum is EXACTLY the sum of the
two independently optimized child absolute pressures. In contrast to the
shared-orientation terminal, it does not become twice paired-width
pressure for equal children.

For any fixed signing, the corresponding zero-temperature objective is

    max_x [|H_P(x)|+|H_R(x)|]
      =max{Q(A),Q(A with all internal P signs reversed)},

where these Q's use the current weighted profile. Thus this is an actual
two-polarity robust cap, not merely a renamed single cap. Its optimized
intermediate value may exceed the ordinary optimized absolute cap.

## Count cleanup survives, but the changed assignment must be tracked

The edge-flip derivative proof works under the joint law of the two bits
and the spins: each edge still has a single +/-1 observable. Vertex-star
contraction works after conditioning on both orientation bits, summing the
vertex spin, and using the resulting cosh(row field). When comparing two
counts, delete the known vertices whose labels differ. On the remainder,
both the edge weights and the assignment of edges to orientation groups
have the same structural form. The diagonal-scaling proof then applies
unchanged. Common-terminal deletion also works because the terminal
grouped optimum depends only on the two child sizes.

Direct cap rounding likewise applies: include both relative orientation
choices in the Bernstein union bound. This changes its logarithmic
parameter by only a fixed multiple of log2. Thus the O_delta(N^(5/4))
one-time expected count cleanup remains valid for the robust cap too.

The reveal's frozen-branch drift must incorporate changes in orientation
assignment, not only changes in variance. A continuous-variance generator
for a fixed edge-group assignment is insufficient by itself to describe
the full discrete reveal.

## Exhaustive actual-sign diagnostic

`computations/flatify_adversary_2026_09_07_grouped_reveal_drift.py` exhausts
all full signings modulo switching at N=4,6,7 and beta=.7,1.5,3. The output
is saved in computations/results with the same stem. It asserts the
terminal factorization against separate full-sign child optimizations,
keeps vertex labels aligned, and checks the drift telescoping identity.
These transcendental calculations use floating point, not certified intervals.

At beta=3:

| N | terminal minus initial | smallest conditional drift | largest |
|---|---:|---:|---:|
|4| -0.2927264031 | -1.0053644288 | 1.0689570386 |
|6| 3.5125163313 | 0 | 1.8426499388 |
|7| 1.6431163246 | -0.2288786012 | 1.0541672360 |

All conditional drifts at N=6 are nonnegative at the three tested beta
values, but both signs survive at N=4 and N=7. Frozen drift and selection
loss are also saved; their separate values depend on which tied current
optimizer the enumerator selects. Only their difference is the optimized
drift. No scalable obstruction or favorable asymptotic drift follows from
this finite test.
