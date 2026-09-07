# Independent audit: scalable exceptional cloud and actual optimal children

2026-09-07. The independent track's
`computations/flatify_independent_2026_09_07_exceptional_cloud_exact.py`
was read and independently replayed: PASS. All base signs, Hadamard
frames, exact bilinear values, Hamming distances, and the q=2,4 architecture
inflations were checked by its integer operations. This audit also
reconstructs the all-scale argument independently of those finite tests.

Let R_q[(s,t),(u,v)]=H_q[t,u]H_q[v,s] and w[(s,t)]=H_q[t,s], with q a
power of two. Direct orthogonality gives R_q^2=q^2 I and R_q w=q w.
The tensor bridge C tensor R_q is again the actual rank-two architecture:
the new fibres are indexed by (old fibre,s), and each old frame is
tensored with H_q, with the two columns for an opposite fibre selected
at indices (old pair bit,opposite new fibre). Thus this is not merely
an arbitrary orthogonal matrix inflation.

The base center value zero remains zero, while the base target 152 becomes
152q^3 and each distance becomes 3q^2. At the specified extra coordinates,
the old signed fields are 0 and 2 and their signed mutual coupling is 1.
Flipping an identical h-subset S of the corresponding inflated coordinates
on both sides changes the target by

    -4qh + 4 w_S^T R_q[S,S] w_S.

A nonnegative subset quadratic always exists, not just at the tested
q=2,4. Indeed K=diag(w)R_qdiag(w) has diagonal one and constant row sum q.
For a uniform subset of size h its expected quadratic is

    h + h(h-1)*(q^3-q^2)/(q^2*(q^2-1)) >= 0.

Take h=floor(q^2/5). The exact distances are then floor(32q^2/10), and
the target is at least (756/5)q^3. Dividing by (32q^2)^(3/2) gives
189/(160sqrt(2))>.83. The centers still have bilinear value exactly zero.

## Actual-minimizer lifting

This bad bridge alignment can coexist with TWO actual optimal children.
Let A be any optimal signing at order n, switch a cap-attaining spin to
the all-one vector, and if necessary reverse the global edge polarity
so that its energy there is +M_n. Randomly permute its vertices and then
switch the all-one vector to the prescribed bridge center x0. For a
prescribed target x at distance j from x0, every edge's expected target
spin product is

    kappa_n=((n-2j)^2-n)/(n(n-1)).

Therefore the expected target child energy is kappa_n M_n, and some
permutation achieves at least that. Make this choice independently on
the other side. Switching, permutation, and global polarity preserve
actual optimality exactly. Both center child energies equal M_n, while
both target child energies are at least kappa_n M_n. No near-minimizer
surrogate or prescribed nonoptimal child has entered the argument.

For j=floor(n/10), kappa_n=16/25+O(1/n). Combining the positive bridge
target with those child energies gives a full parent's actual quadratic
cap at least

    [189/(160sqrt(2)) + (32/25)(M_n/n^(3/2)) - o(1)] n^(3/2).

The active universal upper bound M_n/n^(3/2)<.493608094+o(1) makes this
larger than 2sqrt(2)M_n by more than .07 n^(3/2) asymptotically. Thus
the missing exceptional sector can genuinely break a whole-parent
inequality even with two actual optimal children and zero center bridge.

This establishes EXISTENCE of bad architecture/child alignments. It does
not show every favorable choice fails, does not refute the proved random
profile-sector existence theorem, and does not prove an exponentially
large family of exceptional pairs. The all-scale sequence is sufficient
for its negative scope; no all-order inflation is asserted here.
