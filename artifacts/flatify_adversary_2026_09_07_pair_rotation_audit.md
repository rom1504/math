# Independent check of pair-local rotation rigidity

Status: PASS. Independently reconstructed the director's
`flatify_director_pair_rotation_rigidity_2026_09_07.md`.

For `sigma_ij=a_ij d_ij=+1`, conjugating the tile gives a rotation through
theta_i-theta_j; for sigma=-1, it gives a reflection with angle determined
by theta_i+theta_j. Every tile entry has magnitude one precisely when,
with phi_i=4theta_i/pi,

```math
\phi_j=\sigma_{ij}\phi_i+1\pmod2.
```

Composing around a triangle multiplies phi_i by the three-sigma product and
adds a sum of three signed ones. That sum is an odd integer. A positive
product is impossible modulo 2. Negative product on every triangle is
equivalent to sigma_ij=-t_i t_j: set t_1=1,t_i=-sigma_1i and use triangles
through vertex 1. Choosing phi_i=t_i/2 checks both equal-sign and opposite-sign
cases of the tile relation, proving sufficiency. This verifies the full
if-and-only-if, not just a necessary cycle condition.

The exact row factor is sqrt(2+1/(k-1)), differing from sqrt(2) by O(1/k).
Multiplying a two-child cap O(k^(3/2)) gives O(sqrt(k)) normalization error.
The k empty paired entries cost at most k when filled. For selectable D=-A,
the exact tiles reduce, up to switching/reordering, to the ordinary H2 lift.

The separate unbounded-amplitude example also checks: each hollow bipartite
H2 block has unscaled cap 2 and row squared norm 2. Scaling by
sqrt((N-1)/2) and taking N/4 disjoint blocks gives total cap
N sqrt(2(N-1))/4 and row squared norm N-1. Its limiting constant 1/(2sqrt2)
is below the full-sign lower endpoint, but its entries grow as sqrt(N).
It does not refute fixed-bounded-amplitude selected-child replacement.

No normalization, sign, or scope correction was found.
