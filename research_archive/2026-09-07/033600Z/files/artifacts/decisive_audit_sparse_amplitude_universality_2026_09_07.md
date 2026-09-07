# Independent audit: sparse amplitudes defeat two-sided optimizer universality

The proof in
`decisive_independent_sparse_amplitude_universality_obstruction_2026_09_07.md`
passes fresh reconstruction.  In particular its event is simultaneous
over all signings; no exponentially large union bound is missing.
The conclusion concerns a variance-one **triangular array**, not fixed
Gaussian magnitudes, and does not refute favorable one-sided flatification.

## 1. Exact local wave and the scale of its error

For a regular tree, the sphere sizes are `s_k=d(d-1)^(k-1)`.  The
normalized path eigenvector coefficients in the source satisfy
`sum a_k^2=1`, `a_0=0`, and `sum a_k a_(k+1)=cos(pi/(r+1))`.
For a central edge, a vertex at distances `k,k+1` from its endpoints
lies in one of exactly two side branches, each of size `(d-1)^k`.
The two signed paths differ precisely by the central edge.  Thus every
covariance summand, whatever the signing, has its sign equal to that
edge sign, and the total correlation is exactly

    A_ij 2sqrt(d-1)/d cos(pi/(r+1)).

The radius condition is sufficient: both radius-r neighborhoods of an
edge lie in the radius-(r+1) ball of either endpoint.  A cycle in such
a ball would create a breadth-first-search non-tree edge and a cycle
of length at most `2r+3`.  Hence the stated high-girth assumption
supplies the joint tree computation.

For nearly regular degrees in `[(1-epsilon)d,(1+epsilon)d]`, the
number of vertices in any fixed-depth side branch is
`d^k[1+O_r(epsilon+1/d)]`.  Sphere counts have the same relative
error.  With coefficients `a_k d^(-k/2)`, variance is
`1+O_r(epsilon+1/d)`, while edge covariance is

    A_ij (2/sqrt(d))
      [cos(pi/(r+1))+O_r(epsilon+1/d)].

Dividing by the two standard deviations retains that relative error.
It is not merely an absolute `o(1)` error, which would be insufficient
after summing edges.  Gaussian sign rounding, followed by
`arcsin(t)>=t` for `t>=0`, gives the precise signed-energy factor
`(2/pi)n sqrt(d)` after multiplication by `nd/2` edges.

## 2. The random graph event really has probability tending to one

Take `d=(log n)^3` and `epsilon=(log n)^(-1/2)`.  Chernoff plus the
union bound gives all degree bounds since `epsilon^2 d=(log n)^2`.
Its centering at `(n-1)p` rather than `d=np` changes only `O(1/n)`
relative error, smaller than the stated window.

For fixed `r`, a vertex whose radius-(r+1) ball contains a cycle
has a connected rooted unicyclic witness: two BFS paths and one
non-tree edge.  This witness has `e=v<=2r+3`, **including its root**.
The number of graph embeddings times their edge probabilities is
at most `n^v(d/n)^e=d^e`.  There are finitely many such rooted shapes
at fixed r.  Hence the expected **total** number of bad roots is
`O_r(sum_(e<=2r+3)d^e)`, with no extra factor n.  Markov then makes
the number at most `n^(1/4)` with probability tending to one.  The
degree bound implies only `o(nd)` edges touch those vertices.

At bad vertices the proof uses distinct *fresh* independent Gaussian
coordinates.  A good vertex may still use the old Gaussian associated
with a bad graph vertex; this does not impair independence from that
bad vertex's fresh output coordinate.  Every edge with a bad endpoint
therefore has exactly zero expected signed energy.  No adverse
contribution from those edges is omitted.

All of these events and quantitative counts depend only on the
unsigned support graph.  On their intersection the Gaussian-wave
proof works for **every** edge signing, with the same constants.
The auxiliary Gaussian spin witness may depend on the chosen signing,
which is allowed when lower-bounding its maximum.  Consequently the
bound survives minimization over all signings without any probabilistic
union over them.

## 3. Normalization and the noninterchangeable quantifiers

For `g_e=epsilon_e xi_e/sqrt(p)`, adaptive edge signs absorb the
independent `epsilon_e`.  The support signing cap is multiplied by
`p^(-1/2)=sqrt(n/d)`.  Thus

    (2/pi-o(1))n sqrt(d) cos(pi/(r+1))

becomes `(2/pi-o(1))cos(pi/(r+1))` after division by `n^(3/2)`.
For each desired error first choose a **fixed** r and then let n
grow; only afterward improve r.  The polylogarithmic cycle count is
not claimed for a radius growing with n.

The entries are centered and variance one, and deterministically
`max |g_e|/sqrt(n)<=1/sqrt(d)->0`.  Also
`E|g_e|^3/sqrt(n)=1/sqrt(d)->0`.  For fixed deterministic signing,
a normalized finite-temperature spin pressure has coordinate third
derivative `O_beta(n^(-5/2))`; summing the usual replacement errors
over `O(n^2)` entries is therefore `o(1)`.  The auxiliary global
orientation spin can be included in that smooth fixed-sign pressure
without changing this derivative bound.  This calculation does not
justify replacement after a nonsmooth adaptive minimum over signings.

Hence the example separates optimized sparse triangular-array value
at least `2/pi-o_P(1)` from the deterministic fair-sign-amplitude
value with limsup below `.494515125`.  It proves failure of a
variance-only, two-sided optimized universality principle even under
that ordinary third-moment scale.  It proves nothing comparable for
fixed standard-Gaussian magnitudes; they are not locally tree-like.
It is compatible with `M_n/n^(3/2)<=T_n(g)+o(1)`.  The positive-branch
lower bound for the sparse graph is stronger than needed, but does
not change the fact that the original dense problem is two-sided:
all-negative dense couplings have one-sided cap only `O(n)`.

## 4. Reproducible local checks

`computations/decisive_audit_sparse_gaussian_wave_2026_09_07.py` checks
the signed local variance and covariance identities with exact
rational arbitrary radial coefficients on edge-rooted trees, then
checks the normalized sine wave and independent bad-coordinate repair.
The asymptotic graph-event and all-signings quantifiers are proved
above, not inferred from this finite replay.
