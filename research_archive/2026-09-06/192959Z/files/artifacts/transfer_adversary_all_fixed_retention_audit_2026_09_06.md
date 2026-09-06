# Independent audit: all fixed retentions on near-minimizers

Date: 2026-09-06, approximately 18:58 UTC. Status: PASS for
`transfer_seed_all_fixed_retention_nearmin_loss_2026_09_06.md`.
This is an audit of the stated scoped result, not a new selector
extension or a convergence theorem.

## Weak limit and first exit

Fix `q in (0,1)` and a finite `r` with `q^r<rho`. Along exact parent
minimizers realizing `ell=liminf M_n/n^(3/2)`, nested exact-size uniform
sampling has uniform marginals at every depth. Repeated floors satisfy
`d_(j,s)/D_j -> q^s` for each of the finitely many `s`.

Principal monotonicity bounds all normalized cap coordinates by one
fixed constant. The definition of the global liminf bounds them below
by `ell-o(1)`, uniformly over this finite list of deterministic sizes.
Thus a joint weak-limit subsequence exists. The terminal probability
bound implies that its last coordinate is at least `1/2` almost surely;
it does not need to imply strict inequality there. Since `ell<1/2`,
this already forces a first exit from the exact level `ell`.

The events whose previous coordinates all equal `ell` and whose next
coordinate exceeds `ell` partition almost surely. One has probability
at least `1/r`. Continuity from below then permits a positive `epsilon`
with the strict probability margin `3/(4r)` in the source proof.
The open pair event with parent coordinate below `ell+delta` and child
above `ell+4 epsilon` contains that limiting event. Portmanteau and a
deterministic diagonal choice `delta_j -> 0` therefore give the stated
prelimit margin `2/(3r)`. No sampling depth grows with order.

## Deterministic parent extraction and constants

Condition on the entire selected parent, not on a favorable child.
The next selector is still uniform at exact size `floor(q d_j)`.
If every parent below `ell+delta_j` had conditional child probability
less than `1/(2r)`, the joint pair event could have probability at most
`1/(2r)`, contradicting the larger margin. Selecting one deterministic
parent realization per order is consequently legitimate.

These selected parent caps converge to `ell`. Their orders themselves
realize the liminf: the global minimum at such an order is at most the
selected cap and is at least `ell-o(1)` after normalization. Hence the
normalized optimality gap also tends to zero. On the selected child
event, the increase is greater than `3 epsilon` once `delta_j<epsilon`,
so the theorem's weaker increase `epsilon_q=epsilon` and probability
`eta_q=1/(2r)` are valid with room to spare.

The finite first-exit variant is valid as well. Earlier no-exit history
only weights the distribution of good selected parents. Conditional
on the current parent, the next-step law remains uniform by the
Markov property of the nested sampling process. Thus a largest
first-exit event yields the stated deterministic parent and conditional
probability without assuming that a restriction of a minimizer is
itself a minimizer.

## Essential scope

The theorem produces some offending liminf-realizing NEAR-minimizer
sequence for every fixed retention. It does not produce exact optimizing
parents, does not assert failure on every near-minimizer sequence, and
does not rule out a well-chosen seed sequence or exceptional selector.
Those distinctions are correctly stated in the source. The theorem
does not settle the original convergence problem.

## Discrete-entropy note: scope audit

The scope of
`transfer_reconstruction_discrete_threshold_entropy_target_2026_09_06.md`
also passes. Its exact fixed-order entropy curvature is
`eta^2(N-1)/2`, obtained from orthonormal pair products under the
independent threshold law. It explicitly does not claim a uniform
Taylor remainder at fixed nonzero covariance parameter as order grows.

Its independent two-coordinate block example has covariance entries
of order one, not `s/sqrt(N)`. It is therefore a valid obstruction to
an operator-only rare-threshold contraction, but not a counterexample
to the actual flat-sign target; the note states this correctly.

One minor sign convention was reported: when thresholding at positive
`z`, use `X=1-2*1{G>=z}`, or explicitly reverse both bits, to retain
the product mean `m=1-2 delta`. Literal `sign(G-z)` has mean `-m`.
The mutual information and entropy comparison are unchanged by this
simultaneous relabeling. No mathematical conclusion depends on it.

No numerical test is needed for the first-exit theorem. No temporary
research files were produced in this audit. Further selector extensions
are frozen as directed; the remaining original comparison obligation
is still a positive cap-preserving construction or another genuine
cross-order comparison.
