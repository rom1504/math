# Independent audit of the bounded-cap ramp rank/coherence counterexample

Verdict: PASS after reading the complete deterministic construction in
`decisive_independent_ramp_rank_coherence_counterexample_2026_09_07.md`.
This is not a counterexample to sparse sign repair or to a condition on
selected actual width minimizers.

With the author's notation b=4^s, l=16b, n=bl, one has l=4sqrt(n).
Both H_b and H_l are symmetric Hadamards. The explicit switched order-4
factor in H_l has total sum zero, so 1^T H_l 1=0 exactly. Hence every
block-constant Boolean vector has cross energy zero, despite the cross
matrix generally not annihilating block-constant vectors.

Writing E0=n(l-1)/2 and Bn=n(sqrt(n)+sqrt(l))/2, spectral norm gives

    E0<=P<=E0+Bn,       0<=R<=Bn+n/2.

Consequently

    E0-w >= n[l-2-2sqrt(n)-2sqrt(l)]/4,
    I <= n[l-1+sqrt(n)+sqrt(l)]/4.

These are the precise numerator and denominator needed to derive the
author's c_n from an actual proportional ramp, including the term
-4e/n. The numerator is eventually positive and c_n tends to 2/5.
No assertion about the exact two endpoints is needed.

For random fair block signs the covariance is I_b tensor J_l, equal to
l times the block-constant projector. Thus the ramp forces

    Tr(Pi)>=c_n*n/l=(1/10-o(1))sqrt(n).

Separating diagonal and within-block off-diagonal entries also gives

    max_(i!=j)|Pi_ij| >= [c_n-Tr(Pi)/n]/(l-1).

Both statements apply to positive contractions as well as projectors.
They refute both the o(sqrt(n)) trace regime and the simultaneous
o(n) trace / o(n^(-1/2)) coherence regime for arbitrary bounded-cap
signings. The original matrix has bounded normalized operator norm:
the clique block norm is l-1 and the cross norm is at most
sqrt(n)+sqrt(l).

Finally, replacing only within-block clique edges by the corresponding
Hadamard tensor signs is an actual repair with at most n(l-1)/2 edges
changed. The resulting full tensor has operator norm sqrt(n), so its
hollow cap is at most (1/2+o(1))n^(3/2). The old width is at least
E0/2=(1-o(1))n^(3/2). This verifies the essential scope distinction:
the family is easily repairable and far from globally width-minimizing.
Its block projector is sign-compatible, so the separate exact-mean
sign-compatible surgery is not excluded by this counterexample.
