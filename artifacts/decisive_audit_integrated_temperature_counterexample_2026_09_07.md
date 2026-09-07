# Independent audit: full-sign integrated temperature-payment obstruction

The proof in
`decisive_bridge_integrated_temperature_payment_counterexample_2026_09_07.md`,
including its stronger Section 5, passes fresh algebraic audit.  It concerns
cap-bounded nonminimizers, not actual optimized pressure or original minima.

For `k=2^(2r+1)`, `ell=4k`, `N=k ell`, the construction is hollow and
every off-diagonal entry is a sign.  The inner matrix
`B=[[H,-H],[-H,H]]` kills constant block vectors and has operator norm
`sqrt(2ell)`.  Thus the off-block operator bound is
`Loff<=sqrt(2N)(1+1/sqrt(k))`.  Decomposing each clique spin vector into
its mean and its orthogonal part gives exactly the source's loss
`Delta=half sum(ell^2-m_i^2)` and gain at most `(Loff/ell)Delta`.
For large `k`, the positive cap is therefore `N ell/2-N/2`.

The first-bit Sylvester split has outer children `+H_(k/2)` and
`-H_(k/2)`, with their diagonals removed.  Both `k/2` and `ell/2` are
powers of four.  The stated two order-four sign vectors have eigenvalues
`+2` and `-2`; their tensor products give both signs at every required
order.  Zero trace means hollowing preserves their quadratic energies.
An outer vector with positive quadratic value, tensored with an inner
negative eigenvector `(v,-v)`, gives exactly

    H_child=-d/2-sqrt(2)d^(3/2)/2,       d=N/2.

The outer vector need not remain an eigenvector after diagonal removal;
only its quadratic value is used, and that value is unchanged.  The
inner vector does remain an exact eigenvector and has zero clique sum.

Combining this negative witness with the coherent positive witness
cancels their opposite `d/2` corrections in the paired width.  After
division by the child normalization and summing both children, the
lower bound is `beta N(3sqrt(2))/4`.  The parent's two log-sum upper
bounds give `Nlog2+beta N(2+Loff/sqrt(N))/4`.  At `beta=4` their
difference is consequently at most

    N[log2-2(sqrt(2)-1)+sqrt(2/k)].

The limiting coefficient is strictly negative.  Meanwhile the parent's
positive cap tends to `N^(3/2)` from below and its opposite cap is at
most `N/2+N Loff/2`; hence `Q(A)<=N^(3/2)` eventually.

For a fixed parent, the entire cross-edge deletion path is a convex
combination of the parent and a block-switching conjugate.  Its cap
never exceeds the parent's.  The exact integrated second derivative
of each branch includes its full cross-statistic variance, and the
first derivative at the decoupled endpoint vanishes by child spin
reversal.  Thus the negative pressure difference really is a failure
of full integrated gain to pay child reheating on a uniformly
cap-bounded full-sign path.  No variance term has been selectively
dropped.  Global signing optimality remains an essential untested
hypothesis for any stronger payment theorem.
