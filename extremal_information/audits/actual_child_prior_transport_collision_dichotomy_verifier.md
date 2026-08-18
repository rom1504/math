# Independent verification of the child-prior transport conversion

**Verdict:** mathematical pass, with the operational scope narrowed in the
audited note.

For `Z=||M||_F^2/d in[0,1]`, the hypotheses give `E_qZ>=a`.  If
`p=q{Z>=a/2}`, then `a<=p+(1-p)a/2`, so `p>=a/(2-a)`.  Entropy duality and
the transport gap yield `log K_0>=ca sqrt(d)/2` on that event, and GC.3
gives the deleted inequality simultaneously for every edge.  Theorem
37.56's integrated overlap is exactly the overlap under the declared path
mixture, so the actual-path specialization and its fixed positive mass are
valid.

Failure of a uniform gap supplies only an existential MGF-violating query.
It need not have macroscopic norm, concise description, or target relevance
and may require the complete convex hull to find.  The audit now states
this limitation and does not call the witness an operational coherent
phase.  No proof or constant correction is otherwise needed.

The sector-subgaussian refinement also passes.  Conditioning, Gaussian
linearization, and the second factor bound produce the determinant with
coefficient `kappa_*`; `||A||op^2<=d` gives denominator
`1-kappa_*/4`, and comparison with `1/(4sqrt(d))` yields exactly
`sqrt(d)>=2kappa_*/(4-kappa_*)`.  The rare-spike example has covariance
norm `1+(m-1)e^(-km)` but contributes probability `exp(-kN)/2` to the
aligned rank-one atom, proving the stated nonlinear-MGF violation.  It is
correctly labeled nonactual and is used only to refute covariance-only
closure.

The optimized-tilt refinement PT.1c also passes.  With
`K=max{1,kappa_*}` and `theta=1/(2sqrt(Kd))`, the determinant argument costs
at most `||M||_F^2/(6d)`.  Entropy duality yields
`[1/(2sqrt(K))-1/(6sqrt(d))]||M||_F^2/sqrt(d)`, which is at least
`||M||_F^2/(4sqrt(Kd))` exactly when `d>=4K/9`.  At comparable splits,
`kappa_*=O(N^(2alpha))`, `0<=alpha<1`, therefore gives the stated
`Omega(N^(1-alpha))` positive-mass collision tail.
