# Independent audit: quenched orientation cavity

Date: 2026-09-07. Verdict: PASS for the corrected positive-edge version of
`decisive_independent_quenched_orientation_cavity_2026_09_07.md`.
No thermodynamic limit or cross-order comparison is concluded.

The added scalar field costs at most `h E|g|` in log partition, uniformly
over A and the magnitude profile. For two replicas at the same g, the
orientation mismatch probability is `sech^2(v)/2`. Integrating against
a translated normal density bounded by `1/(h sqrt(2pi))`, and using
`integral sech^2(v) dv=2`, gives exactly `1/(h sqrt(2pi))`. Replacing
`s1 s2` by 1 in a bounded replica observable costs at most twice that.
The signing is selected before g throughout; this quantifier is retained.

For a positive edge, optimality gives `E artanh(t A r)<=0`. Dividing
by t and summing the absolutely bounded remainder gives
`A E r<=t^2/[3(1-t^2)]`. Consequently the sign-selection defect is
either zero or `2|E r|`, and is at most twice that bound. Expanding
`log(1+t A r)` retains `-t^2 E r^2/2`; it cannot be replaced by a
pointwise sign condition. At t<=1/2 the defect plus logarithm remainder
is at most `(8/9+2/3)t^3<2t^3`.

The exact full-cavity identity is

`A R=(t+y)/(1+t y)=t+y-t y^2+t^2(y^3-y)/(1+t y)`.

Combining the same defect bound with `|y^3-y|<=1` gives a remainder
below `(8/9+2)t^2<3t^2`. Thus the stated differential formula and its
integrated O(n^(-1/2)) error at fixed beta are correct. Replacing
cavity squared correlations by full squared correlations has another
O(t) error per edge before multiplication by t. The latter squares
are tested in one common full Gibbs measure, so synchronization is
applicable without a cavity-measure mismatch.

For the weighted radial budget, even convexity in simultaneous radial
scaling gives `sum lambda_e A_e E R_e>=0`. Substitution of

`A R=t+(1-t^2)[y-t y^2/(1+t y)]`

and `A E r=-|E r|+d_e` gives the left side of (14), while the defect
contribution is at most
`sum lambda_e(1-t_e^2)2 kappa_(t_e)=(2/3)sum lambda_e t_e^2`.
This checks its exact constant as well as the unweighted specialization.

## Endpoint issue found and corrected

At t=0 an arbitrary unused sign has no optimality constraint. An exact
counterexample to the original unqualified statement is a four-vertex
path with three equal positive magnitudes lambda and zero magnitude
on its endpoint edge. Every path signing is gauge equivalent and is
globally optimal for this profile. Given the orientation s, the
endpoint correlation is `s tanh(lambda)^3`, so the signed cavity
correlation is `r(g)=tanh(lambda)^3` for EVERY g. Choosing the unused
endpoint sign positive makes `A E r+|E r|>0`, despite kappa_0=0.

The author patched the statements to positive edge magnitudes. At an
interior zero of a nonnegative C1 magnitude path its derivative is
zero, so omission does not change the summed differential formula.
Square-root endpoints can be treated from the positive interior with
the already stated integrable bounds. Zero-weight terms also vanish
from the radial budget. The corrected theorem is unaffected.

Neither the l1 budget for the absolute-MEAN kernel nor synchronization
proves a cross-versus-within comparison for that kernel. The signed
one-replica matrix need not be PSD. The common orientation at a
block-diagonal endpoint also remains common: it does not factor into
independent absolute child partitions. These two gaps are correctly
kept explicit by the canonical artifact.
