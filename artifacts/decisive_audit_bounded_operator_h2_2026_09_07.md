# Independent audit: bounded-operator H2 instability

Date: 2026-09-07. PASS for
`decisive_independent_bounded_operator_h2_counterexample_2026_09_07.md`.

The exact operator class is closed under H2 because, up to coordinate
ordering, `L(A)+I=H2 tensor(A+I)+diag(0,2I)` and
`5sqrt(2)-2>=5`. The banked restricted-weave representation yields
`||A+I||/sqrt(N)<=sqrt(32/31)+o(1)` at all orders, so the claimed strict
subhalf class upper bound has a concrete dependency, not an assumed
regularization of arbitrary original minimizers.

For a lift maximizer let T be the smaller equal/opposite support. The
support estimate follows exactly from
`q_L(x,y)=2q_A(x)-4q_A(x_T)+|S|-|T|`, or the exchanged version, and
`|q_A(x_T)|<=||A|| |T|/2`. Thus both supports are macroscopic using the
stated universal lower constant.

The real intermediate two-block target has cap at most
`(1+rho^2)k^2/2`: its diagonal terms cancel between the equal-size blocks,
and maximizing `.5a^2-.5b^2+rho ab` on the square gives this value.
Its old lift witness has energy `2(1+rho)k^2-2k`, so the gain ratio is
exactly `2(1+sqrt(2))`. The final matrix is an actual signing; the target
is used only for edge-flip probabilities.

Scalar Bernstein with deviation parameter 4k and the union bound over
the 2k-dimensional cube gives the stated `sqrt(epsilon) k^(3/2)+k`
noise cap. Chernoff controls every changed degree by `4epsilon k`, whose
failure probability is negligible even for `theta>=c/log N`.
Consequently the perturbation operator norm is at most `8epsilon k`.
The old selected-block cap is `O(k sqrt(N))`; after multiplication by
epsilon it is `O(N)`. Polarization bounds the linear lift noise by six
times its cap. All errors therefore have the displayed
`O(sqrt(theta) N^(5/4)+(1+theta)N)` size.

Finally bounded dyadic restricted-class minima have increments arbitrarily
close to nonnegative. The positive gadget gain survives on those orders.
The vanishing-theta variant is only near-minimal relative to the fixed
operator class. It does not assert that this class attains the original
M_N asymptotically, and hence does not construct bounded-op original
near-minimizers with that instability. The unconditional conclusion of
bounded-op, strict-subhalf H2 counterexamples is valid.
