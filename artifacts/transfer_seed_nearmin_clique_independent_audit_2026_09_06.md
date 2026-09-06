# Independent reconstruction of the near-minimizer clique theorem

Date: 2026-09-06. Audit of
`transfer_adversary_nearmin_clique_orientation_gap_2026_09_06.md`.
Verdict: **PASS**, with the exact-versus-asymptotic and covariance-metric
scope specified below. No prior positive audit verdict is used as a step.

## 1. Spectral core reconstructed from the original cap

For hollow symmetric A, write `H_A(x)=x^T A x/2` and
`beta(A)=max_(x,y signs)|x^T A y|`. Independent rounding extends the
cap bound to the cube. With `u=(x+y)/2`, `v=(x-y)/2`, symmetry gives
`x^T A y=2(H_A(u)-H_A(v))`, so `beta(A)<=4Q(A)`.

Put `c=asinh(1)`, `K_G=pi/(2c)`. The odd tensor-series factorization
of `sin(ct)` gives two families of unit vectors with cross products
`sin(c<u_i,v_j>)`: the sum of the absolute tensor coefficients is
`sinh(c)=1`. Common Gaussian hyperplane rounding has correlation
`(2/pi)arcsin(sin(c<u_i,v_j>))=(2c/pi)<u_i,v_j>` because `c<pi/2`.
Thus every unit-vector bilinear objective is at most `K_G beta(A)`.

The simultaneous diagonal-majorant SDP has primal
`min tr D`, diagonal `D>=A` and `D>=-A`, and dual
`max tr A(X-Y)`, `X,Y>=0`, `diag(X+Y)=1`.
For Gram vectors `X_ij=<a_i,a_j>`, `Y_ij=<b_i,b_j>`, use the two UNIT
families `(a_i,b_i)` and `(a_i,-b_i)`. Their cross products are `X-Y`.
The preceding rounding therefore bounds the dual by `K_G beta(A)`.
Strict feasibility holds on both sides, so a majorant exists with
`tr D<=4K_G Q(A)`. Hollowing gives `D_ii>=0`.

At most half the coordinates can have `D_ii>8K_G Q(A)/N`.
Retain the others. Restricting the two PSD inequalities shows that
the retained principal operator norm is at most this threshold.
Every r-subset of that core consequently has cap
`u<=r ||A_S||op/2<=4K_G(r/N)Q(A)`.
This uses neither a prior lower bound nor a random-restriction result.
The finite proof was checked directly against Section 1 of
`resumed_bound_audit_minimal_proof_2026_09_06.md`.

## 2. Finite replacement and its asymptotic scale

Switch and, if needed, reverse an absolute maximizing state so that
`H_A(1)=M=Q(A)`. Overwrite the chosen r-block by a positive clique.
The energy identity is exact:
`H_(A')=H_A-H_(A_S)+((sum_S x_i)^2-r)/2`.

If P and R are the positive and negative extrema of A', respectively,
testing the OLD all-one ground gives
`P>=M+r(r-1)/2-u`. Pointwise, the clique energy is at least `-r/2`,
so `R<=M+r/2+u`. Hence
`P-R>=r^2/2-r-2u`. The triangle inequality and the tested witness give
`M+r(r-1)/2-u<=Q(A')<=M+r(r-1)/2+u`.

Now start with an exact minimizer at EVERY order and use
`r=floor(N^(2/3))`. The established upper bound `M_N=O(N^(3/2))`
implies `u=O(N^(7/6))=o(r^2)`. Therefore
`Q(A')-M_N=(1/2+o(1))r^2=o(N^(3/2))`, while
`P-R>=(1/2-o(1))r^2` and eventually `Q(A')=P`.
The modified matrices are asymptotically minimizing but NOT exactly
minimizing; indeed their unnormalized gaps diverge like `r^2/2`.
Any sequence `sqrt(N)<<r<<N^(3/4)` gives the same scope.

## 3. Isotropy and the exact normalization

For `tau=o(r^2)`, every absolute near-ground state has positive energy
at least `P-tau`; negative states stop at absolute value R, strictly
below this shell. An isotropic law `E xx^T=I` has mean hollow energy
zero and cannot be supported there.

The chord of the convex function `|h|` on `[-R,P]` gives the exact
bound `(P+R)|h|<=2PR+(P-R)h`. Under an isotropic law,
`E[P-|H|]>=P(P-R)/(P+R)>=(P-R)/2`.
Thus a small exceptional mass far from the shell does not permit
mean slack `o(r^2)` either.

For an arbitrary law on the positive shell, put `K=E xx^T`.
Frobenius Cauchy--Schwarz and flatness give

```math
\frac{\|K-I\|_F}{\sqrt N}
\ge\frac{2(P-\tau)}{N\sqrt{N-1}},
\qquad
\|K-I\|_{op}\ge\frac{\|K-I\|_F}{\sqrt N}.
```

This is the RMS deviation of the N eigenvalues, NOT the entrywise
RMS `||K-I||F/N`. The latter could still tend to zero, and the proof
does not exclude that weaker metric. K is the second-moment matrix;
if a centered covariance is desired, independent global-sign
symmetrization makes the mean zero while preserving the shell and K.

Because `P=M_N+o(N^(3/2))` and `tau=o(N^(3/2))`, every sequence of
such shell-supported laws has lower limit of both displayed metrics
at least `2 liminf M_N/N^(3/2)`, hence at least
`0.8666442233281614` by the banked lower theorem. No bound on the
operator norm of the clique-modified parent was inserted: its exact
Frobenius norm `sqrt(N(N-1))` suffices.

## 4. Replay and quantifier boundary

`computations/transfer_seed_clique_orientation_audit_2026_09_06.py`
independently checks finite replacement inequalities, the pointwise
chord, and the squared rational shell covariance inequality on random
actual signings of orders three through nine. All arithmetic in its
assertions is integer exact; it is a supplementary finite replay,
not a proof of the spectral core or of asymptotic minimality.

The theorem negates an isotropy assertion quantified over ALL
asymptotically minimizing sequences. It supplies no exact minimizing
counterexample at unbounded orders and excludes neither a specially
selected minimizing family nor a specially selected near-minimizing
family with appropriate ground laws. It is not an upward construction
from a small seed and proves neither convergence nor nonconvergence.
