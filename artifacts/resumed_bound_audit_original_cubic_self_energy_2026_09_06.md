# Independent audit: cubic old-response self-energy on bounded-op signings

Date: 2026-09-06. Outcome: PASS for Sections 1--2 of
`resumed_response_original_bounded_op_cubic_and_transport_audit_2026_09_06.md`,
with bounded functions understood as Gaussian-a.e.-continuous or through
their expressly ordered L2 realizations.

Let B=A/sqrt(n-1) be hollow symmetric sign-flat with ||B||op<=L fixed,
and Q=B^2. Then Q_ii=1, |Q_ij|<=1, Tr Q=n, and Tr Q^2<=L^2 n.
For bounded odd old F and r^2=||P1F||_2^2, the claimed formula is

    E[F(X)^T B F(X)]/(2n)
        =(r^2/2) Tr(B^3)/n+o(1).

Its proof does not assume convergence of the cubic moment.

For a fixed odd polynomial, write F=L_F+R. The exact injective input
Gram theorem plus the small output-root collision covariance gives
Cov(L_F)=r^2 Q+o_op(1), so its trace is exactly the cubic term plus o(n).
The full nonlinear covariance theorem gives a normalized-nuclear error
and main sum of odd Schur powers Q^(circ k), k>=3. Hollowness and
flatness of B give

    |Tr(B Q^(circ k))|/n
      <=[1/(n sqrt(n-1))] sum_(i!=j)|Q_ij|^k
      <=L^2/sqrt(n-1).

The fixed-op bound controls the nuclear error. Thus the pure nonlinear
self-energy is o(n).

The cross term uses only the endpoint Frobenius lemma, not an unproved
cross-covariance nuclear estimate. Independently reconstructing that
lemma: equal spin sets force the own-input root a to be a marked vertex
of the forest and exclude its unmarked output root j from both sets.
There are q-1 free labels and 2q-1 matrix factors, giving prefactor
n^(-1/2). A free-free parity edge gains another n^(-1/2) by the bounded
bilinear norm, so its covariance matrix is O_F(1). Otherwise the only
parity graph is an optional a-j edge and p paths a-u-j. The forest
root has exact parity degree k>=3, hence p=k-e>=2. The matrices are
proportional to

    n^(-1/2)(Q^(circ p)-I), or B circ Q^(circ p).

Each has operator norm O(n^(-1/2)), hence Frobenius norm O(1). The
display explicitly deletes the impossible diagonal a=j. This corrects
a harmless diagonal precision in the earlier endpoint note; without
deleting it the error is O_F(1), not o_F(1), and the needed bound still
holds. That earlier display has now been corrected.

The fully injective forest root-map has bounded operator norm; its
marked-collision deletion is a common input-tensor projection. Combining
this with Cov(output-root collision)=O_op(1/n) yields the remaining
O_F(1) cross-covariance term. Therefore Cov(BR_main,X_T)=O_F(1).
Its diagonal trace is O(sqrt(n)), and raw normalized-L2 forest errors
are removed by Cauchy--Schwarz. This proves the first/nonlinear cross
energy is o(n).

Bounded F follows by a fixed odd-polynomial approximation after the
matrix limit; |Tr(B^3)|/n<=L controls its first-chaos coefficient error.
The same ordered approximation for bounded even H gives
E Q_B(SH)/n->0, from the own-input covariance ||H||_2^2 I+o_op(1).
Its exact creation field B[SH]=UH+o_L2 follows from the input/output
relation and requires no involution hypothesis.

Finally the marked cross identity j=E[H U^*F] and the two feasible
means +/-F+SH yield

    Q(B)/n >= |j|+(r^2/2)|Tr(B^3)/n|-o(1),

since max(|a+j|,|a-j|)=|a|+|j|. This is an expectation comparison;
it does not identify new energy variances, nonlinear transports of BF,
or the return Q[SH] with SH. The latter replacement still requires a
separate involution-type hypothesis.

## Additional audit: a genuine fractional-slack gain

`resumed_response_original_fractional_rounding_gain_2026_09_06.md`
also passes independent reconstruction. Condition on an arbitrary mean
u in the cube and independently round to v. For d=1-||u||^2/n, each
row's centered field has exact variance
[nd-(1-u_i^2)]/(n-1). Smoothing absolute value at scale n^(-1/6) and
second-order Lindeberg replacement gives its stated uniform O(n^(-1/6))
expectation error: summed third moments are O(n^(-1/2)), and the smooth
third derivative is O(epsilon^(-2)). Variance-zero rows cause no issue.

Hollowness makes v_i independent of its own field conditional on u.
The Gaussian excess Jensen bound and ||Bu||^2/n<=L^2(1-d) therefore
give a positive fixed gap g_0 whenever conditional slack d>=d_0>0.
The step g_0/(4L) is feasible, since L>=1 and g_0<=sqrt(2/pi), and
gains g_0^2/(8L)-o(1). Both initial and final independent rounding
preserve the quadratic energy expectation exactly.

The stated limitation is also correct. Damping a Boolean endpoint
of positive normalized energy e by a common factor creates slack d
but loses d e. This particular Gaussian-only certificate gains at most
d/(4pi L), which cannot pay the loss at e=.433322... and L>=1.
No claim is made that a more informative field calculation cannot help.

## Additional audit: mean-preserving correlated Gaussian rounding

The director's subsequent stronger rounding theorem also passes, with
the following fully explicit error bound. For any fixed cube mean u,
sample G with covariance I+B/L and threshold its coordinates to means
u_i. Put phi_i=phi(Phi^(-1)((1+u_i)/2)), with value zero at endpoints.
The first Hermite coefficient is -2phi_i under the lower-threshold
convention, so its pair product is positive 4phi_i phi_j. All higher
chaos contributes a remainder bounded by

    |R_ij| <= (B_ij/L)^2 sqrt((1-u_i^2)(1-u_j^2)).

Consequently its absolute normalized signed energy contribution is at
most 1/(2L^2 sqrt(n-1)). The exact first-chaos normalized payoff is

    2/[L n(n-1)] ((sum_i phi_i)^2-sum_i phi_i^2)
      >=(2/L)(average phi_i)^2-1/[pi L(n-1)].

Concavity of the Gaussian isoperimetric profile gives
phi_i>=phi(0)(1-|u_i|)>=phi(0)(1-u_i^2)/2. Hence, writing
d=1-||u||^2/n, this fresh correlated Boolean rounding v obeys

    E Q_B(v)/n >= Q_B(u)/n+d^2/(4pi L)
          -1/[pi L(n-1)]-1/[2L^2 sqrt(n-1)].

The covariance is PSD and has unit diagonal because B is hollow and
||B||op<=L. Singular covariance and deterministic threshold endpoints
are allowed by the L2 Hermite expansion or continuity. Conditioning on
any prior random mean u is harmless when this Gaussian draw is fresh.
No claim about an adaptively transported matrix field is invoked.
The pure banked endpoint has d=0; this theorem by itself does not
improve its universal constant.
