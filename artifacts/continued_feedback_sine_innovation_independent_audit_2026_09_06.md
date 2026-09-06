# Independent audit of the actual sine innovation floor

Date: 2026-09-06. Full fresh read of
`continued_audit_actual_sine_innovation_floor_2026_09_06.md`, together
with the full zero-first nuclear covariance theorem and the full ordered
bounded next-gain proof. Verdict: the repaired special-case theorem passes.

This audit does not turn the result into a floor for arbitrary responses
or a new universal extremal constant.

## 1. Covariance extension and cancellation control

The needed extension of the zero-first theorem is genuinely licensed by
its proof: a finite catalog of deterministic row responses can be used
at fixed coefficient quantization and fixed clipping threshold. The
responses `h sin(tz)-a_i clip(z,-K,K)` are bounded, odd, and uniformly
Lipschitz at that stage. Refine the coefficient catalog after n, then
remove clipping using uniform integrability of averaged actual Z squares.
That integrability follows by fixed Hermite-polynomial approximation,
bounded row moments of the transported squarefree polynomial, and the
vanishing limiting averaged L2 approximation error.

The covariance comparison of C, D_aZ, and their difference identifies
the SYMMETRIC weighted cross matrix. Scalar polarization by itself would
not identify the full nonsymmetric cross covariance when a_i varies;
the repaired proof correctly makes no such claim. The diagonal after
transport by B sees only the symmetric part. Normalized nuclear norm
duality against diagonal signs therefore gives

`average_i |E eta_i (L_0)_i| -> 0`.

This absolute conclusion, not just a signed trace cancellation, is exactly
what the fixed-variance-cutoff regression quotient requires.

## 2. Positive trace and the correct row variance cap

For R_f and T as in the source, the correlation Schur multiplier gives

`R_f <= L^2 tau^2 I`, `T <= L^2 tau^2 Q`,
`max_i T_ii <= L^2 tau^2`, and `||T|| <= L^4 tau^2`.

The third inequality uses unit ROW norms of B and is stronger than an
operator-norm-only estimate. Since all source degrees p are odd,

`Tr T = sum_p f_p^2 sum_ij Q_ij^(p+1) >= n tau^2`.

For the sine residual,

`K_res=h^2 D_e [sinh(t^2 T)-t^2 T] D_e`

is an entrywise positive odd-power series in the PSD matrix T. Keep its
cubic term and use `Q >= T/(L^2 tau^2)`. The crucial trace is exactly

`Tr(T D_e T^{circ3} D_e)=sum_ij e_i e_j T_ij^4`.

Every summand is nonnegative. Keeping the diagonal and using Jensen
gives the claimed `h^2 t^6 tau^6 exp(-w)/(6L^2)` mean variance floor,
where `w=t^2 L^2 tau^2`. In particular, the tau exponent is six after
division by the tau squared from the Loewner comparison.

The Schur-power bound
`||T^{circ p}|| <= (max_i T_ii)^(p-1) ||T||`
gives `||K_res|| <= h^2 L^2[sinh(w)-w]`. Unit row norms of B make this
a bound for EVERY ideal returned row variance. Combining this upper
cap with the positive average variance and normalized nuclear actual
comparison yields a lower bound on average square roots. The proof
does not make the invalid inference that an average variance floor
alone implies an average standard-deviation floor.

## 3. Gain and feasibility

At every fixed positive variance cutoff, the averaged absolute rho
estimate makes `|sqrt(v)+rho/sqrt(v)|` asymptotically interchangeable
with sqrt(v). The missing small-variance square roots cost at most
sqrt(epsilon). The ordered scalar gain theorem then applies with the
ACTUAL mask `J=1-|f+C|`, bounded below by the declared positive slack s.

For the explicit bounded zero-first source in the source note, choose
h=1/4, s=1/2, and t=1/(L tau). Substituting w=1 gives exactly

`sqrt(2/pi) exp(-1) / [48 sqrt(sinh(1)-1)] * L^(-9)`.

The two continuous-cube endpoints have means
`u +/- (1-|u|) sign(Bu)`. Hollow multilinear rounding makes their
cross gain a valid lower certificate for
`max_x |x^T Bx|/(2n)=Qabs(A)/(n sqrt(n-1))`.

This is a verified positive actual second-query gain on every sequence
with the fixed operator cap. It specifically demonstrates that total
innovation cancellation is avoidable for a concrete admissible source
and response; it does not claim a positive pointwise variance, a positive
floor for all first stages, or an improvement of the preserved bound.
