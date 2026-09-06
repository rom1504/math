# Audit: untruncated old local rate and cap-only energy transport

2026-09-06. **PASS** after full independent reading of
`transfer_fresh_untruncated_local_rate_and_energy_transport_2026_09_06.md`.
The archived qualitative local theorem is correctly credited. The new
rate and energy consequences do not establish untruncated feedback
covariance or a Gaussian law for the literal returns.

## 1. Cubic kernel and all normalization factors

Expanding `Y_i=sum_a B_ia S_a h_2(G_a)` gives coefficient sqrt(2) times
the three possible center terms at each unordered triple. Dividing by
six orders gives exactly `1/(3sqrt(2))` in the displayed symmetric kernel.

For `F_i(a,b,c)=B_ia B_ab B_ac`, the a-versus-(b,c) kernel has norm at
most L/sqrt(m), since the unweighted row-product Gram is Q circ Q and
Schur multiplication by the correlation matrix Q is contractive.
The other cut factors as B times the disjoint-row lift of norm at most
1/sqrt(m). The a=b and a=c exclusions already vanish; b=c is a column
projection for the first cut and a rectangular pinching complement
for the other cuts, costing at most two. Symmetrizing proves
`kappa_i<=sqrt(2)L/sqrt(m)` with no hidden growing-L factor.

The exact covariance gives
`v_i=1-3/m+(2/m)(B^4)_ii`. As `(B^4)_ii>=1` and is at most L^2<=m,
the stated `1-1/m<=v_i<=3` and error bound pass.

## 2. Gaussian brackets and interpolation

With `Y=I_3(k)`, its derivative is `3 I_2(k_a)`. Therefore

```
Gamma_YY=(1/3) sum_a [3 I_2(k_a)]^2
 =3 I_4(sym(k contract_1 k))
  +12 I_2(k contract_2 k)+6||k||HS^2.
```

The constant equals v. Orthogonality gives the coefficients 216 and
288 in the two variance contributions. Both contraction squared norms
are at most kappa^2||k||HS^2<=kappa^2/2, giving 252 kappa^2. The cross
brackets are `Gamma_YG=I_2(k contract_1 b)` and three times that in the
opposite order, so their variances are bounded by 2 kappa^2 and
18 kappa^2, respectively.

The inverse-generator bracket is not symmetric, but the stated
interpolation differentiates with precisely that bracket and the
symmetric Hessian. The sum of its standard deviations is
`sqrt(252)+4sqrt(2)`. Including the variance adjustment from v to one,
the proposed `32 M2 L/sqrt(m)` bound is conservative. Singular covariance
at a small order causes no inverse-covariance issue in this interpolation.

## 3. Boolean replacement and the rectangle rate

The cubic coefficient bound gives max influence at most 9/m for Y;
its total influence is 3v<=9. Adding G gives max at most 10/m and
total at most ten. A degree-two derivative vector on every
Boolean/Gaussian hybrid obeys the stated L4 bound `3sqrt(I_a)`.
Third-order Taylor replacement therefore costs at most

```
(27/6)(1+2sqrt(2/pi)) * 10sqrt(10) * M3/sqrt(m)
 <400 M3/sqrt(m).
```

The elementary spectral bootstrap `L^2<=maxentry(B) beta(B)` follows
directly from an eigenvector at its largest coordinate, and yields
`L^2<=4 C n/sqrt(m)` here. Thus the cap-only smooth-test rate is as
displayed. At mollification width eta=n^(-1/12), the boundary and
second-derivative terms are O(n^-1/12), while the third-derivative term
is smaller. Rectangle conventions are handled by the same uniform
boundary strips.

## 4. Uniform row approximation and energy-only transport

Fixed moments of G and Y are uniformly bounded independently of L:
their degrees are one and three, and variances at most one and three.
The local law therefore extends to the fixed polynomial-growth test
`|f-p|^2` by uniform integrability, retaining the SUPREMUM over output
roots. No degree cutoff grows with n.

For random Hilbert vectors on their actual common probability space,
Grothendieck and beta(B)<=4Cn give

```
|E X^T B Z|/n <=4 K_G C (sup_i||X_i||2)(sup_j||Z_j||2).
```

Choosing `X=sign(B(F-P))` proves the full-parent L1 transport. Choosing
the two actual sign-weighted mask error vectors gives the two-source
bound `4 K_G C delta(2+delta)+o(1)`. Neither step needs independence or
centering. The SUPREMUM row L2 approximation is essential; an averaged
L2 statement could not be substituted here.

These are ENERGY/L1 conclusions. They do not show small averaged
squared L2 norm of the transported error, nor nuclear covariance
convergence, nor a small-ball law for BF.

## 5. Forced-root core cross-check

Forcing the designated root into a principal core adds a row/column
matrix of operator norm at most sqrt(n). On that core the old linear
and cubic Fourier projections scale EXACTLY as sqrt(r) and r^(3/2),
where `r=(|R|-1)/(n-1)`. The inner restricted quadratic centering is r,
not one. The omitted-label errors are bounded by the sum of the removed
seed influences. This independently gives the already archived
qualitative uniform-root law, and creates no retention factor in the
subsequent full-parent energy transport.
