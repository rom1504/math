# Fixed energy-power tilts do not lower the leading physical response

2026-09-17. Bernoulli-track derivation, prompted by the director's
actual-energy-to-cheap-response question. This is a scoped negative
theorem for genuine full signings, not a statement about all
energy-dependent sampling laws. No external novelty claim is made.

## 1. Exact square-tilt Walsh identity

Let A be a full symmetric hollow signing on n>=5 vertices and

```
H(h)=sum_{i<j} A_ij h_i h_j,
d=binom(n,2),
mu_n=E|epsilon_1+...+epsilon_n|,
m=n-1 if n is even, and m=n if n is odd.
```

Under the actual physical sign law with density H(h)^2/d relative to
the uniform cube, its mean absolute response to a Boolean x is EXACTLY

```
E_tilt |h dot x|
 =mu_n {1+[(m-1)W(x)-H(x)^2+d]/[d m(m-2)]},          (1)
W(x)=||Ax||_2^2-n(n-1).
```

To prove it, the square H^2 has only Fourier levels 0,2,4. Its
constant is d and its level-two polynomial is

```
(H^2)_2(x)=2sum_{i<j}(A^2)_ij x_i x_j
         =||Ax||_2^2-n(n-1)=W(x).
```

The remaining level-four polynomial is H(x)^2-d-W(x). For the
absolute-overlap kernel K(h)=|sum_i h_i|, its Fourier coefficients on
levels zero, two and four are respectively

```
mu_n,    mu_n/m,    -mu_n/[m(m-2)].                  (2)
```

Here is an elementary check of the two nonconstant coefficients.
Conditioning on n-2 unused signs, the second difference of absolute
value is 1 at a remaining sum zero (even n), or 1/2 at remaining
sum +-1 (odd n), and zero otherwise. This gives the level-two
coefficient by a central binomial probability. Taking one further
pair difference gives the level-four coefficient: its ratio to
level two is -1/(n-3) for even n and -1/(n-2) for odd n. These are
exactly (2), using the elementary central-binomial formula for mu_n.
Convolution on the cube multiplies equal Fourier coefficients, giving
(1). No unspecified higher levels enter.

If ||A||_op<=L sqrt(n) and Q(A)<=c n^(3/2), formula (1) is uniformly

```
E_tilt |h dot x|/sqrt(n)=sqrt(2/pi)+O_(L,c)(1/n).     (3)
```

In particular the H(x)^2=Q(A)^2 information on an absolute ground
word produces only an O(n^(-1/2)) change in its UNNORMALIZED response.
The separate row-square term W(x) is retained exactly, not assumed
small from the cap alone. Section 2 handles bounded caps without any
operator-norm hypothesis and all fixed even powers.

## 2. Uniform theorem for every fixed positive power

Fix an integer k>=1 and c<infinity. Suppose A is any full signing
with Q(A)<=c n^(3/2). Define the genuine sign-column probability law

```
nu_(A,k)(h)=2^(-n) H_A(h)^(2k)/E H_A(epsilon)^(2k).
```

Its denominator is positive. Uniformly over EVERY Boolean x,

```
|E_(nu_(A,k)) |h dot x|/sqrt(n)-sqrt(2/pi)|
 <=C_(k,c) n^(-1/4).                                (4)
```

Thus these laws cannot yield any fixed leading response discount,
even on the actual absolute ground code. No flat spectral assumption
and no optimizer-specific randomness is used. Constants may depend
on k; (4) makes NO claim for k growing with n, sharp conditioning on
near-ground energy, or exponential Gibbs tilts. Section 3 extends it
to arbitrary fixed-degree positive polynomial densities, including
n-dependent coefficients, while keeping degree fixed.

### 2.1 The actual cap bounds linear--quadratic interaction

For Boolean x,y put u=(x+y)/2 and v=(x-y)/2. Multilinearity and
averaging over independent signs with means u_i or v_i show
|H_A(u)|,|H_A(v)|<=Q(A). Polarization therefore gives

```
|y^T A x|=2|H_A(u)-H_A(v)|<=4Q(A),
||Ax||_1<=4Q(A),
||Ax||_2^2<=||Ax||_infinity ||Ax||_1
          <=4(n-1)Q(A)<=4c n^(5/2).                 (5)
```

This uses the actual Boolean cap, not an inferred spectral estimate.

### 2.2 Moment control throughout the mixed sign/Gaussian path

Let U have independent coordinates, each either a fair sign or a
standard Gaussian. Write

```
Z(U)=x dot U/sqrt(n),    Y(U)=H_A(U)/n.
```

For every fixed r>=2, uniformly along every such mixed path,

```
||Z(U)||_r<=C_r,   ||Y(U)||_r<=C_r.                  (6)
```

For a coordinate removed from the two forms, write

```
Z=z+a U_i,    Y=y+b U_i,
|a|=n^(-1/2),    b=n^(-1)sum_{j!=i}A_ij U_j.
```

Then z,y,b are independent of U_i, and

```
||z||_r+||y||_r<=C_r,    ||b||_r<=C_r n^(-1/2).      (7)
```

These are the ordinary degree-one/two hypercontractive moment bounds:
the squared coefficient sums for Y and y are at most d/n^2<=1/2,
and that for b is (n-1)/n^2. Tensorization applies to the mixed
sign/Gaussian product law. The direct primary source checked for this
ingredient is Mossel--O'Donnell--Oleszkiewicz,
[*Noise stability of functions with low influences: invariance and optimality*](https://annals.math.princeton.edu/wp-content/uploads/annals-v171-n1-p05-p.pdf),
Propositions 3.11--3.12 and Theorem 3.13. Their Theorem 3.18 proof
under Hypothesis 4 was also read: it replaces one coordinate, cancels
three matching moments, and bounds the fourth-order remainder by
hypercontractive moments. We use that mechanism below, with a fully
written joint test-function estimate rather than a black-box claim
about powers of a Boolean multilinearization.

### 2.3 Joint weighted-absolute replacement, with its remainder paid

Put phi_delta(z)=sqrt(z^2+delta^2), 0<delta<=1. Its difference from
|z| is at most delta. Its derivatives of orders 1,...,4 are bounded
by C delta^(1-r), while phi_delta itself is at most |z|+1.
For one replaced coordinate consider

```
F(t)=phi_delta(z+a t)(y+b t)^(2k).
```

Its fourth derivative is a sum, for 0<=j<=min(4,2k), of constant
multiples of

```
a^(4-j) b^j phi_delta^(4-j)(z+a t)
                       (y+b t)^(2k-j).              (8)
```

The first three moments of a fair sign and a Gaussian agree.
Taylor expansion at t=0 therefore cancels all terms through degree
three. For the fourth-order remainders, (7), Holder's inequality,
and the fixed moments of the replaced sign or Gaussian imply

```
E[ |U_i|^4 sup_(|t|<=|U_i|)|F''''(t)| ]
 <=C_k delta^(-3)n^(-2).                            (9)
```

To see the n power explicitly, each term in (8) contains
`a^(4-j)b^j`, whose relevant Holder norm is O_k(n^(-2)); the
remaining polynomial factors have fixed bounded moments by (7).
When 4-j=0, phi_delta contributes at most |z|+|aU_i|+1, again with
bounded moments. This also justifies the supremum over the Taylor
segment for an unbounded Gaussian driver.

Summing n replacements and paying the two smoothing errors gives

```
|E_Rad |Z|Y^(2k)-E_Gauss |Z|Y^(2k)|
 <=C_k[delta+n^(-1)delta^(-3)].                      (10)
```

Choose delta=n^(-1/4). Applying the same argument to Y^(2k) alone
(no absolute-value factor or smoothing) also gives

```
|E_Rad Y^(2k)-E_Gauss Y^(2k)|<=C_k/n.                (11)
```

The polynomial Y^(2k) here is the literal power evaluated at the
mixed variables. It is NOT first reduced modulo U_i^2=1 and then
evaluated at Gaussians; that invalid substitution would change the
Gaussian law and is not part of the proof.

### 2.4 Gaussian separation using only the cap

Let e=x/sqrt(n) and decompose a standard Gaussian vector as
G=Z e+V, with Z standard normal and V independent in e-perpendicular.
Then

```
Y(G)=D+Z L+a_0 Z^2,
D=V^T A V/(2n),   L=e^T A V/n,
a_0=e^T A e/(2n)=H_A(x)/n^2.                        (12)
```

By (5),

```
|a_0|<=c n^(-1/2),
Var(L)<=||Ax||_2^2/n^3<=4c n^(-1/2),
||L||_r<=C_(r,c)n^(-1/4).                            (13)
```

Every fixed moment of D is bounded uniformly: its centered quadratic
coefficient Frobenius norm divided by n is at most
||A||_F/n<=1, and its mean is -a_0 because tr A=0.
In (12), D and L need not be independent, but BOTH are independent
of Z. Expand the fixed power 2k and use Holder with (13). It follows
that

```
E |Z|Y(G)^(2k)=kappa E D^(2k)+O_(k,c)(n^(-1/4)),
E Y(G)^(2k)=E D^(2k)+O_(k,c)(n^(-1/4)).              (14)
```

Together (10)--(14) show

```
|E_Rad |Z|Y^(2k)-kappa E_Rad Y^(2k)|
 <=C_(k,c)n^(-1/4).
```

Finally Jensen gives the dimension-free denominator bound

```
E_Rad Y^(2k)>=(E_Rad Y^2)^k
             =[(n-1)/(2n)]^k>=4^(-k)  (n>=2).
```

Divide to obtain (4). Uniformity in the tested Boolean x is explicit
in every step.

## 3. Every fixed-degree positive energy-polynomial density

There is a stronger, still fixed-degree formulation. Fix D and c.
Let P_n be ANY real polynomial of degree at most D, whose coefficients
may depend on n and A and may have either sign. Assume only

```
P_n(H_A(h)/n)>=0 for EVERY Boolean h,
E_Rad P_n(H_A(h)/n)=1,
Q(A)<=c n^(3/2).
```

The resulting genuine sign law with density P_n(H_A/n) satisfies

```
sup_(Boolean x)
 |E_tilt |h dot x|/sqrt(n)-kappa|
 <=C_(D,c)n^(-1/4).                                 (15)
```

Nonnegativity is needed only on the actual finite energy spectrum,
not on the real line. Thus (15) includes n-dependent normalized
mixtures, and it is not restricted to nonnegative monomial coefficients.

### 3.1 Bounded cap forces a typical-spin quadratic Gaussian limit

For an extremal unit eigenvector e of A with eigenvalue lambda, the
continuous cube point e/||e||_infinity and multilinear averaging give

```
Q(A)>=|lambda|/(2||e||_infinity^2).
```

The eigenvector equation also gives
`|lambda| ||e||_infinity<=sum_i|e_i|<=sqrt(n)`. Hence

```
||A||_op^3<=2n Q(A),
||A||_op/n<= (2c)^(1/3)n^(-1/6).                    (16)
```

This weak spectral estimate is sufficient here. It does NOT assert
||A||_op=O(sqrt(n)). The product-rounding spectral estimate itself is
already in the [August 21 direct-attack record](blank_slate_direct_attack_2026_08_21.md),
in its eigenvector-rounding passage (with the slightly sharper n-1
in place of n). It is reused here, not claimed as a new spectral result.
For a standard Gaussian G, diagonalization gives

```
Y_G=H_A(G)/n=sum_i [lambda_i/(2n)](G_i^2-1),
Var(Y_G)=(n-1)/(2n).
```

For every fixed r>=3 its cumulant is

```
cum_r(Y_G)=[(r-1)!/(2n^r)]sum_i lambda_i^r,
|cum_r(Y_G)|<=C_r (||A||_op/n)^(r-2).                (17)
```

The formula follows by expanding the logarithm of the product of
centered chi-square MGFs at zero. It proves convergence of every
fixed moment to that of N(0,1/2), uniformly over this cap-bounded
class. Applying the polynomial-only replacement argument from (11)
to ANY fixed integer power gives

```
|E_Rad Y^r-E_Gauss Y^r|<=C_r/n.                     (18)
```

Thus the same uniform fixed-moment limit holds for H_A(epsilon)/n.
This is a typical-spin limit, NOT a statement about the extremal tail
at energy n^(3/2). The latter lives at a diverging argument of Y.

### 3.2 Uniform conditioning of the coefficient moment matrix

Let M_(n,D) be the matrix indexed by 0<=i,j<=D with entries
E_Rad Y^(i+j). By (17)--(18) it converges uniformly to the corresponding
moment matrix M_D of a nondegenerate N(0,1/2). The latter is positive
definite: a nonzero real polynomial cannot vanish Gaussian-almost
everywhere. Therefore there are constants a_D>0 and n_0(D,c) such
that

```
v^T M_(n,D) v>=a_D ||v||_2^2   for n>=n_0(D,c).     (19)
```

This finite-dimensional statement is uniform even when the coefficient
vector v=v_n depends on the actual signing and the order.

On the Boolean cube, P_n(H_A/n) has multilinear degree at most 2D.
The fourth-moment hypercontractive inequality and Holder interpolation
give

```
||P_n(Y)||_4<=3^D ||P_n(Y)||_2,
||P_n(Y)||_2<=||P_n(Y)||_1^(1/3)||P_n(Y)||_4^(2/3),
||P_n(Y)||_2<=3^(2D)||P_n(Y)||_1=3^(2D).              (20)
```

The final equality uses nonnegativity on the cube and normalization.
Combining (19)--(20) bounds the Euclidean norm, and hence the sum of
absolute values, of all polynomial coefficients by a constant
depending only on D (for n>=n_0). Boolean multilinearization is used
ONLY to invoke its Boolean norm inequality (20); Gaussian powers in
(17)--(18) remain their literal polynomial evaluations.

### 3.3 Finish by coefficientwise joint replacement and separation

The proof of Section 2.3--2.4 did not require an even exponent except
for the monomial probability normalization. For each fixed integer
0<=r<=D it proves

```
|E_Rad |Z|Y^r-kappa E_Rad Y^r|
 <=C_(r,c)n^(-1/4).                                 (21)
```

Multiply (21) by the bounded coefficients and sum. Normalization
E P_n(Y)=1 gives (15) for n>=n_0. For the finitely many smaller
orders, a larger C_(D,c) suffices because every physical response is
at most n, independently of the polynomial representation of its
density. Thus the theorem is uniformly quantified at every order.

The fixed degree D is essential to this proof. No uniform control
of its constants as D grows is asserted, and no low-temperature or
near-ground conditioning theorem is inferred from it.

## 4. Relation to the campaign's other negative results

The older high-power/Bohnenblust--Hille analyses concerned coefficient
tests, moment recovery, or growing-degree filtering. The theorem here
concerns a DIFFERENT object: an actual positive physical sign law
obtained by fixed-degree positive energy reweighting, tested on every scalar Boolean
response. It supplies no new moment method or coefficient inequality.

It narrows the actual-energy-to-cheap-response search: a fixed power
of H cannot exploit the shared high absolute energy of the nearcode
at leading order, even though that energy forces a large quadratic-
feature covariance eigenvalue. Neither that necessary feature condition
nor this scoped failure rules out nonlinear, sharply conditioned,
or degree-growing laws.

Reproducibility: `computations/paper_bernoulli_2026_09_17_fixed_power_tilts.py`
passed 81,600 exact all-query Walsh and cap/polarization checks on 80
full signings of orders 5--12, plus 16 exact kernel-coefficient checks.
It also records finite floating diagnostics for powers 1--4 on three
explicit Hadamard principal signings; those small-order diagnostics
are not evidence of an asymptotic rate. The script passed Python
compilation, and its complete output is preserved at
`tmp/paper_portfolio_2026_09_17/bernoulli/fixed_power_tilts_audit.json`.
