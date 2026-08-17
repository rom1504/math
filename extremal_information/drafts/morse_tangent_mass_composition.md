# Morse tangent mass closes the logarithmic response scale

Status: rigorous fixed-query discrete-Laplace theorem, plus a strict
finite-parameter Gaussian tangent semigroup.  The general functional-carrier
extension is not proved.  The analytic estimate is classical; the project
contribution is the response-state diagnosis and the closed Gaussian
benchmark.

## 1. A response carrier with one subleading speed

Let `D_A,D_B` be compact subsets of `R^d`. For lattice points
`k in nD_A intersect Z^d`, let `A_n(k)>0`, and define `B_n` similarly,
extending both arrays by zero outside their declared supports.
Assume continuous functions `f,g`, real exponents `alpha,beta`, and a
constant `K>=1` such that, at every feasible lattice point in the full
convolution fibre,

```math
K^(-1)n^alpha e^(n f(k/n))
 <=A_n(k)<=K n^alpha e^(n f(k/n)),

K^(-1)n^beta e^(n g(k/n))
 <=B_n(k)<=K n^beta e^(n g(k/n)).                     \tag{TM.1}
```

For `z` with `nz in Z^d`, set `D_z=D_A intersect(z-D_B)` and

```math
C_n(nz)=\sum_k A_n(k)B_n(nz-k),                       \tag{TM.2}
```

where the sum is over feasible lattice points. Put

```math
F_z(x)=f(x)+g(z-x),
\qquad h(z)=\max_x F_z(x).                            \tag{TM.3}
```

Say the exposed fibre at `z` is uniformly Morse if `F_z` has one maximizer
`x_z` at distance at least `eta` from the boundary of `D_z` and there are
`c,C>0` such that

```math
F_z(x)<=h(z)-c||x-x_z||_2^2
```

on `D_z`, while

```math
F_z(x)>=h(z)-C||x-x_z||_2^2
\quad\hbox{when }||x-x_z||_2<=eta.                    \tag{TM.4}
```

These inequalities follow locally and globally from a unique interior
nondegenerate maximum on a compact smooth domain; stating them directly
avoids irrelevant differentiability assumptions.

### Theorem TM-A (Morse tangent-mass composition)

Under (TM.1)--(TM.4), there are constants `0<c_z<=C_z<infinity`, independent
of `n`, such that, for every sufficiently large admissible integer `n` with
`nz in Z^d`,

```math
c_z n^(alpha+beta+d/2)e^(nh(z))
 <=C_n(nz)<=
C_z n^(alpha+beta+d/2)e^(nh(z)).                      \tag{TM.5}
```

Equivalently, at the two declared speeds `(n,log n)`,

```math
nu(C_n(nz))=(h(z),alpha+beta+d/2).                    \tag{TM.6}
```

Thus the pair consisting of the leading roof and its polynomial prefactor
exponent closes under a `d`-dimensional nondegenerate convolution saddle.
The missing `d/2` is precisely the tangent counting mass of the exposed
fibre.

#### Proof

By (TM.1), `C_n(nz)` lies between `K^(-2)n^(alpha+beta)` and
`K^2n^(alpha+beta)` times

```math
S_n=\sum_k e^(nF_z(k/n)).                              \tag{TM.7}
```

The first inequality in (TM.4) gives

```math
S_n<=e^(nh(z))\sum_(k in Z^d)
 exp(-c||k-nx_z||_2^2/n)
<=C'e^(nh(z))n^(d/2),                                 \tag{TM.8}
```

where the last bound is the product of `d` elementary Gaussian lattice-sum
bounds. For the reverse inequality, take every lattice point satisfying
`||k-nx_z||_2<=sqrt(n)`. There are at least `c'n^(d/2)` such points for all
large `n`; they remain feasible by interiority, and the second inequality in
(TM.4) makes every summand at least `e^(nh(z)-C)`. This proves (TM.5) and
hence (TM.6). `square`

If (TM.1)--(TM.4) hold with common constants for grid queries `z=m/n` in a
compact set, the constants in (TM.5) are uniform.  This gives uniform
logarithmic-order control.  It does not by itself give a reusable exact-
amplitude functional carrier; that would require uniform remainders,
derivative bounds, a smooth saddle map, and closure of the output class.

### Theorem TM-B (the tangent-amplitude composition law)

Assume in addition that `f,g` are `C^3` near `x_z,z-x_z`, that

```math
J_z=-D^2F_z(x_z)
```

is positive definite, and that positive continuous amplitudes `a,b` satisfy,
uniformly near the saddle,

```math
A_n(k)=n^alpha e^(nf(k/n))(a(k/n)+o(1)),
\qquad
B_n(k)=n^beta e^(ng(k/n))(b(k/n)+o(1)).                \tag{TM.9}
```

where both `o(1)` errors are uniform relative errors in a fixed saddle
neighbourhood.  Then, along the admissible integers, the complete tangent
mass composes by

```math
C_n(nz)\sim
{(2pi)^(d/2)a(x_z)b(z-x_z)\over sqrt(det J_z)}
n^(alpha+beta+d/2)e^(nh(z)).                           \tag{TM.10}
```

Hence one fixed-query Morse saddle has the explicit one-step rule

```math
\begin{aligned}
h(z)&=\max_x\{f(x)+g(z-x)\},\\
gamma&=alpha+beta+d/2,\\
c(z)&={(2pi)^(d/2)a(x_z)b(z-x_z)\over sqrt(det J_z)}.
\end{aligned}                                         \tag{TM.11}
```

#### Proof

The quadratic upper bound makes the contribution outside any fixed saddle
neighborhood exponentially negligible. Inside it, Taylor expansion gives,
for `k=nx_z+sqrt(n)y`,

```math
n(F_z(k/n)-h(z))=-{1\over2}y^TJ_zy+o(1)               \tag{TM.12}
```

uniformly for `y` in a growing bounded window, while the tails are dominated
by a Gaussian. The shifted `n^(-1/2)Z^d` Riemann sum therefore converges:

```math
n^(-d/2)\sum_k e^(n(F_z(k/n)-h(z)))a(k/n)b(z-k/n)
\longrightarrow
a(x_z)b(z-x_z)\int_(R^d)e^(-y^TJ_zy/2)dy.             \tag{TM.13}
```

The Gaussian integral is `(2pi)^(d/2)/sqrt(det J_z)`, proving (TM.10).
`square`

For finitely many uniformly separated nondegenerate maximizers, the same
proof sums their amplitudes.  On a full-rank tangent lattice `L`, (TM.10)
has the additional factor `covol(L)^(-1)`.  Thus exact order-one asymptotics
need tangent mass and lattice density, not just the logarithmic exponent.
No finite-state or minimality claim is made for arbitrary continuum fields
`(f,a)`.

## 2. Localization and known counting models

### Lemma TM-C (local saddle asymptotics suffice)

The conclusions of Theorems TM-A--TM-B remain valid if their two-sided input
asymptotics are assumed only on a fixed neighbourhood of the saddle, provided
that:

1. the input arrays have a global upper bound of the form polynomial times
   `exp(nf)` and `exp(ng)`; and
2. outside that neighbourhood, `F_z<=h(z)-eta` for some `eta>0`.

Indeed the polynomially many lattice points outside the neighbourhood have
total mass `exp(n(h(z)-eta/2))` for all large `n`, while the local saddle
contributes `Theta(n^(alpha+beta+d/2)e^(nh(z)))`.  The statement is uniform
when the polynomial degree, gap, and local asymptotics are uniform.

### Corollary TM-D (Vandermonde and multinomial shells)

For `p` in a compact subinterval of `(0,1)` and admissible integers with
`np in Z`, Stirling's uniform bounds near the saddle give

```math
{n choose pn}=Theta(n^(-1/2)e^(nh(p))).                \tag{TM.14}
```

The elementary global entropy bound and strict concavity give the
off-neighbourhood exponential gap in Lemma TM-C.  The binary entropy saddle
in Vandermonde's convolution is one-dimensional, so Theorem TM-A returns

```math
-1/2-1/2+1/2=-1/2,                                   \tag{TM.15}
```

exactly repairing the pointwise lexicographic prediction `-1`.

More generally, empirical-type multiplicities for a fixed `q`-symbol
alphabet, in the first `q-1` integer count coordinates and on admissible
types, live in dimension `d=q-1` and obey

```math
{n!\over\prod_i(np_i)!}
=Theta(n^(-d/2)e^(nH(p)))                             \tag{TM.16}
```

uniformly on compact subsets of the simplex interior. Concatenating two
fragments is a `d`-dimensional convolution, and Theorem TM-A again returns exponent
`-d/2`. Thus the same carrier closes both binary code-shell counts and
finite-alphabet mean-field type counts.

Boundary summands are handled by the method-of-types upper bound and the same
strict-concavity localization.  This does not reprove the sharp Stirling
constants. It explains why the
polynomial entropy exponent is stable under composition and identifies the
dimension term that a pointwise rate roof omits.

## 3. A strict finite-parameter tangent semigroup

The arbitrary functional tuple `(f,alpha,a)` has uncontrolled description
complexity.  There is, however, a natural closed finite-parameter class.

For a positive-definite `d by d` matrix `P`, vector `mu`, and parameters
`c,alpha in R`, `a>0`, define the full-lattice array

```math
G_n^(c,mu,P,alpha,a)(k)
=n^alpha a\exp\left{nc-{1\over2n}(k-nmu)^TP(k-nmu)\right},
\qquad k in Z^d.                                      \tag{TM.17}
```

Let `Sigma=P^(-1)` and replace `a` by its total-mass amplitude

```math
m=a{(2pi)^(d/2)\over\sqrt{det P}}.                    \tag{TM.18}
```

### Theorem TM-E (Gaussian tangent semigroup)

For two arrays of the form (TM.17), their convolution at `ell in Z^d` is,
uniformly when the parameters range over a compact set with the precision
eigenvalues bounded away from zero and infinity and `ell/n` stays compact,

```math
(G_n^theta*G_n^phi)(ell)
=(1+o(1))G_n^(theta star phi)(ell),                    \tag{TM.19}
```

where in the coordinates `(c,mu,Sigma,alpha,m)` the product is

```math
(c,mu,Sigma,alpha,m) star (c',mu',Sigma',alpha',m')
=\left(c+c',mu+mu',Sigma+Sigma',
       alpha+alpha'+{d\over2},mm'\right).             \tag{TM.20}
```

The operation is associative and commutative.  It is a strict
finite-dimensional composable response state of dimension
`d(d+1)/2+d+3`; on a bounded parameter region its total log covering number
is `O_d(log(1/epsilon))` at coordinate precision `epsilon`.

#### Proof

Completing the square in the exponent gives the output leading profile

```math
c+c'-{1\over2}(z-mu-mu')^T
 (P^(-1)+P'^(-1))^(-1)(z-mu-mu'),                    \tag{TM.21}
```

with saddle Hessian `P+P'`.  The parameter-uniform Gaussian lattice sum, or
Theorem TM-B with its uniform quadratic tails gives output power
`alpha+alpha'+d/2` and local amplitude

```math
a_out={(2pi)^(d/2)aa'\over\sqrt{det(P+P')}}.           \tag{TM.22}
```

The determinant identity

```math
det(P+P')=det(P)det(P')det(P^(-1)+P'^(-1))             \tag{TM.23}
```

turns (TM.22) exactly into `m_out=mm'`.  Formula (TM.20) is therefore the
claimed update.  Its associativity and commutativity are now coordinatewise.
Uniformity follows from Gaussian domination and compactness of the declared
parameter region.  The dimension and grid-cover bound are immediate.
`square`

This theorem is not just a rephrasing of a transfer table: an infinite
lattice array and every future convolution are represented by finitely many
parameters, and the representation is closed without enumerating the
descriptor grid.  It is nevertheless a deliberately narrow benchmark.  It
does not claim that arbitrary Morse profiles, codes, or dense quadratic
landscapes synchronize to a Gaussian tangent state.

### Corollary TM-F (finite integer recovery on compact query sets)

Fix compact parameter and query sets as in Theorem TM-E.  There are a box radius
`R<infinity` and a common `C<infinity` such that the finite integer arrays

```math
\widehat G_n^theta(k)
=1_(||k-nmu||_infinity<=Rn)
 floor(e^(nC)G_n^theta(k))                            \tag{TM.24}
```

obey, uniformly on the query set,

```math
(widehat G_n^theta*widehat G_n^phi)(ell)
=(1+o(1))G_n^((theta+C) star (phi+C))(ell),             \tag{TM.24a}
```

where `theta+C` means shifting only the leading `c` coordinate by `C`.
Thus the finite-parameter semigroup has an all-order realization by finite
integer-weight landscapes.

Indeed compactness keeps every relevant saddle a linear distance inside one
common pair of boxes.  Increase `C` so both unrounded factors are at least
`e^(delta n)` throughout one common saddle neighbourhood.  There

```math
0\le XY-floor(X)floor(Y)\le X+Y,
```

so the relative loss is at most `2e^(-delta n)`.  Outside that neighbourhood,
the Gaussian quadratic gap makes the entire unrounded product sum
exponentially smaller than the main saddle mass; the rounded product is no
larger.  The truncated tails have the same bound.  This proves (TM.24a).
It is abstract integer recovery, not realization inside a constrained code
or graph class.

## 4. Falsifiers and the boundary of the class

The Morse hypothesis is structural information. It cannot be discarded.
For a one-dimensional quartic saddle,

```math
\sum_(|k|<=n)exp(-n(k/n)^4)
=Theta(n^(3/4)),                                      \tag{TM.25}
```

not `Theta(n^(1/2))`; a flat interval gives `Theta(n)`. Multiple isolated
Morse saddles preserve the `d/2` exponent but require their amplitudes if one
wants the order-one constant. Saddles whose Hessian rank changes with the
query can therefore force a stratified state.

The finite semigroup theorem is generative only on its Gaussian class.  For
general Morse profiles `(h,alpha)` is insufficient, while arbitrary
functional amplitudes have uncontrolled information.  A finite
stratification by saddle type is a research target, not a proved universal
repair.

## 5. Literature boundary

The estimate is a finite-dimensional discrete Laplace method; local lattice
central-limit methods have long been used for coefficients and convolutions
(for example [Moran](https://doi.org/10.2307/3213083)). Recent discrete
Laplace formulations make the same nondegenerate Hessian and lattice-density
factor explicit. The present statement is recorded because the earlier
multi-speed response carrier isolated exactly the missing `d/2` term and
because Corollary TM-F is a scoped compositional repair on a natural model class,
not because the Gaussian-sum estimate itself is new.
