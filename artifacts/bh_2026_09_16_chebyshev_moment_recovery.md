# Chebyshev recovery: archived upper theorem, BH overhead, generic sharpness

2026-09-16. **No new original signing bound or convergence claim.** The
exact moment-to-cap recovery theorem below ALREADY appears in
`blank_slate_direct_attack_2026_08_21.md`, Section 2.1, Theorem 2.1.
That section and its complete proof were reconstructed, including its
parity-lattice moment LP and the warning that the remaining asymptotic
moment problem is equivalent to the original convergence question.

The present note makes the high-degree filter comparison with the new
BH paper explicit and adds a self-contained matching lower example for
a GENERIC finite atomic moment oracle. The latter is not realized by
full sign quadratics or by dyadic energy histograms.

## 1. Exact recovery from endpoint mass and moments

Let Y be a finite real random variable, `R=||Y||infinity>0`, and suppose
`P(|Y|=R)>=p>0`. If another probability measure nu supported on [-c,c]
has the same even moments through order 2d, then

```
R sech(arcosh(p^(-1/2))/d) <= c.                       (1)
```

Indeed `T_d(y/c)^2` is an even polynomial of degree 2d. Its expectation
under nu is at most one. Moment matching identifies that expectation
with its expectation under Y, which is at least

```
p cosh^2(d arcosh(R/c))                 when R>c.
```

Solving proves (1); the case c>=R is immediate. This is the archived
Chebyshev proof, not a new consequence of polynomial BH.

For the actual signing energy `Y=H_A/n^(3/2)`, global spin reversal
gives endpoint mass at least `p_n=2^(-(n-1))`. Therefore the minimum
central support radius compatible with its moments satisfies

```
R sech(L_n/d) <= L_d <= R,
L_n=arcosh(2^((n-1)/2)).                               (2)
```

This remains valid if the competing measure is restricted to the correct
energy parity lattice, as in the archived LP. A cheaper scalar test can
also be defined without solving that LP:

```
F_d(c)=E T_d(Y/c)^2,
C_d=inf{c>0:F_d(c)<=1}.
```

For every feasible c, the same argument gives `c>=R sech(L_n/d)`;
every c>=R is feasible. Thus C_d obeys (2) as well. Monotonicity of
F_d is NOT assumed: it can oscillate as c varies. The infimum definition
is important. Everything here uses exact moments, with no claim of
efficient or numerically well-conditioned moment computation.

When `d/n -> infinity`, uniformly for caps R in a fixed bounded range,

```
0<=R-L_d <=R[1-sech(L_n/d)]=O((n/d)^2).                (3)
```

For `d=alpha n`, the multiplicative lower factor tends to
`sech(log(2)/(2alpha))`. In comparison the ordinary power proxy
`P_d=(E|Y|^(2d))^(1/(2d))` only satisfies

```
R exp[-(n-1)log(2)/(2d)] <=P_d<=R,
R-P_d=O(n/d)                      for d/n -> infinity. (4)
```

This is the classical Chebyshev/Lanczos square-root acceleration, not
the new paper's mechanism. A primary reference for the analogous
power-versus-Lanczos logarithmic/squared-logarithmic behavior is
[Kuczynski--Wozniakowski, SIAM J. Matrix Anal. Appl. 13 (1992)](https://epubs.siam.org/doi/10.1137/0613066).
The elementary proof above establishes the precise statement used here;
no numerical-linear-algebra theorem is imported into the signing problem.

## 2. Linear-degree filters detect a fixed false cap without BH

Fix c>0 and suppose the actual cap obeys `R>=c+epsilon`, epsilon>0.
Put `alpha=arcosh(1+epsilon/c)`. Then

```
||T_d(Y/c)||_2 >=sqrt(p_n) cosh(d alpha).              (5)
```

Consequently the hypothetical cap R<=c is contradicted already by
Parseval/L2 whenever

```
d > arcosh(p_n^(-1/2))/alpha.                         (6)
```

At fixed c,epsilon this is

```
d = n log(2)/(2 alpha)+O(1/alpha).
```

For small epsilon/c, `alpha~sqrt(2epsilon/c)`. Thus the degree scale
is n/sqrt(epsilon), with the indicated fixed-c interpretation. Neither
BH nor the coefficient signs are required for detection at this scale.

### Exact comparison of old and new BH overhead

The scalar filter has degree d in Y, but its Walsh degree is at most
`m=min(2d,n)`. This distinction is essential, especially for d>n.
Under the hypothetical cap R<=c, the full BH coefficient test would
require `||hat(T_d(Y/c))||_(2m/(m+1))<=B_m`. Since its coefficient
ell_p norm dominates its ell_2 norm, a sufficient contradiction is

```
sqrt(p_n) cosh(d alpha)>Bbar_n,
```

where Bbar_n bounds all B_m for m<=n. Equivalently it suffices that

```
d>arcosh(p_n^(-1/2) Bbar_n)/alpha.                    (7)
```

The earlier subexponential degree bound gives
`log Bbar_n=O(sqrt(n log(n+1)))`; the new polynomial bound gives
`log Bbar_n=O(log(n+1))`; these are (1.2) and Theorem 1.1 of
[Ivanisvili's paper](https://arxiv.org/html/2609.12427v1).
Therefore they add, respectively,

```
O(sqrt(n log n)/alpha)       and       O(log n/alpha)
```

to the same leading `n log(2)/(2alpha)` scalar degree. The elementary
L2 test adds only O(1/alpha). Hence the new BH bound improves an
unnecessary overhead; it is not what makes linear-degree filtering
complete. These are contradictions to a CAP-CONDITIONED coefficient
test, never violations of the actual BH inequality itself.

The weighted square-function version has the same conclusion after
removing its constant mode. Its lower bound is
`S(f)>=m^(-5) sqrt(Var f)`. For actual Y,
`E Y^2=(n-1)/(2n^2)`, so at least half its mass lies in [-c,c] for
large n. Comparing that mass with an endpoint atom gives

```
Var(T_d(Y/c)) >=(p_n/2)[cosh(d alpha)-1]^2.
```

Thus the paper's `S(f)<=K m^22` test is contradicted once the last
square root exceeds K m^27, again an O(log n/alpha) overhead over
the L2 scale. No missing constant coefficient is being smuggled into
the weighted seminorm.

## 3. Matching generic atomic moment ambiguity

**This section concerns arbitrary finite positive measures. It does NOT
construct two signings, satisfy all quadratic-chaos identities, or
realize the measures as dyadic Boolean energy histograms.** It proves
that the scale in (3) cannot be improved using only generic moment and
minimum-atom information.

Let nu be the arcsine probability measure on [-1,1]. The polynomials
`1,sqrt(2)T_1,...,sqrt(2)T_d` are orthonormal, and their kernel is

```
K_d(z,x)=1+2 sum_(j=1)^d T_j(z)T_j(x).
```

Fix R>1. There are two finite positive probability measures, each with
d+1 atoms, whose moments agree through degree 2d:

* mu_0 is ordinary Gauss--Chebyshev quadrature: nodes
  `cos((2j-1)pi/[2(d+1)])`, j=1,...,d+1, each of weight 1/(d+1).
  Its support is inside [-1,1].
* mu_R has one node at R and d nodes inside (-1,1), namely the roots
  of `K_d(R,x)`. Every quadrature weight at a node x_j is
  `1/K_d(x_j,x_j)`.

In particular,

```
mu_R({R})=1/K_d(R,R),
every interior weight >=1/(2d+1).                    (8)
```

Here is a self-contained proof of the second construction. Kernel
reproduction shows, for every polynomial q of degree at most d-1,

```
integral (R-x)K_d(R,x)q(x) dnu(x)=R q(R)-(xq)(R)=0.
```

Thus K_d(R,x) is an orthogonal polynomial for the positive measure
`(R-x)dnu(x)`. It has d simple roots in (-1,1). Use these roots and R
as interpolation nodes and define the quadrature weights by integrating
their degree-d Lagrange polynomials. Division by
`(x-R)K_d(R,x)` proves exactness through degree 2d: the quotient has
degree at most d-1, whose product integrates to zero by the identity
above, and the remainder has degree at most d.

Exactness on each squared Lagrange polynomial gives positive weight
`w_j=integral L_j(x)^2 dnu(x)>0`. Moreover, for any polynomial q of
degree at most d, exactness of `L_j q` gives
`integral L_j q=w_j q(x_j)`. Kernel reproduction then identifies
`L_j=w_j K_d(x_j,.)`, hence `w_j=1/K_d(x_j,x_j)`. On [-1,1],
`K_d(x,x)<=2d+1`, proving (8). This is a prescribed-node positive
quadrature/Christoffel-kernel construction; no new BH input is involved.

### Exact endpoint scale

For every `0<p<1/(2d+1)` there is a unique R=cosh(alpha)>1 satisfying

```
K_d(R,R)=1/p.
```

The kernel is strictly increasing on [1,infinity), equals 2d+1 at one,
and diverges. Both measures above then have EVERY atom at least p,
identical moments through 2d, and respective caps at most one and R.
Since

```
(1/2)exp(2d alpha)
 <=K_d(cosh(alpha),cosh(alpha))
 <=(2d+1)exp(2d alpha),
```

we have, with `L=log(1/p)`,

```
[L-log(2d+1)]/(2d) <=alpha<=[L+log2]/(2d).            (9)
```

If `log d=o(L)` and `L=o(d)`, then

```
alpha~L/(2d),        R-1~L^2/(8d^2).
```

Thus two allowed atomic laws have exactly the same first 2d moments
and a cap gap on the same scale as the Chebyshev upper theorem. In
fact the first cap is `cos(pi/[2(d+1)])=1-O(d^(-2))`, so their cap
gap is also asymptotic to `L^2/(8d^2)`. This agrees with the leading
one-sided recovery loss `1-sech(arcosh(p^(-1/2))/d)`.
Taking `p=2^(-(n-1))`, with polynomially bounded d>>n, gives the
matching generic order `(n/d)^2`. Any moment-only point estimator must
err by at least half their endpoint gap on one of these two measures.
This is a sharpness statement for the stated oracle class, not for the
much narrower class of actual signing energies.

### Symmetric-energy variant

The generic ambiguity persists even if both energy laws are required
to be invariant under `y -> -y`. Assume `2p<1/(2d+1)`, first choose
the outer quadrature atom to have mass `2p`, and then replace each law
mu by `(mu+mu reflected)/2`. The moments still agree through degree
2d, the absolute caps are unchanged, and every atom has mass at least
p: reflection only halves an original atom or merges such masses.
There are at most `2(d+1)` nodes in each law. Replacing L by
`log(1/(2p))=L-log2` leaves the cap-gap asymptotic unchanged.

This extra energy symmetry is NOT the spin-antipodal identity
`H_A(-x)=H_A(x)`, and neither symmetry realizes the quadratures as
dyadic Boolean histograms or as full sign quadratics.

## 4. Status and replay

`computations/bh_mechanism_2026_09_16_chebyshev_recovery.py` checks exact
integer/Fraction Chebyshev recurrences, moment-versus-value equality, and
endpoint amplification on actual small full signings. It separately
checks the generic quadrature construction numerically, explicitly
labeling that regression as non-certified rather than an exact theorem.
The mathematical proof in Section 3 supplies its exact justification.
The canonical replay is
`computations/results/bh_mechanism_2026_09_16_chebyshev_recovery_replay.json`:
40 exact actual-signing tests and five floating-point generic quadrature
regressions, all PASS. The source also passes `py_compile`.

The existing archive already supplies the stronger parity-lattice
moment-LP formulation and its convergence equivalence. This note adds
classification and a generic sharpness example, not a new minimax
lower/upper bound or a new solution of the original problem.
