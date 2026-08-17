# Dense sign bridges can carry exponential response information at the
# `n^(3/2)` scale

Status: main-agent probabilistic theorem draft awaiting independent audit.

The degree-one packing obstruction does not decide whether density itself
regularizes a bridge.  It does not: a single dense sign bridge can expose
exponentially many independent response coordinates with macroscopic
`n^(3/2)` margins, provided the child landscape class is unrestricted.

## Theorem RD.1 (random dense bridge packing)

There are universal constants `a,gamma,L>0` and, for every sufficiently
large `n`, a sign matrix

```math
B_n\in\{-1,1\}^{n\times n},                                  \tag{RD.1}
```

spin states `x_1,...,x_N` and spin queries `y_1,...,y_N`, with

```math
N\ge\exp(gamma n),                                           \tag{RD.2}
```

such that, for `sigma in {0,1}^N`, the `2^N` Boolean-indexed landscapes

```math
h_sigma(x_c)=a n^(3/2)sigma_c,
\qquad
h_sigma(x)=-L n^(3/2)\quad(x\notin\{x_1,...,x_N\})            \tag{RD.3}
```

obey

```math
(P_Bh_sigma)(y_c)
=\max_x\{h_sigma(x)+x^TBy_c\}
=D_c+a n^(3/2)sigma_c,                                      \tag{RD.4}
```

where `D_c=x_c^TBy_c` is independent of `sigma`.  Consequently their
bridge-response functions are pairwise separated by exactly `a n^(3/2)` in
sup norm.  Every summary that answers all these bridge continuations to
error `<a n^(3/2)/2` therefore needs at least

```math
2^N=2^(exp(Omega(n)))
```

states, or `N=exp(Omega(n))` bits.

### Proof

Take `B` with independent Rademacher entries and take independent uniform
`y_1,...,y_N in {-1,1}^n`.  Put

```math
x_c=sign(By_c),\qquad D_c=||By_c||_1,                         \tag{RD.5}
```

with either sign at a zero coordinate.

We use three standard Rademacher estimates.  There are universal positive
constants `c_i` such that:

1. for a universal `C_0`,
   `P(||B||_(2->2)>C_0sqrt(n))<=2exp(-c_1n)`;
2. for fixed `y`,

   ```math
   P(||By||_1<d_0n^(3/2))<=2exp(-c_2n)                       \tag{RD.6}
   ```

   for one fixed `d_0>0`; and
3. conditional on `B,y_d` and on `||B||_(2->2)<=C_0sqrt(n)`, for
   `c!=d`,

   ```math
   P(x_d^TBy_c>d_1n^(3/2))
   <=exp(-d_1^2n/(2C_0^2)).                                  \tag{RD.7}
   ```

For completeness, (RD.6) follows because the rows of `By` are independent
Rademacher sums, `E|sum_(j<=n)epsilon_j|>=sqrt(n/2)` by
Khintchine's inequality, and each centered absolute sum is subgaussian with
scale `O(sqrt n)`; take any `d_0<1/sqrt2`.  The operator estimate follows
from the usual two-net argument and the scalar Rademacher tail.  For (RD.7),
conditional on `B,y_d`, the coefficient vector is `B^Tx_d`, whose Euclidean
norm is at most `C_0n`.  Hoeffding's inequality for the independent signs of
`y_c` gives (RD.7).

Choose `0<d_1<d_0`, then choose `gamma_0>0` smaller than one quarter of all
three exponential rates, and set `N=floor(exp(gamma_0 n))`.  For large `n`,
this is at least `exp(gamma_0n/2)`; rename `gamma_0/2` as `gamma` in (RD.2).
A union bound over
the `N` diagonal events and `N(N-1)` ordered off-diagonal events shows that
with positive probability

```math
D_c>=d_0n^(3/2),
\qquad x_d^TBy_c<=d_1n^(3/2)\quad(d!=c),                     \tag{RD.8}
```

and `||B||_(2->2)<=C_0sqrt n` simultaneously.  In particular the `x_c` are
distinct.

Fix such a realization and choose

```math
0<a<(d_0-d_1)/2,
\qquad L>C_0+d_0+a.                                           \tag{RD.9}
```

At query `y_c`, every listed competitor `x_d`, `d!=c`, scores at most
`(d_1+a)n^(3/2)<d_0n^(3/2)`, while `x_c` scores at least the latter amount
even when `sigma_c=0`.  Every unlisted state scores at most

```math
-Ln^(3/2)+||B||_(infinity->1)
\le(-L+C_0)n^(3/2),                                           \tag{RD.10}
```

because `||B||_(infinity->1)<=n||B||_(2->2)`.  It also loses to `x_c`.
This proves (RD.4).  Conversely, max-plus transfer is sup-norm
nonexpansive, so
`||P_Bh_sigma-P_Bh_tau||_infinity<=||h_sigma-h_tau||_infinity
=a n^(3/2)`.  A differing exposed coordinate attains equality.  Coordinate
separation and the standard packing argument prove the information lower
bound. `square`

## Scope

This is an intrinsic response packing, not a lower bound against one proof
architecture: the query transform itself contains `N` independently exposed
coordinates at the target scale.  It proves that no all-landscape theory can
compress arbitrary dense sign bridges merely from density, spectral norm, or
the `n^(3/2)` normalization.

It does **not** prove that the motivating near-minimizers are incompressible.
The constructed child landscapes are deliberately programmable and need not
be quadratic sign landscapes or near-minimizers.  A route back to the
original problem must therefore exploit rigidity of that special child
class, a restricted query family, or a mechanism other than uniform bridge
response compression.

## Corollary RD.2 (extensive information already for linear children)

There are universal `c,gamma>0`, dense sign matrices `B_n` with
`||B_n||_(2->2)=O(sqrt n)`, and `N>=exp(gamma n)` linear landscapes

```math
h_c(x)=-x^TBy_c,
\qquad ||By_c||_2=O(n),                                      \tag{RD.11}
```

whose bridge responses are pairwise at least `c n^(3/2)` apart.  Hence
uniform error below `c n^(3/2)/2` needs `N` states, or `Omega(n)` bits, even
on this linear subclass.

### Proof

In the same random experiment, decrease `gamma` so that all query pairs have
Hamming distance at least `delta n` for one fixed `delta>0`.  Conditional on
such a pair, the rows of `B(y_c-y_d)` are independent twice-Rademacher sums
of length at least `delta n`.  Khintchine plus subgaussian concentration and
a union bound over `N^2` pairs give, simultaneously,

```math
||B(y_c-y_d)||_1>=c n^(3/2)\quad(c!=d).                      \tag{RD.12}
```

The operator event gives the norm bound in (RD.11).  But

```math
(P_Bh_c)(y_c)=0,
\qquad
(P_Bh_d)(y_c)=||B(y_c-y_d)||_1.                              \tag{RD.13}
```

At the reverse query `y_d` the two displayed values are interchanged.  Their
response difference therefore takes both signs with magnitude at least
`c n^(3/2)`.  The packing holds even projectively, after arbitrary scalar
calibration, and the packing argument finishes the proof. `square`

This corollary rules out the objection that Theorem RD.1 is hard only because
its children are arbitrary lookup tables.  Its lower bound is smaller--linear
rather than exponential bits--but already extensive in interface size.  The
fields in (RD.11) have the natural root-mean-square scale `sqrt n`; they are
not fixed-magnitude sign-quadratic landscapes.

## Why this differs from the fixed-rank packing

The fixed-rank lower bound charges exponential information in the feature
rank.  Here the bridge has dense sign entries and operator norm `O(sqrt n)`;
the exposed margins themselves are order `n^(3/2)`.  Thus this is the sharp
negative endpoint of the bridge hierarchy:

```math
fixed rank -> finite-dimensional roof,
structured full rank -> synchronized count quotient,
generic dense sign bridge + unrestricted children -> exponential bits.
```
