# A sharp energy-only counterexample: centered Curie--Weiss ground code

2026-09-17. This is a scope test for the actual-ground-code minimax
question, NOT a counterexample in the original hollow full-sign class.
Both its nonzero constant and its non-unit edge coefficients are explicit.

## 1. Exact finite response game

Let n>=2 be even, let B_n={x in {+-1}^n: sum_i x_i=0}, and put

```
C_n=B_n union {1,-1},
f_n=E_(X uniform B_n) |b dot X|,   b ANY fixed element of B_n.
```

Permutation symmetry makes f_n independent of b. The exact minimax
over ALL probability laws nu on physical sign columns h is

```
inf_nu max_(x in C_n) E_nu |h dot x|
       = n f_n/(n+f_n).                            (1)
```

No isotropy, support-size, independence, or Gaussian approximation is
assumed in the infimum.

### The elementary permutation inequality

For a fixed physical h write t=|sum_i h_i| and
F(h)=E_(X uniform B_n)|h dot X|. Then

```
F(h)>=f_n(1-t/n).                                  (2)
```

To prove it, replace h by -h if necessary, so its sum is t>=0. Fix a
balanced b. Construct a random permutation H of h by putting +1 in
EVERY position where b=+1, and placing the remaining t/2 plus signs
uniformly among the n/2 positions where b=-1. This is possible even
at t=0,n, and

```
E H = (t/n)1+(1-t/n)b.
```

For each balanced X, Jensen gives E_H|H dot X|>=(1-t/n)|b dot X|.
Average over uniform X. Every permutation of h has the same F, proving
(2). This uses no asymptotic concentration statement.

For an arbitrary nu put a=E_nu|sum_i h_i|. Testing the all-ones word
and averaging the balanced test words gives

```
max_C E_nu|h dot x| >= max{a, f_n(1-a/n)}
                       >= n f_n/(n+f_n).           (3)
```

For the reverse inequality, use the column law that assigns mass
alpha=f_n/(n+f_n) to the all-ones word and mass 1-alpha to a uniform
balanced word. At the all-ones query its response is alpha n; at every
balanced query it is (1-alpha)f_n. Both equal (1). Global sign
symmetrization gives a centered law without changing any response.
The SAME mixture on query words is an exact minimax dual certificate.

## 2. Exact formula and asymptotics

Write r=n/2. The intersection size K of the plus coordinates of two
independent balanced words has law

```
P(K=k)=binom(r,k)^2/binom(2r,r),
their overlap=4K-2r.
```

Consequently

```
f_n=4r binom(r-1,floor(r/2))^2/binom(2r,r).          (4)
```

For a direct finite proof, use

```
(2k-r)binom(r,k)^2
 =r[binom(r-1,k-1)^2-binom(r-1,k)^2]
```

and telescope the positive half of the symmetric sum. The usual
factorial Stirling estimates, separately for even and odd r, give

```
f_n=kappa sqrt(n)+O(n^(-1/2)),  kappa=sqrt(2/pi),
n f_n/(n+f_n)=kappa sqrt(n)-kappa^2+O(n^(-1/2)).    (5)
```

Thus even the unrestricted physical-law game has NO leading discount
from the independent-sign response constant.

## 3. It is an exact absolute quadratic ground code -- outside the original class

Consider the centered weighted quadratic

```
F_n(x)=[(sum_i x_i)^2-n^2/2]/(2sqrt(n)).            (6)
```

Its absolute cap is exactly n^(3/2)/4, and its COMPLETE absolute
ground code is precisely C_n: the endpoints of the possible squared
magnetization interval are 0 and n^2. In the original edge convention,

```
F_n(x)=sum_(i<j) (1/sqrt(n))x_i x_j
       +(n-n^2/2)/(2sqrt(n)).
```

The edge coefficients are therefore 1/sqrt(n), not full signs, and
the constant is nonzero for n>=4. Removing the constant destroys the balanced
absolute-ground sector; replacing the coefficients by unit full signs
changes the cap to order n^2. Neither modification is innocuous.

Equations (1)--(6) rigorously rule out any proposed leading-discount
theorem based ONLY on bounded n^(3/2)-scale quadratic cap and membership
in a complete absolute ground code, if that theorem allows arbitrary
weighted coefficients and offsets. They do NOT rule out the desired
theorem for hollow full signings with zero constant. They also do not
construct a minimizing sequence for the original problem.

## 4. Replay

The accompanying uniquely named script checks the rational formula,
every finite inequality (2) by magnetization class, and the primal/dual
value. For small n it enumerates all physical columns and all balanced
queries. Larger checks use exact hypergeometric sums. Numerical limits
are diagnostics only; (5) follows analytically from the displayed
binomial formula.

Script: `computations/paper_bernoulli_2026_09_17_centered_energy_scope.py`.
Output: `tmp/paper_portfolio_2026_09_17/bernoulli/centered_energy_scope.json`.

## 5. Removing the constant alone does not rescue the general implication

There is also a zero-constant, bounded-coefficient weighted example.
Take even integers 4<=m<=l, write n=m+l, and set

```
a=m(m-2)/[l(l-2)]<=1,
H(x,y)=[(sum x)^2-m-a((sum y)^2-l)]/2,
x in {+-1}^m, y in {+-1}^l.                        (7)
```

This is literally a hollow quadratic with edge coefficient +1 inside
the m block, -a inside the l block, and zero across the blocks. There
is no constant in its edge-polynomial representation. Its positive
and negative extrema have equal magnitude

```
Q=[m(m-1)+a l]/2=[m+a l(l-1)]/2,                  (8)
```

because m(m-2)=a l(l-2). Its COMPLETE absolute ground code is

```
C=( {+-1_m} times B_l ) union ( B_m times {+-1_l} ). (9)
```

For any physical column law nu, write b=E_nu|sum_(l block) h_i|.
Averaging queries in the SECOND sector of (9), and using the symmetry
of the balanced first block, shows its maximum response is at least b.
Averaging the FIRST sector shows it is at least
E_nu F_l(h_l)>=f_l(1-b/l), by (2). Hence

```
inf_nu max_C E_nu |h dot z| >= l f_l/(l+f_l).       (10)
```

Conversely, take independently uniform balanced columns in BOTH
blocks. Their overlaps with the respective ferro queries are exactly
zero. The two sector responses are respectively f_l and f_m, proving

```
inf_nu max_C E_nu |h dot z| <= max(f_m,f_l).         (11)
```

Choose m to be an even integer asymptotic to n^(3/4), and l=n-m.
Equations (5), (8), (10), and (11) give

```
Q=(1/2+o(1))n^(3/2),
inf_nu max_C E_nu |h dot z|=(kappa+o(1))sqrt(n).    (12)
```

Multiplying (7) by a fixed positive 2c changes the cap to
(c+o(1))n^(3/2) without changing (9) or (12). For 0<c<=1/2, all edge
magnitudes remain at most one. Thus ZERO CONSTANT and bounded
coefficient magnitudes, even with a prescribed positive cap scale,
are not enough for a uniform leading response discount.

This still does not approach an original-model counterexample. The
cross coefficients are zero and the large block has coefficient
a~n^(-1/2), not unit magnitude. In particular

```
||A||_F^2=m(m-1)+a^2 l(l-1)~2n^(3/2),
```

whereas a genuine full signing has Frobenius square n(n-1). Replacing
these weak/zero coefficients by unit signs cannot be treated as a
small perturbation uniformly on all the ground-code words. No such
rounding or transfer is claimed. The decisive surviving restriction
in this scope test is the original full-unit dense coefficient class.

The localization track independently audited Sections 1--3 PASS,
including the exact finite minimax certificate. The n=2 zero-constant
exception to the first example's wording was corrected in that audit.
