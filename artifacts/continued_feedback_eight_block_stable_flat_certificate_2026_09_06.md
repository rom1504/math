# A certified eight-block bound removes a stable near-flat class

Date: 2026-09-06. Director's bounded eight-block proposal, derived and
certified here. This is a rigorous class-specific one-row upper bound,
not an upper cap for all rows and not a convergence theorem.

Throughout m=2^d, k=15m/16, p=15/16, q=1/16. Let H be the Walsh matrix.
A supported spin is xi in {0,±1}^m, with support T and weight k. Its
normalized magnitudes are a_u=|(H xi)_u|/sqrt(k). The epsilon-flat class
is `sum_u(a_u-1)^2<=epsilon²m`.

## 1. The theorem and operational gain

There are fixed epsilon0>0 and c>0, and selector events G_m with
P(G_m)->1, such that at t=4

```math
\limsup_m\frac1m\log\mathbb E_T\left[
 1_{G_m}\sum_{x\text{ epsilon0-flat}}L_4(|H[T,:]^Tx|)
 \right]
\le-4(1-\sqrt p)-c.                                    (1)
```

The same bound holds when T is conditioned on G_m. Thus the nearly flat
class is genuinely eliminated at a tilt also left open by the known
typical-Gaussian lower obstruction. This does NOT bound the complementary
rows. The certified zero-width tilted margin is at least .00874170.
The existence of a positive fixed width follows from the explicit finite
dual and the vanishing error terms proved below; no n-dependent width is
substituted for a stable neighborhood.

## 2. Four-wise marginals of a random affine eight-block

Take independent directions a,b,c in F_2^d and a uniformly random origin
x. The eight input points are x+u1 a+u2 b+u3 c, u in F_2^3. At most four
distinct vertices are either affinely independent, or, in the four-point
case only, form a parallelogram. Affinely independent images are uniform
affinely independent input points; compared with independent sampling,
the discrepancy for bounded tests is O(1/m).

Write g=1_T-p. Require G_m to satisfy `||g||_{U²}=o(1)`, uniformly over
its selectors. Such events have probability tending to one. For example,
the fourth Fourier moment of a fixed-weight random selector gives
E||g||_{U²}^4=O(1/m); Markov supplies any fixed slower vanishing cutoff.
This fourth-moment estimate follows equally by the four-sample
hypergeometric formula. It does not involve the spin choices.

With the normalized convention
`||f||_{U²}^4=m^{-4} sum_u |(Hf)_u|^4`, epsilon-flatness gives

```math
\|\xi\|_{U^2}^4
\le8p^2(\epsilon^4+1/m),\qquad
|\mathbb E_x\xi_x|\le\sqrt p(\epsilon+m^{-1/2}).         (2)
```

Indeed, putting a_u=1+e_u, use
`sum e_u^4 <= (sum e_u²)^2` and `(1+e)^4<=8(1+e^4)`.
The mean estimate follows from the single zero-frequency coefficient.

For a parallelogram, the Fourier identity and Hölder give

```math
\left|\mathbb E_{x,a,b}\prod_{j=1}^4 f_j(x+v_j)\right|
\le\prod_j\|f_j\|_{U^2},                               (3)
```

after the appropriate character signs; in F_2 these signs coincide.
Expand the positive/negative symbol indicators as `(p+g±xi)/2`.
Equations (2)--(3), and exact independence of three distinct
parallelogram vertices before the negligible degeneracy restriction,
show that every <=4 symbol marginal of the random eight-block differs
by `O(epsilon)+o_m(1)` from the product law P, where

```math
P(-1)=P(+1)=15/32,\qquad P(0)=2/32.                      (4)
```

The estimate is uniform over ALL epsilon-flat xi on G_m. No U³ bound
or independent eight-symbol law is asserted. The first four marginals
are precisely the data used in the finite optimization below.

## 3. The finite four-wise-independent inequality

Index the eight symbols by F_2^3, and let

```math
\mathcal V=\{000,100,010,001\}.
```

For X in {-1,0,1}^8, let Y be the FOUR unnormalized Walsh coefficients
with frequencies in V. If X has product law P, write P_Y for the induced
law of Y. This is an exactly enumerable finite distribution: there are
6561 input patterns, each of mass
`prod_j w(X_j)/32^8`, with w(-1)=w(+1)=15 and w(0)=2. There are 1537
distinct Y values.

For EVERY law Q on these 6561 patterns having the same <=4-coordinate
marginals as P^8, the following inequality is certified:

```math
\boxed{\mathbb E_Q\log P_Y(Y)
\le A:=-\frac{65899554490727779}{10485760000000000}
<-6.28467125804212.}                                    (5)
```

Here is a reproducible exact dual certificate, independent of the
floating optimizer that discovered it.

* The marginal basis consists of the constant and the products of
  specified ±1 indicators on up to four distinct coordinates. In the
  script's degree/coordinate/symbol order there are 1697 rows, with
  exact expectations `(15/32)^r` for a degree-r row.
* The affine-frequency stabilizer of V has 24 elements. Its induced
  signed coordinate permutations, together with global sign reversal,
  act in 48 ways and partition the patterns into 255 orbits. They
  preserve P^8, all the constraints, and P_Y(Y). The script verifies
  the last statement by EXACT integer masses, not a floating tolerance.
* Average the marginal basis over each orbit. The fixed sparse vector
  `CERTIFIED_DUAL` in the script, divided by 10^10, majorizes
  `log P_Y(Y)` on every orbit. All averaged-basis values are exact
  integers divided by orbit sizes. Its expectation is exactly A.
* Each logarithm is enclosed by rational arithmetic: range-reduce to
  [1,2), put z=(x-1)/(x+1), and use the first 35 terms of
  `log x=2 sum_{j>=0} z^(2j+1)/(2j+1)`. The remainder is at most
  `2 z^71/[71(1-z²)]`. Range-reduction multiples of log2 use the
  appropriate endpoint according to their sign. Every one of the 255
  inequalities is checked against this rational UPPER enclosure.

The checker runs WITHOUT calling an optimizer:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python computations/continued_feedback_eight_block_profile_lp_2026_09_06.py --verify-only
```

The floating LP discovery option is kept separate. Exact optimality of
the candidate dual is unnecessary. For approximately correct four-wise
marginals, the same fixed dual gives (5) with error
`O(epsilon)+o_m(1)`: there are finitely many fixed coefficients, and the
maximum marginal error in Section 2 is uniform. Large but finite dual
coefficients shrink the admissible epsilon0; they do not invalidate a
fixed positive-width conclusion.

## 4. Stable spectral encoding and exact conditional recovery

Partition the input cube into affine three-dimensional cosets. In each
block keep the four Walsh sums Y at V. Globally these are equivalent,
by inverse Walsh transforms on the quotient, to exactly m/2 full
spectral coefficients. Encode only their signs and replace magnitudes
by sqrt(k). This costs `(m/2)log2` nats. The direct-sum basis and the
chosen three directions can be encoded with o(m) extra nats.

There are M=m/8 blocks. Orthogonality of the quotient transform implies
that the squared error across all four block sums is at most

```math
M^{-1}\sum_u(|(H\xi)_u|-\sqrt k)^2\le8\epsilon^2 k.     (6)
```

Every block-sum coordinate is an integer in [-8,8]. Rounding to that
17-symbol alphabet gives at most `32epsilon² k` incorrect coordinates.
Their fraction among the m/2 stored block sums is
`e=64p epsilon²`. Encode the error locations and corrected values,
at total cost at most

```math
\tfrac m2[h(e)+e\log17]+o(m).                           (7)
```

Now the block sums Y are exact. Recover the input block using its exact
product-P conditional probability `P^8(X)/P_Y(Y)`. Conditional coding
does not assume that the actual input blocks are independent: for ANY
fixed Y list, the product of these conditional probabilities sums to
one over all compatible inputs. Hence the number of inputs having
total negative log conditional probability at most R is at most e^R.

For any actual xi of weight k, its negative log product-P probability
is EXACTLY mH3, where `H3=h(p)+p log2`; the split between + and - does not
matter. The recovery length is therefore

```math
mH_3+\sum_{\text{blocks}}\log P_Y(Y_{\text{block}}).      (8)
```

Averaging (8) over affine three-block partitions and using Section 2
and (5), some partition has length at most
`m[H3+A/8+O(epsilon)+o(1)]`. The origin chosen inside a block does not
matter: translating it changes only known signs of Y and leaves P_Y
unchanged. Thus averaging over cosets is the same as averaging over all
affine origins, as used in Section 2.

Combine (7)--(8) with the m/2 spectral sign bits and divide by the
number binomial(m,k) of selectors. This proves

```math
\limsup_m\frac1m\log\mathbb E_T[1_{G_m}N_T(\epsilon)]
\le c_8+\tfrac12[h(64p\epsilon^2)
                    +64p\epsilon^2\log17]+O(\epsilon),          (9)
```

```math
c_8=(p+\tfrac12)\log2+A/8<.210815164800.                 (10)
```

This is a genuine stable DENSE-profile count. The O(epsilon) constant
is determined by the finite dual and the elementary moment estimates,
not by an unproved limiting distribution.

## 5. Certified negative permanent exponent at t=4

The deterministic near-flat permanent bound from
`continued_feedback_stable_flat_profile_count_2026_09_06.md`, Section 3,
is uniform over every removed coordinate:

```math
\frac1m\log L_t\le\tfrac12\log[(1+e^{-4t})/2]
                         +4t\epsilon+o(1).             (11)
```

At epsilon=0 and t=4, adding the weave tilt to (9)--(11) gives

```math
p\log2+A/8+\tfrac12\log(1+e^{-16})+4-\sqrt{15}
<-.00874170918.                                         (12)
```

The displayed strict sign is certified by rational bounds in the checker:
the same rational log2 upper bound, `sqrt(15)>3.872983346207416` checked
by squaring integers, and `exp(-16)<1/8000000`. The last follows from
the first 50 positive terms of the exact rational exponential series
for exp(16). Use `log(1+z)<=z`.

All the error terms in (9) and (11) tend to zero as epsilon decreases,
after m tends to infinity on G_m. Choose a fixed epsilon0 sufficiently
small and then a fixed c>0 below the remaining margin. This proves (1).
Intersections with the previous high-probability selector filters are
allowed without changing the proof or its exponential upper rate.

## 6. Limits

The theorem excludes one stable spectral neighborhood, not all profiles
with low entropy. A near-flat count cannot be applied to every arbitrary
profile, and a finite positive dual for this 6561-pattern problem is not
a certificate for the complete one-row partition sum. The complementary
profile sum is the remaining original-signing upper-construction task.
