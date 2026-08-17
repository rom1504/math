# Independent audit of `short_seed_gram_broadcast.md`

**Verdict: PASS, with interpretation and verifier scope qualifications.**

I independently reconstructed the small-bias sampler, the alternating-form
distance estimate, the finite-field `6k`-wise-independent family, the moment
and net argument, and the response packing.  Every displayed theorem
constant is valid.  No mathematical repair is required before
canonicalization.

Two qualifications should remain conspicuous.  First, “explicit” here means
uniformly computable by finite search, not polynomial-time deterministic:
the proof gives a high-success-probability short random seed, but certifying a
good seed by the stated method checks all `2^h` modulations.  Second, the
`O(k log k)` count is **shared public description**, excluding the `h` hidden
bits and the exponentially large declared same-support query language.  The
result still does not give a disjoint exact-sign composition or a compact
query decoder.

## 1. Small-bias existence and repetition

For a fixed nonzero character, Hoeffding for independent variables in
`[-1,1]` gives

```math
\Pr\{|s^{-1}\sum_jX_j|>1/8\}\le 2e^{-s/128}.
```

There are `2^r-1` characters, so the failure probability is bounded by
`exp((r+1)log 2-s/128)`.  With `s=256r` its exponent is
`(r+1)log 2-2r<0`, including the endpoint `r=2` (indeed already `r=1`).
Thus SG.1 is correct.

Storing the ordered list of `s` points takes `sr=256r^2=k` bits.  Repeating
each indexed occurrence `r` times changes neither its empirical measure nor
the stored list, so it does not multiply this cost by another factor of
`r`.  Multiplicities and equal labels cause no issue below because an
alternating form vanishes on every diagonal pair, including distinct
occurrences carrying the same point.

## 2. Radical mass and unordered support

If nonzero `B` has rank `d`, then `d>=2` and its radical has codimension
`d`.  Fourier expansion of the radical indicator gives exactly

```math
\mu(\operatorname{rad}B)
\le2^{-d}+(1-2^{-d})/8\le11/32.
```

For `p` outside the radical, `q -> B(p,q)` is a nonzero character, whence

```math
\mathbb E_{p,q}(-1)^{B(p,q)}
\le {11\over32}+{1\over8}{21\over32}
={109\over256}.
```

The exact consequence is

```math
\Pr\{B(p,q)=1\}\ge147/512>1/4.
```

The draft deliberately rounds this down to `1/4`.  In the repeated ordered
list the diagonal contribution is zero and every supported unordered pair
is counted twice, hence

```math
2|\{i<j:B(p_i,p_j)=1\}|\ge k^2/4.
```

This proves SG.7.  Applying it to the difference of two forms proves
injectivity and relative distance at least
`(k^2/8)/binom(k,2)=k/(4(k-1))>1/4`.  There is no ordered/unordered factor
missing.

## 3. Finite-field independence and seed accounting

At the theorem's smallest order, `k=1024`; in general
`6k<=binom(k,2)<=q`.  Evaluations of a uniformly random polynomial of degree
below `t=6k` at any at most `t` distinct points are independent uniform
elements of `F_q`, by the Vandermonde map.  The absolute trace
`F_(2^d)->F_2` is a nonzero linear functional (also when `d` is even), so
the traced values are independent fair bits.  Multiplication by a fixed
word `chi_B` preserves this property.

The seed has exactly `td=6k ceil(log_2 E)` bits.  A canonical edge ordering,
field representation, and repetition rule need at most the displayed
`O(log k)` auxiliary description (and can instead be fixed uniformly).

The orthogonal-character lower bound SG.14a is also correct.  Characters
indexed by subsets of size at most `floor(t/2)` are mutually orthogonal,
because the symmetric difference of two such subsets has size at most `t`.
The support of any exact `t`-wise-independent sign distribution therefore
has at least `sum_(j<=floor(t/2)) binom(E,j)` points.  For `E=Theta(k^2)` and
`t=6k`, its logarithm is `Omega(k log k)`.  This establishes optimality for
fixed-length exact-independent sample spaces, not for every conceivable
spectral pseudorandom generator; the draft states this distinction.

## 4. Moment, net, and union-bound constants

For fixed `B,z`, write `X=sum_(i<j)c_(ij)epsilon_(ij)` with
`c_(ij)=2z_i z_j chi_B(i,j)`.  Every monomial of `X^(2m)`, `m=3k`, uses at
most `2m=6k` distinct edge variables, so its expectation agrees exactly
with full independence.  Rademacher-to-Gaussian moment comparison gives

```math
\mathbb E X^{2m}
\le(2m-1)!!\left(\sum c_e^2\right)^m,
\qquad
\sum c_e^2=2(1-\sum_i z_i^4)\le2.
```

Thus

```math
\mathbb E X^{2m}\le{(2m)!\over m!}\le(2m)^m,
```

and Markov at `4 sqrt(k)` yields `(3/8)^(3k)`.  A `1/4` Euclidean net has
size at most `9^k`; for a symmetric matrix the quadratic-form net lemma
multiplies the threshold by `(1-2/4)^(-1)=2`.  Therefore the simultaneous
operator threshold is exactly `8 sqrt(k)`.

Finally,

```math
{1\over k}\log\left(9^k2^h(3/8)^{3k}\right)
\le\log9+{\log2\over512}+3\log(3/8)
=-0.743909\ldots,
```

so the advertised `<exp(-0.74k)` failure probability has room.  Dependence
between different `B` events is irrelevant to the union bound.

## 5. Response-gap normalization

For `B!=T`, the difference `D=A_B-A_T` has magnitude two on at least
`m=k^2/8` unordered edges.  A vertex bipartition cuts at least `m/2` of
them.  If `d_i` are the supported cut degrees on one side, sharp real
Khintchine at `p=1`, followed by `d_i<=k`, gives

```math
\max_{x,y}\left|\sum_{i,j}D_{ij}x_i y_j\right|
\ge \sqrt2\sum_i\sqrt{d_i}
\ge {\sqrt2\over\sqrt k}\sum_i d_i
\ge {m\over\sqrt{2k}}.
```

Flipping all spins on one side reverses the cross term while leaving the
two internal terms fixed, so one of the two full quadratic values has
absolute value at least this cross value.  Hence

```math
Q(D)\ge {\sqrt2\over16}k^{3/2}.
```

This step incurs neither an omitted factor two nor ordinary polarization.
On the upper side, `||A_B||_op,||A_T||_op<=8 sqrt(k)`, so

```math
Q(A_B)\le(k/2)8\sqrt k=4k^{3/2},
\qquad
Q(D)\le(k/2)16\sqrt k=8k^{3/2}.
```

The already-declared query indexed by `B` has response zero on child `B`
and at least the lower gap on child `T`.  Thus the `2^h` response vectors are
pairwise separated.  A deterministic worst-case summary with uniform error
strictly below half the gap has at least `2^h` states, hence at least `h`
bits.  Since

```math
h/k=(r-1)/(512r)\ge1/1024
```

for `r>=2`, SG.25 is exact.

## 6. Uniformity, scope, and verifier

Enumerating all candidate `S` costs at most `2^k`; enumerating all seeds,
modulations, and exact integer-matrix spectral tests costs
`2^{O(k log k)}`.  The analytic existence proof guarantees termination.
Thus “uniformly computable” is justified.  In canonical prose I recommend
using that phrase rather than unqualified “explicit,” whose conventional
meaning sometimes includes deterministic polynomial time.

The bundled verifier passes.  It checks the numerical exponents, full-space
small-bias examples through `r=4`, and a four-wise polynomial-trace toy over
`F_8`.  It does **not** construct a theorem-scale `1/8`-biased list, search a
good `6k`-wise seed, compute the simultaneous spectral event, or instantiate
the response packing at `k>=1024`.  Those parts are proved analytically and
should not be described as computationally certified by the script.

The theorem remains a same-support additive-overlay statement.  Its
response queries may have coefficients in `{0,+-2}` after addition, and its
query language has `2^h` members.  It neither gives an appended disjoint
exact-sign parent nor proves closure for arbitrary dense quadratic systems.
These are scope limitations, not defects in SG.1--SG.4.
