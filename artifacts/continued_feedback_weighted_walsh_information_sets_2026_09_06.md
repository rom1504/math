# Weight-aware Walsh information sets and the remaining profile gap

Date: 2026-09-06. An elementary averaged supported-spin count, followed by
an explicit obstruction to using only the currently controlled profile
statistics. The obstruction is to a RELAXATION, not to the actual weave.

Let m=2^d, let the Walsh frequencies and input coordinates be F_2^d,
and let T be uniform among k-subsets, k/m->p>2/3.

## 1. A weight-aware replacement for the 3^s span count

For a fixed Fourier support U of size s, write E_U for its real Walsh
span. There is a coordinate set I, |I|=s, such that restriction to I is
injective on E_U, by row pivots of a basis matrix. Translation of the
input variable preserves E_U, since it multiplies each Fourier character
by a sign. Therefore EVERY translate I+v is also an information set.

If xi is ternary and has weight k, the average number of zero coordinates
on I+v, over all m translations v, is `(1-k/m)s`. For at least one v this
number is at most its floor. The restriction to that information set has
at most that many zeros and otherwise two sign choices. Injectivity gives
the exact bound

```math
\#\{\xi\in E_U\cap\{0,\pm1\}^m:|\operatorname{supp}\xi|=k\}
\le m\sum_{j\le\lfloor(1-k/m)s\rfloor}\binom sj2^{s-j}.   (1)
```

This is a union bound over translates, not an assertion that each vector
has the same favorable information set. Since p>2/3, the relevant ternary
binomial summands are increasing up to their endpoint and their rate is

```math
s\{h(p)+p\log2\}+o(m).                                  (2)
```

For comparison, `h(p)+p log2<log3` strictly at this p. Unlike the general
3^s argument, (1) uses translation invariance of Walsh spectral subspaces.

## 2. An actual averaged count of supported spin rows

Let N_T(s) count all spin rows on T whose Walsh spectrum has at most s
nonzero coefficients. Each full ternary vector contributes to exactly
one selector, and every support of size at most s is contained in an
s-element frequency set. Consequently

```math
\mathbb E_T N_T(s)
\le\frac{m\binom ms}{\binom mk}
       \sum_{j\le\lfloor(1-k/m)s\rfloor}\binom sj2^{s-j}.  (3)
```

At s/m->q this has exponent at most

```math
c_p(q)=h(q)+q\{h(p)+p\log2\}-h(p),                       (4)
```

also bounded by the trivial p log2. Thus (3) is a genuine profile COUNT,
not merely a typical-spin or good-selector probability statement.
Markov's inequality also makes (4) an improved selector anti-sparsity
criterion wherever c_p(q)<0. Conditioning on any good-selector event
of probability tending to one preserves this exponential upper bound,
since the denominator changes by a factor 1+o(1).

## 3. Its direct one-row use, and why that use is insufficient

For exact support-size classes with a deterministic permanent estimate
`L_t<=exp{m ell_t(q)+o(m)}`, (3) gives

```math
\mathbb E_T\sum_{x:\,|\operatorname{supp}\hat x|\le qm}
 L_t(|\hat x|)
\le\exp\{m[c_p(q)+\ell_t(q)]+o(m)\}.                    (5)
```

One must use a uniform ell bound on the declared class; merely evaluating
a plateaued member is not such an upper bound. Formula (5) isolates an
explicit operational use and its exact counting normalization.

At q=1, however, c_p(1)=p log2: the bound has no savings at all for dense
spectra. The uniform approximate anti-sparsity condition from
`continued_feedback_selector_antisparsity_2026_09_06.md` is compatible
with magnitude profiles converging to the single atom at one. Their
permanent pressure tends, at fixed t, to

```math
\log[(1+e^{-4t})/2].                                     (6)
```

Thus the RELAXATION which knows only support-count bounds (4), the second
moment, and those uniform anti-sparsity inequalities still permits
`exp{mp log2+o(m)}` rows all in this nearly flat class. Its allowed tilted
exponent is

```math
p\log2+\tfrac12\log[(1+e^{-4t})/2]+t(1-\sqrt p)
\ge(p-\tfrac12)\log2>0.                                 (7)
```

Therefore no manipulation of these constraints alone can certify the
desired negative one-row exponent. This is an explicit feasible point
of the relaxed scalar profile/count problem; it is NOT a construction
of that many actual spin rows or a falsification of the full criterion.

Exact flatness might be forbidden by integrality or a high-gcd selector
filter, but this does not repair the relaxation: take two adjacent
admissible integer Fourier magnitudes surrounding sqrt(k), with the
proper even parity. Their normalized magnitudes both tend to one, their
gcd can remain only two, and frequencies can all be nonzero. Choosing
the mixture weights to make the second moment tend to one produces
(6) unchanged. These formal profiles satisfy the approximate geometric
constraints and do not assert that the inverse Walsh vector is ternary.
That last feasibility condition is precisely the missing information.

An upper theorem must therefore control counts of near-flat or otherwise
low-entropy DENSE spectra, or provide a stronger averaged inequality.
The selector filters and sparse-span counts rigorously remove important
exceptions, but they are not substitutes for that dense-profile input.
