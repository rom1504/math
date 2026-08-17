# Independent audit: Boolean-port Fourier feature algebra

**Verdict: PASS.**  PF.1, the multiplier formula (PF.7), the exact quotient
count, and the `p<=3`/`p=4` boundary are correct at the stated **labelled
response-table** level.  I found no mathematical error.  I replaced one
floating-determinant verifier check by exact character diagonalization; the
strengthened verifier passes 256 exact checks.

## 1. Normalization and multiplier audit

The projective group has order

```math
K=|G_p|=2^{p-1}.
```

Because both `g_p(z)=|sum_i z_i|` and every even character are invariant
under `z -> -z`, normalized Fourier expectation on `G_p` equals normalized
expectation on the full Boolean cube.  Thus the multiplier attached to an
even set `S`, `|S|=2k`, really is

```math
\widehat g_p(S)=2^{-p}\sum_z
 \left|\sum_i z_i\right|\prod_{i\in S}z_i.          \tag{PA.1}
```

For odd `p`, the draft's passage to `P=p+1` is exact: the old sum `a` is
odd, and

```math
{|a+1|+|a-1|\over2}=|a|.                           \tag{PA.2}
```

For even `P`, independent binomial evaluation through `P=12` and exact
character diagonalization through `p=7` give

```math
\widehat g_p(\varnothing)={(P-1)!!\over(P-2)!!},
```

and

```math
\widehat g_p(S)=(-1)^{k-1}
{(2k-3)!!(P-2k-1)!!\over(P-2)!!}.                  \tag{PA.3}
```

The endpoint convention is sound: when `2k=P`, the second numerator factor
is `(-1)!!=1`.  Every factor in absolute value is positive, so no Fourier
channel vanishes.

The convolution normalization in PF.10 is also correct.  PF.4 uses an
unnormalized group sum, whereas normalized convolution contains a factor
`1/K`; hence

```math
\widehat L(S)=K\widehat\mu(S)\widehat g_p(S).       \tag{PA.4}
```

## 2. Invertibility and minimal quotient

Let `R` be the `K by K` response matrix

```math
R_{s,\epsilon}=g_p(s\epsilon).
```

If `C` is the projective character table, exact integer arithmetic gives

```math
CRC^T=K^2\operatorname {diag}
 \bigl(\widehat g_p(S):S\text{ even}\bigr),
\qquad CC^T=KI.                                    \tag{PA.5}
```

Formula (PA.3) therefore makes `R` invertible.  Equality of all labelled
responses is equivalent to equality of every projective row count.  All
weak compositions of `n` into `K` parts are realizable by choosing the
corresponding rows, so the exact number of response classes is indeed

```math
{n+K-1\choose K-1}.                                \tag{PA.6}
```

There is no hidden realizability or integrality restriction.  Equality of
histograms is exactly equality up to row permutation and independent row
signs.  The draft assumes the same row count; this is harmless and slightly
redundant, since the nonzero constant Fourier channel already recovers total
mass.

The minimality claim is correctly scoped.  It is a cardinality/minimal
contextual-quotient statement for the entire labelled endpoint table.  It
does not mean that a particular coordinate system literally has to store
each histogram bin, nor that the scalar maximum over endpoint words has the
same number of equivalence classes.

## 3. The first non-Gram channel

The Gram matrix contains the constant channel and all degree-two even
Fourier coefficients of the projective histogram.  These exhaust the even
characters for `p=1,2,3`, so for `p<=3` equal Gram matrices imply equal
**labelled response tables**, and therefore equal joint maxima as a
corollary.

At `p=4` there is one additional character, the degree-four parity.  The
uniform/even-parity example changes exactly this invisible coordinate and
has joint supports `3n/2` and `2n`.  It therefore proves failure not only of
labelled-table recovery but also of the coarser endpoint-maximized response.
This establishes four as the first possible port count for the pure-linear
equal-Gram collision.

The distinction from the interacting trust problem is correctly maintained
in the draft.  PF.1 says nothing by itself about how a quadratic child
correlates with row types.  The separate common-Hadamard theorem EG.2 supplies
that strengthening.

My independently derived theorem TC.3 in
`three_port_gram_closure.md` corroborates EG.2 with the same order-16 seed
and an independent tensor argument: two four-port tuples of Boolean top
eigenvectors have equal `(G,R)`, supports `2n` and `7n/4`, and at width
`m=r` their Boolean trust responses differ by at least `rn/8`.  Arbitrary
public exact-sign completion costs only `O(n)`.  Thus the degree-four
feature exposed abstractly by PF.1 remains leading-order observable under a
common dense regular-Hadamard quadratic child.

## 4. What is and is not new

Fourier invertibility of a finite group convolution is classical.  The
useful theory-level conclusion here is its exact identification of the
minimal composable feature algebra for Boolean port contexts, together with
the sharp arity boundary and the interacting EG.2/TC.3 realization.  It
does not solve approximate compression: small Fourier multipliers may make
some channels expensive to decode, but uniform response distortion depends
on the admissible nonnegative histograms and cannot be read off from
invertibility alone.

## 5. Verifier repair

The original verifier used a rounded floating determinant for response
matrices.  That is numerically fragile and is not an exact certificate of
invertibility.  It now constructs the full projective character table and
checks (PA.5) in integer arithmetic, including the predicted rational
multiplier on every diagonal entry.  Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_boolean_port_fourier_feature_algebra.py
```

to reproduce the audit.

## 6. Addendum: audit of the sparse coreset theorem PF.2

**Verdict: PASS.**  The Hoeffding constant in PF.16 is sufficient (in fact
it has one extra `log 2` of slack), the projective-Hamming metric calculation
is exact, and the packing exponent is correctly scoped to fixed distortion.

For one endpoint query the two-sided Hoeffding failure probability is at
most `2 exp(-2K eta^2)`.  There are `2^(p-1)` distinct projective queries,
so the complete union bound is

```math
2\,2^{p-1}e^{-2K\eta^2}=2^p e^{-2K\eta^2}.          \tag{PA.7}
```

PF.16 chooses

```math
K\ge{(p+1)\log2\over2\eta^2},
```

which makes (PA.7) at most `1/2`, hence gives strictly positive success
probability.  The nearly sharp union-bound threshold would be just above
`p log 2/(2 eta^2)`; retaining the displayed slack avoids an equality
edge case.  Sampling with repetitions produces an empirical probability
measure supported on at most `K` types, and ordered-sample encoding gives
exactly the claimed upper count `|G_p|^K`.

For point masses, choose representatives with projective Hamming distance
`h(s,t)<=p/2`.  Reverse triangle inequality gives the upper bound
`2h/p`, while the query `epsilon=s` has responses `1` and
`1-2h/p`.  Therefore

```math
d_p(\delta_s,\delta_t)={2h(s,t)\over p}             \tag{PA.8}
```

exactly.  Greedy deletion of projective balls of radius `theta p` has ball
volume `2^{(H_2(theta)+o(1))p}` and ambient size `2^{p-1}`, proving the
stated packing exponent.  More precisely, this packing forces the same
linear-bit lower bound for covering radius strictly below `theta` (or with
an arbitrarily small constant slack); pairwise distance equal to `2theta`
alone should not be read as a closed-ball lower bound exactly at the
boundary radius `theta`.

The static/dynamic caveat in PF.2 is essential and correctly stated.
Convex combination of existing approximants preserves their worst error,
whereas independent resparsification at successive levels can add error;
PF.2 by itself is not an all-depth mergeable summary theorem.
