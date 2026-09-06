# Independent audit of the stable flat-profile count

Date: 2026-09-06. Full read and reconstruction of
`continued_feedback_stable_flat_profile_count_2026_09_06.md`. The theorem
and its quantified limitation pass. It is a genuine stable dense-profile
bound, not a consequence of exact Walsh divisibility.

## 1. Recovery from exact pair sums

For a direction a!=0 the unordered input-pair types are 00, ++, --,
+-, 0+, 0-. Given the pair sums and their six counts, zero sums require
choosing which are 00 and orienting the opposite-sign pairs; sums plus
or minus one each require one orientation bit. Thus the exact recovery
count is bounded by

```math
2^{n_{0+}+n_{0-}+n_{+-}}
\binom{n_{00}+n_{+-}}{n_{00}}.
```

Its logarithmic entropy bound is concave and homogeneous in the counts.
Every unordered pair of distinct inputs occurs for exactly one direction.
Therefore the average pair count is N_rN_s/(m-1) for distinct symbols,
and N_r(N_r-1)/(2(m-1)) for equal symbols. Jensen supplies a direction
whose recovery cost is at most the function of these average counts.

For fixed zero and total nonzero counts, the dependence on n_(+-)
has derivative log2+log[(n00+n+-)/n+-]>0. Maximizing N+N- at a balanced
global sign split is consequently a valid upper bound even for an
unbalanced actual row. This recovers precisely the displayed B(p).
Direction and six histogram counts cost only O(log m) extra bits.

## 2. Approximate spectral signs and the repair factor

On U=a-perp, of size M=m/2, the restricted Walsh transform is exactly
the M-point transform of the pair sums, with inverse factor 1/M.
Encoding the signs of these M coefficients costs M bits. Replacing
their magnitudes by sqrt(k) yields squared pair-sum error at most

```math
M^{-1}(\epsilon^2 k m)=2\epsilon^2 k.
```

Nearest rounding to the five pair-sum symbols makes at most
8 epsilon^2 k errors, since every wrong rounding costs at least 1/4
in squared distance. Dividing by M gives e=16p epsilon^2. The code
for changed locations and their correct values has size at most
sum_(j<=eM) binomial(M,j)5^j. Its entropy exponent is h(e)+e log5
when e<5/6. Using five instead of four replacement symbols is a valid
deliberately loose upper count.

The resulting description determines the input row injectively. Dividing
the count of weight-k ternary rows by binomial(m,k) counts exactly the
expected number of supported spin rows. Conditioning the selector on a
probability-tending-to-one event can only increase this upper bound by
the reciprocal event probability, with no exponential change.

For independent X1,X2 with probabilities (p/2,1-p,p/2), conditioning
their ordered pair on its sum gives conditional entropy exactly 2B(p).
Since H(X)=h(p)+p log2, the two stated formulas for c_flat agree:

```math
c_{\rm flat}=\tfrac12\log2+B(p)-h(p)
=(p+\tfrac12)\log2-\tfrac12H(X_1+X_2).
```

## 3. Permanent continuity is uniform over deletions

Let a=|h|/sqrt(k), so ||a||2=sqrt(m) and ||a-1||2<=epsilon sqrt(m).
Cauchy--Schwarz gives sum|a^2-1|<=2 epsilon m. On any retained subset
and any permutation pi,

```math
\sum|a_u a_{\pi(u)}-1|
\le\|a-1\|_2\|a\|_2+\sqrt m\|a-1\|_2
\le2\epsilon m.
```

The quadratic part of log K therefore changes a permutation product's
logarithm by at most 4t epsilon m. The 1-Lipschitz log-cosh part costs
another 4t epsilon m. After averaging permutation products and taking
the square root the error is 4t epsilon per coordinate, as stated.
This uses no coordinatewise closeness or bound on the largest amplitude.
The maximum over the removed diagonal coordinate is harmless because
all these inequalities are uniform in that deletion.

## 4. What the numerical margin does and does not show

At p=15/16, c_flat=0.3168851694. For epsilon=0 the derivative of the
tilted upper expression is

```math
1-\sqrt p-\frac2{e^{4t}+1}.
```

Its unique zero is t=(1/4)log[2/(1-sqrt(p))-1]=1.0317185344,
and the corresponding minimum is still +0.0110751745. At t=6 the
upper expression is about +0.16084. Positive epsilon only increases
the displayed bound.

Therefore this theorem does not eliminate the nearly flat class at a
tilt useful for the whole row sum, despite its substantial improvement
over the trivial p log2 row count. It neither assumes nor proves that
approximately flat rows lie near exact bent functions. A stronger
recovery argument or additional actual profile counting is still needed.
