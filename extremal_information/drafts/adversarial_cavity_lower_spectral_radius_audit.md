# Independent audit: adversarial cavity and lower spectral radius

Audited file:
[`adversarial_cavity_lower_spectral_radius.md`](adversarial_cavity_lower_spectral_radius.md)

Verdict: **ACCEPT with minor clarification.**  The partition-function/norm
comparison, lower-spectral-radius limit, zero-temperature convergence,
contractive projective carrier, additive eigenfunction, finite-net
cycle-mean estimate, diagonal nonperiodic example, and dense-bridge lower
bound are mathematically correct.  The cited theorems match their uses after
the already-applied correction removing an inaccurate numbered-example
attribution.  The only material wording caution is that (AC.20)--(AC.21)
compress the **asymptotic mean pressure**; they do not by themselves give a
depth-independent additive approximation to every unnormalized
finite-horizon response.

## 1. Transfer product and lower spectral radius

With the transpose convention in (AC.4),

```math
Z_{n,\beta}(w)=\mathbf1^T
T_{d_n}\cdots T_{d_1}\mathbf1.
```

For a nonnegative `q by q` matrix `P`, the total entry sum lies between its
largest column sum and `q` times that sum, proving (AC.8).  Taking minima
over the identical word set preserves the factor-`q` comparison.

If

```math
a_n=\min_{|w|=n}\|P_w\|_1,
```

concatenating a minimizing length-`m` product with a minimizing length-`n`
product gives `a_(m+n)<=a_m a_n`.  Hence `log a_n` is subadditive and
Fekete proves existence of `lim a_n^(1/n)`.  The factor `q` in (AC.8)
vanishes after division by `n`, so (AC.5) follows with no convexification of
the disorder alphabet.

The linked Guglielmi--Zennaro review does support the norm-independent lower
spectral-radius definition and limit in its Section 2.  Its attribution is
slightly more nuanced than “the Gurvits limit”: it follows Gurvits for the
definition and spectral-radius formula, while citing later sources for some
equivalent limit statements.  Since the norm-limit used here is proved
directly by Fekete, this is at most an attributional wording issue, not a
dependency.

## 2. Zero temperature

There are exactly `q^(n+1)` spin paths.  Therefore, for each fixed disorder
word,

```math
\max H_w\le\beta^{-1}\log Z_{n,\beta}(w)
\le\max H_w+(n+1)\beta^{-1}\log q.
```

Taking minima over the same words proves (AC.11).  Since
`F_(n,beta)/(beta n)` converges, these inequalities imply

```math
\limsup_n G_n/n-\liminf_nG_n/n\le\log(q)/\beta
```

for every fixed `beta`.  The sequence is bounded by the local reward bound;
sending `beta` to infinity forces its limsup and liminf to agree.  Taking
the resulting limit in (AC.11) then gives the full sandwich in (AC.7) for
every `beta`, and its zero-temperature limit.  The order of limits is valid
and no uniform-in-`beta` transfer theorem is assumed.

## 3. Compact projective carrier

Under (AC.14), one may define explicitly

```math
X=\left\{p\in\Delta_{q-1}:
 {a\over qb}\le p_i\le {b\over qa}\text{ for every }i
 \right\}.
```

Every `tau_d(p)` lies in this compact convex subset of the simplex interior,
even when the input is on the boundary, so `X` is common and invariant.
This explicit definition would improve the draft's phrase “a compact
interior simplex.”

For a positive matrix with entries in `[a,b]`, its projective diameter is at
most

```math
2\log(b/a),
```

because every cross ratio is at most `(b/a)^2`.  Birkhoff's coefficient
`tanh(Delta/4)` is therefore at most

```math
\tanh\left({1\over2}\log(b/a)\right)={b-a\over b+a},
```

which verifies (AC.15).  The Birkhoff and Carroll citations state precisely
the projective contraction formula being used.

For probability vectors `p,p'`, let

```math
m=\min_i p_i/p_i',
\qquad M=\max_i p_i/p_i'.
```

Normalization gives `m<=1<=M`, and every positive linear functional ratio
lies in `[m,M]`.  Hence its absolute log is at most
`log(M/m)=d_H(p,p')`, proving (AC.16).

## 4. Additive eigenfunction and identification of the eigenvalue

Let `L=1/(1-kappa)` and fix `p_* in X`.  The family

```math
\{f\in C(X): f(p_*)=0,\ \operatorname{Lip}_{d_H}(f)\le L\}
```

is convex, closed, uniformly bounded, and equicontinuous, hence compact in
the sup norm by Arzela--Ascoli.  The normalized map

```math
f\longmapsto\mathcal Vf-(\mathcal Vf)(p_*)
```

is continuous and maps this set into itself because
`1+kappa L=L`.  Schauder gives `V u=u+lambda`.  Additive homogeneity and
sup-norm nonexpansiveness then give (AC.19).

Telescoping (AC.12) shows

```math
\mathcal V^n0(p)=
\log\min_{|w|=n}\|P_wp\|_1.
```

For fixed interior `p`,

```math
(\min_i p_i)\|P\|_1\le\|Pp\|_1\le\|P\|_1
```

for every nonnegative `P`.  The factor is independent of product length.
Thus the linear growth rate in (AC.19) is exactly
`log check-rho`, as claimed.

## 5. Finite net and cycle mean

On the compact interior set `X`, Hilbert and Euclidean metrics are uniformly
bi-Lipschitz.  Its dimension is `q-1`, so a Hilbert `delta`-net has
`O_(q,a,b)(delta^(-(q-1)))` points.

For each net point `p` and symbol `d`, direct an edge to a nearest net point
`p'` to `tau_d(p)` and give it reward `r_d(p)`.  Since `u` is
`L`-Lipschitz, the additive eigen-equation implies for every such edge

```math
r_d(p)+u(p')\ge\lambda+u(p)-L\delta.                    \tag{A.1}
```

Summing around any directed cycle gives mean at least
`lambda-L delta`.  Conversely, at each net point choose a symbol attaining
the minimum in the eigen-equation.  Following these selected edges in the
finite graph eventually reaches a cycle, and along each selected edge the
reverse estimate holds:

```math
r_d(p)+u(p')\le\lambda+u(p)+L\delta.                    \tag{A.2}
```

That cycle has mean at most `lambda+L delta`.  This proves (AC.21) with the
displayed constant and does not accumulate an extra reward-rounding error.

What is compressed here is the minimum asymptotic cycle mean.  A fresh
rounding error at every step may create an `O(n delta/(1-kappa))` total
finite-horizon error, although its error per step is uniformly controlled.
Accordingly, “asymptotic response-rate conclusion” or “mean-pressure
carrier” would be more precise than an unqualified depth-independent
response carrier.

## 6. Nonperiodic diagonal example

The matrices in (AC.22) commute.  A word with `k` copies of `A` and `n-k`
copies of `B` is

```math
\operatorname{diag}
\left(3^{-k}2^{n-k},\ 3^k2^{-(n-k)}\right),
```

whose logarithmic spectral radius divided by `n` is exactly (AC.23).
Rational frequencies can approach `log2/log6`, so the lower spectral radius
is one.  Equality for a finite word would require
`3^k=2^(n-k)`, impossible for positive length by unique factorization.

The revised source wording is accurate: this elementary example is merely
analogous to the lower-finiteness failures studied by Bochi--Morris, not a
numbered example imported from that paper.  Bousch--Mairesse does concern
aperiodic matrix-product optimization, but no result from it is used in the
proof.  The draft correctly notes that zero off-diagonal entries make this a
hard-constraint limit outside the uniformly positive finite-temperature
hypothesis.

## 7. Dense-bridge falsifier

For a column sum `S_j=(B^Tx)_j` under uniform Boolean `x`,

```math
\mathbb ES_j^2=n,
\qquad
\mathbb ES_j^4=3n^2-2n\le3n^2.
```

Log-convexity of `L^p` norms gives

```math
\|S_j\|_2\le\|S_j\|_1^{1/3}\|S_j\|_4^{2/3},
```

and therefore `E|S_j|>=sqrt(n/3)`.  Summing expectations over columns gives
an `x` with `||B^Tx||_1>=n^(3/2)/sqrt(3)`; optimizing `y` coordinatewise
proves (AC.25).  Column dependence is irrelevant.

This rigorously shows that ordinary balanced-cut transfer pays a leading
bridge term.  The subsequent claims about separator width `n`, transfer
dimension `2^n`, and projective contraction tending to one correctly
describe the standard local transfer representation, but should not be read
as an information-theoretic no-go for every possible nonlocal quotient.  The
draft already makes that limitation explicit.

## 8. Recommended repairs

Before canonical promotion, I recommend only:

1. define one explicit invariant compact set `X` as above;
2. insert the two cycle inequalities (A.1)--(A.2), which make (AC.21)
   immediate and distinguish mean-pressure error from total-horizon error;
3. qualify “response-rate conclusion” as an asymptotic mean-pressure result;
4. optionally replace “standard Gurvits ... limit” by “the standard lower
   spectral-radius limit” while retaining the review citation.

No theorem statement, constant, or portfolio judgment needs to change.
