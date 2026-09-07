# Independent audit: bounded-cap Maiorana--McFarland permutation packing

Audited files:

- `extremal_information/drafts/bounded_cap_maiorana_permutation_packing.md`
- `extremal_information/experiments/verify_bounded_cap_maiorana_permutations.py`

Verdict: **PASS, with minor verifier/documentation repairs only.**  I found no
proof-critical error in Theorem BMP.1, Lemma BMP.2, the permutation
conditioning, the Turan extraction, or the response separation.  The claimed
packing rate is genuinely `exp(Omega(q log q))` for `n=q^2`, hence the stated
`Theta(sqrt(n) log n)` response-information rate for this explicit family.

## 1. Walsh and Maiorana--McFarland identities

With coordinates ordered as `(a,b),(u,v) in V^2`,

```math
(Ws_\pi)(a,b)
=\sum_{u,v}(-1)^{a\cdot u+b\cdot v+u\cdot\pi(v)}
=q(-1)^{b\cdot\pi^{-1}(a)}.
```

Thus `y_pi=q^{-1}Ws_pi` is Boolean and `Wy_pi=q s_pi`, because
`W^2=nI` and `n=q^2`.  For `b=s_id`, the displayed formula gives `Wb=qb`.
Consequently `H=D_bWD_b` is symmetric, `H^2=nI`, and `H 1=q 1`.
Moreover

```math
tr(H)=tr(W)
=\sum_{a,b}(-1)^{a\cdot a+b\cdot b}=0
```

for `m>=1`.  Hence subtracting the diagonal produces a hollow signing `A`
while, on Boolean vectors,

```math
H_A(x)=\tfrac12x^TAx=\tfrac12x^THx.
```

The operator bound `|x^THx| <= q n` and the Boolean eigenvector `1` prove
`Q(A)=qn/2=n^{3/2}/2`.  Diagonal switching preserves this exactly.  All
normalizations in BMP.3--BMP.6 are correct.

## 2. Pair Rayleigh identity

For `w=s_pi odot s_sigma` and
`tau(v)=v+pi(v)+sigma(v)`, one has

```math
(b\odot w)(u,v)=(-1)^{u\cdot\tau(v)}.
```

Writing `g=b odot w` and summing the `u` variable in `g^TWg` gives

```math
g^TWg
=q\sum_{x,y}(-1)^{x\cdot y+\tau(x)\cdot\tau(y)}.
```

Since `w^THw=g^TWg` and `qn=q^3`, division by `qn` gives exactly the
average over `q^2` pairs in BMP.9.  There is no missing factor of `q` or
`2`.

## 3. Random-function Fourier-rank tail

The Fourier identity

```math
rho(tau)=q^{-1}\sum_x \widehat f_{\tau(x)}(x)
```

is exact.  If the average is greater than `1/4`, more than `q/7` summands
are greater than `1/8`.  Parseval permits fewer than `64` such frequencies
for each fixed output value `a`, so `|tau(G)|>q/448`.  A set this large
cannot lie in a subspace of dimension at most `m-9`, whose size is at most
`q/512`; hence its span contains `r=m-8` independent values.

For fixed independent `a_1,...,a_r` and fixed witnesses `x_i`, the map

```math
z -> (a_1\cdot z,...,a_r\cdot z)
```

is surjective.  Therefore, for each `y`, its value at the uniform random
`tau(y)` is uniform on `F_2^r`; these vectors are independent over `y`.
Equivalently, all `rq` signs entering the `r` Fourier coefficients are
jointly independent Rademachers (the deterministic `x_i dot y` signs do
not change this).  The coefficients are therefore independent, not just
pairwise independent.  Hoeffding gives `exp(-q/128)` per coefficient and
`exp(-rq/128)` jointly.  The union bound over at most `q^r q^r` fixed
choices proves BMP.10.  Selecting the witnesses adaptively causes no gap,
because the final union is over every deterministic choice; the extra
conditions `tau(x_i)=a_i` are safely discarded in the upper bound.

## 4. Conditioning and Turan extraction

For two independent uniform functions `pi,sigma`, the values
`v+pi(v)+sigma(v)` are independent uniform values, so BMP.10 applies.
Conditioning both functions to be permutations multiplies an event
probability by at most

```math
(q^q/q!)^2 <= e^{2q}.
```

This proves BMP.14.  The bad-pair relation is symmetric.  Every diagonal
pair is bad (`tau(v)=v` and `rho=1`), so if `N=q!` and `e` is the number of
unordered off-diagonal bad pairs, the ordered bad-pair density is exactly

```math
(N+2e)/N^2.
```

Thus Caro--Wei/Turan in the form
`alpha >= N^2/(N+2e)` gives the reciprocal-density lower bound claimed in
BMP.16.  Since `r=m-8=(log q)/(log 2)-8`,

```math
rq/128-2r log q-2q = Omega(q log q).
```

The constants only become positive at a very large absolute threshold;
this is compatible with the theorem's unspecified `q_0`.

## 5. Response bound and projective factor

At its own query, changing variables by `D_{s_pi}` gives the field `q 1`;
the spectral upper bound is attained at `u=1`, so the response is
`3qn/2`.  At child `sigma` and query `y_pi`, the field is `qw`.  With
`K=2qI-H`,

```math
K^{-1}=(2qI+H)/(3q^2),
```

and completing the square on `||u||_2^2=n` gives

```math
R_sigma(y_pi)
 <= qn + \tfrac12(qw)^T K^{-1}(qw)
 = qn\left(1+\frac{2+rho(pi,sigma)}6\right).
```

For `rho<=1/4` this is `11qn/8`, so the signed response difference is at
least `qn/8` at `y_pi` and at most `-qn/8` at `y_sigma`.  Its oscillation
is at least `qn/4`; the definition `d_proj=osc/2` therefore yields exactly
`qn/8=n^{3/2}/8`.  The factor is correct.

If two children shared a summary whose decoded responses both had uniform
(or projective) error strictly below `qn/16`, the triangle inequality
would put their projective distance strictly below `qn/8`.  Hence all code
members require distinct summary states.  The lower bound is
`log |C|=Omega(q log q)` bits.  Conversely, explicitly listing a
permutation costs `log(q!)=Theta(q log q)` bits and reconstructs its child
and response, so the family-level `Theta(q log q)` description claim is
valid.  This does not assert the code is efficiently constructible.

## 6. Verifier audit

Executed:

```text
./.venv/bin/python \
  extremal_information/experiments/verify_bounded_cap_maiorana_permutations.py
```

Result:

```text
q=2: maximum code 2, 4 Rayleigh checks, 4 response checks
q=4: maximum code 20, 576 Rayleigh checks, 36 response checks
first positive crude margin m=265
total exact identities=620
```

The script correctly checks all ordered Rayleigh identities at `q=2,4`,
computes the exact maximum compatibility clique at `q=4`, and exhaustively
enumerates Boolean responses for a six-member audited subcode.

### Minor repairs recommended

1. Section 6 says the verifier uses *integer arithmetic*, but two
   assertions convert exact integer quantities to floats (`rho` and the
   final `/2`).  Either change that wording to "exact finite enumeration"
   or cross-multiply the inequalities.  For example, BMP.18 can be checked
   exactly as

   ```python
   assert 6 * q * value <= 8 * q * q * n + n * bias_sum
   ```

   and the projective gap as

   ```python
   assert 4 * (at_i - at_j) >= q * n
   ```

   This is not a theorem defect; the present numbers are far below any
   floating-point precision issue.

2. Add cheap structural regression assertions for `H == H.T`,
   `H @ H == n I`, `H 1 == q 1`, `trace(H)==0`, and `y_id==b`.  The proof
   establishes them, and the exhaustive energy check indirectly catches
   the key consequence, but the current verifier does not explicitly test
   each displayed structural identity.

3. Keep the wording that only six representatives receive exhaustive
   Boolean-response checks.  The maximum code has twenty members; the
   all-pair response theorem for those twenty follows from the proved
   square-completion inequality plus the fully checked pair biases, not
   from 400 exhaustive response enumerations.

These are auditability improvements only.  The rigorous status of BMP.1
is justified as written.
