# A direct cap comparison from diagonal contraction and variance budget

Date: 2026-09-07. Status: **proved elementary lemma; pending independent audit**.
This is a quantitative operation on weighted signings, not an endpoint
flatification theorem. In particular it pays a leading error when the excess
variance is of order N squared.

For a nonnegative hollow symmetric profile v on at most N vertices, let

```math
q(v)=\min_A\max_{x\in\{\pm1\}^d}
       \left|\sum_{i<j}A_{ij}\sqrt{v_{ij}}x_ix_j\right|.
```

## 1. Variance-budget comparison

Suppose 0<=s_i<=1 and w_ij>=s_i^2 s_j^2 v_ij. Set

```math
T=\sum_{i<j}(w_{ij}-s_i^2s_j^2v_{ij}),\qquad
K=\max_{i<j}\sqrt{w_{ij}},\qquad L=(N+2)\log2.
```

Then

```math
q(w)\le q(v)+\sqrt{2TL}+\frac{4K L}{3}.                 (1)
```

**Proof.** Select an optimizing signing A for v. Independently choose each
new sign C_ij so that

```math
E[\sqrt{w_{ij}}C_{ij}]=s_i s_j\sqrt{v_{ij}}A_{ij}.
```

This is feasible, including zero entries. The centered edge variable has
variance w_ij-s_i^2s_j^2v_ij and absolute value at most 2K. For any fixed spin
vector, Bernstein's exponential inequality gives a two-sided tail at most
2 exp(-L) at t=sqrt(2TL)+4KL/3. A union bound over at most 2^N spin vectors
has failure probability at most 1/2. Some actual choice of the signs has
all its centered quadratic responses bounded by t.

The mean quadratic form is the original form at the point (s_i x_i)_i in
the cube [-1,1]^d. A hollow quadratic form is affine separately in each
coordinate, so its absolute value on the cube is bounded by its maximum
on the vertices. Its cap is consequently at most q(v). This proves (1).
The Bernstein inequality here follows directly by bounding centered moments
E|X|^j <= (2K)^(j-2) E X^2 in the exponential series; no asymptotic
concentration hypothesis is used.

Taking all s_i=1 gives a genuine rounding operation whenever the extra
variance T is small. Taking s_i<1 permits an anisotropic contraction before
rounding. The bound depends on the variance introduced, not on the number
of coordinates whose values change.

## 2. Insertion and deletion

If v is a principal restriction of w, every row of w has variance sum at
most N-1, and r vertices were removed, then

```math
0\le q(w)-q(v)\le\sqrt{2r(N-1)L}+\frac{4K L}{3}.        (2)
```

For the lower bound, average a full quadratic form over all spins of the
removed vertices. Its remaining principal form cannot have larger cap.
For the upper bound, keep an optimal remainder and independently assign
unbiased signs to all newly incident edges. The same proof as (1) uses
T<=r(N-1). There is no claim that independently inserting r times and
adding r one-vertex bounds is as sharp as this simultaneous insertion.

## 3. Application: hidden partition counts, without temperature smoothing

Use the exact row-regular posterior profiles v(k,r) from
`flatify_independent_2026_09_07_reveal_count_concentration.md`, with
delta N<=m,n and u=N-k unrevealed vertices. Their entry amplitudes are
bounded by K_delta for sufficiently large N.

If two feasible counts r,r' differ by D and every count between them has
remaining positive fraction in [delta/2,1-delta/2], delete the D known
vertices whose labels differ. On the common remainder, known-known
entries coincide. The relative changes of nonzero entries touching the
unknown set are at most C_delta D/u, by the displayed exact formulas for
known/unknown and unknown/unknown variances. If D/u is sufficiently small,
scale only unknown vertex amplitudes by sqrt(1-epsilon), where
epsilon=C_delta D/u<1. This gives coordinate domination in either direction.

The variance removed by this contraction is at most
2 epsilon u(N-1)=O_delta(ND). The total edge masses of the two remainders
differ by at most O_delta(D^2): each removed set has D rows of mass N-1,
and only its internally counted edge mass differs. Thus the T in (1) is
O_delta(ND). Adding back the D deleted vertices with (2) yields

```math
|q(v(k,r))-q(v(k,r'))|\le C_delta N(\sqrt D+1).           (3)
```

For larger central D/u, divide the interval into a bounded number
(depending only on delta) of shorter intervals on which the contraction
above has epsilon<1/2. The sum of their square roots is at most a
delta-dependent constant times sqrt D. This proves the same bound (3).

Let K_k be the actual hypergeometric count and r_0 the nearest feasible
integer to mk/N. For u>=sqrt N, use the central event
|K_k-r_0|<=delta u/4. Its complement has probability O_delta(1/u).
The exact variance p(1-p)ku/(N-1)<=u/2 gives
E sqrt(|K_k-r_0|)<=O(u^(1/4)). Equation (3) therefore contributes at most
O_delta(N^(5/4)). The complement contributes at most
O_delta(N^(3/2)/u), since the unbiased-sign construction in (1) gives
q(v)<=O_delta(N^(3/2)) for every row-regular profile.

For u<sqrt N, delete all unknown vertices. Each profile agrees there with
a terminal partition, and terminal optimized caps depend only on the two
terminal sizes, not the labels. Applying (2) to the two profiles and their
respective terminal completions gives O_delta(N sqrt u+N), again at most
O_delta(N^(5/4)). We have proved, uniformly in reveal time,

```math
\sup_k E|q(v(k,K_k))-q(v(k,r_0(k)))|
                 =O_delta(N^(5/4)).                       (4)
```

This independently recovers the pressure argument's strongest exponent
without a growing-temperature passage. It does **not** compare the uniform
initial profile to the split terminal profile: neither endpoint has a count
fluctuation in the first place.

## 4. Exact remaining gap

Applying (1) directly to favorable flatification still has T=Theta(N^2)
for the obvious diagonal contractions, hence an O(N^(3/2)) error. The new
operation becomes useful for the convergence target only if a construction
first makes the required variance addition o(N^2), or if another argument
controls its actual cap effect more sharply. No such construction is
proved here. The top-level cross-order obligation is unchanged.
