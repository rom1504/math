# Exact quadratic energy variance on full rank-one dictionaries

2026-09-17. This is a sharp finite scope theorem for the independent
rank-one dictionary used by the Bernoulli response construction. Its
whole-dictionary response barrier is an overcover barrier, not an
automatic obstruction on an actual positive-cap nearcode.

## 1. Exact character expansion and sharp maximum

Partition N coordinates into blocks of sizes n_a=k_a p_a. In block a
sample independently uniform signs c^(a) in {+-1}^k_a and v^(a) in
{+-1}^p_a, and set X^(a)_ij=c^(a)_i v^(a)_j. All blocks are independent.
For any hollow symmetric full signing A, write H_A(X)=sum_(u<v) A_uv X_u X_v.
Then

```
E H_A(X)=0,
E H_A(X)^2 <= V(k,p)
 := (N^2+3 sum_a n_a^2)/2
       -(3/2) sum_a n_a(k_a+p_a)+N.                 (1)
```

Equality is attained by the genuine full signing A=J-I. Thus (1) is
the exact maximum variance over every possible full signing, not just
an operator-norm bound.

Proof: distinct characters in the independent underlying c,v signs
are orthogonal. Within a block, the character classes of physical
edges are:

- Edges between two columns in the same row: each column pair has
  one character v_j v_l, with k contributing edges.
- Edges between two rows in the same column: each row pair has one
  character c_i c_l, with p contributing edges.
- Edges changing both row and column: each rectangle has one
  character c_i c_l v_j v_t, with its two opposite diagonals.

Consequently one block contributes at most

```
k^2 binom(p,2)+p^2 binom(k,2)+4 binom(k,2)binom(p,2)
 =2n^2-(3/2)n(k+p)+n.                              (2)
```

Every cross-block physical edge has its own distinct character, so
its contribution to the variance is exactly one; there are
sum_(a<b) n_a n_b such edges. Adding (2) gives (1). At A=J-I every
coefficient sum in each class has its largest possible absolute
value, proving sharpness. Empty character classes when k=1 or p=1
cause no difficulty and the same formulas remain exact.

For one block, (1) reads `2N^2-(3/2)N(k+p)+N`. For any number and
sizes of blocks, V<=2N^2. The previously recorded quadratic-feature
covariance norm is at most max_a{k_a,p_a,2}; (1) is sharper for full
sign coefficients because it uses every individual character-class
size rather than only their largest size.

## 2. Consequence for actual high-energy subsets

Let F carry the uniform independent rank-one law above. For every
full signing A and threshold u>0,

```
P_F{|H_A(X)|>=u} <= V/u^2.                          (3)
```

In particular, if Q(A)>=c N^(3/2) and T<=eta N^(3/2), eta<c,

```
P_F{X in E_A(T)} <= 2/[(c-eta)^2 N].                (4)
```

Thus the WHOLE dictionary cannot be a nearcode at a fixed positive
normalized cap. Any theorem forcing its minimax response beta(lambda)
is a theorem about that sufficient overcover and may protect many
low-energy words which the actual parent maximum need not use.

More generally, if a law pi on F has density w relative to this
uniform law and is entirely supported on |H_A|>=u, Cauchy--Schwarz
gives

```
E_F w^2 >= u^2/V.                                 (5)
```

Indeed u<=E_pi|H_A|<=sqrt(E_F w^2 * E_F H_A^2).
A genuine positive-cap nearcode law therefore requires density
L2 norm of order at least sqrt(N), even if it has many support words.

Neither (3) nor (5) implies that conditioning on high energy lowers
the minimax absolute-overlap response. A subset with only O(1/N) of
the dictionary can still be exponentially large in k+p. Conversely,
the whole-dictionary minimax theorem does not establish a matching
lower bound on that conditioned subset. A relation between this
necessary energy concentration and favorable scalar column response
is a distinct missing value theorem; no such implication is asserted.

The response theorem and its actual sign realization are in
`artifacts/paper_bernoulli_rank_one_dictionary_2026_09_17.md`.
The fixed-degree energy-tilt failure on the full physical cube is in
`artifacts/paper_bernoulli_fixed_power_tilts_2026_09_17.md`; it concerns
column laws, whereas (3)--(5) here concern high-energy query subsets.

The exact finite character and sharpness replay is
`computations/paper_discrepancy_2026_09_17_rank_one_energy_scope.py`:
155 full signings on five single/multiple-block configurations, with
all underlying cube words enumerated and all variances checked as integers.
