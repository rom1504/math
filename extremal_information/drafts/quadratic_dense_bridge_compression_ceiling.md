# A strict compression bound, and a ceiling for coefficientwise methods

Status: rigorous bridge-independent upper bound.  This note treats only the
compression side for genuine bounded/sign quadratic Boolean children.  It
does not import the arbitrary-landscape packing theorem.

## 1. Setup

Put

```math
m={n\choose 2},\qquad
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,qquad x\in\{-1,1\}^n,
```

and, for an arbitrary fixed bridge `B`,

```math
F_A(y):=(P_BH_A)(y)=\max_x\{H_A(x)+x^TBy\}.
```

The basic contraction is

```math
\|F_A-F_{A'}\|_\infty
\le \|H_A-H_{A'}\|_\infty
\le \sum_{i<j}|a_{ij}-a'_{ij}|.                 \tag{QC.1}
```

The first inequality follows by comparing the same maximands at each `y`.
It is valid for every bridge; density, randomness, and an operator-norm
bound on `B` are irrelevant to the upper bounds below.

## 2. Sign quadratics: a constant-fraction sparse code

The strongest summary bound in this note allows the decoder to use a sparse
weighted quadratic as its surrogate.  The input child remains a genuine sign
quadratic.

### Theorem QC.1 (universal sparse-mask compression)

Let `0<epsilon<=1` and `n>=64/epsilon^2`.  There is a family of at most `m+2`
masks, depending only on `(n,epsilon)`, such that every sign matrix `A` can be
encoded by one mask and the signs of `A` on that mask.  The decoder constructs
a quadratic `H_Ahat` with coefficients in `{-1/q,0,1/q}`, where

```math
p={\epsilon^2\over2},\qquad q=1-p\ge {1\over2},
```

and obtains

```math
\|H_A-H_{\widehat A}\|_\infty\le\epsilon n^{3/2},
\qquad
\|P_BH_A-P_BH_{\widehat A}\|_\infty\le\epsilon n^{3/2}.    \tag{QC.S1}
```

The number of stored bits is at most

```math
b_{\rm sparse}(n,\epsilon)
\le
\left\lceil(1-\epsilon^2/4)m\right\rceil
+\left\lceil\log_2(m+2)\right\rceil.                       \tag{QC.S2}
```

In particular, for every fixed positive `epsilon`, this is a strict
constant-fraction saving from the exact `m` sign bits, while remaining
`Theta(n^2)`.

#### Proof and arithmetic audit

Choose a random mask by taking independent `Z_e~Bernoulli(q)` on the `m`
edges and set

```math
\widehat a_e={a_eZ_e\over q}.
```

For fixed `A` and `x`, the error summands

```math
X_e=a_ex_ix_j(1-Z_e/q)
```

are independent and centered.  Since `p<=q`, they satisfy `|X_e|<=1`, and

```math
\sum_e E X_e^2=m{p\over q}
=m{\epsilon^2\over2-\epsilon^2}
\le m\epsilon^2\le {\epsilon^2n^2\over2}.                 \tag{QC.S3}
```

Bernstein's inequality, with `E=epsilon n^(3/2)`, gives

```math
Pr\{|H_A(x)-H_{\widehat A}(x)|>E\}
\le2\exp\left\{-{E^2\over
2(\epsilon^2n^2/2+E/3)}\right\}.
```

The boundary assumption gives `epsilon sqrt(n)>=8`, hence
`E/3<=epsilon^2n^2/24`.  The last display is at most
`2 exp(-12n/13)`.  A union bound over all `2^n` spins therefore bounds the
uniform-error failure probability by

```math
u_n\le2\exp\{-(12/13-\log 2)n\}<1/4.                      \tag{QC.S4}
```

Here `n>=64`, so the final numerical inequality has ample slack.

Let `D=sum_e(1-Z_e)` be the number of erased edges.  It has mean `pm`, and a
Chernoff bound gives

```math
Pr\{D<pm/2\}\le e^{-pm/8}<1/4.                             \tag{QC.S5}
```

For the last boundary check, `m>=n^2/3` and

```math
{pm\over8}\ge{\epsilon^2n^2\over48}\ge{4096\over48}>85.
```

Thus, for each fixed `A`, a random mask is simultaneously accurate and has
at most `(1-p/2)m=(1-epsilon^2/4)m` retained edges with probability at least
`1/2`.

Now sample `m+2` masks independently.  A fixed `A` is missed by all of them
with probability at most `2^{-(m+2)}`.  The expected number of uncovered
members of `{-1,1}^m` is at most `1/4`, so some fixed family covers every
`A`.  Discard masks with too many retained edges; they were not good for any
`A`.  Store the chosen mask index and the at most `(1-epsilon^2/4)m` retained
signs.  This proves (QC.S2).  The decoded nonzero coefficients have magnitude
`1/q<=2`, and (QC.1) proves the response bound. `square`

For `n<64/epsilon^2`, the exact `m`-bit code remains valid.  The theorem is a
summary/ambient-cover statement: its decoded centers are sparse
`2`-bounded quadratics, not members of the original sign family.  If centers
must themselves be sign quadratics, use the next theorem.

## 3. Internal sign covers: a strict Hamming code

Let `A` range over `{-1,1}^m`, fix an error budget

```math
E=\epsilon n^{3/2},
```

and suppose `1 <= r := floor(E/2) <= m/2`.  Write

```math
V(m,r)=\sum_{j=0}^r {m\choose j}.
```

### Theorem QC.2 (universal Hamming-cover compression)

There is a fixed codebook of sign matrices, independent of `B`, with at most

```math
K\le
\left\lceil {2^m\over V(m,r)}(m\log 2+1)\right\rceil       \tag{QC.2}
```

codewords such that every sign quadratic `A` has a codeword `C(A)` satisfying

```math
\|P_BH_A-P_BH_{C(A)}\|_\infty\le E.                        \tag{QC.3}
```

Consequently it suffices to store

```math
b_{\rm sign}(n,\epsilon)
\le m-\log_2 V(m,r)+O(\log m)                              \tag{QC.4}
```

bits.  For every fixed `epsilon>0`, as `n` tends to infinity,

```math
b_{\rm sign}(n,\epsilon)
\le {n(n-1)\over2}
-{\epsilon\over4}n^{3/2}\log_2 n
+O_\epsilon(n^{3/2}).                                      \tag{QC.5}
```

Thus the sign family admits a rigorous strict saving of
`Theta_epsilon(n^(3/2) log n)` bits from its exact `m`-bit description.  Its
proved leading rate is nevertheless still quadratic.

#### Proof

Sample `K` independent uniform centers in the Hamming cube `{-1,1}^m`.  A
fixed coefficient vector is missed by all their radius-`r` balls with
probability

```math
\left(1-{V(m,r)\over2^m}\right)^K
\le \exp\{-K V(m,r)/2^m\}.
```

For the value of `K` in (QC.2), the expected number of missed vectors is less
than one.  Hence some choice of centers covers the entire cube.

If `A` and `C` differ on at most `r` edges, then

```math
\|H_A-H_C\|_\infty\le 2r\le E.
```

Equation (QC.1) proves (QC.3), and taking the binary logarithm of (QC.2)
proves (QC.4).  Finally,

```math
V(m,r)\ge {m\choose r}\ge (m/r)^r.
```

Since `m=(1/2+o(1))n^2`, `r=(epsilon/2+o(1))n^(3/2)`, and
`log_2(m/r)=(1/2)log_2 n+O_epsilon(1)`, (QC.5) follows. `square`

### Exact ceiling for this proof architecture

The sphere-covering inequality gives the reverse bound

```math
K_{\rm Ham}(m,r)\ge {2^m\over V(m,r)}.                     \tag{QC.6}
```

Therefore the minimum number of radius-`r` coefficient Hamming balls obeys

```math
\log_2 K_{\rm Ham}(m,r)
=m-\log_2V(m,r)+O(\log m)
=m-{\epsilon\over4}n^{3/2}\log_2n
 +O_\epsilon(n^{3/2}).                                     \tag{QC.7}
```

This is a decisive ceiling only for the coefficientwise Lipschitz method:
it cannot produce `o(n^2)` bits.  It is **not** a response-entropy lower
bound, because a response ball may contain coefficient vectors that are far
apart in Hamming distance but cancel after Boolean optimization.

## 4. Bounded real quadratics: remove the spurious `log n`

Deterministically rounding every coefficient finely enough in `l_1` costs
`O(n^2 log(n/epsilon))` bits.  An elementary simultaneous randomized
rounding argument improves this to `O(n^2 log(1/epsilon))`.

### Theorem QC.3 (uniform discrepancy rounding)

Let every `a_ij` lie in `[-1,1]`, let `n>=2`, and put

```math
L=\lceil 2/\epsilon\rceil,
\qquad
G_L=\{-1+2j/L:0\le j\le L\}.
```

For every `A` there is an `Ahat in G_L^m` for which

```math
\|H_A-H_{\widehat A}\|_\infty
\le \epsilon n^{3/2},
\qquad
\|P_BH_A-P_BH_{\widehat A}\|_\infty
\le \epsilon n^{3/2}.                                     \tag{QC.8}
```

Hence bounded real quadratics have a response code using at most

```math
b_{\rm bd}(n,\epsilon)
\le m\,\lceil\log_2(L+1)\rceil
=O\!\left(n^2\log(1+1/\epsilon)\right)                    \tag{QC.9}
```

bits, uniformly over the bridge.

#### Proof

Independently round each coefficient to one of its two neighboring grid
points so that its rounded value has expectation `a_ij`.  If
`xi_ij=ahat_ij-a_ij`, then the variables are independent and mean zero, and
the range of each `xi_ij x_i x_j` has length at most `delta=2/L<=epsilon`.
For a fixed `x`, Hoeffding's inequality gives

```math
Pr\left\{\left|\sum_{i<j}\xi_{ij}x_ix_j\right|\ge t\right\}
\le 2\exp\{-2t^2/(m\delta^2)\}.
```

Take

```math
t=\delta\sqrt{m(n+2)\log 2/2}.
```

A union bound over the `2^n` spin assignments has failure probability at
most `1/2`.  Some rounding therefore works simultaneously for all `x`.
Moreover, for `n>=2`,

```math
t\le \sqrt{(\log 2)/2}\,\delta n^{3/2}
<\epsilon n^{3/2}.
```

The response conclusion follows from (QC.1), and the grid count gives
(QC.9). `square`

The theorem is existential as stated.  It gives a codebook/entropy bound;
it does not claim that finding the successful simultaneous rounding is
polynomial time.

## 5. Why switching and elementary response constraints do not change the rate

For `s in {-1,1}^n`, let `(A^s)_ij=s_i s_j a_ij`.  Although

```math
H_{A^s}(x)=H_A(s\mathbin\odot x),
```

changing variables in the response also changes the bridge from `B` to
`D_sB`.  Thus switching is not a symmetry of a fixed generic bridge.  Even
if one grants the largest usual switching-permutation quotient, its acting
group has size at most

```math
2^{n-1}n!,
```

so the number of sign-coefficient orbits is at least

```math
{2^m\over 2^{n-1}n!}
=2^{m-O(n\log n)}.                                         \tag{QC.10}
```

Gauge identification alone therefore cannot make the description
subextensive in `n^2`.  For a labeled fixed bridge,
usually not even this quotient is available.

Every quadratic response also satisfies

```math
F_A(-y)=F_A(y),                                             \tag{QC.11}
```

and flipping one query bit changes `F_A` by at most `2n`.  The first fact
only halves the response domain.  At accuracy `epsilon n^(3/2)`, the second
permits interpolation only from a Hamming net of radius
`O(epsilon sqrt(n))`; such a net still has exponentially many query points.
Neither constraint yields a competitive response-table code.

## 6. What remains genuinely open

The following common routes all stop at a quadratic leading bit count.

1. Sparse unbiased sampling proves only a constant-factor rate reduction;
   its Bernstein balance erases `Theta(epsilon^2)` of the edges.
2. Coefficient Hamming/`l_1` control is exhausted by its matching
   sphere-cover bound.
3. Switching and vertex permutations save at most `O(n log n)` bits, and are
   generally broken by the fixed labeled bridge.
4. Operator-norm approximation needs error `O(epsilon sqrt(n))`, because
   `|x^T\Delta A x|<=n||\Delta A||_op`; generic volumetric covers of the
   resulting `Theta(n^2)`-dimensional matrix body still have
   `Theta(n^2)` logarithmic size at fixed `epsilon`.
5. A Frieze--Kannan cut decomposition at error
   `epsilon n^(3/2)` uses the normalized accuracy
   `delta=epsilon/sqrt(n)`.  Its standard `O(delta^(-2))` cut terms are then
   `O(n/epsilon^2)` in number, and recording their vertex subsets already
   costs `O(n^2/epsilon^2)` bits.

These are ceilings on proof methods, not an impossibility theorem for the
response family.  A true `o(n^2)` upper bound could still exist, but it must
show that the maximization map identifies large, Hamming-distant sets of
quadratic coefficients.  Equivalently, it needs a new constraint on the
joint optimizer profile

```math
y\longmapsto\arg\max_x\{H_A(x)+x^TBy\}
```

that is special to quadratic `H_A` and a fixed dense sign bridge.  Neither
density of `B`, evenness (QC.11), query Lipschitzness, nor the usual gauge
group supplies such a constraint.

## 7. Verdict

For genuine sign-quadratic inputs, the strongest universal summary bound
established here is

```math
\boxed{
b_{\rm sign}(n,\epsilon)
\le \left(1-\epsilon^2/4\right){n(n-1)\over2}+O(\log n) .}
```

It is a strict constant-fraction compression theorem at additive
`epsilon n^(3/2)`, valid for every dense or nondense bridge.  The decoded
surrogate is sparse and `2`-bounded.  If code centers are required to remain
sign-valued, the internal Hamming cover instead gives

```math
{n(n-1)\over2}
-{\epsilon\over4}n^{3/2}\log_2n+O_\epsilon(n^{3/2})
```

bits.  Both rates are still quadratic.  No rigorous lower bound ruling out a
more structural `o(n^2)` response code is proved here.
