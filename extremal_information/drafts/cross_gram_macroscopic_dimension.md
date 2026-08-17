# The macroscopic dimension of a cross-Gram table

**Status.** Rigorous task-local upper theorem plus a matching-order Walsh
example.  This tests whether the `O(p^2)` table in SA.3 represents
`Theta(p^2)` independently variable macroscopic information.

The answer is negative for the strongest natural notion of independence: a
full affine bit cube with fixed one-port self states.  Positivity of the two
spectral Gram sectors gives a total squared-amplitude budget `O(p)`.  Hence
only `O(p)` raw cross coordinates can vary by a fixed amount.  The exact
table can still have `Theta(p^2)` small entries, and this theorem does not
rule out a nonlinear code which aggregates many `o(1)` entries through a
genuinely state-dependent broadcast.

The total-system normalization and the archived state-local flux ceiling are
addressed explicitly in Sections 4--5.

## 1. A matrix-cube lemma inside the PSD cone

For an unordered pair `e={i,j}`, let

```math
D_e=e_ie_j^T+e_je_i^T.                              \tag{CG.1}
```

### Lemma CG.1 (fixed-diagonal PSD cubes have linear energy)

Let `K_0` be a positive semidefinite `p by p` matrix with diagonal entries
in `[0,1]`.  Let `e_1,...,e_h` be distinct off-diagonal pairs and let
`eta_1,...,eta_h` be real.  Suppose

```math
K_\sigma=K_0+\sum_{a=1}^h\sigma_a\eta_aD_{e_a}
\succeq0
\qquad\hbox{for every }\sigma\in\{+-1\}^h.         \tag{CG.2}
```

Then

```math
\boxed{\sum_{a=1}^h\eta_a^2\le {p\over2}.}         \tag{CG.3}
```

In particular, if every `|eta_a|>=delta`, then

```math
h\le {p\over2\delta^2}.                            \tag{CG.4}
```

#### Proof

First suppose `K_0` is positive definite with diagonal one.  Put

```math
S_\sigma=\sum_a\sigma_a\eta_aD_{e_a},
\qquad
T_\sigma=K_0^{-1/2}S_\sigma K_0^{-1/2}.            \tag{CG.5}
```

Both sign words `sigma` and `-sigma` occur in (CG.2), so

```math
I+-T_\sigma\succeq0.
```

Every eigenvalue of `T_sigma` lies in `[-1,1]`; hence

```math
\operatorname{tr}T_\sigma^2\le p.                 \tag{CG.6}
```

Let `M=K_0^{-1}`.  Averaging (CG.6) over independent signs cancels all cross
terms and gives

```math
\begin{aligned}
p
&\ge\sum_a\eta_a^2
 \operatorname{tr}(MD_{e_a}MD_{e_a})\\
&=2\sum_{a:e_a=\{i,j\}}\eta_a^2
  (M_{ii}M_{jj}+M_{ij}^2).                         \tag{CG.7}
\end{aligned}
```

For a positive definite correlation matrix, `M_ii>=1`: the Schur
complement formula gives
`M_ii=1/(1-c^TB^{-1}c)>=1`.  Thus (CG.7) proves (CG.3).

For general `K_0`, replace it by `K_0+epsilon I` and normalize by its
diagonal.  The normalized center is positive definite with diagonal one;
the coefficient on edge `{i,j}` becomes

```math
{\eta_a\over
 \sqrt{(K_{0,ii}+\epsilon)(K_{0,jj}+\epsilon)}}.
```

The denominators are at most `1+epsilon`.  The positive-definite result
therefore gives

```math
{1\over(1+\epsilon)^2}\sum_a\eta_a^2\le {p\over2}.
```

Let `epsilon` decrease to zero. `square`

The use of every sign word is essential.  A single codebook of PSD matrices
can have many correlated directions; CG.1 concerns genuinely independently
toggleable coordinates.

## 2. The Gram--Rayleigh amplitude law

Let `J` be a real symmetric involution on `R^n`, and let
`w_1,...,w_p` be vectors of squared norm `n`.  Define

```math
G_{ij}={w_i^Tw_j\over n},
\qquad
R_{ij}={w_i^TJw_j\over n}.                         \tag{CG.8}
```

In SA.3 one takes `J=H/r`.  Put

```math
K^+={G+R\over2},
\qquad
K^-={G-R\over2}.                                   \tag{CG.9}
```

These are positive semidefinite, because with
`P_+-=(I+-J)/2`,

```math
K^+_{ij}={w_i^TP_+w_j\over n},
\qquad
K^-_{ij}={w_i^TP_-w_j\over n}.                    \tag{CG.10}
```

Their diagonal entries lie in `[0,1]` and sum to one.

### Theorem CG.2 (only linear many independent macroscopic cross bits)

Suppose a family of `p`-port states is indexed by every
`sigma in {+-1}^h`, has fixed one-port self data, and has an affine
independent cross table

```math
\begin{aligned}
G_\sigma&=G_0+\sum_{a=1}^h\sigma_ag_aD_{e_a},\\
R_\sigma&=R_0+\sum_{a=1}^h\sigma_ar_aD_{e_a},      \tag{CG.11}
\end{aligned}
```

where the `e_a` are distinct off-diagonal pairs.  If every pair
`(G_sigma,R_sigma)` is realizable by vectors as in (CG.8), then

```math
\boxed{\sum_{a=1}^h(g_a^2+r_a^2)\le2p.}            \tag{CG.12}
```

Consequently, if every bit changes its raw Gram--Rayleigh coordinate by

```math
\sqrt{g_a^2+r_a^2}\ge\delta,                       \tag{CG.13}
```

then

```math
h\le {2p\over\delta^2}.                            \tag{CG.14}
```

#### Proof

The two PSD sectors have fixed diagonals because the one-port self data are
fixed.  Their edge amplitudes are

```math
\eta_a^+={g_a+r_a\over2},
\qquad
\eta_a^-={g_a-r_a\over2}.                          \tag{CG.15}
```

Apply Lemma CG.1 to each sector:

```math
\sum_a(\eta_a^+)^2\le p/2,
\qquad
\sum_a(\eta_a^-)^2\le p/2.                        \tag{CG.16}
```

Adding and using

```math
(\eta_a^+)^2+(\eta_a^-)^2={g_a^2+r_a^2\over2}
```

proves (CG.12), then (CG.14). `square`

For Boolean ports this theorem remains valid verbatim; it used only their
Euclidean Gram representation.  In particular it applies to Boolean top
eigenvectors, where `R=G` and the entire budget lies in the positive sector.

## 3. Linear order is attainable in tensor Walsh systems

The `O(p)` conclusion is of the right order.  Let `H_16` and the orthogonal
Boolean top eigenvectors `1,v_0` be as in SA.4.  In

```math
H=H_{16}^{\otimes j},
\qquad n=16^j,                                     \tag{CG.17}
```

the `2^j=n^(1/4)` tensor words

```math
e_c=\bigotimes_{t=1}^j
 \begin{cases}\mathbf1,&c_t=0,\\v_0,&c_t=1,
 \end{cases}
\qquad(c\in\{0,1\}^j)                             \tag{CG.18}
```

are pairwise orthogonal Boolean `+sqrt(n)` eigenvectors.

For even `p<=n^(1/4)`, allocate two distinct basis words `(a_i,b_i)` to
each `i<=p/2`.  Define ports

```math
w_{2i-1}=a_i,
\qquad
w_{2i}=\begin{cases}
b_i,&\sigma_i=0,\\
a_i,&\sigma_i=1.
\end{cases}                                        \tag{CG.19}
```

All cross entries outside the matching vanish.  On matching edge `i`, both
`G` and `R` equal `sigma_i`.  Thus this is a Boolean top-eigenvector cube of
`p/2` independent constant-amplitude cross bits with identical one-port
self states.

Each bit is exposed at one-block scale by the exact SA.4 two-port
composition: the repeated pair has cap `5rn/2`, while the orthogonal pair
has cap at most `(1/2+sqrt(2))rn`.  Hence the local gap is at least

```math
(2-\sqrt2)rn=(2-\sqrt2)n^{3/2}.                    \tag{CG.20}
```

This construction attains `Theta(p)` independent SA.4-scale bits, matching
CG.2 in order.  It does not by itself turn those `p/2` local pair queries
into one total-scale scalar packing; that normalization is addressed next.

## 4. Total-system normalization

With `p` SA ports of width `m=r=sqrt(n)`, the full composed order is

```math
N=n+p\sqrt n.                                      \tag{CG.21}
```

A fixed-pair gap of size `c n^(3/2)` has normalized size

```math
{c n^{3/2}\over N^{3/2}}
=c\left({n\over n+p\sqrt n}\right)^{3/2}.          \tag{CG.22}
```

It remains a fixed positive total-scale gap only for `p=O(sqrt n)`; it
vanishes when `p/sqrt n` diverges.  The tensor construction in Section 3 has
`p<=n^(1/4)`, so order normalization alone does not kill the local pair
scale.  This is only a kinematic check: activating all other ports may move
the optimizer and erase a local SA.4 gap, so it is not a simultaneous
`p`-port packing theorem.

More generally, consider a fixed-arity Gram-based declared composition whose
response changes by at most

```math
Lrn\sqrt{g_a^2+r_a^2}                               \tag{CG.23}
```

when only independent coordinate `a` changes.  This includes any uniformly
Lipschitz fixed-arity function of the normalized SA.3 features; SA.4 has an
absolute constant scale.  A bit separated by `epsilon N^(3/2)` must have

```math
\sqrt{g_a^2+r_a^2}
\ge {\epsilon\over L}\left({N\over n}\right)^{3/2}. \tag{CG.24}
```

Theorem CG.2 then gives

```math
h\le {2L^2\over\epsilon^2}
p\left({n\over N}\right)^3
\le {2L^2\over\epsilon^2}p.                        \tag{CG.25}
```

Thus a fixed-arity feature algebra cannot turn the quadratic table into
`Theta(p^2)` independently toggleable total-scale bits.  When
`p>>sqrt(n)`, even a unit-amplitude single pair cannot meet (CG.24).

Equation CG.25 is deliberately scoped to independently toggleable features
and fixed-arity Lipschitz exposure.  It does not prohibit a code whose
states differ in `Theta(p^2)` microscopic entries and whose query aggregates
them collectively.

## 5. Compatibility with the archived total-scale no-go

Theorems 21.23/TC.1--TC.2 rule out a different tempting construction:
place one hidden Gram/flux bit in each of many disjoint Walsh child blocks
and add public connectors.  With `k` blocks, the entire state-local response
diameter is only `O(N^(3/2)/sqrt(k))`, and fixed-distortion packing entropy
is bounded independently of `k`.  Coding the local bits does not restore a
positive total-scale rate.

The port cube above does not contradict that theorem.  A port vector changes
`n sqrt(n)` old--new sign coefficients and is therefore a state-dependent
cross-block object, outside TC.1's onsite hypothesis.  Conversely, the local
SA.4 exposure of `p/2` matching bits is not automatically a total-system
packing: all `p` port variables and their completion must be included in
`N`, exactly as in CG.21--CG.22.

This yields a sharp research boundary:

1. `Theta(p^2)` exact entries may be needed to reconstruct an arbitrary
   one-layer spherical certificate.
2. Only `O(p)` of those entries can be independently toggled at constant
   raw amplitude with fixed self data.
3. A proposed `Theta(p^2)` macroscopic information lower bound must therefore
   use many vanishing-amplitude entries coherently, via a nonlocal
   state-dependent broadcast.  Counting the entries or querying them one at
   a time is invalid at total scale.

The last possibility is not closed by CG.2.  It is the precise surviving
route rather than a hidden claim of quadratic minimality.

## 6. Research judgment

The `O(p^2)` table in SA.3 is an exact reconstruction state, not an intrinsic
count of macroscopic independent degrees of freedom.  At the SA.4 scale its
coordinatewise independently toggleable constant-amplitude dimension is
`Theta(p)`, with matching upper and lower orders.  This says nothing by
itself about dense affine directions, nonlinear codes, or collective metric
entropy.  It is a real compression signal, but not yet an `O(p)` sufficient
state: small cross entries can influence collective high-arity queries.

The next discriminating theorem is a metric-entropy bound for the elliptope
pair `(G,R)` under the **collective** multi-port cap pseudometric.  Either its
fixed-distortion entropy is `exp(O(p))`, yielding a genuine approximate
quotient, or a dense small-amplitude code plus one exact-sign broadcast query
will provide the missing quadratic lower bound.
