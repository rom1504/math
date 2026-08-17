# Solution-hidden benchmark: Potts-tree separator response

Status: independently derived before literature lookup; exact quotient,
composition, and rate--distortion proofs audited algebraically and accompanied
by a finite exact verifier.

## 1. Frozen operational problem

Fix an alphabet `Q={0,...,q-1}`, with `q>=2`.  A finite rooted graphical-model
fragment has root variable `x` and is attached to an exposed separator variable
`s` by the ferromagnetic Potts score

```math
\psi_K(x,s)=K\mathbf 1\{x=s\},\qquad K>0.                 \tag{PT.1}
```

All scores are in max-sum convention.  After maximizing every private variable
of the fragment while holding `x=i`, write the resulting conditional optimum as
`u_i`.  The message exposed to the future is therefore

```math
m_j=\max_i\{u_i+K\mathbf 1\{i=j\}\}.                       \tag{PT.2}
```

The derivation below used only (PT.1)--(PT.2), arbitrary future factors at `s`,
and recursive tree composition.  The state and all theorem statements were
frozen before the literature comparison in Section 9.

The image theorem itself does not require the private fragment to be a tree:
only its conditional vector `u` matters.  The explicit composition rule later
specializes to a pairwise Potts tree.  A separator tuple can also be treated as
one super-variable; then `q` is the number of joint separator assignments.

## 2. Exact image of the Potts shield

Let

```math
M=\max_i u_i,\qquad b=\max_jm_j,\qquad r=m-b\mathbf 1.       \tag{PT.3}
```

### Theorem PT.1 (reachable normalized state)

The outgoing state satisfies

```math
b=M+K,
\qquad
r_j=\max\{-K,u_j-M\}.                                      \tag{PT.4}
```

Its exact reachable projective carrier is

```math
\mathcal C_K
=\{r\in[-K,0]^q:\max_jr_j=0\}.                             \tag{PT.5}
```

Moreover every pair `(b,r)` in `R x C_K` is realizable.

#### Proof

Equation (PT.2) can be rewritten as

```math
m_j=\max\{M,u_j+K\}.                                       \tag{PT.6}
```

It is at most `M+K`.  If `i` maximizes `u`, choosing `j=i` attains `M+K`,
which proves the formula for `b`; subtracting it proves (PT.4).  Hence every
normalized image lies in `C_K`.

Conversely, given `r in C_K` and any real `b`, take

```math
u_j=b-K+r_j.                                                 \tag{PT.7}
```

Now `M=b-K`, and because `r_j>=-K`, (PT.6) gives `m_j=b+r_j`.
Thus (PT.5) is the exact image, not merely an outer bound. `square`

The raw `q`-entry message has therefore split canonically into one additive
scalar toll `b` and a bounded `(q-1)`-dimensional control state `r`.

## 3. Arbitrary futures and exact contextual minimality

Maximize an arbitrary future fragment conditional on `s=j`, and call its
conditional optimum `g_j`.  Conversely every finite vector `g in R^q` is
realized by a unary factor on `s`.  The completed response is

```math
R_{b,r}(g)=b+\max_j(r_j+g_j).                                \tag{PT.8}
```

### Theorem PT.2 (coarsest exact state)

1. Two fragments have identical absolute responses for every future if and
   only if their messages `m` are equal.
2. Their response functions differ by a future-independent constant `c` if
   and only if `m-m'=c 1`.
3. Consequently `(b,r)` is the coarsest exact absolute-response state, while
   `r` is the coarsest exact projective state.

#### Proof

The reverse implications are immediate.  For the forward implications, fix
coordinate `j` and choose a finite pinning future

```math
g_j=0,\qquad g_\ell=-L\quad(\ell\ne j),                     \tag{PT.9}
```

where `L` exceeds every score advantage that another coordinate has over `j`
in either message.  Coordinate `j` then maximizes both responses, so the
response difference is exactly `m_j-m'_j`.  Applying this pin for every `j`
recovers all coordinates, either identically or up to the same `c`. `square`

In particular, the formal permutation symmetry of the Potts factor does not
give another quotient under arbitrary futures: a future unary factor can name
and distinguish every label.  The additive gauge is the only universal exact
reduction.

## 4. Exact tree composition

At a tree vertex `v`, suppose child `c` sends

```math
m_c(i)=b_c+r_c(i)                                           \tag{PT.10}
```

to `v`, and let `a_i` be the unary score at `v`.  Define

```math
U_i=a_i+\sum_c r_c(i),\qquad M=\max_iU_i.                   \tag{PT.11}
```

If the edge from `v` to its parent has Potts strength `K_v`, then

```math
b_{v\to p}=\sum_cb_c+M+K_v,
\qquad
r_{v\to p}(j)=\max\{-K_v,U_j-M\}.                          \tag{PT.12}
```

Indeed, the conditional score before the parent edge is
`sum_c b_c+U_i`, so Theorem PT.1 applies.  Formula (PT.12) proves closure for
arbitrary depth and degree, including edge-dependent strengths.

It also proves congruence: replacing any child fragment by another fragment
with the same `(b_c,r_c)` leaves the parent state unchanged.  Induction carries
this equality through every enclosing tree and every future.

In max-plus notation, (PT.12) is the usual combination by addition followed by
max-marginalization.  Normalization merely chooses a section of the additive
projective gauge, so it is compatible with every bracketing of the contraction.

### Binary corollary

For `q=2`, the projective message is the single gap

```math
d=m_1-m_0\in[-K,K].                                         \tag{PT.13}
```

If `d_c` are the child gaps, then

```math
h=(a_1-a_0)+\sum_cd_c,
\qquad
d_{v\to p}=\operatorname{clip}(h,-K_v,K_v).                 \tag{PT.14}
```

Together with `b=max(m_0,m_1)`, the message is reconstructed by

```math
m_0=b-\max(d,0),\qquad m_1=b-\max(-d,0).                    \tag{PT.15}
```

Thus the exact binary control state is a bounded scalar, not a raw boundary
table.

## 5. The operational response metrics

For full messages define

```math
d_{\rm abs}(m,m')
=\sup_g\left|\max_j(m_j+g_j)-\max_j(m'_j+g_j)\right|.       \tag{PT.16}
```

### Theorem PT.3 (contextual isometry)

```math
d_{\rm abs}(m,m')=\|m-m'\|_\infty.                          \tag{PT.17}
```

After optimizing the additive calibration,

```math
d_{\rm proj}([m],[m'])
=\inf_c\|m-m'-c\mathbf1\|_\infty
=\frac12\operatorname{osc}(m-m').                          \tag{PT.18}
```

#### Proof

Put `Delta=m-m'`.  For every vector `z`,

```math
\min_j\Delta_j
\le \max_j(z_j+\Delta_j)-\max_jz_j
\le \max_j\Delta_j.                                       \tag{PT.19}
```

This gives the upper bound in (PT.17), while the pinning futures (PT.9) attain
every coordinate difference.  The best uniform approximation of `Delta` by a
constant is its midrange, with error half its oscillation, proving (PT.18).
`square`

The transfer `T_Ku=m` is also projectively nonexpansive.  If

```math
\alpha\mathbf1+v\le u\le\beta\mathbf1+v,                   \tag{PT.20}
```

monotonicity and additive homogeneity imply the same inequalities after
applying `T_K`.  Taking the smallest `alpha,beta` yields

```math
\operatorname{osc}(T_Ku-T_Kv)
\le\operatorname{osc}(u-v).                                 \tag{PT.21}
```

This is sharp: in the unsaturated binary region, (PT.14) is the identity.
Consequently there is no depth-uniform strict contraction under the stated
assumptions.

## 6. Sharp semantic rate--distortion law

Let `N_abs(epsilon)` be the least number of normalized decoder states needed
to approximate every `r in C_K` within contextual absolute error `epsilon`
while carrying `b` separately.  Define `N_proj(epsilon)` analogously using
(PT.18).

### Theorem PT.4 (minimax code size)

For fixed `q` and `0<epsilon<K/20`,

```math
N_{\rm abs}(\epsilon),N_{\rm proj}(\epsilon)
=\Theta_q\!\left((K/\epsilon)^{q-1}\right),                 \tag{PT.22}
```

and hence

```math
\log_2N(\epsilon)
=(q-1)\log_2(K/\epsilon)+O_q(1).                            \tag{PT.23}
```

#### Upper bound

Write

```math
\mathcal C_K=\bigcup_{j=0}^{q-1}
\{r\in[-K,0]^q:r_j=0\}.                                    \tag{PT.24}
```

Each member lies on at least one of these `q` cube faces.  On a selected face,
snap the remaining `q-1` coordinates to an `epsilon`-mesh.  This gives at most

```math
q(\lceil K/\epsilon\rceil+1)^{q-1}                          \tag{PT.25}
```

states and sup-norm error at most `epsilon`.  Equation (PT.17) gives the
absolute guarantee, and `d_proj<=||.||_infinity` gives the projective one.

#### Lower bound

Fix the face `r_{q-1}=0`, restrict all other coordinates to `[-K,-K/2]`, and
use a `5epsilon`-spaced Cartesian grid.  There are at least

```math
\lfloor K/(10\epsilon)\rfloor^{q-1}                         \tag{PT.26}
```

points up to an inessential endpoint convention.  For two distinct points,
one coordinate difference has magnitude at least `5epsilon`, while the last
coordinate difference is zero.  Their sup distance is at least `5epsilon`
and their projective distance is at least `5epsilon/2>2epsilon`.  By the
triangle inequality, one radius-`epsilon` code cell cannot contain both.
This proves the converse. `square`

For `q=2`, (PT.18) becomes

```math
d_{\rm proj}(d,d')=\tfrac12|d-d'|,                           \tag{PT.27}
```

so the optimal rate is `log_2(K/epsilon)+O(1)` bits.

This is genuinely lossy rather than a restatement of exact dynamic
programming: a continuum of exact response classes is mapped to a finite
codebook, the error is defined by all possible futures, and the packing gives
a matching semantic converse.

If the absolute baseline is also to be encoded and is known to lie in an
interval of length `B`, a scalar grid adds `log_2(1+B/epsilon)+O(1)` bits.
Without any bound on `b`, no finite absolute-response code exists.

## 7. A closed lattice implementation

Assume `K=M eta`, all edge strengths are multiples of `eta`, and all local
score-table entries lie in `eta Z`.  Formula (PT.12) then shows inductively that
every normalized message lies in

```math
\mathcal C_{K,\eta}=\mathcal C_K\cap(\eta\mathbb Z)^q.       \tag{PT.28}
```

There are exactly

```math
|\mathcal C_{K,\eta}|=(M+1)^q-M^q.                          \tag{PT.29}
```

Indeed, there are `(M+1)^q` vectors with coordinates in
`{-M eta,...,0}`, and exactly `M^q` have no zero coordinate.  For `q=2`, this
is the `2M+1`-state gap grid.

Now start from arbitrary real local factors and round each microscopic factor
table once to the nearest `eta`-grid.  If `F` tables are rounded, every complete
assignment score changes by at most `F eta/2`.  Taking a conditional maximum
cannot increase a uniform error, so pointwise

```math
\|m-\widetilde m\|_\infty\le F\eta/2.                       \tag{PT.30}
```

Theorem PT.3 then gives, for every unrounded future `g`,

```math
|R_m(g)-R_{\widetilde m}(g)|\le F\eta/2.                    \tag{PT.31}
```

For total target error `epsilon`, taking

```math
M=\left\lceil {KF\over2\epsilon}\right\rceil,
\qquad \eta=K/M                                             \tag{PT.32}
```

gives a compositional normalized carrier of order
`(KF/epsilon)^(q-1)`.  For an extensive allowance `F delta`, the carrier is
instead of order `(K/delta)^(q-1)`, independent of fragment size.

There are two distinct tradeoffs:

1. Computing an exact real message and quantizing only once at its final
   boundary achieves the optimal static carrier (PT.22), independent of `F`,
   but does not reduce internal exact computation.
2. Rounding microscopic factors once makes every intermediate update exactly
   closed on the lattice.  Its error is charged once per factor rather than
   once per message update.

Repeatedly quantizing messages is not silently harmless.  Child projective
errors add before the nonexpansive transfer, and the binary identity region
allows same-signed rounding errors to accumulate linearly until saturation.

## 8. Exact verifier

Run

```bash
python3 extremal_information/experiments/verify_potts_tree_response.py
```

The verifier uses only integer and rational arithmetic.  It checks:

- the exact clipped image and its surjectivity;
- the lattice count `(M+1)^q-M^q`;
- the binary clamp recurrence;
- absolute and projective contextual metrics using pinning futures;
- exact state composition against brute-force enumeration of random trees;
- projective nonexpansiveness and a sharp equality witness;
- finite covering and packing constructions; and
- the one-time factor-rounding bound under arbitrary rational futures.

Finite enumeration is a consistency check, not a substitute for the proofs.

## 9. Post-freeze literature comparison

The combination/max-marginalization rule is classical valuation-algebra,
bucket-elimination, and junction-tree dynamic programming.  Max-product
message updates are max-plus/tropical contractions.  Approximate bucket
schemes give an important comparison, but use structural factor partitioning
rather than the arbitrary-future response metric above.

Primary sources consulted only after freezing the state and theorems:

- [Shenoy, *Valuation-Based Systems for Discrete
  Optimization*](https://arxiv.org/abs/1304.1121)
- [Dechter, *Bucket elimination: A unifying framework for
  reasoning*](https://ics.uci.edu/~dechter/publications/r76A.pdf)
- [McAuley--Caetano, *Faster Algorithms for Max-Product
  Message-Passing*](https://www.jmlr.org/papers/v12/mcauley11a.html)
- [Wainwright--Jaakkola--Willsky, *MAP estimation via agreement on
  (hyper)trees*](https://arxiv.org/abs/cs/0508070)
- [Dechter--Rish, *Mini-Buckets: A General Scheme for Bounded
  Inference*](https://ics.uci.edu/~dechter/publications/r62.pdf)

The reviewed sources support the classical exact contraction and
reparameterization context.  No reviewed primary source stated the combined
exact Potts image (PT.5), arbitrary-future minimality, sharp contextual entropy
(PT.22), and lattice closure (PT.29).  This is a scoped literature-search
observation, not a claim of priority.

## 10. Benchmark verdict

**Pass, independently predicted.**  Starting only from the fragment,
separator, and arbitrary future factors recovers the classical max-sum
message, identifies its exact additive gauge quotient, proves that no stronger
exact quotient survives arbitrary futures, and derives a sharp lossy
rate--distortion law.  The Potts shield adds the decisive structure: its
projective image is the bounded clipped carrier `C_K`, with a closed rational
grid and a binary saturating-gap realization.
