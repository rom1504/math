# A landmark quotient for all future syndrome-radius responses

**Status.** Proved theorem-builder draft, independently audited.  The proof
is self-contained.  This resolves the positive side of Open
Question 1 more strongly than requested: error one already has `o(2^w)`
description bits.  Together with Theorem 8.2, it gives a lattice-scale
transition between error below `1/2` and error one.  It does **not** determine
the optimal `o(2^w)` rate.

## 1. Result

Let `G=F_2^w`.  For a spanning support

```math
S\subseteq G\setminus\{0\},
```

write

```math
\lambda_S(x)=\min\{k:x=s_1+\cdots+s_k,\ s_i\in S\},
\qquad
\rho(S)=\max_x\lambda_S(x).                       \tag{SL.1}
```

For an arbitrary appended support `T`, the future response is

```math
\mathcal R_S(T)=\rho(S\cup T),                    \tag{SL.2}
```

and the complete future-response metric is

```math
d_{\rm resp}(S,S')
=\sup_T|\mathcal R_S(T)-\mathcal R_{S'}(T)|.      \tag{SL.3}
```

The supremum can be taken over all supports `T`; allowing the empty support
causes no problem because `S` and `S'` already span.

Put

```math
V(w,r)=\sum_{j=0}^r\binom wj.                     \tag{SL.4}
```

### Theorem SL.1 (all-future landmark sketch)

For every `0<=r<=w`, there is a deterministic summary `Z_(w,r)(S)` and a
decoder `D_(w,r)` such that

```math
\sup_{S,T}
|D_{w,r}(Z_{w,r}(S),T)-\rho(S\cup T)|\le r,        \tag{SL.5}
```

where `S` ranges over all spanning supports.  The number of bits in the
summary is at most

```math
w^2+K(w,r)\lceil\log_2(w+1)\rceil,                \tag{SL.6}
```

where

```math
K(w,r)\le
\left\lceil{2^w\over V(w,r)}(w\log 2+1)\right\rceil. \tag{SL.7}
```

At radius one one may replace (SL.7) by the sharper explicit bound

```math
K(w,1)\le 2^{w-\lfloor\log_2(w+1)\rfloor}.       \tag{SL.7a}
```

In particular, every future response can be answered to additive error one
using

```math
O\left({2^w\log w\over w}+w^2\right)=o(2^w)      \tag{SL.7b}
```

bits.  Since Theorem 8.2 requires `Theta(2^w)` bits at every uniform error
strictly below `1/2`, allowing one unit of error makes the exact response
information rate, normalized by `2^w`, collapse from one to zero.

Consequently, for fixed `0<delta<1/2` and `r=floor(delta*w)`, the description
length is

```math
2^{(1-h_2(\delta)+o(1))w}=o(2^w).                 \tag{SL.8}
```

Equivalently, the range of the approximate quotient has size

```math
\exp\{2^{(1-h_2(\delta)+o(1))w}\}
=\exp\{o(2^w)\}.                                  \tag{SL.9}
```

Thus arbitrary appended-fragment covering-radius queries admit a
subexponential-in-`2^w` quotient at every fixed positive relative error.

### Corollary SL.2 (an actual response net)

Let `N_resp(w,epsilon*w)` be the least number of actual future-response maps
needed to cover all maps `(R_S)_S` in (SL.3) to radius `epsilon*w`.  For fixed
`0<epsilon<1`, set

```math
r=\left\lfloor{\epsilon w\over2}\right\rfloor.
```

Then

```math
\log_2 N_{\rm resp}(w,\epsilon w)
\le w^2+K(w,r)\lceil\log_2(w+1)\rceil             \tag{SL.10}
```

and, for `epsilon<1`, the right side is

```math
2^{(1-h_2(\epsilon/2)+o(1))w}=o(2^w).             \tag{SL.11}
```

The endpoint `epsilon>=1` is trivial because every radius lies in `[0,w]`.
In fact, using the radius-one summary and one representative per nonempty
summary cell gives a radius-two net with `exp(o(2^w))` elements.

## 2. Proof

### Step 1: a support supplies its own low-entropy geometry

Choose, by a fixed deterministic rule, an ordered basis

```math
B(S)=(b_1,\ldots,b_w)\subseteq S.                 \tag{SL.12}
```

Let `d_B` be Hamming distance in this basis.  Since every basis vector is an
allowed generator,

```math
\lambda_S(x-y)\le d_B(x,y).                       \tag{SL.13}
```

The word metric triangle inequality therefore gives

```math
|\lambda_S(x)-\lambda_S(y)|
\le\lambda_S(x-y)
\le d_B(x,y).                                     \tag{SL.14}
```

Thus the full `2^w`-entry coset-leader table is a one-Lipschitz function on
a Hamming cube whose coordinate system costs only `w^2` bits to declare.

### Step 2: cover the cube by landmarks

There is a binary covering code `C subseteq F_2^w` of Hamming covering radius
`r` and size at most (SL.7).  Here is a self-contained probabilistic proof.
Choose `K` independent uniform centers.  A fixed point is uncovered with
probability at most

```math
\exp\{-K V(w,r)/2^w\}.
```

For the integer in (SL.7), the expected number of uncovered points is less
than one, so some choice covers the cube.

For the sharper radius-one statement, put

```math
m=\lfloor\log_2(w+1)\rfloor.
```

Take an `m by w` binary parity-check matrix whose first `2^m-1` columns are
all nonzero vectors of `F_2^m`, and fill any remaining columns by repetitions.
Its kernel has size `2^(w-m)` and covering radius one: the syndrome of any
word is either zero or is itself one of the displayed columns.  This proves
(SL.7a) without an asymptotic covering-code theorem.

Map this one fixed code into `G` using the ordered basis `B(S)`, obtaining
`C_B`.  Store the basis and the exact integers

```math
(\lambda_S(c))_{c\in C_B}.                        \tag{SL.15}
```

For each `x`, choose a nearest landmark `pi_B(x)` by fixed tie breaking and
define

```math
\widehat\lambda_S(x)=\lambda_S(\pi_B(x)).         \tag{SL.16}
```

Equations (SL.14) and the covering property give

```math
\|\widehat\lambda_S-\lambda_S\|_\infty\le r.     \tag{SL.17}
```

There are fewer than `2^(w^2)` ordered ambient bases, and every stored value
lies in `{0,...,w}`.  This proves (SL.6).  The standard binomial estimate

```math
V(w,\lfloor\delta w\rfloor)
=2^{(h_2(\delta)+o(1))w}                          \tag{SL.18}
```

proves (SL.8)--(SL.9).  Only the elementary existence bound (SL.7) is used;
modern constructions of covering codes near the sphere-covering rate are not
needed for this information theorem.

### Step 3: future composition does not amplify landmark error

Allow `lambda_T` to take value `+infinity` when `T` does not span.  Binary
concatenation is min-plus convolution:

```math
\lambda_{S\cup T}(x)
=\min_{u\in G}\{\lambda_S(u)+\lambda_T(x-u)\}.    \tag{SL.19}
```

If two finite functions differ uniformly by at most `r`, their min-plus
convolutions with the same second function differ uniformly by at most `r`.
Taking a maximum also has Lipschitz constant one.  Hence the decoder

```math
D_{w,r}(Z_{w,r}(S),T)
=\max_x\min_u
 \{\widehat\lambda_S(u)+\lambda_T(x-u)\}          \tag{SL.20}
```

obeys (SL.5).

Finally, two supports with the same summary use the same stored basis and
landmark values.  Both word profiles lie within `r` of the common extension
(SL.16), so their response distance is at most `2r`.  Choosing one actual
support from each nonempty summary cell proves Corollary SL.2. `square`

## 3. Abstract mechanism

The proof isolates a reusable sufficient condition that is stronger than
merely saying that a dynamic program has a boundary state.

### Theorem SL.3 (Lipschitz-interface response compression)

Let `(X,d)` be finite, let `F` be a class of real functions on `X`, each
`L`-Lipschitz and taking values in an alphabet of size `A`.  Suppose every
future environment `e` acts through an operator `P_e` satisfying

```math
\|P_e f-P_e g\|_\infty\le\|f-g\|_\infty,          \tag{SL.21}
```

and the declared response `Q` is also one-Lipschitz in sup norm.  If `C` is
an `r`-net of `X`, then storing `f|_C` gives every response `Q(P_e f)` to
error at most `Lr`, using at most

```math
|C|\log_2 A                                      \tag{SL.22}
```

bits, in addition to whatever fixed metric chart is needed to specify `d`.

#### Proof

Nearest-landmark extension is uniformly within `Lr` of `f`; apply the two
nonexpansive maps. `square`

For syndrome supports, the state-dependent chart is a basis inside `S`, its
cost is `O(w^2)`, `X` is the binary Hamming cube, `L=1`, `P_T` is min-plus
convolution by `lambda_T`, and `Q` is maximum.

This theorem also applies outside linear-code covering radius.  For example,
consider a shortest-path module whose boundary cost profile is Lipschitz on
a product-state interface and whose composition glues modules by min-plus
convolution.  A metric landmark set compresses every future terminal-cost
query with the same non-amplification guarantee.  The result is useful only
when the interface covering number is smaller than the complete state space;
an arbitrary unstructured boundary has no such conclusion.

## 4. Why this is more than a restatement

Each ingredient is classical in isolation:

1. covering codes are Hamming nets;
2. word lengths are Lipschitz in any contained basis; and
3. min-plus convolution is sup-norm nonexpansive.

The theorem combines them into a **future-stable response sketch**.  It is
not the exact dynamic-programming state `lambda_S`, which has `2^w` entries,
and it is not ordinary one-shot source coding, because the same sketch answers
an adversarially chosen continuation `T`.  The resulting description exponent
is strictly below the exact support exponent already at constant additive
error one.

There is an important limitation.  If two approximate landmark profiles are
convolved, both sup-norm errors can add.  Therefore the landmark sketch is not
by itself a closed algebra with error independent of the number of
sketch--sketch compositions.  The hard-core quotient below supplies that
separate property at a larger error scale.

It also explains why the direct-sum lower bound and the positive theorem do
not conflict.  Theorem 8.3 forces `Omega_epsilon(w)` bits; (SL.7b) still uses
`2^{w-o(w)}` bits.  Closing this enormous gap requires exploiting further
algebraic restrictions on word metrics or proving that some such metrics
have large landmark complexity even though they arise from supports.

## 5. Falsification and sharpness tests

The theorem itself has three immediate decisive tests.

1. **Basis test.**  If the stored chart is not generated by edges already in
   `S`, (SL.14) can fail.  The contained-basis hypothesis is essential.
2. **Composition test.**  Any proposed extension to a different query algebra
   must verify the nonexpansive inequality (SL.21).  A Lipschitz profile alone
   is insufficient if future composition amplifies sup error.
3. **Counting test.**  At relative radius `delta`, any Hamming landmark set
   has at least `2^w/V(w,r)` points by the sphere-covering bound.  Thus the
   landmark method itself cannot improve the exponent `1-h_2(delta)`, up to
   subexponential factors.  A substantially smaller syndrome quotient would
   have to use algebraic restrictions on the class `(lambda_S)_S`, not a
   better generic landmark cover.

## 6. Primary literature coordinate

The only external coordinate needed is the classical Hamming covering-code
problem.  Potukuchi and Zhang, *Improved efficiency for covering codes
matching the sphere-covering bound*, arXiv:1902.07408, define a covering code
as a Hamming net and prove linear covering codes with rate approaching
`1-h_2(delta)` at relative radius `delta`; see their Theorem I.1 and the
discussion of the sphere-covering bound.  The present proof deliberately uses
the weaker elementary random-center estimate (SL.7), so no unverified
external hypothesis enters Theorem SL.1.

The primary diameter--density coordinate for the closed quotient below is
Klopsch and Lev, *Generating abelian groups by addition only*,
arXiv:0911.2966.  Their Theorem 2.9 states, for every finite abelian group and
`R>=4`, exactly the bound `|A|<=2|G|/(R+1)` for a generating set of
positive diameter at least `R`.  Their positive diameter agrees with our
word diameter in `F_2^w`.  The short Kneser derivation is included below, so
the imported theorem is a verified literature coordinate rather than a
hidden proof obligation.

## 7. A complementary closed algebra

The preceding sketch optimizes response accuracy.  A different quotient
optimizes algebraic closure.

For an integer `R>=4`, define

```math
q_R(S)=
\begin{cases}
\bot,&\rho(S)<R,\\
S,&\rho(S)\ge R.
\end{cases}                                      \tag{SL.23}
```

### Theorem SL.4 (sparse hard-core quotient)

On spanning supports, `q_R` is an exact quotient of the union semilattice.
Explicitly, define

```math
x\odot_R y=
\begin{cases}
\bot,&x=\bot\text{ or }y=\bot,\\
q_R(x\cup y),&x,y\ne\bot.
\end{cases}                                      \tag{SL.24}
```

Then

```math
q_R(S\cup T)=q_R(S)\odot_R q_R(T)                \tag{SL.25}
```

for all spanning `S,T`; in particular, repeated composition causes no error
accumulation.  Decode the final radius exactly on a retained state and as
`R/2` on `bot`.  Uniform response error is at most

```math
{R\over2}-1.                                      \tag{SL.26}
```

Moreover, for `G=F_2^w`, the number of bits needed to encode the quotient is

```math
O\left({2^w\log R\over R}\right).                \tag{SL.27}
```

Consequently, any threshold `R=R(w)` with

```math
R(w)\longrightarrow\infty,
\qquad R(w)=o(w)                                  \tag{SL.28}
```

gives a closed `exp(o(2^w))`-state composition algebra with `o(w)` radius
error.  For fixed `0<epsilon<1/2` and all sufficiently large `w`, taking
`R=floor(2*epsilon*w+2)` gives error at most `epsilon*w` and
`O_epsilon(2^w*log(w)/w)` bits.

The same state also answers a query by an arbitrary exact appended support,
even if that query support does not span.  On `bot` return `R/2`; on a retained
state use the explicitly retained `S` to compute `rho(S union T)`.  Thus no
full support is hidden in the bottom decoder.  Exact retention of the sparse
high-radius supports is intentional and is paid for by (SL.27).  The closed
binary operation (SL.24) itself is asserted only when both summarized
fragments span, as in the full-rank fragment algebra.

#### Proof and constant audit

If `rho(S)<R`, then

```math
\rho(S\cup T)\le\rho(S)<R                         \tag{SL.29}
```

for every `T`.  Thus the low-radius class is an absorbing ideal under union,
which proves (SL.25), associativity, and closure.  Every spanning support has
radius in `{1,...,w}`.  Conditional on the bottom state, the final radius lies
in `{1,...,R-1}`; its midpoint `R/2` has maximum error `R/2-1`, proving
(SL.26).

It remains to count retained states.  We give the needed diameter--density
bound directly.  Put `A=S union {0}`, `k=R-1`, and let `H` be the period of
the proper sumset `kA`.  Kneser's theorem and `H`-periodicity give

```math
k|A+H|-(k-1)|H|
\le |kA|
\le |G|-|H|.                                    \tag{SL.30}
```

The support spans `G`, while `H` is proper, so `A+H` contains at least two
`H`-cosets and `|H|<=|A+H|/2`.  Since `k>=3`, substituting this into
(SL.30) yields

```math
{k+2\over2}|A+H|\le |G|,
\qquad
|S|+1=|A|\le {2|G|\over R+1}.                    \tag{SL.31}
```

With `N=2^w` and `m=floor(2N/(R+1))`, the quotient therefore has at most

```math
1+\sum_{j=0}^{m}\binom{N-1}{j}                   \tag{SL.32}
```

states.  Since `m/N<=2/5`, the elementary binomial bound gives

```math
\log_2\left(1+\sum_{j=0}^{m}\binom{N-1}{j}\right)
\le 1+m\log_2{eN\over m}
=O\left({N\log R\over R}\right),                 \tag{SL.33}
```

which is (SL.27). `square`

### Abstract hard-core principle

Let `(X,vee)` be a join semilattice and let an integer cost `c` be antitone:

```math
c(x\vee y)\le\min\{c(x),c(y)\}.                \tag{SL.34}
```

Then collapsing the sublevel ideal `{x:c(x)<R}` to one bottom state and
retaining its complement exactly is always a semilattice quotient.  If the
high-cost objects have small metric entropy, this becomes a nontrivial
approximate response algebra.  This applies beyond codes to any monotone
resource system in which adding allowed actions can only lower the extremal
cost; for example, connected graph modules under edge union with diameter as
the query.  The difficult, model-specific input is not the quotient algebra
but a density or counting theorem for the high-cost objects.  Here that input
is the Kneser consequence (SL.31).

This theorem is not merely dynamic programming or convex duality.  Its key
step is a structural rarity law: the objects on which an accurate extremal
response is hard are sparse, while the easy sublevel set is an algebraic
ideal.  It gives a closed quotient without approximating the complete
coset-leader landscape.

## 8. Director-level lesson

This is a positive answer to the stated unrestricted-support dichotomy, but
not yet a universal law of extremal information.  The generative principle
is:

> Future extremal responses compress when the exact interface profile is
> Lipschitz on a domain of low metric entropy and the entire continuation
> algebra is nonexpansive in the profile norm.

It identifies three separate sources of feature growth: geometric entropy of
the interface, amplification by composition, and the entropy of the
high-response hard core.  The syndrome model now escapes the first at error
one and escapes the second at every diverging sublinear error scale.  Neither
mechanism determines the optimal response entropy between the linear block
lower bound and the `O(2^w log w/w)` upper bound.

The genuinely general law surviving this example is a dichotomy:

> Future extremal responses compress either because their exact profile is
> regular on a low-entropy interface and continuation is contractive, or
> because the difficult response states form a sparse complement of an
> absorbing sublevel ideal.

The next nontrivial theorem should test whether these two mechanisms can be
combined into a query-sensitive bound—retain a sparse exposed hard core and
landmark-compress its internal response profiles—or whether one can construct
a model where both interface entropy and hard-core entropy are necessarily
large.
