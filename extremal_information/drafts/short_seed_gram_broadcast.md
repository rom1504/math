# Short-seed public advice for alternating-form Gram broadcast

Status: task-local rigorous draft.  This strengthens the public-advice part
of `state_dependent_gram_broadcast.md`; it does not change the declared
same-support overlay language or make that language closed in the exact-sign
class.

The old proof chose the public base signing with one fully independent bit
per edge.  The point here is that simultaneous spectral flatness of the
whole alternating-form coset is only a linear-moment event.  A
`6k`-wise-independent edge family therefore suffices.  A second, elementary
small-bias argument compresses the public label list.  The resulting shared
data have `O(k log k)` bits rather than `Theta(k^2)` bits.

The sample spaces are completely explicit and a uniformly computable good
seed exists.  The proof does **not** supply a deterministic polynomial-time
algorithm for finding that seed: exhaustive search takes
`2^{O(k log k)}` time, while a random short seed succeeds with exponentially
high probability.

## 1. A linear-description alternating-form sampler

Let `V=F_2^r`, `r>=2`, and put

```math
h=\dim\operatorname{Alt}(V)={r(r-1)\over2}.          \tag{SG.1}
```

A multiset `S` in `V`, with its uniform probability measure `mu`, is
`delta`-biased when

```math
\left|\mathbb E_{p\sim\mu}(-1)^{\ell(p)}\right|
\le\delta                                               \tag{SG.2}
```

for every nonzero linear form `ell` on `V`.

### Lemma SG.1 (a stored linear-size small-bias multiset suffices)

There is a `1/8`-biased multiset `S` of size

```math
s=256r.                                                \tag{SG.3}
```

It can be stored in `sr=256r^2` bits.

#### Proof

Sample `s` independent uniform points of `V`.  For a fixed nonzero `ell`,
Hoeffding's inequality gives

```math
\Pr\left\{\left|{1\over s}\sum_{p\in S}(-1)^{\ell(p)}
                   \right|>{1\over8}\right\}
\le 2\exp(-s/128).                                    \tag{SG.4}
```

A union bound over fewer than `2^r` nonzero forms is strictly below one,
because

```math
(r+1)\log2-s/128=(r+1)\log2-2r<0                     \tag{SG.5}
```

for every `r>=1`.  Hence such a multiset exists. `square`

Repeat every indexed occurrence of `S` exactly `r` times and call the
resulting ordered list `P=(p_1,...,p_k)`.  Thus

```math
k=sr=256r^2.                                          \tag{SG.6}
```

The repetition rule is implicit, so the stored labels still cost exactly
`k` bits, not `kr` bits.

### Lemma SG.2 (every alternating form is seen on a constant edge density)

For every nonzero `B in Alt(V)`,

```math
#\{i<j:B(p_i,p_j)=1\}\ge {k^2\over8}.                \tag{SG.7}
```

Consequently `B -> (B(p_i,p_j))_(i<j)` is an injective
binary linear code of dimension `h` and relative distance greater than
`1/4`.

#### Proof

Let `R=rad(B)`, and let the rank of `B` be `d>=2`.  Fourier expansion of the
indicator of the codimension-`d` space `R`, followed by (SG.2), gives

```math
\mu(R)
\le 2^{-d}+(1-2^{-d})/8
\le {11\over32}.                                      \tag{SG.8}
```

For `p notin R`, the map `q -> B(p,q)` is a nonzero linear form.  Therefore

```math
\begin{aligned}
\mathbb E_{p,q\sim\mu}(-1)^{B(p,q)}
&\le \mu(R)+{1\over8}(1-\mu(R))\\
&\le {109\over256}<{1\over2}.                       \tag{SG.9}
\end{aligned}
```

In particular

```math
\Pr_{p,q\sim\mu}\{B(p,q)=1\}
={1-\mathbb E(-1)^{B(p,q)}\over2}\ge {1\over4}.      \tag{SG.10}
```

The empirical law of `P` is exactly `mu`.  Because `B` is alternating and
symmetric in characteristic two, the `k^2` ordered pairs include zero
diagonal and every nonzero unordered pair twice.  Hence (SG.10) gives

```math
2#\{i<j:B(p_i,p_j)=1\}
=k^2\Pr_{p,q\sim\mu}\{B(p,q)=1\}\ge{k^2\over4},
```

which is (SG.7).  Positive support for every nonzero `B` gives injectivity;
dividing (SG.7) by `binom(k,2)` gives relative distance
`k/(4(k-1))>1/4`. `square`

The sharper constants in (SG.8)--(SG.9) actually give support at least
`147k^2/1024`; the coarser (SG.7) is convenient below.

## 2. An explicit `6k`-wise-independent edge sample space

Let

```math
E={k\choose2},\qquad
d=\lceil\log_2E\rceil,\qquad
q=2^d,
\qquad t=6k.                                          \tag{SG.11}
```

Since `k>=1024`, one has `t<=E<=q`.  Fix distinct field elements
`alpha_e in F_q`, one for every unordered edge `e`.  A seed is a coefficient
vector

```math
u=(u_0,...,u_{t-1})\in F_q^t.                        \tag{SG.12}
```

It produces the hollow symmetric signing

```math
A_u(e)=(-1)^{\operatorname{Tr}_{F_q/F_2}
                    (\sum_{j=0}^{t-1}u_j\alpha_e^j)}. \tag{SG.13}
```

The family is explicit once a canonical irreducible polynomial of degree
`d` and a canonical edge ordering are fixed.  Its seed length is

```math
td=6k\lceil\log_2E\rceil=O(k\log k).                 \tag{SG.14}
```

This logarithmic factor is intrinsic to the chosen derandomization
principle.  Indeed, if a distribution on `E` signs is `t`-wise independent,
then the characters indexed by subsets of at most `floor(t/2)` coordinates
are mutually orthogonal in `L^2` of that distribution.  Consequently

```math
|\operatorname{supp}\mathcal D|
\ge\sum_{j=0}^{\lfloor t/2\rfloor}{E\choose j}.       \tag{SG.14a}
```

For `E=Theta(k^2)` and `t=6k`, (SG.14a) requires seed length
`Omega(k log k)`.  Thus (SG.14) is order-optimal among exact
`6k`-wise-independent sample spaces, though another spectral
derandomization principle could conceivably use fewer bits.

### Lemma SG.3 (short-seed simultaneous spectral flatness)

For the list `P` in Lemma SG.2, some seed `u in F_q^t` has the property that
every `B in Alt(V)` satisfies

```math
\|A_u\odot\chi_B\|_{2\to2}\le8\sqrt k,
\qquad
\chi_B(i,j)=(-1)^{B(p_i,p_j)}.                       \tag{SG.15}
```

Moreover a uniform random seed fails (SG.15) with probability at most

```math
9^k2^h(3/8)^{3k}<\exp(-0.74k).                       \tag{SG.16}
```

#### Proof

Polynomial interpolation shows that the values of a uniform polynomial of
degree below `t` at any at most `t` distinct field points are independent
and uniform in `F_q`.  Their traces are consequently independent fair
bits.  Thus the edge signs in (SG.13) are `t=6k`-wise independent.
Multiplication by any fixed character word `chi_B` preserves this property.

Fix `B` and a unit vector `z in R^k`, and write

```math
X=z^T(A_u\odot\chi_B)z
=\sum_{i<j}2z_iz_j\chi_B(i,j)A_u(i,j).               \tag{SG.17}
```

Put `m=3k`, so `2m=t`.  Expanding `X^{2m}` involves at most `t` distinct
edge signs in every monomial.  Its expectation under (SG.13) therefore
equals its expectation for fully independent Rademacher signs.  The
Khintchine moment bound and

```math
\sum_{i<j}(2z_iz_j)^2
=2\left(1-\sum_i z_i^4\right)\le2                  \tag{SG.18}
```

give

```math
\mathbb E X^{2m}
\le(2m-1)!!\,2^m
={(2m)!\over m!}
\le(2m)^m.                                           \tag{SG.19}
```

Markov's inequality now yields

```math
\Pr\{|X|>4\sqrt k\}
\le\left({2m\over16k}\right)^m
=(3/8)^{3k}.                                         \tag{SG.20}
```

Take a `1/4`-net of the unit sphere with at most `9^k` points and union-bound
over this net and the `2^h` forms.  If no scalar event fails, the standard
symmetric quadratic-form net inequality gives operator norm at most
`2(4sqrt(k))=8sqrt(k)`.  Finally

```math
\begin{aligned}
{1\over k}\log\left(9^k2^h(3/8)^{3k}\right)
&\le \log9+{\log2\over512}+3\log(3/8)\\
&<-0.743,
\end{aligned}                                        \tag{SG.21}
```

because `h/k=(r-1)/(512r)<1/512`.  This proves both existence and
(SG.16). `square`

## 3. Response packing with `O(k log k)` shared data

Fix a good seed from Lemma SG.3 and put

```math
A_B=A_u\odot\chi_B,
\qquad H_B(x)={1\over2}x^TA_Bx.                      \tag{SG.22}
```

As in the original Gram-broadcast construction, declare the complete
same-support additive overlay language `{-H_T:T in Alt(V)}` and write

```math
R_T(B)=\max_x|H_B(x)-H_T(x)|.                        \tag{SG.23}
```

### Theorem SG.4 (short-advice Gram broadcast)

For every `r>=2`, at order `k=256r^2` there are `2^h` exact hollow sign
children, with `h=r(r-1)/2`, such that

```math
R_B(B)=0,
\qquad
{\sqrt2\over16}k^{3/2}
\le R_T(B)\le8k^{3/2}\quad(B\ne T).                 \tag{SG.24}
```

Every child has Boolean cap at most `4k^{3/2}`.  Hence uniform response
error below `(sqrt(2)/32)k^{3/2}` requires at least

```math
h\ge{k\over1024}                                    \tag{SG.25}
```

bits.  All shared public data can be encoded in

```math
k+6k\lceil\log_2\tbinom{k}{2}\rceil+O(\log k)
=O(k\log k)                                          \tag{SG.26}
```

bits.

#### Proof

If `B!=T`, Lemma SG.2 says that `D=A_B-A_T` has magnitude two on at least
`k^2/8` unordered entries.  A random vertex bipartition cuts at least half
of these entries for some partition.  Optimizing signs on one side after
averaging over the other and using the sharp `p=1` Khintchine inequality
gives

```math
Q(D)\ge {\sqrt2\over2}{k^2/8\over\sqrt k}
={\sqrt2\over16}k^{3/2}.                             \tag{SG.27}
```

On the other hand, (SG.15) gives

```math
Q(D)\le {k\over2}\|D\|_{op}\le8k^{3/2},
\qquad
Q(A_B)\le4k^{3/2}.                                   \tag{SG.28}
```

The query indexed by one member of any pair distinguishes their response
vectors by the lower bound in (SG.24), proving the packing claim.  Since
`r>=2`, `h>=r^2/4=k/1024`.  The small-bias multiset costs `sr=k` stored bits;
its repetition rule costs only `O(log k)` bits.  The base signing is
specified by the polynomial seed (SG.14), proving (SG.26). `square`

## 4. What is and is not explicit

1. **Explicit family.**  The map from `(r,S,u,B)` to every coefficient of
   `A_B` is a finite-field formula.  It uses `O(k log k)` shared bits and
   `h=Theta(k)` hidden bits.  No quadratic edge table is stored.
2. **Randomized construction.**  A random `S` has positive success
   probability, and after fixing a good `S`, a random polynomial seed is
   good with probability at least `1-exp(-0.74k)`.  Generating and evaluating
   a candidate is polynomial-time; certifying it by checking all `2^h`
   modulations is exponential-time.
3. **Deterministic uniform construction.**  Enumerate multisets until the
   finite bias inequalities hold, then enumerate polynomial seeds until all
   `2^h` integer matrices pass the exact spectral test.  This halts by
   Lemmas SG.1 and SG.3 and takes `2^{O(k log k)}` time.  Thus the result is
   algorithmic and uniform, but not known to be polynomial-time explicit.
4. **Remaining scope costs.**  The query family still has `2^h` members and
   acts by same-support additive overlay.  The public description is short,
   but the construction does not prove an exact-sign parent closure, a
   compact decoder, or relevance to arbitrary dense near-minimizers.

The finite arithmetic checks are in
[`../experiments/verify_short_seed_gram_broadcast.py`](../experiments/verify_short_seed_gram_broadcast.py).
