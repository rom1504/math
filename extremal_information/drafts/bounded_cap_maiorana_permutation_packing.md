# A cap-`1/2` response packing of size `exp(Omega(sqrt(n) log n))`

Status: rigorous theorem.  The algebra is checked by an exact finite
verifier.  The asymptotic counting argument is self-contained.  This
strengthens the fixed-permutation Maiorana--McFarland packing from
`Omega(sqrt(n))` to `Omega(sqrt(n) log n)` response bits.  It is not a
statement about exact minimizers and does not imply convergence of the
original signing problem.

## 1. The result

For a hollow symmetric signing `A`, write

```math
H_A(x)=\frac12x^TAx,
\qquad
Q(A)=\max_x|H_A(x)|,
```

and, for a sign bridge `B`, put

```math
(P_BH_A)(y)=\max_x\{H_A(x)+x^TBy\}.
```

Responses are compared projectively:

```math
d_{\rm proj}(R,S)=\frac12\operatorname{osc}_y(R(y)-S(y)).
```

### Theorem BMP.1 (permutation Maiorana--McFarland packing)

There are absolute constants `c>0` and `q_0` such that, whenever

```math
q=2^m\ge q_0,
\qquad n=q^2,
```

there are a common bridge `B in {+-1}^{n times n}` and a collection
`{A_pi: pi in C}` of hollow symmetric signings with

```math
|C|\ge \exp(cq\log q)
```

such that

```math
\|B\|_{2\to2}=\sqrt n,
\qquad
Q(A_\pi)=\frac12n^{3/2},
```

and, for distinct `pi,sigma in C`,

```math
d_{\rm proj}(P_BH_{A_\pi},P_BH_{A_\sigma})
\ge\frac18n^{3/2}.                                    \tag{BMP.1}
```

Consequently any summary answering all continuations through this one
bridge to uniform error less than `n^(3/2)/16` has at least

```math
\exp(\Omega(\sqrt n\log n))\text{ states},
\qquad
\Omega(\sqrt n\log n)\text{ bits}.                    \tag{BMP.2}
```

The permutation itself gives a matching `O(q log q)`-bit description for
this explicit family.  Thus its response information is
`Theta(sqrt(n) log n)` at the declared scale.

The new ingredient is a counting lemma for approximate self-isometries of
the Walsh bilinear form.  It uses all `q!` permutation parameters before a
single Turan selection, rather than paying the `q` freely chosen Boolean
values one pair at a time.

## 2. Walsh regularization and the permutation family

Let `V=F_2^m`, and index the `n=q^2` coordinates by `(u,v) in V^2`.  Set

```math
W_{(a,b),(u,v)}=(-1)^{a\cdot u+b\cdot v}.
```

For a permutation `pi:V->V`, define the Maiorana--McFarland sign vector

```math
s_\pi(u,v)=(-1)^{u\cdot\pi(v)}.                         \tag{BMP.3}
```

Summing first over `u` gives

```math
(Ws_\pi)(a,b)
=q(-1)^{b\cdot\pi^{-1}(a)}.
```

In particular

```math
y_\pi:=q^{-1}Ws_\pi\in\{-1,1\}^n,
\qquad Wy_\pi=qs_\pi.                                  \tag{BMP.4}
```

Let `b=s_id` and

```math
\mathcal H=D_bWD_b,
\qquad
A=\mathcal H-\operatorname{diag}(\mathcal H).
```

The direct transform identity for `b` says `Wb=qb`.  Hence

```math
\mathcal H\mathbf1=q\mathbf1,
\qquad
\mathcal H^2=nI.
```

Also `tr(mathcal H)=tr(W)=0`.  Therefore `A` is a hollow symmetric signing
and, on Boolean vectors,

```math
H_A(x)=\frac12x^T\mathcal Hx,
\qquad
Q(A)=\frac12qn=\frac12n^{3/2}.                         \tag{BMP.5}
```

Switch this one child by the full permutation family:

```math
A_\pi=D_{s_\pi}AD_{s_\pi}.                             \tag{BMP.6}
```

All children retain the exact cap in (BMP.5).  We use the common bridge
`B=W`.

## 3. The pair coordinate is an inner-product-preservation bias

For two permutations `pi,sigma`, put

```math
w=s_\pi\odot s_\sigma,
\qquad
\tau(v)=v+\pi(v)+\sigma(v),                             \tag{BMP.7}
```

where addition is in `V`.  Since

```math
(b\odot w)(u,v)=(-1)^{u\cdot\tau(v)},
```

summing over `u` once and then taking the remaining inner product gives
the exact identity

```math
\begin{aligned}
w^T\mathcal Hw
&=(b\odot w)^TW(b\odot w)\\
&=q\sum_{x,y\in V}(-1)^{x\cdot y+\tau(x)\cdot\tau(y)}.
\end{aligned}                                           \tag{BMP.8}
```

Thus the normalized Rayleigh coordinate is

```math
\rho(\pi,\sigma)
:=\frac{w^T\mathcal Hw}{qn}
=\mathbb E_{x,y\in V}
 (-1)^{x\cdot y+\tau(x)\cdot\tau(y)}.                 \tag{BMP.9}
```

The right side is the signed advantage with which `tau` preserves the
standard binary inner product.  This identity turns the permutation
parameter into a countable obstruction: a bad pair is an approximate
self-isometry with bias greater than `1/4`.

## 4. Approximate self-isometries are superexponentially sparse

### Lemma BMP.2 (random-function inner-product-bias tail)

Let `tau:V->V` be a uniformly random function (independent uniform values),
and let

```math
\rho(\tau)=\mathbb E_{x,y}
 (-1)^{x\cdot y+\tau(x)\cdot\tau(y)}.
```

For `m>=9`, set `r=m-8`.  Then

```math
\Pr\{\rho(\tau)>1/4\}
\le q^{2r}\exp(-rq/128).                                \tag{BMP.10}
```

#### Proof

For `a,x in V`, let

```math
f_a(y)=(-1)^{a\cdot\tau(y)},
\qquad
\widehat f_a(x)=q^{-1}\sum_y(-1)^{x\cdot y}f_a(y).
```

Then

```math
\rho(\tau)=q^{-1}\sum_x\widehat f_{\tau(x)}(x).        \tag{BMP.11}
```

If this average exceeds `1/4`, more than `q/7` values of `x` satisfy

```math
\widehat f_{\tau(x)}(x)>1/8.                            \tag{BMP.12}
```

Indeed, if their fraction is `p`, the average is at most
`p+(1-p)/8`; an average above `1/4` forces `p>1/7`.

For fixed `a`, normalized Fourier Parseval gives
`sum_x widehat f_a(x)^2=1`.  Hence fewer than `64` values of `x` can satisfy
(BMP.12) with this same `a`.  The image under `tau` of the good set
therefore has more than `q/448` elements.  Its linear span has dimension at
least `m-8`: a subspace of dimension at most `m-9` has at most `q/512`
elements.

It follows that there are `r=m-8` linearly independent values
`a_1,...,a_r` and witnesses `x_1,...,x_r` for which

```math
\widehat f_{a_i}(x_i)>1/8
\quad(1\le i\le r).                                    \tag{BMP.13}
```

Fix such independent `a_i` and arbitrary `x_i`.  For every `y`, the vector

```math
(a_1\cdot\tau(y),\ldots,a_r\cdot\tau(y))
```

is uniform on `F_2^r`; these vectors are independent as `y` varies.
Consequently the `r` Fourier coefficients in (BMP.13) are independent
averages of `q` independent Rademacher signs.  Hoeffding gives probability
at most `exp(-rq/128)` that all of (BMP.13) hold.  There are at most `q^r`
ordered choices of the `a_i` and `q^r` choices of the `x_i`.  The union
bound proves (BMP.10). `square`

### Corollary BMP.2a (permutation-pair tail)

If `pi,sigma` are independent uniform permutations of `V` and `tau` is
defined by (BMP.7), then

```math
\Pr\{\rho(\pi,\sigma)>1/4\}
\le
\exp(2q+2r\log q-rq/128).                               \tag{BMP.14}
```

To see this, first choose `pi,sigma` as independent uniform functions.
Then the values `v+pi(v)+sigma(v)` are independent and uniform, so Lemma
BMP.2 applies.  Conditioning both functions to be permutations costs at
most

```math
\left(\frac{q^q}{q!}\right)^2\le e^{2q},               \tag{BMP.15}
```

using `q! >= (q/e)^q`.

### Corollary BMP.2b (large low-bias permutation code)

For all sufficiently large `q=2^m`, there is a set `C` of permutations
such that

```math
|C|\ge
\exp(rq/128-2r\log q-2q)
=\exp(\Omega(q\log q)),                                \tag{BMP.16}
```

and `rho(pi,sigma)<=1/4` for all distinct `pi,sigma in C`.

Indeed, make a graph on the `q!` permutations, joining a bad pair.  Its
ordered edge density including the diagonal is bounded by (BMP.14).
Turan's bound in the form

```math
\alpha(G)\ge\frac{|G|^2}{|G|+2e(G)}
```

gives an independent set of size at least the reciprocal of this density,
which is (BMP.16).

The conditioning step is important.  A bounded-differences estimate sees
only `q` function values and yields an `exp(-O(q))` tail, too weak for the
`q!` reservoir.  The Fourier-rank argument extracts `m-O(1)` independent
output characters, each carrying a speed-`q` deviation, and obtains the
necessary speed `q log q`.

## 5. Response separation

At query `y_pi`, changing variables `u=D_{s_pi}x` and using (BMP.4) gives

```math
(P_WH_{A_\pi})(y_\pi)
=\max_u\left\{\frac12u^T\mathcal Hu+q\mathbf1^Tu\right\}
=\frac32qn.                                             \tag{BMP.17}
```

For child `sigma` at the same query, the field becomes `qw`, with `w` from
(BMP.7).  Put `K=2qI-mathcal H`.  Since `mathcal H^2=q^2I`,

```math
K^{-1}=\frac{2qI+\mathcal H}{3q^2}.
```

Completing the square on the sphere containing the Boolean cube gives

```math
\begin{aligned}
\frac12u^T\mathcal Hu+qw^Tu
&\le qn+\frac12(qw)^TK^{-1}(qw)\\
&=qn\left(1+\frac{2+\rho(\pi,\sigma)}6\right).
\end{aligned}                                           \tag{BMP.18}
```

For a pair in `C`, this is at most `11qn/8`.  Hence at `y_pi` the response
of child `pi` exceeds that of child `sigma` by at least `qn/8`; at
`y_sigma` the reverse difference holds.  The oscillation is at least
`qn/4`, proving (BMP.1).

## 6. Exact finite verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_bounded_cap_maiorana_permutations.py
```

The verifier uses integer arithmetic to check:

1. the permutation Maiorana--McFarland transform for every permutation at
   `q=2,4`;
2. the Rayleigh/inner-product-bias identity for every ordered pair;
3. `H^2=nI`, `H1=q1`, symmetry, zero trace, and the self-dual identity
   query;
4. the exact bad-pair graph and a maximum low-bias code at `q=4`;
5. exact Boolean responses and the projective `1/8` gap on a six-member
   audited subcode at `n=16`.

The exact graph computation finds a maximum code of `20` of the `24`
permutations at `q=4`; the response enumeration uses six representatives
from it.  These are only identity/regression checks.  The theorem's
asymptotic code comes from the Fourier-rank tail and Turan's theorem.

## 7. What changed and what remains open

This result answers the first bounded-cap upgrade discriminator.

- The `sqrt(n)` barrier of the fixed-permutation Boolean-table family is
  not intrinsic.  Permutation-valued bent data carry an additional
  `log n` factor while preserving the exact cap and the same fixed dense
  bridge.
- The proof identifies a new intermediate invariant: approximate
  inner-product-preservation bias.  Its high-bias set is sparse at speed
  `q log q`, whereas a scalar random-label Hoeffding argument only sees
  speed `q`.
- The construction is still sublinear in `n`.  It neither proves an
  `Omega(n)` lower bound nor a ceiling for all bounded-cap families.
- These cap-`1/2` children are regular-Hadamard switchings, not known
  near-minimizers below `1/2`.  Any claim about the latter still needs new
  rigidity or a different packing.

The next honest question is whether another bent class has description
entropy `Omega(n)` while pair products obey an analogous list-size theorem,
or whether the algebraic maximum for a single regular-Hadamard switching
orbit is `exp(O(sqrt(n) polylog n))`.  The weighted-neighborhood route in
BC.1 remains a genuinely different possible path to a linear rate.
