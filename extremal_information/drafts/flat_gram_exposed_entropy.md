# Exposed entropy of flat Gram-difference landscapes

Status: rigorous task-local draft.  This note does not modify the canonical
theorem list.  It answers the exposed-set question needed by Lemma SC.5a and
records an important scope condition for applying that lemma to a disjoint
compiler.

## 1. A sign-balance lemma for quadratic Boolean polynomials

Let

```math
P(x)=\sum_{i<j}d_{ij}x_ix_j,
\qquad x\in\{\pm1\}^k,
```

be nonzero.  Uniform expectation on the cube is denoted by `E`.

### Lemma FE.1 (both signs occupy constant cube mass)

Every nonzero homogeneous quadratic Boolean polynomial satisfies

```math
\Pr\{P>0\}\ge {1\over324},
\qquad
\Pr\{P<0\}\ge {1\over324}.                         \tag{FE.1}
```

#### Proof

Put `sigma=||P||_2`.  The degree-two Bonami--Beckner inequality at
`q=4` gives

```math
||P||_4\le3||P||_2=3\sigma.                        \tag{FE.2}
```

Interpolation between `L^1` and `L^4` gives

```math
\sigma=||P||_2
\le ||P||_1^{1/3}||P||_4^{2/3},
```

and hence

```math
||P||_1\ge {\sigma\over9}.                         \tag{FE.3}
```

Because `P` has no constant term, `E P=0`.  Therefore

```math
\mathbb E P_+=\mathbb E P_-={1\over2}||P||_1
\ge {\sigma\over18}.                              \tag{FE.4}
```

Cauchy--Schwarz now gives

```math
{\sigma\over18}
\le\mathbb E P_+
\le\sigma\Pr(P>0)^{1/2},
```

and the same estimate for `P_-`.  This proves (FE.1). `square`

The constants are deliberately elementary.  The important feature is that
the entropy deficit is `O(1)`, rather than a positive multiple of `k`.

## 2. Signed differences have a full-entropy exposed bulk

For a symmetric hollow real matrix `D`, write

```math
H_D(x)={1\over2}x^TDx,
\qquad q(D)=\max_x|H_D(x)|.
```

### Theorem FE.2 (uniform signed exposed bulk)

If `D` is nonzero, there are a sign `s in {+-1}`, a point
`p in {+-1}^k`, and a set `X subset {+-1}^k` such that

```math
|X|\ge {2^k\over324},
\qquad
sH_D(x)\ge sH_D(p)+q(D)\quad(x\in X).              \tag{FE.5}
```

Moreover `sH_D(p)=sH_D(-p)=-q(D)`, so this is exactly an antipodal
exposed-set hypothesis.

#### Proof

Choose `s` so that the minimum of `sH_D` equals `-q(D)`, and choose a
minimizer `p`.  Apply Lemma FE.1 to `sH_D` and set

```math
X=\{x:sH_D(x)\ge0\}.
```

The set has at least `2^k/324` elements.  Quadratic homogeneity gives
`H_D(-p)=H_D(p)`, and every `x in X` lies at least `q(D)` above the two
antipodal baseline points. `square`

Combining Theorem FE.2 with Lemma SC.5a gives the following immediate
selector lower bound.  If

```math
F(x)=\max_{q\in[K]}(c_q+b_q\mathbin\cdot x),
\qquad ||b_q||_2\le D_0,
\qquad ||F-sH_D||_\infty\le\eta<q(D)/2,
```

then

```math
\log K
\ge-\log324+{(q(D)-2\eta)^2\over2D_0^2}.           \tag{FE.6}
```

This is stronger for Lemma SC.5a than merely having an exponentially large
near-maximal Hamming ball: the entropy deficit in (FE.6) is constant.

## 3. Consequence for the short-seed alternating-form family

Use the notation of Theorem 21.26.  For `B!=T`, put

```math
D_{B,T}=A_B-A_T,
\qquad c_0={\sqrt2\over16}.
```

The theorem gives, simultaneously for every ordered distinct pair,

```math
q(D_{B,T})\ge c_0k^{3/2},
\qquad ||D_{B,T}||_{2\to2}\le16\sqrt k.             \tag{FE.7}
```

### Corollary FE.3 (simultaneous exposed entropy)

For every `B!=T`, one of the two signed difference landscapes
`+-H_{D_{B,T}}` has an antipodally exposed set of at least `2^k/324`
points and exposed gap at least

```math
c_0k^{3/2}.                                        \tag{FE.8}
```

Consequently any affine-selector realization of that orientation with
uniform error `eta<c_0k^(3/2)/2` and slope norm at most `D_0` obeys

```math
\log K\ge-\log324+
 { (c_0k^{3/2}-2\eta)^2\over2D_0^2}.               \tag{FE.9}
```

If the selectors are auxiliary spins `y in {+-1}^m`, with slopes
`b_y=Cy`, then `D_0<=||C||_op sqrt(m)` and `K<=2^m`.  In particular:

1. For an arbitrary complete sign cross-block, the elementary bound
   `||C||_op<=sqrt(km)` gives, when `eta=o(k^(3/2))`,

   ```math
   m\ge
   \left({1\over256\log2}\right)^{1/3}k^{2/3}(1-o(1)).
                                                               \tag{FE.10}
   ```

2. Under the genuinely flat cross-block hypothesis
   `||C||_op<=Lsqrt(k)`, the same estimate gives

   ```math
   m\ge {1\over16L\sqrt{\log2}}k-o(k).              \tag{FE.11}
   ```

Indeed, (FE.10) follows from
`m log2+log324 >= c_0^2 k^2/(2m^2)` and
`c_0^2=1/128`; (FE.11) follows from
`m log2+log324 >= c_0^2 k^2/(2L^2m)`.

These conclusions are deterministic and uniform after any public seed
satisfying Theorem 21.26 has been fixed.  In the randomized construction of
that theorem, the simultaneous statement therefore holds on the same event
of probability at least `1-exp(-0.74k)` (conditional on the chosen
small-bias label list).  No additional union bound over `(B,T)` is needed.

## 4. What can be said for the absolute landscape

There is also a uniform, but much weaker, absolute-value statement.  It is
useful mainly for showing why the signed orientation matters.

### Proposition FE.4 (spectral Hamming-ball exposure)

Suppose

```math
q(D)\ge ck^{3/2},
\qquad ||D||_op\le L\sqrt k,
\qquad |D_{ij}|\le2.
```

Put `delta=(c/(4L))^2`.  There are `p,z in {+-1}^k` and the Hamming ball

```math
X=\{x:d_H(x,z)\le\lfloor\delta k\rfloor\}
```

such that

```math
|H_D(p)|\le\sqrt2k,
\qquad
|H_D(x)|\ge {1\over2}q(D)\quad(x\in X).            \tag{FE.12}
```

Thus, for fixed `c,L`, the absolute landscape has an antipodally exposed
set of size

```math
\exp\{h(\delta)k-O(\log k)\},                       \tag{FE.13}
```

with gap `q(D)/2-O(k)`, where
`h(delta)=-delta log(delta)-(1-delta)log(1-delta)`.

#### Proof

Parseval gives

```math
\mathbb E H_D(x)^2=\sum_{i<j}D_{ij}^2<2k^2,
```

so some `p` has the first property.  Choose `z` attaining `q(D)`.  If
`d_H(x,z)=d`, then

```math
|H_D(x)-H_D(z)|
\le2||D||_op\sqrt{dk}
\le2Lk^{3/2}\sqrt{d/k}.                             \tag{FE.14}
```

For `d<=delta k`, this is at most `ck^(3/2)/2<=q(D)/2`.
The binomial estimate for a Hamming ball proves (FE.13). `square`

For the short-seed constants `c=c_0` and `L=16`, one may take

```math
\delta={1\over524288}.                              \tag{FE.15}
```

This is positive-rate entropy, but its deficit from the whole cube is
`(log2-h(delta))k`, rather than `O(1)`.  Lemma SC.5a therefore need not give
a useful selector lower bound for the absolute landscape at natural slope
constants.  The full-entropy signed result (FE.9) is the robust conclusion.

This loss is not merely an artifact of the Hamming-ball proof.

### Proposition FE.5 (flat absolute landscapes can have a linear entropy deficit)

Let `k=s^2`, with `s` even, partition the coordinates into `s` blocks of
size `s`, and let `D_ij=2` within a block and `D_ij=0` between blocks.  Then

```math
||D||_op=2(s-1)<2\sqrt k,
\qquad
q(D)=k^{3/2}-k,                                    \tag{FE.16}
```

but, for every fixed `alpha>0`,

```math
2^{-k}\#\{x:|H_D(x)|\ge\alpha k^{3/2}\}
\le\exp\{-\alpha k/4+O(\sqrt k)\}.                \tag{FE.17}
```

Thus spectral flatness and a `Theta(k^(3/2))` maximum do not force an
absolute exposed set with only `O(1)` entropy deficit.

#### Proof

If `M_b` is the magnetization in block `b`, then

```math
H_D(x)=\sum_{b=1}^s(M_b^2-s).                       \tag{FE.18}
```

This proves (FE.16), and also `H_D>=-k`; hence the negative absolute tail in
(FE.17) is empty for all sufficiently large `s`.  For `0<lambda<1/2`, the
Gaussian identity and `cosh u<=exp(u^2/2)` give

```math
\mathbb E\exp\{\lambda M_b^2/s\}
\le(1-2\lambda)^{-1/2}.                             \tag{FE.19}
```

The blocks are independent under uniform cube measure.  Chernoff's bound,
with `lambda=1/4`, applied to
`sum_b M_b^2/s`, proves

```math
\Pr\{H_D\ge\alpha s^3\}
\le
\exp\{-\alpha s^2/4+O(s)\},
```

which is (FE.17). `square`

The matrix in FE.5 is itself the difference of two complete sign matrices
(choose opposite signs inside each block and equal signs between blocks),
although those two individual sign matrices need not be spectrally flat.
It is a counterexample to an argument based only on flatness of the
difference, not a counterexample specific to the short-seed family.

## 5. Scope: when SC.5a may actually be applied

Theorem FE.2 concerns an affine envelope which approximates the signed
difference itself.  It must not be applied blindly when a fixed quadratic
child remains outside the selector maximum.  More generally, if

```math
G(x)=h(x)+\max_q(c_q+b_q\mathbin\cdot x)
```

approximates `f(x)`, then the comparison proof of SC.5a applies to the
**residual** `f-h`, not to `f`.  In a compiler written as

```math
H_{A_B}(x)+\max_y G_T(x,y)
\approx H_{A_B}(x)-H_{A_T}(x),
```

the residual seen by the affine selectors is `-H_{A_T}`, not the difference
`H_{A_B-A_T}`.  Therefore (FE.9) directly obstructs:

- a realization which presents the selected signed difference as one affine
  envelope; or
- a compiler language which exposes both signed orientations and to which
  the appropriate orientation can be applied.

It does **not**, by itself, turn a one-sided compiler of the fixed overlay
`-H_{A_T}` into a lower bound when the large absolute extremum of that
particular orientation occurs on the wrong side.  Likewise, a pointwise
realization only of `|H_{A_B-A_T}|` is governed by Proposition FE.4, not by
the full-entropy signed conclusion.  This orientation/common-baseline issue
is the remaining logical boundary, not lack of exposed entropy in the
signed Gram-difference landscapes.

## 6. Research judgment

The requested uniform deterministic statement is true, and in a form much
stronger than random-base evidence would suggest: signed quadratic
differences always have a constant-fraction exposed bulk.  The argument is
not special to alternating forms.  The short-seed construction supplies the
missing `Theta(k^(3/2))` gap simultaneously for all pairs.

What it does not prove is an unconditional linear lower bound for every
exact-sign disjoint compiler.  That conclusion additionally needs either a
spectrally flat old--new block and access to the correct signed orientation,
or a new two-sided/orientation-free selector inequality.  Under no spectral
assumption, the same argument still gives the scalable `Omega(k^(2/3))`
auxiliary lower bound (FE.10).
