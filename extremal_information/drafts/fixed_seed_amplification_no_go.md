# A parity-pole obstruction to diluting the five-port coherence seed

**Status.** Rigorous no-go theorem for a broad, explicitly defined
tensor/replica/direct-mixture amplification class, with a quantitative
perturbation version and finite verifier. The theorem concerns prescribed
selector certificates. It neither compares full trust optima nor rules out
an exact-sign construction which changes the seed embedding nonlocally.

## 1. The seed and four logically different notions

Let `H_0` and `w_1,...,w_5` be the order-16 exact-sign seed of Theorem
ES.2. Normalize

```math
T_0={H_0\over4},
\qquad
\langle f,g\rangle={1\over16}f^Tg.
```

For every odd `S subseteq[5]`, put `z_S=product_(i in S)w_i`, and let
`x=Maj_5(w_1,...,w_5)`. The exact seed identities are

```math
\langle z_S,T_0z_S\rangle={13\over16}
\quad(|S|\text{ odd}),
\qquad
\langle x,T_0x\rangle={1\over2}.                \tag{PA.1}
```

It is important not to conflate the following statements.

1. **PSD defect geometry.** Since `T_0` is a self-adjoint contraction,
   `I-T_0` is positive semidefinite. The proof below first takes place in
   this normalized Hilbert-space category.
2. **Exact signs.** `H_0` is entrywise signed, and tensoring it with another
   exact-sign matrix remains exact-sign.
3. **Tightness of the spectral roof.** If an auxiliary Boolean vector is a
   positive top eigenvector, the tensor roof is Boolean-attained. For a
   general auxiliary vector it need not be.
4. **The full trust optimum.** Equations such as (PA.1) evaluate one
   prescribed same-field selector. A different Boolean spin may have a
   larger quadratic value. Nothing below asserts separated full optima.

The question addressed here is narrow: can the fixed five-port gap be
diluted into vanishing deficits for **every** active product pole while its
prescribed majority witness retains constant defect, using the usual
tensor, odd-monomial, replica, or direct-mixture operations? The answer is
no for the class below.

## 2. Odd-monomial selector lifts retain an exposed bad pole

Let `p` be odd. On an auxiliary probability space let `u` be Boolean and
let `S` be a self-adjoint contraction. For each `j in[p]`, choose an odd
set `A_j subseteq[5]` and a sign `sigma_j`, and define the lifted port

```math
v_j=\sigma_j z_{A_j}\otimes u.                  \tag{PA.2}
```

Assume only that the declared selector presentation preserves the bad
witness projectively:

```math
Maj_p(v_1,...,v_p)=\sigma x\otimes u             \tag{PA.3}
```

for one global sign `sigma`. Repeating seed ports, adding cancelling pairs,
permuting them, and replacing a port by any odd product of seed ports are
all included. In particular, no balanced-replication hypothesis is needed.

### Lemma PA.1 (parity pole)

The full active product pole of (PA.2) is

```math
\prod_{j=1}^pv_j=\sigma' z_A\otimes u,
\qquad
A=A_1\mathbin\triangle\cdots\mathbin\triangle A_p,          \tag{PA.4}
```

where `|A|` is odd. Consequently `A` is nonempty and (PA.1) applies.

#### Proof

Coordinate products of Boolean monomials add their exponent vectors over
`F_2`, which gives (PA.4). Parity is a linear functional over `F_2`.
Every `A_j` has parity one and `p` is odd, hence `A` has parity one. The
auxiliary factor is `u^p=u`. `square`

The full set `[p]` has a nonzero Fourier coefficient for odd majority:

```math
\widehat {Maj_p}([p])
=(-1)^{(p-1)/2}{\binom{p-1}{(p-1)/2}\over2^{p-1}}\ne0.       \tag{PA.4a}
```

Thus (PA.4) is one of the active product poles inspected by the robust-selector
certificate. Its coefficient may be small, but the proposed
vanishing-**every**-marginal premise cannot omit it.

### Proposition PA.1b (the top-monomial subgroup cannot represent the witness)

There is an even broader one-block obstruction on the **actual** 16-row
seed support. Allow arbitrary subsets `A_j subseteq[5]`, not necessarily
odd, and port-dependent auxiliary Boolean factors:

```math
v_j=\sigma_jz_{A_j}\otimes u_j.                  \tag{PA.4b}
```

The seed monomials attaining the positive spectral roof are exactly the
subgroup

```math
K=\{0,5,9,12,17,20,24,29\}\subseteq\mathbb F_2^5,           \tag{PA.4c}
```

where masks use the port order `w_1,...,w_5`. Every monomial outside `K`
has normalized Rayleigh value at most `13/16` (the next shell is exactly
`52/64`). More explicitly, the 32 exponent masks split into energies
`64,52,48,32` with multiplicities `8,16,4,4`, and
`K=span_{F_2}{5,9,17}` is the energy-64 class.

Suppose the output of (PA.4b) agrees projectively with the prescribed
selector `x` on every seed row, for each fixed auxiliary coordinate. Then
some `A_j` lies outside `K`. Since every singleton coefficient of odd
majority is nonzero, this `v_j` is an active product pole. Under
`T_0 tensor S` its ratio is at most

```math
{13\over16},                                                  \tag{PA.4d}
```

so its deficit is at least `3/16`.

Indeed, at zero-based seed rows `3` and `10`, every `K`-monomial takes the
same value at the two rows, whereas `x(3)=-1` and `x(10)=1`. If every
`A_j` belonged to `K`, every input in (PA.4b), and hence their pointwise
majority, would agree at those two rows for each auxiliary coordinate.
This contradicts the premise. The Rayleigh bound follows by tensor
factorization and `||S||<=1`; direct seed enumeration shows that all
non-top monomial ratios lie in `[1/2,13/16]`.

Thus **no monomial tensor lift preserving the concrete bad seed response
can have all active product deficits tend to zero**, even with
port-dependent decorations and even-monomial voters.

This proposition is a one-block/tensor statement. In an abstract direct
mixture, different blocks can put their forced non-top singleton at different
port labels. Averaging all odd Fourier channels closes that apparent escape.

### Theorem PA.1c (top-coset packing survives arbitrary direct mixtures)

Let `mathcal O_p` be the `2^(p-1)` odd subsets of `[p]`. Consider an
arbitrary weighted direct mixture of blocks. On block `a`, take

```math
T_a=T_0\otimes S_a,
\qquad
v_(a,j)=\sigma_(a,j)z_(A_(a,j))\otimes u_(a,j),               \tag{PA.4f}
```

where the seed subsets and auxiliary Boolean factors may depend on both
the block and the port. Assume only that the majority output agrees
projectively with `x` on every concrete seed row on every block. Then
the global product-pole deficits satisfy

```math
\boxed{
{1\over2^{p-1}}\sum_{B\in\mathcal O_p}d_B\ge {3\over32},
\qquad
\max_{B\in\mathcal O_p}d_B\ge {3\over32}.}       \tag{PA.4g}
```

Thus even blockwise relabelling, port-dependent decorations, and arbitrary
monomial parity patterns cannot dilute this fixed seed into vanishing
individual deficits.

#### Proof

On block `a`, let `L_a:F_2^p->F_2^5` be the exponent map

```math
L_a(1_B)=\mathbin\triangle_{j\in B}A_(a,j).
```

Proposition PA.1b says the image of `L_a` is not contained in `K`. Choose a
linear functional on `F_2^5/K` which is nonzero on that image. Its pullback
is a nonzero vector `q_a in F_2^p`. If `B` is uniform on the affine
hyperplane `mathcal O_p`, the functional `q_a dot 1_B` is identically one
when `q_a=1`, and is balanced otherwise. It is therefore one with
probability at least `1/2`, and whenever it is one, `L_a(1_B)` lies outside
`K`.

The corresponding `B`-product in (PA.4f) has a non-top seed monomial and
an arbitrary Boolean auxiliary factor. Its Rayleigh ratio is at most
`13/16`, so its deficit is at least `3/16`. Every other deficit is
nonnegative because `T_a` is a contraction. Hence the average over odd
`B` is at least `(1/2)(3/16)=3/32` on every block. Average over the direct
mixture and exchange the two finite averages. The maximum dominates the
average. `square`

The constant `3/32` is the universal quotient-counting guarantee, not claimed
optimal for selector representations of `x`. When every port has odd seed
degree, Lemma PA.1 gives the sharper exposed full pole and the quantitative
selector relation below.

## 3. Tensor no-go and its sharp defect relation

Use the tensor contraction

```math
T=T_0\otimes S,
\qquad
s=\langle u,Su\rangle\in[-1,1].                 \tag{PA.5}
```

Let `d_sel` be the deficit of (PA.3), and `d_full` that of (PA.4), both
relative to the normalized spectral roof one.

### Theorem PA.2 (one-factor tensor/selector dilution is impossible)

Every odd-monomial selector lift satisfying (PA.2)--(PA.3) obeys

```math
\begin{aligned}
d_{sel}&=1-{s\over2},\\
d_{full}&=1-{13s\over16},                        \tag{PA.6}\\
\boxed{d_{full}\ge {3\over8}d_{sel}}.
\end{aligned}
```

In fact

```math
d_{full}-{3\over8}d_{sel}={5\over8}(1-s)\ge0.  \tag{PA.7}
```

Therefore

```math
\max_{B:\widehat{Maj_p}(B)\ne0}d_B\longrightarrow0
\quad\Longrightarrow\quad d_{sel}\longrightarrow0         \tag{PA.8}
```

throughout this class. A constant selector defect cannot coexist with
vanishing individual product deficits.

#### Proof

Tensor quadratic values multiply. Equations (PA.1), (PA.3), and (PA.4)
give (PA.6). Subtraction gives (PA.7). Since the full product is active,
its deficit is bounded by the maximum in (PA.8). `square`

The constant `3/8` is sharp within this class: if `Su=u`, then
`d_full=3/16` and `d_sel=1/2`, giving equality. Thus attaching arbitrarily
many positive top tensor factors, or common-top tensor amplification of
the exact seed, preserves rather than dilutes the obstruction.

### Exact-sign specialization

If `S=K/r_K`, where `K` is a symmetric exact-sign matrix and
`r_K>=||K||_op`, then `H_0 tensor K` is again entrywise signed and (PA.6)
holds at roof `4r_K`. Moreover

```math
\operatorname {tr}(H_0\otimes K)=0.
```

Deleting its diagonal therefore leaves **every** Boolean quadratic energy
unchanged. The operator norm of the hollowed signing can differ from the
declared roof by at most one; this is a relative `o(1)` issue whenever
`r_K->infinity`. If `Ku=r_Ku`, the declared roof is Boolean-tight and the
two deficits are exactly `3/16` and `1/2` at every tensor scale.

Thus PA.2 is an exact-sign no-go, not merely a weighted-contraction one,
for genuine tensor decoration. It is not an all-purpose exact-sign no-go:
a nonlocal completion that is not a tensor decoration may change the
relevant coherence at leading order.

### Corollary PA.3 (recursive selector products make a marginal worse)

Take `L` independent copies of the seed and child

```math
T_L=T_0^{\otimes L}.
```

On each copy first apply five-input majority. Compose the resulting `L`
bits by any Boolean outer selector `F` for which
`widehat F([L])` is nonzero. (Odd-input majority is the principal example.)
Then the Fourier expansion of the composed selector contains the full
`5L`-leaf product with coefficient

```math
\widehat F([L])\widehat {Maj_5}([5])^L\ne0.       \tag{PA.8a}
```

Its product pole is `z_[5]^(tensor L)`, whose normalized Rayleigh value and
deficit are

```math
\left({13\over16}\right)^L,
\qquad
1-\left({13\over16}\right)^L.                   \tag{PA.8b}
```

The same induction applies to a finite tree of odd-majority gates: the
full leaf coefficient is the product of the nonzero top-degree gate
coefficients. Consequently recursive selector products do not dilute the
seed. They force at least one individually inspected active product pole
*away* from the top, with deficit tending to one as the number of seed
factors grows.

This is an exact-sign and Boolean-roof-tight statement: `H_0^(tensor L)`
is exact-sign and has the Boolean top pole `1^(tensor L)`. Its trace is
zero, so hollowing again preserves all Boolean energies.

## 4. Direct mixtures cannot hide seed-supported loss

Direct sum is naturally an operation on probability-space contractions,
not on dense exact-sign matrices. Let blocks have weights `lambda_a`, and
on each block use an odd-monomial lift of the form above, possibly with a
different `S_a,u_a` and a different odd seed pole in (PA.4). A common
global port list may have different blockwise presentations; only (PA.2)
and (PA.3) are required on each block. Direct-sum quadratic values average,
so PA.7 gives the following.

### Corollary PA.4 (sharper common-factor convex dilution law)

For every direct mixture of seed-generated blocks,

```math
\boxed{d_{full}\ge {3\over8}d_{sel}.}            \tag{PA.9}
```

One may also mix in perfectly synchronized blocks on which both declared
deficits vanish, without changing (PA.9). Hence reducing the total seed
weight to zero reduces its joint selector loss to zero as well; retaining
constant seed-supported joint loss leaves a constant exposed full-product
deficit.

Theorem PA.1c already rules out arbitrary blockwise monomial dilution;
PA.9 additionally relates the exposed marginal quantitatively to the
prescribed selector defect in the common-factor odd-monomial class. These
cover the idealized block-diagonal/dilution operation. A block
diagonal contraction has zero cross-block entries and is not a dense
exact-sign matrix. Filling those entries arbitrarily is outside the
corollary unless the completion has `o(1)` normalized quadratic effect on
the two declared witnesses.

## 5. Stable version

The obstruction is not destroyed by a spectrally negligible completion or
an approximate selector identity. Suppose in one block

```math
\left|\langle X,TX\rangle-{s\over2}\right|\le\eta,
\qquad
\left|\langle Z,TZ\rangle-{13s\over16}\right|\le\eta.       \tag{PA.10}
```

Then the corresponding deficits satisfy

```math
\boxed{
d_Z\ge {3\over8}d_X-{11\over8}\eta.}            \tag{PA.11}
```

The same inequality holds for a weighted direct mixture with `eta`
replaced by the weighted average error. In particular, `eta=o(1)` still
forbids a constant joint defect with vanishing full-product defect.

#### Proof

Write the two errors in (PA.10) as `e_X,e_Z`. The exact calculation becomes

```math
d_Z-{3\over8}d_X
={5\over8}(1-s)-e_Z+{3\over8}e_X
\ge-{11\over8}\eta.
```

Average this identity for the direct-mixture statement. `square`

An additive matrix perturbation `E` with
`||E||_op` at most `eta` times the declared roof satisfies (PA.10) for
every Boolean witness. Thus low-operator-norm exact-sign repairs are
covered. Frobenius-smallness alone is insufficient at the relevant
quadratic scale unless it implies this operator bound.

## 6. What the theorem rules out—and the remaining escape

The five-port seed cannot be upgraded to the vanishing-marginal exact-sign
falsifier by any combination of:

1. common-factor tensor decoration;
2. arbitrary port-dependent monomial tensor decoration which preserves the
   concrete bad seed response;
3. odd replication or cancelling-pair padding of the selector;
4. replacing lifted ports by odd seed monomials;
5. arbitrary convex/direct-sum monomial dilution (with the stronger
   selector-relative constant in the common-factor odd-monomial class); or
6. perturbations negligible in normalized operator norm.

The reason is not a loose Cauchy--Schwarz estimate. It is the exact parity
pole (PA.4), which the declared majority Fourier algebra necessarily tests.

A positive construction must therefore break at least one premise. For
example, it could use a nonlocal exact-sign lift in which no global active
product restricts to one fixed bad seed factor, or a growing seed whose odd
product shell itself approaches the top eigenspace while joint coherence
does not. Either mechanism is genuinely new; tensoring or diluting the
present seed is not enough.

## 7. Verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_fixed_seed_amplification_no_go.py
```

The verifier checks the complete odd seed shell, a nontrivial
odd-monomial/cancelling-pair lift, exact tensor formulas, convex mixtures,
the sharp constant, and the perturbation inequality on an exact rational
grid.
