# Second-order augmented cut-code amalgamation

Status: exact decomposition and a histogram-level switching theorem.  The
theorem identifies a nonlocal state strictly smaller than the full bridge
response table.  Present finite data certify only constant gains, not the
leading-order gain required for convergence.

## 1. Derivation in the quadratic objective

Let `A` and `D` be signings on disjoint vertex sets of sizes `m` and `n`,
and let `W` be an `m` by `n` sign bridge.  Write

```math
H_A(x)=\sum_{i<j}A_{ij}x_ix_j,
\qquad
H_D(y)=\sum_{i<j}D_{ij}y_iy_j.
```

The block signing `P` has energy

```math
H_P(x,y)=H_A(x)+H_D(y)+x^{\mathsf T}Wy.              \tag{SO1}
```

Replacing `y` by `-y` preserves both internal energies and reverses the
bridge energy.  Hence the absolute parent cap has the exact form

```math
\boxed{
\operatorname{cap}(P)
=\max_{x,y}\left(
 |H_A(x)+H_D(y)|+|x^{\mathsf T}Wy|
 \right).}                                          \tag{SO2}
```

This identity must precede the coding interpretation: the bridge cannot be
treated as a signed term that is allowed to cancel the internal energy,
because the opposite bridge orientation is always another parent spin.

Put

```math
M_A=\operatorname{cap}(A),\qquad
M_D=\operatorname{cap}(D),\qquad
L_W=\max_{u,v}|u^{\mathsf T}Wv|,
```

and define the two nonnegative deficits

```math
\eta(x,y)=M_A+M_D-|H_A(x)+H_D(y)|,
\qquad
\zeta_W(u,v)=L_W-|u^{\mathsf T}Wv|.                 \tag{SO3}
```

For a row switching `p` and column switching `q`, let
`W^(p,q)=diag(p) W diag(q)`.  Equations (SO2)--(SO3) give

```math
\operatorname{cap}(P^{p,q})
=M_A+M_D+L_W-\Gamma(p,q),                            \tag{SO4}
```

where

```math
\Gamma(p,q)=
\min_{x,y}\bigl(
 \eta(x,y)+\zeta_W(p\mathbin\odot x,q\mathbin\odot y)
 \bigr).                                            \tag{SO5}
```

Thus `Gamma` is exactly the state-dependent amount by which amalgamation
beats the independent rectangular term.  It is not a reformulation with an
unspecified error.

## 2. Exact augmented-code meaning

Use the absolute sign version of the augmented cut code

```math
\mathcal C_k^{\pm}
=\{(\sigma z_i z_j)_{i<j}:\sigma,z_i\in\{\pm1\}\}. \tag{SO6}
```

If `E_k=binom(k,2)`, then the covering distance of a signing `A` is

```math
r(A)=d(A,\mathcal C_k^{\pm})={E_k-M_A\over2}.        \tag{SO7}
```

The analogous absolute code on the rectangular bridge gives
`r(W)=(mn-L_W)/2`.  Combining (SO4) and (SO7) yields the exact second-order
covering-radius decomposition

```math
\boxed{
r(P^{p,q})=r(A)+r(D)+r(W)+{\Gamma(p,q)\over2}.}      \tag{SO8}
```

The usual independent deep-hole amalgamation retains only `Gamma>=0` and
therefore pays all of `L_W`.  Any improved covering-radius amalgamation must
produce a positive uniform `Gamma`; (SO8) fixes its normalization and sign.

## 3. Histogram switching theorem

Work projectively, so that

```math
\mathcal X_k=\{\pm1\}^k/\{\pm\mathbf1\},
\qquad
S=|\mathcal X_m\mathbin\times\mathcal X_n|
=2^{m+n-2}.
```

Let

```math
N_\eta(e)=\#\{(x,y):\eta(x,y)=e\},
\qquad
N_\zeta(z)=\#\{(u,v):\zeta_W(u,v)=z\}.             \tag{SO9}
```

> **Histogram amalgamation theorem.** If an integer `g>=0` satisfies
>
> ```math
> \boxed{
> \sum_{e+z<g}N_\eta(e)N_\zeta(z)<S,}               \tag{SO10}
> ```
>
> then some row/column switching of `W` obeys
>
> ```math
> \Gamma(p,q)\ge g,
> \qquad
> \operatorname{cap}(P^{p,q})
> \le M_A+M_D+L_W-g.                                \tag{SO11}
> ```

In particular, when `A` and `D` attain `M_m` and `M_n`, respectively,

```math
M_{m+n}\le M_m+M_n+L_W-g.                           \tag{SO11a}
```

To prove the theorem, choose `p` and `q` uniformly and projectively.  For
each fixed child pair `(x,y)`, the pair `(p odot x,q odot y)` is uniform on
`X_m times X_n`.  The expected number of pairs violating the first inequality
in (SO11) is exactly the left side of (SO10) divided by `S`.  If it is less
than one, a switching with no violation exists.

This state is genuinely smaller than full bridge optimization.  The signed
energy histograms of `A` and `D` determine `N_eta`, and the magnitude
histogram of `W` determines `N_zeta`.  Each has only polynomially many
integer bins, while the complete alignment table has `2^(m+n-2)` entries.
The theorem uses random switching to supply the missing alignment.

It is nevertheless a nonlocal state: proving sharp shell counts may be hard,
and exact computation of a histogram need not be algorithmically easy merely
because its output has few bins.

## 4. Exact leading-order obligation

For a desired residual bridge cost `r_(m,n)`, set

```math
g=L_W-r_{m,n}.
```

Then (SO10) becomes the concrete shell-separation condition

```math
\sum_{e+z<L_W-r_{m,n}}N_\eta(e)N_\zeta(z)<2^{m+n-2}, \tag{SO12}
```

and proves

```math
\operatorname{cap}(P^{p,q})\le M_A+M_D+r_{m,n}.     \tag{SO13}
```

Thus a family with `r_(m,n)=o((m+n)^(3/2))` for comparable blocks would
remove the leading rectangular Gale--Berlekamp term.  More generally, one
can insert the exact residual required by a proposed geometric recurrence in
(SO12).  This is the precise new lemma; it concerns convolution of shell
counts rather than full pointwise bridge responses.

## 5. Balanced shell-entropy diagnostic

Take `m=n`, and normalize both histograms in (SO9) by

```math
S=2^{2n-2}.
```

Let `(X,Y)` and `(U,V)` be independent uniform projective pairs, the first
sampled from the child-pair state space and the second from the bridge-pair
state space.  Setting `g=L_W-r_n`, condition (SO12) is exactly

```math
\Pr\left{
 |H_A(X)+H_D(Y)|+|U^{\mathsf T}WV|
 >M_A+M_D+r_n
\right}<2^{-2n+2}.                                 \tag{SO14}
```

Thus the required large-deviation rate is

```math
\liminf_{n\to\infty}-{1\over n}\log\Pr\{\text{event in (SO14)}\}
>2\log2.                                            \tag{SO15}
```

This is the exact balanced shell-entropy criterion.  To remove the leading
rectangular term one needs (SO15) with `r_n=o(n^(3/2))`.

### 5.1 What universal degree-two bounds provide

For every signing, `H_A(X)` is a centered Rademacher polynomial of degree two
with

```math
\|H_A\|_2^2={n\choose2}.
```

Bonami hypercontractivity gives, for every real `p>=2`,

```math
\|H_A\|_p\le(p-1)\sqrt{{n\choose2}},
```

and therefore

```math
\Pr\{|H_A|\ge t\}
\le\inf_{p\ge2}
 \left({(p-1)\sqrt{\binom n2}\over t}\right)^p.     \tag{SO16}
```

At `t=alpha n^(3/2)`, the optimizing moment has `p=Theta(sqrt(n))`,
so (SO16) is only `exp(-Theta(sqrt(n)))`.  It has zero rate on the left
side of (SO15).

The general Hanson--Wright inequality does not repair this.  The symmetric
sign matrix has Frobenius norm squared `Theta(n^2)` and, without additional
structure, operator norm as large as `Theta(n)`.  Hence

```math
\Pr\{|H_A|\ge t\}
\le2\exp\left[-c\min\left(
 {t^2\over n^2},{t\over n}
 \right)\right].                                   \tag{SO17}
```

At the project scale the second term again gives only
`exp(-Theta(sqrt(n)))`.  Near-optimality of the Boolean cap does not itself
give the missing `O(sqrt(n))` operator-norm bound, so (SO17) is the strongest
uniform conclusion available from these norms.

### 5.2 A tractable Hadamard bridge still does not close the entropy

Let `W` be a Hadamard matrix.  Conditional on `U`, the bridge energy is a
Rademacher sum whose coefficient vector has squared norm exactly `n^2`.
Thus the elementary moment generating function bound gives

```math
\Pr\{|U^{\mathsf T}WV|\ge t\}
\le2\exp\left(-{t^2\over2n^2}\right),               \tag{SO18}
```

while Cauchy--Schwarz gives `L_W<=n^(3/2)`.

This bridge histogram is exactly the binary Hadamard **glow** studied by
[Banica, *The glow of Fourier matrices: universality and fluctuations*,
Theorem 1.6](https://arxiv.org/abs/1403.2108).  Banica's random variable
`Omega` is

```math
\Omega=U^{\mathsf T}WV.
```

Passing from all sign pairs to projective pairs only divides every
multiplicity by the same factor, so `N_zeta` is precisely the magnitude-glow
histogram reflected about `L_W`:

```math
N_\zeta(z)=\#\{(U,V):|\Omega|=L_W-z\}.              \tag{SO18a}
```

Theorem 1.6 gives Gaussian limiting moments for `Omega/n` at every fixed
moment order.  This correctly describes the central `Omega(n)` scale.  It
does not control (SO14), where the relevant bridge values can be
`Theta(n^(3/2))`, equivalently `Omega/n=Theta(sqrt(n))`.  Recovering a
Gaussian large-deviation estimate there by the moment method requires moment
orders `p=Theta(n)`, while the theorem takes `n` to infinity with `p` fixed.
There is no uniformity statement at the required growing order.  The paper
also singles out the Walsh glow as a further conjectural computation, so its
fixed-moment Gaussianity cannot be promoted to the needed Walsh extreme-tail
enumerator.

It is useful to grant the children a bound stronger than the rigorous
universal estimates and see whether even that would suffice.  Suppose,
optimistically, that each child energy were sharply subgaussian with its
exact variance proxy `binom(n,2)`, for all moment parameters.  Combining
those two hypothetical child bounds with the rigorous bridge moment bound
and taking the four possible signs gives

```math
\Pr\{|H_A+H_D|+|U^{\mathsf T}WV|>T\}
\le4\exp\left(-{T^2\over4n^2}+o(n)\right).           \tag{SO19}
```

Write

```math
M_A=M_D=(c+o(1))n^{3/2},
\qquad r_n=(s+o(1))n^{3/2}.
```

For (SO19) to reach the entropy threshold (SO15), it would require

```math
2c+s>\sqrt{8\log2}=2.35482\ldots.                   \tag{SO20}
```

For the strongest known near-optimal structured scale `c<=1/2`, this asks
for

```math
s>1.35482\ldots.
```

But the entire Hadamard rectangular budget has `s<=1`.  Therefore even this
optimistic sharp-subgaussian calculation cannot certify a positive
leading-order `g=L_W-r_n`; the rigorous universal bounds are strictly weaker.
The miss is by a leading constant, not a logarithm or a lower-order term.

This does not falsify (SO14).  Near the exact cap, a structured landscape can
have much thinner shells than any variance bound detects.  It proves that a
successful use of the histogram theorem needs an exponentially sharp
near-cap shell enumerator specific to the chosen family, not
hypercontractivity, Hanson--Wright, or second moments.

### 5.3 Is the histogram state genuinely smaller?

Yes at the level of information: (SO10) depends only on polynomially many
marginal shell counts and is invariant under every permutation of states
within a shell.  The full bridge response records an exponentially large
alignment table, and the finite audit below shows that its best switching
gain can be strictly larger than the histogram-certified gain.

At theorem scale the price of this compression is now explicit.  The
marginal counts must be known to exponential accuracy at rate `2 log 2` near
their extreme shells.  Such a shell theorem is not logically equivalent to
cap optimization—it neither identifies a maximizing spin nor retains its
alignment with a bridge, and algebraic families can in principle have
explicit weight enumerators.  For an arbitrary near-optimal signing,
however, obtaining precisely those extreme-shell counts is at least a new
counting obligation not implied by near-optimality.  If one restores state
identities to avoid (SO15), one has returned to the full-response bridge
problem.

## 6. Finite exact audit

The reproducer evaluates (SO10) and exhausts every switching of each bridge.
The first three bridges were originally found by state-dependent CP-SAT, so
they are checks of the identity, not independent structured-family evidence.
The last two use fixed Sylvester bridges and are the relevant clean test.

| children | child caps | bridge | `L_W` | histogram-certified `g` | best switching `g` | certified parent cap |
|---:|---:|:---|---:|---:|---:|---:|
| `5+5` | `4,4` | CP-SAT witness | 13 | 2 | 8 | 19 |
| `6+6` | `5,5` | CP-SAT witness | 16 | 2 | 8 | 24 |
| `6+7` | `5,9` | CP-SAT witness | 16 | 2 | 10 | 28 |
| `4+4` | `4,4` | Sylvester | 8 | 4 | 4 | 12 |
| `8+8` | `10,10` | Sylvester | 20 | 4 | 8 | 36 |

The criterion therefore certifies a strict second-order gain, including for
a bridge chosen without the child response table.  It does **not** yet show a
gain proportional to `L_W`; the certified Sylvester gain remains four from
order four to order eight.  The larger exhaustive switching gains are finite
evidence only.

## 7. Research judgment and stopping boundary

The histogram theorem is a valid nonlocal compression mechanism and avoids
the earlier exact-linear response-rank obstruction: it preserves shell
counts, not every optimized response.  Its useful asymptotic target is
exactly (SO12) with a recurrence-scale residual.

No present estimate proves that target.  Replacing the histograms by only
their second moments gives ordinary subgaussian union bounds and restores a
leading rectangular cost.  Tracking the identity of the low-deficit pairs,
rather than proving the marginal convolution bound (SO12), is the full bridge
optimization problem and is outside this route.  The next defensible step is
therefore a shell-enumerator theorem for a specified structured child and
bridge family; without such a theorem, further exact alignment work would be
the stopped full-response route.

## Reproduction

```bash
.venv/bin/python \
  computations/audit_second_order_amalgamation_histograms.py \
  --output computations/results/second_order_amalgamation_histograms.json
```

The script uses exact integer enumeration throughout and verifies that the
histogram-certified gain never exceeds the exhaustively optimized switching
gain.
