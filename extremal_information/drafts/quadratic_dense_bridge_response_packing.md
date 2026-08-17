# Quadratic children retain extensive information across a dense sign bridge

Status: rigorous theorem derived from the audited RD.1/RD.2 random-bridge
events.  The lower bound uses an explicit complete sign-quadratic subclass,
not programmable lookup landscapes.  A finite exact verifier accompanies the
proof.

## 1. Question and answer

For a sign bridge `B in {-1,1}^{n times n}` and a quadratic Boolean child

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad x\in\{-1,1\}^n,                                     \tag{QD.1}
```

define its response

```math
(P_BH_A)(y)=\max_x\{H_A(x)+x^TBy\}.                         \tag{QD.2}
```

The question is whether restricting the arbitrary RD.1 child landscapes to
`|a_ij|<=1`, or even to `a_ij in {-1,1}`, destroys every scalable packing at
error `epsilon n^(3/2)`.

It does not.  There is an explicit family of complete `+-1` quadratic forms
whose response functions contain `Omega(n)` bits at that scale, even modulo
additive constants.  On the other hand, the quadratic restriction does give a
large strict quotient compared with unrestricted RD.1 landscapes: a sign
quadratic has only `Theta(n^2)` input bits, and bounded real coefficients have
an `O(n^2 log(n/epsilon))`-bit coefficient quantization at the target error.

Thus the proved rate bracket, after the independent compression audit in
`quadratic_dense_bridge_compression_ceiling.md`, is

```math
\boxed{\Omega(n)\ \text{bits}
\ \le\ R_{\rm quadratic}(\epsilon n^{3/2})
\ \le\ O(n^2\log(1+1/\epsilon))\ \text{bits}}.             \tag{QD.3}
```

For the sign subclass the upper bound improves to `n(n-1)/2` exact bits.  The
gap between linear and quadratic bits remains open.

## 2. Contextual response metric

For an arbitrary future landscape `G(y)`, gluing returns

```math
\operatorname{Opt}_B(A,G)
=\max_y\{(P_BH_A)(y)+G(y)\}.                                \tag{QD.4}
```

As in the audited bridge-query isometry, coordinate-pinning futures give

```math
\sup_G|\operatorname{Opt}_B(A,G)-\operatorname{Opt}_B(A',G)|
=\|P_BH_A-P_BH_{A'}\|_\infty.                              \tag{QD.5}
```

Modulo a global score offset, the metric is

```math
d_{\rm proj}(A,A')
=\frac12\operatorname{osc}(P_BH_A-P_BH_{A'}).               \tag{QD.6}
```

Consequently a response packing is an intrinsic information lower bound, not
a failure of one proposed algorithm.

## 3. An exact pole-locking lemma

For `z in {-1,1}^n`, define the complete gauge-ferromagnetic quadratic

```math
H_z(x)=\sum_{i<j}z_iz_jx_ix_j
=\frac{(x^Tz)^2-n}{2}.                                      \tag{QD.7}
```

Every coefficient is exactly `+-1`; the only parameter is the planted pole
pair `{z,-z}`.

### Lemma QD.1 (pole locking under a sublinear coordinate field)

If `h in R^n` satisfies

```math
\|h\|_\infty<n/2,                                           \tag{QD.8}
```

then

```math
\max_x\{H_z(x)+h^Tx\}
=\frac{n^2-n}{2}+|h^Tz|,                                    \tag{QD.9}
```

and every optimizer is one of the two poles `+-z`.

#### Proof

For an arbitrary `x`, choose `s in {-1,1}` so that

```math
r=d_H(x,sz)\le n/2.                                         \tag{QD.10}
```

The quadratic loss relative to that pole is

```math
H_z(sz)-H_z(x)
={n^2-(n-2r)^2\over2}
=2r(n-r)\ge nr.                                             \tag{QD.11}
```

The field can favor `x` over `sz` by at most

```math
h^Tx-h^T(sz)\le2r\|h\|_\infty<nr.                          \tag{QD.12}
```

Thus every nonpole loses strictly to a pole.  Maximizing between `z` and
`-z` gives the absolute value in (QD.9). `square`

This lemma is the mechanism missing from a direct attempt to replace the
programmable RD.1 landscapes by arbitrary quadratics: it proves exposure for
the whole Boolean cube, not merely against a selected list of competitors.

## 4. Exponential packing of sign-quadratic children

### Theorem QD.2 (dense-bridge packing for `+-1` quadratics)

There are universal constants `gamma,g>0` such that, for all sufficiently
large `n`, there exist

```math
B_n\in\{-1,1\}^{n\times n},
\qquad N\ge\exp(\gamma n),                                  \tag{QD.13}
```

queries `y_1,...,y_N in {-1,1}^n`, and pairwise distinct complete sign
quadratics `H_1,...,H_N` of the form (QD.7), such that for every `c!=d`,

```math
\|P_BH_c-P_BH_d\|_\infty\ge g n^{3/2},
\qquad
d_{\rm proj}(H_c,H_d)\ge g n^{3/2}.                         \tag{QD.14}
```

The bridge may additionally be chosen with
`||B_n||_(2->2)=O(sqrt(n))`.

#### Proof

Use the audited RD.1 random experiment.  Fix constants
`C_0,d_0,c_op,c_diag>0` for which

```math
\Pr\{\|B\|_{2\to2}>C_0\sqrt n\}\le2e^{-c_{op}n},
\qquad
\Pr\{\|By\|_1<d_0n^{3/2}\}\le2e^{-c_{diag}n}.               \tag{QD.15a}
```

Take `B` with independent Rademacher entries and independent uniform sign
queries `y_c`.  Put

```math
h_c=By_c,\qquad z_c=\operatorname{sign}(h_c),\qquad
D_c=\|h_c\|_1=z_c^Th_c.                                    \tag{QD.15}
```

Set

```math
d_1=d_0/2,
\qquad c_\times={d_0^2\over8C_0^2},
\qquad
0<\gamma<{1\over8}
\min\{c_{op},c_{diag},1/8,c_\times\}.                       \tag{QD.15b}
```

Take `N=floor(exp(2 gamma n))`, which is at least `exp(gamma n)` for all
sufficiently large `n`.  The following events hold simultaneously with
positive probability:

```math
\begin{aligned}
\|B\|_{2\to2}&\le C_0\sqrt n,\\
D_c&\ge d_0n^{3/2} &&(c\le N),\\
\|h_c\|_\infty&<n/2 &&(c\le N),\\
|z_d^Th_c|&\le d_1n^{3/2} &&(c\ne d).
\end{aligned}                                               \tag{QD.16}
```

Here is the complete union-bound audit.

1. The first two lines are exactly the operator and diagonal estimates used
   in RD.1.  Each diagonal failure has probability `exp(-Omega(n))`.
2. For fixed `(i,c)`, `(By_c)_i` is a length-`n` Rademacher sum, so

   ```math
   \Pr\{|(By_c)_i|\ge n/2\}\le2e^{-n/8}.                    \tag{QD.17}
   ```

   A union bound over `nN` coordinates proves the third line.
3. Conditional on `B,y_d`, the vector `z_d` is fixed and `y_c` remains an
   independent sign vector.  On the operator event,

   ```math
   \|B^Tz_d\|_2\le C_0n.                                    \tag{QD.18}
   ```

   The two-sided Rademacher tail therefore gives

   ```math
   \Pr\{|z_d^TBy_c|>d_1n^{3/2}\mid B,y_d\}
   \le2\exp\{-d_1^2n/(2C_0^2)\}
   =2e^{-c_\times n}.                                       \tag{QD.19}
   ```

   Union-bound over the fewer than `N^2` ordered pairs.

With (QD.15b), the total failure probability is at most

```math
2e^{-c_{op}n}+2Ne^{-c_{diag}n}
+2nNe^{-n/8}+2N^2e^{-c_\times n}<1                          \tag{QD.19a}
```

for all sufficiently large `n`.  Notice that (QD.16) also forces
the `z_c` to be distinct even modulo sign: otherwise an off-diagonal absolute
inner product would equal `D_c`.

Fix a realization satisfying (QD.16), and let `H_d=H_{z_d}`.  Lemma QD.1
applied to the field `h_c` gives the exact response matrix

```math
(P_BH_d)(y_c)
=\frac{n^2-n}{2}+|z_d^Th_c|.                                \tag{QD.20}
```

At the diagonal coordinate `y_c`,

```math
(P_BH_c)(y_c)-(P_BH_d)(y_c)
\ge {d_0\over2}n^{3/2}.                                     \tag{QD.21}
```

At `y_d` the same response difference is at most
`-(d_0/2)n^(3/2)`.  Thus its sup norm is at least this amount and its
oscillation is at least twice this amount.  Set `g=d_0/2` and use
(QD.5)--(QD.6). `square`

### Corollary QD.2a (extensive semantic information)

For every fixed `epsilon<g/2=d_0/4`, a summary answering all future
continuations to absolute error at most `epsilon n^(3/2)` on this quadratic
subclass needs at least

```math
N\ge\exp(\gamma n)                                          \tag{QD.22}
```

states, or `Omega(n)` bits.  The same conclusion holds projectively, even if
the decoder may choose an arbitrary additive calibration for each response.
Indeed, if two children shared one decoded state, the projective triangle
inequality would put their distance at most `2 epsilon n^(3/2)`, contradicting
(QD.14).  The bit count is at least `(gamma/log(2))n-O(1)`.

The common quadratic maximum `(n^2-n)/2` is irrelevant: it cancels between
children and may be removed as a projective baseline.  The separated signal
in (QD.21) is genuinely at the requested `n^(3/2)` scale.

Since `a_ij in {-1,1}`, Theorem QD.2 simultaneously answers the weaker
`|a_ij|<=1` question.

## 5. The upper direction: quadratic children are still a strict quotient

The lower bound is extensive, but it does not reproduce RD.1's
`exp(Omega(n))` **bits**.  Such a conclusion is impossible for sign
quadratics: there are only

```math
p={n(n-1)\over2}                                             \tag{QD.23}
```

coefficients, so storing their signs is an exact `p`-bit composable state.

For bounded real coefficients there is a uniform approximate analogue.  The
coefficient grid below is the simplest deterministic construction; the
independent unbiased-rounding theorem QC.2 improves its bit count to
`O(n^2 log(1+1/epsilon))` while preserving the same uniform response error.

### Proposition QD.3 (coefficient-grid upper bound)

Suppose `|a_ij|<=1`, and fix a target error `epsilon n^(3/2)`.  Round every
coefficient to a grid of mesh

```math
\Delta={2\epsilon n^{3/2}\over p}.                           \tag{QD.24}
```

Then, for every spin assignment,

```math
|H_A(x)-H_{\widehat A}(x)|
\le {p\Delta\over2}=\epsilon n^{3/2}.                       \tag{QD.25}
```

Consequently

```math
\|P_BH_A-P_BH_{\widehat A}\|_\infty
\le\epsilon n^{3/2}                                        \tag{QD.26}
```

for every bridge `B`, and the same bound survives every exact future.  The
code uses at most

```math
p\log_2\!\left(2+\left\lceil {p\over\epsilon n^{3/2}}\right\rceil\right)
=O\!\left(n^2\log(n/\epsilon)\right)                         \tag{QD.27}
```

bits.

#### Proof

Each of the `p` summands changes by at most `Delta/2`, proving (QD.25).
Taking a maximum is sup-norm nonexpansive, first over `x` and then over the
future `y`, proving (QD.26).  The interval `[-1,1]` has at most
`2+ceil(2/Delta)` grid values, which is the count in (QD.27) up to a harmless
factor of two inside the logarithm. `square`

This is a true response approximation, not a convergence statement.  The
quantized quadratic is an exact replacement child, and its error is uniform
over all bridge queries and all later landscapes.

## 6. What both attacks establish, and what remains open

The two directions give a sharp qualitative answer.  Together with the
Hamming-cover theorem QC.1 and unbiased-rounding theorem QC.2, they give the
strongest current quantitative bracket.

1. **No subextensive-in-`n` response state.**  Complete `+-1` quadratic
   children already force `Omega(n)` bits at the target scale.
2. **A strict quotient relative to arbitrary children.**  Quadratic
   coefficient data give polynomially many bits, versus the
   `exp(Omega(n))` bits required by the programmable RD.1 landscapes.
3. **The optimal quadratic rate is unresolved.**  The present bounds leave
   `Omega(n)` versus `O(n^2)` bits for signs and for bounded real
   coefficients at every fixed positive error.

### Scope: arbitrary sign quadratics versus bounded-cap near-minimizers

Theorem QD.2 concerns the full class of sign-quadratic children and an
explicit highly ordered subclass inside it.  It does **not** prove the same
packing for the bounded-cap or near-minimizer children that motivate the
original signing problem.  Indeed,

```math
\max_xH_z(x)={n^2-n\over2},
\qquad
\operatorname{osc}(H_z)=\Theta(n^2).                        \tag{QD.28}
```

Thus these planted-pole children have quadratic-scale internal caps and
spreads.  They lie outside any hypothesis requiring a cap, or a projective
energy spread, of order `O(n^(3/2))`; they are also far from signings that
minimize such a cap.  Subtracting their common maximum as a response baseline
does not change the `Theta(n^2)` internal spread that supplies pole locking.

The exact conclusion is therefore:

- arbitrary complete `+-1` quadratic children do not admit an `o(n)`-bit
  uniform response quotient at `epsilon n^(3/2)` error;
- the smaller bounded-cap/near-minimizer class remains open and may still
  possess additional rigidity.

Several tempting routes do not close this gap.

- Interpolating independent RD.1 bonuses on `exp(Omega(n))` exposed states
  asks a quadratic polynomial with only `Theta(n^2)` coefficients to satisfy
  exponentially many independent conditions.
- Encoding RD.2's linear fields with a single anchor spin requires quadratic
  coefficients of their natural `Theta(sqrt n)` size; rescaling them into
  `[-1,1]` reduces the obvious separation below `n^(3/2)`.
- Random sign quadratics have the correct `n^(3/2)` energy scale, but their
  optimizers move with both the child and the query.  Finite random separation
  is not, by itself, a uniform exposure proof.  Lemma QD.1 succeeds because
  the entire optimizer set is rigorously reduced to two poles.

An `Omega(n^2)` lower bound would require a new family in which quadratically
many coefficient choices remain observable through the same `2^n` bridge
queries with target-scale margins.  No such exposure mechanism is proved
here.

## 7. Exact finite verification

Run

```bash
python3 extremal_information/experiments/verify_quadratic_dense_bridge_response.py
```

The verifier uses integer and rational arithmetic to check:

- Lemma QD.1 against exhaustive Boolean maximization for small `n`;
- a seeded finite dense-sign certificate and its absolute/projective response
  gaps using the exact formula (QD.20);
- that every child coefficient in the witness is `+-1`; and
- the coefficient-rounding and arbitrary-future bounds in QD.3.

The finite search distinguishes the precise pole-locking claim; it is not
used in the asymptotic probabilistic proof.

## 8. Verdict

**Scalable counterexample found within an explicit quadratic subclass.**
Dense sign bridges retain extensive `n^(3/2)`-scale response information even
for complete `+-1` quadratic children.  Quadratic restriction nevertheless
eliminates the doubly exponential state demand of arbitrary programmable
landscapes.  The defensible conclusion is the rate bracket (QD.3), not either
extreme conjecture.
