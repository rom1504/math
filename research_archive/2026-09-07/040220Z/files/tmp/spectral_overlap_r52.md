# Wave 52A: exact core-load calculus and the small-core circularity boundary

## Status and research judgment

This attack does **not** prove an asymptotic overlap lemma or convergence.
It gives three exact, independently checked conclusions.

1. There is a simple core-load formula for every common-core retention and
   an exact pair-collision formula for `P_2-lambda_2`.
2. The formula exposes a decisive circularity boundary.  Level two, and more
   generally every core scale `ell=O(L)` for
   `L=n^(3/4-c)`, has no spectral amplification at the target exponent.  A
   saved retention estimate at such a scale already implies the saved center
   degree directly, before using the spectral-excess theorem.
3. A genuinely different spectral implementation must use a core with
   `-log rho_ell >> L` (in particular `ell >> L`, typically a linear-sized
   core) and prove saved *positive* excess above `lambda_2`.  No such theorem
   is proved here.  One-port selector exchange is finitely false even for
   exact minimizers and row-optimal box completions.

All algebra and finite calculations are reproduced by
`tmp/spectral_overlap_r52_check.py`; its saved output is
`tmp/spectral_overlap_r52.out`.

## 1. Exact core-load formula

Let `N=C(n,m)`, let `F_z` be the favorable selector family of a center `z`,
and put

```math
r_z=|F_z|,\qquad a_z=r_z/N,\qquad
c_z(R)=|\{S\in F_z:R\subseteq S\}|.
```

For an `ell`-set `R`, set

```math
\rho_\ell=\frac{\binom m\ell}{\binom n\ell},qquad
\mu_z(R)=\frac{c_z(R)}{r_z},qquad
u_z(R)=\frac{\mu_z(R)}{\rho_\ell}.
```

Then `E_R u_z(R)=1`.  Exact common-core counting gives the **Verified
identity**

```math
\boxed{
p_\ell(z)
=\frac{\sum_{|R|=\ell}c_z(R)^2}
       {r_z\binom m\ell\binom{n-\ell}{m-\ell}}
=a_z\,\mathbb E_R u_z(R)^2.}
\tag{R52A.1}
```

Since `0<=u_z(R)<=1/rho_ell`, Jensen and `u^2<=u/rho_ell`
give the sharp elementary bounds

```math
\boxed{a_z\le p_\ell(z)\le\frac{a_z}{\rho_\ell}.}
\tag{R52A.2}
```

For a center class `C`, define as in Wave 51

```math
P_\ell=\frac{\sum_{z\in C}a_z^2p_\ell(z)}
                 {\sum_{z\in C}a_z^2},qquad
r_C=\frac{\sum_{z\in C}a_z^3}{\sum_{z\in C}a_z^2},qquad
M=\max_{z\in C}a_z.
```

Averaging (R52A.2) proves the **Verified aggregate sandwich**

```math
\boxed{r_C\le P_\ell\le\frac{r_C}{\rho_\ell}
\le\frac M{\rho_\ell}.}
\tag{R52A.3}
```

Combining this with (10.1261), where

```math
D_\ell=(1-\lambda_1(\ell))
       +n(\lambda_1(\ell)-\lambda_2(\ell)),
```

gives the exact **Verified hybrid extraction**

```math
\boxed{
M\ge\max\left\{
\rho_\ell P_\ell,
\frac{[P_\ell-\lambda_2(\ell)]_+}{D_\ell}
\right\}.}
\tag{R52A.4}
```

The first term is direct collision counting; only the second is spectral
amplification.

## 2. Exact higher-core and level-two surpluses

Put

```math
d_\ell=\binom m\ell\binom{n-\ell}{m-\ell},qquad
b_\ell=\binom{m-2}{\ell-2}
       \binom{n-\ell-2}{m-\ell}.
```

Direct factorial cancellation gives

```math
\lambda_2(\ell)d_\ell=b_\ell.
```

Therefore the aggregate excess has the **Verified exact load form**

```math
\boxed{
P_\ell-\lambda_2(\ell)
=\frac{\displaystyle
 \sum_{z\in C}r_z
 \left\{\sum_{|R|=\ell}c_z(R)^2-b_\ell r_z\right\}}
 {\displaystyle d_\ell\sum_{z\in C}r_z^2}.}
\tag{R52A.5}
```

Equivalently, because

```math
\sum_{|R|=\ell}c_z(R)^2
=\binom m\ell r_z
 +2\sum_{|R|=\ell}\binom{c_z(R)}2,
```

the required positive term is collision between *distinct* favorable
selectors through a common core; the diagonal contribution is explicit.

At level two, (R52A.5) becomes the particularly concrete formula

```math
\boxed{
P_2-\lambda_2(2)
=\frac{\displaystyle
 \sum_{z\in C}r_z
 \left\{\sum_{i<j}c_z(ij)^2
             -\binom{n-4}{m-2}r_z\right\}}
 {\displaystyle
 \binom m2\binom{n-2}{m-2}\sum_{z\in C}r_z^2}.}
\tag{R52A.6}
```

Thus `P_2>=lambda_2(2)` is exactly a weighted pair-load collision
inequality.  But it already forces, without (10.1261),

```math
\boxed{
M\ge r_C\ge\rho_2\lambda_2(2)
=\frac{2(n-m)(n-m-1)}
       {n(n-1)(n-2)(n-3)}=\Theta(n^{-2})}
\tag{R52A.7}
```

on a fixed-density window.  Equality `P_2=lambda_2(2)` is allowed in
(R52A.7).  This is much stronger than the saved degree needed for
convergence, but it is not an overlap amplification: at fixed density
`rho_2` is a constant and (R52A.3) says `P_2` and `r_C` are equivalent up
to a constant factor.

The diagonal/self-loop must not be mistaken for signing evidence.  Every
nonempty family has

```math
p_2(z)\ge h_2:=\binom{n-2}{m-2}^{-1},qquad
\frac{h_2}{\lambda_2(2)}
=\frac{\binom m2}{\binom{n-4}{m-2}}.
\tag{R52A.8}
```

All stored examples with `m>n/2` and `n<=10` satisfy `h_2>=lambda_2(2)`.
Their nonnegative level-two excess is therefore contaminated by the
finite self-loop and is **not asymptotic evidence**.  At fixed density the
binomial in the denominator of (R52A.8) is exponential, so the inequality
reverses and the diagonal is negligible.

## 3. Precise circular and genuinely spectral regimes

Let `L=n^(3/4-c)`, with the same `0<c<1/4` as the restriction iteration.
Equation (R52A.4) gives an exact test:

- if `-log(rho_ell lambda_2(ell))=O(L)`, even the non-strict inequality
  `P_ell>=lambda_2(ell)` directly gives a saved center;
- if `-log rho_ell=O(L)`, any saved lower bound on `P_ell` (hence any saved
  positive spectral excess) directly gives a saved center via the first
  term of (R52A.4);
- spectral amplification can be genuinely useful only when
  `rho_ell lambda_2(ell)=exp{-omega(L)}` and one proves
  `P_ell-lambda_2(ell)>=exp{-O(L)}`.

Uniformly on a compact fixed-density window and for `ell=o(n)`,

```math
\log\rho_\ell=\ell\log(m/n)+O(\ell^2/n),
\qquad -\log\lambda_2(\ell)=O(\log n).
\tag{R52A.9}
```

Consequently every `ell=O(L)` lies in the direct/circular regime at the
project's exponential resolution.  A non-vocabulary-only spectral program
must use `ell >> L` (or, most cleanly, a linear-sized core), where
`rho_ell=exp{-omega(L)}`.  The refined **Open target** is therefore:

> At a controlled project-row cap, find a core scale with
> `-log rho_ell >> L` and prove
> `P_ell-lambda_2(ell)>=exp{-O(L)}`.

This is strictly narrower than the Wave 51 statement allowing arbitrary
`2<=ell<m`.  A level-two pair-load theorem could still prove convergence,
but (R52A.3) shows that it would simply be a direct saved-degree theorem in
equivalent notation.

## 4. Adaptive scales do not remove the missing mechanism

Let `H_z(j)` count ordered pairs `(S,T) in F_z^2` with
`|S\cap T|=j`.  Then (R52A.5) is the binomial-transform identity

```math
P_\ell-\lambda_2(\ell)
=\frac{\displaystyle
 \sum_z r_z\left\{
   \sum_jH_z(j)\binom j\ell-b_\ell r_z\right\}}
 {\displaystyle d_\ell\sum_zr_z^2}.
\tag{R52A.10}
```

This gives an exact adaptive-scale diagnostic but no positive average over
scales.  A singleton family has `p_ell=h_ell` and, at fixed density for all
sufficiently large `n`, `h_ell<=lambda_2(ell)` for every `2<=ell<m`.
Thus nonemptiness, cap nesting, and positive mixtures over core sizes cannot
force a winning scale.  Exact-minimizer structure must supply a genuinely
global overlap statement at a large core.

## 5. Selector exchange: exact finite walls and scale obstruction

For adjacent selectors `U+v` and `U+w`, the exact port-regret formula
(10.1211) says that favorability at `v` zeroes one orientation's three
nonnegative regrets at that port.  It puts no sign or size constraint on
the corresponding regrets at `w`.  Keeping the same center `z` while
exchanging the selector leaves `R_2(z)` *exactly unchanged*, so row
bookkeeping does not repair this missing regret transfer.

Two exact-minimizer checks make the obstruction concrete:

- in the order-ten minimizer at `m=6`, the two row-ten centers have the same
  five-selector family.  Distinct ordered intersections have histogram
  `{4:10, 2:10}` and never size `5`; hence **every one-port exchange out of
  each favorable selector is unfavorable**, despite the optimal box row;
- in the order-nine minimizer at `m=7`, the minimum active row class has
  family sizes `[1,3,1]`.  Two box-optimal centers have singleton favorable
  families.  Thus exact minimality plus a row-optimal witness does not even
  force a second favorable selector pointwise.

These are **finite falsifications of a pointwise or one-exchange theorem**,
not asymptotic counterexamples to aggregate large-core excess.  They also
show why an unchecked exchange iteration cannot be used.

There is a separate scalable count obstruction.  Even if one could make
every selector in a Johnson ball of radius `s=o(n)` favorable, that ball has
only

```math
\sum_{j\le s}\binom mj\binom{n-m}j=\exp\{o(n)\}
```

members.  Since `p_2(z)<=r_z h_2` and `h_2=exp{-Theta(n)}` at fixed density,
such local propagation still has `p_2=exp{-Theta(n)}<<lambda_2(2)`.
Level-two success requires global, linear-depth propagation or directly
the saved selector degree.  The exact large-core route instead needs a
global retention/excess theorem, not a bounded selector-exchange lemma.

## 6. What remains live

- **Verified:** (R52A.1)--(R52A.10), including the exact pair/higher-core
  load formulas, the hybrid extraction, and the small-core circularity
  boundary.
- **Falsified:** one-port or pointwise exchange from a low-row box witness;
  any claim that the stored small examples support asymptotic level-two
  excess (their threshold is forced by the self-loop).
- **Numerical only:** the displayed order-nine/order-ten family data.
- **Open:** minimizer-specific saved excess at a genuinely large core
  (`-log rho_ell >> L`), or a direct pair-load/degree theorem strong enough
  to bypass spectral amplification entirely.
