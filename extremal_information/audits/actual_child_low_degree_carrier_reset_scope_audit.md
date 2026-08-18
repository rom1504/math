# Reset-scope audit of the normed low-degree actual-child carrier

Status: **rigorous scope audit and coherent-secant theorem**.  The normed
low-degree carrier is not merely a relabelling of an arbitrary row-factor
table: it gives a uniformly value-complete, diffuse, polynomial-description
class of product competitors, and every extensive coherent-retuning phase has
one fixed-degree witness in this class.  In fact the witness can be organized
as one coherent product path with a linearly favorable directional score.

It is not yet a child-side branch *decision* statistic.  The restricted
optimum still pairs the carriers with the exact actual-child pressure, and no
proved state smaller than that response oracle evaluates or globally
optimizes those pairings.  Exact response closure is impossible: even normed
affine positive-part carriers span the full pressure table.  Therefore the
proper campaign classification is a **factor/certificate RESET**, not a
Level-6 closure and not a proof that the balanced product phase has been
decided.

## 1. Setup and the exact branch gain

Fix two actual contracted-temperature minimizing children, one orientation
and one row direction.  There are `m` bridge rows of width `n`, `N=m+n`, and
the physical amplitude is `u=beta/sqrt(N)`.  Let `L` be the exact bridge
pressure and define, for a row product `P=\otimes_i P_i`,

```math
\mathcal F(P)=\mathbb E_PL+{1\over\lambda}
 \sum_{i=1}^mD(P_i\Vert U_n).                     \tag{RS.1}
```

Let `r=\otimes_i r_i` be the canonical product and let `p^*` minimize
`mathcal F` over all row products.  The exact variational identities give

```math
\boxed{
\mathcal J-\mathcal I^{\leftarrow}
=\lambda\{\mathcal F(r)-\mathcal F(p^*)\}.}       \tag{RS.2}
```

For `d>=0` and `K>=1`, write `C_(n,d,K)` for the carrier densities

```math
q_g={g_+\over\mathbb E_Ug_+},\qquad
\deg g\le d,\qquad \mathbb E_Ug=1,\qquad
\|g\|_2\le K.                                    \tag{RS.3}
```

Put `K_0=exp(lambda^2 beta^2/2)`.  Theorem LDC.3 supplies a sequence
`eta_d=O_(beta,lambda)(d^(-1/8) log d)` tending to zero such that

```math
0\le \inf_{P_i\in\mathcal C_{n,d,K_0}}\mathcal F(P)
       -\mathcal F(p^*)
\le m\eta_d                                      \tag{RS.4}
```

uniformly in the orders, minimizing children, orientation, and row
direction.

## 2. One coherent finite-degree retuning path

The following makes precise the sense in which (RS.4) is stronger than a
collection of separately paid scalar channels.

**Theorem RS.1 (coherent finite-degree secant and tangent certificate).**
Fix `alpha,beta,lambda>0`.  There is a finite degree
`d=d(alpha,beta,lambda)`, independent of `m,n` and of the actual children,
with the following property.  Whenever

```math
\mathcal J-\mathcal I^{\leftarrow}\ge\alpha N,    \tag{RS.5}
```

there are row laws `q_i in C_(n,d,K_0)` such that, for
`q=\otimes_iq_i`,

```math
\boxed{
\lambda\{\mathcal F(r)-\mathcal F(q)\}
\ge {\alpha N\over2}.}                            \tag{RS.6}
```

Moreover, put

```math
p_i(s)=(1-s)r_i+s q_i,\qquad P_s=\bigotimes_i p_i(s),
\qquad 0\le s\le1.                                \tag{RS.7}
```

Then there is `s in (0,1)` at which the single coherent directional score
obeys

```math
\boxed{
-\sum_{i=1}^m\left[
 \lambda\,\mathbb E_{(q_i-r_i)\otimes P_{s,-i}}L
 +\mathbb E_{U_n}(q_i-r_i)\log p_i(s)
 \right]
\ge {\alpha N\over2}.}                            \tag{RS.8}
```

Here a density and its Radon--Nikodym derivative with respect to `U_n` are
identified.  Formula (RS.8) retains the joint conditional pressure before
taking the final scalar sum; it does not bound the row channels separately.

*Proof.* Choose `d` so that `lambda eta_d<=alpha/2`.  By (RS.4), some
carrier product `q` has

```math
\mathcal F(q)\le\mathcal F(p^*)+m\eta_d.
```

Equations (RS.2), `m<=N`, and the choice of `d` prove (RS.6).

Let `G(s)=\lambda\{\mathcal F(r)-\mathcal F(P_s)\}`.  The canonical factors are
strictly positive.  Hence, for `0<s<1`, every `p_i(s)` is strictly positive
and ordinary differentiation under the finite sums gives

```math
G'(s)=-\sum_i\left[
 \lambda\,\mathbb E_{(q_i-r_i)\otimes P_{s,-i}}L
 +\mathbb E_U(q_i-r_i)\{1+\log p_i(s)\}
 \right].                                        \tag{RS.9}
```

Both densities have mass one, so the terms containing `1` vanish.  The
function `G` is continuous on `[0,1]` and absolutely continuous there (the
finite-dimensional entropy `x log x` has an integrable one-sided derivative
along a segment).  Since `G(0)=0` and (RS.6) gives
`G(1)>=alpha N/2`, the integral of `G'` is at least `alpha N/2`.
Thus (RS.8) holds at some interior differentiability point. `square`

The endpoint `q` uses

```math
m\sum_{a=0}^d{n\choose a}=N^{O_{\alpha,\beta,\lambda}(1)}       \tag{RS.10}
```

bounded real coefficients.  Every factor is diffuse (here `q_i(E)` denotes
probability mass):

```math
q_i(E)\le K_0\sqrt{U_n(E)},\qquad
\max_b q_i(\{b\})\le {K_0\sqrt{\sum_{a=0}^d{n\choose a}}\over2^n}.
                                                               \tag{RS.11}
```

Thus (RS.6)--(RS.8) cannot hide a deterministic bridge word or an
exponential row table in one coefficient.  A fixed-accuracy quantization
still has polynomial description length.  Given the coefficient list,
point-evaluation access to the actual pressure, and a separately certified
value (or sampling oracle) for the common endpoint `mathcal F(r)`,
Proposition LDC.4 gives randomized one-sided verification of (RS.6) from
carrier samples.  This explicit endpoint qualification matters: the carrier
theorem does not itself compress the canonical response.  Subject to it,
(RS.6) is a genuine lower-information *certificate* for alternative (iii).

## 3. Why the restricted optimum is not yet a branch observable

Define the restricted reverse projection

```math
\mathcal I_{d,K_0}^{\leftarrow}
=\lambda\left\{
 \inf_{P_i\in\mathcal C_{n,d,K_0}}\mathcal F(P)-V_\lambda
 \right\}.                                       \tag{RS.12}
```

The proved estimate

```math
0\le\mathcal I_{d,K_0}^{\leftarrow}
       -\mathcal I^{\leftarrow}
\le\lambda m\eta_d                               \tag{RS.13}
```

means that, *as a variational value*, one fixed degree distinguishes a
fixed extensive gap.  It does not supply a two-sided operational decision:

1. a low value has a carrier witness, but proving that the global restricted
   minimum is high still needs a lower-bound mechanism;
2. evaluating the objective couples all rows through the exact pressure
   `L`; the coefficient list compresses the factors, not this query;
3. optimizing a polynomial-dimensional black-box nonconvex objective is not
   made easy merely by naming its minimum.

There is also an important information-accounting caveat.  The two child
sign matrices themselves contain only `O(N^2)` bits and determine every
pressure value by (possibly exponential) computation.  For `d>1`, the
formal carrier coefficient count in (RS.10) can even exceed `N^2`.
Therefore polynomial coefficient dimension is not, by itself, an absolute
information lower bound or an algorithmic improvement over the raw child
input.  It is a rigorous compression relative to the exponentially long
optimal-factor tables.  Turning that relative compression into an
optimizer-specific child observable requires a query, precision, or
computational closure theorem.

The exact obstruction is substantive.  Theorem AB.1 shows that, for every
fixed `K>1` and large enough row width, normed degree-one positive-part
carriers already contain a basis of all row functions.  Tensor products of
those carriers therefore make the complete exact response surface
table-complete.  The high-degree eigenvalues are exponentially attenuated,
so this does not refute macroscopic approximation of the optimized scalar;
it does refute any claim that the carrier's presenting Walsh degree alone
gives an exact child-response closure.

Consequently `I_(d,K_0)^leftarrow` is not yet the requested
optimizer-specific low-information *observable*.  It is a restricted
optimization problem with a genuinely smaller competitor representation.
Calling it a completed branch decider would simply move the missing response
oracle behind a scalar symbol.

## 4. RESET/STRIKE judgment

The campaign result passes four nontrivial tests.

- **Not a renamed factor table:** fixed degree and norm give polynomial
  description and quantitative diffuseness, uniformly for actual optimizing
  children.
- **Not independently paid channels:** the carrier coefficients enter one
  product expectation, and Theorem RS.1 yields one coherent product path.
- **Macroscopic completeness:** every fixed extensive branch-(iii) gain has
  a fixed-degree witness with only a fixed-factor loss.
- **Strictly narrower remaining obligation:** arbitrary optimal row tables
  are no longer needed on the competitor side.

It fails two stronger tests.

- It does not decide alternatives (i)--(iii) from a proved child statistic.
- It does not close, approximate, or lower-bound the carrier response using
  a state smaller than the actual pressure landscape.

Accordingly the correct classification is:

> **RESET, narrowly, by a materially weaker SML and a coherent
> lower-information certificate; not a branch decision and not Level 6.**

The new smallest missing lemma is not “find the full product shadow.”  It is:

> **Macroscopic finite-degree child closure.**  For one fixed `d`, evaluate,
> upper/lower bound, or generate a linearly favorable normed-carrier response
> from an actual-child state strictly smaller than the complete pressure
> response surface, to `o(N)` or fixed-density accuracy; alternatively give
> a scalable actual-minimizer obstruction to every such robust closure.

The exact affine-basis result shows why “retain the degree-`d` Walsh
coefficients of `L`” is not that lemma.  The required result must be about
the optimized value at macroscopic precision, not exact algebraic response.
