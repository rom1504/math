# Adversarial audit of the low-degree density carrier

Status: **PASS for the one-shot variational approximation, subject to an
essential normalization constraint; FAIL as a child-only closure theorem.**

The proposed construction starts from a row density `f=dP/dU_n`, projects
it to Walsh degree at most `d`, and repairs positivity:

```math
g_d=\Pi_{\le d}f,
\qquad
\widetilde f_d={ (g_d)_+\over E_U(g_d)_+}.          \tag{LDA.1}
```

Under the actual-factor hypotheses, this does approximate the complete
row-product variational objective by `o_d(1)N`, uniformly in the child
orders.  The proof is sound and is recorded below because several plausible
failure modes do not occur.  However, the coefficient list describes only
the competing **factors**.  It does not give a finite response state for the
actual collision--cavity interaction.  In particular, positivity clipping
restores high Walsh degrees, and evaluating or optimizing the product
expectation can still require a large cross-row response surface.

## 1. The uniform Fourier estimate is valid

Let `f` be a probability density on the `n`-cube and suppose

```math
\|f\|_2\le K,
\qquad
|\log f(x)-\log f(x^{(a)})|\le {A\over\sqrt n}     \tag{LDA.2}
```

for every one-bit flip.  Put `delta=A/sqrt(n)`.  Pointwise,

```math
|f(x)-f(x^{(a)})|
\le(e^\delta-1)f(x).
```

The Walsh Dirichlet identity consequently gives

```math
\sum_S |S|\widehat f(S)^2
={1\over4}\sum_aE_U(f-f^{(a)})^2
\le {K^2A^2e^{2A}\over4}.                          \tag{LDA.3}
```

Thus, with

```math
\tau_d={KAe^A\over2\sqrt{d+1}},
```

one has

```math
\|f-g_d\|_2\le\tau_d.                              \tag{LDA.4}
```

This estimate is uniform in `n`.  Bounded `D_2(P||U)` is exactly the
uniform `L^2` bound needed here; the weak-coordinate assumption supplies
the bounded Dirichlet energy.

## 2. Positive-part repair does not destroy the approximation

Since `E_Ug_d=E_Uf=1`, if `a_d=E_U(g_d)_-`, then

```math
Z_d:=E_U(g_d)_+=1+a_d\ge1.                         \tag{LDA.5}
```

On the set where `g_d<0`, positivity of `f` gives
`(g_d)_-<=|f-g_d|`.  Hence

```math
\|(g_d)_-\|_2\le\tau_d,
\qquad a_d\le\tau_d.
```

Using

```math
{(g_d)_+\over Z_d}-g_d
={(g_d)_--a_dg_d\over1+a_d}
```

and `||g_d||_2<=K` yields

```math
\boxed{
\|\widetilde f_d-f\|_2
\le(2+K)\tau_d=:\varepsilon_d,
\qquad
\|\widetilde f_d\|_2\le K.}                      \tag{LDA.6}
```

Therefore clipping is harmless for this one-shot `L^2` approximation.
It is not harmless algebraically: `(g_d)_+` generally has Walsh coefficients
of all degrees.

## 3. The full product-objective telescope is dimension-free

Let `P=\bigotimes_{i=1}^mP_i` and replace its factors successively by the
laws `\widetilde P_i` from (LDA.1).  Suppose a bridge interaction `h` has
one-bit oscillation at most `B/sqrt(N)`, with `m,n<=N`.  At the `i`th
telescoping step, average over all other rows and call the resulting row
function `H_i`.  It has the same bit oscillation.  Efron--Stein under the
fair row law gives

```math
\operatorname {Var}_{U_n}H_i\le {B^2n\over4N}\le {B^2\over4}. \tag{LDA.7}
```

The two row densities have equal mass, so Cauchy--Schwarz and (LDA.6)
give

```math
|E_{P_i}H_i-E_{\widetilde P_i}H_i|
\le {B\over2}\varepsilon_d.                       \tag{LDA.8}
```

After all rows,

```math
|E_Ph-E_{\widetilde P}h|
\le {B\over2}m\varepsilon_d.                       \tag{LDA.9}
```

There is no hidden `sqrt(N)` loss.  The relevant comparison is in
`L^2(U_n)` against the variance of the conditional row query, rather than
a total-variation comparison times its full range.

The entropy term is also stable.  On any probability space, the map

```math
f\longmapsto E[f\log f]
```

is uniformly continuous on an `L^2` ball.  Explicitly, truncating at
`R=\varepsilon^{-1/2}` and at `\varepsilon` gives a modulus

```math
|E[f\log f]-E[q\log q]|
\le\omega_K(\|f-q\|_2),
\qquad
\omega_K(s)=O_K(\sqrt s\log(e/s)).                 \tag{LDA.10}
```

The precise modulus is unimportant here; it tends to zero independently
of the cube dimension.  If a reference row density `r_i` has log-density
bit oscillation at most `C/sqrt(N)`, then, after subtracting its fair mean,
Efron--Stein similarly gives

```math
|E_{P_i}\log r_i-E_{\widetilde P_i}\log r_i|
\le {C\over2}\varepsilon_d.                       \tag{LDA.11}
```

Consequently, for

```math
\mathcal G(P)=E_Ph+{1\over\lambda}
              \sum_iD(P_i\Vert r_i),
```

```math
\boxed{
|\mathcal G(P)-\mathcal G(\widetilde P)|
\le m\left[
 {B\over2}\varepsilon_d
 +{1\over\lambda}
   \left\{\omega_K(\varepsilon_d)
          +{C\over2}\varepsilon_d\right\}
 \right]
=o_d(1)N.}                                        \tag{LDA.12}
```

This proves the advertised analytic statement for comparable splits (and
with the evident aspect-ratio constants otherwise).

## 4. Applicability to the actual optimizing factors

The hypotheses are available for the actual globally optimal row-product
shadow.  The coordinate best-response equation AC.17 gives

```math
D_2(p_i^*\Vert U_n)\le\lambda^2\beta^2{n\over N}
```

and one-bit log-density oscillation `2lambda beta/sqrt(N)`.  The canonical
reference factor has the same two bounds.  The actual collision--cavity
interaction `h` has one-bit oscillation at most `4beta/sqrt(N)`.  Thus all
constants in (LDA.12) depend only on the fixed thermodynamic parameters and
the split window.

Define a repaired carrier class by

```math
\mathcal A_{d,K}
=\left\{
 {g_+\over E_Ug_+}:
 \deg g\le d, E_Ug=1, \|g\|_2\le K
 \right\}.                                        \tag{LDA.13}
```

If `G_d` is the product variational minimum restricted to factors in
`A_(d,K)`, while `G_*` is the unrestricted product minimum, then (LDA.12)
proves

```math
\boxed{0\le G_d-G_*\le o_d(1)N.}                  \tag{LDA.14}
```

This is stronger than merely extracting a low-degree separator from the
unknown optimizer: it gives a genuine one-shot variational restriction.

## 5. Essential falsifier: norm control without mean control is circular

The conditions `E_Ug=1` and hence `E_Ug_+>=1` in (LDA.13) are essential.
A fixed polynomial norm by itself does **not** prevent point queries.
Fix a vertex `v` and put

```math
g_v(x)={1\over n}\{
       \langle v,x\rangle-(n-1)\}.
```

This affine polynomial has uniformly bounded `L^2(U_n)` norm.  It is
positive only at `x=v`, so

```math
{(g_v)_+\over E_U(g_v)_+}=2^n\mathbf1_{\{v\}},     \tag{LDA.15}
```

the point-mass density.  Thus an unrestricted positive-part affine carrier
already reconstructs every row point mass.

The inherited unit mean repairs this defect.  For every member `q` of
`A_(d,K)`,

```math
\|q\|_2\le K,
\qquad
U_n(\operatorname {supp}q)\ge K^{-2}.             \tag{LDA.16}
```

Moreover, the evaluation norm of the degree-`d` Walsh space gives

```math
2^{-n}q(x)
\le 2^{-n}K
 \sqrt{\sum_{j=0}^d{n\choose j}}.                 \tag{LDA.17}
```

For fixed `d`, no carrier law can approximate a point mass.  Quantizing
the coefficient ball produces a finite polynomial-dimensional net; the
positive-part normalization is Lipschitz on the unit-mean slice because
its denominator is at least one.

## 6. Why this is not yet a child-only response carrier

Equation (LDA.14) compresses the decision variables but not automatically
the query being optimized.  Three distinctions are mandatory.

1. The positive part of a degree-`d` polynomial generally has full Walsh
   degree.  Hence its expectation against `h` is not determined merely by
   degree-`d` Walsh coefficients of `h`.
2. Even without clipping, products of row-degree-`d` densities query all
   cross-row tensors with row degree at most `d`.  There are

   ```math
   \left(\sum_{j=0}^d{n\choose j}\right)^m
   ```

   such formal response coefficients.  A polynomial number of factor
   parameters must not be confused with a polynomial-size representation
   of the objective surface.
3. The repaired density can have zeros and therefore need not retain any
   finite log-oscillation bound.  The construction is a one-shot recovery
   class, not a carrier known to close under another interaction step.

Thus (LDA.14) is a strict **factor-side** information reduction and a valid
new representation of the retuning branch.  It is not by itself an
optimizer-specific child observable deciding that branch, a summable
almost-subadditivity error, or a reusable Level-6 state.  To obtain such a
result one still needs a bounded-degree child-response closure theorem:
evaluate the restricted product objective from an actual-child state that
is itself smaller than the collision--cavity landscape, or prove that the
actual interaction synchronizes the exponentially many cross-row tensors.

## 7. Audit judgment

- **Fourier-tail falsifier:** none; (LDA.3)--(LDA.4) are valid.
- **Positivity/entropy falsifier:** none after enforcing `E g=1`; the
  dimension-free bounds (LDA.6) and (LDA.10) suffice.
- **Telescoping-energy falsifier:** none; conditional-row variance removes
  the apparent extra `sqrt(N)` loss.
- **Real-bit falsifier:** point masses occur if only polynomial norm is
  imposed; the unit-mean slice plus fixed `L^2` bound repairs this exactly.
- **Actual applicability:** pass for the optimal factors and the one-shot
  product objective.
- **SML effect:** a real but partial weakening.  Full row-factor tables can
  be replaced, at `o(N)` objective loss, by a slowly growing bounded-degree
  recovery class.  The actual-child response needed to optimize that class
  remains uncompressed, so claiming a child-only branch decision or Level 6
  from this theorem alone would be circular.

