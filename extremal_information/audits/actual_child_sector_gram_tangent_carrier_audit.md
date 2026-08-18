# Audit of SQ.3: the sector--Gram tangent carrier

## Verdict

**PASS for the boxed carrier formula (SQ.23), with required corrections to
the displayed proof and scope.**  The sector--Gram state determines the
exact first nonzero interaction coefficient at fixed child orders.  There
is no missing factor of two in (SQ.23), and EO.2 gives the two separate
`O(u^6)` expansions claimed at tangent scale.  As currently written,
however:

1. the index ranges in (SQ.26) and the first factor of (SQ.28) are wrong;
2. the chain `T_u=J_u=...` should be written as two asymptotic equalities,
   since `J_u-T_u=O(u^8)` rather than zero identically;
3. EO.2 proves `O(u^8)` for the **canonical marginal-retuning term**, not
   for a globally best row-product retuning;
4. the state is a constant-dimensional exact-real response carrier, but no
   finite-precision information bound or repeated-composition closure has
   yet been proved.

These corrections do not change (SQ.23) or its genuine conclusion.

## 1. Normalization in (SQ.23)

For a child `C`, write

```math
C_C^a(r,s)=E_{\mu_{C,a,t}}X_rX_s,
\qquad
G_C(a,b)=\sum_{r<s}C_C^a(r,s)C_C^b(r,s).
```

For the zero-bridge sector `a`, the children are independent and

```math
\Gamma_{ik;j\ell}^{\epsilon}
=\sum_{a=\pm1}\pi_a^\epsilon
 C_A^a(i,k)C_D^{\epsilon a}(j,\ell),
```

with exactly the weights in (SQ.24).  EO.2 defines

```math
K_\epsilon
=\sum_{i<k}\sum_{j,\ell}
 (\Gamma_{ik;j\ell}^{\epsilon})^2.                 \tag{A.1}
```

Expanding (A.1) gives

```math
K_\epsilon
=\sum_{a,b}\pi_a^\epsilon\pi_b^\epsilon
 \left(\sum_{i<k}C_A^a(i,k)C_A^b(i,k)\right)
 \left(\sum_{j,\ell}C_D^{\epsilon a}(j,\ell)
                       C_D^{\epsilon b}(j,\ell)\right).
```

The first parenthesis is `G_A(a,b)`.  In the second, the `n` diagonal
terms are one and every unordered off-diagonal pair occurs twice, so it is

```math
n+2G_D(\epsilon a,\epsilon b).
```

This proves (SQ.23) exactly.  In particular, the asymmetric factors are
correct: there is no `m+2G_A` because EO.2 sums over `i<k`, whereas its
right indices `(j,ell)` are ordered and may coincide.

As a normalization check, when `D` has order two,
`G_D(a,b)=ab\tanh^2t`.  Formula (SQ.23) becomes

```math
2\left\|\sum_a\pi_a v_A^a\right\|^2
+2\tanh^2t\left\|\sum_a a\pi_a v_A^a\right\|^2,
```

which is exactly the `2 sum_(i<k)(a_ik^2+b_ik^2)` formula in EO.4.

The displayed (SQ.26) must therefore read `sum_(i<k) sum_(j,l)`, not
`sum_(i,k=1)^m sum_(j,l=1)^n`.  Likewise the first parenthesis in (SQ.28)
must be `sum_(i<k)`.  The prose immediately after (SQ.29) uses the correct
ranges, so these are display errors rather than an error in the boxed
formula.

## 2. What state is sufficient, and how large is it?

The defined state is sufficient for `K_epsilon`: it contains the child
orders, the sector weights needed to form `pi`, and every entry of the two
sector Gram matrix.  It does not require the vectors `v_C^a`, their
coordinate labels, or a row/bridge response table.

The dimension claim can be sharpened.  `G_C` is symmetric, so it has three
independent entries.  Moreover (SQ.24) is unchanged when both `Z_C^+` and
`Z_C^-` are multiplied by a common factor.  Thus for this theorem one may
replace the two partitions by the single bias `gamma_C`.  A reduced carrier
is

```math
(d,\gamma_C,G_C),
```

with four independent real coordinates apart from the discrete order (one
bias and three Gram entries).  The larger state in (SQ.22) remains a valid
constant-dimensional upper carrier; it is not a minimality theorem.

Two scope qualifications are essential:

- a constant number of exact real coordinates is not yet a bound on the
  number of bits required at a prescribed error.  No conditioning or
  Lipschitz theorem for approximate `S_2` is proved here;
- (SQ.23) evaluates the tangent response of one pair of children.  It does
  not compute the sector--Gram state of the composed parent, and therefore
  is not yet a closed algebra under repeated composition.

Accordingly, “composes” is correct only in the restricted sense that a
fixed bilinear/rational map of the two child carriers returns the scalar
tangent coefficient.

## 3. The tangent expansion

EO.2 proves, at fixed finite `A,D,t,lambda` as `u -> 0`,

```math
\mathsf T_u
={\lambda^2u^4\over2}K_\epsilon+O(u^6),
\qquad
\mathsf M_u=O(u^8),
\qquad
\mathcal J_u
={\lambda^2u^4\over2}K_\epsilon+O(u^6).           \tag{A.2}
```

Substitution of (SQ.23) in the first and third equations proves the intended
(SQ.25).  The conversion from `rho=tanh u` creates no worse error because
`rho^4=u^4+O(u^6)`.

The source should not write `T_u=mathcal J_u` as an exact chain.  What is
true is that they have the same expansion through order four, while

```math
\mathcal J_u-\mathsf T_u=\mathsf M_u=O(u^8).
```

Also, `mathsf M_u` in EO.2 is the integrated drift of the canonical
one-row marginals relative to the canonical factors.  EO.2 does not
optimize over every row-product law.  Thus “best-product retuning term”
must be replaced by “canonical marginal-retuning term” unless a separate
local variational proof is supplied.

Finally, every remainder in (A.2) is a fixed-system remainder.  Nothing in
SQ.3 makes it uniform when the orders grow and `u=t=beta/sqrt(N)`.  The
draft's subsequent warning about this limitation is correct.

## Final judgment

After the display and terminology corrections above, SQ.3 is a rigorous
new finite response carrier: a genuinely sub-landscape child statistic
computes the first nonzero cross-row interaction response.  It is not yet a
physical-scale carrier, a finite-bit compression theorem, or a reusable
composition state.

---

# Addendum: audit of Corollary SQ.4 and Theorem SQ.5

## Addendum verdict

**SQ.4 passes with its stated tangent-only qualification.  The measure
identities and KL directions in SQ.5 also pass.**  SQ.5 needs a scope
correction: conditioning on `c` exhausts the orientation dependence of the
**canonical row product**, not the orientation dependence of the full joint
escort.  Moreover, neither `c` nor the conditional kernel is presently a
proved low-information state.  Thus “orientation-exhaustive” is valid only
with “canonical-row” inserted, and “strictly narrower” is not yet justified
in the project's information-footprint sense.

## 4. Check of Corollary SQ.4

For each sector, `C_C^a=E_(mu_(C,a,t))XX^T` is positive semidefinite.  Hence

```math
|G_C(a,b)|
\le\sqrt{G_C(a,a)G_C(b,b)}\le g_C,
```

and

```math
d+2G_C(a,b)=\operatorname {tr}(C_C^aC_C^b)
```

is nonnegative.  It is at most `d+2g_C`.  Applying these facts to the
correct formula (SQ.23), whose sector weights sum to one, gives exactly

```math
0\le K_\epsilon\le g_A(n+2g_D).
```

There is no missing symmetric factor.  The asymmetry comes from the
`i<k` versus ordered `(j,l)` convention already checked above.

Let `g=max(g_A,g_D)`.  Since `n<=N`,

```math
K_\epsilon\le g(N+2g).
```

Thus `K_epsilon>=eta N^3` forces `g>=c_eta N^(3/2)` for all sufficiently
large `N`.  In a witnessing sector,

```math
\operatorname {tr}C^2=d+2g_C,
\qquad
\lambda_{\max}(C)
\ge{\operatorname {tr}C^2\over\operatorname {tr}C},
```

and `d<=N` yields an order-`sqrt(N)` eigenvalue.  Global spin-flip symmetry
gives `E X=0`, so this eigenvalue is indeed the variance of a unit-vector
aggregate spin observable.  The constant need not depend on the aspect
ratio for this bare implication, although retaining a comparability
parameter is harmless.

At `u=beta/sqrt(N)`, `u^4K=o(N)` whenever `K=o(N^3)`.  This controls only
the formal leading tangent term.  The `O(u^6)` remainder is not uniform in
order, so “harmless” must continue to be read as “harmless at tangent
coefficient level,” exactly as the source's final paragraph warns.

## 5. Check of the fibre conditional law

On a nonempty fibre `c(b)=z`, (SQ.5) gives

```math
{dr_\epsilon\over dU_n}(b)
\propto s(b)^{-\lambda}(1+\theta_\epsilon z)^{-\lambda}.
```

The second factor is positive and constant on the fibre.  It therefore
cancels from the conditional normalization, leaving exactly `kappa_z` in
(SQ.35).  Since the canonical bridge law is the row product
`r_epsilon^(tensor m)` and the event `C_i=z_i` is rowwise, conditioning
preserves the product and proves

```math
r_\epsilon^{\otimes m}(dB\mid\mathbf C)
=\bigotimes_i\kappa_{C_i}(dB_i).
```

The kernel depends on `D,t,u,lambda` through `(s,c)`, but not on the
orientation `epsilon` or on the left child `A`.  This orientation
independence is exact.

All finite-temperature likelihoods are strictly positive.  Consequently
every contrast vector with positive canonical mass also has positive full
escort mass, so no hidden support convention is needed in the conditional
KL below.

## 6. Check of the KL directions and scale alternative

CR.6 defines the canonical error in the direction

```math
\mathcal J_\epsilon
=D(r_\epsilon^{\otimes m}\Vert q_\epsilon).
```

Apply the KL chain rule to the deterministic map
`B -> boldsymbol C`.  With `R_epsilon` the image of the first argument and
`Q_epsilon` the image of the second, the exact direction and averaging law
are

```math
D(R_\epsilon\Vert Q_\epsilon)
+E_{R_\epsilon}
 D(r_\epsilon^{\otimes m}(\cdot\mid\mathbf C)
   \Vert q_\epsilon(\cdot\mid\mathbf C)).
```

Substitution of the conditional product kernel proves (SQ.37).  The source
has not reversed either KL, and the expectation correctly uses
`R_epsilon`, not `Q_epsilon`.  Both summands are nonnegative, so the
`eta N/2` alternative follows immediately and has the correct scale.

The first term is a data-processing lower bound on the **canonical**
reverse-product error `J_epsilon`.  It should not be confused with the
globally optimized reverse product projection `I^leftarrow`.

## 7. Required information-footprint correction

The exact disintegration does not yet prove that the residual is a strict
low-information reduction.

1. SQ.2 proves that `c` is minimal for the binary row experiment only when
   `\gamma_A\ne0`.  When `gamma_A=0`, the two canonical row escorts are
   identical and the minimal orientation statistic is trivial, even though
   SQ.5 may still condition on the nonconstant function `c`.
2. A scalar-valued statistic need not have a small alphabet.  Here
   `c(b)=c(-b)`, but it may still have as many as `2^(n-1)` distinct values
   and can in principle identify a projective row word.  Across `m` rows it
   may retain order `mn` bits, the leading bridge information scale.
3. The reference kernel `kappa_z` requires the values of `s(b)` inside
   every fibre.  Thus an exact model description for the second term uses
   the full row tables `(s,c)`, not only the range of `c`.
4. The conditional full escort
   `q_epsilon(.|boldsymbol C)` can remain orientation-dependent through the
   joint likelihood ratio `p_+/p_-` in (SQ.19).  The theorem exhausts
   rowwise/canonical orientation, but it does not exhaust joint sector
   orientation.
5. Computing `Q_epsilon` or its conditional law can still require the full
   bridge escort.  No compressed evaluation theorem is supplied by the KL
   chain rule itself.

Therefore the rigorous conclusion is:

> SQ.5 is an exact canonical-row sector disintegration.  It gives a useful
> high-transport data-processing certificate or a conditional residual,
> but it does not by itself prove that either branch has a smaller
> information footprint than the original joint bridge response.

With this scope correction, the theorem is mathematically sound and the
linear-gap dichotomy is exact.

## 8. Audit of Theorem SQ.6

**The constants and conditional max-divergence claims in SQ.6 pass.**  The
nonconstant part of the row log-likelihood ratio is

```math
f(c)=\lambda\{\log(1+\theta_-c)-\log(1+\theta_+c)\}.
```

It is monotone.  Evaluating its limiting endpoint difference and using
`atanh(theta_epsilon)=gamma_D+epsilon gamma_A` gives

```math
f(1)-f(-1)
=2\lambda\{\operatorname {atanh}\theta_-
            -\operatorname {atanh}\theta_+\}
=-4\lambda\gamma_A.
```

The actual contrast range lies inside `(-1,1)`, so

```math
\operatorname {osc}\ell\le4\lambda|\gamma_A|.
```

There is no missing factor two.

For every signing, positivity of the sector sums gives
`|gamma_A|<=tQ(A)`.  If `A` is an exact pressure minimizer of order `m`, a
spin attaining `|H_A|=Q(A)` and its distinct global negative contribute

```math
F_A(t)=\log E_x\cosh(tH_A(x))\ge tQ(A)-m\log2.
```

Annealed averaging and `log cosh t<=t^2/2` give

```math
F_A(t)\le{m\choose2}\log\cosh t
        \le {m(m-1)t^2\over4}.
```

At `t=beta/sqrt(N)` and `m<=N`, these inequalities prove exactly

```math
|\gamma_A|
\le N\left(\log2+{\beta^2\over4}\right).
```

Thus (SQ.40) and the safe state-count upper bound (SQ.41) are correct.

On a quantization cell of width at most `eta`, let
`L=dr_+/dr_-`.  Then `max L/min L<=e^eta`, while

```math
{d r_+(\cdot\mid T)\over d r_-(\cdot\mid T)}
={L\over E_{r_-(\cdot\mid T)}L}.
```

The denominator lies between the cellwise minimum and maximum.  Hence this
conditional density ratio lies in `[e^(-eta),e^eta]`, proving both
directions of (SQ.42).  Conditioning a row product on its full row-label
vector again leaves a product, so tensorization gives the stated two-sided
`D_infty` bound `m eta`.  For `eta_N=o(1)` this is `o(N)` because `m<=N`.

The theorem's information conclusion needs the following precise reading.

- It constructs a polynomial-size **per-row** alphabet when `beta,lambda`
  are fixed and `1/eta` is polynomial.  The full `m`-row label vector can
  still have `(O(N/eta))^m` values and `O(m log(N/eta))` description length.
- Constructing `T_eta` requires the row likelihood-ratio/contrast table;
  SQ.6 does not compute that map from a smaller child state.
- Label distributions can retain linear response divergence.  SQ.6 makes
  the two canonical conditional reference kernels close after the labels
  are revealed; it does not make the two full conditional escorts close.
- Consequently it rules out an exponentially large **per-row canonical
  orientation alphabet** as necessary for conditional accuracy, not
  orientation-dependent joint complexity or the contrast-image dynamics.

With those qualifications, the source's explicit disclaimers are
substantially correct.  Phrases such as “orientation-blind within-cell
residual” should still be replaced by “canonical-orientation-insensitive
reference kernel,” because `q_epsilon(.|T_eta)` may retain joint sector
orientation.
