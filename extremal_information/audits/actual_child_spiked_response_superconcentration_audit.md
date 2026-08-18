# Natural-channel audit of the actual-child spiked response

Status: **rigorous recoupling theorem and exact optimizer-specific
blockade**.  The degree-two spiked response has a strict tensor contraction
when it is written as a conditional expectation in its natural probability
geometry.  Hence every *individual* high-row-order coefficient is
exponentially attenuated, despite the unit `L^2(U)` norm of the fully active
test in ST.10.  This does not control the Boolean minimum: exponentially many
attenuated coefficients may still create a linear rare excursion.

The argument below uses the exact pressure of the actual contracted-
temperature minimizing children throughout.  It identifies the missing
optimizer input as external-disorder superconcentration of that pressure
under one explicit degree-two row law.  The isolated child-minimality
identities are zero-field identities and do not imply this statement.

## 1. The carrier is an exact binary channel

Retain ST.1--ST.7.  Thus, for `y in {+-1}^n`,

```math
 z(b)={\langle y,b\rangle\over\sqrt n},
 \qquad e(b)={1+z(b)^2\over2},
 \qquad q_v(b)=e(b)+vz(b).                         \tag{SC.1}
```

Let

```math
 d\mu_y=e\,dU_n,
 \qquad a(b)={z(b)\over e(b)}={2z(b)\over1+z(b)^2}.
                                                               \tag{SC.2}
```

Then `mu_y` is a probability law, `|a|<=1`, and

```math
 {q_v\,dU_n\over d\mu_y}=1+va.                    \tag{SC.3}
```

Equivalently, if `V` is a fair sign and, conditionally on `V=v`,
`B` has law `q_vU_n`, then `B` has marginal `mu_y` and

```math
 \mathbb E[V\mid B]=a(B).                          \tag{SC.4}
```

Put

```math
 \kappa_n^2=\mathbb E_{\mu_y}a^2
 =\mathbb E_{U_n}{2z^2\over1+z^2}.                 \tag{SC.5}
```

This number does not depend on `y`.  For `n>=2`, `kappa_n<1`, because

```math
 1-\kappa_n^2
 =\mathbb E_{U_n}{(z^2-1)^2\over2(1+z^2)}>0.       \tag{SC.6}
```

Moreover, by the central limit theorem and bounded convergence,

```math
 \kappa_n^2\longrightarrow
 \kappa_\infty^2
 =\mathbb E_{G\sim N(0,1)}{2G^2\over1+G^2}<1.     \tag{SC.7}
```

Numerically `kappa_infty^2=0.6886409151...`.  In particular the row
channel has a uniform strict contraction for all sufficiently large row
widths.

## 2. Exact weighted Parseval inequality

Let `L` be the exact bridge pressure of the actual children and let

```math
 R(v)=\mathbb E_{\otimes_iq_{v_i}U_n}L(B),
 \qquad v\in\{\pm1\}^m.                            \tag{SC.8}
```

Under the joint experiment consisting of independent copies of (SC.3),

```math
 \boxed{R(V)=\mathbb E[L(B)\mid V].}               \tag{SC.9}
```

Define `phi=a/kappa_n`.  The functions

```math
 \phi_S(B)=\prod_{i\in S}\phi(B_i),
 \qquad S\subseteq[m],                             \tag{SC.10}
```

are orthonormal in `L^2(mu_y^(otimes m))`.  Consequently the Boolean Walsh
coefficients of `R` obey the exact identity

```math
 \widehat R(S)
 =\kappa_n^{|S|}
   \langle L,\phi_S\rangle_{L^2(\mu_y^{\otimes m})},             \tag{SC.11}
```

and Bessel's inequality gives

```math
 \boxed{
 \sum_{S\ne\varnothing}
   \kappa_n^{-2|S|}\widehat R(S)^2
 \le \operatorname {Var}_{\mu_y^{\otimes m}}L.}   \tag{SC.12}
```

*Proof.*  Expanding (SC.3) row by row gives

```math
 R(v)=\mathbb E_{\mu_y^{\otimes m}}
       L(B)\prod_i\{1+v_i a(B_i)\}.
```

Multiplication by `v_S` and averaging over `v` proves (SC.11).  For
nonempty `S`, `phi_S` is orthogonal to constants.  Bessel's inequality for
the orthonormal family (SC.10) proves (SC.12). `square`

This resolves an ambiguity in ST.10.  In `L^2(U_B)`, the fully active test
`prod_i z(B_i)` indeed has norm one.  In the natural Markov geometry of the
query, however,

```math
 |\widehat R([m])|
 \le\kappa_n^m
    \sqrt{\operatorname {Var}_{\mu_y^{\otimes m}}L}.             \tag{SC.13}
```

Since every augmented signing pressure satisfies
`0<=L<=beta N^(3/2)/2`, the right side is exponentially small at comparable
splits, up to a polynomial factor.  More generally,

```math
 \sum_{|S|\ge k}\widehat R(S)^2
 \le \kappa_n^{2k}
      \operatorname {Var}_{\mu_y^{\otimes m}}L.     \tag{SC.14}
```

Thus row degrees above `C log N` are negligible in `L^2(v)` for a suitable
constant `C`.  This is a genuine typical-query compression.

It is not an extremal compression.  The general conversion
`||f||_infty<=2^(m/2)||f||_2` loses the full query entropy.  At the limiting
value in (SC.7), the raw estimate obtained from (SC.14) would require more
than `1.85m` row degrees to compensate this loss, which is impossible.
Therefore weighted Parseval by itself gives no nontrivial uniform bound on
`min_v R(v)` or on the range.

## 3. The precise superconcentration implication

The conditional-expectation representation gives a clean sufficient
theorem.  Suppose that, for a number `sigma_N^2`, the actual pressure obeys

```math
 \log\mathbb E_{\mu_y^{\otimes m}}
  \exp\{s(L-\mathbb E L)\}
 \le {s^2\sigma_N^2\over2}
 \qquad(s\in\mathbb R).                            \tag{SC.15}
```

Conditional Jensen and (SC.9) give the same bound for
`R(V)-E_VR(V)`.  Since every query word has mass `2^(-m)`, optimizing the
exponential Markov bound separately at each word yields

```math
 \boxed{
 \max_v|R(v)-\mathbb E_VR(V)|
 \le\sqrt{2m\log2\,\sigma_N^2},
 \qquad
 \operatorname {range}R
 \le2\sqrt{2m\log2\,\sigma_N^2}.}                 \tag{SC.16}
```

Hence `sigma_N^2=o(N)` is sufficient for `range(R)=o(N)` at comparable
splits, and a power saving in `sigma_N^2` gives a power saving in the range.
This is stronger than (ST.16), whose universal proxy is `Theta(N)` and
therefore gives only the leading `O(N)` scale after paying for all `2^m`
queries.

Equation (SC.15) is not the full bridge table.  It is one scalar centered
MGF for each declared spike direction `y`.  It is therefore an admissible
low-information optimizer target.  But it is not currently proved.

## 4. Why isolated child minimality does not supply (SC.15)

For a left-child switching `A^v_(ik)=v_iv_kA_(ik)`, exact gauge covariance
gives

```math
 L_{A,D}(\operatorname {diag}(v)B)
 =L_{A^v,D}(B),
 \qquad
 R_{A,D,y}(v)
 =\mathbb E_{q_{+,y}^{\otimes m}}L_{A^v,D}(B).     \tag{SC.17}
```

Every `A^v` is again an exact contracted-temperature minimizing child.
Thus `R` measures how the external bridge splits an exactly degenerate
switching orbit.  Isolated minimality asserts equality of the zero-bridge
child pressures along this orbit; it gives no comparison after the common
external bridge ensemble in (SC.17) is applied.

The real extension of the pressure is convex and has gradient

```math
 \nabla_BL(B)
 =t\,\mathbb E_{\nu_B}[\tau XY^{\mathsf T}].        \tag{SC.18}
```

Therefore a standard convex-concentration proof of (SC.15) would require a
power-saving two-replica response estimate under the *spiked external
disorder* `mu_y^(otimes m)`, for example a suitable high-probability form of

```math
 \boxed{
 \sup_y
 \left\|\mathbb E_{\nu_B}[\tau XY^{\mathsf T}]\right\|_F
 =o(N).}                                           \tag{SC.19}
```

The exact flip and contraction inequalities for a minimizing child are
zero-external-field lower-MGF inequalities.  They neither upper-bound
(SC.18) nor remain available after the spiked bridge field is turned on.
The deletion--reinsertion identity controls an existing row relative to its
own deletion; it does not control all switching representatives `A^v`
against the common bridge ensemble.  Sector balancing and sector--Gram data
only control the zero-field tangent and likewise do not imply (SC.19).

This is the concrete blockade.  Fourier energy plus exact log-partition
structure proves exponential attenuation of each high-order channel and
typical-query compression, but the uniform minimum still asks for an
optimizer-to-external-field superconcentration theorem.  Without (SC.15),
(SC.19), or an equivalent power-saving tilted Dirichlet estimate, the
current identities stop at the same `Theta(N)` entropy-versus-fluctuation
balance as arbitrary bounded-difference functions.

## 5. Scope decision

The strict Markov contraction corrects the claim that the fixed-degree
carrier has no high-row-order attenuation in any natural norm.  It does not
resolve `L_balanced-product-phase`: the surviving obstruction is the rare
minimum of a strict-noise image, not one unattenuated tensor coefficient.
Accordingly this is a theorem-level sharpening and a precise missing
optimizer identity, but not a Level-6 recurrence and not by itself a reset.
