# Adversarial audit of the signed sector-mode product certificate

## Verdict

**PASS WITH SCOPE QUALIFICATIONS.**  The product-gain identity, Kronecker
carrier, stable-mode theorem, binary retuning theorem, and all displayed
constants are correct.  The stable theorem also extends to every orientation
by OU.1, with exactly the proxy and cutoff stated in the revised draft.

The signed carrier is a real improvement over the scalar sector--Gram mass:
it has only `O(m^2+n^2)` child covariance coordinates and can certify either
zero product retuning or one extensive coherent retuning direction.  It does
not, however, decide the actual-minimizer product phase uniformly.  The
remainder hypothesis is still an all-order condition, and the exact
directional observable is low-query only in a sampling/point-evaluation oracle
model.  No theorem yet forces actual optimizing children into either signed
certificate regime.

There are two wording qualifications worth retaining in the canonical
statement:

1. polynomially many **real coordinates** do not by themselves give a proved
   finite-precision or polynomial-time representation; and
2. the exact directional score avoids materializing a Gibbs table, but each
   point evaluation of `h` may still require an expensive child partition
   computation.

These qualifications do not affect the mathematical certificates.

## 1. Exact audit of the product-gain identity

From

```math
{dq\over dr}={e^{-\lambda h}\over Z},
\qquad Z=E_re^{-\lambda h},
```

one has, for every row product `P`,

```math
D(P\Vert q)=D(P\Vert r)+\lambda E_Ph+\log Z.
```

At `P=r`, this is

```math
\mathcal J=D(r\Vert q)=\lambda E_rh+\log Z.
```

Therefore

```math
\mathcal J-D(P\Vert q)
=\lambda(E_rh-E_Ph)-D(P\Vert r)=\mathscr G_h(P).
```

Taking the supremum over row products proves (SM.2), including its sign and
normalization.  Because `r` is itself a row product, the supremum is
nonnegative.  No attainment or coordinate fixed-point premise is being used.

## 2. Kronecker carrier and Frobenius normalization

Stack the bridge by left rows.  Conditional on sector `a`, the two child
spin laws are independent, so for `i\ne k`,

```math
\Gamma_{ik;j\ell}
=\sum_{a=\pm1}\pi_a^\epsilon
 (C_A^a)_{ik}(C_D^{\epsilon a})_{j\ell}.
```

Replacing `C_A^a` by `\widehat C_A^a=C_A^a-I_m` changes only the diagonal
left blocks, which are absent from `H_2`.  Hence, in row-major ordering,

```math
M=\sum_a\pi_a^\epsilon
 \widehat C_A^a\otimes C_D^{\epsilon a}
```

has diagonal blocks zero and off-diagonal block `M_(ik)=Gamma_(ik)`.  Thus

```math
{1\over2}B^{\mathsf T}MB
=\sum_{i<k}B_i^{\mathsf T}\Gamma_{ik}B_k=H_2(B).
```

Every off-diagonal block occurs twice in the symmetric matrix, so

```math
\boxed{\|M\|_F^2
=2\sum_{i<k}\|\Gamma_{ik}\|_F^2=2K_\epsilon.}
```

There is no missing factor of two in (SM.3)--(SM.5).  A matrix-vector product
does use only the two sector covariance matrices of each child and the two
sector weights.  This is `O(m^2+n^2)` stored real coordinates, although the
implicit matrix acts on an `mn`-dimensional space.

## 3. Stable-mode theorem and orientation-uniform constants

Suppose

```math
\log E_{r_{\rm row}}e^{\langle v,B\rangle}
\le {\sigma^2\over2}\|v\|_2^2.
```

Entropy duality gives, for `m_i=E_(P_i)B_i`,

```math
D(P_i\Vert r_{\rm row})
\ge\sup_v\left\{\langle v,m_i\rangle
                 -{\sigma^2\over2}\|v\|_2^2\right\}
={1\over2\sigma^2}\|m_i\|_2^2.                 \tag{A.SM.1}
```

For `P=\otimes_iP_i`, row independence and the zero diagonal blocks imply

```math
E_PH_2={1\over2}m^{\mathsf T}Mm,
\qquad E_rH_2=0.                                  \tag{A.SM.2}
```

Since `m^TMm>=-rho_-(M)||m||^2` and
`E_rR-E_PR<=osc R=Omega`, (A.SM.1)--(A.SM.2) yield

```math
\mathscr G_h(P)
\le {1\over2}\{\lambda t^2\rho_-(M)-\sigma^{-2}\}
       \|m\|_2^2+\lambda\Omega.                 \tag{A.SM.3}
```

This proves (SM.12)--(SM.14) exactly.

The orientation-uniform specialization is also correct.  CR.8 gives, in
either orientation and either row direction,

```math
D_2(r_{\epsilon,t}\Vert U_n)
\le\lambda^2t^2n\le\lambda^2\beta^2.
```

OU.1 therefore gives

```math
\sigma^2
=4\exp\{\lambda^2t^2n/2\}
\le4\exp\{\lambda^2\beta^2/2\}.                 \tag{A.SM.4}
```

Using the conservative right-hand side in (A.SM.4), the exact
orientation-uniform stability condition is

```math
\boxed{
\lambda t^2\rho_-(M)
\le {1\over4e^{\lambda^2\beta^2/2}}.}             \tag{A.SM.5}
```

Thus the theorem can be applied directly in a target-reaching orientation;
the bias-canceling orientation need not first be proved target-relevant.  In
the balanced presentation, density domination gives the alternate proxy
`sigma^2=e^(lambda(beta^2/2+log 2))`.  These two valid cutoffs do not
uniformly dominate one another.

## 4. Binary retuning theorem and constants

Central symmetry makes every odd `{-1,1}`-valued feature fair.  Hence the
tilt in (SM.20) satisfies

```math
E_(P_(i,a))B=\tanh(a)w_i,
\qquad
D(P_(i,a)\Vert r_{\rm row})
=a\tanh a-\log\cosh a=d(a).                       \tag{A.SM.6}
```

Substituting (A.SM.6) into (A.SM.2) proves (SM.22), and changing the law in
the remainder changes its expectation by at most `osc R`, proving (SM.23).
For `0<=a<=1`,

```math
d(a)\le {a^2\over2},
\qquad
\tanh^2a\ge a^2(1-2a^2/3).                       \tag{A.SM.7}
```

Under (SM.24), (A.SM.7) gives

```math
\mathscr G_(t^2H_2)(P_a)
\ge {a^2k\over2}
 \{(1+\delta)(1-2a^2/3)-1\}.
```

For `a^2=3delta/(4(1+delta))`, the braced term is `delta/2`, so the result is

```math
{3\delta^2\over16(1+\delta)}k.
```

This verifies (SM.25)--(SM.26), including the factor `3/16`.  The selected
amplitude obeys `a<1` when `0<delta<=1`.

An odd tie rule in (SM.28) always exists, for example a fixed coordinate bit
on the tie set.  What is **not** proved is a rounding theorem saying that a
negative eigenvector necessarily yields (SM.24): the row correlations
`w_i=E[B phi_i(B)]` may be small, localized, or poorly aligned after
hyperplane rounding.  The condition is correctly retained as an explicit
falsifiable hypothesis.

## 5. Remainder and sampling scope

Writing `h=t^2H_2+R+c` is exact and invariant under the arbitrary additive
constant `c` because only `osc R` is used.  Under the already stated
all-word cumulant convergence premise,

```math
\operatorname {osc}R\le2\mathfrak C_{\ge4}(t)
```

is the previously audited bound.  The signed theorem does not make that
premise easier: a linear absolute cluster tail supplies no sign and does not
force either spectral stability or rounded retuning.

The row bounded-difference calculation is also correct.  Flipping one bridge
bit changes each relevant log likelihood by at most `2t`; changing an entire
row changes `h` by at most `4tn`.  Hoeffding's lemma therefore gives the
centered log-MGF coefficient

```math
{1\over8}\sum_{i=1}^m(4tn)^2=2mt^2n^2.            \tag{A.SM.8}
```

In fact, since `m+n=N` and `t^2=beta^2/N`, (A.SM.8) is `O_beta(N^2)` for
all splits, not only comparable ones.  Independent samples consequently
estimate either expectation in (SM.29) to additive error `epsilon N` with
`O_beta(epsilon^(-2)log(1/zeta))` point evaluations.

This is a legitimate finite-query certificate **given** sampling access to
the row law and point-evaluation access to `h`.  It does not prove that those
oracles are computationally cheap, nor that an exact scalar expectation has
low finite-bit description complexity.  The safe conclusion is that no full
bridge table need be materialized, not that no hard Gibbs computation occurs.

## 6. Scalable signed counterexample to unsigned sufficiency

There is an exact physical-scale product model showing that `K` and the
absolute higher-cluster tail cannot decide the signed product branch.
This is an algebraic carrier counterexample, **not** a claim about actual
optimizing children.

Let `n` be odd, `m=n+1`, `N=m+n`, and set `c=1/n`.  With `J_s` denoting the
all-one `s by s` matrix, define

```math
M_+=c(J_m-I_m)\otimes J_n,
\qquad
M_-=-c(J_m-I_m)\otimes J_n.                       \tag{A.SM.9}
```

Both are valid sector-covariance carriers.  Indeed

```math
C_A^\pm=I_m\mathbin\pm c(J_m-I_m)
```

are sign covariance matrices: `C_A^-` is obtained from the uniform law on
balanced `m`-bit words, while `C_A^+` is obtained by mixing the independent
fair law with weight `1-c` and the all-equal fair law with weight `c`.
Also `C_D=J_n` is the covariance of an all-equal fair word.  Thus
`M_pm=(C_A^pm-I) tensor C_D` has exactly the form (SM.5).

The unsigned data agree exactly:

```math
K_+=K_-={1\over2}\|M_\pm\|_F^2
={n(n+1)\over2}.                                  \tag{A.SM.10}
```

Take the fair row product `r=U_n^(tensor m)`, let the interaction be exactly
quadratic (`R=0`), and put

```math
\lambda t^2={2\pi\over N}.                        \tag{A.SM.11}
```

The nonzero eigenvalues relevant here are

```math
\rho_-(M_+)=1,
\qquad
\rho_-(M_-)=n.                                    \tag{A.SM.12}
```

For `M_+`, the fair-row proxy is `sigma^2=1`, and
`2pi/N<1`; SM.1 with `Omega=0` gives

```math
\mathcal J_+-\mathcal I_+^\leftarrow=0.           \tag{A.SM.13}
```

For `M_-`, choose on every row the majority bit
`phi(b)=sgn(sum_j b_j)`.  There are no ties because `n` is odd.  If

```math
\alpha_n={E|\sum_(j=1)^nB_j|\over n},
```

then `w_i=alpha_n 1_n` and

```math
-{\lambda t^2\over m}w^{\mathsf T}M_-w
={2\pi\over N}\alpha_n^2n^2
\longrightarrow2,                                \tag{A.SM.14}
```

because `(E|sum_j B_j|)^2/n -> 2/pi`.  Hence (SM.24) holds, for example
with `delta=1/2`, for all sufficiently large odd `n`.  SM.2 gives

```math
\mathcal J_- -\mathcal I_-^\leftarrow
\ge {m\over32}=\Omega(N).                         \tag{A.SM.15}
```

The two models have the same `K`, the same absolute quadratic coefficients,
and zero order-at-least-four remainder, yet one has zero product retuning
and the other has extensive coherent retuning.  Thus no theorem depending
only on unsigned `K` and absolute cluster mass can decide the signed branches
in the general carrier algebra.  An optimizer-specific restriction could
still rule out one member of this pair; that is exactly the missing actual-
child input.

## 7. Information and frontier judgment

The covariance factorization (SM.5), negative edge, and rounded response
profile are demonstrably smaller than the full `2^(mn)` bridge response
landscape at the level of coordinate count.  They therefore provide
concrete low-information **sufficient certificates**:

```text
small signed negative edge + sublinear remainder
    => no extensive product retuning;

one extensive rounded negative direction + sublinear remainder
    => extensive coherent product retuning.
```

OU.1 makes the first certificate available in the target-reaching
orientation, which removes a separate orientation comparison from this
subroute.  That is a material narrowing of the certificate architecture.

It is not yet a decision theorem for `L_balanced-product-phase`.  The gap
between the two sufficient regimes is real, the rounded eigenmode need not
capture the optimal product direction, and neither `Omega=o(N)` nor a
signed-mode dichotomy has been proved for actual optimizing children.
Consequently the signed carrier supports a narrower new SML--prove that
actual children satisfy spectral stability or exhibit a macroscopic rounded
direction, with a controlled remainder--but does not by itself justify a
Level-6 recurrence or claim that the full product phase has been resolved.
