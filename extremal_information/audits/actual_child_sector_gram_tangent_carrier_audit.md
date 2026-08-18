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
