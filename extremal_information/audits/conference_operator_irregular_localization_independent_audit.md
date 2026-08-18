# Independent audit: operator-irregular conference localization

**Frozen source:**
`extremal_information/drafts/conference_operator_irregular_localization.md`

**SHA-256:**
`414f1c2c4ea9f2d4f230872ba30e4d289a233e261b7a669d9b9fa97e011747d9`

**Verdict: PASS.**  The transport inequality, effective conditional-support
conclusion, cylinder corollary, twin-row stratum, universal-double
counterexample, and stated scope are correct at the frozen hash.  I found no
required repair.  The final diffuse-product example should be read with
`u,v in {+-1}^r`; making that implicit convention explicit would improve the
exposition but does not affect any theorem.

## 1. One-coordinate pressure normalization

With

```math
H_B(x,y)={1\over2}(x,y)^TS_{\epsilon,B}(x,y)
=H_A(x)+\epsilon H_A(y)+x^TBy,
```

flipping one bridge entry changes `H_B` by exactly `+-2`.  Since
`|log cosh(a)-log cosh(b)|<=|a-b|`, the normalized log partition changes by
at most `2t`.  Thus OL.3 has the correct factor; there is no extra factor two
from the two symmetric off-diagonal blocks of `S`.

## 2. Chain-rule effective-support inequality

For any ordering `pi`, relative-entropy chain rule gives exactly

```math
D(q\|U_r)=\sum_jd_j^\pi.
```

In the hybrid from `q^(j-1)` to `q^j`, the prefix has the same `q` law in
both measures and the suffix is uniform.  Conditional on the prefix, the
suffix-averaged function of bit `pi(j)` has oscillation at most `2t`.
For a binary conditional law `p` and the uniform bit law,

```math
|E_pg-E_{U_1}g|
\le2t\,\|p-U_1\|_{TV}
\le2t\sqrt{D(p\|U_1)/2}.
```

Jensen over the prefix consequently gives

```math
|E_{q^j}f-E_{q^{j-1}}f|\le\sqrt2t\sqrt{d_j^\pi}.
```

Summing and using `sqrt(2)t=beta/sqrt(r)` proves

```math
|E_qf-E_Uf|
\le {\beta\over\sqrt r}\inf_\pi\sum_j\sqrt{d_j^\pi}
=\beta\sqrt{D(q\|U_r)s_*(q)/r}.
```

For nonzero entropy,

```math
1\le { (\sum_j\sqrt{d_j^\pi})^2\over\sum_jd_j^\pi}
\le r^2,
```

so the asserted range of `s_*` is also exact.  Combining a pressure drop
`eta r-o(r)` with `D<=Cr` yields

```math
s_*(q)\ge
\left({\eta^2\over\beta^2C}-o(1)\right)r^2,
```

with the constant in OL.8 correct.

For the negative tilt,

```math
D(q_\lambda\|U)
=-\lambda E_{q_\lambda}f-\log E_Ue^{-\lambda f}
\le\lambda E_Uf,
```

because `f>=0` and Jensen gives
`-log E exp(-lambda f)<=lambda E f`.  This verifies OL.10b.  If
`phi(lambda)=-log E exp(-lambda f)`, then `phi` is concave,
`phi'=E_(q_lambda)f`, and

```math
\phi'(\lambda)\le\phi(\lambda)/\lambda.
```

Thus a drop of the negative-moment pressure itself indeed implies the mean
drop needed for OL.10c.

## 3. Fixed-support and row-cylinder corollaries

If `q=q_C tensor U_(C^c)` and the `m` coordinates in `C` are exposed first,
all later entropy increments vanish.  Cauchy--Schwarz gives

```math
\sum_{j\in C}\sqrt{d_j}\le\sqrt{mD(q\|U)},
```

and hence the exact bound

```math
|E_qf-E_Uf|\le\beta\sqrt{mD(q\|U)/r}.
```

For conditioning on a row event `E_I`, the complementary signs remain
independent uniform and
`D(U(.|E_I)||U)=log(1/U(E_I))`.  Substituting `m=kr` and
`D<=Cr` gives `beta sqrt(Ckr)`, which is `o(r)` precisely when `k=o(r)`.
This proves OL.2, including its column analogue.

## 4. Twin-row stratum

For `k` specified rows, the first common row is free and every remaining row
matches it with probability `2^{-r}`.  Therefore

```math
U(E_I^{twin})=2^{-(k-1)r}.
```

The conditioned `k by r` submatrix is `1_k w^T`, whose only nonzero
singular value is `sqrt(kr)`.  Since `B` is an off-diagonal compression of
`S`, `||S||_op>=||B||_op`; multiplying by
`t=beta/sqrt(2r)` gives the exact lower bound
`beta sqrt(k/2)` in OL.15.

The conditional-mean estimate follows from the cylinder result with
`m=kr` and `D=(k-1)r log 2`:

```math
|E[f|E_I^{twin}]-E_Uf|
\le\beta\sqrt{k(k-1)r\log2}=O_{\beta,k}(\sqrt r).
```

For conditional bounded differences, there are `r` common-row bits with
oscillation at most `2tk` and `(r-k)r` ordinary bits with oscillation at
most `2t`.  The sum of squared oscillations is

```math
r(2tk)^2+(r-k)r(2t)^2
=2\beta^2(r+k^2-k)=O_{\beta,k}(r).
```

McDiarmid at deviation `eta r` therefore yields only an
`exp(-c r)` upper bound, exactly as OL.17 states.  The note does not claim
that this bound is sharp inside every twin stratum.

## 5. Universal-double overwrite counterexample

For `epsilon=-1` and `B^0=A+I`, the archived conference normalization is

```math
f_{-,r}(B^0)/r\longrightarrow\tau_\beta=2\psi(\beta).
```

Overwriting `k` rows with the same sign row `w` produces a `k by r`
repeated-row submatrix, so `||B^w||_op>=sqrt(kr)`.  At most `kr` bridge
bits differ from `B^0`; applying the verified one-bit oscillation gives

```math
|f(B^w)-f(B^0)|
\le2tkr
=\sqrt2\,\beta k\sqrt r.
```

For fixed `k` this is `o(r)`, proving OL.20 with the correct scaling.  Once
`k>2kappa^2/beta^2`, the parent is outside the operator-regular set while
retaining the conference pressure rate.  Different `w` give exactly `2^r`
distinct bridges, only probability `2^{r-r^2}` under the full cube.  The
source therefore correctly labels this a pointwise falsifier rather than a
speed-`r` favorable basin.

## 6. Scope and surviving stratum

The rectangular subgaussian norm tail in OL.22 has the standard
`exp(-c s^2)` form.  A fixed number of twin rows supplies a genuine
speed-`r` localized operator event, while OL.2 shows only that its
conditional **mean** remains typical.  OL.17 cannot exclude a further
speed-`r` favorable subset.  Conversely OL.3 proves that no deterministic
statement “localized spike implies typical pressure” can hold.

Thus the conclusion is exactly as narrow as advertised:

```text
entropy-O(r) plus a linear pressure gain
    forces effective conditional support Theta(r^2),
```

but this neither proves nor disproves the remaining diffuse weak-tilt phase.
It supplies no unconditioned lower-pressure LDP and makes no claim that
operator localization alone characterizes the basin.

For the final illustrative product law, taking sign vectors
`u,v in {+-1}^r` gives entry means
`alpha u_i v_j/sqrt(r)`, total entropy
`(alpha^2/2+o(1))r`, and mean singular value `alpha sqrt(r)`, as stated.
Writing the sign-vector convention explicitly would remove the only minor
notational ambiguity.

