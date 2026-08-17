# Adversarial audit: bounded-cap Hadamard response packing

**Verdict: PROMOTE.**  The Walsh construction, exact cap, resolvent bound,
response packing, and information conclusion are correct.  I found no hidden
sign, indexing, or normalization error.

## 1. Algebraic construction

With `q=2^m`, `n=q^2`, and coordinates `(u,v)`, summing first over `u`
gives

```math
(Ws_g)(a,b)=q(-1)^{a\cdot b+g(a)}.
```

Thus the stated Maiorana--McFarland transform and `Wy_g=q s_g` are exact.
For `b=s_0`, the conjugate `mathcal H=D_bWD_b` is symmetric and satisfies

```math
\mathcal H^2=nI,
\qquad
\mathcal H\mathbf1=q\mathbf1.
```

Diagonal conjugation preserves trace, while the `2m`-fold Walsh tensor has
trace zero.  Hence, after removing the diagonal,

```math
H_A(x)=\frac12x^T\mathcal Hx
```

for every Boolean `x`.  The spectral bound gives absolute value at most
`qn/2`, and `x=1` attains the positive endpoint.  Therefore

```math
Q(A)=\frac12qn=\frac12n^{3/2}
```

exactly; the assertion concerns the absolute cap, not only the maximum.

## 2. Code size and switching classes

For two independent Boolean tables the bias is a sum of `q` independent
Rademacher variables, so

```math
Pr\{|S|>q/2\}\le 2e^{-q/8}.
```

Taking `N=floor(exp(q/32))`, the union bound is at most
`2 exp(-q/16)<1` for all sufficiently large `q`.  This yields the claimed
`exp(Omega(q))` family.

There is no hidden collision among the switched children.  If
`D_sAD_s=D_tAD_t`, then every nonzero off-diagonal coefficient gives
`(s_i t_i)(s_j t_j)=1`; for `n>=3`, all `s_i t_i` are equal.  Thus the only
collisions are `s=t` and `s=-t`, corresponding to equal or complementary
tables.  Both have absolute bias `q` and are excluded by the code condition.

## 3. Rayleigh and resolvent constants

For `w=s_g odot s_h`, one has `b odot w=s_{g+h}`.  Taking the inner product
of its Walsh transform gives

```math
w^T\mathcal Hw=qS(g,h)^2,
\qquad
\rho=\frac{w^T\mathcal Hw}{qn}=\frac{S(g,h)^2}{q^2}\le\frac14.
```

The proposed stronger completion-of-square estimate is correct.  On the
sphere `||u||_2^2=n`, put `K=2qI-\mathcal H`.  Then

```math
K^{-1}=\frac{2qI+\mathcal H}{3q^2},
```

and therefore

```math
\begin{aligned}
\frac12u^T\mathcal Hu+qw^Tu
&\le qn+\frac12(qw)^TK^{-1}(qw)\\
&=qn+\frac{2qn+w^T\mathcal Hw}{6}\\
&=qn\left(1+\frac{2+\rho}{6}\right)\\
&\le \frac{11}{8}qn.
\end{aligned}
```

No Boolean relaxation is being asserted as exact: enlarging the Boolean cube
to its containing sphere is used only for this valid upper bound.

## 4. Response and information separation

At its matched query, each child has response `3qn/2`.  At every unmatched
code query the preceding estimate gives at most `11qn/8`.  For children
`g,h`, their response difference is consequently at least `qn/8` at `y_g`
and at most `-qn/8` at `y_h`.  Its oscillation is at least `qn/4`, so

```math
d_proj(A_g,A_h)\ge\frac18qn=\frac18n^{3/2}.
```

The same probes give the stated sup-norm lower bound.  If one decoded state
approximated two responses, even with a separate additive calibration for
each child, the projective triangle inequality would put their distance at
most `2 epsilon n^(3/2)`.  Hence `epsilon<1/16` forces all packed children
into distinct summary states, yielding `exp(Omega(sqrt(n)))` states and
`Omega(sqrt(n))` bits.

## 5. Verification and scope

Running

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_bounded_cap_hadamard_response.py
```

completed successfully.  It checked 200 exact Walsh, Rayleigh, cap, and
response identities at the declared small orders.  The script does not by
itself certify the asymptotic random-code size or explicitly enumerate all
switching classes, but the analytic arguments above establish both.

The contextual information conclusion assumes the declared response query
family, including coordinate-pinning futures.  It should not be silently
reinterpreted as a lower bound when allowed futures are restricted, for
example, only to hollow sign-quadratic continuations.  The theorem also does
not establish the same packing for exact minimizers or for cap strictly below
`1/2`.
