# Wave 19: Gibbs edge and row stability

This note audits whether discrete signing stability can force the lower
cavity reward required in (10.617).  The conclusion is a scoped no-go: exact
edge, row, and Hamming-`p` row variations of a **pressure minimizer** all put
an upper bound on the cavity reward.  An exact `Q`-minimizer supplies no sign
for these pressure variations.  This does not disprove (10.617); it closes
this direct variational mechanism.

## 1. Exact edge flip

Let `A` have order `r`, let `mu_A` be its oriented Gibbs law at raw inverse
temperature `beta`, and put

```math
s_e(\sigma,x)=\sigma a_{ij}x_ix_j\in\{\pm1\}
```

for `e={i,j}`.  Flipping `a_ij` changes the oriented quadratic energy by
`-4s_e`, so, exactly,

```math
P_\beta(A^e)-P_\beta(A)
=\frac1\beta\log\mathbb E_{\mu_A}e^{-4\beta s_e}.
```

If `p_e=mu_A(s_e=-1)`, this is

```math
\frac1\beta\log\left[(1-p_e)e^{-4\beta}+p_e e^{4\beta}\right].
```

Consequently, if `A` minimizes `P_beta` over all order-`r` signings,

```math
p_e\ge\frac1{1+e^{4\beta}},
\qquad
\mathbb E_{\mu_A}s_e\le\tanh(2\beta).
```

This is a one-edge marginal condition, not a cavity lower bound.

For an exact `Q`-minimizer the displayed variation still holds, but it has
no sign: a `Q`-minimizer need not minimize `P_beta`.  Writing
`L_beta=Q-P_beta`, `Q(A^e)>=Q(A)` yields only

```math
P_\beta(A^e)-P_\beta(A)
=[Q(A^e)-Q(A)]-[L_\beta(A^e)-L_\beta(A)]
\ge-\frac{r\log2}{\beta},
```

which is far below the local square-root scale.

## 2. Exact row replacement

Fix coordinate `i`, write `C=A[-i]`, and let `b` be the deleted row.  For
any replacement sign row `c in {+-1}^{r-1}`, let `A^{i to c}` retain `C`
and replace `b` by `c`.  Under the oriented child Gibbs law
`mu_{beta,C}`, define

```math
F(c)=\mathbb E_{\mu_{\beta,C}}\cosh(2\beta c^{\mathsf T}y),
\qquad
\kappa(c)=\frac1\beta\log F(c).
```

The common child partition function cancels, giving

```math
\boxed{
P_\beta(A^{i\to c})-P_\beta(A)
=\kappa(c)-\kappa(b)
=\frac1\beta\log\frac{F(c)}{F(b)}.
}
```

Thus a pressure-minimizing signing makes its actual row a **minimizer** of
the cavity reward:

```math
F(b)\le F(c)\quad\hbox{for every sign row }c.
```

Averaging uniform `c` and using independence coordinate by coordinate gives

```math
\mathbb E_cF(c)=(\cosh 2\beta)^{r-1}.
```

Therefore pressure minimality proves the wrong-direction estimate

```math
\boxed{
0\le\kappa_{\beta,i}(A)
\le\frac{r-1}{\beta}\log\cosh(2\beta).
}
```

At fixed `beta` this is linear in `r`.  At
`beta_r=r^{c-1/2}`, `0<c<1/2`, it is asymptotic to
`2r^{c+1/2}`, still larger than the required `Theta(sqrt(r))` reward.

## 3. Exact Hamming-`p` star transform

The full biased row perturbation contains no hidden reverse inequality.
Independently flip every incident edge with probability `p`, and write the
resulting row as `c=b odot t`.  In a parent state put

```math
s_j=\sigma b_jx_ix_j,
\qquad h_i=\sum_{j\ne i}s_j.
```

For a fixed `t`, the partition-function ratio is

```math
R_t=\frac{\mathcal Z_\beta(A^{i\to c})}{\mathcal Z_\beta(A)}
=\mathbb E_{\mu_A}
\exp\left(-4\beta\sum_{j:t_j=-1}s_j\right).
```

Set

```math
a_p=\sqrt{1+2p(1-p)(\cosh4\beta-1)},
```

and

```math
\lambda_p=\frac12\log
\frac{1-p+pe^{4\beta}}{1-p+pe^{-4\beta}}.
```

Factoring each coordinate gives the exact transform

```math
\boxed{
\mathbb E_tR_t
=a_p^{r-1}\mathbb E_{\mu_A}e^{-\lambda_ph_i}.
}
```

Under the child law, the last moment is

```math
\mathbb E_{\mu_A}e^{-\lambda_ph_i}
=\frac{
\mathbb E_{\mu_{\beta,C}}
\cosh((2\beta-\lambda_p)b^{\mathsf T}y)}
{
\mathbb E_{\mu_{\beta,C}}
\cosh(2\beta b^{\mathsf T}y)}.
```

If `A` minimizes pressure, every `R_t>=1`; hence

```math
\boxed{
\mathbb E_{\mu_{\beta,C}}
\cosh(2\beta b^{\mathsf T}y)
\le
a_p^{r-1}
\mathbb E_{\mu_{\beta,C}}
\cosh((2\beta-\lambda_p)b^{\mathsf T}y).
}
```

For `p=1/2`, `a_p=cosh(2 beta)`, `lambda_p=2 beta`, and the rightmost
expectation is one, recovering the uniform-row upper bound.  At `p=0` and
`p=1` the inequality is equality (no change and full-star switching,
respectively).  For every intermediate `p`, it still upper-bounds the
denominator defining `kappa`; no choice of `p` reverses the direction.

## 4. A useful near-`Q` pressure-minimizing surrogate

Although it does not solve the cavity bound, pressure minimization can be
used without losing the leading scale.  Let `A_beta` minimize `P_beta` over
order-`n` signings and let `A_*` be an exact `Q`-minimizer.  The pressure
approximation gives

```math
P_\beta(A_\beta)
\le P_\beta(A_*)\le q_n,
```

and hence

```math
\boxed{
Q(A_\beta)\le q_n+\frac{n\log2}{\beta}.
}
```

With the raw temperature

```math
\beta_n=n^{c-1/2},\qquad0<c<1/2,
```

the loss is `O(n^(3/2-c))`, a genuine power saving.  Such a
target-specific near-minimizer could therefore serve as a root in a
power-saving endpoint argument, provided a separate theorem supplies a
good principal restriction.  Its temperature changes with `n`, so it is
not the literal fixed-`beta`, every-exact-minimizer assertion (10.617).

## 5. Scope and quantifiers

1. **Exact `Q`-minimizers.** Edge and row pressure variations have no sign.
   Zero-temperature flip stability controls top deficit layers, but does not
   order finite-temperature partition functions.  The entropy-gap fallback
   is much too large.
2. **Pressure minimizers.** They are power-saving near-`Q` minimizers at the
   stated raw temperature, but discrete signing optimality makes every
   actual row minimize its cavity reward.  Uniform and Hamming-`p` averaging
   consequently yield reward upper bounds.
3. **No-go scope.** This closes the direct argument “signing replacement
   stability implies `kappa>=T-K`.”  It is not a counterexample to (10.617),
   and it does not exclude a genuinely quadratic theorem using joint
   near-ground structure, triangle/cut constraints, ancestry, or optimized
   restriction endpoints.
4. **Convergence quantifier.** A target-specific pressure-minimizing root
   with the power-saving `Q` error is admissible for a weakened endpoint
   comparison.  It does not establish the literal every-minimizer path
   criterion, and no favorable endpoint follows from the variations above.

The formulas are independently checked by
`gibbs_edge_stability_r19.py` on exhaustive small sign/state spaces.
