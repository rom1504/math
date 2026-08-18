# KL information can lower conference pressure only through irregular mass

**Status.** Independent task-local verification.  This note does not modify
the frozen arbitrary-conditioning source.  It verifies the proposed
one-sided theorem for laws with `D(q||U)=O(r)` and records the centering
needed when the entropy-transport parameter is of order `sqrt(r)`.

## Proposition KI.1 (regular/irregular entropy decomposition)

Let `U_r` be uniform on bridges and let `K_r` be a regular event such that,
for fixed constants `a>0,L<infinity`,

```math
U_r(K_r^c)\le e^{-ar+o(r)}.
\tag{KI.1}
```

Suppose there is a convex extension `g_r` agreeing with the conference
pressure `f_r` on `K_r`, with exact mean `mu_r=E_U g_r`, such that

```math
\mu_r=h_\beta r+o(r),
\qquad
\log\mathbb E_U e^{-s(g_r-\mu_r)}\le Ls^2
\quad(s\ge0).
\tag{KI.2}
```

Then every sequence of bridge laws `q_r` satisfying

```math
D(q_r\|U_r)\le Cr+o(r)
\tag{KI.3}
```

obeys

```math
\boxed{
\liminf_{r\to\infty}{\mathbb E_{q_r}f_r\over r}
\ge\left(1-{C\over a}\right)_+h_\beta.}
\tag{KI.4}
```

**Proof.**  Put `p_r=U_r(K_r^c)` and
`delta_r=q_r(K_r^c)`.  Binary data processing gives

```math
D(q_r\|U_r)
\ge d_{bin}(\delta_r\|p_r)
\ge\delta_r\log(1/p_r)-\log2.
\tag{KI.5}
```

Consequently

```math
\limsup\delta_r\le\min\{1,C/a\}.
\tag{KI.6}
```

If `C<a`, this also keeps `1-delta_r` uniformly positive.  The exact KL
chain rule across the partition `{K_r,K_r^c}` gives

```math
D(q_r\|U_r)
=d_{bin}(\delta_r\|p_r)
 +(1-\delta_r)D(q_r(\cdot\mid K_r)\|U_r(\cdot\mid K_r))
 +\delta_rD(q_r(\cdot\mid K_r^c)\|U_r(\cdot\mid K_r^c)).
\tag{KI.7}
```

Thus

```math
D(q_r(\cdot\mid K_r)\|U_r(\cdot\mid K_r))=O(r).
\tag{KI.8}
```

Because `f_r=g_r` on `K_r`, (KI.2) implies

```math
\log\mathbb E_{U_r(\cdot\mid K_r)}
 e^{-s(f_r-\mu_r)}
\le Ls^2+\log{1\over U_r(K_r)}.
\tag{KI.9}
```

Apply the entropy variational inequality under `U_r(.|K_r)`.  Optimizing
at `s=Theta(sqrt(r))` using (KI.8) yields

```math
\mathbb E_{q_r(\cdot\mid K_r)}f_r
\ge\mu_r-O(\sqrt r)=h_\beta r-o(r).
\tag{KI.10}
```

This step must be centered at the exact `mu_r`: inserting an unspecified
`o(r)` error inside the MGF before taking `s=Theta(sqrt r)` would not be
justified.

Since `f_r>=0`, equations (KI.6) and (KI.10) give

```math
\mathbb E_{q_r}f_r
\ge(1-\delta_r)(h_\beta r-o(r)),
\tag{KI.11}
```

which proves (KI.4) when `C<a`.  For `C>=a`, the displayed lower bound is
zero and follows directly from `f_r>=0`. `square`

## Corollary KI.2 (entropy cost of reaching the child target)

Let `tau_beta=h_beta-gamma(beta)`.  If a law sequence satisfies

```math
\limsup {\mathbb E_{q_r}f_r\over r}\le\tau_\beta,
\tag{KI.12}
```

then every regular event satisfying (KI.1)--(KI.2) forces

```math
\boxed{
\liminf {D(q_r\|U_r)\over r}
\ge {a\gamma(\beta)\over h_\beta}.}
\tag{KI.13}
```

Indeed, apply KI.1 with any `C` strictly above the entropy liminf along a
subsequence and rearrange
`(1-C/a)h_beta<=h_beta-gamma(beta)`.  If operator thresholds provide every
`a<C_*(beta)`, taking `a` upward gives the optimized bound

```math
\liminf {D(q_r\|U_r)\over r}
\ge {C_*(\beta)\gamma(\beta)\over h_\beta}.
\tag{KI.14}
```

Strict inequalities in the definition of `C_*` cause no loss after taking
the supremum.

## Verdict and scope

The proposed one-sided theorem is **valid**.  It is compatible with the
high-pressure KL-mixture counterexample: KI.1 controls only how far the
mean can move downward, whereas a KL budget can put fixed mass on a much
rarer high-pressure event.

The theorem requires:

1. the exact conditional KL decomposition, not the unconditional transport
   bound alone;
2. `C<a` to make the regular conditional entropy `O(r)` (the other regime
   is only the trivial nonnegative bound);
3. the all-`s>=0` negative MGF for the convex extension, evaluated around
   its exact mean;
4. no upper control on pressure in the irregular sector, because
   nonnegativity alone is used there.

No hidden ESD or minimizer hypothesis enters the argument.
