# Uniform row-product no-gain passes to arbitrary latent mixtures

**Status.** Task-local theorem note; no canonical edits.  This isolates the
exact logical consequence of a uniform noncentral mean-peel row-product
theorem.  The latent space may grow with the order or be uncountable, and the
mixture may have linear row total correlation.

## 1. Uniform component input

Let `U_r` be uniform on the row cube and, for a row law `nu`, write

```math
K_2(\nu)=E_{U_r}\left({d\nu\over dU_r}\right)^2.
\tag{LM.1}
```

No centrality is assumed.  Fix the conference children, orientation, and
`0<beta<sqrt(2)/6`, and put

```math
\mathcal S_r(B)=\left(h_\beta-{f_r(B)\over r}\right)_+.
\tag{LM.2}
```

The required component input is the following uniform form of the
noncentral mean-peel theorem:

```math
\epsilon_r(K):=
\sup_{\nu:\,K_2(\nu)\le K}
E_{B\sim\nu^{\otimes r}}\mathcal S_r(B)
\longrightarrow0
\qquad(K<\infty\text{ fixed}).
\tag{LM.3}
```

It is enough that the row-product theorem be stated sequentially: if it
holds for every sequence `nu_r` with `sup_r K_2(nu_r)<=K`, then (LM.3)
follows by contradiction, choosing an almost-maximizing bad component at
each violating order.  Thus “uniform” is not an extra compactness
assumption; it is the exact quantifier content of an all-sequences theorem.

This input is supplied by the independently audited companion theorem:

```text
extremal_information/drafts/bounded_l2_noncentral_row_extension.md
sha256 e4e9a9e83e369bafabe3896e98efb8a95e8e9d49f4f70778d5ddf5b57568e282

extremal_information/audits/noncentral_mean_peel_extension_adversarial_audit.md
sha256 8e0e22d0146146843de328414aae398e50cf10eaa3b5913ab56af769fc4ccbc3
```

Its statement is uniform over all triangular row-law sequences with the
same fixed `K`, so the contradiction argument above verifies (LM.3).

The noncentrality matters only inside the component theorem.  Its mean-peel
removes the deterministic rank-one row mean before the regular covariance
argument and restores it one-sidedly by convexity.  Nothing in the mixture
step below symmetrizes the component laws.

## 2. Arbitrary latent-mixture theorem

For every `r`, let `(Z_r,pi_r)` be an arbitrary probability space and let
`z -> nu_(z,r)` be a measurable probability kernel on the finite row cube.
Define the bridge law

```math
q_r=\int \nu_{z,r}^{\otimes r}\,\pi_r(dz).
\tag{LM.4}
```

The latent space, its support cardinality, and `pi_r` may all depend on `r`.

### Theorem LM.1 (bounded-component latent mixtures cannot lower pressure)

If, for one fixed finite `K`,

```math
K_2(\nu_{z,r})\le K
\quad\text{for }\pi_r\text{-almost every }z\text{ and every }r,
\tag{LM.5}
```

then

```math
\boxed{
E_{C_r\sim q_r}\mathcal S_r(C_r)
\le\epsilon_r(K)\longrightarrow0.}
\tag{LM.6}
```

**Proof.**  The row cube is finite, so the component expectation in (LM.3)
is a measurable polynomial in the row probabilities.  Tonelli's theorem
and (LM.4) give the exact identity

```math
E_{q_r}\mathcal S_r
=\int E_{\nu_{z,r}^{\otimes r}}\mathcal S_r\,\pi_r(dz).
\tag{LM.7}
```

Every integrand is at most `epsilon_r(K)` by (LM.3), independently of `z`.
Integration proves (LM.6). `square`

There is no union bound over latent states.  Consequently a continuum of
components or a support growing faster than exponentially causes no loss.

## 3. Why this genuinely goes beyond sublinear total correlation

Let `R_1,...,R_r` have law (LM.4), and let

```math
\bar\nu_r=\int\nu_{z,r}\,\pi_r(dz)
```

be their common row marginal.  The conditional-iid information identity is

```math
D(q_r\|\bar\nu_r^{\otimes r})
=rI(Z_r;R_1)-I(Z_r;R_1,\ldots,R_r)
\le rI(Z_r;R_1).
\tag{LM.8}
```

Also,

```math
I(Z_r;R_1)+D(\bar\nu_r\|U_r)
=\int D(\nu_{z,r}\|U_r)\,\pi_r(dz).
\tag{LM.9}
```

Rényi monotonicity and (LM.5) give

```math
D(\nu_{z,r}\|U_r)\le\log K_2(\nu_{z,r})\le\log K.
```

Therefore

```math
D(q_r\|\bar\nu_r^{\otimes r})\le r\log K.
\tag{LM.10}
```

This upper bound can have linear order.  For example, split the row cube
into two uniform halves and let one latent bit choose which half supplies
all rows.  The components have `K_2=2`, while the row total correlation is
`(r-1)log 2`.  Theorem LM.1 nevertheless rules out a downward pressure
phase.  Hence linear total correlation by itself is not the obstruction;
what matters here is whether it decomposes into uniformly regular product
components.

## 4. The exact divergent-density loophole

An essential-supremum bound is stronger than necessary.  Put

```math
K_{z,r}=K_2(\nu_{z,r}),
\qquad
\alpha_r(K)=\pi_r\{z:K_{z,r}>K\}.
```

Since pressure is nonnegative,

```math
0\le\mathcal S_r\le h_\beta.
```

Splitting (LM.7) at level `K` gives the quantitative truncation bound

```math
\boxed{
E_{q_r}\mathcal S_r
\le\epsilon_r(K)+h_\beta\alpha_r(K).}
\tag{LM.11}
```

Thus LM.6 still holds whenever the component density constants are uniformly
tight under the latent laws:

```math
\lim_{K\to\infty}\limsup_{r\to\infty}\alpha_r(K)=0.
\tag{LM.12}
```

The remaining loophole is therefore precise: a nonvanishing amount of
latent mass must escape every fixed `L^2` density bound.  Mere divergence of
the essential supremum is not enough if the divergent components have
vanishing latent mass.

Without (LM.12), the fixed-`K` mean-peel theorem supplies no uniform error.
The covariance-peel rank, row transport cost, and restoration cost may all
become macroscopic as `K_{z,r}` grows.  An extreme illustration takes `z`
uniform on the row cube and `nu_(z,r)=delta_z`.  Then the mixture has a
uniform one-row marginal but all rows are identical, while

```math
K_2(\delta_z)=2^r.
```

Theorem LM.1 makes no assertion about such a law.  This example is a scope
falsifier, not evidence that its pressure is favorable.

## 5. Exact conclusion

With the audited uniform noncentral row-product theorem (LM.3), the
latent-mixture no-gain statement is **proved**, even for growing latent
support and `Theta(r)` total correlation.  The next noncovered class is not
“all latent mixtures,” but mixtures carrying non-tight component
`L^2` complexity.  Progress there requires a quantitative component theorem
with an explicit admissible growth rate for `K_2`, or a counterexample using
the escaping high-density latent mass.
