# Latent-mixture no-gain from tight Renyi-2 component complexity

**Status.** Task-local theorem note; no canonical edits. This note records the
sharp compactness consequence of the fixed-density component theorem and its
quantifier limitation. It does **not** assert a quantitative growing-density
component theorem.

## 1. Setup and the master truncation inequality

Use the notation of
`latent_mixture_row_product_no_gain.md`. Thus

```math
q_r=\int \nu_{z,r}^{\otimes r}\,\pi_r(dz),
\qquad
0\le \mathcal S_r(B)
=\left(h_\beta-{f_r(B)\over r}\right)_+
\le h_\beta,
\tag{RT.1}
```

and the audited noncentral row-product theorem supplies, for every fixed
finite `K`,

```math
\epsilon_r(K):=
\sup_{\nu:\,K_2(\nu)\le K}
E_{\nu^{\otimes r}}\mathcal S_r
\longrightarrow0.
\tag{RT.2}
```

Here

```math
K_2(\nu)=E_{U_r}\left({d\nu\over dU_r}\right)^2,
\qquad
Y_r(z):=\log K_2(\nu_{z,r})
=D_2(\nu_{z,r}\Vert U_r).
\tag{RT.3}
```

All row laws are absolutely continuous because the row cube is finite. In
particular `0\le Y_r\le r\log 2`. Splitting the exact mixture identity over
`{Y_r\le t}` and its complement gives, for every deterministic `t\ge0`,

```math
\boxed{
E_{q_r}\mathcal S_r
\le \epsilon_r(e^t)+h_\beta\,\pi_r\{Y_r>t\}.}
\tag{RT.4}
```

This is (LM.11) in logarithmic, information-theoretic coordinates. More
generally, (RT.4) holds with an order-dependent threshold `t=t_r`; the issue
is then whether `\epsilon_r(e^{t_r})` is controlled.

## 2. The compactness theorem

### Theorem RT.1 (tight component Renyi-2 complexity implies no gain)

Suppose the component complexities are asymptotically tight under their
latent laws:

```math
\lim_{T\to\infty}\limsup_{r\to\infty}
\pi_r\{Y_r>T\}=0.
\tag{RT.5}
```

Then

```math
\boxed{E_{q_r}\mathcal S_r\longrightarrow0.}
\tag{RT.6}
```

**Proof.** Fix `T<\infty`. Since `e^T` is fixed, (RT.2) and (RT.4) imply

```math
\limsup_{r\to\infty}E_{q_r}\mathcal S_r
\le h_\beta\limsup_{r\to\infty}\pi_r\{Y_r>T\}.
```

Now let `T\to\infty` and use (RT.5). `\square`

This is strictly more general than a uniform essential-supremum bound on
`K_2`. Arbitrarily large component densities are allowed, provided their
latent mass escapes to infinity only vanishingly.

### Corollary RT.2 (moment and Orlicz criteria)

Each of the following is sufficient for (RT.6):

1. `\sup_r E_{\pi_r}Y_r<\infty`;
2. `\sup_r E_{\pi_r}Y_r^p<\infty` for some `p>0`;
3. more generally, `\sup_r E_{\pi_r}\Psi(Y_r)<\infty` for some increasing
   function `\Psi:[0,\infty)\to[0,\infty)` with `\Psi(t)\to\infty`.

**Proof.** Markov's inequality gives respectively a uniform tail bound
`M/T`, `M/T^p`, or `M/\Psi(T)`, and hence (RT.5). `\square`

The first criterion is a bounded mean component `D_2` cost. It allows, for
example, components with `Y_r=\Theta(r)` on latent mass `O(1/r)`; neither
the essential supremum of `K_2` nor the largest component density need be
bounded.

This criterion is not a disguised sublinear-total-correlation assumption.
Writing `R_1,\ldots,R_r` for the rows and
`\bar\nu_r=\int\nu_{z,r}\pi_r(dz)`, the conditional-iid identities from
(LM.8)--(LM.9) give

```math
\begin{aligned}
D(q_r\Vert\bar\nu_r^{\otimes r})
&\le r I(Z_r;R_1)\\
&\le r\int D(\nu_{z,r}\Vert U_r)\,\pi_r(dz)\\
&\le r E_{\pi_r}Y_r.
\end{aligned}
\tag{RT.7}
```

Thus bounded mean `Y_r` permits `\Theta(r)` row total correlation, and the
two-half mixture in the companion note attains that order while satisfying
`Y_r=\log2` identically.

## 3. The strongest diagonal conclusion available without a rate

Pointwise convergence in (RT.2) does imply the existence of *some* slowly
growing density window, but it does not identify its growth rate.

### Lemma RT.3 (nonconstructive diagonal density window)

There is a deterministic nondecreasing sequence `b_r\to\infty` such that

```math
\epsilon_r(e^{b_r})\longrightarrow0.
\tag{RT.8}
```

Consequently, no gain follows under either of the conditions

```math
\pi_r\{Y_r>b_r\}\longrightarrow0,
\tag{RT.9}
```

or, more strongly but sometimes conveniently,

```math
E_{\pi_r}Y_r=o(b_r).
\tag{RT.10}
```

**Proof.** For each positive integer `j`, (RT.2) permits an integer `R_j`
such that

```math
r\ge R_j\quad\Longrightarrow\quad
\epsilon_r(e^j)\le {1\over j}.
```

Choose the `R_j` strictly increasing and define
`b_r=\max\{j:R_j\le r\}`, with any harmless convention before `R_1`.
Then `b_r\to\infty` and
`\epsilon_r(e^{b_r})\le1/b_r`. Apply (RT.4) with `t=b_r` to obtain (RT.9).
Condition (RT.10) implies (RT.9) by Markov's inequality. `\square`

The window `b_r` depends on the unknown convergence moduli in the component
theorem. The lemma therefore gives a qualitative compactness statement, not
an explicit admissible growth regime.

For a supplied tail envelope

```math
\pi_r\{Y_r>t\}\le T_r(t),
\tag{RT.11}
```

the exact consequence of the present proof is

```math
\boxed{
E_{q_r}\mathcal S_r
\le\inf_{0\le t\le r\log2}
\bigl\{\epsilon_r(e^t)+h_\beta T_r(t)\bigr\}.}
\tag{RT.12}
```

In particular, a quantitative no-gain conclusion follows whenever one can
choose `t_r` for which both terms in braces tend to zero.

## 4. Why no prescribed growing-density rate follows

### Proposition RT.4 (logical no-rate obstruction)

The two properties

```math
0\le\epsilon_r(K)\le h_\beta,
\qquad
K\mapsto\epsilon_r(K)\text{ nondecreasing},
\tag{RT.13}
```

together with fixed-`K` convergence (RT.2), imply no bound of the form

```math
\epsilon_r(e^{a_r})\longrightarrow0
\tag{RT.14}
```

for an arbitrarily prescribed divergent sequence `a_r`.

**Proof.** Let `a_r\to\infty` be any prescribed sequence (it may be
truncated at `r\log2`). The abstract profile

```math
\widetilde\epsilon_r(K)
=h_\beta\,\mathbf1\{\log K\ge a_r\}
\tag{RT.15}
```

satisfies (RT.13), and for every fixed `K` it is eventually zero. Yet
`\widetilde\epsilon_r(e^{a_r})=h_\beta` for every `r`. `\square`

This is a quantifier obstruction: (RT.15) is not asserted to be the actual
conference-pressure error profile. It proves that the qualitative input
(RT.2), by itself, cannot justify a named growing regime such as
`Y_r=o(r)`, `Y_r=O(\log r)`, or `K_2\le r^c`. Any such regime requires an
explicit rate in the component theorem or a new argument exploiting the
high-density components directly. For example, deterministic
`Y_r=\sqrt r` is sublinear but not tight, and (RT.2) alone says nothing
about it.

## 5. Research conclusion

The correct information-theoretic compactness variable for the current
latent-product argument is the random component complexity

```math
D_2(\nu_{Z_r,r}\Vert U_r)=\log K_2(\nu_{Z_r,r}),
```

not latent support size, row total correlation, or the essential supremum
of component density. Tightness of this variable is a rigorous sufficient
condition for no gain and bounded mean component `D_2` is an elementary
checkable criterion.

Beyond tightness, the exact frontier is (RT.12). Qualitative fixed-density
control yields only a nonconstructive diagonal window. A theorem for any
explicit divergent component-complexity scale must obtain quantitative
dependence on `K_2`, or use structural information absent from the
fixed-`K` theorem.
