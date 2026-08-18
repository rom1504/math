# Independent audit: latent-mixture row-product no-gain theorem

**Frozen source:**
`extremal_information/drafts/latent_mixture_row_product_no_gain.md`

**SHA-256:**
`b0fc78f7b4d1de219e2b813ea86fcbc0fbfcffa94beb2c24dbf983595cbc1ec3`

**Verdict:** **PASS.**  The uniform-supremum quantifier argument, measurable
mixture identity, conditional-iid total-correlation formulas, two-half
example, and density-tightness truncation are all correct.  No source repair
is required.

This verdict takes LM.3 from the cited audited noncentral row-product theorem
as the component input, exactly as the source declares.  The passage from
that input to arbitrary latent mixtures introduces no cardinality, entropy,
or measurability loss.

## 1. Sequential uniformity is equivalent to LM.3

Suppose the sequential component theorem holds for every triangular
sequence `(nu_r)` satisfying

```math
\sup_r K_2(\nu_r)\le K.
```

If LM.3 failed, there would be `eta>0`, orders `r_j -> infinity`, and row
laws `nu_(r_j)` with `K_2(nu_(r_j))<=K` such that

```math
\mathbb E_{\nu_{r_j}^{\otimes r_j}}\mathcal S_{r_j}\ge\eta.
```

Choose such a law within an arbitrarily vanishing error of the supremum if
the supremum is not attained, and fill the unselected orders with the
uniform row law.  This produces one bounded-`K_2` triangular sequence that
violates the sequential theorem.  Hence the sequential theorem implies

```math
\epsilon_r(K)=
\sup_{K_2(\nu)\le K}
\mathbb E_{\nu^{\otimes r}}\mathcal S_r\longrightarrow0.
```

Conversely, LM.3 plainly implies every such sequential statement.  Thus no
compactness or rate assumption is hidden in the word “uniform.”

## 2. Measurability and Tonelli

The row cube is finite.  Writing `p_z(a)=nu_(z,r)({a})`, measurability of the
kernel means each `z -> p_z(a)` is measurable.  For a fixed bridge
`B=(b_1,...,b_r)`,

```math
z\longmapsto \nu_{z,r}^{\otimes r}(B)
=\prod_{i=1}^r p_z(b_i)
```

is measurable.  Therefore

```math
z\longmapsto
\mathbb E_{\nu_{z,r}^{\otimes r}}\mathcal S_r
=\sum_B\mathcal S_r(B)
  \prod_{i=1}^r p_z(b_i)
```

is a finite measurable sum (indeed, a polynomial in the row probabilities).
Since `mathcal S_r>=0`, Tonelli applies on an arbitrary latent probability
space.  Expanding the definition of `q_r` gives exactly

```math
\mathbb E_{q_r}\mathcal S_r
=\int\mathbb E_{\nu_{z,r}^{\otimes r}}\mathcal S_r\,\pi_r(dz).
```

Every integrand is bounded by the same `epsilon_r(K)`, proving LM.6.  There
is no union bound and hence no dependence on latent support cardinality.

The same finite-sum argument shows that
`z -> K_2(nu_(z,r))` is measurable, as required for LM.11.

## 3. Exact information identities LM.8--LM.10

The cleanest verification uses KL chain rules rather than entropy notation.
Let `P` be the joint law

```math
P(dz,dr_1\cdots dr_r)
=\pi_r(dz)\prod_{i=1}^r\nu_{z,r}(dr_i),
```

let `q_r` be its row-vector marginal, and let `bar nu_r` be the one-row
marginal.  Conditional iid structure gives

```math
\begin{aligned}
D(P\|\pi_r\otimes\bar\nu_r^{\otimes r})
&=\int D(\nu_{z,r}^{\otimes r}
          \|\bar\nu_r^{\otimes r})\,\pi_r(dz)\\
&=r I(Z_r;R_1).
\end{aligned}
```

Decomposing the same relative entropy through the `R_1,...,R_r` marginal
gives

```math
D(P\|\pi_r\otimes\bar\nu_r^{\otimes r})
=I(Z_r;R_1,\ldots,R_r)
 +D(q_r\|\bar\nu_r^{\otimes r}).
```

Subtracting proves the equality in LM.8, including for an uncountable latent
space.  Its inequality follows from nonnegativity of mutual information.

Similarly,

```math
\begin{aligned}
\int D(\nu_{z,r}\|U_r)\,\pi_r(dz)
&=D(P_{Z_rR_1}\|\pi_r\otimes U_r)\\
&=I(Z_r;R_1)+D(\bar\nu_r\|U_r),
\end{aligned}
```

which is LM.9.  Since the row cube is finite and `U_r` has full support,
there are no absolute-continuity pathologies.

For a row density `g`, Jensen under `nu` gives

```math
D(\nu\|U)=\mathbb E_\nu\log g
\le\log\mathbb E_\nu g
=\log\mathbb E_Ug^2
=\log K_2(\nu).
```

Thus LM.5 and LM.9 imply `I(Z;R_1)<=log K`; inserting this into LM.8 proves
LM.10.  All logarithms are consistently in nats.

## 4. Two-half example

Let the row cube be partitioned into two equal disjoint halves `A_+` and
`A_-`, let `nu_+` and `nu_-` be their uniform laws, and choose one latent
sign for all `r` rows.  Each component density relative to `U_r` is

```math
g_\pm=2\mathbf1_{A_\pm},
```

so `K_2(nu_+)=K_2(nu_-)=2`, and the common row marginal is `U_r`.
On either all-plus or all-minus component support, the density of `q_r`
relative to `U_r^{otimes r}` is the constant `2^(r-1)`.  Hence

```math
D(q_r\|U_r^{\otimes r})=(r-1)\log2.
```

Equivalently, LM.8 has `I(Z;R_1)=I(Z;R_1,...,R_r)=log2`.  The example and
its claimed exact total correlation are correct.

## 5. Truncation and tightness

Because `f_r(B)>=0`,

```math
0\le\mathcal S_r(B)\le h_\beta.
```

On `{K_(z,r)<=K}`, the component expectation is at most `epsilon_r(K)`;
on its complement it is at most `h_beta`.  Splitting the exact mixture
identity therefore gives

```math
\mathbb E_{q_r}\mathcal S_r
\le(1-\alpha_r(K))\epsilon_r(K)
   +h_\beta\alpha_r(K)
\le\epsilon_r(K)+h_\beta\alpha_r(K),
```

which is LM.11.  For each fixed `K`, LM.3 makes the first term vanish, so

```math
\limsup_r\mathbb E_{q_r}\mathcal S_r
\le h_\beta\limsup_r\alpha_r(K).
```

Sending `K` to infinity under LM.12 proves the desired limit.  The order of
limits in LM.12 is therefore exactly the one needed; no rate relating `K`
to `r` is assumed.

The point-mass illustration also has the stated scale.  If
`nu_(z,r)=delta_z`, then its density is `2^r 1_{\{z\}}` and

```math
K_2(\delta_z)=2^{-r}(2^r)^2=2^r.
```

Uniformly mixing these components gives a uniform individual-row marginal
but identical rows, so it lies precisely outside the tight bounded-density
regime.

## 6. Scope

LM.1 is a theorem about mixtures of conditionally iid row laws.  It does not
cover residual row dependence after conditioning on the latent state, nor
does LM.12 control a non-tight positive mass of increasingly singular
components.  Those exclusions are explicit in the source.  Subject to them,
the latent-mixture consequence is rigorous and strictly stronger than the
sublinear-total-correlation corollary: the two-half example has linear total
correlation but remains covered.

## Final-hash confirmation

The source was subsequently updated only to replace the two provenance
hashes for its companion theorem and audit by their final frozen values:

```text
bounded_l2_noncentral_row_extension.md
e4e9a9e83e369bafabe3896e98efb8a95e8e9d49f4f70778d5ddf5b57568e282

noncentral_mean_peel_extension_adversarial_audit.md
8e0e22d0146146843de328414aae398e50cf10eaa3b5913ab56af769fc4ccbc3
```

The resulting final latent-mixture source has SHA-256

```text
720c87a8d61d8a72f317334cc9a4bd578bca19f088e8c907bd38315698ce1430
```

No theorem statement, proof step, constant, or scope condition changed.
The **PASS** verdict therefore applies unchanged to this final hash.
