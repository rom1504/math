# Independent audit of ES.3: extensive row retuning

**Object audited:** Section 5, Theorem ES.3, of
[`../drafts/actual_child_effective_row_support_dichotomy.md`](../drafts/actual_child_effective_row_support_dichotomy.md).

**Verdict:** **PASS, with two scope qualifications.**  The product-functional
identity, telescoping estimate, chi-square/Renyi constants, entropy change,
and positive-density threshold are correct.  There is no missing absolute
value in the proof.  The two qualifications are:

1. the displayed three-way asymptotic classification is exhaustive after
   passage to subsequences (or if the relevant normalized quantities have
   limits), but not literally as a trichotomy of arbitrary oscillating full
   sequences under the usual eventual meaning of `Omega(N)`; and
2. ES.3 is a nontrivial localization theorem for the gain obtained by
   retuning the canonical product, but is not yet a strict reduction of the
   optimal row-product variational problem, because the factors `p_i^*` are
   still defined by that global product oracle.

Neither qualification changes ES.33--ES.35.

## 1. Variational identity and signs

For

```math
q=q_\lambda={e^{-\lambda L}\over E_Ue^{-\lambda L}}U,
\qquad
\mathcal F(P)=E_PL+{1\over\lambda}D(P\Vert U),
```

the Gibbs identity is

```math
D(P\Vert q)=\lambda\{\mathcal F(P)-V_\lambda\}.
```

Since the canonical law `r=\bigotimes_i r_i` is an admissible row product
and `p^*` minimizes over row products,

```math
\mathcal J=D(r\Vert q),
\qquad
\mathcal I^{\leftarrow}=D(p^*\Vert q),
```

and hence exactly

```math
{\mathcal J-\mathcal I^{\leftarrow}\over\lambda}
=\mathcal F(r)-\mathcal F(p^*)\ge0.
```

Thus the equality and sign in ES.33 are correct.  Existence of `p^*` follows
from compactness of the finite product of row simplices.  Every global
minimizer is coordinatewise minimizing; the entropy term then makes its
coordinate best response the strictly positive Gibbs law AC.17, so the
Renyi estimate invoked by ES.3 applies to any chosen global minimizer.

## 2. Relative Renyi and chi-square constants

The preceding ES.3f gives

```math
D_2(p_i^*\Vert r_i)\le5C,
\qquad C=\lambda^2u^2n.
```

For `a_i^2=\chi^2(p_i^*\Vert r_i)`, the identity

```math
D_2(p_i^*\Vert r_i)=\log(1+a_i^2)
```

therefore yields

```math
a_i^2\le e^{5C}-1=X,
\qquad a_i\le\sqrt X.                       \tag{A.1}
```

The factor `5` is correct.  Both `p_i^*` and `r_i`, viewed as densities with
respect to the fair row law, have one-bit log oscillation at most
`2\lambda u`.  The logarithmic moment lemma with
`V_0=C/2` and Holder exponents `3/2,3` gives

```math
\log E_U{(p_i^*)^2\over r_i}
\le {2\over3}(9V_0)+{1\over3}(12V_0)
=10V_0=5C.
```

This also confirms that ES.3 uses chi-square in the direction
`p_i^*\Vert r_i`, which is precisely the direction needed by the subsequent
Cauchy--Schwarz estimates.

## 3. Product telescoping and energy term

Choose any row order and define intermediate products by replacing the
factors of `r` by the corresponding factors of `p^*` one at a time.  At the
step changing row `i`, average `L` over all other factors of that intermediate
product and call the resulting row function `f_i`.  The surrounding factors
do not affect the pointwise estimate

```math
|f_i(b)-f_i(b^{(j)})|\le2u.
```

For `Z=f_i-E_{U_n}f_i`, bounded differences gives, with
`v=n(2u)^2=4u^2n`,

```math
E_{U_n}Z^4\le v^2=16u^4n^2.
```

As `D_2(r_i\Vert U_n)\le C`, Cauchy--Schwarz gives

```math
\operatorname {Var}_{r_i}(f_i)
\le E_{r_i}Z^2
\le e^{C/2}(E_{U_n}Z^4)^{1/2}
=4e^{C/2}u^2n.                              \tag{A.2}
```

Applying Cauchy--Schwarz once more, now under `r_i`, proves

```math
|E_{p_i^*}f_i-E_{r_i}f_i|
\le \sqrt{\operatorname {Var}_{r_i}(f_i)}
       \sqrt{\chi^2(p_i^*\Vert r_i)}
\le2e^{C/4}u\sqrt n\,a_i.                    \tag{A.3}
```

This reconstruction is independent of which mixture of old and new factors
appears in the intermediate product.  In particular, the telescope does not
silently assume that `L` is row additive.

## 4. Entropy change identity and absolute values

Put `g_i=\log(dr_i/dU_n)`.  The exact reference-change identity is

```math
D(p_i^*\Vert U_n)-D(r_i\Vert U_n)
=D(p_i^*\Vert r_i)+(E_{p_i^*}-E_{r_i})g_i.          \tag{A.4}
```

There is no sign assumption on the second term, so taking an absolute value
is necessary.  ES.40 does so.  Since `g_i` has one-bit oscillation at most
`2\lambda u`, the same fourth-moment calculation as above gives

```math
\operatorname {Var}_{r_i}(g_i)
\le4e^{C/2}\lambda^2u^2n.                           \tag{A.5}
```

Also, monotonicity of Renyi divergence and
`\log(1+x)\le x` give

```math
D(p_i^*\Vert r_i)
\le D_2(p_i^*\Vert r_i)
=\log(1+a_i^2)\le a_i^2.                            \tag{A.6}
```

Equations A.4--A.6 yield exactly

```math
|D(p_i^*\Vert U_n)-D(r_i\Vert U_n)|
\le a_i^2+2e^{C/4}\lambda u\sqrt n\,a_i.
```

After division by `\lambda`, the two linear contributions from A.3 and the
entropy term total `4e^{C/4}u\sqrt n\,a_i`.  Using A.1 to bound
`a_i^2\le\sqrt X a_i` gives

```math
\mathcal F(r)-\mathcal F(p^*)
\le\left(4e^{C/4}u\sqrt n+{\sqrt X\over\lambda}\right)
   \sum_i a_i,
```

with exactly the constant `K` in ES.32.

## 5. Positive-density threshold

If `\mathcal J-\mathcal I^{\leftarrow}\ge\eta N`, ES.33 first gives

```math
\sum_i a_i\ge {\eta\over\lambda K}N.               \tag{A.7}
```

Let `t=\eta/(2\lambda K)` and let `k` rows have `a_i\ge t`.  From A.1,
`m\le N`, and the strict inequality below threshold,

```math
\sum_i a_i
\le k\sqrt X+(m-k)t
\le k\sqrt X+tN.
```

Comparison with A.7 yields

```math
k\ge {\eta\over2\lambda K\sqrt X}N,
```

which is ES.35.  If the displayed lower bound exceeds `m`, the antecedent is
impossible; this is consistent rather than a contradiction.  At the
physical scaling and a comparable split, fixed positive
`\beta,\lambda` make `C,X,K` order-one constants, so the conclusion is
genuinely positive density.  At the degenerate endpoint `u=0`, `X=K=0` and
ES.35 should simply be read through its vacuous antecedent (or stated for
`u>0`); the formula itself has a zero denominator there.

## 6. Counterexample search and exhaustivity

Three extremal sanity checks are consistent with the theorem:

- If `p^*=r`, then `\mathcal J-\mathcal I^{\leftarrow}=0`, so the retuning
  antecedent is false.
- If `q` is itself a row product different from `r`, then
  `\mathcal I^{\leftarrow}=0`; a linear canonical mismatch must indeed be
  carried by extensively many bounded row changes, as ES.35 says.
- If the best product remains `r` while `q` has irreducible dependence,
  then `\mathcal J=\mathcal I^{\leftarrow}` and ES.3 correctly assigns no
  retuning gain.

No counterexample to ES.33--ES.35 results.

The final asymptotic list is best interpreted subsequentially.  Under the
standard full-sequence meanings of `o(N)` and `Omega(N)`, the numerical
pattern

```math
\mathcal J_N=N,
\qquad
\mathcal I_N^{\leftarrow}=
\begin{cases}N,&N\text{ odd},\\0,&N\text{ even},\end{cases}
```

belongs to none of the three displayed cases: the retuning gap is neither
`o(N)` nor eventually `Omega(N)`.  This does not undermine the structural
dichotomy.  Because both normalized resources lie in a compact interval,
every sequence has a subsequence on which the canonical error has a limit,
and then a further subsequence on which the normalized retuning gap has a
limit.  Along that subsequence the three cases are exhaustive; oscillation
is exactly the stated possibility of mixtures of branches.

## 7. Reduction status

ES.3 is more than a restatement of the Gibbs identity: local Renyi
regularity plus the telescope proves the new fact that an extensive
improvement over the canonical row product cannot be achieved by retuning
only `o(N)` factors.  It therefore eliminates a sparse-row mechanism.

It does **not**, by itself, give a strict reduction of the optimal-product
oracle.  The quantities `p_i^*` solve the coupled mean-field equations
AC.17, and determining their common retuning direction may require the full
row-product variational problem.  A strict reduction will occur only if an
optimizer-specific statistic, strictly coarser than that oracle, either
controls `\mathcal I^{\leftarrow}` or exposes the positive-density row
shift.  The theorem source identifies this remaining obligation accurately.
