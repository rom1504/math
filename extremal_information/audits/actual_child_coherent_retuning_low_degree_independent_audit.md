# Independent audit: bounded-degree coherent retuning witness

**Object audited:**
[`../drafts/actual_child_coherent_retuning_low_degree.md`](../drafts/actual_child_coherent_retuning_low_degree.md)

**Verdict:** **PASS as a representation theorem, with one material scope
correction and one minor omitted estimate.**  CRW.1, CRW.12--CRW.15,
CRW.2, the quantized coefficient count, and CRW.3 are mathematically
sound.  The polynomial coefficient list is a strict finite-description
reduction for an *existential separator*, even though it depends on
`p^*`.  It is not a strict information, query, or optimization reduction
for deciding the retuning branch.  The final sentence should therefore say
`strict certificate-description reduction`, not unqualified `strict
information reduction`.

## 1. CRW.1 centering and exponential moments

Let `f=log(dq/dU_k)` and suppose each bit flip changes `f` by at most
`A/sqrt(k)`.  The bounded-differences form of Hoeffding's lemma gives

```math
\log\mathbb E_Ue^{s(f-\mathbb E_Uf)}\le{s^2A^2\over8}.
```

At `s=1`, normalization and Jensen imply, with the claimed constants,

```math
-{A^2\over8}\le\mathbb E_Uf\le0.
```

Consequently `|E_U ell|<=A^2/8` for `ell=log(dp/dr)`.  Since the one-bit
oscillation of `ell` is at most `2A/sqrt(k)`, another application gives

```math
\log\mathbb E_Ue^{s(\ell-\mathbb E_U\ell)}
\le{s^2A^2\over2}.
```

Applying this to both signs proves a dimension-free bound on
`E_U exp(s|ell|)`.  If `rho_q=dq/dU`, then the Renyi hypothesis and
Cauchy--Schwarz give

```math
\mathbb E_qe^{s|\ell|}
\le\|\rho_q\|_{L^2(U)}
   \{\mathbb E_Ue^{2s|\ell|}\}^{1/2}
\le e^{C/2}\{\mathbb E_Ue^{2s|\ell|}\}^{1/2}.
```

This verifies CRW.5a uniformly under both `r` and `p`.  There is no
centering or reference-measure error.

The proof's sentence treating the negative tail of
`Z^2=(e^ell-1)^2` is slightly abbreviated.  The fact that `Z^2<=1` there
does not by itself give uniform integrability; one also uses CRW.5a, for
example

```math
\mathbb E_r[Z^2\mathbf1_{\{\ell<-L\}}]
\le r(\ell<-L)\le K_s e^{-sL}.
```

The already proved exponential-moment estimate supplies this missing line,
so this is an exposition omission rather than a defect in CRW.1.

## 2. Chi-square-to-Jeffreys and the Walsh tail

Under `r`, put `Z=e^ell-1`.  Uniform integrability permits an `L` for
which the truncated chi-square mass is at least `eta/2`.  On
`|ell|<=L`,

```math
(e^\ell-1)\ell\ge e^{-L}(e^\ell-1)^2,
```

and hence

```math
D(p\Vert r)+D(r\Vert p)
=\mathbb E_r(e^\ell-1)\ell\ge e^{-L}\eta/2.
```

The Walsh normalization is also consistent.  With characters taking
values in `{+-1}`,

```math
\sum_{a=1}^k\mathbb E_U
 \{\ell(B)-\ell(B^{(a)})\}^2
=4\sum_S|S|\widehat\ell(S)^2,
```

so the degree tail has squared `L^2(U)` norm at most `A^2/(d+1)`.
Moreover

```math
\left\|{dp\over dU}-{dr\over dU}\right\|_2
\le2e^{C/2}.
```

The constant coefficient cancels against `p-r`, and CRW.9 follows with no
missing direction or absolute value.

## 3. CRW.12--CRW.15 for the actual optimizing product

Using `dr_i/dU_n proportional p_i^{-lambda}`, one has, up to a constant
independent of the product law,

```math
\mathbb E_Ph+{1\over\lambda}\sum_iD(P_i\Vert r_i)
=\mathbb E_P\log p_{\rm for}
 +{1\over\lambda}D(P\Vert U_B).
```

Thus a global product minimizer is a coordinate minimizer and entropy makes
that coordinate minimizer strictly positive.  Its Euler equation is
exactly CRW.12.  Expanding `h` against the definition of `r_i` also recovers
AC.17, so this is not a stationary-point assumption weaker than global
optimality.

For the canonical inverse factor, the erased-row likelihood has one-bit
log oscillation at most `2u`, hence

```math
\operatorname {osc}_{b_a}\log{dr_i\over dU_n}\le2\lambda u.
```

For `p_i^*`, AC.17 is the Gibbs law for an average of the actual bridge
landscape; flipping one bridge bit changes that averaged landscape by at
most `2u`.  Therefore

```math
\operatorname {osc}_{b_a}\log{dp_i^*\over dU_n}\le2\lambda u.
```

The conditional-cube lemma applies separately to both laws and gives

```math
D_2(r_i\Vert U_n),D_2(p_i^*\Vert U_n)
\le n\log(1+\tanh^2(\lambda u))
\le\lambda^2u^2n.
```

Hence CRW.15 is correct for every selected global product minimizer.  It
does not need the weaker relative estimate
`D_2(p_i^*||r_i)<=5lambda^2u^2n` from ES.3.

At physical scale one may take fixed
`A=4lambda beta` in CRW.1 because CRW.14 is at most
`4lambda beta/sqrt(n)`, and fixed `C=lambda^2 beta^2` in CRW.15.  Thus the
uniform application in CRW.2 has no hidden `N` dependence.

## 4. What the coefficient count does and does not reduce

Once ES.35 supplies a positive-density set of rows, CRW.1 gives fixed
degree and fixed rowwise separation.  Summing the row polynomials is valid
because both `p^*` and `r` are products.  For fixed `d`, the number of
nonconstant coefficients is

```math
m\sum_{a=1}^d{n\choose a}=O(N^{d+1}),
```

whereas explicit optimal row tables contain `m2^n` entries.  This is a
strict asymptotic reduction in the description dimension of the separating
observable.

Arbitrary real coefficients alone would not constitute a finite-bit
statement, but the quantization paragraph repairs that issue.  CRW.6 puts
each row coefficient vector in a fixed `ell_2` ball.  A mesh
`epsilon/sqrt(M)`, with `M=sum_(a<=d)binom(n,a)`, changes each row function
by at most `epsilon` in `L^2(U_n)`.  CRW.15 changes its directional gap by
at most `2e^{C/2}epsilon`; choosing fixed small `epsilon` preserves a linear
fraction of CRW.19.  Encoding the `O(N^{d+1})` rounded coordinates uses
`O(N^{d+1}log N)` bits, consistent with the displayed `N^{O(d)}log N`
bound.  Encoding the row set `I` costs only `O(N)` additional bits.

The dependence on `p^*` nevertheless matters operationally.  CRW.13
requires expectations against all the other globally optimal factors, and
the midpoint threshold in CRW.3 likewise uses optimizer-dependent
expectations.  Producing or validating those data from the children can
still require solving the full coupled product variational problem.  A
short encoding of the *answer* is not a low-query method for finding it.
Nor do the coefficients represent or reconstruct either the optimal row
tables or the full bridge landscape; they only represent one robust
separator extracted from them.

Accordingly, these two phrases need different readings:

- `strict representation reduction` in CRW.2 is correct if it means the
  finite description size of an existential separating certificate;
- `strict information reduction` in the last paragraph is too strong
  without the proposed bounded-degree child-closure lemma and should be
  replaced by `strict certificate-description reduction`.

## 5. Variance certificate and campaign consequence

CRW.6 bounds `||g_i||_2` uniformly.  Fixed-degree hypercontractivity and
CRW.15 then bound the second moment of each `g_i` under both row laws.
Independence makes variances add, so CRW.22 and the `O(1/N)` Chebyshev error
are correct.  This remains an optimizer-relative distinguishability
statement, not a child-only branch test.

Thus the draft proves a genuine polynomial-size representation theorem but
does not provide the optimizer-specific low-information observable required
to decide alternative (iii).  It does not weaken the SML and, by itself,
does not qualify as a RESET or make Level 6 credible.
