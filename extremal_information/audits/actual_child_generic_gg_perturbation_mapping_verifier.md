# Independent verifier: generic GG perturbation mapping

Status: **PASS after two local qualifications; no route-level conclusion is
strengthened**.  I checked the normalizations against Panchenko's primary
sources and reconstructed the pressure and bridge-response estimates.  The
quantitative rate and the stated bridge/escort boundary are correct.  One
coefficient-domain mismatch should be repaired in the mapping note, and the
finite-to-limit quantifiers should be stated more narrowly.

Primary sources used in this verification:

- D. Panchenko, [*Introduction to the SK
  model*](https://arxiv.org/abs/1412.0170), Section 10, especially the
  perturbation comparison and the quantitative proof of the
  Ghirlanda--Guerra identities;
- D. Panchenko, [*The Parisi ultrametricity
  conjecture*](https://arxiv.org/abs/1112.1003), for GG implying
  ultrametricity;
- D. Panchenko, [*The free energy in a multi-species
  Sherrington--Kirkpatrick model*](https://arxiv.org/abs/1310.6679),
  Theorems 2--4.

## 1. Finite normalization audit

The ordered-tuple convention in (GG.2) is exactly the source convention.
For

```math
g_{p,k}(\sigma)=k^{-p/2}\sum_{i_1,\ldots,i_p}
g_{i_1\cdots i_p}\sigma_{i_1}\cdots\sigma_{i_p},
```

independence gives

```math
E g_{p,k}(\sigma)g_{p,k}(\tau)
=k^{-p}\left(\sum_i\sigma_i\tau_i\right)^p
=R(\sigma,\tau)^p.
```

Thus (GG.2)--(GG.3) have no missing factorial or power of `k`.  For
`x_p in [1,2]`, (GG.5) is also exact:

```math
\sum_{p\ge1}4^{-p}x_p^2\le4\sum_{p\ge1}4^{-p}=4/3.
```

For a deterministic base Hamiltonian, convexity of log-sum-exp and then
`E exp(sg_x(sigma))=exp(s^2 C_x/2)` give

```math
0\le E_g\log Z_s-\log Z_0\le s^2C_x/2.
```

So (GG.6) is an **unnormalized** pressure estimate; divided by child order it
is `o(1)` whenever `s_k^2/k -> 0`.

### Required local correction to (GG.7)

Panchenko's concentration envelope `v_k(s)` is a supremum over

```math
0\le x_p\le3\qquad(p\ge1),
```

not only over the `[1,2]` cube.  The proof varies one coefficient by
`y in [0,1]`, hence needs values between `0` and `3`.  On this enlarged cube,

```math
\sup_\sigma E g_x(\sigma)^2
\le9\sum_{p\ge1}4^{-p}=3,
```

so Gaussian concentration gives, uniformly in the deterministic signing,

```math
\sup_{0\le x_p\le3}E|\log Z_s-E\log Z_s|
\le \sqrt3\,s
```

(and any larger universal constant is harmless).  Therefore the displayed
`v_k(s_k)<=Cs_k` remains correct, but its supremum should explicitly be over
`[0,3]^{\mathbb N}`.  The `4/3` constant remains correct for the actual perturbation
coefficients in `[1,2]`; it simply cannot be reused for this enlarged
concentration envelope.

With that correction, the source's exact estimate is

```math
E_x\Delta_k(f,n,p)
\le {2^p\over n}\left({2\over s_k}
+48{\sqrt{v_k(s_k)}\over s_k}\right),
```

under `v_k(s_k)/s_k^2<=4^(-p)`.  Hence `s_k=k^gamma`, with any
`0<gamma<1/2`, gives (GG.12), namely `O_{p,n}(k^(-gamma/2))`.  The lower
restriction `gamma>1/4` in the random SK example is indeed absent for a
deterministic base.

One wording qualification: the quantitative GG theorem needs
`s_k -> infinity` and `v_k(s_k)/s_k^2 -> 0`; `s_k^2/k -> 0` is additionally
needed to make the pressure perturbation negligible.  Grouping all three as
the assumptions used in the application is fine, but the third is not an
assumption of the GG identity itself.

## 2. Diagonalization and limiting GG conclusions

The limit claims are correct with the standard quantifiers.  One first
diagonalizes over a **countable convergence-determining algebra** of overlap
tests (and over integer `p,n`).  Polynomial approximation and a monotone-class
extension then give the full bounded-measurable GG identities for a
subsequential limiting weakly exchangeable Gram array.  It would be too strong
to say that one finite diagonal choice directly controls every bounded
measurable test simultaneously.

For such a limiting Gram array, the cited ultrametricity theorem applies:
positive semidefiniteness and weak exchangeability come from the finite
replica arrays, and the full GG identities also yield positivity of the
off-diagonal overlaps.  Thus (GG.13) is correctly oriented.  The statement
that the whole overlap-array law is determined by the one-overlap law is also
correct **at this limiting, annealed level**.  It is neither a finite-`k`
statement nor a quenched assertion for one frozen auxiliary perturbation.

## 3. Uniform bridge-response cost

The union bound (GG.14) is exact: every `g_x(sigma)` is centered Gaussian of
variance at most `C_x`, so a two-sided Gaussian tail followed by a union over
`2^k` configurations gives the factor `2^(k+1)`.  Consequently
`s_k||g_x||_infty=O_P(s_k sqrt(k))`.

The oscillation estimate (GG.15) is also correct **because `L(B)` is the
normalized child-Gibbs response in (GG.18)**.  Indeed, if
`d mu_h/d mu=e^h/E_mu e^h`, then

```math
e^{-osc(h)}\le {d\mu_h\over d\mu}\le e^{osc(h)}.
```

For every positive bridge weight `F_B`, this implies

```math
|\log E_{\mu_h}F_B-\log E_\mu F_B|\le osc(h).
```

Applying this to the two product factors gives (GG.15), uniformly in `B`, and
for comparable children its size is `O_P(N^(1/2+gamma))=o_P(N)`.  If `L`
instead meant an **unnormalized parent log partition function**, oscillation
alone would be false (an additive constant field is the immediate
counterexample).  The note should therefore cross-reference the normalized
definition when first stating (GG.15).

This is a scalar-response stability estimate only.  It does not say that the
replica overlap law computes the response, and it does not produce a bridge
quotient.

## 4. Multi-species hypotheses

The normalization in (GG.16) agrees with the source.  The source's doubly
indexed perturbation has conditional variance at most `4`, hence unnormalized
expected-pressure cost at most `2s_N^2`.

The synchronization conclusion (GG.17) is correct provided the following
hypotheses, implicit in the source, are retained:

1. there are finitely many species with limiting proportions
   `lambda_s>0`;
2. one has a subsequential limiting joint species-overlap array;
3. that array satisfies the **full multi-species GG identity** for every
   bounded measurable function of the species-overlap vector, obtained from
   a dense set of weight vectors and all `p`;
4. the conclusion is for the limiting annealed array.

Under these hypotheses Theorem 4 gives nondecreasing
`1/lambda_s`-Lipschitz maps `L_s:[0,1]->[0,1]`.  Positivity of the domain and
range is not an extra assumption: it follows from positive semidefiniteness,
the species GG identities, and Talagrand's positivity principle.

The mapping note correctly identifies the product-law obstruction.  A joint
weighted-species field has covariance

```math
(\lambda_1w_1R_1+\lambda_2w_2R_2)^p,
```

whose terms for `p>=2` couple the species and generally cannot be written as
an additive perturbation of the two children.  Conversely, independent
child perturbations preserve the product law but do not imply the full joint
multi-species GG identity required by synchronization.

## 5. Scope verdict

The boundary in Sections 5--6 is accurate.  GG identities compress the
annealed replica Gram law up to Hilbert-space isometry.  The actual bridge
escort queries coordinate-labelled cavity responses, is not invariant under
arbitrary Hilbert isometries, ranges over exponentially many bridge
environments, and includes a negative partition-function power.  None of the
cited theorems transfers through those changes of quantifier or measure.

Accordingly, this audit validates the mapping as a useful **strict scope
statement**, not as a RESET for the actual-child overlap SML.  The only
mathematical correction is the `[0,3]` concentration domain; the other points
are quantifier/definition qualifications.  No factor in (GG.2)--(GG.12), no
pressure exponent, and no bridge-response exponent needs changing.
