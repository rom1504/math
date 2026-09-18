# Boolean BH boundedness: positive-mechanism investigation

Research start: 2026-09-18 04:29 UTC. This is the independent positive-route
track of the bounded three-hour campaign. The target is the unrestricted
complex Walsh inequality

```math
B_m=\sup_{n\ge1}\sup_{0\ne f:\deg f\le m}
 \frac{\|\widehat f\|_{2m/(m+1)}}{\|f\|_\infty},
\qquad \sup_m B_m<\infty.
```

No claim in this note proves that target or disproves it. In particular,
failure of an approximation mechanism is not a lower bound tending to
infinity for `B_m`.

## Orientation and candidate selection

The README and the September 16 independent paper audit, fresh refinement
check, and final synthesis were read. The old polynomial exponent is not a
target here. The Slote--Volberg v2 abstract/introduction was checked directly:
its polynomial weighted bootstrap does not assert uniform boundedness.

Ten candidate mechanisms considered before selecting the present route:

1. A bounded-cost Boolean dilation/decomposition, exploiting the sharp
   constant for Boolean-valued functions. Concrete advantage: an exact
   endpoint class with constant below 2 is already known.
2. A geometric-residual, approximate Boolean decomposition, allowing
   dimension-free degree inflation. Concrete advantage: exact decomposition
   has an elementary spectral-l1 obstruction, while approximation could
   bypass it.
3. A near-unimodular defect inequality with a summable error term.
4. An entropy versus `L-infinity/L2`-gap inequality at the Renyi exponent
   `m/(m+1)`.
5. A dimension reduction theorem for coefficient-norm extremizers.
6. A dual `L1` extension constructed from Riesz products of low-weight Walsh
   characters.
7. A tensor-power self-improvement of polynomial BH to bounded BH.
8. A contractive nonlinear composition that amplifies a BH violation without
   proportionately increasing degree.
9. A permutation-symmetric/Krawtchouk model theorem to locate a possible
   extremizer obstruction.
10. A completely bounded-norm interpolation or dilation with controlled
    scalar loss.

The first two were selected for immediate falsification because they supply
a literal closing argument, not only a renamed norm. Root coordination
reserved entropy, complex unimodular geometry, and counterexamples for other
tracks initially. The initial obstruction was derived independently here;
later sections explicitly label cross-track extensions and reciprocal
audits.

## 1. The tempting closing inequality

For fixed `m`, suppose every real polynomial `f` of degree at most `m`,
bounded by 1, belonged to the closed absolutely convex hull of Boolean-valued
polynomials of degree at most `m`, with coefficient mass at most an absolute
constant `C`. The sharp Boolean-valued BH inequality would immediately give
`B_m(real)<=2C`, and real/imaginary decomposition would handle complex
scalars with at most another factor 2.

This exact statement is false, even if the atom degree is an arbitrary
finite function `D(m)` and the coefficient mass is an arbitrary finite
function `C(m)`. Every Boolean polynomial of degree `D` has Fourier
coefficients in `2^(1-D) Z`; Parseval implies its Fourier l1 norm is at most
`2^(D-1)`. A bounded-total-mass combination of such atoms has bounded Fourier
l1 norm. The degree-two forms in Section 2 have Fourier l1 norm `sqrt(N)`.

One could try to repair the argument by approximating in supremum norm with
a relative residual contraction `delta<1`, then iterating. The following
theorem defeats **every degree-independent contraction**, even with arbitrary
finite atom degree and atom mass as functions of the input degree.

## 2. Verified route obstruction: no uniform Boolean-atom contraction

### Statement with quantifiers

Let `D(m)` be any positive-integer-valued function and `C(m)<infinity` any
finite nonnegative function. For each `k>=1` put

```math
m=4k,\qquad \gamma_k=1-(1-e^{-2})^k.
```

There is a sequence of real cube polynomials `f_N` of degree at most `m`,
with supremum norm at most 1, such that

```math
\liminf_{N\to\infty}\ \inf_g \|f_N-g\|_\infty\ge\gamma_k,
```

where `N` runs through powers of 2, and the infimum is over every finite sum

```math
g=\sum_\nu c_\nu b_\nu,\qquad
\sum_\nu |c_\nu|\le C(m),
```

with Boolean-valued `b_nu` of degree at most `D(m)` on the same cube. The
coefficients `c_nu` may even be complex. The conclusion also holds for
supremum-norm closures of this class.

Since `gamma_k -> 1`, there is no absolute `delta<1` for which every bounded
degree-`m` polynomial admits such an approximation with error at most
`delta`, uniformly over `m` and ambient dimension. The theorem concerns a
specific candidate proof route, **not** the BH constants themselves.

### 2.1 A bounded quadratic Hadamard form

Let `U=H_N/sqrt(N)` be a normalized real Sylvester Hadamard matrix, so
`U^T U=I` and every entry is `+/-N^(-1/2)`. On the `2N`-cube define

```math
F_N(x,y)=\frac1N x^T U y.
```

Cauchy--Schwarz gives `|F_N|<=1`. The polynomial has degree 2. Its `N^2`
Fourier coefficients all have absolute value `N^(-3/2)`, so

```math
\|\widehat F_N\|_1=\sqrt N,\qquad
\|\widehat F_N\|_{4/3}=1,\qquad
\|F_N\|_2=N^{-1/2}.
```

Under independent uniform signs, `E F_N=0` and `E F_N^2=1/N`.

### 2.2 A locally uniform but globally correlated distribution

Let `G` be an `N`-dimensional standard Gaussian vector and set `H=U^T G`.
Conditionally on `G`, draw all `2N` signs independently with means

```math
\mathbb E[X_i\mid G]=\sin G_i,\qquad
\mathbb E[Y_j\mid G]=\sin H_j.
```

Write `mu_N` for the resulting law on the cube. The elementary Gaussian
characteristic function gives, for standard Gaussians of correlation `r`,

```math
\mathbb E\sin A\sin B=e^{-1}\sinh r.
```

Here `r=U_ij`. Consequently

```math
\mathbb E_{\mu_N} F_N
=\frac{e^{-1}}N\sum_{i,j}U_{ij}\sinh(U_{ij})
=e^{-1}\sqrt N\sinh(N^{-1/2})\longrightarrow e^{-1}.
```

The stochastic sine rounding is useful because its concentration is
elementary; no four-sign Gaussian identity is being assumed.

### 2.3 A quantitative variance estimate

Condition on `G`, and abbreviate `a=sin G`, `b=sin(U^T G)`, `X=a+xi`,
`Y=b+eta`. The noises `xi,eta` are independent and centered, and both
coordinate covariance matrices are diagonal and at most `I`. The three
terms

```math
N^{-1}\xi^T U b,\qquad N^{-1}a^T U\eta,
\qquad N^{-1}\xi^T U\eta
```

are pairwise uncorrelated conditionally on `G`. Each conditional variance
is at most `1/N`: for the first two use `||a||,||b||<=sqrt(N)` and
orthogonality of `U`; for the third use `||U||_F^2=N`. Thus

```math
\mathbb E\operatorname{Var}(F_N\mid G)\le3/N.
```

The conditional mean is

```math
M(G)=N^{-1}\sin(G)^T U\sin(U^T G).
```

Its gradient is the sum of
`N^(-1) diag(cos G) U sin(U^T G)` and
`N^(-1) U diag(cos(U^T G)) U^T sin G`. Each has Euclidean norm at most
`N^(-1/2)`. Gaussian Poincare therefore gives `Var(M)<=4/N`, and the
conditional-variance identity proves

```math
\operatorname{Var}_{\mu_N}(F_N)\le7/N.
```

In particular `F_N -> e^(-1)` in `L2(mu_N)`.

### 2.4 Uniform local indistinguishability

Select any `r` of the `X` coordinates and any `s` of the `Y` coordinates,
where `r+s<=K` and `K` is fixed. Before rounding, the corresponding Gaussian
covariance matrix is

```math
\Sigma=\begin{pmatrix}I_r&A\\ A^T&I_s\end{pmatrix},
\qquad |A_{ij}|=N^{-1/2}.
```

Its off-diagonal block has squared Frobenius norm `rs/N<=K^2/(4N)`, and
operator norm at most `K/(2sqrt(N))`. For sufficiently large `N` it is
nonsingular. The Gaussian relative entropy from the independent standard
law is

```math
\tfrac12(\operatorname{tr}\Sigma-(r+s)-\log\det\Sigma)
=-\tfrac12\log\det(I-A^T A)
\le \|A\|_F^2
\le K^2/(4N),
```

where the inequality uses all eigenvalues of `A^T A` at most `1/2`.
Pinsker gives total variation at most `K/sqrt(8N)`. Applying the same
independent sine-rounding kernels cannot increase total variation.
Under independent Gaussian inputs those rounded signs are independent
uniform signs because sine is odd. Therefore every set of at most `K`
coordinates has marginal distance at most `K/sqrt(8N)` from uniform,
uniformly over which coordinates are selected.

Every Boolean-valued function of degree at most `D` depends on at most
`D 2^(D-1)` coordinates. An elementary proof suffices here: if variable `i`
is relevant, its derivative `(b(x_i=1)-b(x_i=-1))/2` is a nonzero polynomial
of degree at most `D-1`. A nonzero degree-`d` multilinear polynomial is
nonzero on at least a fraction `2^(-d)` of the cube, by induction on a
variable appearing in a top-degree monomial. Since the derivative is in
`{0,+/-1}`, its influence is at least `2^(1-D)`. Total influence is at most
`D` by Parseval. This gives the asserted number of relevant variables.

Thus, for fixed `D`, every such Boolean function satisfies

```math
|\mathbb E_{\mu_N}b-\mathbb E_{\mathrm{unif}}b|
\le \frac{D2^{D-1}}{\sqrt{2N}}.
```

A sum with total coefficient mass at most `C` has expectation gap at most
`C D2^(D-1)/sqrt(2N)`, which tends to zero. This is uniform in the number
of atoms and their identities.

### 2.5 Polynomial amplification and closing the obstruction

For fixed `k`, use the real polynomial

```math
P_k(t)=1-2(1-t^2)^k,
\qquad f_N=P_k(F_N).
```

It is bounded by 1 on `[-1,1]`, so `||f_N||_infinity<=1`; its reduced
Walsh degree is at most `4k`. Continuity and boundedness of `P_k`, together
with the two `L2` limits above, imply

```math
\mathbb E_{\mathrm{unif}} f_N\longrightarrow-1,
\qquad
\mathbb E_{\mu_N} f_N\longrightarrow
1-2(1-e^{-2})^k.
```

The expectation gap of `f_N` tends to `2 gamma_k`. If `g` is any admissible
Boolean-atom sum, its gap tends uniformly to zero. Since the difference
between expectations under two probability measures is bounded in
absolute value by twice the supremum norm, this proves

```math
2\|f_N-g\|_\infty\ge2\gamma_k-o(1),
```

uniformly over `g`, as claimed.

An entirely finite version follows from `|P_k'|<=4k`, the variance bound,
`|E_mu F_N-e^(-1)|<=1/N`, and
`0<=E_unif P_k(F_N)+1<=2k/N`. If `K=D2^(D-1)`, then every atom sum of mass
at most `C` obeys

```math
\|f_N-g\|_\infty
\ge\gamma_k-2k\sqrt{7/N}-3k/N-CK/(2\sqrt{2N}),
```

once `K^2/(4N)<=1/2`. The mean-bias estimate uses
`sinh(z)/z-1<=z^2` for `0<=z<=1`, directly from its positive power series.

### 2.6 Audit notes and scope

- All examples are real, hence also belong to the requested complex class.
- The input degree `m=4k` is fixed before `N->infinity`; only afterward is
  `k` chosen large to defeat a proposed absolute contraction.
- Atom degree and total coefficient mass may depend arbitrarily on `m`, but
  not on ambient dimension. No restriction on the number of atoms is used.
- The amplified examples do not contradict bounded BH. In fact the base
  form has BH coefficient norm exactly 1. Its purpose is to separate a
  boundedness proof route from its proposed Boolean endpoint class.
- Gaussian Poincare, the Gaussian KL formula, Pinsker, and data processing
  for total variation are the standard analytic inputs. All task-specific
  calculations are supplied above.

Monte Carlo sanity replay:
`computations/bh_boundedness_positive_2026_09_18_rounding.py`, default seed
`20260918`, `k=8`, and 8192 samples at each `N=64,256,1024`. The observed
correlated means were respectively `0.3681622,0.3684225,0.3681625`; exact
means are `0.3688382,0.3681190,0.3679393`. The observed `N Var(F_N)` values
were `0.72239,0.72585,0.74724`, consistent with (but not proving) the upper
bound 7. The amplified gaps were `1.10061,1.29916,1.35579`, tending toward
the proved limit `1.3750980` for `k=8`. The script writes no files.

### 2.7 Extension to every finite-order phase atom

The same obstruction holds with Boolean atoms replaced by functions of
the form `eta b`, where `|eta|=1` and `b` takes values in roots of unity.
The order of those roots may vary between atoms and grow with ambient
dimension. The only restrictions remain degree at most a finite `D(m)`
and total coefficient mass at most a finite `C(m)`, both independent of
ambient dimension. In particular, arbitrary finite degree inflation does
not rescue a uniform residual contraction `delta<1` through this endpoint
class.

Here is the additional ingredient, independently checked from the Galois
argument in `bh_boundedness_barrier_2026_09_18.md`, Section 7.2. If `b` has
degree at most `D` and values in `mu_L`, every `2^D bhat(S)` is an algebraic
integer in `Q(zeta_L)`. Indeed, multilinear interpolation on `{0,1}^n`
expresses its coefficients as integral combinations of vertex values;
the change to Walsh monomials introduces denominators dividing `2^D`.
Every field embedding `sigma` preserves degree and sends all vertex
values to roots of unity, so Parseval gives

```math
\sum_S |\sigma(\widehat b(S))|^2=1.
```

All embeddings have the same nonzero Fourier support. Write `r` for the
number of embeddings. For a nonzero coefficient, the nonzero integral
field norm gives

```math
\prod_\sigma |\sigma(\widehat b(S))|^2\ge 2^{-2Dr}.
```

Generalized Holder therefore implies

```math
|\operatorname{supp}\widehat b|4^{-D}
\le\sum_{S\in\operatorname{supp}\widehat b}
       \prod_\sigma |\sigma(\widehat b(S))|^{2/r}
\le\prod_\sigma\left(\sum_S
       |\sigma(\widehat b(S))|^2\right)^{1/r}=1.
```

Thus `b` depends on at most `K=D 4^D` coordinates, uniformly in `L`.
Multiplication by a global phase changes no dependence. The local total
variation estimate in Section 2.4 already applies to complex functions
of modulus at most one. Replacing `K=D 2^(D-1)` by `K=D 4^D` in the proof
and its finite-error estimate gives exactly the same lower limit
`gamma_k=1-(1-exp(-2))^k`. Degree-zero atoms are constant and have zero
expectation discrepancy, so that case causes no exception.

This extension was suggested by the barrier track after its independent
audit of Section 2; the algebraic-integer, Holder, junta, and quantifier
steps above were rechecked here. It strengthens the approximation no-go,
not the unrestricted BH lower bound.

### 2.8 Final extension: all unimodular atoms, with arbitrary phases

The qualitative junta theorem proved in Section 9.4 removes the torsion
assumption altogether. For each fixed `D`, every degree-at-most-`D`
complex polynomial of modulus identically one depends on at most the
finite number `J_D` of variables, independently of ambient dimension and
of its phase alphabet.

Replace the atoms in Section 2's main statement by **arbitrary unimodular
polynomials** of degree at most `D(m)`, and set `K=J_(D(m))` in the local
indistinguishability estimate. The proof, finite-error estimate, and all
quantifiers are unchanged. In particular, for every finite
dimension-independent `D(m)` and `C(m)`, the same `f_N` satisfies

```math
\liminf_{N\to\infty}\inf_{
 \substack{g=\sum_\nu c_\nu b_\nu,
             \ \sum_\nu|c_\nu|\le C(m)\\
             |b_\nu|\equiv1,\ \deg b_\nu\le D(m)}}
\|f_N-g\|_\infty\ge\gamma_k,
\quad m=4k,\quad\gamma_k\longrightarrow1.
```

Supremum-norm closures are included, and arbitrary global phases require
no separate treatment. This rules out every absolute residual contraction
`delta<1` by bounded-total-mass, dimension-independent-degree unitary
atoms. It does not rule out reductions whose atom degree grows with
ambient dimension, or approximate rather than exact unit modulus. It is
a barrier to a uniform boundedness mechanism, not a divergent BH example.

## 3. Verified positive class theorem: permutation-symmetric polynomials

For every `n,m>=1` and every complex-valued permutation-symmetric function
`f` on the `n`-cube of Walsh degree at most `m`,

```math
\boxed{\|\widehat f\|_{2m/(m+1)}
 \le\sqrt{6e}\,\|f\|_\infty.}
```

The constant is not optimized. This eliminates this symmetric class as a
source of diverging unrestricted BH constants. It is not a solution for
general coefficient phases or magnitudes.

### 3.1 A spin-rotation influence bound

Assume `1<=m<=n`; the case `m>n` follows afterward by coefficient-norm
monotonicity. Write `p(k)` for the value of `f` on inputs with exactly `k`
negative coordinates. Symmetry and the elementary-symmetric Walsh expansion
show that `p` is an ordinary polynomial in `k` of degree at most `m`.

On `(C^2)^(tensor n)` let

```math
J_x=\frac12\sum_{i=1}^n X_i,\qquad
J_z=\frac12\sum_{i=1}^n Z_i,\qquad
F=p(n/2-J_z).
```

Here `X,Z` are the usual Pauli matrices. The diagonal operator `F` has norm
`||f||_infinity`. Under rotation,

```math
A(t)=e^{itJ_x}F e^{-itJ_x}
```

is an operator-valued trigonometric polynomial of degree at most `m`:
conjugation replaces `J_z` by a linear combination of `J_z cos t` and
`J_y sin t`, and expansion of each power introduces frequencies of
absolute value at most its degree. Noncommutativity changes the matrix
coefficients, not this frequency bound. Also `||A(t)||=||F||` for real `t`.

Apply the scalar trigonometric Bernstein inequality to
`<u,A(t)v>` for arbitrary unit vectors `u,v`, then take the supremum.
It follows that

```math
\|[J_x,F]\|\le m\|f\|_\infty.
```

Restrict to the symmetric tensor subspace, with normalized Dicke vectors
`|k>` (`k` negative coordinates). The only nonzero off-diagonal entries of
`J_x` connect successive `k`, and

```math
\langle k+1|J_x|k\rangle
=\tfrac12\sqrt{(k+1)(n-k)}.
```

Thus, for every `k`, the squared norm of the `k`th commutator column is

```math
\frac14\bigl[(k+1)(n-k)|p(k+1)-p(k)|^2
+k(n-k+1)|p(k)-p(k-1)|^2\bigr]
\le m^2\|f\|_\infty^2,
```

with nonexistent endpoint terms omitted. Average this with
`w_n(k)=2^(-n) binom(n,k)`. The coefficient of
`|p(k+1)-p(k)|^2` becomes

```math
\frac14(k+1)(n-k)[w_n(k)+w_n(k+1)]
=\frac{n(n+1)}8 w_{n-1}(k).
```

The total influence is exactly

```math
I(f)=\sum_S |S||\widehat f(S)|^2
=\frac n4\sum_{k=0}^{n-1}w_{n-1}(k)|p(k+1)-p(k)|^2.
```

Therefore

```math
\boxed{I(f)\le\frac{2m^2}{n+1}\|f\|_\infty^2.}
```

This proof is valid directly for complex `p`; no real/imaginary factor is
needed.

### 3.2 Converting influence to the critical coefficient norm

Let `a_r` be the common coefficient on level `r`, and put
`v_r=binom(n,r)|a_r|^2`. Then `sum v_r=||f||_2^2` and
`sum r v_r=I(f)`. At `q=2m/(m+1)`, the norm on each coefficient level
satisfies

```math
\|\widehat f^{=r}\|_q^2=\binom nr^{1/m}v_r.
```

Counting-measure comparison across the at most `m+1` levels gives

```math
\|\widehat f\|_q^2
\le(m+1)^{1/m}\sum_{r=0}^m\binom nr^{1/m}v_r.
```

For `1<=r<=m`, write `s=r/m` and `t=n/m>=1`. The elementary binomial bound,
the inequality `s(1-log s)<=1`, and convexity of `t^s` on `[0,1]` yield

```math
\binom nr^{1/m}
\le(et/s)^s\le e\,t^s
\le e(1+s(t-1))\le e(1+nr/m^2).
```

The last upper bound also holds for `r=0`. Since
`(m+1)^(1/m)<=2`, the spin-rotation estimate now gives

```math
\|\widehat f\|_q^2
\le 2e\left(\|f\|_2^2+\frac n{m^2}I(f)\right)
\le2e\left(1+\frac{2n}{n+1}\right)\|f\|_\infty^2
\le6e\|f\|_\infty^2.
```

If `m>n`, apply the proved statement at degree bound `n` and use
`q_m>=q_n` and monotonicity of counting-measure coefficient norms. This
establishes all advertised quantifiers.

### 3.3 Why the spin argument does not prove the unrestricted theorem

For arbitrary `f`, simultaneous spin rotation still gives
`||[J_x,f(Z)]||<=m||f||_infinity`. But the Dicke subspace need not be
invariant under `f(Z)`. Without that invariant symmetric subspace the
hypercube edges do not combine into the amplified coefficients
`sqrt((k+1)(n-k)) Delta p(k)`. The normalized Hilbert--Schmidt estimate only
recovers `I(f)<=m^2||f||_infinity^2`, which is insufficient. This is the
exact missing structural input, not a suppressed symmetry assumption.
Ordinary Parseval separately gives the stronger general estimate
`I(f)<=m||f||_infinity^2`; neither supplies the needed factor of order
`m^2/n` when the ambient dimension is large.

### 3.4 Numerical falsification attempts before the proof

`computations/bh_boundedness_positive_2026_09_18_symmetry.py` represents
degree-`m` symmetric polynomials in the Chebyshev basis on the **discrete**
Hamming-weight grid. It uses linear programs for central differences and
iterated supporting linear programs to seek large total influence. These
are floating-point local searches, not global certificates.

The first replay used seed `20260918`, 8 starts, and 6 supporting-plane
rounds, for `(n,m)` equal to `(4,2),(8,2),(16,2),(16,4),(32,4),(64,4),
(64,8),(128,8),(256,8),(256,16),(512,16),(1024,16)`. The largest found
`n I(f)/m^2` in this list was `0.75`; the tested values were all below
`0.751`. The proof above gives the rigorous upper bound `2n/(n+1)` for
all dimensions/degrees and does not rely on that computation. The script
also prints the candidate Chebyshev coefficients so each numerical witness
can be replayed without relying on an unstored optimizer state.

## 4. Verified closure bounds and a degree-compression warning

### 4.1 Tensor products cannot amplify a fixed critical ratio

Let `f=product_i f_i` on disjoint coordinate blocks, with nonzero factors of
degrees bounded by positive `d_i`, and put `M=sum_i d_i`. Interpolation of
counting-measure sequence norms uses the exact identity

```math
\frac1{q_M}=\frac{d_i}{M}\frac1{q_{d_i}}
+\left(1-\frac{d_i}{M}\right)\frac12.
```

Writing `R_i=||fhat_i||_(q_di)/||f_i||_infinity`, Parseval gives

```math
\frac{\|\widehat f\|_{q_M}}{\|f\|_\infty}
\le\prod_i R_i^{d_i/M}
\le\max_i R_i.
```

Thus no disjoint tensor-product family built from bounded-degree seeds can
produce a diverging critical BH ratio. This assertion does not assume flat
coefficients or unimodularity. Constant factors can simply be removed.

### 4.2 Addressing/selector gates contract constants toward 2

Let `f,g` have caps at most 1 and degree at most `m-1`; they may share
coordinates. Add a fresh selector variable `s` and define

```math
F(s,x)=\frac{1+s}{2}f(x)+\frac{1-s}{2}g(x).
```

Suppose each nonconstant branch, at its own positive degree bound, has
critical coefficient norm at most `C`, where `C>=1`. Constants have norm
at most 1. Interpolation with Parseval then gives
`||fhat||_(q_m),||ghat||_(q_m)<=C^((m-1)/m)`.

For complex scalars `a,b` and `1<=q<=2`, the elementary inequality

```math
|a+b|^q+|a-b|^q\le2(|a|^q+|b|^q)
```

follows by concavity of the `q/2` power and the parallelogram identity.
Applying this to each common Fourier index of the two branches yields

```math
\|\widehat F\|_{q_m}
\le 2^{1/m} C^{(m-1)/m}.
```

In particular a constant `C>=2` is preserved by such gates. Together with
the tensor bound, this excludes amplification by the usual product and
addressing operations at their natural degree bounds.

**Essential warning:** the selector estimate uses the degree bound
`1+max(deg f,deg g)`. If the highest-degree parts of `f` and `g` agree,
the actual degree of `F` can be smaller. Replacing the natural degree bound
by that smaller actual degree is not licensed by the calculation. An
operation exploiting systematic degree cancellation/compression remains
outside this closure theorem.

## 5. Tested stronger shortcuts that do not close boundedness

- A global `l1`/`L2` interpolation estimate of the form
  `||fhat||_1 ||f||_2^(m-1)<=C^m||f||_infinity^m` would imply the desired
  BH inequality but is false. Add a fixed nonzero constant (or a dominant
  monomial on disjoint variables) to a small multiple of the Hadamard
  quadratic form: its cap and L2 norm stay bounded while its Fourier l1
  norm grows as `sqrt(N)`.
- A uniform finite-degree Boolean approximation with any fixed error below
  1 is defeated by Section 2. Allowing enormous but dimension-independent
  atom degree or coefficient mass does not repair it.
- The symmetric spin argument is genuinely structural, not an unrestricted
  degree argument; Section 3.3 identifies exactly where the symmetry is
  used. Symmetrizing arbitrary coefficient phases cannot simply be assumed
  to preserve the target coefficient norm.

## 6. Independent audits and campaign coordination

At 2026-09-18 05:01--05:04 UTC, the independent barrier-track agent
reconstructed Sections 2 and 3 in full and reported both proofs correct.
Its audit checked the three conditional noise variances, Gaussian gradient
bound, determinant/Pinsker constant, elementary junta bound, all quantifier
orders, the finite approximation-error bound, the spin-rotation frequency
bound for complex matrices, the exact Dicke-column average, the coefficient
norm conversion, and the `m>n` case. This is an independent mathematical
reconstruction, not just a successful numerical replay.

The same agent correctly noted that for unrestricted functions Parseval
already gives `I(f)<=m||f||_infinity^2`, stronger than the generic
commutator argument's `m^2` bound. Neither is the needed `m^2/n` estimate;
the full permutation symmetry remains essential to Section 3.

### 6.1 Reciprocal audit of the torsion non-density theorem

I independently reconstructed Section 8 of
`artifacts/bh_boundedness_barrier_2026_09_18.md` and its exact replay at
2026-09-18 05:04 UTC. Verdict: **passes**.

The explicit unimodular eight-tuple has only the full vanishing subsum:
the imaginary radical coefficients force complete conjugate pairs, and
the remaining real equation is `15h-6b-10c-14d=0`. Multiplication by the
top parity removes the cubic Fourier coefficient and gives a degree-two
unimodular function on three variables.

The supplied special-case proof of Mann's theorem is valid. For a prime
power `p^a` in the minimal common order, coprime cyclotomic degrees give
the asserted minimal polynomial over the remaining cyclotomic field.
Dividing the grouped polynomial by that cyclotomic polynomial forces
equal coefficients within each residue class modulo `p^(a-1)`.
Irreducibility and the normalized term 1 force one such residue, so
`a>=2` is impossible. When `p>k`, an empty exponent class forces every
grouped coefficient to vanish; irreducibility then removes `p`. Thus the
term ratios lie in `mu_210` for an irreducible eight-term relation.

The target ratio has trace `-2/5`, and therefore is not a root of unity:
that rational trace would have to be an integer if it were a sum of
algebraic integers. The code's quantitative checks are stronger. I
reviewed its rational square-root enclosures, outward-rounded Machin
formula, sine/cosine Taylor remainders of orders 47 and 46, respectively,
and all interval-distance inequalities. The dedicated replay

```sh
.venv/bin/python -c 'import json; from computations.bh_boundedness_barrier_2026_09_18 import torsion_nondensity_exact; print(json.dumps(torsion_nondensity_exact(), indent=2))'
```

passed all 254 proper-subset and 210 root comparisons. The certified
proper-subsum modulus lower bound is `1/15`; the ratio is at distance
strictly greater than `1/200` from `mu_210` (diagnostic lower bound
`0.006878361991474618`). A sup-norm approximation error at most `1/400`
would preserve irreducibility and move every ratio by at most `1/200`,
giving a contradiction. Fixing additional variables preserves this
obstruction.

Finally, the asserted convex-hull gap is numerically and algebraically
correct: on the eight-point cube, sup-distance at least `1/400` implies
`L2` squared distance at least `1/1280000`, hence real inner product at
most `1-1/2560000`. Unimodularity converts uniform convex approximation
into an average `L2` squared error, contradicting this gap. This only
defeats norm-one torsion reductions; an arbitrary absolute factor is not
ruled out by that theorem.

### 6.2 Platform and scope note

During this campaign the root agent and later this agent temporarily
reported an external API error concerning an unsupported
`access_programs.cyber` parameter. A peer follow-up resumed this track;
no permission changes, altered tool schemas, or access bypass was used.
All research remained in the assigned artifact/script paths. No ledger,
steering file, or commit was changed by this track.

### 6.3 Independent audit of the genuine degree-four lower example

At 2026-09-18 05:10 UTC, I independently reconstructed Section 7 of
`artifacts/bh_boundedness_counterexamples_2026_09_18.md` and reviewed and
reran its integer Eisenstein replay. Verdict: **passes**. The result is

```math
B_4>2+2^{-18}.
```

This is a genuine new campaign lower bound for the requested unrestricted
complex constants, though it does not show divergence. No external novelty
claim is made here.

The five-variable root-valued seed has all 16 degree-at-most-two
coefficients of magnitude `1/4`. Its exact 32-value computation puts all
values among sixth roots, proving cap 1 and critical ratio 2. The local
function `a=conjugate(zeta)+zeta x_3+x_4` has degree 1, while the supplied
complete coefficient table and independent XOR convolution show that
`b=f^2 conjugate(a)` has reduced degree exactly 3.

For the two-block functions `F=f tensor f` and `h=a tensor b-b tensor a`,
the antilinear involution `J_F u=F^2 conjugate(u)` satisfies `J_F h=-h`.
Thus `Re(conjugate(F)h)=0` pointwise. The exact quotient histogram is
`0` at 640 vertices and each of `+/-4i sqrt(3)` at 192 vertices, giving
cap exactly `sqrt(1+48 epsilon^2)`.

There are 36 genuinely new coefficients of `h` outside the 256-term
support of `F`: 24 have modulus `sqrt(3)/4`, and 12 have modulus `3/4`.
Their squared mass is exactly `45/4`, a lower bound for their `8/5` mass.
On the old support, the convex tangent inequality for `|z|^(8/5)` has a
common modulus prefactor because `|Fhat|=1/16`. Its summed linear term is
zero by Parseval and the pointwise tangency above. Therefore no old
coefficient loss is hidden when the new `epsilon^(8/5)` mass is added.

For `epsilon=2^(-10)`, denominator concavity gives the increase
`3/81920`, while the rational coefficient-mass gain is at least
`33/1310768>2^(-16)`. The comparison `2^(8/5)<7/2` is certified by fifth
powers. On the interval needed for the final `5/8` power, its derivative
is greater than `1/4`. Thus the final claimed improvement by `2^(-18)`
is rigorous. The decimal `2.0000467013452567` is only a diagnostic.

I also checked the focused exact-normalization escape test suggested by
this track. The radial correction requires
`C=-h^2 conjugate(F)`, whose exact coefficient at mask 503 is `3/4` and
has weight 8. The integer replay verifies 751 nonzero coefficients and
degree exactly 8. Since `F,h` have degree 4, the top correction cannot
cancel for any nonzero perturbation. This normalization therefore does
not produce a degree-four unimodular lower example.

The later SU(2)-completion diagnostic has the same precise obstruction:
the off-diagonal entry `f P_0`, with
`P_0=(1+x_3)(1-x_4)/4`, has degree 4, while the desired diagonal entries
have degrees 1 and 3. The counterexample track certified a nonzero
degree-four coefficient `(-1+zeta)/16` at mask 15. This blocks the naive
balanced degree-two-per-block unitary completion; no general
nonexistence theorem for other completions is claimed.

### 6.4 Independent audit of the persistent asymptotic gap

At 2026-09-18 05:16 UTC, I independently checked all arguments and explicit
constants in Section 8 of the counterexample artifact. Verdict: **passes**.
The stronger scoped result is

```math
\liminf_{m\to\infty}B_m>2+2^{-207}.
```

This remains a constant lower bound and does not decide boundedness. Its
mechanism does go beyond plain tensoring: joint polynomial normalization
of many nearly unitary blocks costs much less degree than exact
normalization of each block separately.

For the exact `F,h` above, put `epsilon=2^(-104)`,
`R=|h|^2/48 in {0,1}`, `lambda=log(1+48 epsilon^2)/2`, and
`G=(F+epsilon h) exp(-lambda R)`. Then `|G|=1`, while `R` has degree at
most 8. The proof establishes

```math
H(|\widehat G|^2)\ge 8\log 2+
 \varepsilon^2[(45/4)\log(1/\varepsilon)-256].
```

I checked the finite entropy estimate without relying on an unspecified
asymptotic remainder. In `G=F+epsilon h+w`, the exact norms
`||h||_2^2=18`, `||48FR||_2^2=864` yield
`||w||_2<=24 epsilon^2(epsilon sqrt(18)+sqrt(864)/48)<16 epsilon^2`.
Thus `||G-F||_2<4.5 epsilon`, and old probabilities stay above half their
initial value `1/256`. The summed Taylor penalty is bounded by
`[7(81/4)+153^2/256] epsilon^2<256 epsilon^2`. The 36 new coefficients
have initial magnitudes between `1/4` and `3/4`, so their relative
amplitude error is at most `64 epsilon`; their squared mass loses at
most a factor `128 epsilon`. This gives the displayed logarithmic
entropy gain, with all other summands nonnegative.

On `k` independent blocks, let `J=sum R_i` and multiply
`P_k=product(F_i+epsilon h_i)` by the Taylor polynomial of `exp(-lambda J)`
having `L_k=ceil(4 lambda k)` terms. The resulting polynomial `Q_k` has
degree at most `M_k=4k+4(L_k-1)`, and

```math
\|Q_k-G^{\otimes k}\|_\infty
 \le e^{\lambda k}(e\lambda k/L_k)^{L_k}
 \le e^{-\lambda k/3}.
```

The degree bound uses a further exact reduction, independently checked
after the initial audit: `R_i^2=R_i`, so every `J^j` reduces to squarefree
products of at most `j` distinct block indicators. Moreover,
`P_i R_i=F_i R_i+epsilon h_i=C_i/48+epsilon h_i`, since `h_i R_i=h_i`.
Consequently, each selected indicator raises its block's degree from at
most 4 to at most 8, at cost 4, rather than paying an additional 8. A
squarefree term indexed by `A` therefore has degree at most `4k+4|A|`.

The approximation inequality follows from monotonicity in `L_k>=4 lambda k` and
`log 2>=2/3`. The approximation transfers to the changing critical norm:
on `10k` variables, `||Ehat||_(q_m)<=2^(5k/m)||E||_infinity`. This factor
stays bounded when `m>=M_k`, so no ambient error is discarded. Choosing
the largest `k` with `M_k<=m` handles every sufficiently large integer
`m`, since the gaps of `M_k` are bounded and `m/k -> d=4+16 lambda`.
The exact product Renyi identity then gives the lower limit
`exp(H(|Ghat|^2)/(2d))`.

Finally, `lambda<=24 epsilon^2`, `d<5`, and the entropy bound give
`H-2d log 2 >=[402 log 2-256]epsilon^2>=12 epsilon^2`, hence
`H/(2d)>log 2+epsilon^2`. Exponentiating proves the stated strict gap
`2 epsilon^2=2^(-207)`.
The exceedingly large degrees needed by this explicit tiny perturbation
do not affect the asymptotic quantifiers. A fixed ambient ratio remains,
so this construction alone cannot diverge.

## 7. Extending spin rotation to orbit quotients: exact obstruction

This section attacks the following proposed closing mechanism directly:
perhaps transitivity, or a sufficiently large permutation orbit, supplies
the same dimension gain as full permutation symmetry in Section 3. The
quotient calculation is exact, and the transitive influence inequality it
suggests is false.

### 7.1 Exact orbit-quotient commutator energy

Let a permutation group `G<=S_n` leave `f` invariant. Its vertex orbits are
`O`, with normalized orbit vectors `|O>` and probabilities
`w_O=|O|/2^n`. For adjacent distinct orbits, let `a_OP` be the number of
single-bit neighbors in `P` of any vertex in `O`. This is well defined;
edge counting gives `w_O a_OP=w_P a_PO=:e_OP`. In the orbit basis,

```math
\langle P|J_x|O\rangle=\tfrac12\sqrt{a_{OP}a_{PO}},
\qquad F|O\rangle=f_O|O\rangle.
```

The same rotation/Bernstein argument as Section 3 gives operator norm
`||[J_x,F]||<=m||f||_infinity` on this invariant subspace. Averaging its
squared columns therefore proves

```math
\mathcal E_G(f):=
\frac14\sum_{\{O,P\}}e_{OP}(a_{OP}+a_{PO})|f_O-f_P|^2
\le m^2\|f\|_\infty^2.
```

The sum is over unordered orbit pairs. In the same notation,

```math
I(f)=\frac12\sum_{\{O,P\}}e_{OP}|f_O-f_P|^2.
```

Thus an orbit edge receives an amplification factor
`(a_OP+a_PO)/2`, not the group size or orbit size. For full symmetric
orbits the sum is always `n+1`, recovering Section 3 exactly. For weak
symmetries it can remain 2. A sufficient condition
`a_OP+a_PO>=kappa` only yields `I<=2m^2||f||_infinity^2/kappa`.
Replacing `kappa` by `n` without a separate argument is invalid.

Fourier coefficients are constant on subset orbits `Omega`. If `v_Omega`
is the squared coefficient mass on such an orbit, and there are `K`
orbits intersecting levels at most `m`, the elementary coefficient step is

```math
\|\widehat f\|_{q_m}^2
\le K^{1/m}\sum_\Omega |\Omega|^{1/m}v_\Omega.
```

This identifies a second independent cost: a large number of Fourier
orbits. Transitivity on coordinate positions controls neither this number
nor the needed edge amplification. The spin argument alone does not close
either cost.

### 7.2 A leaf-transitive Boolean family defeats dimension amplification

Define the symmetric three-input Boolean polynomial

```math
E(x,y,z)=\frac{xy+xz+yz-1}{2}.
```

It is `+1` when the three inputs agree and `-1` otherwise. Put `f_0(x)=x`
and recursively apply `E` to three disjoint copies of `f_t`. Then
`f_t` is Boolean, is invariant under the leaf-transitive automorphism
group of the depth-`t` ternary tree, and has

```math
n_t=3^t,\qquad m_t=\deg f_t=2^t,\qquad \|f_t\|_\infty=1.
```

The degree equality follows inductively because products from distinct
pairs of child blocks have disjoint nonempty top Fourier supports; they
cannot cancel. Let `mu_t=E f_t` and `I_t=I(f_t)`. Direct differentiation,
using independence and Booleanity of the children, gives the exact
recurrences

```math
\mu_{t+1}=\frac{3\mu_t^2-1}{2},\qquad
I_{t+1}=\frac32(1+\mu_t^2)I_t,
\qquad \mu_0=0, I_0=1.
```

For example, a leaf derivative in the first child is its child derivative
times `(f_t(second)+f_t(third))/2`, whose squared expectation is
`(1+mu_t^2)/2`. Summing all three child contributions gives the formula.
Consequently,

```math
\frac{n_t I_t}{m_t^2}\ge(9/8)^t\longrightarrow\infty.
```

There is therefore no absolute `C` such that every transitive-invariant
bounded degree-`m` polynomial satisfies `n I(f)<=C m^2||f||_infinity^2`.
This falsifies the precise influence inequality needed to reuse the
symmetric coefficient envelope. These examples are Boolean-valued, so
their critical BH ratios remain below 2: this is a mechanism obstruction,
not evidence for divergent BH constants.

Exact replay is in
`computations/bh_boundedness_positive_2026_09_18_orbits.py`. Exhaustive
Walsh transforms and orbit-edge counts on `n=1,3,9` check the degree,
influence, and both formulas in Section 7.1. For depth two it finds 20
vertex orbits, `I=45/16`, and averaged commutator energy `117/16`.
Exact rational recurrences are additionally checked through depth 12;
the diagnostic `n I/m^2` there is approximately `15.25549`. The limiting
failure follows from the proved lower bound, not from this finite replay.
The barrier track independently reconstructed all orbit factors and tree
recurrences and confirmed the theorem during this campaign.

### 7.3 Independent audit of the fixed-block positive extension

The barrier track's Section 16 proves that symmetry within each of
`ell` disjoint variable blocks gives

```math
\|\widehat f\|_{q_m}^2
\le 2e(1+2\ell)\binom{m+\ell}{\ell}^{1/m}\|f\|_\infty^2
\le2e(1+2\ell)(\ell+1)\|f\|_\infty^2.
```

I independently checked the complete proof: the joint Dicke commutator
has distinct block transitions, giving exactly
`sum_j(n_j+1)I_j/2<=m^2||f||_infinity^2`; its weighted AGM uses weights
`r_j/m` and `1-sum_j r_j/m` with terms `sqrt(n_jr_j)/m` and the same
remaining weight. Squaring the AGM inequality gives the required
multilevel multiplicity envelope. The multilevel count, final Holder
step, and all complex-scalar and degree quantifiers are correct. Verdict:
**passes**. The dependence on `ell` is explicit and does not imply a
uniform result for unrestricted functions.

## 8. Scalarizing dense quantum Hadamard constructions does not diverge

The new primary preprint [Slote, arXiv:2608.01424](https://arxiv.org/html/2608.01424)
constructs flat Pauli unitaries using an anticommuting core and a
Hadamard-generated commuting suffix. Its Section 4 motivates a concrete
scalarization test. The following is an application of the root track's
Hadamard-chain rigidity mechanism, with the degree of the scalar core
explicitly counted.

Let `U` be a real normalized `L`-by-`L` Hadamard matrix, and let
`p_l(z)=theta_l chi_(S_l)(z)` be distinct scalar Walsh characters times
phases `|theta_l|=1`, with maximum degree `r>=1`. On `k>=2` further
independent blocks of `L` signs, put

```math
T(z,x)=L^{-1/2}p(z)^T U D(x^{(1)})U\cdots U D(x^{(k)})\mathbf1.
```

There are exactly `L^(k+1)` distinct Fourier terms, all of magnitude
`L^(-(k+1)/2)`, and the actual degree is `d=r+k`. Yet

```math
\|T\|_\infty\ge\sqrt L/2,
\qquad
\frac{\|\widehat T\|_{q_d}}{\|T\|_\infty}
\le2L^{(1-r)/(2d)}\le2.
```

To prove the cap lower bound, fix `z` and all but the last two sign
layers. The row vector before the penultimate diagonal has squared norm
`L`, since all preceding matrices are orthogonal or diagonal unitary.
Choose a global phase so its real part has squared norm at least `L/2`.
Randomize the penultimate signs. Flatness of the final Hadamard matrix
makes the variance of each final real coordinate at least `1/2`.
The real Rademacher Khintchine inequality therefore gives expected
absolute value at least `1/2` for each coordinate. Some sign choice has
real-part l1 norm at least `L/2`; choose the final signs to align these
real parts and divide by `sqrt L`. If every core phase is real, the same
argument gives cap at least `sqrt(L/2)` and ratio at most `sqrt2`.

The `k>=2` hypothesis is essential to this particular proof. No assertion
is made here about arbitrary scalar replacements, single-layer suffixes,
non-Hadamard networks, or general noncommutative-to-classical transfers.

## 9. Circle-valued structure: affine blocks and local rigidity

The unrestricted unit-circle-valued support bound `|supp fhat|<=4^m`
has not been proved here. In particular, neither the torsion argument
nor numerical searches establish it. The following two exact structural
results do apply without a torsion assumption.

### 9.1 At most `2^m` simultaneously affine data variables

Suppose `|F(x,z)|=1`, `deg F<=m`, and no Fourier monomial involves two
of the designated data variables `z_1,...,z_N`. Write

```math
F(x,z)=a_0(x)+\sum_{j=1}^N z_j a_j(x),
\qquad\deg a_j\le m-1.
```

For each fixed `x`, the coefficients of `z_j` and `z_jz_k` in `|F|^2`
show that the nonzero complex vectors among `a_0(x),a_1(x),...,a_N(x)`
are pairwise orthogonal in the real plane. At most two are nonzero.
Every nonzero complex polynomial `a_j` of degree at most `m-1` is nonzero
on a fraction at least `2^(1-m)` of its cube, by the elementary support
induction already used in Section 2. Consequently,

```math
N\,2^{1-m}\le\mathbb E_x\#\{j:a_j(x)\ne0\}\le2,
\qquad N\le2^m.
```

Only relevant data variables are counted. If `a_0(x)` never vanishes,
the right side is at most one, improving the bound to `2^(m-1)`.
The same proof for real `s`-dimensional sphere-valued outputs gives
`N<=s 2^(m-1)`. For degree two, every independent set in the Fourier
interaction graph of relevant variables therefore has size at most four.
This does not bound the whole graph: dense interactions remain the
missing case.

### 9.2 A general tangent-gap rigidity lemma

Let `f` be unimodular of degree at most `m` on a fixed base cube. Suppose
there is a positive `tau` such that every complex polynomial `h` on that
cube of degree at most `m-1` satisfies

```math
\|\operatorname{Re}(\overline f h)\|_2\ge\tau\|h\|_2.       \tag{TG}
```

Then, for any number of additional variables, every unimodular polynomial
`g` of degree at most `m` satisfying `||g-f||_infinity<tau/2` is independent
of all the additional variables.

Indeed, write `g=a(x)+r(x,z)` where `a=E_z g` and `E_z r=0`. If
`delta=||g-f||_infinity`, then `||a-f||_infinity<=delta` and
`||r||_infinity<=2delta`. Every nonconstant fresh Fourier component of
`r` has base degree at most `m-1`; orthogonality of fresh characters and
(TG) give `||Re(bar f r)||_2>=tau||r||_2`. Projecting `|g|^2-|f|^2=0`
onto nonconstant fresh characters yields

```math
2\operatorname{Re}(\overline f r)
=-2\operatorname{Re}(\overline{a-f}\,r)
  -\bigl(|r|^2-\mathbb E_z|r|^2\bigr).
```

The projection is L2-contractive, so the right side has L2 norm at most
`2delta||r||_2+||r||_infinity||r||_2<=4delta||r||_2`. Thus
`2tau||r||_2<=4delta||r||_2`, forcing `r=0` at the asserted radius.
This argument is dimension-uniform and requires no analytic path or
algebraic curve-selection assumption.

The same calculation proves a robust L2 statement if `g` is not exactly
unimodular: for `delta<tau/2`,

```math
\|g-\mathbb E_z g\|_2
\le\frac{\||g|^2-1\|_2}{2\tau-4\delta}.
```

This is only L2 control, not dimension-uniform critical-coefficient-norm
control; the distinction is important.

### 9.3 An exact uniform radius for the five-bit flat seed

Let `f` be the counterexample track's five-bit, 16-term unimodular seed,
with coefficient masks and sixth-root powers as given in Section 6.3 and
the replay script below. For every affine complex `h` on these five bits,

```math
\|\operatorname{Re}(\overline f h)\|_2^2\ge\tfrac14\|h\|_2^2.
```

Therefore **every unimodular degree-at-most-two polynomial on any larger
cube within strict uniform distance `1/4` of this seed depends only on the
original five bits**. Arbitrary global phases, input sign changes, and
permutations give the corresponding equivalent statements.

Here is a small exact certificate. Write `f=u+iv` and express affine
`h=c+id` in the ordered real Walsh basis `(1,x_1,...,x_5)`. The Gram
matrix of the real map `h -> Re(bar f h)=u c+v d` is

```math
\frac1{16}\begin{pmatrix}7I&\sqrt3 C\\\sqrt3 C&9I\end{pmatrix},
\qquad
C=\begin{pmatrix}
0&0&0&1&-1&0\\
0&0&-1&0&0&1\\
0&-1&0&0&0&1\\
1&0&0&0&-1&0\\
-1&0&0&-1&0&0\\
0&1&1&0&0&0
\end{pmatrix}.
```

The matrix `C` is symmetric and every absolute row sum is 2, so
`||C||<=2`. After subtracting `I/4`, its Schur complement is positive:
equivalently the scalar 2-by-2 lower envelope has diagonal `3,5` and
off-diagonal magnitude `2sqrt3`, with determinant `15-12>0`.
The exact smallest Gram eigenvalue is `(8-sqrt13)/16>1/4`, though the
weaker bound is enough.

`computations/bh_boundedness_positive_2026_09_18_unitary.py` independently
evaluates all 32 seed values, reconstructs this Gram matrix over the
rationals by using `h=c+i sqrt3 d`, and checks positive rational pivots
after subtracting one quarter of the corresponding norm metric. The
replay passes. No floating eigenvalue is used as a certificate.

This is local rigidity, not a classification of all degree-two unitary
functions. Other components or distant examples on six or more variables
are not excluded. The general higher-degree support/entropy target also
remains open in this investigation.

### 9.4 Every fixed-degree circle-valued polynomial is a finite junta

**Theorem.** For each integer `m>=0` there is a finite `J_m` such that,
for every ambient dimension `n`, every function
`f:{-1,1}^n -> {z in C:|z|=1}` of Walsh degree at most `m` depends on
at most `J_m` input coordinates. No restriction on the phase alphabet is
imposed. The bound here is qualitative and very large; it is not an
exponential-in-degree support bound.

Here is an explicit finite recurrence in terms of multicolor uniform
hypergraph Ramsey numbers. Set `J_0=0`, `J_1=2`. For `m>=2` let

```math
L_m=\max\{2m,\,2^{m-1}J_{m-1}+1\},\qquad
J_m=R_m(L_m,2m,2m,2m,2m,2m,2m,2m,2m)-1.
```

The notation `R_m(t_0,...,t_8)` denotes the least finite `N` such that
every nine-coloring of the `m`-subsets of an `N`-point set has a
monochromatic `t_j`-point clique in color `j` for some `j`. Existence is
the finite hypergraph Ramsey theorem.

The degree-one base case follows by writing `f=c+sum a_i x_i` and
comparing coefficients of `|f|^2=1`: the nonzero complex vectors among
`c,a_1,...,a_n` are pairwise orthogonal as real two-dimensional vectors.
There are at most two, so at most two relevant coordinates.

For the induction, discard irrelevant coordinates and suppose their
remaining number is at least `J_m+1`. Color an `m`-set `S` with color 0
if `ahat(S)=0`. Otherwise color it by the one of eight half-open argument
arcs of width `pi/4` containing `ahat(S)/|ahat(S)|`.

A `2m`-set `U` cannot be a clique of any nonzero color. Indeed, because
`deg f<=m`, the coefficient of `|f|^2` at `U` has exactly the terms

```math
0=\widehat{|f|^2}(U)
 =\sum_{\substack{S\subset U\\|S|=m}}
       \widehat f(S)\overline{\widehat f(U\setminus S)}.
```

To see why there are no omitted terms, `S triangle T=U` with
`|S|,|T|<=m` forces `|S|=|T|=m`, `S intersect T` empty, and
`S union T=U`; outside coordinates cannot occur. Pair the displayed terms
under complementation. Within a common argument arc every resulting
real part is strictly positive, contradicting the equality to zero.

Ramsey therefore supplies a set `A` of size `L_m` all of whose `m`-subset
coefficients vanish. Fix any assignment `y` to the coordinates outside
`A`. The restricted polynomial `f_y` is still unimodular and has degree
at most `m-1` in `A`: a monomial containing `m` coordinates of `A` could
have no outside factor, and its coefficient was zero. By induction,
each `f_y` has at most `J_(m-1)` relevant coordinates.

On the other hand, each original coordinate `i in A` is relevant. Choose
one nonzero coefficient `ahat(T)` with `i in T`, and put `B=T intersect A`.
The coefficient of `chi_B` in `f_y` is the nonzero outside polynomial

```math
c_B(y)=\sum_{C\subseteq A^c}\widehat f(B\cup C)\chi_C(y),
\qquad\deg c_B\le m-|B|\le m-1.
```

It is not identically zero because its `C=T\setminus A` coefficient is
nonzero. Hence it is nonzero on at least a fraction `2^(1-m)` of outside
assignments. On all those assignments, `i` remains relevant to `f_y`.
Linearity of expectation, with no independence assumption between these
events, gives

```math
J_{m-1}\ge\mathbb E_y\#\{i\in A:i\text{ relevant to }f_y\}
\ge |A|2^{1-m}>J_{m-1},
```

the desired contradiction.

The root and barrier tracks independently reconstructed and adversarially
checked every step, including the top-degree autocorrelation and the
survival probability after restriction. This theorem implies a finite
dimension-independent Fourier support bound at each fixed degree, but
the Ramsey growth is much too large to imply `C^m` support or linear
entropy bounds. Its direct campaign consequence is the all-unimodular
atom obstruction in Section 2.8.

The same argument works for polynomials valued in any fixed-dimensional
Euclidean unit sphere: cover the sphere of coefficient directions by
finitely many sets with positive pairwise inner products, and replace
the eight phase colors by those direction colors. The junta bound may
depend on the output dimension. Thus this does not settle a quantum junta
question whose matrix dimension grows with the input dimension.

Targeted primary-literature searches for sphere-valued low-degree juntas,
unimodular polynomial Ramsey arguments, and unit-modulus juntas did not
identify a prior statement during this campaign. That limited search is
not a novelty certificate, and no external priority claim is made.

## 10. Quantitative entropy target and precise failures of shortcuts

The strongest remaining targeted inequality in this branch is

```math
H_f=-\sum_S|\widehat f(S)|^2\log|\widehat f(S)|^2
\stackrel{?}{\le} C v(f),\qquad
v(f)=\|F^*J_xF-J_x\|,
```

for arbitrary unit-circle-valued `f`. Here all logarithms are natural.
The generic rotation Bernstein inequality gives `v(f)<=deg f`, so this
would in particular imply linear Fourier entropy in degree. It remains
unproved. The five-variable flat seed has entropy `log 16=4 log 2`, not
`8 log 2`; its numerical speed is approximately `1.8660061245` and the
ratio is approximately `1.48584117`.

The diagnostic script
`computations/bh_boundedness_positive_2026_09_18_entropy.py` optimizes
`H_f/v(f)` over unrestricted phase vectors, using analytic derivatives
checked against finite differences. Eight five-variable starts reached
ratios between `1.5985` and `1.7323`; six six-variable starts reached
ratios between `1.7046` and `1.7830`; four seven-variable starts reached
`1.764175052, 1.816253190, 1.810064153, 1.801808921`. All four of the
seven-variable starts stopped at their 800-iteration limit. Most of the
smaller starts also stopped at the iteration limit,
and the largest-eigenvalue objective is nonsmooth at eigenvalue
collisions. These are neither certified maxima nor low-degree examples,
and no upper bound is inferred from them.

### 10.1 Fourier marginal entropies cannot close the speed inequality

The always-valid subadditivity estimate
`H_f<=sum_i h(I_i(f))` cannot be followed by a dimension-free comparison
of that sum with `v(f)`, even for homogeneous Boolean functions.

Let

```math
g(x)=\frac{x_1x_3+x_1x_4+x_2x_3-x_2x_4}{2},
```

and recursively compose copies on four disjoint blocks, starting with
`g_0(x)=x`. The depth-`t` function is Boolean, homogeneous of degree
`m_t=2^t`, uses `n_t=4^t` coordinates, and every coordinate influence is
`2^-t`. Because a homogeneous Boolean degree-`m_t` function has exactly
`m_t` sensitive neighbors at every vertex, its sensitivity graph is
`m_t`-regular; hence `v(g_t)=m_t`. Consequently

```math
\frac{\sum_i h(I_i(g_t))}{v(g_t)}
=2^t h(2^{-t})=t\log2+1+o(1)\longrightarrow\infty.
```

The actual Fourier entropy does not diverge relative to degree. Distinct
ancestries cannot collide because all inner functions are balanced.
The coefficients are flat on `4^(2^t-1)` sets, giving exactly
`H(g_t)=(2^t-1)log4`. Thus the failed shortcut loses crucial correlations
between Fourier-coordinate indicators; it is not a counterexample to
the targeted entropy inequality.

### 10.2 Scalar quantum walks: exact identification, no imported bound

The Fourier convolution operator
`U=sum_S fhat(S) translation_S` on `Z_2^n` is unitary exactly when `f`
is unimodular, and `deg f<=m` is propagation radius at most `m`.
This identifies the problem with homogeneous scalar quantum walks, but
the primary classifications inspected do not supply the missing support
or entropy estimate.

Lopez Acevedo--Roland--Cerf,
[Exploring scalar quantum walks on Cayley graphs](https://arxiv.org/html/quant-ph/0609234),
Sections II--IV, classify two/three-generator cases and selected cyclic
and Johnson-graph examples. Their basic quadrangularity condition is
automatic for our exponent-two group: `(S,T)` and `(T,S)` have the same
difference. A triangle of uniquely represented unordered differences
would force three pairwise real-orthogonal nonzero complex coefficients,
which is impossible, but this only recovers a local obstruction and does
not control a full Hamming ball's many representations.

Bisio--D'Ariano--Erba--Perinotti--Tosini,
[Quantum walks with a one-dimensional coin](https://arxiv.org/html/1603.07666),
Proposition 1, treat `F x Z^d` with `d>=1`. They diagonalize the finite
factor and classify the infinite directions. Arbitrary phases on the
finite factor remain, so the theorem does not constrain our finite-cube
Fourier support.

Bu--Garcia--Jaffe--Koh--Li,
[Complexity of quantum circuits via sensitivity, magic, and coherence](https://arxiv.org/html/2204.12051),
Theorem 23, give a quantum entropy/influence estimate with a logarithmic
dimension factor. It is not the scalar entropy/speed inequality above.
Their general quantum FEI conjecture also cannot be imported as a
theorem; the noncommutative examples discussed in Section 8 refute it.

## 11. Sharp support bound for bipartite quadratic unitaries

**Theorem.** If `f` is a complex unit-circle-valued degree-at-most-two
polynomial and its quadratic interaction graph is bipartite, then

```math
|\operatorname{supp}\widehat f|\le8,
\qquad \|\widehat f\|_{4/3}\le 2^{3/4}.
```

Both bounds are sharp. Linear and constant terms are allowed, and the
number of original variables is unrestricted.

Choose a bipartition of the variables, putting isolated vertices on
either side. The polynomial is separately affine in the two blocks.
With one additional sign variable on each side, homogenize it as

```math
B(x_0,x,y_0,y)=x_0y_0 f(x_0x,y_0y)=X^T A Y.
```

This preserves unimodularity and the number of nonzero Fourier
coefficients. It suffices to bound the number of nonzero entries of `A`.

Fix `Y`. The homogeneous linear polynomial `X^T AY` has constant modulus
one. Comparing its squared-modulus coefficients shows that its nonzero
coefficients are pairwise orthogonal as real vectors in `C=R^2`.
Thus at most two row sums `(AY)_i` are nonzero at every `Y`.

Every nonzero homogeneous complex linear form on a sign cube is nonzero
on at least half the cube, by conditioning on one nonzero coefficient.
Hence there are at most four nonzero rows of `A`; by symmetry there are
at most four nonzero columns.

We use two elementary sharpenings of that half-cube fact. First, equality
of the zero probability with `1/2` forces exactly two nonzero
coefficients, equal or opposite. Indeed, if any two coefficients are
neither equal nor opposite, their four signed sums are distinct, and
conditioning on the other signs gives zero probability at most `1/4`.
Otherwise all nonzero coefficients are equal up to signs, and the
central binomial probability is `1/2` only for two coefficients.
Second, a homogeneous linear form with exactly three nonzero complex
coefficients has zero probability at most `1/4`: either the preceding
distinct-sums argument applies, or all three are equal up to signs and
zero is impossible.

If there are four rows, their nonzero probabilities sum to at most two
and each is at least `1/2`. Thus every row has exactly two entries, and
there are eight entries total. The same applies if there are four
columns. In the remaining case the matrix has at most three rows and
three columns. More than eight entries would require a full `3 x 3`
matrix; but its three row forms are each nonzero with probability at
least `3/4`, giving expected active-row count at least `9/4`, a
contradiction. This proves the support bound. Parseval and the elementary
support-size norm comparison give the stated critical norm bound.

For sharpness on six variables, take

```math
B=\frac{(x_1+x_2)(y_1+i y_2)
             +(x_1-x_2)(y_3+i y_4)}{2\sqrt2}.
```

If `x_1=x_2`, only the first branch remains, and otherwise only the
second remains; in either case the modulus is one. Its eight homogeneous
quadratic coefficients all have modulus `1/sqrt8`, proving sharpness.
The Gaussian-integer replay
`computations/bh_boundedness_positive_2026_09_18_bilinear.py` checks
`|sqrt8 B|^2=8` by exact Walsh convolution and at all 64 vertices.

The barrier track independently audited the complete argument. The
sixth-root 16-coefficient extremizer from Section 6 is therefore
necessarily outside this bipartite class. No general quadratic
classification or all-degree support theorem is claimed.

## 12. Sharp all-degree theorem for unitary block-affine polynomials

**Theorem.** Let `k>=1`, let the Boolean coordinates be partitioned into
`k` arbitrary disjoint blocks, and let `f` be complex-valued, unimodular,
and affine separately in each block. Then

```math
|\operatorname{supp}\widehat f|\le 2^{2k-1},
\qquad
\|\widehat f\|_{2k/(k+1)}\le 2^{1-1/(2k)}.
```

Both bounds are sharp for every `k`. This is a genuine all-degree
dimension-uniform theorem for the stated unitary subclass, not for all
bounded block-multilinear polynomials or all unimodular degree-`k`
polynomials.

### 12.1 Odd Boolean postprocessing preserves every affine block

Fix all coordinates outside one block. The resulting affine function is
`c+sum_i a_i x_i` and has constant modulus one. As in the degree-one base
case of Section 9, all nonzero vectors among `c,a_1,...,a_n` are pairwise
orthogonal over the reals. In particular, there are at most two relevant
variables. If there are two, the constant term is zero and the slice is
odd under simultaneous negation of those two variables.

Let `phi` be ANY odd map from `range(f) union -range(f)` to `{-1,1}`.
On a slice depending on zero or one variable, `phi(f)` is automatically
affine. On a slice depending on two variables, it is an odd function of
those two signs. Its constant and quadratic coefficients vanish, so it
is again affine. This holds for every outside restriction; consequently
the global Fourier expansion of `phi(f)` has at most one coordinate in
each block. Thus every such odd Boolean postprocessing has degree at
most `k`.

The active variables of a slice may depend on the outside restriction.
This causes no gap: a Fourier coefficient containing two coordinates
from that block is the outside average of a slice coefficient which
vanishes for every restriction.

### 12.2 A parity lattice converts postprocessing into sharp sparsity

Partition `range(f) union -range(f)` into antipodal pairs and choose one
unit representative `rho_j` from each pair. Define the real signed
indicator `g_j` to equal `1` where `f=rho_j`, `-1` where `f=-rho_j`, and
zero elsewhere. A pair may have only one of its values actually
observed; the definition still applies. Pointwise,

```math
f=\sum_j\rho_jg_j,\qquad \sum_jg_j^2=1.
```

For every independent choice of signs `sigma_j`, the function
`h_sigma=sum_j sigma_j g_j` is an odd Boolean postprocessing of `f`.
It therefore has degree at most `k`. The usual Boolean Fourier
granularity says

```math
\widehat h_\sigma(S)\in 2^{1-k}\mathbb Z
\quad\text{for every }S\text{ and every }\sigma.
```

For completeness, write a Boolean function in the `0,1` variables
`y_i=(1+x_i)/2`. Its multilinear coefficients are integer finite
differences of its values. All nonconstant coefficients are even,
because every value is odd; only terms of degree at most `k` occur.
Expanding each product of `y_i` back into signs proves the displayed
granularity, including the constant coefficient when `k>=1`.

Fix `S` and put `c_j=ghat_j(S)`. Comparing two sign vectors differing
only at `j` shows that `w_j=2^k c_j` is an integer. Taking all signs
positive further shows that `sum_j w_j` is even. A nonzero integer vector
with even coordinate sum has squared Euclidean norm at least two.
Therefore every Fourier position at which at least one `ghat_j(S)` is
nonzero contributes at least `2*4^-k` to `sum_j |ghat_j(S)|^2`.

Parseval and the signed partition identity give

```math
\sum_S\sum_j|\widehat g_j(S)|^2
=\sum_j\mathbb E g_j^2=1.
```

There are consequently at most `4^k/2` such positions. The Fourier
support of `f=sum rho_j g_j` is contained in their union, proving the
sharp support bound. Since `||f||_2=1`, the critical norm bound follows
by the support-size comparison.

An alternative, weaker proof is now transparent: odd phase quantizations
of `f` into even-order roots of unity preserve all affine blocks and
converge uniformly to `f`. The torsion theorem would give support at
most `4^k`. The parity lattice above gains the exact factor two and
needs no algebraic-number argument.

### 12.3 Sharp examples and scope

Let `r=k-1`, `N=2^r`, and let `g` be the Boolean address function with
`r` address bits and `N` independent data bits. Its `N^2=4^(k-1)` Fourier
coefficients all have magnitude `1/N`. Put each address bit in its own
block and all data bits in the final block. Let `h` be an independent
copy, grouping its corresponding blocks with those of `g`. Then

```math
f=(g+i h)/\sqrt2
```

is unimodular and affine in those `k` blocks, has actual degree `k`, and
is flat on `2*4^(k-1)` Fourier coefficients. It attains the asserted
critical ratio exactly. The exact Gaussian-integer script
`computations/bh_boundedness_positive_2026_09_18_blockunit.py` verifies
the construction for `k=1,...,5` by Walsh convolution, including the
block-affine support condition.

The barrier researcher independently checked the complete odd-map,
granularity, even-lattice, Parseval and sharpness arguments: PASS.
The bipartite quadratic theorem in Section 11 is the case `k=2` and has
an independent elementary matrix proof.

The related real-valued block-rounding argument of Arunachalam,
Chakraborty, Koucky, Saurabh and de Wolf,
[Improved bounds on Fourier entropy and Min-entropy, v2](https://arxiv.org/html/1809.09819v2),
Section 5.1, Lemma 5.7, says that a real `d`-block polynomial uniformly
within `1/3` of a Boolean function forces that Boolean function to have
degree at most `d`. Their Theorem 5.5 rules out flat Boolean approximants
of superexponential sparsity. The complete rounding proof was read and
reread during the bounded final literature check. It does not state the
sharp complex circle-valued support bound above. The exact odd-map and
even-lattice argument was derived directly here; no external novelty
claim is made.

This route does not extend just by coloring an arbitrary interaction
hypergraph into `O(deg f)` affine blocks. For example, recursive ternary
equality from Section 7 has degree `2^t` but its co-occurrence graph is
complete on `3^t` vertices: distinct top child groups are joined at every
level. Thus the required block count can exceed degree by an unbounded
factor even for Boolean-valued functions, whose BH ratios are already
uniformly bounded for a different reason.
To check the within-child edges in that recursion as well, a coefficient
supported in just one child is multiplied by the other children's mean
`mu`; the means satisfy `mu_1=-1/2` and
`mu_(t+1)=(3mu_t^2-1)/2`, and are nonzero for every `t>=1` because they
are rational and cannot equal `+/-1/sqrt(3)`. Cross-child edges come from
the unique product of those two nonempty child coefficients.

### 12.4 An exact bound on the number of output phases

The same argument gives at most `2^k` antipodal output pairs, hence at
most `2^(k+1)` distinct values. Indeed, each nonzero signed indicator
`g_j` is a difference of two degree-`k` Boolean postprocessings divided
by two, so has degree at most `k`. The elementary polynomial support
lemma gives `Pr[g_j != 0] >= 2^-k`. These events partition the cube.
This output bound is sharp: take a product of `k` two-variable affine
unitaries with values `+/- exp(+/- i theta_j)`, choosing the angles so
that all `2^(k+1)` resulting phases are distinct.

### 12.5 Odd postprocessing does not preserve general Walsh degree

The degree-preserving property in Section 12.1 needs the affine-block
hypothesis. For an explicit failure, use the eight unit complex numbers

```math
z=\left(1,1,-\frac15\pm\frac{2i\sqrt6}{5},
 -\frac13\pm\frac{2i\sqrt2}{3},
 -\frac7{15}\pm\frac{4i\sqrt{11}}{15}\right).
```

Their sum is zero and no two are antipodal. Assign them to the eight
vertices arbitrarily and put `f(x)=x_1x_2x_3 z_x`. Then `|f|=1` and its
only potentially degree-three Fourier coefficient is zero, so `deg f<=2`.
Define the odd map on `range(f) union -range(f)` by `phi(z_x)=1` and
`phi(-z_x)=-1`; absence of antipodal pairs makes this consistent. Then
`phi(f(x))=x_1x_2x_3`, of degree three. These are the same phases used in
the barrier researcher's Section 8.1, but this conclusion needs only
their zero sum and absence of antipodal pairs, not the Mann argument.

Consequently the proof of Section 12 does not establish the conjectural
support bound `|supp fhat|<=4^m` for arbitrary degree-`m` unimodular
polynomials. That broader support statement remains unproved here.

The final targeted primary-source search did not locate a resolution of
that broader arbitrary-phase support question; this is a limited search
report, not a literature-wide priority or open-problem certification.
Terminology matters: Volberg's
[unimodular-coefficient Sidon paper](https://arxiv.org/abs/2205.04936)
and the new [cb-Hadamard rigidity paper](https://arxiv.org/abs/2609.16329)
use unimodularity for coefficients or coefficient paths, not the scalar
condition `|f(x)|=1` at every Boolean input. Their statements therefore
cannot be imported as a classification of the present unitary-valued
polynomials.

The barrier researcher supplied an even stronger example using the
five-bit flat sixth-root seed: its odd Boolean postprocessing `f^3`
has degree five. Independent exact Eisenstein convolution here gives
all 32 Fourier positions nonzero and top coefficient `9/16` at mask 31.
Thus even this particular torsion seed cannot be Booleanized while
preserving its degree two.

## 13. Dimension-uniform stability for one affine block

Here is a quantitative extension which does not assume exact constant
modulus. If a complex affine polynomial `p` satisfies

```math
\big\|\,|p|^2-1\,\big\|_\infty\le\eta\le\tfrac14,
```

then there is a complex affine unimodular polynomial `u` with

```math
\|p-u\|_\infty\le9\sqrt\eta.
```

The constant is independent of the number of variables. The square-root
order cannot be improved uniformly in dimension.

### 13.1 A direct proof

Homogenize the constant coefficient with a new sign `x_0`, obtaining
`z=sum_j a_j x_j`. Its values, up to sign, are those of `p`, so the same
annulus condition holds. Write `S=sum_j |a_j|` and
`M=||z||_infty <= sqrt(1+eta)`. Averaging the maximum real projection over
rotations of the complex plane gives `(2/pi)S <= M`. Parseval gives
`sum |a_j|^2 >= 1-eta`; hence the largest coefficient has magnitude

```math
r\ge\frac{1-\eta}{S}\ge\frac14.
```

Rotate that coefficient to the positive real number `r`, label its sign
`x_1`, and write the remaining affine form as `U+iV`, where `U,V` are
real homogeneous linear forms in the remaining signs. Comparing the
two choices of `x_1` gives

```math
4r|U|\le2\eta,\qquad \|U\|_\infty\le2\eta.
```

Averaging those two squared moduli yields
`|r^2+U^2+V^2-1| <= eta`. Therefore
`|V^2-(1-r^2)| <= eta+4 eta^2 <= 2 eta`, and the oscillation of `V^2`
on its cube is at most `4 eta`.

Put `W=sum |v_j|`, where `V=sum v_j x_j`, and let `a=max |v_j|`.
If `a>W/2`, the minimum possible absolute value of `V` is exactly
`2a-W`. Its maximum is `W`, so the oscillation of `V^2` is
`4a(W-a) >= 4(W-a)^2`; consequently `W-a <= sqrt(eta)`.
If `a<=W/2`, greedy choice of signs gives a value of absolute value at
most `a`. Thus the oscillation is at least `W^2-a^2 >= 3W^2/4`, and
`W <= 4 sqrt(eta)/sqrt(3) < 3 sqrt(eta)`. The case `V=0` is automatic.
In either case, retaining one largest imaginary coefficient discards
at most `3 sqrt(eta)` in supremum norm.

The retained form `u_0=r x_1+i v_j x_j` has constant modulus `R_0>0`,
and

```math
\|z-u_0\|_\infty\le2\eta+3\sqrt\eta\le4\sqrt\eta.
```

Also `|| |z|-1 ||_infty <= eta`. It follows that
`|R_0-1| <= 4 sqrt(eta)+eta`. Normalize `u=u_0/R_0` and undo the scalar
rotation. The total error is at most `8 sqrt(eta)+eta <= 9 sqrt(eta)`.
Restricting the homogenizing sign to one proves the affine assertion.

### 13.2 The square-root rate is necessary

For `t>0`, take

```math
p_n(x)=1+\frac{it}{n}\sum_{j=1}^n x_j.
```

Then `|| |p_n|^2-1 ||_infty=t^2`. Every exact affine unimodular function
depends on at most two signs. Fix those signs and set all the others
first to one and then to minus one. The two values of `p_n` differ by
at least `2t(n-2)/n`, while that exact function is unchanged. Hence its
uniform distance from `p_n` is at least `t(1-2/n)`.

For multiple affine blocks, applying this argument to slices does not
yet produce a global exact unitary: the retained coordinates and phases
may vary with the other blocks. No dimension-uniform multi-block
stability theorem is claimed here.

The barrier researcher independently reconstructed Sections 13.1--13.2
and the output-phase bound in Section 12.4: PASS, including the greedy
discrepancy step, normalization cost, and matching square-root obstruction.

## 14. Final adversarial audit: the genuine BH lower bounds

At approximately 06:37 UTC, this track reread and reconstructed the counterexample
track's Sections 5, 7 and 8 afresh. The most important campaign result
is a genuine constant gap above two, NOT divergence:

```math
B_4>2+2^{-18},\qquad
\liminf_{m\to\infty}B_m>2+2^{-207}.
```

The finite seed and proof pass the following checks independently of
the earlier audits in Section 6:

1. The five-variable phase table has exactly sixteen degree-at-most-two
   coefficients and all 32 values are sixth roots. The transfer partner
   `b=f^2 bar(a)` has reduced degree three, not the nominal degree five.
   Thus both the two-block product `F` and tangent `h` have degree four.
2. The quotient `h/F` takes exactly `0,+4i sqrt(3),-4i sqrt(3)` with
   multiplicities `640,192,192`. Consequently the cap is exactly
   `sqrt(1+48 epsilon^2)`. The 36 coefficients outside the old support
   have squared mass exactly `45/4`; no numerical threshold is used.
3. Convexity retains the complete old-support critical mass because its
   first variation vanishes exactly by Parseval. At `epsilon=2^-10`,
   the rational new-mass surplus is at least `33/(16*81923)>2^-16`.
   Taking the `5/8` power yields the strict ratio improvement `2^-18`.
4. For the normalized unitary `G=(F+epsilon h)exp(-lambda R)`, where
   `R=|h|^2/48` and `lambda=log(1+48epsilon^2)/2`, the correction has
   coefficient `ell_2` norm at most `16epsilon^2`. On the old support,
   the entropy Taylor loss is strictly below `256epsilon^2`; on the 36
   new positions, the entropy gain is at least
   `(45/4)epsilon^2 log(1/epsilon)`. These are finite inequalities for
   every `0<epsilon<=2^-10`, not an unspecified asymptotic expansion.
5. The joint normalization polynomial has degree at most
   `M_k=4k+4(ceil(4lambda k)-1)`: after reducing indicator powers, a
   selected block costs degree eight instead of four. Its uniform error
   is at most `exp(-lambda k/3)`, from the ordinary negative-exponential
   Taylor remainder. The final coefficient error is at most this error
   times `2^(5k/m)`, bounded for `m>=M_k`.
6. Bounded gaps of `M_k` justify every sufficiently large degree, not
   just a subsequence. At the fixed explicit `epsilon=2^-104`, the
   entropy margin is `H(G)-2(4+16lambda)log2>=12epsilon^2`, which gives
   the asserted strict asymptotic gap by the exact changing-exponent
   Renyi identity.

The newly discovered failure of a velocity lower bound for cap-bounded
`L2` approximation does not affect this theorem: its normalization
converges in the uniform norm and its full coefficient error is paid
explicitly. Ordinary tensoring is not the construction; the joint
normalization uses additional degree at rate `16lambda`. The fixed
ambient-to-degree ratio still imposes a constant ceiling, so this is
not evidence establishing divergence or resolving uniform boundedness.

The exact tangent replay and the normalization replay were rerun at this
checkpoint: PASS. The latter checks eleven rational inequalities; its
500-digit entropy values are diagnostics only, not proof premises.

## 15. Final reciprocal checks of approximation scope

These are audits of other tracks' scoped statements, not additional
claims resolving BH boundedness.

- Counterexample track Section 17: independent full audit PASS for
  `H(|Ghat|^2)<=64 d_eff` for a fixed independent-phase target and
  cap-bounded tensor approximants with `limsup L2 error<1`. Conditional
  fixing of outside coordinates preserves cap and improves squared
  error. The correlation identity evaluates the restricted polynomial
  at `-i tan(theta_j)`. The reciprocal Joukowski-polynomial argument
  proves Bernstein--Walsh for complex coefficients, and the finite
  logarithmic correction tends to zero under the stated hypotheses.
  The exact dyadic sum is `(pi/2)[1+2log(8/pi)]<5`; large angles are paid
  by ordinary Fourier-degree concentration. The 400-case numerical
  evaluation replay passed, but the audit is analytic. No unstated
  simultaneous moving-seed limit is used.
- Counterexample track Section 12: the four-index moment expansion gives
  `(69/8)||K||_F^2+(3/8)||K1||_2^2` exactly. The two-phase cut and random
  cut argument yields cap at least the sum of row norms. Polynomial
  approximation of the bounded radial phase gives nested Fourier
  supports with geometric tail and logarithmic cardinality, establishing
  entropy `O_c(log N)` without ambient-dimension continuity loss. The
  spectral degree lower bound remains uniform over any number of tensor
  copies, by bounding the expected fraction of low-degree blocks. The
  stronger exact address-parity restriction in Section 12.3 also passes.
- Barrier track Sections 28--28.1: independent full audit PASS for the
  controlled disjoint-Boolean phase entropy bounds. The biased-Boolean
  spectral mixture is exact even at its empty coefficient, because the
  real and imaginary contributions have zero cross term. Parseval gives
  a common data marginal for the two classical joint laws, validating
  the conditional-entropy combination. The rare-event bounds use a
  shared control partition for every control, including the scalar
  offset; the stated shared-decision-tree extension is correct.
- Barrier track Section 23: the exact eight-by-eight matrix replay passed
  `A^3=3A`, `tr(A^2)=12`, and the displayed characteristic polynomial.
  Identifying two coordinates produces parity on two signs, proving
  uniform tensor degree rate two despite velocity `sqrt(3)`.
- Barrier track Section 29: the cap-bounded `L2` counterexample was
  checked against [de Wolf's primary paper](https://arxiv.org/html/0802.1816),
  Appendix A and its polynomial-method statement. Exact search uses
  `ceil((pi/4)sqrt(k/j))` queries as stated there; adding one verification
  query per round gives the safe budget below `4sqrt(kt)`. The invariant
  for success whenever the number of marked inputs is at most `t`,
  rare-AND substitution, Chernoff estimate, and rate `1/sqrt(2)<4` all
  check. Thus uniform approximation cannot be replaced by `L2`
  convergence in the velocity lower bound, even with cap one.

At the same final-check checkpoint, this track reran its own exact
block-unitary sharpness certificates for `k=1,...,5`, affine tangent-gap
certificate, orbit-quotient certificates, and bipartite sharpness
certificate: all PASS. The symmetry and sine-rounding scripts are
diagnostics, whereas their corresponding theorems have analytic proofs
and independent audits.

## 16. Independent audit of the flat-diagonal compiler certificate

The director's Section 9 and the barrier researcher's Sections 31--32
were reconstructed independently here: PASS. This is a useful stronger
positive closure statement discovered near the end of the campaign.

For a unitary scalar seed `G` on `n` signs, define `V(G)` as the supremum
of the absolute expectation of `M_G^* J M_G-J`, with
`J=(1/2)sum_i X_i`, over positive trace-one matrices whose physical cube
diagonal is exactly uniform. If `Q` on `nk` signs has degree at most `D`,
cap `C>0`, and `delta=||Q-G^[tensor k]||_2`, the exact finite statement is

```math
\frac Dk\ge\frac{V(G)-n\delta}{C}.
```

The two error terms really do each cost at most `nk delta/2`. For the
second, `Tr rho U^* E K`, use state Cauchy--Schwarz with `A=E^*U`,
`B=K`. Diagonality and commutation of `E,U` give `A^*A=|E|^2`, so the
flat diagonal yields the ordinary uniform-cube squared `L2` error.
No commutation of `rho` with `K` is assumed. The first error uses
`A=KU`, `B=E`. The main term is bounded by the operator commutator
norm, at most `D C` by simultaneous-rotation Bernstein. Tensoring the
seed density gives exactly `k` times the seed expectation.

The quadratic-phase entropy proof has admissible witnesses, not just
arbitrary operator-norm witnesses. The constant X eigenvector is flat.
The other states have Bloch `Z_j=z_j/2`, `X_j=0`, and selected signs of
`Y_j=sqrt(3)/2`; averaging over independent uniform `z_j` makes their
physical diagonal exactly uniform even when the selected Y signs depend
on the entire vector `z`. The linear-term extension's fictitious sign
preserves the needed absolute first and second moments. Its small-row
bound and large-row influence bound give the fully general quadratic
statement

```math
H\bigl(|\widehat {e^{iP}}|^2\bigr)\le64V(e^{iP})
\qquad(P\text{ real},\ \deg P\le2).
```

### 16.1 The actual ratios of fixed-seed compilers are bounded

The director's strengthening is also valid, and is stronger than merely
bounding the entropy-transfer lower-bound criterion. Fix a NONCONSTANT
quadratic phase `G` on exactly `n` signs. Let `Q_k` be on exactly `nk`
signs, have actual degree `D_k`, and satisfy
`||Q_k-G^[tensor k]||_2->0` and `||Q_k||_infty->1`. Since
`V(G)>=I(G)>0`, the finite inequality forces
`liminf D_k/k>=V(G)`. Thus

```math
\|\widehat {Q_k-G^{\otimes k}}\|_{q_{D_k}}
\le2^{nk/(2D_k)}\|Q_k-G^{\otimes k}\|_2\longrightarrow0.
```

The factor is bounded because the seed is fixed and the degree rate is
bounded away from zero. Along any subsequence with `D_k/k->d`, the
finite-seed Renyi entropy tends to Shannon entropy, and the ACTUAL ratio
therefore tends to `exp(H(G)/(2d))`. Actual degrees are at most `nk`, so
all subsequences have further subsequences with
`d in [V(G),n]`. Consequently

```math
\limsup_k\frac{\|\widehat Q_k\|_{q_{D_k}}}
                   {\|Q_k\|_\infty}\le e^{32}.
```

This does not cover constant targets, arbitrary auxiliary dimensions,
or moving seeds without quantitative error control. Those exclusions
are substantive: the coefficient-error factor would otherwise not be
uniformly bounded. It also does not cover arbitrary BH polynomials not
arising as these fixed-seed tensor approximants. The counterexample
researcher independently checked this actual-ratio passage as well.

For arbitrary unitary seeds the same actual-ratio passage holds along
positive degree-rate subsequences, but the entropy bound by `64V` is
special to quadratic phases. Indeed, for Boolean-valued seeds one has
`V=I` exactly, so a universal entropy bound by `V` would include the
unresolved Boolean Fourier entropy/influence conjecture.

### 16.2 Closing synthesis scope audit

The closing synthesis was read in full at approximately 06:57--07:00 UTC.
Its sharp block-affine and symmetric statements, all-unitary junta scope,
flat-diagonal inequality, and fixed-quadratic actual-ratio corollary pass.
The requested precision changes distinguish the normalized input unit
ball in the atomic obstruction, require `C>0` before division, and refer
to fixed quadratic-phase compiler families rather than claiming a
single-polynomial critical-Renyi bound for all quadratic phases. The
thirteen-record closing replay report was also inspected directly:
all return codes are zero. No unrestricted boundedness or divergence
claim is justified or made.

## Sources and provenance

Primary references inspected:

- Slote--Volberg, [Polynomial Bohnenblust--Hille bounds for product of cyclic
  groups, v2](https://arxiv.org/html/2609.07758v2), for current problem context.
- Arunachalam--Dutt--Escudero Gutierrez--Palazuelos,
  [A cb-Bohnenblust--Hille inequality with constant one and its applications
  in learning theory](https://link.springer.com/article/10.1007/s00208-025-03142-5),
  Proposition 1.5 for the already-known sharp Boolean-valued endpoint.

The stochastic sine-rounding obstruction and its detailed calculations
above were derived in this campaign. No external novelty claim is made.
