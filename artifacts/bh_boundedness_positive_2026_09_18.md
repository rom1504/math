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
tracks. All conclusions below are independent of those tracks.

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
