# Own-free transport and Gaussian Sobolev approximation

Date: 2026-09-05. Status: independent proof of the transport and approximation
steps. Application to injective rooted-tree fields additionally uses their
moment theorem, audited separately. No recursive AMP assertion is assumed.

## 1. Hypotheses and notation

Let `A` be symmetric hollow with off-diagonal entries `±1`, put `m=n-1` and
`B=A/sqrt(m)`, and let `S` be uniform independent signs. Let
`U_i(S)∈R^d` be a fixed finite collection of polynomial fields satisfying:

1. `U_i` is independent of `S_i` exactly.
2. Their Walsh degrees are bounded by a constant independent of `n`, and
   `sup_i E|U_i|²=O(1)`.
3. For each fixed finite `p`,
   `sup_(i≠k)||D_kU_i||_p≤C_p/sqrt(n)`, where `D_k` is the central half
   difference with respect to `S_k`.
4. Uniformly over roots, `U_i` converges in law to a fixed standard Gaussian
   vector `Z∈R^d`.

The fixed-degree hypothesis and hypercontractivity give uniform moments of
every fixed order. Thus hypothesis 4 also gives uniform moment convergence
for continuous functions of polynomial growth. A fixed nondegenerate
Gaussian covariance can equivalently be standardized.

Injective rooted-tree fields have exact own-freeness. Their uniform
influence estimate follows directly from coefficient counting: if a tree
has `D` nonroot vertices, each degree-`D` Walsh coefficient is bounded by
`C n^(-D/2)`, each fixed spin belongs to at most `O(n^(D-1))` supports, and
orthogonality gives `||D_kU_i||₂²=O(1/n)`. Fixed-degree hypercontractivity
upgrades this to every fixed `L^p` norm. Multiple fixed tree shapes only
change constants.

## 2. Exact two-endpoint identity

Let `g_i=g(U_i)` for any measurable `g` with the required second moments.
Since `g_j` is independent of `S_j`, and `g_k` is independent of `S_k`,
conditioning on all spins except those two gives

\[
g_j=a+S_k b,\qquad g_k=c+S_j d,
\]

where all four coefficients are independent of both endpoint spins.
Consequently, for `j≠k`,

\[
\boxed{\quad E[S_jS_k g_jg_k]
       =E[(D_k g_j)(D_j g_k)].\quad}                  \tag{1}
\]

There is no omitted argument dependence: own-freeness is precisely what
makes each central derivative independent of the other endpoint as well.
The identity would not be justified for arbitrary recursively defined
fields that retain their own spin.

Expanding the transported square, applying (1), and then Cauchy--Schwarz
gives, for every root `i`,

\[
\begin{aligned}
E|(B[Sg])_i|^2
&\le \sum_j B_{ij}^2 E g_j^2
 +\sum_{j\ne k}|B_{ij}B_{ik}|
       \sqrt{E|D_k g_j|^2 E|D_j g_k|^2}\\
&\le \max_j E g_j^2+
       n\max_{j\ne k}E|D_k g_j|^2.                  \tag{2}
\end{aligned}
\]

Here `Σ_j B_ij²=1`; the exact off-diagonal coefficient count is at most
`n-2`, so the simpler factor `n` is valid.

## 3. Sobolev bound for a fixed smooth remainder

Let `g` be fixed smooth with polynomially bounded derivatives through
second order. In particular `g=H-P`, where `H` is fixed bounded smooth with
bounded derivatives and `P` is any fixed polynomial, is allowed.
Write `Delta=D_kU_j`. Taylor expansion of the exact two-endpoint difference,
evaluating the first derivative at the actual field, gives

\[
D_k g(U_j)=\nabla g(U_j)\cdot\Delta+R_{jk},\qquad
\sup_{j,k}\|R_{jk}\|_2=O_g(n^{-1}).                  \tag{3}
\]

Indeed a pointwise remainder bound is
`C_g |Delta|²(1+|U_j|^a+|Delta|^a)` for a fixed exponent `a`. Hölder,
the influence estimates in every fixed `L^p`, and the uniform field moments
give (3). The fact that the chosen first derivative still depends on the
endpoint spin is harmless: no independence is used in this estimate.

Another Hölder estimate gives

\[
E|\nabla g(U_j)\cdot\Delta|^2
\le (E|\nabla g(U_j)|^4)^{1/2}
     (E|\Delta|^4)^{1/2}.
\]

Multiplying by `n`, the remainder square and mixed term in (3) vanish.
Uniform Gaussian moment convergence therefore yields

\[
\limsup_n n\max_{j\ne k}E|D_k g(U_j)|^2
       \le C_U\|\nabla g(Z)\|_{L^4}^2.              \tag{4}
\]

Combining (2) and (4),

\[
\boxed{\quad
\limsup_n\max_i E|(B[Sg(U)])_i|^2
\le \|g(Z)\|_{L^2}^2+C_U\|\nabla g(Z)\|_{L^4}^2.
\quad}                                               \tag{5}
\]

The constant depends only on the fixed field family, not on the particular
polynomial approximation remainder. The vanishing error terms may depend
on that remainder; this is why the polynomial must be fixed before taking
the dimension limit.

## 4. Continuity of the normalized energy functional

For field functions `f,g`, Cauchy--Schwarz over both the spins and the root
index gives

\[
\frac1n|E f(U)^T B[Sg(U)]|
\le \left(\frac1n\sum_i E f(U_i)^2\right)^{1/2}
     \left(\frac1n\sum_i E|(B[Sg(U)])_i|^2\right)^{1/2}.
\]

Thus the normalized bilinear energy is continuous, in the ordered limit,
in the `f` argument in Gaussian `L²` and in the `g` argument in the norm

\[
\|g\|_{L^2(\gamma_d)}+\|\nabla g\|_{L^4(\gamma_d)}.
\]

In particular, a polynomial energy identity proved for all fixed
polynomials extends to fixed smooth bounded response functions by first
taking `n→∞` and then taking polynomial approximations. No polynomial
degree is permitted to grow during a dimension limit. The same inequality
also controls the error from approximating the first response function;
the approximating sequence has bounded Gaussian `L²` norm.

## 5. Self-contained polynomial density argument

The required polynomial approximation need not rely on a bare assertion
that ordinary Hermite partial sums converge in `L^4`.

For a smooth function with polynomial growth, multiply by a smooth cutoff
on an expanding ball. Gaussian tails show convergence in `W^{1,4}` (hence
also in `L²`); the cutoff-gradient term tends to zero as well. It is
therefore enough to approximate a compactly supported smooth function.

Its Fourier transform is rapidly decreasing. The Fourier integral is
absolutely integrable as an integral of functions taking values in Gaussian
`W^{1,4}`, since a plane wave `exp(i t·x)` has this norm at most `C(1+|t|)`.
Truncate the frequency integral and approximate it by finite Riemann sums.
Strong continuity in frequency, dominated by Gaussian moments of `|x|`,
shows convergence in the required norm.

Each plane wave can in turn be approximated by its Taylor polynomials.
Both the Taylor polynomials and their first derivatives are dominated by
`C_t exp(|t||x|)`, which belongs to every finite Gaussian `L^p` space.
Dominated convergence therefore proves convergence in Gaussian `W^{1,4}`.
Taking real parts gives real polynomials. This proves the required density.

If the response function is jointly even or jointly odd, replace each
polynomial by its even or odd part under `x→-x`. This projection is
contractive in the stated Gaussian norms and preserves the desired parity.

The proof applies directly to the smooth bounded response functions needed
for Boolean rounding. Extension to arbitrary nonsmooth masks would require
their own appropriate smoothing argument; it is not silently included in
the Sobolev estimate.
