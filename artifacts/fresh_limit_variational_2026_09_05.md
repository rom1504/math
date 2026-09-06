# Independent variational attack — 2026-09-05

This note was developed without reading the historical route assessments or
ledger. It concerns an auxiliary relaxation, not a proof that the Boolean
sequence converges.

## Normalizations and the first proof obligation

Let `A` be a real symmetric hollow `n` by `n` matrix with off-diagonal entries
in `{−1,+1}`, let `d=n−1`, and write

\[
H_A(x)=\tfrac12x^TAx,\qquad Q(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|.
\]

The original quantity is `M_n=min_A Q(A)` and `c_n=M_n/n^(3/2)`.
Define the absolute vector-spin relaxation

\[
V(A)=\frac12\max\{|\operatorname{Tr}(AX)|:X\succeq0,\ X_{ii}=1\}.
\]

Thus `Q(A)≤V(A)`. Hilbert-space or tensor amplification arguments that silently
replace Boolean spins by unit vectors must confront this relaxation.

**First explicit proof obligation for a vector-based convergence argument.**
One must either establish scalar-preserving realizability at arbitrary large
orders, or establish `V(A_n)−Q(A_n)=o(n^(3/2))` for an appropriate sequence of
exact Boolean minimizers. The theorem below does not establish either claim.
In fact, if the latter gap-vanishing assertion held for all exact minimizers,
the known upper bound `limsup c_n≤1/2` would force `c_n→1/2`.

## Theorem 1: exact minimizers of the relaxation are conference matrices

For every `n≥2`,

\[
V(A)\ge\frac{n\sqrt{n-1}}2.
\]

Equality holds if and only if `A²=(n−1)I`.

### Proof

Put `B=A/√d`. Then `B_ii=0`, `(B²)_ii=1`, and `Tr(B²)=n`.
The two matrices

\[
X_\pm=\frac12(I\pm B)^2
\]

are positive semidefinite and have diagonal one. Consequently, writing
`t=Tr(B³)/2`,

\[
\operatorname{Tr}(BX_+)=n+t,\qquad
\operatorname{Tr}(BX_-)=-n+t.
\]

The maximum of the absolute values is `n+|t|≥n`, proving the lower bound.
This is also the Gram construction
`u_i^±=(e_i±B_i)/√2`.

Suppose equality holds. The preceding formulas force `Tr(B³)=0`, and `X_+`
is optimal for the SDP `max Tr(BX)` with `X` a correlation matrix. Its dual is

\[
\min\left\{\sum_i z_i:\operatorname{diag}(z)-B\succeq0\right\}.
\]

Strict primal and dual feasibility give an attained dual optimum `Z=diag(z)`
with `Tr Z=n`. Complementary slackness gives

\[
(Z-B)(I+B)=0.
\]

The diagonal of this identity is `z_i−(B²)_ii=z_i−1`, so `Z=I`.
The full identity then gives `I−B²=0`. Conversely, if `B²=I`, then
`−I≤B≤I`, hence `|Tr(BX)|≤Tr X=n` for every correlation matrix `X`.
This proves equality and the classification. □

The converse can equivalently be seen from the spectral bound
`Q(A)≤V(A)≤n ||A||op/2`.

## Theorem 2: quantitative spectral stability

Suppose, for some `ε≥0`,

\[
V(A)\le(1+\varepsilon)\frac{n\sqrt{n-1}}2.
\]

Then

\[
\frac1n\sum_{j=1}^n\bigl(|\lambda_j(B)|-1\bigr)^2
\le 2(1+\sqrt2)\varepsilon.
\]

Equivalently, there is a real symmetric orthogonal matrix `U` for which

\[
\frac{\|A-\sqrt{n-1}\,U\|_F^2}{n(n-1)}
\le 2(1+\sqrt2)\varepsilon.
\]

This does not assert entrywise closeness to a sign conference matrix, nor
existence of a conference matrix of the same order.

### Proof

Let `Z=diag(z)` and `W=diag(w)` be optimal duals for the positive and negative
SDPs, so

\[
Z-B\succeq0,\quad W+B\succeq0,\quad
\operatorname{Tr}Z,\operatorname{Tr}W\le(1+\varepsilon)n.
\]

All `z_i,w_i` are strictly positive: each row of `B` has a nonzero
off-diagonal entry, while a positive semidefinite matrix with a zero diagonal
entry has a zero corresponding row. Let

\[
\delta_+=\operatorname{Tr}Z-n-\tfrac12\operatorname{Tr}B^3,\qquad
\delta_-=\operatorname{Tr}W-n+\tfrac12\operatorname{Tr}B^3.
\]

Feasibility of `X_±` gives `δ_±≥0`, and
`δ_++δ_-≤2εn`. Moreover,

\[
\|(Z-B)^{1/2}(I+B)\|_F^2=2\delta_+.
\]

The diagonal of `(Z−B)(I+B)` is `z_i−1`. Cauchy–Schwarz, applied to the
`i`th column and `(Z−B)^(1/2)e_i`, therefore yields

\[
\sum_i\frac{(z_i-1)^2}{z_i}\le2\delta_+.
\]

The analogous minus statement holds for `w`. Set `S=(Z+W)/2`, with diagonal
`s_i`. Convexity of `f(s)=(s−1)²/s` now gives

\[
\sum_i\frac{(s_i-1)^2}{s_i}\le2\varepsilon n,
\qquad \sum_i s_i\le(1+\varepsilon)n. \tag{1}
\]

Write `B=B_+−B_-` for its positive/negative spectral parts. Multiplying the
two dual inequalities by `B_+` and `B_-`, respectively, and taking traces,
gives

\[
n=\operatorname{Tr}B^2
\le\operatorname{Tr}(ZB_++WB_-)
=\operatorname{Tr}(S|B|).
\]

The last equality holds because `B_ii=0` implies
`(B_+)_ii=(B_-)_ii=|B|_ii/2`. Put `h_i=|B|_ii`.
Cauchy–Schwarz in the spectral measure of `e_i` gives
`0≤h_i≤√((B²)_ii)=1`. Therefore

\[
0\le\sum_i s_i(1-h_i)
\le\sum_i s_i-n\le\varepsilon n. \tag{2}
\]

Using (1), (2), and `(1−h_i)²≤1−h_i`,

\[
\begin{aligned}
n-\|B\|_*
&=\sum_i(1-h_i)\\
&=\sum_i s_i(1-h_i)+\sum_i(1-s_i)(1-h_i)\\
&\le\varepsilon n+
\sqrt{\left(\sum_i\frac{(1-s_i)^2}{s_i}\right)
\left(\sum_i s_i(1-h_i)^2\right)}\\
&\le(1+\sqrt2)\varepsilon n.
\end{aligned}
\]

Finally, `Tr B²=n`, so

\[
\sum_j(|\lambda_j(B)|-1)^2=2n-2\|B\|_*
\le2(1+\sqrt2)\varepsilon n.
\]

Choosing `U` to have the same eigenvectors as `B` and eigenvalues
`sign(λ_j(B))` (either sign at zero) proves the equivalent statement. □

## Status and next falsifier

The two proofs above were frozen before archive comparison. The exact
classification and stability need independent audit. Next, test the scalar
gap in conference families. Even a fixed asymptotic gap for a single family
would refute a universal assertion that conference/vector extremality implies
Boolean/vector agreement; it would not by itself refute a gap-vanishing
assertion restricted to exact Boolean minimizers.

For any block argument, the actual identity remains

\[
Q\!\begin{pmatrix}A&B\\B^T&D\end{pmatrix}
=\max_{x,y}\bigl(|H_A(x)+H_D(y)|+|x^TBy|\bigr).
\]

No signed-cancellation shortcut has been used here.

## Sharper polar-Gram result (root contribution, independently checked)

After the preceding proofs were frozen, the root agent supplied the following
shorter argument and stronger stability bound. The derivation below has been
independently checked by this agent.

### Theorem 3: exact vector/nuclear tradeoff

\[
\boxed{V(A)\,\|A\|_*\ge\frac{n^2(n-1)}2.}
\]

Let `h_i=|A|_ii`, `D=diag(h)`, and `S=||A||_*=Σ_i h_i`.
The `h_i` are strictly positive, and

\[
X_\pm=D^{-1/2}(|A|\pm A)D^{-1/2}
\]

are correlation matrices. Their objective difference gives

\[
V(A)\ge\frac12\operatorname{Tr}(AD^{-1/2}AD^{-1/2})
=\frac12\sum_{i\ne j}\frac1{\sqrt{h_ih_j}}.
\]

The function `(u,v)↦(uv)^(−1/2)` is jointly convex on the positive quadrant.
Jensen over ordered distinct pairs gives

\[
\frac1{n(n-1)}\sum_{i\ne j}\frac1{\sqrt{h_ih_j}}
\ge\frac nS.
\]

This proves the tradeoff. Since `S≤n√(n−1)` by Cauchy–Schwarz on singular
values, it also gives Theorem 1 immediately. Equality in the universal bound
forces equality in that Cauchy–Schwarz inequality, hence `A²=(n−1)I`.

Under the hypothesis of Theorem 2 it gives `||B||_*≥n/(1+ε)`, and therefore
the sharper bound

\[
\boxed{
\frac1n\sum_j(|\lambda_j(B)|-1)^2
\le\frac{2\varepsilon}{1+\varepsilon}.
}
\]

The longer dual proof above is retained as the independently frozen initial
derivation, but this is the preferred statement and proof.

### Theorem 4: Boolean/nuclear tradeoff by paired Gaussian rounding

\[
\boxed{
Q(A)\ge\frac{n(n-1)}\pi
\arcsin\!\left(\frac n{\|A\|_*}\right)
\ge\frac{n^2(n-1)}{\pi\|A\|_*}.
}
\]

Here `n≤||A||_*≤n√(n−1)`. The lower bound follows, for example, from
`||A||op≤n−1` and `||A||_*≥||A||_F²/||A||op`.

Use the polar correlation matrices `X_±` from Theorem 3. For an edge `ij`,
put

\[
p_{ij}=\frac{|A|_{ij}}{\sqrt{h_ih_j}},\qquad
q_{ij}=\frac{a_{ij}}{\sqrt{h_ih_j}}.
\]

Let `Y_±` be signs of centered Gaussian vectors with covariance `X_±`.
The Gaussian arcsine identity gives

\[
\mathbb E H_A(Y_\pm)=\frac2\pi
\sum_{i<j}a_{ij}\arcsin(p_{ij}\pm q_{ij}).
\]

For `|p|+t≤1`, `t≥0`,

\[
\arcsin(p+t)-\arcsin(p-t)\ge2\arcsin t. \tag{3}
\]

Indeed, the left side is even in `p` and nondecreasing for `p≥0`, since
`arcsin′(u)=(1−u²)^(−1/2)` is even and increases with `|u|`. The
correlation constraints `p±q∈[−1,1]` ensure `|p|+|q|≤1`.
Since `sign(q_ij)=a_ij`, (3) and the fact that each expected energy lies in
`[−Q(A),Q(A)]` imply

\[
Q(A)\ge\frac2\pi\sum_{i<j}
\arcsin\!\left(\frac1{\sqrt{h_ih_j}}\right).
\]

Jensen first for the convex function `arcsin` on `[0,1]`, and then for
`(uv)^(−1/2)`, proves the boxed inequality. □

### Interpretation and limits

The theorem is a genuine Boolean lower bound but does not improve the
currently known universal constant `0.336493364431`: taking only
`||A||_*≤n√(n−1)` gives asymptotic constant `1/π`.
It does quantitatively exclude large nuclear deficiency in a low-energy
signing. If `Q(A_n)/n^(3/2)→c`, then

\[
\liminf\frac{\|A_n\|_*}{n\sqrt{n-1}}\ge\frac1{\pi c}.
\]

At `c=0.336493364431` the right side is about `0.946`. This suggests a
specific additional obligation: extend the stronger exact-conference
rounding lower bounds robustly in the normalized spectral defect
`2−2||A||_*/(n√(n−1))`, with a useful explicit error. Neither such a robust
bound nor a resulting improvement is claimed here.

## Numerical smoke test and exact finite conference caps

`computations/fresh_limit_variational_audit.py` checked all displayed
inequalities on 33 seeded random signings (three at every order 2 through 12)
and the four conference matrices below. The Boolean caps were obtained by
exhausting all projective Boolean assignments. The SDP values were only
numerical diagnostics, computed with CLARABEL; their exact values for
conference matrices follow from Theorem 1.

| Conference field | Order | Exact Boolean cap | Exact vector cap | Boolean/vector ratio |
|---|---:|---:|---:|---:|
| GF(5) | 6 | 5 | `3√5` | 0.745356 |
| GF(9) | 10 | 15 | 15 | 1 |
| GF(13) | 14 | 21 | `7√13` | 0.832050 |
| GF(17) | 18 | 33 | `9√17` | 0.889297 |

The order-6 and order-14 examples are also exact Boolean minimizers, according
to the exact values supplied in the task. Thus even simultaneous Boolean and
vector extremality does not imply a zero finite scalar/vector gap. These
finite computations do not decide whether the gap is subleading.

## Additional candidate: a sign-specific direct nuclear bound

The root proposed testing the following much stronger assertion:

\[
Q(A)\stackrel{?}{\ge}\frac12\|A\|_*-O(n),
\]

uniformly over hollow symmetric sign matrices. This is not asserted as a
theorem. Together with Theorem 4 it would imply
`liminf c_n≥1/√(2π)≈0.398942`, because with
`s=||A||_*/n^(3/2)` the two bounds would give
`c≥max(s/2,1/(πs))` asymptotically.

The exact finite conference deficits `||A||_*/2−Q(A)` above are approximately
`1.7082,0,4.2389,4.1080`; all are compatible with an `O(n)` remainder.
No scalable sign counterexample has been obtained.

An important amplification obstruction is elementary. For every symmetric
Hadamard matrix `H` of order `m`, the vector `vec(H)` is Boolean and

\[
(H\otimes H)\operatorname{vec}(H)=m\operatorname{vec}(H).
\]

Thus the tensor square saturates the spectral bound regardless of a scalar
gap in the seed. More generally, the partial Boolean witness `vec(A)` for
`A⊗A` has quadratic half-energy `Tr(A⁴)/2`, which is at least
`n(n−1)²/2`. Its zeros are precisely on the `n` diagonal coordinates of
`vec(A)`, and can be randomly filled without decreasing its expected
half-energy. The product matrix itself still has zero edge entries and is
not an admissible signing; no claim about an arbitrary sign completion is
made. This explains why naive tensor-square tests favor a coefficient
`1/2` instead of preserving a low seed ratio.

As a literature sanity check, the maximum-excess problem distinguishes
maximizing over all Hadamard matrices from maximizing within a fixed
equivalence class. [Hirasaka, Momihara, and Suda's primary paper](https://www.numdam.org/item/ALCO_2018__1_5_697_0.pdf)
constructs particular regular or biregular families; it does not give a
uniform same-spin near-regularization theorem for every symmetric sign
matrix. Therefore these constructions cannot be substituted for the
displayed direct nuclear assertion.

## Attempted robust extension: a simpler best-response second step

The archived exact-conference two-step calculation was inspected only after
the independent proofs above. A direct `B↔sign(B)` coupling gives an `L²`
error proportional to `sqrt(1−||B||_*/n)`. Passing through a sharp final
threshold requires either anti-concentration or a fixed dither; the naive
dither-Lipschitz bound worsens the constant. No useful linear-defect
cancellation was established.

For completeness, there is a fully deterministic random-update inequality.
For any hollow symmetric `B`, a Boolean vector `x`, and
`y=sign(Bx)`, independently replace each coordinate of `x` by the
corresponding coordinate of `y` with probability `u`. Conditional independence
and the absence of a diagonal give

\[
\mathbb E H_B(x')=(1-u)^2H_B(x)
+u(1-u)\|Bx\|_1+u^2H_B(y).
\]

Consequently, for any distribution of initial `x`,

\[
(1+u^2)Q(B)\ge (1-u)^2\mathbb E H_B(x)
+u(1-u)\mathbb E\|Bx\|_1. \tag{4}
\]

Applied formally to the exact-conference one-probe population, set
`a=2φ(t)`, `b=2Φ(t)−1`, where `a=tb` at the one-probe optimum
`t≈0.876900985553`. Then the initial half-energy is `ab≈0.336493364432`.
The full conference local field has the law
`aS+sqrt(1−a²)G`, whose absolute mean is about `0.805416113371`.
Substitution in (4), followed by optimizing `u`, gives only about
`0.341900688967`, at `u≈0.0816623`.

This is a population calculation conditional on the conference field law,
not an independently proved universal theorem. It exposes the obstruction:
if one retains just the conditional linear part `aS+bG`, its absolute mean
is exactly `2ab`, and (4) gives no improvement at all. All gain is supplied
by the nonlinear residual field; the simpler update does not eliminate the
missing residual-universality proof.

## Checkpoint conclusion

Proved and independently checked: vector/conference classification, the
sharp inverse-nuclear vector tradeoff and stability, and the exact paired
Gaussian Boolean/nuclear inequality. Computations corroborate the formulas
and separate finite scalar and vector extrema. No convergence result,
strict liminf/limsup separation, improved universal lower constant, or
scalable counterexample to the proposed direct nuclear bound is claimed.

## Second assignment: parity covariance and the true regularization

Historical-stage note: later in the same campaign Schmidt's primary theorem
settled the H2 saturation question below. See
`fresh_schmidt_odd_walsh_regularization_2026_09_05.md`. The finite plateaus
and restricted-certificate obstructions remain valid; they do not imply a
positive uniform global defect. The final original-problem frontier is in
`ACTIVE_STATE.md`, not the earlier numerical targets in this chronological report.

This section concerns full symmetric real seeds, including their diagonals.
Use `q(B)=max_x |x^TBx|/2` here, and let `H_s` range over the fixed regular
Hadamard semigroup generated by the order-4 and order-144 matrices in the
algebra agent's note. Define

\[
R(B)=\sup_s\frac{q(H_s\otimes B)}{s^{3/2}},\qquad
T(B)=\frac12\min_{D\succeq\pm B}\max_xx^TDx.
\]

For full sign seeds the diagonal can be removed after amplification at a
subleading cost, so a certified `R(B)<m^(3/2)/2` would be materially useful.
The known `T(B)≥m^(3/2)/2` for every full symmetric sign seed means that such
a certificate must distinguish `R` from `T`.

### Exact reformulation 1: regularized bilinear norm

The algebra agent's one-factor conversion and tensor stability imply the
stronger identity

\[
\boxed{2R(B)=\sup_s
\frac{\|H_s\otimes B\|_{\infty\to1}}{s^{3/2}}.} \tag{5}
\]

Indeed `q(C)≤||C||_(∞→1)/2` gives one inequality. Conversely, their explicit
order-4 conversion gives `R(C)≥||C||_(∞→1)/2` for every symmetric `C`.
Apply this to `C=H_s⊗B` and use `R(H_s⊗B)=s^(3/2)R(B)`.

### Exact reformulation 2: `2T` is a Hilbert-factorization ideal norm

Let

\[
\Gamma_2(B:\ell_\infty^m\to\ell_1^m)
=\inf_{B=P^TQ}\|P\|_{\infty\to2}\|Q\|_{\infty\to2},
\]

allowing arbitrary finite Hilbert-space dimension for the two factors.
This convention is **not** the frequently used `γ₂` norm for maps from
`ℓ₁` to `ℓ∞`, and is not its ordinary vector-SDP dual.
For symmetric `B`,

\[
\boxed{2T(B)=\Gamma_2(B:\ell_\infty^m\to\ell_1^m).} \tag{6}
\]

To prove the nontrivial direction, balance a factorization `B=P^TQ` by a
scalar rescaling so `||P||_(∞→2)=||Q||_(∞→2)=sqrt(c)`.
Symmetry gives `B=(P^TQ+Q^TP)/2`. Therefore

\[
D=\frac12(P^TP+Q^TQ),\qquad
D\pm B=\frac12(P\pm Q)^T(P\pm Q)\succeq0,
\]

and `max_x x^TDx≤c`, proving `2T≤Γ₂`. In the reverse direction,
`D⪰±B` first implies `D⪰0` and `ker D⊆ker B`. On the support of `D`,
set `C=D^(-1/2)BD^(-1/2)`; then `−I⪯C⪯I`. Extending `C` by zero on
the kernel gives `B=D^(1/2) C D^(1/2)` with `C` a self-adjoint contraction.
The factors `P=D^(1/2)` and `Q=CD^(1/2)` have norm product at most
`||D^(1/2)||_(∞→2)²=max_xx^TDx`.

Consequently, `R=T` asks whether the fixed regular-Hadamard tensor tests in
(5) compute this Hilbert-factorization ideal norm. It is not the already
falsified identification with the same-spin vector SDP.

### The parity-covariance formulation and what convexification does give

For the order-4 generator alone,

\[
U_4=H_4/2=2E-I,
\]

where `E` is averaging on a four-point probability space. Its `t`th tensor
power acts as `+1` on even ANOVA levels and `−1` on odd levels.
For Boolean vector-valued `f=(f_1,...,f_m)`, put

\[
G_{ij}=\mathbb E f_if_j,\qquad
K_{ij}=\mathbb E f_i(U_4^{\otimes t}f_j).
\]

Then `G∈CUT_m` and `−G⪯K⪯G`. The latter follows by splitting each function
into its even and odd components. The `T` dual relaxes to all pairs satisfying
these conditions; equality with `R` requires a realizability theorem, not
just the displayed necessary conditions.

The closure of the actual parity-covariance pairs is convex. To see this,
use a selector depending only on differences of fresh pairs of four-point
coordinates. Every function of such differences is in the `+1` eigenspace:
each nonconstant character occupies both coordinates of its pair. A selector
of probability approaching any desired `λ` can therefore combine two
independent constructions without changing their parity action. Disjointness
of selector events makes both `G` and `K` the corresponding convex mixtures.

A fresh balanced `+1` eigenfunction makes all individual means zero without
changing `(G,K)`. A balanced `−1` eigenfunction changes `(G,K)` to `(G,−K)`.
For rank-one `G=ss^T`, these constructions realize every
`K=κss^T`, `|κ|≤1`; this is the entire relaxed slice in that case. Any
separation must therefore exploit coherence between distinct cut states.

These claims are for the order-4 parity family. A bound for this smaller
family need not hold after allowing the additional order-144 generator.

### A rigorous obstruction to a broad class of nonquadratic certificates

Suppose a continuous even function `Φ:R^m→R` satisfies `Φ(0)=0`
and is a pointwise separable contraction potential for the order-4
transform:

\[
\sum_{a=1}^4\Phi\bigl((U_4z)_a\bigr)
\le\sum_{a=1}^4\Phi(z_a)
\qquad\text{for every }z_1,z_2,z_3,z_4\in\mathbb R^m. \tag{7}
\]

Then **`Φ` is a quadratic form**. No homogeneity assumption is necessary.
Since `U₄²=I`, applying (7) twice makes it equality. Input `(u,0,0,0)`
and evenness first give `Φ(u)=4Φ(u/2)`. For `z=(u,v,0,0)`, its transform is

\[
\tfrac12(-u+v,\ u-v,\ u+v,\ u+v).
\]

Evenness and this derived dyadic scaling now give the parallelogram identity

\[
\Phi(u+v)+\Phi(u-v)=2\Phi(u)+2\Phi(v).
\]

Polarization and continuity imply `Φ(x)=x^TDx` for a symmetric matrix `D`.
If `Φ≥0`, then `D⪰0`.

Thus a tensor-stable separable `L²` potential, valid on arbitrary real vector
inputs, cannot furnish a genuinely nonquadratic improvement over the
quadratic-majorant paradigm. A successful nonquadratic certificate would
have to retain Boolean reachability, additional joint history, or another
constraint not captured by such a pointwise contraction.

### Minimal full-sign test case

The seed `B=[[1,1],[1,−1]]` has `q(B)=1` and `T(B)=√2`.
Exact enumeration gives `q(H₄⊗B)=10`, hence `R(B)≥1.25`.
Saturation of `T` in the parity formulation requires `G=I₂` and, up to the
objective orientation,

\[
Uf\approx(f+g)/\sqrt2,\qquad
Ug\approx(f-g)/\sqrt2.
\]

This is an approximate Boolean-invariant two-dimensional rotation problem.
Finite-order integrality prevents exact equality, but its elementary gap
vanishes with the outer order and supplies no asymptotic separating
certificate. Neither saturation nor a uniform gap is claimed.

### Exact outer-order 16 audit for the minimal full seed

The integer-only script `computations/fresh_walsh32_exact_certificate.py`
proves

\[
\|W_{32}\|_{\infty\to1}=160,
\qquad q(H_4^{\otimes2}\otimes H_2)=80.
\]

Here `W32` is the Sylvester Walsh matrix. Enumerating the 32768 Boolean
inputs to `W16` up to global complementation gives eight sorted absolute
spectrum types. The Walsh recursion writes a length-32 transform as
`(a+b,a−b)`, so its 1-norm is `2 sum_i max(|a_i|,|b_i|)`.
The bound obtained by assigning any `k` coordinates to the `k` largest
entries of the first type and the remaining coordinates to the second
type is at most 160, except for the pair of types

\[
(0^6,4^8,8^2),\qquad (2^{10},6^6).
\]

These two types have respectively 840 and 448 distinct absolute profiles.
All 376320 profile comparisons have value at most 160; their exact value
histogram is `{136:6720,144:120960,152:208320,160:40320}`.
This is an exhaustive finite integer certificate, with no floating-point
optimizer. Independent signed equivalence of `H4` to `W4` transfers the
bilinear upper bound to the full seed lift. The script supplies a Boolean
vector of quadratic energy −160, proving equality for `q`.

Thus the normalized witness is still `80/16^(3/2)=5/4`, the same value as
outer order four. **This is only a finite-level exact result, not an upper
bound for `R(H2)`.** In particular, algebra's outer-order 64 witness already
has normalized value `676/512=1.3203125`.

### A concrete non-Walsh invariant-partition test

The following is an exact reduction of a proposed realization step, **not
an existence claim**. Let `F` be the order-12 Hadamard used in the square
construction of `H144`. Index the latter by ordered pairs `(i,j)` and put
`D_(i,j)=F_ij`. Its action on a real 12-by-12 array `X` is

\[
U_{144}X=F\circ\left[\frac1{12}F(F\circ X)^TF\right],
\]

where `∘` denotes entrywise multiplication. Thus, writing
`Z=F∘X` and `M(X)=Z^TF/12`, we have `M(U144 X)=M(X)^T`.

An equitable partition of the 144 coordinates into three cells of size 48
with quotient `2E3−I3` is equivalent to three disjoint 0/1 masks `P_c`
summing to the all-one matrix and satisfying

\[
U_{144}P_c=\frac23\mathbf1-P_c,
\quad\text{equivalently}\quad
Z_c^TF+F^TZ_c=8I,
\qquad Z_c=F\circ P_c. \tag{8}
\]

The diagonal equations force exactly four selected entries in every
column of each mask. In the original Hadamard matrix, the required block
row sums are −4 into one's own cell and +8 into either other cell.

A decomposition of this specific `F` into disjoint matrices `Z_1,Z_2,Z_3`
forming an orthogonal design of type `(4,4,4)` would suffice: the identities
`Z_c^TZ_c=4I` and `Z_c^TZ_d+Z_d^TZ_c=0` imply (8).
Existence of such an orthogonal design in the relevant signed class, or of
the weaker masks (8), has not been established here. Algebra is testing
the finite partition separately.

The motivation is to obtain a genuine additional invariant probability
algebra from the order-144 generator, rather than assuming every conclusion
about Walsh parity survives that generator. Even a successful partition
would not by itself prove `R=T`: a coherent composition/realization
argument would still be needed.

## Adversarial audit of the rooted-response lower-bound candidate

This section audits the separate proposed proof in
`artifacts/fresh_limit_rooted_response_2026_09_05.md` (root and algebra
contributions), following the parent's priority change at approximately
19:55 UTC. The earlier `R/T` investigation is paused. The conclusions below
are an independent proof audit, not a claim to have settled convergence.

### Normalization and spectral bootstrap

Use `m=n−1`, `B=A/sqrt(m)`, `Q=B²`, and independent Rademacher `S`.
For a hypothetical low-cap sequence `q(A)≤C n^(3/2)`, we indeed have
`||Q||op=O(sqrt(n))`, `Qii=1`, and `Tr Q²=O(n^(3/2))`.
Besides the eigenvector proof in the main candidate, there is a short
alternative verification: every row of `A` belongs to the cube, hence
`||A²||_(∞→∞)≤β(A)`, where `β(A)=max_{x,y}x^TAy` over sign vectors.
Thus `||A||op²≤β(A)≤4q(A)` by symmetry and polarization. No bounded
operator norm of `B` is being assumed.

### Endpoint covariance and the critical tail estimate

For a fixed smooth even `g`, set `w_j=S_jg(G_j)` with `G=BS`.
For distinct `j,k`, averaging `S_j,S_k` gives exactly

\[
\mathbb E w_jw_k=\tfrac14\mathbb E
 [g(U+B_{jk})-g(U-B_{jk})]
 [g(V+B_{jk})-g(V-B_{jk})].
\]

The centered differences are `2Bjk g'(U)+O(Bjk³)` in every required
fixed moment. Hence the leading term is `m^(-1)E g'(U)g'(V)`.
A replacement of the other input spins by Gaussians costs `O(m^(-1/2))`
*inside this expectation*, uniformly even when the two fields are nearly
perfectly correlated. Indeed one applies the ordinary third-order Taylor
replacement to `g'(u)g'(v)`, not a Berry--Esseen bound requiring an inverse
covariance. Each two-dimensional increment has size `O(m^(-1/2))`;
polynomial-growth derivative remainders have bounded moments from row
subgaussianity. The two common-field variances are `1−1/m` and their
covariance is exactly `Qjk`. Adding independent Gaussian increments of
variance `1/m` changes the smooth expectation by `O(1/m)`.

Since `g'` is odd, its Gaussian covariance kernel has only positive
coefficients on odd Schur powers:

\[
K_g(Q)=\sum_{r\text{ odd}}c_r Q^{\circ r},\qquad c_r\ge0,
\quad \sum c_r=\mathbb E g'(Z)^2.
\]

The PSD correlation-matrix Schur multiplier is a contraction on operator
norm. Therefore `||Kg(Q)||op≤E g'(Z)² ||Q||op`. Sandwiching by one flat
unit row of `B` costs only `O(||Q||op/m)=o(1)`; the entrywise remainder
costs `O(m^(-3/2)) ||B_i||_1²=O(m^(-1/2))`. This verifies

\[
\mathbb E(Bw)_i^2=\mathbb E g(Z)^2+o(1)
\]

uniformly in `i`. This estimate is valid for a smooth bounded function
minus a *fixed* polynomial. Its constants can grow with polynomial degree;
the order of limits must remain `n→∞`, then degree `→∞`.

### Explicit finite-degree comparison

For a row of `m` coefficients `a_j=±m^(-1/2)`, let `P_r` be the ordered
distinct-index homogeneous sum of degree `r`, with `P_0=1`, `P_1=G`.
Multiplying by `G` and separating the new index from the `r` existing
indices gives the exact recurrence

\[
P_{r+1}=G P_r-r\bigl(1-(r-1)/m\bigr)P_{r-1}.
\]

Comparison with the Hermite recurrence shows that for each fixed `r`,
`H_r=P_r+sum_{s<r,s≡r mod2} O_r(1/m)P_s`. For example,

\[
H_4=P_4-\frac8m P_2-\frac2m.
\]

Multiplication by the root spin `S_j` creates no repetition at that root,
because `Bjj=0`. The even tail-transport estimate above, used inductively,
bounds all transported lower-degree row polynomials. Thus their
`O_r(1/m)` coefficients really give an `o(1)` error after transport.

For the resulting degree-`r+1` multilinear polynomial, a fixed support has
at most `r+1` root choices. Its coefficient is bounded by
`C_r m^(-(r+1)/2)`. Counting supports containing one specified input gives
influence `O_r(1/n)`. Fixed-degree hypercontractivity and the ordinary
one-input Taylor replacement therefore give joint invariance with the
local linear field, with error `O_r(n^(-1/2))` for smooth tests.

The Gaussian counterpart is a pure Wick chaos: `Z_j` is independent of
`(BZ)_j`. Terms omitted by multilinearization have at most `r` distinct
indices, coefficients `O_r(m^(-(r+1)/2))`, and Hermite factorial weights
bounded in terms of the fixed degree. Their squared norm is `O_r(1/n)`.
For `r=2`, the Rademacher rooted polynomial is already exactly cubic and
multilinear, which supplies a particularly transparent independent check.

### Directional decoupling

Gaussian integration by parts verifies the candidate's exact covariance

\[
\mathbb E[Z_jH_r(G_j)Z_kH_r(G_k)]
=r!1_{j=k}+\frac{r r!}{m}1_{j\ne k}Q_{jk}^{r-1}.
\]

The derivative in the unit direction `B_i` is the sum of
`m^(-1)sum_{j≠i}H_r(G_j)` and
`r sum_j Bij Qij Z_j H_(r−1)(G_j)`.
The first squared norm is at most `r! Tr Q²/m²=o(1)`.
For the second, when `r=2` the covariance matrix is
`I+(J−I)/m`, hence has bounded operator norm; when `r≥4`, the same Schur
bound applies. Its squared coefficient norm is at most `(Q²)ii/m=o(1)`.
These checks validate the conditional Gaussian Poincaré step: each fixed
rooted chaos differs in `L²` by `o(1)` from a function of the Gaussian
coordinates orthogonal to `G_i`, and is therefore asymptotically independent
of `G_i`. No residual Gaussian limit is needed.

### Moment passage and strictness

The fixed-polynomial approximants have uniform fourth moments by
hypercontractivity. The transported `L²` approximation error tends to zero
as the approximation degree grows, after taking `n→∞`. This implies
uniform integrability of the residual squares, not merely tightness:
split `R=R_D+(R−R_D)`, use small `L²` error for the second term, and the
fourth-moment bound for the fixed approximant on exceptional events.
Consequently the second moments and the covariance with the cubic witness
pass to subsequential limits. Even without this full square-UI observation,
the covariance `RF` is UI from bounded `ER²` and bounded `EF⁴` by splitting
at `|F|>M` and then at `|R|>L/M`.

The candidate's constants are safely conservative. At small smoothing,
`|E RF|≥0.4`, `ER²≤1`, `EF²=2`, and `EF⁴≤2916`. With
`E={|R|>0.1}`,

\[
0.4\le0.1\sqrt2+2916^{1/4}\mathbb P(E)^{1/4},
\]

so the asserted `P(E)≥4e−7` follows. For
`a∈[0.5,0.6]`, `b∈[0.5,0.7]`, and `|r|≤0.1`,

\[
j''(r)\ge\frac2{0.7}\phi(1.4)>0.4,
\qquad j(r)=\mathbb E\max(a,|bZ+r|).
\]

Thus `j(0.1)−j(0)>0.002`, and the response gain `4e−10` in the candidate
has a factor-two safety margin. Since `j` is even and strictly minimized
at zero, this strictness does not actually need the residual mean-zero
condition once independence from `Z` has been established.

### Own-spin removal and final normalization check

Removing `S_i` from other row fields gives
`(Bu)_i=aS_i+W_i+o_(L²)(1)`, with `W_i` independent of `S_i`.
Concentration of the empirical even derivative follows from
`Tr Q²/n²=o(1)`. The analogous change in `(Bv)_i` has leading term
`S_i m^(-1)sum_j S_j h'(G_j^(i))`; endpoint extraction bounds its
off-diagonal covariances by `O(1/m)`, so its `L²` norm vanishes.
The replacement fields are jointly independent of `S_i`; averaging that
spin and using convexity gives `E max(|Bu|,|Bv|)≥E max(a,|Bv|)+o(1)`.
No tightness or uniform moment bound on `W_i` is needed at this step.

Finally, with `e/n→c*` and `ell/n≥2c*+δ`, the exact partial-update formula
gives the limiting improvement

\[
\frac{(1-p)^2c_*+p(1-p)(2c_*+\delta)}{1+p^2}-c_*
=\frac{p[\delta-(2c_*+\delta)p]}{1+p^2}>0
\]

for `0<p<δ/(2c*+δ)`. With the conservative `δ=4e−10`, choosing
`p=1e−10` gives a gain exceeding `3e−20` before any desired further
rounding of the stated constant. This arithmetic is conditional on the
candidate proof and on its baseline estimate; it is not a convergence
claim. The conversion from `q(B)/n` to `q(A)/n^(3/2)` multiplies by
`sqrt((n−1)/n)→1` and loses no asymptotic constant.

**Audit status:** no fatal gap found in Sections 3--7 after the explicit
checks above. The candidate's finite-degree and limiting mechanisms are
substantive and appear valid; they should be retained explicitly in any
final proof rather than replaced by a general AMP/universality assertion.

### Independent derivation of the baseline energy

The baseline need not be imported from an earlier one-probe theorem.
Conditional independence of the dithered coordinates and hollowness give

\[
e=\mathbb E u^TBv
  =\sum_{i\ne j}B_{ij}\mathbb E[f(G_i)S_jh(G_j)].
\]

Write `G_i=U+B_ij S_j`. The pair `(U,G_j)` is independent of `S_j`,
so averaging that spin gives

\[
\mathbb E[f(G_i)S_jh(G_j)]
=B_{ij}\mathbb E[f'(U)h(G_j)]+O_\tau(m^{-3/2}).
\]

The covariance of `(U,G_j)` is exactly `Qij`; their variances are
`1−1/m` and `1`. The same smooth two-dimensional replacement as above
therefore yields

\[
e=\frac1m\sum_{i\ne j} L(Q_{ij})+O_\tau(\sqrt n),
\quad L(q)=\mathbb E[f'(Z_1)h(Z_2)],\quad
\operatorname{Corr}(Z_1,Z_2)=q.
\]

Both functions in this last expectation are even. Their centered Hermite
expansions start at degree two, and Cauchy--Schwarz gives

\[
|L(q)-ab|\le q^2
 \sqrt{\operatorname{Var}(f'(Z))\operatorname{Var}(h(Z))}.
\]

Consequently the total centered term is at most
`C_tau Tr Q²/m=O_tau(sqrt(n))`. Since there are exactly `nm`
off-diagonal terms,

\[
\boxed{e=ab\,n+O_\tau(\sqrt n).}
\]

This closes the baseline obligation under the same low-cap hypothesis
already used by the response-field argument. All smoothing constants are
fixed before taking `n→∞`.

## Independent audit of the stronger rooted Gaussian limit

The parent and algebra subsequently strengthened the response argument by
proving that the entire *rooted marginal* converges to a Gaussian, using
explicit chaos contractions. I independently checked the proof in
`artifacts/fresh_limit_rooted_gaussian_2026_09_05.md`, including the two
cross-root cases and completeness under symmetrization.

For fixed even `r≥2`, fix a root coordinate `i`, write `a_j=B_ij`,
`b_j=B_j`, and use the unsymmetrized kernel

\[
T=\sum_j a_j e_j\otimes b_j^{\otimes r}.
\]

Let `1≤ell≤r` be a contraction order and `q=r−ell+1≥1`.
Expanding contractions of `Sym(T)` reduces to contractions of two permuted
copies of `T`. Permutations of uncontracted coordinates do not affect the
norm. The four possibilities are exhaustive:

1. The roots are paired to each other. Their inner product identifies
   their summation indices. The squared norm is
   `sum_jk a_j²a_k² Q_jk^(2q)≤Tr Q²/m²`.
2. Neither root is contracted. The two surviving basis vectors identify
   the corresponding indices when squaring. The squared norm is
   `sum_jk a_j²a_k² Q_jk^(2ell)≤Tr Q²/m²`.
3. Exactly one root meets an opposite branch. For fixed surviving root `k`,
   put `(v_k)_j=a_j B_kj Q_jk^(ell−1)`. The squared norm is exactly
   `sum_k a_k² v_k^T P v_k`, where `P=Q^{∘q}`. Since
   `||P||op≤||Q||op` and `||v_k||²≤1/m`, this is at most
   `||Q||op/m`.
4. Both roots meet opposite branches (necessarily `ell≥2`). The remaining
   coefficient matrix is `C_jk=a_j a_k B_jk² Q_jk^(ell−2)`, with zero
   diagonal. The squared norm is exactly `Tr(C P C^T P)`, and hence is at
   most `||P||op² ||C||F²≤||Q||op²/m²`.

The last step uses `sum_j a_j²=1`; importantly the exponent `q` never
vanishes for a nontrivial contraction. Thus every displayed bound tends
to zero. The number of permutations depends only on the fixed chaos
degree, so the triangle inequality proves the same conclusion for the
fully symmetrized kernel. This check finds no missing mixed root placement.

I directly read the required primary statements:
[Nualart--Peccati (2005), Theorem 1](https://arxiv.org/pdf/math/0503598)
equates fixed-chaos Gaussian convergence with vanishing nontrivial
contractions after variance normalization;
[Nualart--Ortiz-Latorre (2007), Theorem 7](https://arxiv.org/pdf/math/0703240)
gives joint Gaussian convergence from componentwise convergence and limiting
identity covariance. Their hypotheses match the kernels here. The common
Hilbert space may be taken as `ell²`, embedding each finite-dimensional
Gaussian input. Include the degree-one field `G_i` and normalize each
higher chaos by `sqrt(r!)`; distinct chaos degrees give zero covariance.

The already audited multilinear replacement transfers every fixed finite
chaos list to Rademacher inputs. The already audited even tail-transport
bound then permits the truncation degree to grow *after* the dimension
limit. Consequently the stronger marginal conclusion is justified:

\[
\sum_j B_{ij}S_j h((BS)_j)
\Longrightarrow N(0,\mathbb E h(Z)^2).
\]

The bounded second moments suffice for convergence of expectations of
the linearly growing Lipschitz test `z↦max(a,|z|)`. Uniformity over root
coordinates follows by applying the same argument to an arbitrary chosen
sequence of coordinates; every contraction and approximation estimate was
uniform. No Gaussian limit for `B f(BS)` is implied or needed.

Therefore the stronger response floor is

\[
L(t)=a(t)[2\Phi(a(t)/\sqrt{b(t)})-1]
       +2\sqrt{b(t)}\phi(a(t)/\sqrt{b(t)}),
\quad a(t)=2\phi(t),\ b(t)=2\Phi(t)-1,
\]

after taking smoothing to zero. Combining with the independently verified
baseline gives the valid asymptotic lower-bound expression

\[
\frac{(1-p)^2a(t)b(t)+p(1-p)L(t)}{1+p^2}.
\]

The numerical certification of a stated decimal is delegated to the
parent/algebra; the structural audit finds no fatal gap in this stronger
argument. This improves a lower bound only and does not resolve convergence.
