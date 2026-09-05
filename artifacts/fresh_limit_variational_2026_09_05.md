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
