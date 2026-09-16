# Independent audit of arXiv:2609.12427v1

Audit date: 2026-09-16 UTC. This audit was conducted directly from the complete primary paper, without accepting the result on the strength of publication, author attribution, or the paper's AI acknowledgment. The mathematical reconstruction and the parameter refinement below are independent checks, not a claim of novelty to the author.

## Outcome

I found no mathematical gap in the stated proof. In particular, the two mixed-norm inequalities have the correct direction, complex-valued functions cause no hypercontractivity problem, freezing coordinates preserves the class over which the fixed-dimensional supremum is taken, and the displayed contraction constant is actually about `0.674204810873338`, below the claimed `3/4`.

There are two separate outputs:

1. The stated polynomial Boolean Bohnenblust–Hille theorem and its level-flat influence corollary withstand the audit.
2. A fully quantified refinement of the same proof gives a weighted estimate with denominator `r^5`, degree growth `m^3`, and consequently ordinary BH growth `O(m^(11/2))`. Its influence corollary has exponent 11. This is a secondary analytic observation, **not progress on the asymptotic quadratic signing constant**.

Indeed, a concrete six-variable signing independently enumerated below prevents any unrestricted dimension-free quadratic BH constant, however sharpened, from giving a direct asymptotic signing lower bound greater than `0.390057886553453`. The reported project lower bound is already `0.4333221116640807`.

## Source and preserved inputs

The primary source is Paata Ivanisvili, [*Polynomial growth of Bohnenblust–Hille constants on the Hamming cube*, arXiv:2609.12427v1](https://arxiv.org/html/2609.12427v1), submitted 11 September 2026. The entire paper was read, including both remarks, the influence application, acknowledgments, and bibliography. The arXiv page identifies a CC BY 4.0 license.

Preserved primary downloads:

- `tmp/bh_2026_09_16/2609.12427v1_independent.pdf`, SHA-256 `e6d21b8df03c1d7a9f28f9767c8b4b4271b25a14e01d061e07d87bd1439075a3`.
- `tmp/bh_2026_09_16/2609.12427v1_independent.html`, SHA-256 `247a3aa14cdb3a80317280ccbf4418313dc8bb042ffa1b47df3340655a5a29b7`.
- `tmp/bh_2026_09_16/independent_contraction_check.py`, exact rational/integer verification of the refinement and exact enumeration of the six-variable obstruction.

Canonical replay: `computations/bh_2026_09_16_independent_audit_checks.py`; its preserved machine-readable output is `computations/results/bh_2026_09_16_independent_audit_checks.json`. Both the scratch replay and canonical replay were rerun successfully with the workspace's `.venv/bin/python`; all 31 exact assertions passed.

The paper's low-degree input is the existence of fixed-degree Boolean BH constants. Its cited [Defant–Mastyło–Pérez primary preprint](https://arxiv.org/abs/1706.03670) states the real-valued theorem. Passing to complex-valued polynomials is harmless: apply the real theorem separately to the real and imaginary parts and use the triangle inequality in coefficient `ell_q`. A factor of 2 suffices, and only finiteness in degrees 1 through 11 is needed.

## 1. Exact scope of the result

The probability space is the cube `{-1,1}^n` with uniform probability measure. Write

\[
f(x)=\sum_{S\subseteq[n]}a_Sx^S,
\qquad q_r=\frac{2r}{r+1},
\qquad L_r(f)=\left(\sum_{|S|=r}|a_S|^{q_r}\right)^{1/q_r}.
\]

For every positive integer degree bound `m`, uniformly over the number of coordinates, the paper proves

\[
\|a\|_{q_m}\le K_1m^{27}\|f\|_\infty,
\qquad
\left(\sum_{r=1}^m\frac{L_r(f)^2}{r^{10}}\right)^{1/2}
\le K_2m^{22}\|f\|_\infty.
\]

The functions may be complex-valued and inhomogeneous. The first expression includes the constant coefficient; the second does not. The constants are absolute but the paper does **not** provide a useful numerical value for them or a sharp fixed-degree constant. Empty Fourier levels vanish, so the argument also applies after freezing coordinates even if `m` then exceeds the remaining dimension.

The influence corollary concerns real-valued `f` bounded by 1, with the strong hypothesis

\[
|a_S|=a_{|S|}\quad\text{for every subset }S\text{ on each level}.
\]

This means a level is either entirely absent or has all its coefficients of the specified magnitude. It does not mean merely that the *nonzero* coefficients have a common magnitude. Under that hypothesis **every** coordinate has influence at least `c m^(-54) Var(f)^2`.

## 2. Reconstruction of the stated proof

### 2.1 Reduction to the weighted quantity

Normalize `||f||_infinity=1`. Since `q_r <= q_m`, the coefficient `ell_(q_m)` norm within level `r` is no larger than `L_r(f)`. Separating the constant term, comparing the finite sequence of level norms in `ell_(q_m)` and `ell_2`, and then removing the weights gives

\[
\|a\|_{q_m}
\le1+\left(\sum_{r=1}^mL_r(f)^{q_m}\right)^{1/q_m}
\le1+m^{1/(2m)}\left(\sum_{r=1}^mL_r(f)^2\right)^{1/2}
\le1+m^{5+1/(2m)}S(f),
\]

where `S(f)^2 = sum_(r=1)^m L_r(f)^2/r^10`. These are counting-measure coefficient norms and probability-measure function norms; confusing those conventions would introduce dimension factors, but the paper does not do so.

### 2.2 Low levels and the complex Markov point

For a vertex `x`, the one-variable polynomial

\[
p_x(t)=f(tx)=\sum_{r=0}^m f^{=r}(x)t^r
\]

has supremum at most 1 on `[-1,1]`. Indeed, the multilinear extension of `f` is bounded by its vertex supremum: successively maximize the convex function `|A+Bt|` at an endpoint in each coordinate. This is valid for complex `A,B`.

Repeated Markov inequalities give

\[
|f^{=r}(x)|=|p_x^{(r)}(0)|/r!\le m^{2r}/r!.
\]

The real Markov inequality extends with exactly the same constant to complex polynomials. At any point choose a phase making the derivative real and nonnegative, apply the real inequality to the real part of the rotated polynomial, and then take the supremum. Repeating this for derivatives is legitimate.

Let `C_r` be any valid finite homogeneous Boolean BH constant, for complex-valued functions if required. Then

\[
L_r(f)\le C_r\|f^{=r}\|_\infty\le C_rm^{2r}/r!.
\]

For the levels 1 through 11,

\[
H(f)\le K_0m^{22},
\qquad
K_0^2=\sum_{r=1}^{11}\frac{C_r^2}{r^{10}(r!)^2}.
\]

This is an input from finitely many old BH inequalities, not an application of the new theorem being proved. Thus it is not circular. If `m<12`, this already controls the whole weighted sum.

### 2.3 Fixed-dimensional supremum and restriction stability

For fixed finite `m,n`, set

\[
M=\sup\{S(g):\deg g\le m,\ \|g\|_\infty\le1\}.
\]

This is finite without a dimension-free theorem: each coefficient is bounded by 1 and there are finitely many coefficients. Explicitly,

\[
M^2\le\sum_{r=1}^{\min(m,n)}r^{-10}\binom nr^{2/q_r}.
\]

Every coordinate restriction of an eligible polynomial is still bounded by 1 and has degree at most `m`. Re-embed it on the original cube by adding unused coordinates. Its weighted norm is unchanged and bounded by this same `M`. The proof never assumes the all-dimensional supremum is finite.

### 2.4 Random partition and captured coefficient mass

Independently send each coordinate to `I` or `J` with equal probability. On a fixed level of size `r`, the number going to `I` is binomial `(r,1/2)`. Keep only those coefficients with split degree `d` in `[r/6,5r/6]`.

For `r>=12`, the capture probability satisfies

\[
p_r\ge p:=1-2e^{-8/3}=0.861033097554397.
\]

If `C_r(omega)` denotes the `ell_(q_r)` norm of the captured coefficients, linearity gives the exact identity

\[
\mathbb E_\omega C_r(\omega)^{q_r}=p_rL_r(f)^{q_r}.
\]

No independence between the indicators belonging to different subsets is used or needed. Probability-space norm monotonicity gives

\[
L_r(f)\le p_r^{-1/q_r}(\mathbb E C_r^2)^{1/2}.
\]

Because `q_r>=24/13`, one may replace the factor by `p^(-13/24)`. Squaring and summing over `r` commutes with the finite expectation.

### 2.5 Separating bidegrees

For a fixed partition write `c_(A,B)=a_(A union B)` and let `b_(d,e)` be its coefficient `ell_(q_(d+e))` norm restricted to `|A|=d, |B|=e`. The captured norm on level `r` is the `ell_(q_r)` norm of those blocks. There are at most `r+1` blocks, hence

\[
C_r\le(r+1)^{1/(2r)}\left(\sum_{d+e=r\atop\text{balanced}}b_{d,e}^2\right)^{1/2}.
\]

The function `log(r+1)/r` decreases for positive `r`; thus for `r>=12` the prefactor is at most `kappa=13^(1/24)=1.112792798278971`.

### 2.6 Two-block interpolation, with exponents checked

Put `r=d+e`, `theta=d/r`, and

\[
X_{d,e}=\left(\sum_A(\sum_B|c_{A,B}|^2)^{q_d/2}\right)^{1/q_d},
\quad
Y_{d,e}=\left(\sum_B(\sum_A|c_{A,B}|^2)^{q_e/2}\right)^{1/q_e}.
\]

For each fixed row, interpolate its `ell_(q_r)` norm between `ell_2` with weight `theta` and `ell_(q_e)` with weight `1-theta`. Then use Hölder in the row index with exponents `q_d/(theta q_r)` and `2/((1-theta)q_r)`. The necessary identities are

\[
\frac1{q_r}=\frac\theta2+\frac{1-\theta}{q_e}
=\frac\theta{q_d}+\frac{1-\theta}2.
\]

This yields `b_(d,e) <= X_(d,e)^theta Z_(d,e)^(1-theta)`, where `Z` is the `ell_A^2(ell_B^(q_e))` norm. Minkowski gives `Z<=Y`, since `q_e<=2`. Therefore

\[
b_{d,e}\le X_{d,e}^{\theta}Y_{d,e}^{1-\theta}.
\]

The direction of the mixed-norm swap is essential and correct.

### 2.7 Hypercontractive damping has uniformly bounded removal cost

Let

\[
\rho_d=\sqrt{q_d-1}=\sqrt{(d-1)/(d+1)},
\quad \widehat X_{d,e}=\rho_d^eX_{d,e},
\quad \widehat Y_{d,e}=\rho_e^dY_{d,e}.
\]

Balanced pairs satisfy `d,e>=2`, so these numbers are nonzero. Removing the damping costs exactly

\[
D(d,e)=\left[\frac{d+1}{d-1}\frac{e+1}{e-1}\right]^{de/(2(d+e))}.
\]

The function

\[
k\log\frac{k+1}{k-1}
=2\sum_{j\ge0}\frac1{(2j+1)k^{2j}}
\]

decreases for `k>1`. Since `d,e>=2`, its value is at most `2 log 3`. It follows that `log D <= log 3`, so

\[
b_{d,e}\le3\widehat X_{d,e}^{\theta}\widehat Y_{d,e}^{1-\theta}.
\]

### 2.8 The weight creates the contraction

Define `h(theta)=theta^theta(1-theta)^(1-theta)` and

\[
U_{d,e}=\widehat X_{d,e}/d^5,
\qquad V_{d,e}=\widehat Y_{d,e}/e^5.
\]

Dividing by `(d+e)^5` gives exactly

\[
\frac{b_{d,e}}{(d+e)^5}
\le3h(\theta)^5U_{d,e}^{\theta}V_{d,e}^{1-\theta}.
\]

On `[1/6,5/6]`, convexity and symmetry of `log h` place its maximum at the endpoints. That maximum is `(5^5/6^6)^(1/6)=0.637270409443108 < 2/3`. After squaring and weighted AM–GM,

\[
\frac{b_{d,e}^2}{(d+e)^{10}}
\le9(2/3)^{10}\bigl(\theta U_{d,e}^2+(1-\theta)V_{d,e}^2\bigr)
\le9(2/3)^{10}(U_{d,e}^2+V_{d,e}^2).
\]

The last step is deliberately lossy, not a gap.

### 2.9 Full rows, Minkowski, and complex-valued hypercontractivity

Fix `d>=2`. For each `A` of size `d`, collect **all** complementary degrees:

\[
R_A(y)=\sum_{B\subseteq J}c_{A,B}y^B.
\]

With `z_(A,e)^2 = sum_(|B|=e) rho_d^(2e)|c_(A,B)|^2`, Minkowski gives

\[
\|z\|_{\ell_e^2(\ell_A^{q_d})}
\le\|z\|_{\ell_A^{q_d}(\ell_e^2)}.
\]

For clarity, raising this inequality to the power `q_d` turns it into the triangle inequality in `ell_(2/q_d)` applied to the sequences `z_(A,e)^(q_d)`. Since `2/q_d>=1`, its direction is confirmed. Adding the nonnegative complementary degrees 0 and 1 only enlarges the right side.

Parseval and Bonami–Beckner now yield

\[
\left(\sum_{e\ge2}\widehat X_{d,e}^2\right)^{1/2}
\le\left(\sum_{|A|=d}\|T_{\rho_d}R_A\|_2^{q_d}\right)^{1/q_d}
\le\left(\mathbb E_y\sum_{|A|=d}|R_A(y)|^{q_d}\right)^{1/q_d}.
\]

There is no missing real-to-complex factor. For real `0<=rho<=1`, the cube noise operator is a positive probability kernel. Consequently

\[
|T_\rho R|\le T_\rho|R|,
\qquad
\|T_\rho R\|_2
\le\|T_\rho|R|\|_2
\le\||R|\|_q=\|R\|_q
\]

whenever `rho^2<=q-1`. The real inequality may be applied to `|R|`: it holds for all functions on the cube, not only bounded-degree polynomials. This establishes the complex-valued claim without invoking “complex hypercontractivity” with a complex noise parameter.

For fixed `y`, the values `R_A(y)` are exactly the coefficients of the restricted polynomial `f(.,y)` at degree `d`. Thus the last displayed quantity is

\[
\bigl(\mathbb E_yL_d(f(\cdot,y))^{q_d}\bigr)^{1/q_d}
\le\bigl(\mathbb E_yL_d(f(\cdot,y))^2\bigr)^{1/2}.
\]

Square, multiply by `d^(-10)`, sum, and use restriction stability:

\[
\sum_{d,e\ge2}U_{d,e}^2\le M^2.
\]

Interchanging the two coordinate blocks gives `sum V_(d,e)^2 <= M^2`.

Keeping the complete row is indispensable. Replacing `R_A` by a homogeneous component before this step would introduce an unjustified supremum-norm assertion; the paper correctly avoids that replacement.

### 2.10 Closing the fixed-dimensional contraction and constants

Combine the capture factor, the block-count factor, the weight bound, and the two row bounds. The high-level part satisfies

\[
T(f)\le c_0M,
\quad
c_0=(1-2e^{-8/3})^{-13/24}13^{1/24}3\sqrt2(2/3)^5
=0.674204810873338\ldots<3/4.
\]

Therefore

\[
M\le K_0m^{22}+c_0M,
\qquad M\le\frac{K_0}{1-c_0}m^{22}.
\]

The paper replaces the reciprocal by 4. More precise numerical bookkeeping gives `1/(1-c_0)=3.069413034245951`. Dimension disappears only after this finite-dimensional rearrangement, which is valid.

For a wholly rational certification of the original `c_0<3/4`, rather than relying on decimals, the eighth-order exponential Taylor sum at `8/3` exceeds `100/7`; hence `p>43/50`. Also `13*8^24<9^24` implies `kappa<9/8`, and `sqrt(2)<10/7`. Since `p^(-13/24)<p^(-1)`, these give

\[
c_0<\frac{50}{43}\frac98\cdot3\frac{10}{7}\left(\frac23\right)^5
=\frac{2000}{2709}<\frac34.
\]

Returning to the coefficient estimate and using `m^(1/(2m))<=2` gives, for example,

\[
K_2=K_0/(1-c_0),
\qquad K_1=1+2K_0/(1-c_0).
\]

These choices are explicit in terms of the eleven fixed homogeneous constants, not explicit universal numerals independent of unreported input constants.

### 2.11 Influence corollary

Let `v_r = binom(n,r) a_r^2`. The level-flat hypothesis gives

\[
\operatorname{Var}(f)=\sum_rv_r,
\qquad
\operatorname{Inf}_i(f)=\frac1n\sum_rrv_r
\quad\text{for every }i.
\]

Also

\[
L_r(f)^2=\binom nr^{1/r}v_r\ge\frac nr v_r.
\]

The binomial inequality follows directly from
`binom(n,r)=product_(j=0)^(r-1)(n-j)/(r-j)`, whose factors are all at least `n/r`.

Cauchy–Schwarz now gives

\[
\operatorname{Var}(f)^2
\le\left(\sum_rrv_r\right)\left(\sum_rv_r/r\right)
\le\operatorname{Inf}_i(f)\sum_rL_r(f)^2
\le K_2^2m^{54}\operatorname{Inf}_i(f).
\]

Thus the corollary is correct, including its “every coordinate” quantifier and exponent 54. Its coefficient-support hypothesis must not be silently weakened.

## 3. Independent all-degree parameter refinement

This section changes parameters in the already audited proof. It is not necessary for the validity of the paper's stated conclusions.

### 3.1 Better fixed-level extraction without a deep auxiliary theorem

If a one-variable polynomial `p` has degree at most `m` and `||p||_[-1,1]<=1`, then

\[
|[t^r]p(t)|\le(em/r)^r\quad(1\le r\le m).
\]

Here is a self-contained proof. The expression

\[
G(w)=w^mp((w+w^{-1})/2)
\]

is an ordinary polynomial and is bounded by 1 on the unit circle; the maximum principle bounds it by 1 in the disk. For any complex `z`, choose a root of `w+w^(-1)=2z` with `|w|<=1`. Then

\[
|p(z)|\le |w|^{-m}
\le\bigl(|z|+\sqrt{1+|z|^2}\bigr)^m
\le e^{m|z|}.
\]

Cauchy's coefficient estimate on `|z|=R` gives `|[t^r]p|<=e^(mR)R^(-r)`. Set `R=r/m`. Applying this to each `p_x` proves

\[
\|f^{=r}\|_\infty\le(em/r)^r\|f\|_\infty.
\]

Thus retaining only levels 1–3 as initial cases costs `m^3`, not `m^6`.

There is also an independent real-variable route to precisely the low levels needed here. Set `g(theta)=p(cos theta)`, a trigonometric polynomial of degree at most `m`. Trigonometric Bernstein gives `||g^(j)||_infinity<=m^j`. At `theta=pi/2`, the chain rule gives `p'(0)=-g'`, `p''(0)=g''`, and `g'''=p'(0)-p'''(0)`. Consequently the first three coefficients are bounded by `m`, `m^2/2`, and `(m^3+m)/6`, respectively. These sharper constants give the same `m^3` low-level growth. The complex-plane argument above avoids requiring this additional form of Bernstein's inequality.

### 3.2 General weighted contraction formula

Use

\[
S_s(f)^2=\sum_{r=1}^mL_r(f)^2/r^{2s}.
\]

For each high level choose balanced bidegrees, with exact capture probability `p_r`, block count `N_r`, and `theta=d/r`. The preceding proof, with the weighted AM–GM coefficients retained, gives

\[
T_s(f)^2
\le\mathbb E_\omega\sum_{(d,e)}
p_r^{-2/q_r}N_r^{1/r}D(d,e)^2h(\theta)^{2s}
\bigl(\theta U_{d,e}^2+(1-\theta)V_{d,e}^2\bigr),
\]

where `U=widehat X/d^s` and `V=widehat Y/e^s`. The row argument still gives `sum U^2<=M_s^2` and `sum V^2<=M_s^2` for **any** fixed `s`. Hence a sufficient all-level contraction criterion is

\[
p_r^{-1/q_r}N_r^{1/(2r)}D(d,e)h(\theta)^s
\sqrt{2\max(\theta,1-\theta)}\le c<1
\]

for every retained bidegree. It is legitimate to bound this coefficient pairwise before summing; no product of two unrelated suprema is being treated as an equality.

### 3.3 Choice `s=5/2`, high cutoff 4, balance interval `[1/3,2/3]`

Retain `ceil(r/3)<=d<=floor(2r/3)`, set `e=r-d`, and take `r>=4`. All retained degrees are at least 2. Define

\[
p_r=2^{-r}\sum_{d=\lceil r/3\rceil}^{\lfloor2r/3\rfloor}\binom rd,
\qquad N_r=\lfloor2r/3\rfloor-\lceil r/3\rceil+1.
\]

For `r=4`, only `(d,e)=(2,2)` occurs. With `p_4=3/8`, `N_4=1`, `D=3`, and `h=1/2`, the sufficient contraction factor is

\[
c_4=3(8/3)^{5/8}2^{-5/2}=(27/32)^{1/8}
=0.978986545550759\ldots<1.
\]

For the next three degrees it suffices to use `D<=3`, the exact probabilities and block counts, and their actual extreme split ratios. The resulting upper bounds have these exact powers:

| Level | Probability | Block count | Exact power of the upper contraction bound | Upper bound |
| --- | --- | --- | --- | --- |
| 5 | `5/8` | 2 | `c_5^5 = 2^17 3^15 / 5^18` | `0.868107950729711` |
| 6 | `25/32` | 3 | `c_6^12 = 2^67 / (3^23 5^14)` | `0.892901184699018` |
| 7 | `35/64` | 2 | `c_7^14 = 2^110 3^29 / (5^8 7^50)` | `0.862858461812393` |

Each is strictly less than `9/10`, checked by exact rational arithmetic. Also `9/10<c_4` by raising both sides to the eighth power.

For `r>=8`, both split degrees are at least 3. The monotonic damping argument now gives `D<=2^(3/2)`, since `3 log((3+1)/(3-1))=3 log 2`. We also have

\[
h(\theta)\le(4/27)^{1/3},
\quad\max(\theta,1-\theta)\le2/3,
\quad N_r^{1/(2r)}\le9^{1/16},
\quad1/q_r\le9/16.
\]

The all-degree probability bound is

\[
p_r\ge21/32\qquad(r\ge8).
\]

This is not based on checking an unbounded sequence numerically. The exact finite checks for `8<=r<=31` are in the companion script, with the minimum `21/32` attained at `r=10`. For every `r>=32`, Hoeffding gives

\[
p_r\ge1-2e^{-r/18}\ge1-2e^{-16/9}>21/32.
\]

The final strict inequality has an elementary rational certificate:

\[
e^{16/9}>\sum_{j=0}^5\frac{(16/9)^j}{j!}
=\frac{5189207}{885735}>\frac{64}{11}.
\]

Consequently all `r>=8` have contraction factor at most

\[
c_{\rm tail}=(32/21)^{9/16}9^{1/16}2^{3/2}
(4/27)^{5/6}\sqrt{4/3}
=0.967089238922440\ldots<c_4.
\]

The comparison is exactly certified: `c_tail^48 = 2^335/(3^165 7^27)`, and `c_tail<c_4` is equivalent to

\[
2^{365}<3^{183}7^{27}.
\]

Thus `c=c_4` is a rigorously uniform contraction factor over **every** high level.

### 3.4 Refined theorem and constants

Let

\[
D_0^2=\sum_{r=1}^3\frac{C_r^2(e/r)^{2r}}{r^5},
\qquad C_* = D_0/(1-(27/32)^{1/8}).
\]

The low-level estimate and the contraction prove

\[
\boxed{\left(\sum_{r=1}^mL_r(f)^2/r^5\right)^{1/2}
\le C_*m^3\|f\|_\infty.}
\]

The numerical reciprocal in this explicit dependence is `1/(1-c)=47.588558198060895`. Removing the weights as before yields

\[
\boxed{\|\widehat f\|_{2m/(m+1)}
\le(1+2C_*)m^{11/2}\|f\|_\infty.}
\]

The same level-flat corollary calculation gives

\[
\boxed{\operatorname{Inf}_i(f)\ge C_*^{-2}m^{-11}\operatorname{Var}(f)^2
\quad\text{for every }i.}
\]

These statements have the same real/complex and dimension-uniform scope as appropriate in the original theorem. They are consequences of a parameter refinement and better low-level extraction, not numerical conjectures. They are not claimed to be optimal or previously unknown.

## 4. Quadratic-project obstruction independently verified

For a complete signing on six variables, use

\[
H=\begin{pmatrix}
0&1&1&1&1&1\\
1&0&1&1&-1&-1\\
1&1&0&-1&-1&1\\
1&1&-1&0&1&-1\\
1&-1&-1&1&0&1\\
1&-1&1&-1&1&0
\end{pmatrix}.
\]

This is the matrix stored in `computations/results/exact_m6.json`, but the result below was checked independently by summing all 15 unordered edge terms on all 64 spin vectors. The energy histogram is

\[
\{-5:12,\ -3:20,\ 3:20,\ 5:12\}.
\]

Thus `Q(H)=5`, while the coefficient `ell_(4/3)` norm is `15^(3/4)`. No assertion that this matrix is globally optimal is required.

Let `B_2^hom` denote any dimension-free constant valid for **all** homogeneous quadratic Boolean polynomials. The example forces

\[
B_2^{\rm hom}\ge15^{3/4}/5=1.524398244463844\ldots.
\]

For an `n`-vertex complete signing, its coefficient norm is `binom(n,2)^(3/4)`. The ordinary BH inequality alone therefore yields

\[
Q(A)\ge\binom n2^{3/4}/B_2^{\rm hom}.
\]

Even the best possible unrestricted constant in this argument has asymptotic coefficient at most

\[
\frac{2^{-3/4}}{15^{3/4}/5}
=\frac5{30^{3/4}}
=0.390057886553453\ldots.
\]

This rules out improving the reported `0.4333221116640807` lower constant by **merely sharpening the unrestricted degree-two BH constant**. It does not rule out an asymptotic-only result for dense signings, a structural strengthening beyond BH, or a genuinely different high-degree reduction. The distinction between these quantifiers is essential.

In fact the finite witness can be embedded in every ambient dimension `n>=6` by adding zero coefficients. Thus even dimension-dependent constants valid for all quadratic polynomials on the `n`-cube satisfy the same lower bound for every `n>=6`. Merely excluding small ambient dimensions is not an escape. One must restrict to dense/flat coefficient classes, forbid such padding, or introduce other structure outside the paper's unrestricted class.

The level-flat influence corollary does apply to a normalized complete quadratic signing, but its unspecified constants do not circumvent this obstruction. Disjoint tensor products generally cease to be level-flat on the full ambient cube, so they cannot automatically inherit the corollary's hypothesis.

## 5. Independent audit of the reduced-power application

The separate derivation in `artifacts/bh_2026_09_16_power_analysis.md` was read in full. Its central identities, the minimax obstruction, the auxiliary moment argument, and the displayed second/third-power coefficient formulas were independently checked.

Write `g=H_A^k`, let `c` be its reduced Walsh coefficients, and let `L=||H_A||_(2k)`. For `p=4k/(2k+1)` define `F=||c||_p^(1/k)`. If `pi_S=|c_S|^2/sum_T|c_T|^2` and `alpha=p/2`, then direct substitution, with no inequality, gives

\[
\log(F/L)=\frac1{kp}\log\sum_S\pi_S^\alpha
=\frac{H_\alpha(\pi)}{4k^2}.
\]

Consequently `L<=F<=D^(1/(4k^2))L`, where `D` is any upper bound on the Fourier support size. A valid degree bound `m` with exponent `q_m` instead gives entropy denominator `2mk`. The expression `min(2k,2 floor(n/2))` is always such a degree bound; it need not equal the actual degree. For example, the six-variable witness satisfies `H^4=34H^2-225`, so its fourth power has degree 4 rather than the available bound 6. Using a genuinely smaller degree does not remove the support-norm comparison.

There is a sharper support bound for the complete weighted functional than the coarse `sqrt(en)` estimate. For any function `g` invariant under global sign reversal, and any `s>=0`,

\[
S_s(g)^2
\le\sum_{r\ge2\atop r\text{ even}}r^{-2s}\binom nr^{1/r}\|g^{=r}\|_2^2
\le\frac{\binom n2^{1/2}}{2^{2s}}\|g-\mathbb Eg\|_2^2.
\]

The last inequality uses the fact that `binom(n,r)^(1/r)` decreases with `r`: binomial coefficients are log-concave and the zeroth coefficient is 1, so the averages of their successive logarithmic increments decrease. Therefore

\[
S_s(H_A^k)^{1/k}
\le\left(\binom n2^{1/4}/2^s\right)^{1/k}\|H_A\|_{2k}.
\]

This proves directly that when `k/log n -> infinity`, the full weighted coefficient functional cannot exceed the ordinary moment by a nonvanishing multiplicative advantage after the `k`th root. It does not assert a two-sided weighted comparison.

The stronger sublinear minimax obstruction is also correct. Average over independent uniform edge signs. For each fixed spin vector, `H_A(x)` is a sum of `N=binom(n,2)` independent signs, whence

\[
\mathbb E_A\mathbb E_xH_A(x)^{2k}
\le(2k-1)!!N^k\le(2kN)^k.
\]

The moment comparison follows coefficientwise from `cosh(t)<=exp(t^2/2)` as even power series; it is not an unjustified inference from mgf ordering. Some admissible signing therefore has `||H_A||_(2k)<=sqrt(2kN)`. The elementary support bound `D<=(n+1)^(2k)` gives

\[
\min_A F\le(en)^{1/(2k)}\sqrt{2kN}.
\]

The weighted functional satisfies the same weaker upper bound. For `2<=k<=K<=n`, the function `sqrt(k) exp(log(en)/(2k))` decreases and then increases, so its maximum occurs at an endpoint. Separating `K<=log(en)` from `K>=log(en)` gives an absolute bound

\[
\frac{\min_AF}{n^{3/2}},\quad
\frac{\min_AS_s(H_A^k)^{1/k}}{n^{3/2}}
\le C_s\left(n^{-1/4}+\sqrt{K/n}\right).
\]

Thus every sequence `2<=k_n=o(n)` gives zero in this normalization. The witness produced by averaging may depend on `k`; this is allowed under the minimization. It need not minimize the cap `Q`, and this argument makes **no assertion about the moment distribution of actual cap minimizers**. It obstructs universal coefficient lower bounds over all signings, not applications supplied with an additional theorem about minimizers.

At `k` proportional to `n`, the support comparison already supplies a bounded prefactor before the `k`th root. Both the new polynomial BH loss and the older subexponential BH loss become negligible after that root; indeed the elementary dimension-dependent comparison is already enough. At superlinear powers, the global coefficient functional and the cap have asymptotically matching pointwise sandwiches, making the corresponding global coefficient-minimum convergence problem equivalent to the original one. That equivalence should not be inferred for the weighted functional from its one-sided upper bound alone.

The power agent subsequently supplied the following valid extra lemma, which *does* recover weighted equivalence at superlinear powers. Fixed-degree BH gives `Q(A)>=c n^(3/2)` uniformly over complete signings. Put `tau=sqrt(2N)/Q(A)=O(n^(-1/2))`. Markov's inequality and Parseval give `P(|H_A|<=sqrt(2N))>=1/2`, while one maximum-absolute projective spin configuration has probability at least `2^(-(n-1))`. For independent copies `H,H'`, the two ordered cross-events in `Var(H^k)=(1/2)E(H^k-H'^k)^2` consequently imply

\[
\operatorname{Var}(H_A^k)
\ge2^{-n}\bigl(Q(A)^k-(2N)^{k/2}\bigr)^2.
\]

This is valid for odd as well as even `k`, since `|a^k-b^k|>=|a|^k-|b|^k`. For sufficiently large `n`, `tau<1` uniformly, so the displayed difference is positive. Using `S_s(g)>=m^{-s}||g-Eg||_2`, for any degree bound `m<=n`, gives

\[
S_s(H_A^k)^{1/k}
\ge m^{-s/k}2^{-n/(2k)}Q(A)(1-\tau^k)^{1/k}.
\]

Together with the earlier upper bound, this proves the uniform ratio `S_s(H_A^k)^(1/k)/Q(A)->1` when `k/n->infinity`. Hence weighted coefficient-minimum convergence is also equivalent to the original convergence problem in this regime, **with this extra variance argument included**.

The power note's independent-sign averaging, explicit Sylvester-family alternative, and direct decoupling/Gaussian mgf proof all check out. For the cubic coefficient at an edge, the direct count is `(3N-2)a_ij + 6[(A^3)_ij-(2n-3)a_ij]`, which simplifies to the stated formula. The four-vertex cubic formula correctly subtracts twice the star contribution, because the path-plus-edge correlation sum counts each star three times.

## 6. Reproduction and audit limits

Run from the workspace root:

```sh
.venv/bin/python computations/bh_2026_09_16_independent_audit_checks.py --json
```

The check uses standard-library integer and rational arithmetic for every certification; the printed decimals are only explanatory. It has no solver dependency, network access, or filesystem writes. It separately checks the primary proof's displayed numerical contraction and the six-variable energy histogram.

The original proof relies on the established real Bonami–Beckner inequality, Markov's inequality, and the existence of finitely many homogeneous Boolean BH constants. This audit reconstructs their use, supplies the real-to-complex justifications, and verifies all new contraction algebra; it does not re-prove the entire historical hypercontractivity or BH literature. No unsupported application to the unknown asymptotic quadratic signing limit is asserted.
