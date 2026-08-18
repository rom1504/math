# Exact response completeness of normed affine positive-part carriers

Status: **rigorous exact-algebra obstruction, with an essential precision
qualification**.  Even degree-one carriers with unit mean and an arbitrarily
small fixed excess over unit `L^2` norm contain a finite family whose response
values determine the complete bridge-pressure table.  Thus positivity
clipping does not merely restore a few higher Walsh levels: in exact
arithmetic it restores every row level, and tensoring the rows restores every
bridge level.

This theorem applies directly to the actual optimized-child pressure `L`; no
surrogate pressure is introduced.  It concerns the **whole carrier response
surface**, not only its optimized scalar value.  Its inverse has exponentially
small singular values, so the theorem does not rule out an `o(N)`-accurate
evaluation of the optimized value from a smaller actual-child state.

Throughout, `U_n` is the fair law on `{+-1}^n`, and Walsh coefficients use

```math
\widehat f(S)=E_{U_n}f(B)\chi_S(B).
```

## 1. A diffuse affine carrier basis

Fix `K>1`.  For `n>=3`, put

```math
a_n={1\over n-1},
\qquad
\delta_n={1\over n-1},
\qquad
Z_n=1+2^{-n}\delta_n.                              \tag{AB.1}
```

For every `v in {+-1}^n`, define

```math
g_v(b)=1+a_n\langle v,b\rangle,
\qquad
h_v(b)={(g_v(b))_+\over E_{U_n}(g_v)_+}.            \tag{AB.2}
```

**Theorem AB.1 (normed affine response basis).**  For all sufficiently large
`n=n(K)`, every `h_v` belongs to the normed carrier class

```math
h={g_+\over E_Ug_+},
\qquad \deg g\le1,
\qquad E_Ug=1,
\qquad \|g\|_2\le K.                               \tag{AB.3}
```

Moreover the `2^n` functions `{h_v:v in {+-1}^n}` are a basis of the full
space of real functions on the row cube.

*Proof.*  The affine polynomial has

```math
E_Ug_v=1,
\qquad
\|g_v\|_2^2=1+na_n^2=1+{n\over(n-1)^2}.            \tag{AB.4}
```

The last quantity tends to one, so it is at most `K^2` for all sufficiently
large `n`.  If `b=-v`, then `g_v(b)=-delta_n`.  At every other vertex,
`<v,b>` is strictly larger than `-n`, and its smallest possible value is
`-n+2`; hence

```math
g_v(b)\ge1-a_n(n-2)={1\over n-1}>0.                \tag{AB.5}
```

Thus clipping changes exactly one row word:

```math
(g_v)_+=g_v+\delta_n\mathbf1_{\{-v\}},
\qquad E_U(g_v)_+=Z_n.                              \tag{AB.6}
```

Writing `v_S=prod_(j in S)v_j`, its Walsh coefficients are therefore

```math
\widehat h_v(S)=
\begin{cases}
1,&S=\varnothing,\\[2mm]
\displaystyle {a_n-2^{-n}\delta_n\over Z_n}\,v_j,
   &S=\{j\},\\[3mm]
\displaystyle {2^{-n}\delta_n\over Z_n}
 (-1)^{|S|}v_S,&|S|\ge2.
\end{cases}                                        \tag{AB.7}
```

Every coefficient multiplying `v_S` is nonzero.  The matrix
`(h_v(b))_(v,b)` is a group-circulant matrix on the Boolean group.  Its Walsh
eigenvalues are `2^n \widehat h_{\mathbf1}(S)`, all nonzero by (AB.7).
Consequently the matrix is invertible and the translated densities `h_v`
form a basis. `square`

These carriers remain quantitatively diffuse.  Since positive projection
can only decrease the squared norm here and `Z_n>=1`, one also has

```math
\|h_v\|_2\le\|g_v\|_2\le K.                       \tag{AB.8}
```

In particular, the basis property is not obtained by hiding point masses in
the carrier class.

## 2. Tensor completeness for the actual bridge pressure

Let `m,n` be arbitrary, and let

```math
L:\{+-1\}^{m\times n}\longrightarrow\mathbb R       \tag{AB.9}
```

be the exact bridge pressure.  In particular, `L` may be
`L_epsilon(B)=log Zbar_N(A,epsilon D,B;t)` for any two actual
contracted-temperature minimizing children.  For
`boldsymbol v=(v_1,...,v_m)`, put

```math
H_{\boldsymbol v}(B)=\prod_{i=1}^m h_{v_i}(B_i),
\qquad
\mathcal R_L(\boldsymbol v)=E_{U_B}[L(B)H_{\boldsymbol v}(B)]. \tag{AB.10}
```

**Corollary AB.2 (the exact carrier response is table-complete).**  For every
fixed `K>1` and all sufficiently large `n`, the collection

```math
\{\mathcal R_L(\boldsymbol v):
  \boldsymbol v\in(\{+-1\}^n)^m\}                  \tag{AB.11}
```

determines all `2^(mn)` values of `L(B)` by an invertible linear transform.
The same is true if each response value includes the row entropy penalty

```math
{1\over\lambda}\sum_iD(h_{v_i}U_n\Vert U_n).       \tag{AB.12}
```

*Proof.*  The products `H_boldsymbol v` are tensor products of the row basis
in Theorem AB.1, and hence form a basis of the `2^(mn)`-dimensional bridge
function space.  The response vector (AB.11) is the pairing of `L` with this
basis, so its transform matrix is the `m`-fold tensor power of the invertible
row matrix.  Finally, all `h_v` are Boolean translates of one another, so
their entropies are equal.  The term (AB.12) is therefore the same known
constant for every `boldsymbol v` and does not affect invertibility. `square`

For the actual collision--cavity interaction

```math
h(B)=L(B)-\sum_iL_i(B_i)+c,                         \tag{AB.13}
```

all response coefficients involving at least two nonconstant row modes are
unchanged by the row terms.  Hence the obstruction is genuinely a
cross-row obstruction, not an artifact of retaining the erased-row
likelihoods.

## 3. An unclipped exponential subspace

The one-point clipping in Theorem AB.1 is what makes the exact row span full.
There is a separate obstruction that does not use clipping.  Fix any
`0<a<min{1,sqrt(K^2-1)}`.  The `n+1` densities

```math
1,\qquad 1+ab_j\quad(1\le j\le n)                  \tag{AB.14}
```

are admissible degree-one carriers and span the constant-plus-linear row
space.  Their `m`-fold products span

```math
\left(\operatorname {span}\{1,b_1,\ldots,b_n\}\right)^{\otimes m},
\qquad
\dim=(n+1)^m.                                      \tag{AB.15}
```

Thus even before positive-part leakage, fixed row degree does not bound the
number of interacting rows.  Exact evaluation of the full response from
linear pressure moments naturally exposes an `(n+1)^m` cross-row tensor.

When `K>=sqrt2`, the admissible choices `1+sigma b_j` fix one selected bit
of each row.  The corresponding objective values give, without numerical
attenuation, the `2^m` values of the actual pressure averaged over all other
bridge bits.  This is an exact exponentially large slice even at comparable
child orders.

## 4. Precision and value scope

The exact obstruction must not be promoted into an approximate one.  From
(AB.7), every row Walsh eigenvalue of degree at least two has magnitude

```math
{2^{-n}\delta_n\over Z_n}=\Theta(2^{-n}/n).         \tag{AB.16}
```

The inverse row transform therefore amplifies some errors by
`Theta(n2^n)`, and its `m`-fold tensor power is still worse.  This is the
necessary price of reconciling table completeness with the uniform `L^2`
bound: Parseval prevents one diffuse density from having exponentially many
macroscopic Walsh coefficients.

Consequently Corollary AB.2 proves the following and no more.

1. The fixed-degree positive-part carrier restriction does **not** yield an
   exact algebraic response closure.  Its full objective surface contains
   the complete actual pressure table.
2. Retaining only the declared row degree of the presenting polynomial is
   invalid for exact evaluation; clipping can restore every row degree.
3. Even without clipping, cross-row tensor order remains unbounded.

It does **not** prove that the optimized scalar
`V_lambda^(d,K)` needs exponentially many observables at `o(N)` accuracy.
An optimizer-specific synchronization theorem could make the exponentially
attenuated directions irrelevant to that value.  Nor is an unrestricted
information lower bound meaningful when arbitrary real observables are
allowed: the two child sign matrices themselves have polynomial description
and determine `L` by exponential computation.  A stronger no-go theorem
would have to declare a precision/query model and exhibit actual minimizing
children with a robust value gap under that model.

The rigorous conclusion is therefore an exact-surface obstruction and a
scope falsifier for naive finite-Walsh closure, not a resolution of
finite-degree child closure at normalized pressure accuracy.

## 5. The complete basis is macroscopically attenuated

The precision caveat can be made quantitative using the actual pressure's
bit oscillation.  Let `t=beta/sqrt(N)`, `N=m+n`, and suppose that changing
one bridge bit changes `L` by at most `2t`, as it does for the actual
optimized-child pressure.  Cube Poincare gives

```math
\operatorname {Var}_{U_B}L
\le {1\over4}\sum_{i,j}E(L-L^{(ij)})^2
\le mn t^2={\beta^2mn\over N}.                     \tag{AB.17}
```

Write the nonconstant coefficient magnitudes in (AB.7) as `c_S`.  They all
satisfy

```math
|c_S|\le a_n={1\over n-1}.                         \tag{AB.18}
```

Expanding (AB.10) first in the bridge Walsh basis and then in the direction
variables gives

```math
\mathcal R_L(\boldsymbol v)
=\sum_{\boldsymbol S}\widehat L(\boldsymbol S)
 \prod_i c_{S_i}(v_i)_{S_i},                       \tag{AB.19}
```

with the convention `c_emptyset=1`.  Plancherel in the uniformly random
directions therefore proves

```math
\boxed{
\operatorname {Var}_{\boldsymbol v}\mathcal R_L(\boldsymbol v)
\le {1\over(n-1)^2}\operatorname {Var}_{U_B}L
\le {\beta^2mn\over N(n-1)^2}.}                   \tag{AB.20}
```

For comparable splits, the root-mean-square response spread over the entire
table-complete basis is only `O_beta(N^(-1/2))`.  The modes created solely by
clipping, namely row Walsh levels at least two, have the still smaller
multiplier `Theta(2^(-n)/n)` from (AB.16).

There is also a uniform, rather than average, macroscopic bound.  Equations
(AB.4), (AB.6), and `Z_n>=1` give

```math
\|h_v\|_2^2\le1+{n\over(n-1)^2}.
```

Consequently

```math
\begin{aligned}
\|H_{\boldsymbol v}-1\|_2^2
&=\prod_i\|h_{v_i}\|_2^2-1\\
&\le\exp\left\{{mn\over(n-1)^2}\right\}-1,
\end{aligned}                                      \tag{AB.21}
```

and Cauchy--Schwarz with (AB.17) yields

```math
\boxed{
|\mathcal R_L(\boldsymbol v)-E_{U_B}L|
\le
\beta\sqrt{mn\over N}
\left[
 \exp\left\{{mn\over(n-1)^2}\right\}-1
\right]^{1/2}.}                                    \tag{AB.22}
```

This is `O_beta(sqrt(N))=o(N)` uniformly at comparable splits.  Finally,
monotonicity of Renyi divergences gives

```math
\sum_iD(h_{v_i}U_n\Vert U_n)
\le m\log\left(1+{n\over(n-1)^2}\right)
=O(m/n).                                           \tag{AB.23}
```

Thus the entropy charge is only `O(1)` at comparable splits, while none of
these table-complete basis carriers can change the pressure by a linear
amount.  Exact reconstruction uses large signed linear combinations of
their responses; it is not implementable by selecting one admissible
product law.  This proves that AB.2 is a genuine exact-response obstruction
but not a robust obstruction to the normalized optimized value.

## 6. General high-row-degree invisibility

The attenuation is not confined to the special basis in Theorem AB.1.  Let
`F` be any row-effective pressure obtained by averaging `L` over arbitrary
genuine laws on all other rows.  It retains one-bit oscillation at most
`2u`.  Therefore

```math
\sum_{S\subseteq[n]}|S|\widehat F(S)^2
={1\over4}\sum_jE_U(F-F^{(j)})^2
\le nu^2.                                          \tag{AB.24}
```

For the projection onto row Walsh levels at least `k`, this gives

```math
\boxed{
\|\Pi_{\ge k}F\|_2\le u\sqrt{n/k},
\qquad
|E_h\Pi_{\ge k}F|\le K u\sqrt{n/k}}               \tag{AB.25}
```

for every density `h` with `||h||_2<=K`.  In particular, if every factor in
a row product satisfies

```math
\widehat{(h_i-1)}(S)=0\quad(1\le|S|<k),            \tag{AB.26}
```

then replacing its rows sequentially by fair rows and applying (AB.25) at
each genuine intermediate product proves

```math
\boxed{
|E_{\otimes_i h_i}L-E_UL|
\le m\sqrt{K^2-1}\,u\sqrt{n/k}.}                  \tag{AB.27}
```

At physical scale and comparable splits, taking `k` proportional to `n`
makes (AB.27) only `O_(beta,K)(sqrt N)`.  Hence a product phase cannot be
carried solely by row degrees proportional to the row length.

This conclusion is compatible with polynomially visible individual high
modes.  Fix `K>1`, put `R=sqrt(K^2-1)`, and let `S` have size `k`.  For all
sufficiently large `k`, set

```math
\ell=\left\lceil{k\over2}+{\sqrt k\over2R}\right\rceil,
\qquad s=\ell-{k\over2},
\qquad
g_S(b)=1+{1\over2s}\sum_{j\in S}b_j.              \tag{AB.28}
```

Then `Eg_S=1`, `||g_S||_2<=K`, and, for
`h_S=(g_S)_+/E(g_S)_+`, alternating partial-binomial identities give

```math
|\widehat h_S(S)|
={2^{-k}\over zs}{k-\ell\over k-1}
 {k-1\choose\ell-1}
=\Theta_K(k^{-1}),
\qquad z=E(g_S)_+.                                 \tag{AB.29}
```

Indeed `1<=z<=K`, while the binomial coefficient in (AB.29) lies a fixed
`K`-dependent number of standard deviations from the center.  Thus every
large subset can be targeted with only polynomial Fourier attenuation.
Nevertheless (AB.24) gives

```math
|\widehat F(S)|\le u\sqrt{n/k}.                    \tag{AB.30}
```

For `k` proportional to `n`, the isolated response in (AB.29) is at most
`O_(beta,K)(n^(-3/2))` per row.  Isolating it from the other modes also uses
a signed average of translated carrier responses, not one feasible carrier
law.  The exact response surface can therefore contain exponentially many
individually addressable directions while their high-degree-only physical
effect is still sublinear.

Equations (AB.24)--(AB.30) sharpen the remaining SML.  The unresolved
response information is not exact high row degree.  It is coherent
low/mesoscopic row-degree information coupled across a positive density of
rows, together with the entropy cost of the feasible factors.  The formal
cross-row tensor at every fixed row degree remains exponential, and no
actual-minimizer identity in this audit synchronizes it.
