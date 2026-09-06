# A universal positive unmarked first variation under a fixed operator bound

Date: 2026-09-05. Independent proof responding to the parent's proposed
`A/-A`-averaged perturbation. This is a reusable strict-gain theorem for
each fixed smooth one-Gaussian baseline and each fixed operator bound.
It is not a numerical improvement over the banked `0.429` universal bound.

Chronological scope: that numerical reference is an earlier checkpoint.
The final explicit bound is `.4306581794055286`; the universal weighted
projection and normalized-gain files give the later strict improvement.
This supporting module is retained for its fixed-smooth-response argument,
not as the final proof-dependency map.

## 1. Matrix and response hypotheses

Let `A_n` be symmetric hollow sign matrices, `m=n-1`, and `B=A/sqrt(m)`.
Assume `||B||op≤L` for one fixed finite `L`. Let `S` be uniform on the sign
cube and put `G=BS`. All response functions and cutoff parameters below
are fixed before the matrix order tends to infinity.

Choose a nonzero odd smooth bounded function `F`, with bounded derivatives
of every order and

\[
 |F|\le1-3\eta,\qquad 0<\eta<1/4.
\]

For a standard Gaussian `N`, define

\[
 a=E[NF(N)],\quad r^2=EF(N)^2-a^2,
 \quad h(g)=\frac{F(g)-ag}{r}.                         \tag{1}
\]

Then `r>0`: a nonzero bounded function cannot equal a linear Gaussian
function almost everywhere. The response `h` is odd, has at most linear
growth, all its positive-order derivatives are bounded, and

\[
 Eh(N)^2=1,\qquad E[Nh(N)]=Eh'(N)=0.
\]

Define the new unmarked transport

\[
 Z=B h(G).
\]

This definition does not assume that `Z` is locally Gaussian or that
`B²=I`. Neither such assertion is used in this proof.

## 2. An exact positive Gaussian-input direction

First let the input be an independent standard Gaussian vector, so that
`G` is Gaussian with covariance `Q=B²`, whose diagonal is one. Expand
`F=sum_{s odd} f_s He_s/sqrt(s!)`. By (1), the Hermite coefficients of
`h` are `h_1=0` and `h_s=f_s/r` for odd `s≥3`. Exactly,

\[
 \begin{split}
 \frac1n E F(G)^T BZ
 &=\frac1n E F(G)^T Qh(G)\\
 &=\sum_{s\ge3\ {\rm odd}}\frac{f_s^2}{r}
                      \frac1n\sum_{i,j}Q_{ij}^{s+1}
 \ge r.                                                   \tag{2}
 \end{split}
\]

Every power in the last sum is even, and the diagonal alone contributes
one to each normalized trace sum. This proves positivity without a
spectral symmetry or coherence hypothesis.

One may pass the Hermite series by finite `L²` approximation: `||Q||op≤L²`
controls the bilinear energy error. Equivalently use the exact Gaussian
covariance kernel and then monotone convergence of its nonnegative
coefficient expression.

Equation (2) transfers to sign input with an `o(1)` error. Here is a direct
whole-expression Lindeberg proof. The functional is
`n^-1 F(Bs)^T Qh(Bs)`. Differentiating with respect to an input coordinate
uses a flat column `b=B e_k` with
`||b||₂=1`, `||b²||₂=O(n^-1/2)`, and `||b³||₂=O(n^-1)`.
The vector response and its first three derivatives consequently have
`L^p` Euclidean norms `O(sqrt(n)), O(1), O(n^-1/2), O(n^-1)` respectively,
uniformly through Gaussian/sign hybrids, with a fixed polynomial factor
in the coordinate currently being replaced. The third derivative of the
normalized functional is therefore `O_L(n^-3/2)` in expectation. A
second-order Taylor replacement and matching first two input moments give
total error `O_L(n^-1/2)`. The linear growth of `h` causes no difficulty,
since every fixed moment of the hybrid linear fields is uniformly bounded.
Thus for the original sign input,

\[
 \boxed{\qquad \frac1n E F(G)^T BZ\ge r-o(1).\qquad}  \tag{3}
\]

## 3. Moment and influence bounds for the unmarked field

The following bounds hold for every fixed finite `p`, uniformly in `i,k`:

\[
 \|Z_i\|_p\le C_{L,p,h},\qquad
 \|D_k Z_i\|_p\le C_{L,p,h}n^{-1/2}.                  \tag{4}
\]

Here `D_k` is the central sign-cube derivative. In particular the second
estimate includes the own coordinate `k=i`.

For completeness, for a weighted sum `W=sum_j w_j q(G_j)` with bounded
first two derivatives of `q`, its product-space gradient obeys

\[
 \sum_l|D_lW|^2
 \le C_{L,q}\left(\|w\|_2^2+
                              \frac{n\|w\|_1^2}{m^2}\right). \tag{5}
\]

The continuous-gradient contribution is bounded by
`||B||op² ||q'||infinity² ||w||₂²`. The difference between each central
derivative and the corresponding continuous derivative at the actual
field is at most `C ||w||₁/m`, by a second-order scalar Taylor bound.
This proves (5). Tensorized Poincare followed by the usual even-power
induction bounds every fixed centered `L^p` norm by the square root of
the right side, with a constant depending on `p` but not on dimension.

Apply (5) to `Z_i` using `w_j=B_ij`, so `||w||₂=1` and
`||w||₁≤sqrt(m)`. Its mean is exactly zero by oddness. Although `h` itself
has linear growth, its first two derivatives are bounded, so the same
argument applies.

For the derivative use

\[
 D_k Z_i=\sum_j B_{ij}B_{jk}h'(G_j)+O(n^{-1/2})
\]

pointwise, with uniform error. The weights have squared `l²` norm at most
`1/m` and `l¹` norm at most one. Also the common scalar mean
`E h'(G_j)=O(1/n)`: all sign-flat rows have the same normalized
Rademacher-sum law, whose first three moments match a standard Gaussian,
and fourth-order replacement uses `sum_l B_jl^4=1/m`, together with
`Eh'(N)=0`. Equation (5) and power induction now prove (4).

## 4. Truncation and an explicitly cube-feasible curve

Fix an even smooth cutoff `theta_R`, valued in `[0,1]`, equal to one on
`[-R,R]` and zero outside `[-2R,2R]`. Put

\[
 W_i=\theta_R(G_i)\theta_R(Z_i)Z_i.                    \tag{6}
\]

Thus `|W_i|≤2R`. By (4), uniformly in dimension,

\[
 \max_i E|Z_i-W_i|^2\longrightarrow0\quad(R\to\infty). \tag{7}
\]

For example split the two tail events and use bounded fourth moments,
the uniform subgaussian tails of `G_i`, and Markov on `Z_i`. The operator
bound gives

\[
 \left|\frac1nE F(G)^TB(Z-W)\right|
 \le L\sqrt{EF(G_i)^2}\,max_i\|Z_i-W_i\|_2.
\]

Choose `R` once and for all, depending on `L,F`, large enough that this
error is at most `r/4`. Equations (3)--(7) then imply

\[
 \frac1n E F(G)^T BW\ge3r/4-o(1).                    \tag{8}
\]

Choose `0<tau≤eta`. For `|epsilon|≤eta/(2R)`, define

\[
 F_\epsilon=F(G)+\epsilon W,\qquad
 H_\epsilon=1-\sqrt{F_\epsilon^2+\tau^2},\qquad
 \mu_\epsilon=F_\epsilon+S H_\epsilon.                \tag{9}
\]

Here `|F_epsilon|≤1-2eta`, `H_epsilon≥eta`, and
`|F_epsilon|+H_epsilon≤1`. Thus the conditional means in (9) are genuinely
in the cube, for both signs of epsilon. Independent coordinatewise
rounding preserves their expected hollow quadratic energies.

The order is important: first fix `F,eta,tau`; then choose `R` using the
operator-controlled tail estimate; only then choose the step epsilon.

## 5. The dangerous marked-marked derivative is negligible

Put

\[
 H_0(g)=1-\sqrt{F(g)^2+\tau^2},\qquad
 \psi(g)=\frac{F(g)}{\sqrt{F(g)^2+\tau^2}},\qquad
 g_i=\psi(G_i)W_i.
\]

We claim

\[
 \frac1n E[Sg]^T B[S H_0(G)]\longrightarrow0.          \tag{10}
\]

By (4), the smooth cutoff chain rule gives
`max_ik ||D_k g_i||₂=O(n^-1/2)`. Remove the own spin by setting
`g_i^0=E[g_i|S_l,l!=i]`. Then
`||g_i-g_i^0||₂=||D_i g_i||₂=O(n^-1/2)`. The operator bound makes the
normalized energy error from this replacement `O_L(n^-1/2)`.

Both `g_i^0` and `H_0(G_j)` now exclude their own spins exactly. Conditioning
on all other spins and writing the two endpoint expansions proves

\[
 E[S_i g_i^0 S_j H_0(G_j)]
   =E[(D_jg_i^0)(D_iH_0(G_j))].                       \tag{11}
\]

Conditional expectation does not increase the derivative `L²` norm, while
`|D_i H_0(G_j)|≤C|B_ji|=O(n^-1/2)`. Hence the right side of (11) is
`O(1/n)`. Multiplying by `|B_ij|=1/sqrt(m)`, summing the `n(n-1)` terms,
and dividing by `n` proves (10). No cancellation among bridge signs is
being assumed.

## 6. Averaging the two matrix orientations and differentiating

Let

\[
 e_\epsilon(B)=\frac1n E H_B(\mu_\epsilon(B)),\qquad
 L_\epsilon(B)=\frac12[e_\epsilon(B)+e_\epsilon(-B)].
\]

Each term is at most `Q(A)/(n sqrt(m))`, so their average is a legitimate
lower certificate. Under `B -> -B`, the fields satisfy

\[
 G\mapsto-G,\quad Z\mapsto Z,\quad W\mapsto W,
 \quad F_\epsilon\mapsto-F_{-\epsilon},
 \quad H_\epsilon\mapsto H_{-\epsilon}.
\]

Consequently the mixed `F/SH` derivative terms cancel in the average,
while the two self-energy derivatives remain. Exactly,

\[
 L'_0(B)=\frac1nE F(G)^T BW
       -\frac1nE[S\psi(G)W]^T B[S H_0(G)].            \tag{12}
\]

Equations (8) and (10) show `L'_0(B)≥3r/4-o(1)`, and in particular it
is at least `r/2` for all sufficiently large dimensions, uniformly in the
family under consideration.

At epsilon zero, the orientation average retains just the old one-probe
cross energy:

\[
 L_0(B)=\frac1n E F(G)^TB[S H_0(G)]
       \longrightarrow J_0=EF'(N)\,EH_0(N).          \tag{13}
\]

For a direct proof, expand the exact own-free endpoint `D_jF(G_i)`.
Its leading term is `B_ij F'(G_i)` and the accumulated remainder is
`O(n^-1/2)`. The average pair expectation factors because
`Tr Q²≤L²n` implies that all but a vanishing fraction of pairs have
`Q_ij=o(1)`, and the two-row Gaussian replacement is uniform for fixed
smooth responses. This establishes (13) independently of the full tree
hierarchy.

There is also a uniform second-derivative bound. Along (9),
`|mu_epsilon|≤1`, `|mu'_epsilon|≤4R`, and
`|mu''_epsilon|≤4R²/tau`. Therefore

\[
 |L''_\epsilon(B)|\le C_0,
 \qquad C_0=L(16R^2+4R^2/\tau).                      \tag{14}
\]

Choose the positive fixed step

\[
 \epsilon_0=\min\left\{\frac\eta{2R},\frac r{2C_0}\right\}.
\]

Taylor's theorem, (12)--(14), and Boolean rounding give the explicit
qualitative conclusion

\[
 \boxed{\quad
 \liminf_n\frac{Q(A_n)}{n^{3/2}}
 \ge J_0+\frac{r\epsilon_0}{4}>J_0.
 \quad}                                                    \tag{15}
\]

Every parameter is fixed, and the derivative and finite-step errors have
been controlled in the correct order.

## 7. Scope and the next actual bottleneck

The strict gain in (15) depends on the fixed operator bound `L`, the fixed
smooth baseline, and its slack and cutoffs. It is not presently uniform
as `L→infinity`. Therefore the spectral-deletion transfer principle alone
does not convert it into a numerical universal improvement over an
optimized baseline.

Smooth approximations to
`F_t(g)=sign(g)1{|g|>t}`, `H_t(g)=1{|g|≤t}` recover the first-Gaussian
baseline `2phi(t)[2Phi(t)-1]`, whose optimized value is approximately
`0.336493364431`. The gain proved for each fixed smooth approximation
need not dominate its approximation error to that optimum. No claim even
of a numerical improvement over that optimized constant is made here
without choosing and certifying actual parameters.

Most importantly, this one-Gaussian theorem does not yet couple the
unmarked residual to the successful many-tree `0.429` baseline. The exact
positive trace identity (2) uses a response of a *single linear Gaussian
field*. Extending that positivity, or replacing it by an appropriate
operator inequality, is the remaining substantive energy problem for a
many-field gain. Marginal Gaussianity alone is not a substitute.

## 8. Keeping the coherent channel: a stronger canonical direction

There is a useful extension of the same proof. Instead of transporting
only the nonlinear residual, define the own-spin-corrected field

\[
 Z^{\rm can}=B F(G)-aS,\qquad a=E[NF(N)].              \tag{16}
\]

This differs from `r Z` in Section 1 by the coherent channel
`a(Q-I)S`. It need not be Gaussian or have uniformly small individual
off-spin influences. Nevertheless it is a valid positive direction.

For Gaussian input,

\[
 \begin{split}
 \frac1n E F(G)^T BZ^{\rm can}
 &=\frac1nE F(G)^T QF(G)-aEF(N)N\\
 &=\sum_{s\ge1\ {\rm odd}}f_s^2
             \frac1n\sum_{i,j}Q_{ij}^{s+1}-a^2\\
 &\ge r^2+a^2\left(\frac{\operatorname{Tr}Q^2}{n}-1\right)
 \ge r^2.                                                 \tag{17}
 \end{split}
\]

The subtraction is exactly `a²`, since `BS=G`. The same whole-functional
Lindeberg proof transfers (17) to sign input with `o(1)` error. Thus the
fourth spectral moment defect contributes with a nonnegative sign.

The derivative calculation now gives, uniformly in `i,k`,

\[
 D_k Z_i^{\rm can}
   =a(Q-I)_{ik}+O_{L,p}(n^{-1/2})\quad\text{in }L^p.  \tag{18}
\]

Indeed subtract the common mean `E F'(G_j)=a+O(1/n)` from the weighted
derivative sum and use the same product-space estimate (5). Because
`Q_ii=1`, the own-spin derivative in (18) is still `O(n^-1/2)`.
Only the off-spin derivatives can retain a fixed coherent component.

Use the same two cutoffs, now applied to `G,Z^{can}`, and the same
cube-feasible curve (9). The fields still have all fixed moments bounded,
and the own-spin removal still costs `o(1)` in normalized energy. After
removal, the nuisance term in (11) is bounded by

\[
 C|B_{ij}|\bigl(|Q_{ij}|+n^{-1/2}\bigr).
\]

Multiplying by the energy coefficient `|B_ij|` and summing gives the
normalized upper bound

\[
 \frac C{nm}\sum_{i\ne j}
             (|Q_{ij}|+n^{-1/2})=O_L(n^{-1/2}),       \tag{19}
\]

since `sum_ij |Q_ij|≤n||Q||F≤L n^{3/2}`. Thus a maximum-influence bound
is unnecessary: the aggregate coherent influence is small enough.

The `A/-A` symmetry is unchanged, because both `B F(G)` and the correction
`aS` are invariant under reversing the matrix while keeping the input
spin fixed. The first-variation and finite-step theorem therefore hold
with the positive floor `r²` replacing `r`, while retaining the additional
nonnegative term in (17) before truncation. On a conference matrix this
direction reduces exactly to the residual direction already considered.

This is a genuine bounded-operator energy statement accommodating a
possibly non-Gaussian coherent channel. It still does not make the gain
uniform as the operator bound tends to infinity, nor by itself couple
the response to the successful multi-tree baseline.
