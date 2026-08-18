# Bipartite ground-state geometry sharpens the actual-child overlap obstruction

Status: **rigorous finite-order theorem for every actual-child prior**.  The
support-cardinality step in Theorem 37.52 treats the planted rank-one bridge
words as an arbitrary set of `2^(N-1)` cube points.  Their exact bipartite
geometry gives a strictly smaller Rademacher maximum.  A quantitative
Bernoulli-to-Gaussian comparison, Gaussian increment comparison, bounded
differences, and entropy transport lower the balanced strong-channel
threshold from

```math
\sqrt{8\log2}=2.354820\ldots
\quad\hbox{to}\quad
{4\over\sqrt\pi}=2.256758\ldots .
```

No property of a surrogate prior is used.  The theorem holds uniformly over
the exact rank-one prior induced by every pair of actual children and either
orientation.  It is still an obstruction theorem, not a directional product
certificate or a recurrence.

## 1. The Bernoulli bipartite ground state

Let `m+n=N`, `d=mn`, and

```math
\mathcal Q_{m,n}=\{xy^T:x\in\{-1,1\}^m,
                         y\in\{-1,1\}^n\}.
```

There are exactly `2^(N-1)` distinct matrices in `mathcal Q_(m,n)`.  For a
real `m` by `n` matrix `z`, put

```math
X(z)=\max_{Q\in\mathcal Q_{m,n}}\langle z,Q\rangle
    =\max_{x,y}x^Tzy.                               \tag{BG.1}
```

Write

```math
c_0=E|G|=\sqrt{2/\pi},
\qquad
c_{\rm Lin}={1+E|G|^3\over6}
            ={1+2\sqrt{2/\pi}\over6},             \tag{BG.2}
```

where `G` is standard Gaussian, and define the explicit error

```math
R_N=(N-1)(\log2)N^{1/3}
       +c_{\rm Lin}mnN^{-2/3}.                     \tag{BG.3}
```

### Lemma BG.1 (finite Bernoulli bipartite maximum)

If `B` has independent fair sign entries, then

```math
\boxed{
E X(B)
\le c_0\{m\sqrt n+n\sqrt m\}+R_N.}                \tag{BG.4}
```

*Proof.*  For `eta>0`, smooth the maximum by

```math
F_\eta(z)={1\over\eta}
 \log\sum_{Q\in\mathcal Q_{m,n}}
             e^{\eta\langle z,Q\rangle}.           \tag{BG.5}
```

Since `|mathcal Q_(m,n)|=2^(N-1)`,

```math
X(z)\le F_\eta(z)
\le X(z)+{(N-1)\log2\over\eta}.                   \tag{BG.6}
```

For any entry `e`, differentiation of the log partition gives

```math
\partial_e^3F_\eta
=\eta^2E_\eta[(Q_e-E_\eta Q_e)^3].                \tag{BG.7}
```

Because `Q_e` is sign valued,

```math
\left|E_\eta[(Q_e-E_\eta Q_e)^3]\right|
=2|u|(1-u^2)\le {4\over3\sqrt3}<1,
\qquad u=E_\eta Q_e.                              \tag{BG.8}
```

Replace the entries of `B` one at a time by independent standard Gaussians.
The variables match in their first two moments.  Taylor expansion at zero
with integral remainder gives, for one replacement,

```math
|E F_\eta(\ldots,B_e,\ldots)
 -E F_\eta(\ldots,G_e,\ldots)|
\le {E|B_e|^3+E|G_e|^3\over6}\eta^2
=c_{\rm Lin}\eta^2.                               \tag{BG.9}
```

Thus, for a standard Gaussian matrix `G`,

```math
E X(B)\le E X(G)+{(N-1)\log2\over\eta}
                    +c_{\rm Lin}d\eta^2.           \tag{BG.10}
```

It remains to bound the Gaussian maximum.  For

```math
Z_{x,y}=x^TGy,
\qquad
Y_{x,y}=\sqrt n\,g^Tx+\sqrt m\,h^Ty,              \tag{BG.11}
```

where `g,h,G` are independent standard Gaussian arrays, put
`u=x^Tx'` and `v=y^Ty'`.  Direct calculation gives

```math
E(Z_{x,y}-Z_{x',y'})^2
 =2(mn-uv),
```

whereas

```math
E(Y_{x,y}-Y_{x',y'})^2
 =2(2mn-nu-mv).
```

Their difference is

```math
2(m-u)(n-v)\ge0.                                  \tag{BG.12}
```

The finite Sudakov--Fernique comparison therefore yields

```math
E X(G)
\le E\max_{x,y}Y_{x,y}
=\sqrt n\,mE|G|+\sqrt m\,nE|G|
=c_0(m\sqrt n+n\sqrt m).                          \tag{BG.13}
```

Take `eta=N^(-1/3)` in (BG.10).  Equations (BG.3) and (BG.13)
give (BG.4). `square`

The comparison error is `O(N^(4/3))` without hiding its finite-order
origin: `O(N/eta)` is soft-max entropy and `O(mn eta^2)` is Lindeberg
replacement.  Both are `o(N^(3/2))` on comparable splits.

## 2. Entropy transport around the geometric mean

### Lemma BG.2 (Bernoulli transport)

For every law `q` on the bridge cube,

```math
\boxed{
E_qX
\le E_UX+\sqrt{2mnD(q\Vert U)}.}                  \tag{BG.14}
```

*Proof.*  Replacing one bridge bit changes `X` by at most `2`.  The
bounded-differences exponential-moment inequality therefore gives, for
every `theta>=0`,

```math
\log E_Ue^{\theta(X-E_UX)}\le {mn\theta^2\over2}.  \tag{BG.15}
```

Donsker--Varadhan transport gives

```math
E_qX-E_UX
\le {D(q\Vert U)\over\theta}+{mn\theta\over2}.
```

Optimization in `theta`, with the zero-entropy case obtained by continuity,
proves (BG.14). `square`

For the actual negative tilt `dq_s proportional e^(sL)dU`, the exact
single-edge oscillation of the bridge pressure gives

```math
D(q_s\Vert U)\le mn\,\kappa(|s|t),
\qquad
\kappa(a)=a\tanh a-\log\cosh a.                   \tag{BG.16}
```

Combining (BG.4), (BG.14), and (BG.16), uniformly for
`s in [-delta,0]`, gives

```math
\boxed{
E_{q_s}X
\le c_0(m\sqrt n+n\sqrt m)+R_N
       +mn\sqrt{2\kappa(\delta t)}.}              \tag{BG.17}
```

## 3. Uniform actual-child overlap floor

Fix any actual children, either orientation, and let their exact bridge
channel be

```math
L(B)=c+\log E_\mu e^{t\langle B,Q\rangle},
\qquad Q=\tau XY^T,
\qquad t={\beta\over\sqrt N},
\qquad \rho=\tanh t.                              \tag{BG.18}
```

Let `r_e(B_(-e))` be the exact deleted-edge posterior response and put

```math
S_s=E_{q_s}\sum_e r_e^2.                           \tag{BG.19}
```

The one-edge Bayes calculation of Theorem 37.52 is prior independent.  For

```math
A_\rho={1\over1+\rho},
\qquad
C_{\rho,\delta}=A_\rho+{\delta\over1-\rho^2},     \tag{BG.20}
```

it says, for every `s in [-delta,0]`,

```math
{S_s\over mn}
\ge {1\over C_{\rho,\delta}}
 \left\{A_\rho-{E_{q_s}X\over\rho mn}\right\}.  \tag{BG.21}
```

Indeed the full posterior mean `m_e` obeys

```math
m_e={r_e+\rho B_e\over1+\rho B_er_e},
```

the negative tilt has conditional mean
`tanh(s arctanh(rho r_e))`, and summing the resulting conditional lower
bound uses only the pointwise fact

```math
\sum_eB_em_e\le\max_{Q\in\operatorname {supp}\mu}
                       \langle B,Q\rangle\le X(B).
```

Substitution of (BG.17) in (BG.21) proves the finite-order theorem

```math
\boxed{
{S_s\over mn}
\ge {1\over C_{\rho,\delta}}
\left\{
A_\rho
-{c_0(m\sqrt n+n\sqrt m)+R_N\over\rho mn}
-{\sqrt{2\kappa(\delta t)}\over\rho}
\right\}.}                                        \tag{BG.22}
```

No support size or prior weight appears in (BG.22).

### Theorem BG.3 (geometric strong-channel obstruction)

Suppose that along a sequence of comparable splits

```math
{mn\over N^2}\ge\gamma_0>0.                       \tag{BG.23}
```

Define

```math
\beta_{\rm BG}(\gamma_0)
=\sqrt{2\over\pi}
 \sqrt{{1\over\gamma_0}+{2\over\sqrt{\gamma_0}}}. \tag{BG.24}
```

If `beta>beta_BG(gamma_0)`, choose

```math
0<\delta<\min\left\{\lambda,1,
             1-{\beta_{\rm BG}(\gamma_0)\over\beta}\right\}. \tag{BG.25}
```

Then, uniformly over all actual minimizing children and both orientations,

```math
\liminf_{N\to\infty}\inf_{-\delta\le s\le0}
 {S_s\over mn}
\ge
{1-\beta_{\rm BG}(\gamma_0)/\beta-\delta
 \,\over1+\delta}>0.                             \tag{BG.26}
```

Consequently,

```math
\boxed{
\liminf_{N\to\infty}\widehat\rho_N^-(\lambda)
\ge {\delta\over\lambda}
 {1-\beta_{\rm BG}(\gamma_0)/\beta-\delta
  \,\over1+\delta}>0.}                           \tag{BG.27}
```

*Proof.*  If `gamma_N=mn/N^2`, then exactly

```math
\left\{\sqrt{N\over m}+\sqrt{N\over n}\right\}^2
={1\over\gamma_N}+{2\over\sqrt{\gamma_N}}.       \tag{BG.28}
```

The right side decreases with `gamma_N`, so (BG.23) bounds the geometric
term in (BG.22) by `beta_BG(gamma_0)/beta+o(1)`.  Moreover,

```math
{R_N\over\rho mn}=O_{\beta,\gamma_0}(N^{-1/6}),
\qquad
{\sqrt{2\kappa(\delta t)}\over\rho}\longrightarrow\delta,
```

while `A_rho->1` and `C_(rho,delta)->1+delta`.  This proves (BG.26).
Integrating its uniform lower bound over `[-delta,0]` proves (BG.27).
`square`

For a balanced split, `gamma_0=1/4`, and (BG.24) becomes

```math
\boxed{
\beta_{\rm BG}(1/4)={4\over\sqrt\pi}
=2.256758\ldots .}                                 \tag{BG.29}
```

## 4. A lower bound and the ceiling of the maximum-envelope architecture

The Gaussian comparison in (BG.4) is not sharp.  There is nevertheless a
rigorous positive lower bound which shows that replacing it by the exact
bipartite ground-state constant cannot push the present *pointwise-maximum*
argument all the way to an arbitrarily small threshold.

Put

```math
\Psi(a)=E|G+a|
=\sqrt{2\over\pi}e^{-a^2/2}
 +a\{2\Phi(a)-1\},
\qquad a\ge0,                                      \tag{BG.30}
```

where `Phi` is the standard Gaussian distribution function.

### Lemma BG.4 (one alternating step gives a bipartite lower bound)

Suppose `m,n->infinity` and `n/m->r in (0,infinity)`.  Then

```math
\boxed{
\liminf {E_UX(B)\over m\sqrt n}
\ge\Psi\!\left(\sqrt{2\over\pi}\sqrt r\right).}  \tag{BG.31}
```

Interchanging the two shores also gives

```math
\boxed{
\liminf {E_UX(B)\over n\sqrt m}
\ge\Psi\!\left(\sqrt{2\over\pi\,r}\right).}      \tag{BG.32}
```

*Proof.*  Start from the deterministic left spin `x_i^(0)=1`.  For every
column define

```math
S_j=\sum_{i=1}^mB_{ij},
\qquad
y_j=\operatorname {sgn}(S_j),                     \tag{BG.33}
```

with `sgn(0)=1`.  Now perform one maximizing update on the left:

```math
Z_{ij}=B_{ij}y_j,
\qquad
T_i=\sum_{j=1}^nZ_{ij},
\qquad
x_i=\operatorname {sgn}(T_i).                     \tag{BG.34}
```

This is a feasible spin pair, and therefore pointwise

```math
X(B)\ge x^TBy=\sum_{i=1}^m|T_i|.                  \tag{BG.35}
```

For a fixed row `i`, the variables `Z_(ij)` are iid as `j` varies: each is a
function only of the independent `j`th bridge column.  Variables belonging
to different rows in the same column are dependent, but (BG.35) needs no
row independence.  Exchangeability within one column gives their common
mean

```math
\begin{aligned}
\mu_m=E Z_{ij}
 &={1\over m}E\left[S_j\operatorname {sgn}(S_j)\right]\\
 &={E|S_m|\over m}
 ={1\over2^{m-1}}
   {m-1\choose\lfloor(m-1)/2\rfloor},             \tag{BG.36}
\end{aligned}
```

where `S_m` is a sum of `m` fair signs.  In particular, the central binomial
asymptotic gives

```math
\sqrt m\,\mu_m\longrightarrow\sqrt{2\over\pi}.   \tag{BG.37}
```

For each fixed row, `T_i` is thus a sum of `n` iid signs of mean `mu_m`.
The bounded triangular-array central limit theorem and (BG.37) give

```math
{T_i\over\sqrt n}
={\sum_j(Z_{ij}-\mu_m)\over\sqrt n}
 +\sqrt n\,\mu_m
\Longrightarrow
G+\sqrt{2\over\pi}\sqrt r.                       \tag{BG.38}
```

This convergence also holds in first absolute moment.  Indeed,

```math
E\left({T_i\over\sqrt n}\right)^2
=1-\mu_m^2+n\mu_m^2,                              \tag{BG.39}
```

which stays bounded when `n/m` stays bounded.  Hence the variables in
(BG.38) are uniformly integrable.  Taking expectations in (BG.35), using
that all rows have the same marginal law, and applying (BG.30) proves
(BG.31).  Starting on the other shore proves (BG.32). `square`

For balanced shores `m=n=k`, (BG.31) becomes

```math
\liminf_{k\to\infty}{E_UX(B)\over k^{3/2}}
\ge\Psi\!\left(\sqrt{2\over\pi}\right)
=1.039196660145857\ldots .                         \tag{BG.40}
```

The denominator in the pointwise posterior argument is

```math
\rho mn\sim {\beta\over\sqrt{2k}}k^2
={\beta\over\sqrt2}k^{3/2}.                       \tag{BG.41}
```

Consequently, any balanced obstruction proof which seeks a uniform floor on
a tilt interval containing `s=0` and whose only posterior-energy input there
is a universal upper bound on the pointwise envelope `X` must use an upper
constant at least the lower constant in (BG.40).  At `s=0`, such an argument
cannot make its envelope ratio strictly smaller than one unless

```math
\boxed{
\beta>
\sqrt2\,\Psi\!\left(\sqrt{2\over\pi}\right)
=1.469646010751096\ldots .}                        \tag{BG.42}
```

Equation (BG.42) is an **architecture ceiling**, not an actual overlap
threshold.  It neither proves decay below this value nor rules out a
posterior-specific theorem that avoids the maximum `X`, or an argument
localized to negative tilts bounded away from zero.

## 5. Sharp scope

This theorem improves only the universal obstruction threshold.  It does
not prove overlap decay below (BG.24), does not identify the bipartite
spin-glass reconstruction threshold, and does not separate row-additive
cavity response from irreducible reverse-product dependence.  The
`O(N^(4/3))` universality error is negligible at the leading physical scale
but supplies no power-saving recurrence.  Further threshold improvement by
this route requires either a sharper rigorous upper bound on the Bernoulli
bipartite ground-state constant than (BG.4), or a posterior-energy estimate
that uses more than the pointwise maximum.
