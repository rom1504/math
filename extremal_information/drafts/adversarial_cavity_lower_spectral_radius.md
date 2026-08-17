# Adversarial cavity benchmark: lower-spectral thermodynamics

Status: **proved benchmark theorem; discrete-disorder gate passed**.
Portfolio judgment: **keep warm as a semigroup-limit mechanism, but do not
promote to a second full theory**.  The finite-temperature result is not a
Sion relaxation, but its useful cavity compression is already naturally
expressible as a contractive continuous response carrier.  The part which is
genuinely orthogonal--scalar convergence without a finite-period optimizer--
currently fails the dense-bridge transfer test.

## 1. Model

Fix a spin alphabet `[q]`, a finite disorder alphabet `D`, and bounded local
rewards

```math
h_d:[q]\times[q]\longrightarrow\mathbb R,
\qquad d\in D.                                                \tag{AC.1}
```

For a disorder word `w=(d_1,...,d_n)` and a spin path
`sigma=(sigma_0,...,sigma_n)`, put

```math
H_w(\sigma)=\sum_{k=1}^n h_{d_k}(\sigma_{k-1},\sigma_k),

Z_{n,\beta}(w)=\sum_{\sigma\in[q]^{n+1}}
                  \exp\{\beta H_w(\sigma)\}.                \tag{AC.2}
```

The disorder is genuinely adversarial and discrete:

```math
F_{n,\beta}=\min_{w\in D^n}\log Z_{n,\beta}(w),
\qquad
G_n=\min_{w\in D^n}\max_{\sigma}H_w(\sigma).                \tag{AC.3}
```

No probability distribution on `D` and no convex hull of its transfer
matrices is introduced.

## 2. Exact theorem

### Theorem AC (discrete adversarial thermodynamic limit)

Let

```math
T_d(i,j)=e^{\beta h_d(j,i)}>0,
\qquad \mathcal T_\beta=\{T_d:d\in D\}.                      \tag{AC.4}
```

The harmless transpose convention makes column vectors evolve forward.
For every fixed `beta>0`,

```math
p(\beta):=\lim_{n\to\infty}{F_{n,\beta}\over n}
          =\log\check\rho(\mathcal T_\beta),                \tag{AC.5}
```

where the lower spectral radius is

```math
\check\rho(\mathcal T)
=\lim_{n\to\infty}
 \min_{d_1,\ldots,d_n}
 \|T_{d_n}\cdots T_{d_1}\|_1^{1/n}.                         \tag{AC.6}
```

Moreover the adversarial ground-state density exists and is recovered at
zero temperature:

```math
g:=\lim_{n\to\infty}{G_n\over n}
  =\lim_{\beta\to\infty}{p(\beta)\over\beta},

0\le {p(\beta)\over\beta}-g\le {\log q\over\beta}.          \tag{AC.7}
```

Thus this is a genuine minimax thermodynamic-limit theorem with a finite,
nonconvex disorder alphabet.

### Proof

For a nonnegative matrix `P`, with `||.||_1` the maximum column-sum norm,

```math
\|P\|_1\le \mathbf1^TP\mathbf1\le q\|P\|_1.                 \tag{AC.8}
```

If

```math
b_n=\log\min_{d_1,\ldots,d_n}
              \|T_{d_n}\cdots T_{d_1}\|_1,
```

then matrix-norm submultiplicativity gives

```math
b_{m+n}\le b_m+b_n.                                          \tag{AC.9}
```

Fekete's lemma and (AC.8) prove (AC.5).  This is the standard
lower-spectral-radius limit, reviewed with the norm-independent formula in
[Guglielmi--Zennaro, Section 2](https://arts.units.it/retrieve/handle/11368/2972807/338571/GZ2019_finale.pdf).

For every word,

```math
\max_\sigma H_w(\sigma)
\le {1\over\beta}\log Z_{n,\beta}(w)
\le \max_\sigma H_w(\sigma)+{(n+1)\log q\over\beta}.        \tag{AC.10}
```

Taking the minimum over the same discrete word set gives

```math
{G_n\over n}
\le {F_{n,\beta}\over\beta n}
\le {G_n\over n}+\left(1+{1\over n}\right){\log q\over\beta}.
                                                                    \tag{AC.11}
```

After `n -> infinity`, the limsup--liminf gap of `G_n/n` is at most
`log(q)/beta` for every `beta`; sending `beta -> infinity` proves convergence.
Taking the resulting `n`-limit in (AC.11), and then `beta -> infinity`, proves
(AC.7). `square`

## 3. The exact cavity/response state in the positive regime

The transfer theorem has a genuine cavity formulation.  On the probability
simplex let

```math
\tau_d(p)={T_dp\over\|T_dp\|_1},
\qquad
r_d(p)=\log\|T_dp\|_1.                                      \tag{AC.12}
```

Then the log partition function telescopes as an accumulated reward along
the projective orbit.  The finite-horizon adversarial cavity operator is

```math
(\mathcal Vf)(p)=\min_{d\in D}
       \{r_d(p)+f(\tau_d(p))\}.                              \tag{AC.13}
```

Suppose uniformly that

```math
0<a\le T_d(i,j)\le b<\infty.                                \tag{AC.14}
```

After one step every cavity vector lies in the common invariant compact set

```math
X=\left\{p\in\Delta_{q-1}:
 {a\over qb}\le p_i\le {b\over qa}\text{ for every }i\right\}.       \tag{AC.14a}
```

Indeed each coordinate of `T_dp` lies in `[a,b]` and its coordinate sum in
`[qa,qb]`, even when the input lies on the boundary. Thus `X` is contained in
the simplex interior.
In Hilbert's projective metric, Birkhoff's contraction theorem gives

```math
d_H(\tau_d(p),\tau_d(p'))\le\kappa d_H(p,p'),
\qquad
\kappa\le\tanh\!\left({1\over2}\log{b\over a}\right)
       ={b-a\over b+a}<1.                                   \tag{AC.15}
```

The source theorem is
[Birkhoff, *Extensions of Jentzsch's theorem*](https://doi.org/10.1090/S0002-9947-1957-0087058-6);
the explicit `tanh(Delta/4)` coefficient is also stated in
[Carroll](https://doi.org/10.1016/j.laa.2004.02.039).
Order comparison gives

```math
|r_d(p)-r_d(p')|\le d_H(p,p').                              \tag{AC.16}
```

Consequently there are a scalar `lambda` and a Lipschitz function `u` on
`X` such that

```math
\mathcal Vu=u+\lambda,
\qquad
\operatorname{Lip}_{d_H}(u)\le {1\over1-\kappa},
\qquad
\lambda=\log\check\rho(\mathcal T_\beta).                  \tag{AC.17}
```

Here is a short proof, included to avoid importing an inapplicable stochastic
cavity theorem.  Normalize functions by `f(p_*)=0`.  Equations (AC.15)--
(AC.16) imply

```math
\operatorname{Lip}(\mathcal Vf)
\le1+\kappa\operatorname{Lip}(f).                            \tag{AC.18}
```

Thus `f -> V f-(V f)(p_*)` maps the compact convex set of normalized
`1/(1-kappa)`-Lipschitz functions into itself.  Schauder gives a fixed point.
Since `V` is sup-norm nonexpansive,

```math
\|\mathcal V^n0-(u+n\lambda)\|_\infty\le\|u\|_\infty.       \tag{AC.19}
```

Finally, for any fixed interior `p`, `||Pp||_1` and `||P||_1` differ by a
factor bounded independently of the product length, identifying `lambda`
with (AC.6).

This gives a quantitative approximate carrier as well.  A Hilbert-metric
`delta`-net of `X` has

```math
O_{q,a,b}(\delta^{-(q-1)})                                   \tag{AC.20}
```

states.  Round every `tau_d(p)` to a nearest net point `p'` and retain the
exact reward at `p`. Put `L=1/(1-kappa)`. The additive eigen-equation and
Lipschitz bound give, for every rounded edge,

```math
r_d(p)+u(p')\ge\lambda+u(p)-L\delta.                       \tag{AC.20a}
```

Every directed cycle therefore has mean at least `lambda-L delta`. At each
net point choose a symbol attaining the minimum in (AC.17). Along that
selected edge the reverse estimate holds,

```math
r_d(p)+u(p')\le\lambda+u(p)+L\delta.                       \tag{AC.20b}
```

and following selected edges eventually enters a cycle of mean at most
`lambda+L delta`. Hence the minimum cycle mean `lambda_delta` of the rounded
finite controller satisfies

```math
|\lambda_\delta-\lambda|
\le {\delta\over1-\kappa}.                                  \tag{AC.21}
```

This is an **asymptotic mean-pressure** response-rate conclusion, not just
existence of pressure: fixed spin width has polynomial approximate state
complexity, with upper exponent given by the projective cavity dimension
`q-1`. It does not claim depth-independent additive error for the
unnormalized finite-horizon response; fresh rounding may accumulate
`O(n delta/(1-kappa))` total error while retaining the displayed per-step
bound. No matching lower bound is claimed without an exposed-response
hypothesis.

## 4. Why the scalar limit is not merely a periodic finite-state ansatz

The lower-spectral mechanism survives even when projective contraction is
lost.  Consider the two nonnegative transfer matrices

```math
A=\begin{pmatrix}1/3&0\\0&3\end{pmatrix},
\qquad
B=\begin{pmatrix}2&0\\0&1/2\end{pmatrix}.                   \tag{AC.22}
```

A length-`n` word containing `k` copies of `A` has

```math
{1\over n}\log\rho(P)
=\left|{n-k\over n}\log2-{k\over n}\log3\right|.           \tag{AC.23}
```

Rational frequencies approximate

```math
{k\over n}\longrightarrow{\log2\over\log6},                \tag{AC.24}
```

so `check rho({A,B})=1`.  But no finite word attains the lower spectral
radius: equality would require `3^k=2^(n-k)`.  This elementary calculation is
the lower-radius analogue of the broader failure of the lower finiteness
property studied by
[Bochi--Morris](https://arxiv.org/abs/1309.0319), and of aperiodic
matrix-product optimization phenomena in
[Bousch--Mairesse](https://www.imo.universite-paris-saclay.fr/~thierry.bousch/preprints/artetris.pdf).
It is not attributed to a numbered example in either source; equations
(AC.23)--(AC.24) are the complete proof needed here.

The matrices in (AC.22) encode a two-sector chain with forbidden sector
switches, so they are a hard-constraint limit rather than a uniformly
positive finite-temperature model.  The point is appropriately narrow:
scalar adversarial pressure can converge by semigroup subadditivity even
when no bounded-period disorder pattern certifies the limit.

## 5. Decisive transfer falsifier for dense sign disorder

The chain theorem works because gluing two words changes only a bounded
interface, or equivalently multiplies two fixed-dimensional transfer
operators.  Neither fact survives a balanced split of a dense sign
quadratic.

Indeed, for every `n by n` sign matrix `B`,

```math
\max_{x,y\in\{\pm1\}^n}x^TBy\ge {n^{3/2}\over\sqrt3}.        \tag{AC.25}
```

To prove this, choose `x` uniformly.  For each column,
`S_j=(B^Tx)_j` has `E S_j^2=n` and `E S_j^4<=3n^2`.  Interpolation between
`L^1,L^2,L^4` gives `E|S_j|>=sqrt(n/3)`.  Hence some `x` has
`sum_j|S_j|>=n^(3/2)/sqrt(3)`, and choosing
`y_j=sign(S_j)` proves (AC.25).

Thus the cross-interface term is already on the leading `n^(3/2)` scale.
In transfer language, a path decomposition of the complete interaction has
separator width `n`, transfer dimension `2^n`, and local Boltzmann ratios
which grow exponentially with `n`; the Birkhoff coefficient tends to one.
There is no fixed-dimensional product and no subleading boundary correction
to which (AC.9) could be applied.

This is the requested falsifier.  It does **not** say that a hidden
multiplicative norm for dense sign quadratics is impossible.  It says that
the ordinary local transfer/lower-spectral construction does not provide
one: such a norm would need an additional algebraic quotient which removes
a provably leading bridge response.

## 6. Director judgment

This benchmark materially improves the adversarial-statistical-mechanics
card:

1. the discrete disorder constraint is retained exactly;
2. fixed-temperature and ground-state limits both exist;
3. the response state and its approximate complexity are explicit; and
4. the proof never uses Sion or randomized disorder.

It nevertheless does not warrant a new full theory branch.  In the uniformly
positive regime, `(projective cavity, scalar potential)` is precisely a
continuous presented carrier with contraction, so the current contextual-
response theory predicts (AC.20)--(AC.21).  Outside that regime, lower
spectral radius is a genuinely orthogonal scalar semigroup invariant, but
(AC.25) blocks its naive dense transfer.

**Recommendation:** import Theorem AC as a strong adversarial benchmark and
keep lower-spectral/antinorm methods warm.  Promote them only if one finds a
nonlocal product representation for a dense or growing-interface model whose
operator dimension is sub-landscape and whose multiplication defect is
subleading.  The falsifier is failure of either condition.
