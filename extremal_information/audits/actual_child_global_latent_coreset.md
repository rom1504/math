# A global latent coreset for the actual inverse-escort channel

**Status.**  Rigorous integrable coreset theorem and scope audit.  One
global empirical sample from the zero-bridge latent child law defines one
genuine rank-one bridge likelihood.  Its cavity field is therefore exactly
curl-free.  The expected physical cavity error is controlled by the
inverse-escort average of one deleted-posterior chi-square complexity:

```math
t^2E_q\sum_e(r_e-r_e^{(R)})^2
\le 32t^2mn\,{\overline K_{\rm del}\over R}.       \tag{GC.1}
```

The complete likelihood has an exact relative-Monte-Carlo identity, and a
small-ball argument converts it to scalar log-likelihood error.  In
particular,

```math
\log\overline K_{\rm del}=o(N)                    \tag{GC.2}
```

would give a single subexponential-support, integrable coreset with
power-saving physical cavity error and `o(N)` inverse-escort pressure error.

No current actual-child theorem proves (GC.2).  Existing bounded row
Renyi-two complexity concerns the output law `q` relative to a row product;
`overline K_del` is instead the annealed Renyi-two complexity of the
**latent Bayes posterior relative to its child prior**.  The two quantities
have no applicable data-processing comparison.  Thus this theorem sharpens
the generation/integrability problem but does not close it or create a
Level-6 recurrence.

## 1. Actual channel and two posterior complexities

Put

```math
d=mn,qquad N=m+n,qquad t={\beta\over\sqrt N},qquad
\rho=\tanh t.                                     \tag{GC.3}
```

Let `mu` be the exact zero-bridge latent law of the two children on rank-one
sign words `Q`.  The normalized channel kernel and output likelihood are

```math
k_Q(B)=\prod_{f=1}^{d}(1+\rho B_fQ_f),
\qquad P(B)=E_\mu k_Q(B).                          \tag{GC.4}
```

Let the actual inverse escort be

```math
{dq_\lambda\over dU}(B)={P(B)^{-\lambda}\over Z_\lambda}.
                                                               \tag{GC.5}
```

For a deleted edge `e`, define

```math
k_{e,Q}(B_{-e})=\prod_{f\ne e}(1+\rho B_fQ_f),
\qquad D_e(B_{-e})=E_\mu k_{e,Q}(B_{-e}),          \tag{GC.6}
```

and

```math
r_e(B_{-e})={E_\mu Q_ek_{e,Q}(B_{-e})\over D_e(B_{-e})}.
                                                               \tag{GC.7}
```

The deleted posterior is

```math
{d\mu_{e,B}\over d\mu}(Q)={k_{e,Q}(B_{-e})\over D_e(B_{-e})}.
                                                               \tag{GC.8}
```

Its collision complexity and its inverse-escort average are

```math
K_e(B_{-e})
=E_\mu\left({k_{e,Q}\over D_e}\right)^2
=1+\chi^2(\mu_{e,B}\Vert\mu)
=\exp D_2(\mu_{e,B}\Vert\mu),                    \tag{GC.9}
```

and

```math
\boxed{
\overline K_{\rm del}
={1\over d}\sum_{e=1}^{d}E_{q_\lambda}K_e(B_{-e}).}           \tag{GC.10}
```

Similarly, for the complete posterior

```math
{d\mu_B\over d\mu}(Q)={k_Q(B)\over P(B)},
```

put

```math
K_0(B)=E_\mu\left({k_Q(B)\over P(B)}\right)^2
=1+\chi^2(\mu_B\Vert\mu),
\qquad
\overline K_{\rm full}=E_{q_\lambda}K_0(B).       \tag{GC.11}
```

These are annealed collision factors, not merely averages of their
logarithms.

### Lemma GC.0 (exact doubled-channel identity)

Keep the latent child law `mu` fixed, and write

```math
P_u(B)=E_\mu\prod_{f=1}^{d}(1+\tanh(u)B_fQ_f),
```

with `P_{u,e}(B_{-e})` denoting the product with edge `e` deleted.  Then

```math
\boxed{
K_0(B)=(1+\rho^2)^d{P_{2t}(B)\over P_t(B)^2},
\qquad
K_e(B_{-e})=(1+\rho^2)^{d-1}
 {P_{2t,e}(B_{-e})\over P_{t,e}(B_{-e})^2}.}      \tag{GC.11a}
```

*Proof.*  For `z in {+1,-1}`,

```math
(1+\rho z)^2
=(1+\rho^2)\left(1+{2\rho\over1+\rho^2}z\right)
=(1+\rho^2)(1+\tanh(2t)z).
```

Multiplying this identity over all coordinates and averaging against the
same latent law `mu` gives the full formula; deleting one coordinate gives
the second.  `square`

Thus `overline K_full` and `overline K_del` are inverse-escort averages of
explicit scalar two-temperature pressure ratios.  They do **not** conceal
a request for the complete cavity-response table.  Notice, however, that
`d log(1+rho^2)=Theta_beta(N)` at comparable splits.  Establishing
subexponential collision complexity therefore requires a genuine
order-`N` cancellation between this universal prefactor and
`P_{2t}/P_t^2`; it does not follow from a pointwise channel bound.

## 2. One global empirical latent law

Draw, once and for all,

```math
Q^1,\ldots,Q^R\stackrel{\rm iid}{\sim}\mu,
\qquad
\mu_R={1\over R}\sum_{a=1}^{R}\delta_{Q^a}.       \tag{GC.12}
```

This defines the single global likelihood

```math
P_R(B)={1\over R}\sum_{a=1}^{R}k_{Q^a}(B),        \tag{GC.13}
```

and its deleted cavity

```math
r_e^{(R)}(B_{-e})
={\sum_aQ_e^ak_{e,Q^a}(B_{-e})
  \over\sum_ak_{e,Q^a}(B_{-e})}.                  \tag{GC.14}
```

Every denominator is positive.  If exact central symmetry is desired,
replace each sampled word by the pair `{Q^a,-Q^a}`.  This doubles the
support and Rao--Blackwellizes the estimator over the antipodal orbit, so
none of the bounds below worsens.

## 3. Exact complete-likelihood sample complexity

### Theorem GC.1 (relative-MSE identity)

For the random coreset (GC.12),

```math
\boxed{
E_{\mu_R}E_{q_\lambda}
 \left({P_R(B)\over P(B)}-1\right)^2
={\overline K_{\rm full}-1\over R}.}              \tag{GC.15}
```

Thus the iid empirical scheme has exact mean-square sample complexity
`(overline K_full-1)/epsilon^2` for relative likelihood error `epsilon`.

*Proof.*  At fixed `B`, let

```math
W_a(B)={k_{Q^a}(B)\over P(B)}.
```

Then the `W_a` are iid with mean one and second moment `K_0(B)`, while
`P_R/P=R^(-1)sum_aW_a`.  Conditional variance gives

```math
E_{\mu_R}(P_R/P-1)^2={K_0(B)-1\over R}.
```

Average against `q_lambda`.  `square`

## 4. Simultaneous cavity approximation

### Theorem GC.2 (one coreset controls all deleted cavities)

For every edge and deleted bridge word,

```math
\boxed{
E_{\mu_R}(r_e^{(R)}-r_e)^2
\le {32K_e(B_{-e})\over R}.}                      \tag{GC.16}
```

Consequently one deterministic realization of the single coreset exists
for which

```math
\boxed{
t^2E_{q_\lambda}\sum_{e=1}^{d}
 (r_e^{(R)}-r_e)^2
\le32t^2d\,{\overline K_{\rm del}\over R}.}       \tag{GC.17}
```

At comparable splits this is

```math
O_\beta\left(N{\overline K_{\rm del}\over R}\right).          \tag{GC.18}
```

*Proof.*  Fix `(e,B_(-e))` and abbreviate

```math
W_a={k_{e,Q^a}\over D_e},
\qquad Z_a=W_a(Q_e^a-r_e).
```

Then

```math
E W_a=1,qquad E Z_a=0,qquad
E Z_a^2\le4E W_a^2=4K_e.                          \tag{GC.19}
```

Writing bars for sample means gives

```math
r_e^{(R)}-r_e={\overline Z\over\overline W}.
```

On `|bar W-1|<=1/2`, its square is at most `4 bar Z^2`, whose expectation
is at most `16K_e/R`.  Off that event, both responses lie in `[-1,1]`, so
the squared error is at most four, while Chebyshev gives

```math
\Pr\{|\overline W-1|>1/2\}
\le{4(K_e-1)\over R}.                             \tag{GC.20}
```

The off-event contribution is at most `16(K_e-1)/R`; this proves (GC.16).
Sum, average, and use the probabilistic method to obtain one common sample
realization satisfying (GC.17).  `square`

No union bound over edges or bridge words is used.  The guarantee is in the
declared `q_lambda`-weighted aggregate norm, which is the norm entering the
physical overlap and recurrence audits.

## 5. Full and deleted posterior complexity are equivalent at this scale

### Lemma GC.3 (one-edge comparison)

Pointwise, for every edge,

```math
\boxed{e^{-4t}K_e(B_{-e})\le K_0(B)\le e^{4t}K_e(B_{-e}),
\qquad
\overline K_{\rm full}\le e^{4t}\overline K_{\rm del}.}       \tag{GC.21}
```

*Proof.*  The insertion identities give

```math
{k_Q(B)/P(B)\over k_{e,Q}(B_{-e})/D_e(B_{-e})}
={1+\rho B_eQ_e\over1+\rho B_er_e}.              \tag{GC.22}
```

Both numerator and denominator lie in `[1-rho,1+rho]`; hence the ratio is
at most

```math
{1+\rho\over1-\rho}=e^{2t}.
```

The same ratio is at least `e^(-2t)`.  Square and average over `mu` to get
both pointwise assertions.  The upper assertion holds for every `e`; average
it over edges and then over `q_lambda`.  `square`

Thus the single scalar `overline K_del` controls both the complete
likelihood and all deleted cavities.

## 6. Integrability and zero holonomy

Define the half-flip discrete derivative

```math
\nabla_eF(B_{-e})
={F(B_e=+,B_{-e})-F(B_e=-,B_{-e})\over2}.          \tag{GC.23}
```

### Corollary GC.4 (the coreset cavity field is an exact gradient)

For every empirical coreset,

```math
\boxed{
\nabla_e\log P_R
=\operatorname {arctanh}(\rho r_e^{(R)}).}        \tag{GC.24}
```

In particular, every square circulation vanishes exactly.  Moreover,

```math
\boxed{
E_{\mu_R}E_q\sum_e
 |\nabla_e\log P_R-\nabla_e\log P|^2
\le {32\rho^2\over(1-\rho^2)^2}
      d{\overline K_{\rm del}\over R}.}          \tag{GC.25}
```

*Proof.*  Deleted insertion factors `P_R` as

```math
P_R(B)=D_e^{(R)}(B_{-e})
       \{1+\rho B_er_e^{(R)}(B_{-e})\},
```

which proves (GC.24).  The function
`x -> arctanh(rho x)` is `rho/(1-rho^2)`-Lipschitz on `[-1,1]`; apply
Theorem GC.2.  `square`

This resolves the Hodge/path-realization issue for the coreset: unlike an
adaptive pointwise SVD, the approximate vector field comes from one scalar
potential by construction.

## 7. From relative MSE to scalar log-likelihood error

Relative MSE alone does not automatically control a logarithm, because an
empirical mean can have a rare severe lower tail.  The channel has a finite
likelihood range, and a Paley--Zygmund argument gives the needed explicit
correction.

### Lemma GC.5 (positive sample-mean logarithm bound)

Let `W_1,...,W_R` be iid positive variables with

```math
EW=1,qquad EW^2=K,qquad W\ge e^{-L}.
```

Then

```math
\boxed{
E\left|\log{1\over R}\sum_{a=1}^{R}W_a\right|
\le3\sqrt{K/R}
 +{4K\log(2R)\over R}
 +L e^{-R/(4K)}.}                                 \tag{GC.26}
```

*Proof.*  Put `bar W=R^(-1)sum W_a`.  On
`|bar W-1|<=1/2`, use `|log bar W|<=2|bar W-1|` and variance.  The positive
tail is bounded by `log x<=x-1` and the same first absolute moment.  On the
remaining lower tail, Chebyshev gives probability at most `4(K-1)/R`.

Paley--Zygmund gives

```math
\Pr\{W\ge1/2\}\ge{1\over4K}.                     \tag{GC.27}
```

If at least one sample has `W_a>=1/2`, then `bar W>=1/(2R)` and the lower
logarithm is at most `log(2R)`.  The probability that no sample does is at
most `exp{-R/(4K)}`; on that event the deterministic lower bound gives
`-log bar W<=L`.  Combining the pieces proves (GC.26), with harmlessly
enlarged constants.  `square`

For `W=k_Q(B)/P(B)`, every kernel and its mixture lie between
`(1-rho)^d` and `(1+rho)^d`.  Hence

```math
W\ge\left({1-\rho\over1+\rho}\right)^d=e^{-2td}.                \tag{GC.28}
```

Average (GC.26) over `B`.  For any deterministic threshold
`H>=overline K_full`, split according to `K_0(B)<=H` and use Markov on its
complement.  Jensen handles the first two terms.  This proves:

### Theorem GC.6 (averaged scalar-potential bound)

For every `H>0`,

```math
\boxed{
\begin{aligned}
E_{\mu_R}E_{q_\lambda}|\log P_R-\log P|
\le{}&3\sqrt{\overline K_{\rm full}/R}
 +{4\overline K_{\rm full}\log(2R)\over R}\\
&+2td\left{e^{-R/(4H)}
       +{\overline K_{\rm full}\over H}\right\}.
                                                               \tag{GC.29}
\end{aligned}}
```

Again, some single deterministic coreset achieves the displayed averaged
bound.

When both the cavity and scalar-potential conclusions are needed, apply the
probabilistic method to the sum of the two errors normalized by their
expectations.  One realization satisfies both bounds with at most a common
factor two, which is absorbed in the asymptotic constants below.

The threshold term is not cosmetic.  A scalar random variable can have
mean one and bounded second moment while taking an arbitrarily small value
with fixed positive probability; with all `R` samples in that atom, the
sample mean has an arbitrarily large negative logarithm.  Posterior
chi-square controls cavity ratios directly because they remain in
`[-1,1]`; scalar pressure additionally needs the finite channel range or a
stronger lower-tail hypothesis.

## 8. Conditional subexponential-support theorem

### Corollary GC.7 (annealed subexponential posterior complexity is enough)

Assume comparable splits and

```math
\boxed{\log\overline K_{\rm del}=o(N).}           \tag{GC.30}
```

Fix any `zeta in(0,1/2)`, put

```math
\mathcal K=e^{4t}\overline K_{\rm del},
\qquad H=\mathcal K N^{1/2+\zeta},
\qquad R=\left\lceil16H(\log N)^2\right\rceil.   \tag{GC.31}
```

Then `R=exp{o(N)}` and there is one global `R`-atom rank-one channel such
that

```math
\boxed{
t^2E_{q_\lambda}\sum_e(r_e-r_e^{(R)})^2
=O_\beta\left({N^{1/2-\zeta}\over(\log N)^2}\right),}          \tag{GC.32}
```

the same power bound holds for the exact-gradient error up to a fixed
`beta`-dependent constant, and

```math
\boxed{E_{q_\lambda}|\log P_R-log P|=o(N).}       \tag{GC.33}
```

More explicitly, the lower-tail contribution in (GC.29) is
`O_beta(N^(1-zeta))+exp{-Omega((log N)^2)}`, and every other term in
(GC.29) is `o(N)` under (GC.30).

*Proof.*  Combine (GC.17), (GC.21), and (GC.31) to obtain (GC.32).  In
(GC.29),

```math
{\overline K_{\rm full}\over H}\le N^{-1/2-\zeta},
\qquad {R\over4H}\ge4(\log N)^2,
```

while `td=Theta_beta(N^(3/2))`.  This gives the stated lower-tail bound.
The square-root and sample-variance terms are `o(N)` because
`log R=o(N)`.  `square`

This is a genuine global integrability theorem: it replaces the full latent
mixture by one finite mixture before evaluating any bridge word, and all
future cavity queries use the same scalar likelihood.

It is not yet a low-information quotient in the strongest project sense.
Each rank-one atom contains `N-2` projective sign bits, so a raw description
uses `Theta(RN)` bits.  Condition (GC.30) makes the **support size**
subexponential, but does not make the response image `exp{o(N)}` or the
description length `o(N)`.  A recurrence would still need a mergeable or
symmetry-compressed presentation and directional target relevance.

## 9. Collision with current actual-child information

Theorem 37.18 bounds the Renyi-two complexity of each conditional **bridge
row under `q_lambda` relative to the fair row law**.  Theorems 37.25--37.29
likewise control row filtrations or latent decompositions of the output
escort.  In contrast, (GC.9) measures how sharply observing a bridge word
concentrates the **child latent law** relative to `mu`.  Bayes inversion
reverses the channel direction, so the existing row bounds do not upper
bound (GC.9).

The elementary actual Gibbs bounds give only

```math
K_e(B_{-e})\le{1\over\min_{Q\in\operatorname {supp}\mu}\mu(Q)}
=\exp\{O_\beta(N)\}.                              \tag{GC.34}
```

This follows from `sum_Q mu_(e,B)(Q)^2/mu(Q)<=1/min mu` and the actual-child
cap bound on the smallest Gibbs atom.  It is exactly one exponential order
too weak for (GC.30).  Macroscopic conditional spread bounds the largest
prior atom and narrow coordinate basins; it does not bound the smallest
atom or the inverse-escort posterior collision.  Replica-overlap
regularization controls Gram laws, not the labelled likelihood ratios in
(GC.9).

Therefore the exact new conditional target is:

> prove `log overline K_del=o(N)`, equivalently control the inverse-escort
> average of the deleted scalar two-temperature ratio in (GC.11a), or
> prove that
> `overline K_del>=exp{cN}` because an inverse-escort-positive mass of bridge
> words localizes the latent child phase.

Either outcome is structural.  The first gives the global integrable
coreset above; the second proves that prior-sampled latent coresets cannot
be the missing frame-synchronization mechanism.  At present neither branch
is decided, so the theorem narrows but does not discharge
`L_actual-posterior-frame-synchronization`.

The scalar condition is strictly lower-information than retaining the full
cavity table, by GC.0, but is not yet known to be mathematically easier on
the actual optimizing-child class.
