# Generic Ghirlanda--Guerra perturbations: exact actual-child mapping and bridge boundary

Status: **primary-source mapping audit**.  There is a genuine theorem for an
arbitrary deterministic child Hamiltonian: a random mixed-`p`-spin
perturbation with `o(N)` pressure cost forces the Ghirlanda--Guerra identities,
and a multi-species version forces overlap synchronization.  What it compresses
is the **replica Gram law**.  No cited theorem turns that law into a quotient
sufficient for the negative bridge escort, the edge-cavity responses, or
`J-I^leftarrow`.  The last conversion remains a separate missing theorem.

Primary sources checked at theorem/proof level:

- D. Panchenko, [*Introduction to the SK
  model*](https://arxiv.org/abs/1412.0170), especially the perturbation and
  quantitative Ghirlanda--Guerra estimate in Section 10;
- D. Panchenko, [*The Parisi ultrametricity
  conjecture*](https://arxiv.org/abs/1112.1003);
- D. Panchenko, [*The free energy in a multi-species Sherrington--Kirkpatrick
  model*](https://arxiv.org/abs/1310.6679), especially Theorems 2--4.

The normalization below is transcribed from these sources and then specialized
to deterministic quadratic children.

## 1. A deterministic child satisfies the perturbation theorem with a better scale

Let `k` be the child order, let

```math
H_0(\sigma)=tH_A(\sigma),\qquad
H_A(\sigma)=\sum_{i<j}a_{ij}\sigma_i\sigma_j,qquad
t={\beta\over\sqrt N},                                    \tag{GG.1}
```

where `A` is arbitrary and deterministic.  Here `k` and the parent order `N`
are comparable in the intended application; nothing in this section uses
minimality of `A`.

For each `p>=1`, define the Gaussian pure-`p` field

```math
g_{p,k}(\sigma)=k^{-p/2}
 \sum_{i_1,\ldots,i_p\le k}g'_{i_1\cdots i_p}
 \sigma_{i_1}\cdots\sigma_{i_p}.                          \tag{GG.2}
```

Thus, for `R(\sigma,\tau)=k^{-1}\sum_i\sigma_i\tau_i`, exactly

```math
\mathbb E g_{p,k}(\sigma)g_{p,k}(\tau)=R(\sigma,\tau)^p.   \tag{GG.3}
```

For coefficients `x_p in [1,2]`, put

```math
g_x(\sigma)=\sum_{p\ge1}2^{-p}x_pg_{p,k}(\sigma),
\qquad H_s(\sigma)=H_0(\sigma)+s_kg_x(\sigma).             \tag{GG.4}
```

The variance constant

```math
C_x:=\sup_\sigma\mathbb E g_x(\sigma)^2
 =\sum_{p\ge1}4^{-p}x_p^2\le {4\over3}.                  \tag{GG.5}
```

is independent of `k`, `A`, and `t`.  Jensen's inequality gives the exact
expected-pressure comparison

```math
0\le
\mathbb E_g\log Z_s-\log Z_0
\le {C_xs_k^2\over2}.                                    \tag{GG.6}
```

In particular, `s_k^2/k -> 0` makes the perturbation `o(k)` in pressure.

The important optimizer-specific simplification is that `H_0` is
deterministic.  Panchenko's concentration envelope takes a supremum over
`0<=x_p<=3`, because the coefficient interpolation used in the proof can
shift an `[1,2]` coefficient by a number in `[0,1]`.  On this enlarged cube
the variance bound is `3`, rather than the `4/3` in (GG.5).  As a function
of all Gaussian coefficients in (GG.2), `log Z_s` is therefore at most
`sqrt(3)s_k`-Lipschitz.  Gaussian concentration gives, uniformly over
**every** deterministic signing `A`,

```math
v_k(s_k):=\sup_x\mathbb E_g
 |\log Z_s-\mathbb E_g\log Z_s|
 \le C s_k.                                               \tag{GG.7}
```

Here and only here the supremum is over `[0,3]^N`; the actual perturbation
in (GG.4) still has `x_p in [1,2]`, so (GG.5)--(GG.6) retain the sharper
constant `4/3`.

The quantitative GG theorem in the source assumes the first two conditions
below; the third is the additional condition making the perturbation
negligible in normalized pressure:

```math
s_k\longrightarrow\infty,
\qquad {v_k(s_k)\over s_k^2}\longrightarrow0,
\qquad {s_k^2\over k}\longrightarrow0.                  \tag{GG.8}
```

Consequently, for a deterministic child one may take

```math
s_k=k^\gamma\quad\hbox{for any }0<\gamma<1/2.             \tag{GG.9}
```

The familiar restriction `1/4<gamma<1/2` in the SK example is caused by
the `Theta(sqrt(k))` concentration scale of the *random base Hamiltonian*;
it is not needed for a fixed child.

## 2. Exact theorem and quantitative error

Let `G_{k,s}` be the Gibbs measure for (GG.4), let
`R_{\ell,\ell'}` be overlaps of i.i.d. replicas, and for `|f|<=1` define

```math
\begin{split}
\Delta_{k}(f,n,p)=\Big|&\mathbb E\langle
 fR_{1,n+1}^p\rangle
-{1\over n}\mathbb E\langle f\rangle
 \mathbb E\langle R_{1,2}^p\rangle\\
&-{1\over n}\sum_{\ell=2}^n
 \mathbb E\langle fR_{1,\ell}^p\rangle\Big|.
                                                               \tag{GG.10}
\end{split}
```

The expectations include the Gaussian perturbation.  Panchenko's proof gives,
after averaging the auxiliary `x_p` uniformly on `[1,2]`,

```math
\mathbb E_x\Delta_k(f,n,p)
\le {2^p\over n}
 \left({2\over s_k}+48{\sqrt{v_k(s_k)}\over s_k}\right),     \tag{GG.11}
```

provided `v_k(s_k)/s_k^2<=4^{-p}`.  Combining (GG.7) and (GG.9), for each
fixed `p,n`,

```math
\mathbb E_x\Delta_k(f,n,p)=O_{p,n}(k^{-\gamma/2}).             \tag{GG.12}
```

A diagonal choice `x^k` over a countable convergence-determining algebra of
overlap tests makes every fixed test converge along the sequence.  After
polynomial approximation and a monotone-class extension, every
subsequential limiting annealed replica-overlap array satisfies the exact
identities.  Panchenko's ultrametricity theorem implies

```math
R_{2,3}\ge\min(R_{1,2},R_{1,3})\quad\hbox{almost surely}.      \tag{GG.13}
```

Moreover, the exact Ghirlanda--Guerra identities determine the distribution
of the whole limiting overlap array from the one-overlap law `zeta`.

This is a real low-information statement: at fixed overlap resolution
`epsilon`, `zeta` has only `O(1/epsilon)` scalar bins.  It is, however, a
statement about the **law of replica Gram matrices**, not a finite quotient of
physical configurations that answers external bridge queries.

## 3. The perturbation is uniformly sublinear even at the bridge-response level

The process maximum obeys the elementary union bound

```math
\mathbb P\left(\|g_x\|_\infty>u\right)
\le 2^{k+1}\exp\left(-{u^2\over2C_x}\right).                 \tag{GG.14}
```

Hence `||s_kg_x||_infty=O_P(s_k sqrt(k))`.  For two comparable children,
let `L(B)` be the **normalized child-Gibbs response** in (GG.18) and let
`L_s(B)` be that response after separately perturbing the two child factors.
A direct Radon--Nikodym bound gives

```math
\sup_B|L_s(B)-L(B)|
\le \operatorname{osc}(s_mg_x)
   +\operatorname{osc}(s_ng_y)
=O_P(N^{1/2+\gamma}).                                      \tag{GG.15}
```

Thus any `0<gamma<1/2` yields a power-saving uniform perturbation of the
**scalar bridge response**.  This observation prevents an overly broad
no-go claim: generic perturbation is not killed merely by its energetic size.

Equation (GG.15) still does not supply a response quotient.  It says that if
one could compute or approximate the perturbed bridge response from the
overlap object, the answer would transfer back with sublinear error.  The
italicized premise is exactly what is missing.

## 4. What multi-species synchronization really supplies

For finitely many species `s in S`, Panchenko uses a countable dense set
`W subset [0,1]^S` and fields with covariance

```math
\mathbb E h_{N,w,p}(\sigma^1)h_{N,w,p}(\sigma^2)
=R_w(\sigma^1,\sigma^2)^p,
\qquad
R_w=\sum_s\lambda_sw_sR_s.                                \tag{GG.16}
```

The sum over `(w,p)` has bounded variance (at most `4` in the source), and its
pressure cost is at most `2s_N^2`.  The resulting multi-species
Ghirlanda--Guerra identities imply the exact synchronization theorem:

```math
R^s_{\ell,\ell'}=L_s(R_{\ell,\ell'})\quad\hbox{a.s.},
\qquad R=\sum_s\lambda_sR^s,                               \tag{GG.17}
```

where every `L_s:[0,1]->[0,1]` is nondecreasing and
`1/lambda_s`-Lipschitz.

For two children with positive limiting proportions, this offers a scalar
overlap order parameter plus two Lipschitz decoding functions, provided the
limiting joint species-overlap array satisfies the full multi-species GG
identity.  But there is an exact structural tradeoff:

1. the perturbation (GG.16) contains mixed monomials across the two species,
   so the perturbed joint child law is no longer the required product
   `mu_A tensor mu_D`;
2. perturbing the two children separately preserves that product, but supplies
   only two separate Ghirlanda--Guerra structures.  The hypotheses of the
   multi-species synchronization theorem are then absent, and no theorem makes
   the two independent overlap trees share one scalar phase.

Therefore (GG.17) cannot be imported while silently retaining the exact
product bridge channel.

## 5. Why the overlap quotient does not yet control the inverse escort

For the actual children, the object to be controlled is

```math
L(B)=\log\mathbb E_{\mu_A\otimes\mu_D\otimes U_\tau}
 \exp\{t\tau x^TBy\},
\qquad
q_\lambda(B)={e^{-\lambda L(B)}\over
 \mathbb E_Ue^{-\lambda L}}U(B).                           \tag{GG.18}
```

The exact cavity observables are coordinate-labelled quantities such as
`E_{G_B}(tau x_i y_j)`.  Replica overlaps retain only

```math
k^{-1}\langle x,x'\rangle,
```

and determine the asymptotic Gibbs measure only up to Hilbert-space isometry.
The bridge query is not invariant under that equivalence: an orthogonal change
of embedding sends `B` to `U^TBV`, while the allowed bridge set consists of
coordinatewise sign matrices and is not orthogonally invariant.  Thus a Gram
law, even the complete infinite replica Gram law, forgets the physical
coordinate embedding used by (GG.18).

There are two further quantifier failures.

- The Ghirlanda--Guerra conclusion is annealed over the auxiliary Gaussian
  perturbation and concerns finitely many replicas.  It is not a single
  quenched statement uniform over the `2^(mn)` bridge environments selected by
  `q_lambda`.
- Under the inverse escort the bridge environment is reweighted by a negative
  power of its partition function.  The joint density of `(B,spin)` contains
  `Z_B^{-(lambda+1)}`.  It is not the ordinary Gibbs measure to which the
  Gaussian integration-by-parts proof of (GG.10) applies.  Differentiating the
  perturbation also differentiates this negative partition power, producing
  additional escort terms.

Consequently, none of (GG.10)--(GG.17) proves a bound on `I^leftarrow`, `J`,
`J-I^leftarrow`, or the negative edge-cavity overlap.

## 6. Exact theorem/missing-theorem boundary

The literature gives the rigorous implication

```math
\boxed{
\begin{array}{c}
\text{arbitrary deterministic child}\ +\ o(N)\text{-pressure generic perturbation}
\\[2pt]\Longrightarrow\\[-2pt]
\text{GG identities}\ +\ \text{ultrametric replica Gram law},
\end{array}}
                                                               \tag{GG.19}
```

and, if one allows a joint mixed-species perturbation,

```math
\boxed{
\text{multi-species GG}\Longrightarrow
R^s=L_s(R)\text{ for every species}.}                         \tag{GG.20}
```

The genuinely new lemma needed for the bridge program would be a
**coordinate-response lifting theorem** of the following precise kind:

> For separately perturbed product children satisfying quantitative
> Ghirlanda--Guerra identities, construct from a discretized scalar overlap
> law (plus explicitly bounded coordinate data of subexponential complexity)
> an approximation to every inverse-escort row/cavity response with total
> error `o(N)` uniformly over the negatively tilted bridge law.

This is not furnished by ultrametricity or synchronization.  It must also
either synchronize the two independently perturbed child trees without mixed
perturbations, or prove that the mixed-species perturbation's induced
non-product dependence is harmless for the bridge variational problem.

## 7. Verdict

The route earns a **rigorous mapping, not a bridge RESET**.

- Positive result: actual deterministic children admit quantitative GG
  regularization with `O(k^{-gamma/2})` fixed-test error and a uniform
  `O_P(N^(1/2+gamma))` bridge-response perturbation, for any
  `0<gamma<1/2`.
- Genuine quotient: the limiting replica Gram law collapses to one scalar
  overlap distribution; multi-species arrays collapse to one scalar overlap
  plus finitely many Lipschitz functions.
- Decisive boundary: this quotient is not sufficient for the coordinatewise
  rank-one bridge query, is not uniform under the inverse escort, and joint
  synchronization conflicts with preservation of the child-product law.

Any continuation should attack the coordinate-response lifting theorem above,
not re-prove Ghirlanda--Guerra identities or ultrametricity.
