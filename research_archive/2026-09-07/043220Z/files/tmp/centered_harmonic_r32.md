# Wave 32 route 3: canonical centering and the exact screening criterion

## Status

- **Verified:** selector constants are a gauge.  There is a canonical gauge
  in which every selector component is a normalized parent likelihood, but
  no choice of constants changes `Phi`, its posterior, or any cut response.
- **Verified:** a parent-coordinate response is the log moment generating
  function of the exact screening residual.  It splits into a selector-common
  drift plus a posterior KL correction.  The posterior susceptibility sees
  neither an additive common residual nor the corresponding parent response.
- **Verified:** a heat-bath comparison plus the exact screened response energy
  is a sufficient fixed-density parent lemma.  This isolates two genuinely
  minimizer-specific inputs and includes the orientation bit.
- **Verified wall:** for the exact minimizer `A_6`, `m=3`, centered harmonic
  flatness diverges linearly as `beta -> infinity` while the actual parent
  entropy tends to zero exponentially.  Thus (10.921) is much stronger than
  the desired parent estimate even after the optimal selector normalization.
- **Open:** no target-scale comparison or screened response bound is proved
  uniformly for exact minimizers at
  `beta=Theta(n^(-1/2+c))`.  This route does not prove convergence.

All finite identities, the abstract wall, and the exact max-plus table are
checked by `tmp/centered_harmonic_r32.py`.

## 1. Selector constants have a canonical gauge, but no optimization power

Retain

\[
 z_S=\mathbb E_{\nu_\beta}e^{-F_S(D)}.
\]

Directly summing first over outside completions gives

\[
 \boxed{
 z_S=\frac{2^{n-m}Z_S(\beta)}{Z_A(\beta)}.
 }
\tag{R32.1}
\]

Choose

\[
 \boxed{
 \alpha_S^*=-\log z_S,\qquad
 h_S(d)=e^{-[F_S(d)-\alpha_S^*]}=\frac{e^{-F_S(d)}}{z_S},
 \qquad q(S)=\frac{z_S}{\sum_Tz_T}.
 }
\tag{R32.2}
\]

Then `E_nu h_S=1`, `q(S)=Z_S(beta)/sum_T Z_T(beta)`, and

\[
 \boxed{
 f_\beta(d)=\frac{e^{-\Phi(d)}}{\mathbb E_\nu e^{-\Phi}}
 =\mathbb E_{S\sim q}h_S(d).
 }
\tag{R32.3}
\]

Thus `h_S nu` is the normalized lifted selector channel and (R32.3) is its
exact barycenter.  This is the unique gauge which makes every component have
parent mean one (or is unique up to a common constant if only equal means are
required).

It does **not** optimize any fluctuation.  For arbitrary constants `alpha_S`,

\[
 q_\alpha(S)e^{-[F_S(d)-\alpha_S]}
 =\frac{e^{-F_S(d)}}{\sum_Te^{-\alpha_T}}.
\tag{R32.4}
\]

Consequently changing `alpha` changes the soft minimum by an additive
constant only.  In particular its oscillation, Doob energy, normalized
likelihood, and selector posterior are exactly gauge invariant.  The useful
role of (R32.2) is probabilistic normalization, not concentration.

## 2. Exact parent-flip response and centered posterior susceptibility

In the canonical gauge put

\[
 \pi_d(S)=\frac{q(S)h_S(d)}{f_\beta(d)}
 =\frac{e^{-F_S(d)}}{\sum_Te^{-F_T(d)}}.
\tag{R32.5}
\]

Let `T_i d=d^i` be a projective vertex flip and define

\[
 Z_{S,i}(d)=\log\frac{h_S(d^i)}{h_S(d)}
 =F_S(d)-F_S(d^i).
\]

The one-spin recursion (10.916) gives the promised screened form

\[
 \boxed{
 Z_{S,i}(d)=2\mathbf 1_{\{i\in S\}}
 \left[b^V_{S\setminus\{i\},i}(d[S\setminus\{i\}])-\beta a_{S\setminus\{i\},i}(d)\right].
 }
\tag{R32.6}
\]

The parent common response is

\[
 \boxed{
 \chi_i(d):=\log\frac{f_\beta(d^i)}{f_\beta(d)}
 =\log\mathbb E_{S\sim\pi_d}e^{Z_{S,i}(d)}.
 }
\tag{R32.7}
\]

More sharply, the flip exponentially tilts the selector posterior:

\[
 \boxed{
 \pi_{d^i}(S)=\pi_d(S)e^{Z_{S,i}(d)-\chi_i(d)},
 \qquad
 \chi_i(d)=\mathbb E_{\pi_d}Z_{S,i}(d)
 +D(\pi_d\Vert\pi_{d^i}).
 }
\tag{R32.8}
\]

The reverse identity is

\[
 \mathbb E_{\pi_{d^i}}Z_{S,i}(d)-\chi_i(d)
 =D(\pi_{d^i}\Vert\pi_d).
\tag{R32.9}
\]

Hence the natural centered posterior susceptibility is
`Var_(pi_d)(Z_(S,i))`, while the parent response is the log-normalizing
constant `chi_i`.  Adding the same number to every `Z_(S,i)` leaves the
posterior score `Z-chi`, its variance, both posterior divergences, and its
Hellinger distance unchanged, but shifts `chi` by that number.  Therefore no
bound using only selector-centered screening can control the parent common
mode.  One must control the signed combination in (R32.8), including its
cancellation, rather than separately upper-bounding the drift and KL terms.

This posterior change is along a **parent cut-coordinate flip**.  It is not
the adjacent-selector finite-measure Hellinger quantity in (10.800), which
compares `S` and `S'` along the Johnson graph.  The two targets remain
separate.

## 3. Exact local-screening sufficient criterion

Parameterize the oriented-cut cube by one orientation bit and `n-1`
projective vertex bits, with coordinate involutions `T_0,...,T_(n-1)`.  For
the orientation bit define (R32.7) with
`Z_(S,0)=log h_S(T_0d)-log h_S(d)`; (R32.6) is asserted only for vertex
bits.

Let

\[
 \mu_t(d)=\frac{\nu_\beta(d)f_\beta(d)^t}
 {\mathbb E_{\nu_\beta}f_\beta^t},\qquad 0\le t\le1.
\]

For an unordered coordinate edge `e={d,T_jd}`, put

\[
 c_t(e)=\frac{\mu_t(d)\mu_t(T_jd)}
 {\mu_t(d)+\mu_t(T_jd)}.
\]

The heat-bath Dirichlet energy of the endpoint log likelihood is exactly

\[
 \boxed{
 \mathcal E_t(\log f_\beta)
 =\sum_{j=0}^{n-1}\ \sum_{e\in E_j}c_t(e)\chi_j(d)^2.
 }
\tag{R32.10}
\]

The sharp endpoint-specific comparison coefficient is

\[
 C_{\rm scr}(t)=
 \frac{\operatorname{Var}_{\mu_t}(\log f_\beta)}
 {\mathcal E_t(\log f_\beta)},
\]

with the evident zero-energy convention.  Thus (10.900) is equivalently the
exact centered-susceptibility identity

\[
 \boxed{
 \operatorname{Ent}_{\nu_\beta}(f_\beta)
 =\int_0^1t\,C_{\rm scr}(t)
 \sum_{j,e\in E_j}c_t(e)\chi_j(d)^2\,dt.
 }
\tag{R32.10a}
\]

If `C_HB(mu_t)` is the heat-bath Poincare constant in the normalization

\[
 \operatorname{Var}_{\mu_t}(g)
 \le C_{\rm HB}(\mu_t)
 \sum_j\mathbb E_{\mu_t}\operatorname{Var}_{\mu_t}
 (g\mid D_{-j}),
\]

then `C_scr(t)<=C_HB(mu_t)`, so (10.900) and (R32.10) imply the rigorous
sufficient bound

\[
 \boxed{
 \operatorname{Ent}_{\nu_\beta}(f_\beta)
 \le\int_0^1 t\,C_{\rm HB}(\mu_t)
 \sum_{j,e\in E_j}c_t(e)
 \left[\log\mathbb E_{\pi_d}e^{Z_{S,j}(d)}\right]^2dt.
 }
\tag{R32.11}
\]

Thus a target-scale upper bound `O(n^(1/2-2c))` on the right side is an
exact fixed-density parent lemma.  It may be proved by a full Poincare bound,
or more weakly by a comparison valid only for `g=log f_beta`.  Generic
parent comparison is already blocked by (10.813); exact-minimizer structure
must control both the comparison factor and the **combined** screened
response.  Raw `a`, posterior variance alone, and separate absolute bounds on
the two terms in (R32.8) are not substitutes.

## 4. Two centered walls

### 4.1 A selector-common response is invisible to susceptibility

On a two-point parent space with uniform `nu`, take two duplicate selector
components

\[
 e^{-F_1}=e^{-F_2}=(3/4,1/4).
\]

Then `z_1=z_2=1/2`, `h_1=h_2=(3/2,1/2)`, and `q=(1/2,1/2)`.  The selector
posterior is identically `q`, so every posterior KL, Hellinger distance, and
centered screening variance is zero.  Nevertheless

\[
 \chi=\log(1/3),\qquad
 \operatorname{Ent}_\nu(f)
 =\frac12\left[\frac32\log\frac32+\frac12\log\frac12\right]
 =0.1308120359\ldots>0.
\tag{R32.12}
\]

This is an abstract likelihood wall, not a complete-signing endpoint.  It
rigorously rules out deriving the parent estimate from selector-centered
susceptibility by algebra alone.

### 4.2 Exact `A_6`: centered flatness can diverge while entropy vanishes

Take the exact order-six minimizer `A_6` and `m=3`.  For a parent cut `d`,
let

\[
 \lambda(d)=\min_{|S|=3}\max_{z\text{ extending }d[S]}H_{S,d[S]}(z).
\]

If `g_S(d)` is the multiplicity of the maximum outside energy, define the
leading coefficient

\[
 C(d)=\frac{2^3}{\binom63}
 \sum_{S:\,\max H_{S,d[S]}=\lambda(d)}\frac1{g_S(d)}.
\]

Exact enumeration gives

| parent energy | `lambda` | `C` | number of oriented cuts |
|---:|---:|---:|---:|
| `10` | `4` | `2/3` | `12` |
| `6` | `4` | `8/15` | `20` |
| `-6` | `4` | `2/15` | `20` |
| `-10` | `8` | `1` | `12` |

Since the state space is finite,

\[
 U_\beta(d)=C(d)e^{-\beta\lambda(d)}(1+o(1))
\]

uniformly in `d`.  It follows that

\[
 \boxed{
 \operatorname{osc}_d\Phi_\beta(d)
 =4\beta+\log(2/3)+o(1)\longrightarrow\infty.
 }
\tag{R32.13}
\]

On the other hand the parent Gibbs law is supported asymptotically on the
twelve energy-`10` cuts, where the leading coefficient is constant.  The
energy-`6` class has parent relative mass `(5/3)e^(-4 beta)` and endpoint
likelihood ratio tending to `(8/15)/(2/3)=4/5`.  Therefore

\[
 \boxed{
 \operatorname{Ent}_{\nu_\beta}(f_\beta)
 =\frac53\left[\frac15+\frac45\log\frac45\right]e^{-4\beta}
 +o(e^{-4\beta})\longrightarrow0.
 }
\tag{R32.14}
\]

No selector constants can alter (R32.13).  Thus centered harmonic flatness
(10.921) is not a faithful proxy for the parent entropy, even on a fixed
exact minimizer under matched temperatures.  This is a finite-parametric
mechanism wall, not an asymptotic falsifier at the project scaling.

At `beta=1/2,t=0`, the same exact enumeration gives

\[
 \operatorname{Var}_{\mu_0}\Phi=0.00436277919\ldots,
 \qquad \mathcal E_0(\log f)=0.01168882238\ldots,
\]

while the conductance-weighted squares of the two terms in (R32.8) are
`0.1789713344...` and `0.1615779804...`, and the corresponding posterior
variance sum is `0.6648398430...`.  The much smaller exact response is caused
by signed cancellation in (R32.8), not by either term being small.

## 5. Frontier

Optimizing `alpha_S` cannot advance the parent bound: it is only a gauge, and
the canonical choice (R32.2) merely exposes the barycenter likelihood.  The
strong flatness target (10.921) should be downgraded as a construction tool,
because `A_6` makes it diverge while the desired entropy vanishes.  The exact
surviving local route is (R32.11): prove a minimizer-specific comparison for
the particular endpoint log likelihood and control the log-mgf of the
screening residual `b-beta a`, including its selector-common drift and its
KL cancellation.  The adjacent-selector Hellinger estimate remains a
separate requirement.  No such target-scale theorem is presently known.
