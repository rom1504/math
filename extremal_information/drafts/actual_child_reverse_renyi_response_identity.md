# Reverse-Renyi work identity for the actual child channel

Status: **task-local rigorous theorem note, independently audited after the
endpoint-localization correction below**.  The
negative-disorder escort is not an arbitrary law on bridge matrices.  It is
the inverse-power escort of a weak binary channel whose planted word is the
rank-one child-spin matrix.  This note differentiates its reverse Renyi
divergence exactly.  A linear negative-disorder resource forces a macroscopic
density of leave-one-coordinate posterior response at an intermediate bridge
amplitude comparable to the physical amplitude; uniform disappearance of that
density rules out a linear resource.

The result uses the actual contracted-temperature child prior.  It does not
use a conference or Paley surrogate.  It does not by itself decide whether
the response phase occurs for optimizing children: the identity is an exact
classification, not an asymptotic closure.

## 1. Channel representation

Fix actual contracted-temperature children `A,D`, one relative orientation
`epsilon`, and let `mu_epsilon` be their zero-bridge augmented Gibbs law.
Write

```math
Q=(Q_{ij})=(\tau X_iY_j)\in\{-1,1\}^{m\times n},
\qquad d=mn.                                           \tag{RR.1}
```

The law of `Q` is centrally symmetric.  For a bridge channel amplitude
`u>=0`, put `rho=tanh u` and let

```math
W_u(B\mid Q)
=\prod_{e=1}^d{e^{uB_eQ_e}\over2\cosh u}.
                                                               \tag{RR.2}
```

If `Pi_u` is the resulting output law and `U` is the fair bridge law, its
likelihood ratio is

```math
p_u(B):={d\Pi_u\over dU}(B)
=\mathbb E_{\mu_\epsilon}
  \prod_{e=1}^d(1+\rho B_eQ_e).                       \tag{RR.3}
```

At the physical amplitude `u=t=beta/sqrt(N)`, this is exactly the sector
likelihood in (AC.5)--(AC.6), up to the sector normalization already peeled
there.  For `lambda>0`, define

```math
q_{\lambda,u}(B)
={p_u(B)^{-\lambda}\over\mathbb E_Up_u^{-\lambda}}U(B),
\qquad
\mathcal R_\lambda(u)
={1\over\lambda}\log\mathbb E_Up_u(B)^{-\lambda}.
                                                               \tag{RR.4}
```

Thus `R_lambda(u)=D_(1+lambda)(U||Pi_u)` and
`R_lambda(0)=0`.  Moreover, if

```math
L_u(B)=C_u+\log p_u(B),                                \tag{RR.5}
```

then its negative soft value is `C_u-R_lambda(u)`.

Delete coordinate `e` from the channel and define the *extrinsic planted
response*

```math
r_{e,u}(B_{-e})
=\mathbb E[Q_e\mid B_{-e}]
={\mathbb E_Q Q_e\prod_{f\ne e}(1+\rho B_fQ_f)
  \over
  \mathbb E_Q\prod_{f\ne e}(1+\rho B_fQ_f)}.          \tag{RR.6}
```

All denominators are strictly positive for finite `u`.

## 2. Exact reverse-Renyi work identity

**Theorem RR.1 (inverse-channel response identity).**  For every finite
prior on `Q`, every `lambda>0`, and every finite `u>=0`,

```math
\boxed{
\mathcal R_\lambda'(u)
=\sum_{e=1}^d\mathbb E_{q_{\lambda,u}}
 \left[
 {(1-\rho^2)r_{e,u}
  \{\rho r_{e,u}
    +\tanh(\lambda\operatorname{arctanh}(\rho r_{e,u}))\}
  \over1-\rho^2r_{e,u}^2}
 \right].}                                             \tag{RR.7}
```

In particular, if

```math
S_\lambda(u)
=\sum_{e=1}^d\mathbb E_{q_{\lambda,u}}
                  r_{e,u}(B_{-e})^2,                  \tag{RR.8}
```

then, with

```math
c_\lambda=1+\min\{\lambda,1\},
\qquad C_\lambda=1+\max\{\lambda,1\},
```

one has the pointwise sandwich

```math
\boxed{
c_\lambda\rho(1-\rho^2)S_\lambda(u)
\le\mathcal R_\lambda'(u)
\le C_\lambda\rho S_\lambda(u).}                    \tag{RR.9}
```

Consequently,

```math
\boxed{
c_\lambda\int_0^t\tanh u\,\operatorname{sech}^2u\,
                 S_\lambda(u)\,du
\le\mathcal R_\lambda(t)
\le C_\lambda\int_0^t\tanh u\,S_\lambda(u)\,du.}  \tag{RR.10}
```

*Proof.*  Introduce separate channel amplitudes and differentiate on their
diagonal.  With coordinate `e` deleted, (RR.3) factors as

```math
p_u(B)=p_{u,-e}(B_{-e})(1+\rho B_er_{e,u}(B_{-e})).   \tag{RR.11}
```

Therefore, conditionally on `B_(-e)`, the two inverse-escort weights are
proportional to `(1+rho b r)^(-lambda)`.  Put `a=rho r` and
`z=arctanh a`.  Their conditional normalizer is

```math
{(1+a)^{-\lambda}+(1-a)^{-\lambda}\over2}
=(1-a^2)^{-\lambda/2}\cosh(\lambda z).               \tag{RR.12}
```

Since

```math
\mathcal R_\lambda'(u)
=-\mathbb E_{q_{\lambda,u}}\partial_u\log p_u(B),    \tag{RR.13}
```

differentiating (RR.12), dividing by `lambda`, and summing the separate
coordinate derivatives gives (RR.7).  This step keeps the other channel
coordinates fixed, so no derivative of the extrinsic response is omitted.

For `z>=0`, concavity and monotonicity of `tanh` give

```math
\min\{\lambda,1\}\tanh z
\le\tanh(\lambda z)
\le\max\{\lambda,1\}\tanh z.                        \tag{RR.14}
```

Use oddness for `z<0`, `tanh z=a`, and

```math
1-\rho^2\le1-\rho^2r^2\le1                          \tag{RR.15}
```

to obtain (RR.9).  Integration from zero proves (RR.10). `square`

## 3. Extensive reverse Renyi means dense posterior response

For comparable splits, `d=Theta(N^2)` and
`t=beta/sqrt(N)`.  Equation (RR.10) immediately gives a sharp scale
dichotomy.

The identity also contains the ordinary reverse-KL endpoint.  Define

```math
\mathcal R_0(u):=D(U\Vert\Pi_u)=-\mathbb E_U\log p_u(B).
```

Taking `lambda downarrow 0` in (RR.7) gives

```math
\boxed{
\mathcal R_0'(u)
=\sum_e\mathbb E_U
 {\rho(1-\rho^2)r_{e,u}^2
  \over1-\rho^2r_{e,u}^2}.}                          \tag{RR.15a}
```

Consequently the centered negative-disorder gain has the exact work
difference

```math
\boxed{
G_\lambda(t)=\mathcal R_\lambda(t)-\mathcal R_0(t)
=\int_0^t
   \{\mathcal R_\lambda'(u)-\mathcal R_0'(u)\}\,du.} \tag{RR.15b}
```

Thus the part not already paid by ordinary reverse-KL compensation is
precisely the integrated enhancement of planted-coordinate response under
the inverse escort.  This is an identity, not a claim that the integrand is
pointwise nonnegative.

**Corollary RR.2 (posterior-response density dichotomy).**  Fix
`beta,lambda>0` and a balanced split window.  If

```math
\mathcal R_\lambda(\beta/\sqrt N)\ge\eta N,          \tag{RR.16}
```

then for some `0<u_N<=beta/sqrt(N)`,

```math
\boxed{
{S_\lambda(u_N)\over mn}
\ge {\eta N\over
 C_\lambda mn\log\cosh(\beta/\sqrt N)}
=\Omega_{\eta,\beta,\lambda}(1).}                  \tag{RR.17}
```

Moreover, the point can be chosen away from the zero-amplitude endpoint.
With

```math
\alpha=\min\left\{{1\over2},
        \sqrt{{2\eta\over C_\lambda\beta^2}}\right\},
```

one can choose `u_N in [alpha beta/sqrt(N),beta/sqrt(N)]`, with the localized
bound obtained by replacing the numerator on the right side of (RR.17) by
`3 eta N/4`.

Conversely, if

```math
\sup_{0\le u\le\beta/\sqrt N}{S_\lambda(u)\over mn}=o(1),  \tag{RR.18}
```

then `R_lambda(beta/sqrt(N))=o(N)`.

Thus a linear reverse-Renyi resource cannot be carried by a sparse set of
posterior-responsive bridge coordinates.  A positive fraction of all
`Theta(N^2)` coordinates has nonvanishing squared extrinsic response on
average at some intermediate amplitude.

*Endpoint localization.*  Put `t=beta/sqrt(N)`.  Since `mn<=N^2/4` and
`log cosh v<=v^2/2`,

```math
mn\log\cosh(\alpha t)\le {\eta N\over4C_\lambda}.
```

The upper bound in (RR.10) and (RR.16) therefore imply

```math
\int_{\alpha t}^t\tanh u\,S_\lambda(u)\,du
\ge {3\eta N\over4C_\lambda}.
```

Comparison with
`int_(alpha t)^t tanh u du<=log cosh t` proves the strengthened form of
(RR.17).  This localization is needed below: density at an amplitude much
smaller than `N^(-1/2)` would not by itself yield extensive mutual
information.

There is an equivalent statement intrinsic to the inverse escort.  Let

```math
s_{e,u}(B_{-e})
=\mathbb E_{q_{\lambda,u}}[B_e\mid B_{-e}].           \tag{RR.19}
```

The same two-point calculation gives

```math
\boxed{
s_{e,u}
=-\tanh\!\left(\lambda\operatorname{arctanh}
                         (\rho r_{e,u})\right),}      \tag{RR.20}
```

and hence

```math
\min\{\lambda,1\}\rho|r_{e,u}|
\le|s_{e,u}|
\le\max\{\lambda,1\}\rho|r_{e,u}|.                 \tag{RR.21}
```

Because global bridge inversion leaves the actual sector pressure invariant,
every bridge bit is marginally fair under `q_(lambda,u)`.  Therefore

```math
I_{q_{\lambda,u}}(B_e;B_{-e})
=\mathbb E\left[
 {1+s_{e,u}\over2}\log(1+s_{e,u})
 +{1-s_{e,u}\over2}\log(1-s_{e,u})\right].          \tag{RR.22}
```

For `|s|<=tanh(lambda t)=O(N^(-1/2))`, the bracket is
`s^2/2+O(s^4)` uniformly; in fact it is at least `s^2/2` for every
`|s|<=1`.  At the localized point `u_N>=alpha beta/sqrt(N)`, (RR.17) and
(RR.21) therefore force

```math
\boxed{
\sum_{e=1}^{mn}I_{q_{\lambda,u_N}}(B_e;B_{-e})
=\Omega(N).}                                         \tag{RR.23}
```

This is extensive leave-one-out dependence spread over quadratically many
weak coordinates.  It is stronger than the generic effective-support lower
bound (AC.28) because it identifies the dependence with posterior recovery
of the actual rank-one child word.  It is not yet a cross-*row* lower bound:
part of (RR.23) may be stored within rows and hence captured by the
row-product shadow.

The last ambiguity has an exact information-theoretic split.  Write
`R_i=(B_(ij))_(j<=n)` and `R_(i,-j)=R_i\setminus{B_(ij)}` and define

```math
\begin{aligned}
\mathsf W_{\lambda,u}
 &=\sum_{i,j}I(B_{ij};R_{i,-j}),\\
\mathsf C_{\lambda,u}
 &=\sum_{i,j}I(B_{ij};R_{-i}\mid R_{i,-j}).           \tag{RR.23a}
\end{aligned}
```

**Corollary RR.3 (within-row or irreducible cross-row response).**  At the
localized amplitude supplied by Corollary RR.2,

```math
\boxed{
\mathsf W_{\lambda,u_N}+\mathsf C_{\lambda,u_N}
=\sum_{i,j}I(B_{ij};B_{-(ij)})=\Omega(N).}            \tag{RR.23b}
```

Consequently at least one of `W` and `C` is extensive.  The first resource
is entirely internal to the row factors allowed in (AC.15).  The second is
genuinely cross-row: because `q_(lambda,u)` has full support,

```math
\boxed{
\mathsf C_{\lambda,u}=0
\quad\Longleftrightarrow\quad
q_{\lambda,u}=\bigotimes_i(q_{\lambda,u})_{R_i}.}     \tag{RR.23c}
```

*Proof.*  The mutual-information chain rule applied to
`B_(-(ij))=(R_(i,-j),R_(-i))` gives (RR.23b) term by term.  For (RR.23c),
the reverse implication is immediate.  In the forward direction, every
summand in `C` vanishes.  Hence the full conditional law of each bit
`B_(ij)`, given all other bits, depends only on the other bits in row `i`.
For a strictly positive law, expand `log q` in the Walsh basis.  The
conditional log odds of bit `(i,j)` contains every Walsh coefficient whose
support contains `(i,j)`.  Independence from all other rows forces every
coefficient meeting two rows to vanish.  Thus `log q` is a constant plus a
sum of row functions, which is exactly row-product factorization. `square`

This is a qualitative zero-set characterization, not a quantitative
comparison between `C` and the directed projection
`I_lambda^leftarrow`.  Establishing such a comparison requires a row
functional inequality; it does not follow from tight conditional row
Renyi-two alone.

## 4. Relation to the target and the remaining gap

At the physical amplitude, (RR.5) gives the exact identities

```math
V_\lambda=C_t-\mathcal R_\lambda(t),
\qquad
G_\lambda
=\mathcal R_\lambda(t)-D(U\Vert\Pi_t).              \tag{RR.24}
```

Thus every target-reaching inverse-disorder phase has an exact planted-
channel interpretation.  Corollary RR.2 says that a linear Renyi
compensation forces dense posterior response under the inverse escort.
Together with Lemma SH.0, there are now two complementary response tests:

- submacroscopic parent response on one exponentially typical convex
  uniform carrier kills the centered negative gain;
- submacroscopic extrinsic planted response throughout the channel path
  kills the entire reverse-Renyi compensation.

Neither statement decides the row-product target excess `Delta_N`.  The
precise unresolved refinement is:

> **Cross-row posterior-response lemma.**  Split the extensive quantity in
> (RR.23) using (RR.23a).  Prove, for actual optimizing children, either
> that `C_(lambda,u)=o(N)` along the whole channel path together with a
> quantitative row tensorization that transfers the remaining work to the
> bounded-row product shadow, or that `C_(lambda,u)=Omega(N)` on a set of
> amplitudes of the required weight.

This is narrower than classifying the full bridge law: it asks for one
integrated two-replica/leave-one-out response resource.  But a theorem
identifying that split with the directed row-product excess is still
missing, so this note alone does not close the Level-5-to-6 transition.
