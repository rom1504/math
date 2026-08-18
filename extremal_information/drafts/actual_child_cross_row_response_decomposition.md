# Cross-row posterior response: an explicit product certificate and its obstruction

Status: **task-local rigorous theorem and finite falsifier**.  This note uses
the actual optimized-child channel of
[`actual_child_reverse_renyi_response_identity.md`](actual_child_reverse_renyi_response_identity.md).
It gives a nonvariational row-product certificate whose error is exactly a
centered cross-row interaction cumulant, and then writes that interaction as
an integral of full-versus-row-erased posterior responses.  It also proves
that reverse-Renyi work nonadditivity alone cannot represent the directed
row-product gap, even for a centrally symmetric rank-one channel.

The theorem applies verbatim to the actual child prior.  The finite
falsifier is a rank-one channel prior, but is not asserted to arise from
optimized children.

## 1. Full and row-erased channels

Fix the actual prior on

```math
Q=\tau XY^T\in\{\pm1\}^{m\times n},
```

and a channel amplitude `u`.  Write `p_u=dPi_u/dU` for the full output
likelihood and let

```math
p_{i,u}(b_i)=E_{U_{-i}}p_u(b_i,B_{-i})                \tag{CR.1}
```

be its `i`th row marginal likelihood.  Thus `p_(i,u)` is obtained by erasing
every channel row except `i`; it is not an approximation.  Put

```math
\begin{aligned}
 Z_u&=E_Up_u^{-\lambda},
 &\mathcal R_u&={1\over\lambda}\log Z_u,\\
 Z_{i,u}&=E_{U_i}p_{i,u}^{-\lambda},
 &\mathcal R_{i,u}&={1\over\lambda}\log Z_{i,u},      \tag{CR.2}\\
 {dq_u\over dU}&={p_u^{-\lambda}\over Z_u},
 &{dr_{i,u}\over dU_i}&={p_{i,u}^{-\lambda}\over Z_{i,u}},\\
 r_u&=\bigotimes_{i=1}^m r_{i,u}.                    \tag{CR.3}
\end{aligned}
```

Finally define the forward row-interaction information density

```math
h_u(B)=\log p_u(B)-\sum_{i=1}^m\log p_{i,u}(B_i).      \tag{CR.4}
```

It is the logarithm of the full forward output density divided by the
product of its row marginals.  In particular `E_(Pi_u)h_u` is the forward
row total correlation, although the measure relevant below is `r_u`, not
`Pi_u`.

## 2. Exact canonical row-product decomposition

**Theorem CR.1 (capturable row work plus centered interaction).**  At every
finite `u` and `lambda>0`,

```math
\boxed{
\mathcal R_u
=\left(\sum_i\mathcal R_{i,u}-E_{r_u}h_u\right)
 +{1\over\lambda}\mathcal J_u,}                      \tag{CR.5}
```

where

```math
\boxed{
\mathcal J_u
=D(r_u\Vert q_u)
=\log E_{r_u}\exp\{-\lambda(h_u-E_{r_u}h_u)\}.}      \tag{CR.6}
```

The first parenthesis in (CR.5) is exactly the inverse-work value achieved
by the explicit row-product law `r_u`.  Consequently the optimal directed
row-product excess satisfies

```math
\boxed{
\mathcal I_{\lambda,u}^{\leftarrow}
=\inf_{p=\otimes_i p_i}D(p\Vert q_u)
\le\mathcal J_u.}                                    \tag{CR.7}
```

Moreover, every factor of this explicit certificate has the same tight
conditional complexity as the optimal shadow:

```math
\boxed{
D_2(r_{i,u}\Vert U_i)
\le n\log(1+\tanh^2(\lambda u))
\le\lambda^2u^2n.}                                   \tag{CR.8}
```

*Proof.*  Since `p_u=(prod_i p_(i,u))e^(h_u)`,

```math
Z_u=\left(\prod_iZ_{i,u}\right)E_{r_u}e^{-\lambda h_u}.
                                                               \tag{CR.9}
```

Taking logarithms gives the uncentered form of (CR.5).  Also

```math
{dq_u\over dr_u}
={e^{-\lambda h_u}\over E_{r_u}e^{-\lambda h_u}},
```

so direct evaluation of `D(r_u||q_u)` gives (CR.6).  For the potential
`-lambda log p_u`, the Gibbs variational value at `r_u` is

```math
\lambda\left(\sum_i\mathcal R_{i,u}-E_{r_u}h_u\right).
```

The unrestricted value is `lambda R_u`; minimizing KL over row products
therefore proves (CR.7).  Finally, flipping one output bit changes
`log p_(i,u)` by at most `2u`.  Lemma AC.1 applied to its inverse escort
proves (CR.8). `square`

This is not a reformulation of the unknown product optimum.  The law `r_u`
is explicit from the `m` row-erased forward channels, and (CR.7) is an
upper certificate that may be nonoptimal.  In particular,

```math
\mathcal I_{\lambda,u}^{\leftarrow}\ge\eta N
\quad\Longrightarrow\quad
\log E_{r_u}e^{-\lambda(h_u-E_{r_u}h_u)}\ge\eta N.    \tag{CR.10}
```

Thus actual irreducible row dependence forces a linear lower-tail cumulant
of one explicit cross-row observable under a product law with bounded row
Renyi-two complexity.

## 3. Reverse-Renyi cross work is an exact posterior-response integral

For a bridge bit `e=(i,j)`, define two extrinsic planted responses at
amplitude `v`:

```math
\begin{aligned}
a^{\rm all}_{e,v}(B_{-e})
 &=E[Q_e\mid B_{-e}\text{ through all channel rows}],\\
a^{\rm row}_{e,v}(B_{i,-e})
 &=E[Q_e\mid B_{i,-e}\text{ through row }i\text{ only}].       \tag{CR.11}
\end{aligned}
```

Let

```math
\Psi_{\lambda,v}(a)
={(1-\rho^2)a\{\rho a+
 \tanh(\lambda\operatorname{arctanh}(\rho a))\}
  \over1-\rho^2a^2},
\qquad\rho=\tanh v.                                  \tag{CR.12}
```

Applying Theorem RR.1 once to the full channel and once to each row-erased
channel gives the exact signed cross-response rate

```math
\begin{aligned}
\mathcal C_\lambda(v)
:={d\over dv}\left(\mathcal R_v-\sum_i\mathcal R_{i,v}\right)
&=\sum_eE_{q_v}\Psi_{\lambda,v}(a^{\rm all}_{e,v})\\
&\quad-\sum_i\sum_{e\in i}
 E_{r_{i,v}}\Psi_{\lambda,v}(a^{\rm row}_{e,v}),      \tag{CR.13}
\end{aligned}
```

and hence

```math
\boxed{
\mathcal R_t-\sum_i\mathcal R_{i,t}
=\int_0^t\mathcal C_\lambda(v)\,dv.}                 \tag{CR.14}
```

This is the requested exact split of reverse-Renyi work into isolated-row
response and response visible only when the other rows are present.  It is
important that `C_lambda` is signed; Section 5 gives a rank-one example
where its integral is negative.

## 4. The missing mean interaction is also a row-erasure response

The term `E_(r_t)h_t` in (CR.5) is essential.  It too has an exact response
path.  Define the forward coordinate score

```math
\chi_v(b,a)={(1-\rho^2)ba\over1+\rho ba}.             \tag{CR.15}
```

Separate coordinate amplitudes before differentiating.  The usual deleted-
coordinate factorization gives

```math
{d\over dv}\log p_v(B)
=\sum_e\chi_v(B_e,a^{\rm all}_{e,v}),
```

and the analogous identity for each row marginal.  Since `h_0=0`,

```math
\boxed{
h_t(B)=\int_0^t\sum_e
 \left[\chi_v(B_e,a^{\rm all}_{e,v})
       -\chi_v(B_e,a^{\rm row}_{e,v})\right]dv.}      \tag{CR.16}
```

Combining (CR.5), (CR.14), and (CR.16) yields a completely explicit
row-erasure decomposition of the canonical certificate:

```math
\boxed{
\mathcal J_t
=\lambda\int_0^t\left\{\mathcal C_\lambda(v)
 +E_{r_t}\sum_e
 [\chi_v(B_e,a^{\rm all}_{e,v})
  -\chi_v(B_e,a^{\rm row}_{e,v})]\right\}dv.}         \tag{CR.17}
```

Here the expectation in the second line is under the explicit **final-
amplitude** product `r_t`; no derivative of that measure is being omitted.

There is also a usable one-sided criterion.  Since

```math
\left|{\partial\over\partial a}\chi_v(b,a)\right|
={1-\rho^2\over(1+\rho ba)^2}\le e^{2v},              \tag{CR.18}
```

put

```math
\mathcal A_t
=\sup_B\int_0^t e^{2v}\sum_e
 |a^{\rm all}_{e,v}(B_{-e})
  -a^{\rm row}_{e,v}(B_{i,-e})|\,dv.                 \tag{CR.19}
```

Then `|h_t(B)|<=A_t`.  The trivial centered-range bound and Hoeffding's
lemma in (CR.6) give

```math
\boxed{
\mathcal I_{\lambda,t}^{\leftarrow}
\le\mathcal J_t
\le\min\left\{2\lambda\mathcal A_t,
 {\lambda^2\mathcal A_t^2\over2}\right\}.}          \tag{CR.20}
```

At the physical scale `t=beta/sqrt(N)`, the concrete synchronization lemma

```math
\sup_{0\le v\le t}\sup_B
 \sum_e|a^{\rm all}_{e,v}-a^{\rm row}_{e,v}|
=o(N^{3/2})                                           \tag{CR.21}
```

implies `A_t=o(N)` and therefore
`I_(lambda,t)^leftarrow=o(N)`.  Conversely, if
`I_(lambda,t)^leftarrow>=eta N`, then (CR.20) forces
`A_t>=eta N/(2lambda)`, and hence some amplitude and bridge have

```math
\sum_e|a^{\rm all}_{e,v}-a^{\rm row}_{e,v}|
=\Omega_{\eta,\beta,\lambda}(N^{3/2}).               \tag{CR.22}
```

Thus the actual directed row dependence has a genuine posterior-response
dichotomy: either full and row-erased planted responses synchronize at the
`N^(-1/2)` average-coordinate scale, or a linear reverse-product resource
forces a macroscopic witness.  Uniformity in `B` is strong, but it protects
against exactly the exponentially rare bridges that a bounded-temperature
average can miss.

## 5. Falsifier: cross work alone can have the wrong sign

The term `C_lambda` in (CR.13) cannot by itself identify
`I^leftarrow`.  Consider two rows of two bits and the eight rank-one words
`Q=xy^T`.  Give the four central orbits `{Q,-Q}` represented by

```math
\begin{pmatrix}1&1\\1&1\end{pmatrix},\quad
\begin{pmatrix}1&-1\\1&-1\end{pmatrix},\quad
\begin{pmatrix}1&1\\-1&-1\end{pmatrix},\quad
\begin{pmatrix}1&-1\\-1&1\end{pmatrix}              \tag{CR.23}
```

total orbit weights

```math
{3\over80},\quad {17\over400},\quad {1\over25},
\quad {22\over25},                                   \tag{CR.24}
```

split equally between `Q` and `-Q`.  This is a centrally symmetric
rank-one prior.  Pass it through the binary channel with
`rho=4/5` and take `lambda=1`.  In row-state order
`(--),(-+),(+-),(++)`, the exact output likelihood is

```math
p={1\over625}
\begin{pmatrix}
201&369&369&209\\
369&217&2897&369\\
369&2897&217&369\\
209&369&369&201
\end{pmatrix},                                       \tag{CR.25}
```

and both row marginal likelihoods are

```math
p_1=p_2={1\over625}(287,963,963,287).                 \tag{CR.26}
```

Direct rational summation gives

```math
Z={6484767816875\over3248296648443},
\qquad
Z_1=Z_2={390625\over276381},                         \tag{CR.27}
```

so

```math
{Z\over Z_1Z_2}
={306836044316068869\over307023296142578125}<1.       \tag{CR.28}
```

Therefore

```math
\int_0^t\mathcal C_1(v)dv
=\mathcal R_t-\mathcal R_{1,t}-\mathcal R_{2,t}<0.   \tag{CR.29}
```

Nevertheless `q proportional p^(-1)` is not row-product: for example the
upper-left `2`-by-`2` minor in (CR.25) has
`201*217 != 369^2`.  Since row products form a compact closed set,

```math
\mathcal I_{1,t}^{\leftarrow}>0.                     \tag{CR.30}
```

The canonical centered interaction is also positive by strictness of KL
(`J_t=D(r_t||q_t)>0`; numerically it is about `0.16443`).  Thus no theorem
of the form “positive cross reverse-Renyi work equals irreducible row
dependence” is valid for general centrally symmetric rank-one channels.
The mean response term in (CR.17), or equivalently the centered cumulant in
(CR.6), is indispensable.

## 6. New smallest missing lemma

The exact actual-child target is now bounded by a concrete product
certificate rather than by another optimization:

```math
\boxed{
\mathcal I_{\lambda,t}^{\leftarrow}
\le
\log E_{r_t}\exp\{-\lambda(h_t-E_{r_t}h_t)\}.}        \tag{CR.31}
```

The strictly narrower missing statement is an **actual-child cross-response
synchronization/concentration lemma**: prove either (CR.21), or directly
that the centered row-interaction density `h_t` has `o(N)` lower-tail
log-MGF under the explicit bounded-`D_2` product `r_t`.  Conversely, a
linear `I^leftarrow` forces the witness (CR.22).  The finite example shows
that optimizer-specific structure is essential if one tries to use the
signed work difference (CR.14) alone.
