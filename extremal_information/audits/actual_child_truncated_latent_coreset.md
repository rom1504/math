# A truncated global latent coreset and the escaping-collision boundary

**Status.** Rigorous theorem for the exact rank-one bridge channel.  This
note weakens the annealed collision hypothesis in Theorem 37.75 to a
quenched tail condition.  A large value of the mean posterior collision
factor is harmless if it is carried by a sufficiently small set under the
declared bridge law.  The theorem is valid for the actual inverse escort and
for any declared negative-path mixture, but it does not prove the required
tail bound for optimizing children.

## 1. Setup

Use the notation of the global latent-coreset audit.  Thus `d=mn`,
`t=beta/sqrt(N)`, `mu` is the fixed zero-bridge latent child law,
`P=P_t` is its normalized bridge likelihood, and `q` is any declared bridge
law.  Let

```math
K_0(B)=1+\chi^2(\mu_B\Vert\mu)
```

be the complete posterior collision factor.  For a threshold `H>=1`, put

```math
G_H=\{B:K_0(B)\le H\},\qquad
\delta_H=q(G_H^{\mathsf c}).                     \tag{TC.1}
```

Draw one global iid empirical latent law `mu_R` of size `R` from `mu`, and
write `P_R` and `r_e^(R)` for its likelihood and deleted-edge cavities.
The same sample is used at every bridge word and for every edge.

## 2. Tail-sensitive coreset theorem

### Theorem TC.1 (truncated global latent coreset)

For every `H>=1` and integer `R>=1`,

```math
\boxed{
E_{\mu_R}E_q\sum_{e=1}^{d}(r_e^{(R)}-r_e)^2
\le {32e^{4t}dH\over R}+4d\delta_H.}            \tag{TC.2}
```

Moreover,

```math
\boxed{
\begin{aligned}
E_{\mu_R}E_q|\log P_R-\log P|
\le{}&3\sqrt{H/R}+{4H\log(2R)\over R}\\
&+2td\{e^{-R/(4H)}+\delta_H\}.
\end{aligned}}                                  \tag{TC.3}
```

Some single deterministic empirical law satisfies both conclusions with a
common factor at most two.  Its cavity field is exactly integrable:

```math
\nabla_e\log P_R=\operatorname {arctanh}(\rho r_e^{(R)}),
\qquad \rho=\tanh t.                             \tag{TC.4}
```

*Proof.*  On `G_H`, the one-edge comparison from Lemma GC.3 gives

```math
K_e(B_{-e})\le e^{4t}K_0(B)\le e^{4t}H
```

simultaneously for every edge.  The pointwise empirical-posterior estimate
in Theorem GC.2 is therefore at most `32e^(4t)H/R` per edge.  On the
complement, both cavities lie in `[-1,1]`, so their squared difference is at
most four.  Averaging proves (TC.2).

At a fixed `B` in `G_H`, apply the positive sample-mean logarithm bound
GC.5 to

```math
W={k_Q(B)\over P(B)}.
```

It has mean one, second moment `K_0(B)<=H`, and
`W>=exp(-2td)`.  This gives the first three terms in (TC.3), uniformly on
`G_H`.  On its complement, both `P_R` and `P` are mixtures of kernels in
`[(1-rho)^d,(1+rho)^d]`; hence

```math
|\log(P_R/P)|\le2td.
```

This contributes `2td delta_H` and proves (TC.3).  Equation (TC.4) is the
exact deleted-edge insertion identity.  Applying the probabilistic method
to the sum of the two normalized errors gives one common realization.
`square`

### Corollary TC.2 (polynomial escaping mass is enough)

Assume comparable splits.  Suppose that for some fixed `zeta>0` there are
thresholds

```math
\log H_N=o(N),\qquad
\delta_{H_N}=O(N^{-1/2-\zeta}).                 \tag{TC.5}
```

Set

```math
R_N=\left\lceil
16H_NN^{1/2+\zeta}(\log N)^2
\right\rceil.                                   \tag{TC.6}
```

Then `R_N=exp{o(N)}` and one global, curl-free empirical channel obeys

```math
\boxed{
t^2E_q\sum_e(r_e^{(R_N)}-r_e)^2
=O_\beta(N^{1/2-\zeta}),
\qquad
E_q|\log P_{R_N}-\log P|=o(N).}                 \tag{TC.7}
```

*Proof.*  At comparable splits, `t^2d=Theta_beta(N)` and
`td=Theta_beta(N^(3/2))`.  Substitute (TC.5)--(TC.6) into (TC.2): the good
part is `O_beta(N^(1/2-zeta)/(log N)^2)` and the escaping part is
`O_beta(N^(1/2-zeta))`.  In (TC.3), the exponential term is negligible,
the escaping term is `O_beta(N^(1-zeta))`, and the remaining two terms are
`o(N)` because `log H_N=o(N)`.  `square`

## 3. Why this is strictly weaker than the annealed condition

By the two-sided one-edge comparison GC.3, Theorem 37.75's condition on
`overline K_del` is equivalent at exponential scale to

```math
\log E_qK_0=o(N).                                \tag{TC.8}
```

Markov's inequality turns (TC.8) into (TC.5) after multiplying the
threshold by a suitable polynomial.  The converse is false in the ambient
scalar-random-variable class: `K_0` may be `exp(cN)` on a set of
`q`-mass `N^(-1/2-zeta)` and bounded elsewhere.  Then (TC.5) holds with a
constant threshold while `E_qK_0` is exponential.  Thus an extensive
annealed collision exponent, by itself, does not falsify a reusable
typical-response coreset.

This logical strictness has not been realized or proved inside the class of
actual optimizing-child channels.  The tail condition is therefore a
strictly weaker sufficient condition, not yet a theorem that it is
optimizer-operationally easier.

The new sufficient condition has an exact operational interpretation.  It
asks for a subexponential posterior-collision threshold outside only the
amount of bridge mass which the physical `t^2d=Theta(N)` cavity scale and
the deterministic `td=Theta(N^(3/2))` log-likelihood range can afford.  It
does not request the complete cavity table, a generic cumulant hierarchy,
or independent approximations at different bridge words.

## 4. Updated boundary

The smallest collision question relevant to the existing coreset
architecture is therefore no longer the annealed mean.  It is the actual
escort tail:

> Find `H_N=exp{o(N)}` for which
> `q\{K_0>H_N\}=O(N^(-1/2-zeta))` for some `zeta>0`, uniformly on the
> required negative path; or prove that every subexponential threshold is
> exceeded on nonvanishing actual escort mass and identify the localized
> child phase.

The positive branch gives a power-saving physical cavity approximation and
`o(N)` scalar-pressure error.  It still does not provide target reach, a
mergeable `o(N)`-bit state, or the one-row directional cost required by
Theorem 37.76.  Those obligations remain separate.
