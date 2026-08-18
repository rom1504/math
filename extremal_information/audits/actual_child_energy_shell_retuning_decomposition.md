# Actual-child retuning splits into energy-shell and within-shell geometry

Status: **rigorous exact actual-child theorem**.  The quadratic spread
theorem excludes low-rate atomic and narrow-cluster retuning.  This note
identifies what remains without restoring the full child table: every
negative-path latent retuning splits exactly into a polynomial-state
combined-energy response and a conditional geometric entropy deficit inside
equal-energy shells.

The theorem is a strict information reduction for the latent posterior.  It
does not identify the row-product optimizer, bound the row lifetime, or
prove target reach.

## 1. The actual conditional rank-one prior is shell-uniform

Let `A,D` be actual pressure-minimizing children of orders `m,n` at common
raw temperature `t`, and condition on relative orientation `epsilon`.  If
`Q=xy^T`, define

```math
 E_\epsilon(Q)=H_A(x)+\epsilon H_D(y).                 \tag{ES.1}
```

This is well-defined: a rank-one sign matrix determines the projective
classes `[x],[y]`, and quadratic energies are invariant under global spin
sign.  The exact law (PA.11) is

```math
 \mu_\epsilon(Q)
 ={2^{2-m-n}\over\mathcal Z_\epsilon}
   \cosh(tE_\epsilon(Q)).                              \tag{ES.2}
```

For an attainable value `e`, put

```math
 \mathcal S_e=\{Q:E_\epsilon(Q)=e\},
 \qquad p_e=\mu_\epsilon(\mathcal S_e).                \tag{ES.3}
```

Equation (ES.2) shows exactly that

```math
 \boxed{
 \mu_\epsilon(\,\cdot\mid E_\epsilon=e)
 =U_{\mathcal S_e}.}                                  \tag{ES.4}
```

Moreover, each child energy is an integer in an interval of length at most
twice its number of edges.  Hence

```math
 |\operatorname {range}E_\epsilon|
 \le {m\choose2}+{n\choose2}+1=O(N^2).                \tag{ES.5}
```

Thus the shell distribution is polynomial-state data, not an exponential
spin or flip table.  It is determined by the two signed child-energy
histograms and their convolution.

It also has a uniform linear surprise ceiling.  Put

```math
 C_\beta=\log2+{\beta^2\over4}.                       \tag{ES.5a}
```

The sector sums obey

```math
 \mathcal Z_\epsilon
 \le(Z_A^++Z_A^-)(Z_D^++Z_D^-)
 =4Z_A(t)Z_D(t)
 \le4(\cosh t)^{{m\choose2}+{n\choose2}}.             \tag{ES.5b}
```

Using `cosh(tE)>=1` in (ES.2), `log cosh t<=t^2/2`, and
`m^2+n^2<=N^2` gives, for every supported word and hence every nonempty
shell,

```math
 \boxed{
 \mu_\epsilon(Q)\ge e^{-C_\beta N},
 \qquad p_e\ge e^{-C_\beta N}.}                       \tag{ES.5c}
```

There is no hidden superlinear rare shell in this representation.

## 2. Exact shell/geometric chain rule

Let `bar mu` be any law absolutely continuous with respect to
`mu_epsilon`, including the average latent posterior induced by any actual
negative-disorder bridge law.  Write `bar p` for its `E_epsilon` marginal
and `bar mu_e` for its conditional law on `mathcal S_e`.

**Theorem ES.1 (radial/geometric retuning split).**  One has

```math
 \boxed{
 D(\bar\mu\Vert\mu_\epsilon)
 =D(\bar p\Vert p)
  +\sum_e\bar p_e
       D(\bar\mu_e\Vert U_{\mathcal S_e}).}            \tag{ES.6}
```

Equivalently,

```math
 D(\bar\mu_e\Vert U_{\mathcal S_e})
 =\log|\mathcal S_e|-H(\bar\mu_e).                    \tag{ES.7}
```

If the left side of (ES.6) is at least `cN`, then either

```math
 D(\bar p\Vert p)\ge {cN\over2}                       \tag{ES.8}
```

or

```math
 \sum_e\bar p_e
 \{\log|\mathcal S_e|-H(\bar\mu_e)\}
 \ge {cN\over2}.                                      \tag{ES.9}
```

*Proof.*  Apply the relative-entropy chain rule to the deterministic map
`Q -> E_epsilon(Q)` and use (ES.4).  Equation (ES.7) is the KL formula
against a uniform law.  The dichotomy follows by nonnegativity. `square`

The radial alternative has a quantitative rare-event form.

**Corollary ES.2 (linear shell KL exposes rare energy shells).**  If

```math
 D(\bar p\Vert p)\ge cN>0,                             \tag{ES.9a}
```

then, with

```math
 \mathcal R_c=\{e:p_e\le e^{-cN/2}\},                 \tag{ES.9b}
```

one has

```math
 \boxed{
 \bar p(\mathcal R_c)\ge{c\over2C_\beta-c},
 \qquad
 p(\mathcal R_c)
 \le O(N^2)e^{-cN/2}.}                                \tag{ES.9c}
```

Necessarily `c<=C_beta`, so the denominator is positive.

*Proof.*  Let `X(e)=-log p_e`.  Since

```math
 D(\bar p\Vert p)=E_{\bar p}X-H(\bar p),
```

(ES.9a) implies `E_bar p X>=cN`.  Equation (ES.5c) gives
`0<=X<=C_beta N`.  If `a=bar p(R_c)`, then

```math
 cN\le E_{\bar p}X
 \le aC_\beta N+(1-a){cN\over2},
```

which rearranges to the first bound.  The second follows from (ES.5) and
the definition of `R_c`. `square`

Applied after the coarse dichotomy (ES.8)--(ES.9), linear total retuning
`D(bar mu||mu_epsilon)>=cN` therefore yields either

```math
 \bar p\{e:p_e\le e^{-cN/4}\}
 \ge {c\over4C_\beta-c}                               \tag{ES.9d}
```

or the within-shell entropy deficit in (ES.9) is at least `cN/2`.

For the posterior in Theorem LE.1 there is no hidden fibre charge.  The
channel likelihood depends on the augmented latent state only through `Q`,
so both the pointwise and averaged posterior retain the prior conditional
law on every fibre over `Q`.  Therefore

```math
 D(\bar\nu\Vert\nu)
 =D(\bar\mu_Q\Vert\mu_\epsilon),                       \tag{ES.10}
```

and (ES.6) applies to the exact retuning term in (LE.10)--(LE.11).

There is also no hidden information capacity in the shell label.  Under the
joint posterior law `eta(dB,dQ)=q(dB)nu_B(dQ)`,

```math
 \boxed{
 E_qD(\nu_{E\mid B}\Vert p)
 =I_\eta(B;E)+D(\bar p\Vert p),
 \qquad
 I_\eta(B;E)\le\log|\operatorname {range}E|=O(\log N).} \tag{ES.10a}
```

Thus linear shell posterior work is coherent marginal retuning of the
one-dimensional energy law; the genuine information shared with the bridge
through that label is only logarithmic.

## 3. Interaction with the spread theorem

Theorem 37.61 gives

```math
 \max_Q\mu_\epsilon(Q)
 \le e^{-\eta_\beta N+o(N)}.                          \tag{ES.11}
```

Consequently a sublinear-retuning posterior cannot collapse onto a
subcritical-rate catalogue, whether that catalogue crosses many shells or
lies inside one shell.  Equations (ES.6)--(ES.9) sharpen the surviving
possibility:

1. **shell retuning:** the posterior changes the distribution of the
   scalar combined child energy by linear KL, hence puts fixed mass on an
   exponentially rare energy-shell event by Corollary ES.2; or
2. **within-shell retuning:** after conditioning on that scalar, it loses a
   linear amount of entropy relative to a uniform equal-energy shell.

The first branch is visible to an `O(N^2)`-state child statistic.  The
second is the genuinely nonradial, diffuse exponential-rate obstruction.
It is strictly narrower than the complete Gibbs table and than the earlier
phrase “diffuse retuning.”

This is not yet a row-lifetime theorem.  No proved inequality identifies
`J-I^leftarrow` with (or bounds it by) either term in (ES.6), and a shell
entropy deficit need not expose a favorable bridge direction.  The new SML
is to connect the canonical factor-retuning lifetime to the dichotomy
(ES.8)--(ES.9): control the polynomial shell branch directly, and either
compress or lower-bound the within-shell geometric branch using actual
child optimality.
