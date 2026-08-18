# Sector-bias balancing for actual optimizing children

Status: **rigorous task-local theorem**.  This note sharpens the biased
one-vertex extension bound in EE.2.  The scalar bias of either child need not
be bounded.  Nevertheless, after choosing the bridge orientation and the row
filtration direction together, every canonical inverse-row factor has a
dimension-free `D_infinity` bound at the contracted physical scale.  This
eliminates sector-bias escape from the canonical-component support branch in
that selected presentation.  It does not control the joint row interaction or
prove that the full bridge escort is close to a product.

## 1. Exact sector-weight transfer

Let `D` be a signing of order `n`.  Use the sector partitions and sector
output likelihoods

```math
Z_D^a(t)=E_y e^{atH_D(y)},qquad
z_D^a(b;t,u)
={E_{\mu_{D,a,t}}\cosh(u\langle b,Y\rangle)
  \over(\cosh u)^n},qquad a\in\{\pm1\}.
```

Each `z_D^a` is positive and has uniform mean one.  Put

```math
\gamma_D={1\over2}\log{Z_D^+\over Z_D^-}.
```

The unbiased augmented extension response and the response in an external
sector field `g` are respectively

```math
z_D^0
=\sum_{a=\pm1}{e^{a\gamma_D}\over2\cosh\gamma_D}z_D^a,
```

```math
z_D^g
=\sum_{a=\pm1}{e^{a(\gamma_D+g)}
                  \over2\cosh(\gamma_D+g)}z_D^a.       \tag{SB.1}
```

The superscript on `z_D^g` denotes the external field, as in EE.2; the total
sector log-odds in the second mixture is `gamma_D+g`.

**Lemma SB.1 (sharp sector-weight comparison).**  Pointwise on the row cube,

```math
\boxed{
z_D^g(b;t,u)\ge \kappa(\gamma_D,g)z_D^0(b;t,u),
\qquad
\kappa(\gamma_D,g)
=e^{-|g|}{\cosh\gamma_D\over\cosh(\gamma_D+g)}.}       \tag{SB.2}
```

Equivalently, with

```math
\Delta(\gamma_D,g)
=|g|+\log\cosh(\gamma_D+g)-\log\cosh\gamma_D,          \tag{SB.3}
```

one has `z_D^g>=e^(-Delta)z_D^0`.  The penalty is nonnegative.

*Proof.*  Let `w_a` and `q_a` be the two sector weights in (SB.1).  Since all
sector likelihoods are positive,

```math
\sum_aq_az_D^a
\ge\left(\min_a{q_a\over w_a}\right)\sum_aw_az_D^a.
```

Direct calculation gives

```math
\min_a{q_a\over w_a}
={\cosh\gamma_D\over\cosh(\gamma_D+g)}
  \min_a e^{ag}
=e^{-|g|}{\cosh\gamma_D\over\cosh(\gamma_D+g)}.
```

The minimum of two positive numbers whose `w`-weighted mean is one is at
most one, so `Delta>=0`. `square`

The earlier comparison `z_D^g>=e^(-2|g|)z_D^0` follows because the logarithm
of `cosh` is one-Lipschitz.  Formula (SB.2) is strictly sharper when the two
child biases cancel.

## 2. Balancing two actual children

Let `C,D` be the two children, with orders `m,n`, at the common raw
temperature `t`.  In relative orientation `epsilon`, EE.1 identifies the
canonical erased-row likelihood in the direction `C -> D` as

```math
z_{C\to D}^{\epsilon}=z_D^{\epsilon\gamma_C}.          \tag{SB.4}
```

The bridge can be transposed, so either child may be used as the base `D`.
It can also use either relative orientation.

**Theorem SB.2 (orientation--filtration balancing).**  Label the children so
that

```math
|\gamma_C|\le|\gamma_D|,                              \tag{SB.5}
```

and choose `epsilon` so that `gamma_D` and `epsilon gamma_C` have opposite
signs (either choice works when one bias is zero).  Then

```math
\boxed{
z_{C\to D}^{\epsilon}(b;t,u)
\ge {1\over2}z_D^0(b;t,u)
\quad\hbox{for every }b.}                            \tag{SB.6}
```

Thus the conclusion is uniform in the magnitudes of both sector biases.

*Proof.*  Put `a=|gamma_C|` and `d=|gamma_D|`.  After the chosen orientation,
the two arguments in SB.2 have magnitudes `d` and `d-a`, with `0<=a<=d`.
Consequently

```math
\kappa=e^{-a}{\cosh d\over\cosh(d-a)}.
```

The addition formula gives

```math
\cosh d
=\cosh(d-a)\cosh a+\sinh(d-a)\sinh a
\ge\cosh(d-a)\cosh a.
```

Therefore `kappa>=e^(-a)cosh a=(1+e^(-2a))/2>=1/2`, and
SB.1 proves (SB.6). `square`

The chosen orientation is not an artificial sacrifice.  Before adding the
bridge, independence of the two child spin blocks gives the exact identity

```math
\overline Z_{C\oplus\epsilon D}(t)
=\sqrt{Z_C^+Z_C^-Z_D^+Z_D^-}\,
  \cosh(\gamma_C+\epsilon\gamma_D).                 \tag{SB.6a}
```

Thus the same cancellation choice minimizes the zero-bridge augmented
pressure over the two relative orientations.

For identical children this simply chooses opposite orientation.  Their
total sector field is then zero, and the equal-sector response is at least
half the neutral augmented response even when the common bias is arbitrarily
large.

## 3. Dimension-free inverse component complexity

Let

```math
F_k(t)=\min_A\log\overline Z_k(A,t),
\qquad
\delta_n(t)
=n\log\cosh t-\{F_{n+1}(t)-F_n(t)\}.
```

Suppose now that the base child `D` is an exact order-`n` pressure minimizer.
EE.2 proves the optimizer-specific lower envelope

```math
z_D^0(b;t,t)\ge e^{-\delta_n(t)}.                    \tag{SB.7}
```

For `lambda>0`, let `r_{C->D}^epsilon` be the normalized canonical inverse
row escort with density proportional to
`(z_{C->D}^epsilon)^(-lambda)` relative to the fair row law.

**Corollary SB.3 (tight actual-child row factors without a bias bound).**
Under the choice in Theorem SB.2 and at `u=t`, every row satisfies

```math
\boxed{
D_\infty(r_{C\to D}^{\epsilon}\Vert U_n)
\le\lambda\{\delta_n(t)+\log2\}.}                   \tag{SB.8}
```

The same right side bounds every finite Renyi divergence and gives

```math
H_\infty(r_{C\to D}^{\epsilon})
\ge n\log2-\lambda\{\delta_n(t)+\log2\}.            \tag{SB.9}
```

At the contracted scale `t=beta/sqrt(N)`, where `m+n=N`,

```math
\boxed{
D_\infty(r_{C\to D}^{\epsilon}\Vert U_n)
\le\lambda\left({\beta^2n\over2N}+\log2\right)
\le\lambda\left({\beta^2\over2}+\log2\right).}      \tag{SB.10}
```

*Proof.*  Equations (SB.6)--(SB.7) give the pointwise lower bound
`z_(C->D)^epsilon>=exp(-delta_n)/2`.  Its uniform mean is one.  Jensen gives
`E_U z^(-lambda)>=1`, so the normalized inverse density is at most
`exp(lambda(delta_n+log 2))`.  This proves (SB.8), and the Renyi and
min-entropy statements follow.  Finally use
`delta_n(t)<=n log cosh t<=nt^2/2`. `square`

The constant `log 2` is the sharp universal price of comparing an
equal-sector mixture with an arbitrarily biased neutral mixture: as
`a=d -> infinity`, the coefficient in SB.2 tends to `1/2`.  This is only a
sharpness statement for the positive-mixture algebra, not a claim that such
an extreme is attained by a sequence of optimizing signings.

## 4. Implication and remaining obstruction

The theorem gives a complete answer to the **sector-bias component-support**
question for a bridge construction allowed to choose its relative
orientation and row direction:

```math
\boxed{
\text{arbitrary child sector biases}
+\text{ exact pressure minimality}
\Longrightarrow
\text{a canonical row presentation with tight component }D_\infty.}
```

It is stronger than the generic bounded-oscillation `D_2` estimate and than
the biased EE.2 bound, because its constant is independent of both
`gamma_C` and `gamma_D`.  It uses only one scalar comparison to choose among
the two orientations and the two transpose filtrations; it does not store a
row response table or reconstruct the parent optimization.

There are two important scope limits.

1. The statement chooses one relative orientation.  It does not assert the
   same bound for an orientation fixed adversarially in advance.  For the
   parent upper construction this choice is available and minimizes the
   zero-bridge pressure by (SB.6a), but no theorem says that it remains the
   target-reaching orientation after bridge optimization.  An argument
   requiring simultaneous control of both sectors needs more.
2. Tight canonical row factors do not control the interaction remainder
   `h_epsilon`, row total correlation, or the reverse product projection
   `mathcal I_lambda^leftarrow`.  Those joint resources can remain linear.

Accordingly, a standalone bound `gamma_C=O(1)` is no longer needed to rule
out escaping **canonical-product components in the selected presentation**.
The actual missing statement is joint: prove that the balanced presentation
has `o(N)` row dependence and retuning, or exhibit a positive density of
irreducible dependence in that same balanced presentation.

The independent adversarial audit is
[`../audits/actual_child_sector_bias_balancing_adversarial_audit.md`](../audits/actual_child_sector_bias_balancing_adversarial_audit.md).
