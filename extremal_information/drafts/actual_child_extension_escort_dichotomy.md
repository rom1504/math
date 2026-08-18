# Actual child rows as one-vertex extension escorts

Status: **rigorous, independently audited task-local theorem**.  This note concerns the
forward channel and the canonical inverse-row factors induced by the actual
contracted-temperature children.  It identifies those row factors with
one-vertex extension responses of a child, including the exact scalar
sector-bias correction.  Exact child minimality then upgrades the neutral
row from a generic subexponential-density bound to a uniform density bound.
The result sharpens the effective-support classification; it does not
control the remaining joint row interaction.

## 1. Augmented and biased child measures

For a signing `D` of order `n`, put

```math
Z_D^a(t)=E_y e^{atH_D(y)},\qquad a\in\{\pm1\},
```

```math
\overline Z_D(t)=E_y\cosh(tH_D(y))
 ={Z_D^+(t)+Z_D^-(t)\over2},
```

and define its sector bias

```math
\gamma_D(t)={1\over2}\log{Z_D^+(t)\over Z_D^-(t)}. \tag{EE.1}
```

For an arbitrary real `gamma`, let

```math
{d\overline\mu_{D,t}^{\gamma}\over dU_n}(y)
 ={\cosh(tH_D(y)+\gamma)
   \over E_{U_n}\cosh(tH_D+\gamma)}.                \tag{EE.2}
```

This is the spin marginal of the augmented child Gibbs measure after adding
a field `gamma` to the auxiliary sign.  It remains invariant under global
spin flip because the quadratic Hamiltonian satisfies `H_D(-y)=H_D(y)`.

For a bridge amplitude `u`, define the normalized biased extension response

```math
z_{D,t,u}^{\gamma}(b)
 ={E_{\overline\mu_{D,t}^{\gamma}}
       \cosh(u\langle b,Y\rangle)
   \over(\cosh u)^n},
\qquad b\in\{\pm1\}^n.                              \tag{EE.3}
```

It is a likelihood relative to the fair row law, so

```math
E_{U_n}z_{D,t,u}^{\gamma}=1.                        \tag{EE.4}
```

## 2. The exact row-extension identity

Let `A,D` be the left and right children in the actual forward channel, and
fix the relative orientation `epsilon`.  Proposition CR.0 writes the erased
row prior as

```math
\mu_{\rm row}^{(\epsilon)}
=\sum_{s=\pm1}\pi_s^{(\epsilon)}\mu_{D,\epsilon s},
\qquad
\pi_s^{(\epsilon)}\propto Z_A^s(t)Z_D^{\epsilon s}(t).
                                                               \tag{EE.5}
```

**Theorem EE.1 (sector rows are biased extension channels).**  The mixture
in (EE.5) is exactly

```math
\boxed{
\mu_{\rm row}^{(\epsilon)}
=\overline\mu_{D,t}^{\epsilon\gamma_A(t)}.}         \tag{EE.6}
```

Consequently its output likelihood is

```math
\boxed{
p_{{\rm row},u}^{(\epsilon)}(b)
=z_{D,t,u}^{\epsilon\gamma_A(t)}(b).}               \tag{EE.7}
```

When the forward law `Pi` is averaged over `epsilon` (rather than when an
inverse escort is averaged), the erased row prior is the unbiased augmented
child measure and its likelihood is `z_{D,t,u}^0`.

At equal internal and incident amplitudes, the unbiased response is exactly
a one-vertex extension ratio.  If `D\oplus b` denotes the order-`n+1`
signing obtained by adjoining one vertex with incident row `b`, then

```math
\boxed{
z_{D,t,t}^0(b)
={\overline Z_{n+1}(D\oplus b,t)
  \over \overline Z_n(D,t)(\cosh t)^n}.}            \tag{EE.8}
```

*Proof.*  In (EE.5), multiplication of each sector Gibbs law by its weight
cancels `Z_D^{epsilon s}`.  Up to normalization, the resulting density at
`y` is

```math
Z_A^+(t)e^{\epsilon tH_D(y)}
 +Z_A^-(t)e^{-\epsilon tH_D(y)}
=2\sqrt{Z_A^+Z_A^-}
 \cosh(tH_D(y)+\epsilon\gamma_A),                  \tag{EE.9}
```

where evenness of `cosh` was used in the last equality.  This proves
(EE.6), and the binary channel formula proves (EE.7).

For (EE.8), average the new spin `x_0` first:

```math
E_{x_0}\cosh\{tH_D(y)+tx_0\langle b,y\rangle\}
=\cosh(tH_D(y))\cosh(t\langle b,y\rangle).         \tag{EE.10}
```

Divide the resulting partition by `Zbar_D(t)(cosh t)^n`. `square`

Equation (EE.8) is a state reduction specific to the actual channel: the
canonical row factor is the negative escort of the complete one-vertex
extension response, not an arbitrary bounded-density row law.

## 3. Exact minimization forces diffuse canonical rows

Let

```math
F_k(t)=\min_C\log\overline Z_k(C,t)
```

and define the adjacent annealed deficit

```math
\delta_n(t)
=n\log\cosh t-\{F_{n+1}(t)-F_n(t)\}.               \tag{EE.11}
```

The raw-temperature cavity inequalities give

```math
0\le\delta_n(t)\le n\log\cosh t.                  \tag{EE.12}
```

For `lambda>0`, let the canonical inverse row escort be

```math
{dr_{D,t,u}^{\gamma}\over dU_n}(b)
 ={(z_{D,t,u}^{\gamma}(b))^{-\lambda}
   \over E_{U_n}(z_{D,t,u}^{\gamma})^{-\lambda}}.  \tag{EE.13}
```

**Theorem EE.2 (optimizer extension-support bound).**  Suppose `D` is an
exact minimizer defining `F_n(t)` and take `u=t`.  Then for every `b`,

```math
z_{D,t,t}^0(b)\ge e^{-\delta_n(t)}.                 \tag{EE.14}
```

For every real `gamma`,

```math
z_{D,t,t}^{\gamma}(b)
\ge e^{-\delta_n(t)-2|\gamma|},                    \tag{EE.15}
```

and hence

```math
\boxed{
D_\infty(r_{D,t,t}^{\gamma}\Vert U_n)
\le\lambda\{\delta_n(t)+2|\gamma|\}.}            \tag{EE.16}
```

For comparison, the one-bit channel oscillation gives the universal bound

```math
D_\infty(r_{D,t,u}^{\gamma}\Vert U_n)
\le2\lambda un.                                    \tag{EE.16a}
```

Thus at `u=t` the sharp combined statement proved here is

```math
\boxed{
D_\infty(r_{D,t,t}^{\gamma}\Vert U_n)
\le\lambda\min\{\delta_n(t)+2|\gamma|,\,2tn\}.}   \tag{EE.16b}
```

Bound (EE.16a) is only `O(sqrt(N))` at physical scaling and is generic; the
`O(1)` neutral bound from exact minimization is the new part.

Equivalently, its min-entropy obeys

```math
\boxed{
H_\infty(r_{D,t,t}^{\gamma})
\ge n\log2-\lambda\{\delta_n(t)+2|\gamma|\}.}    \tag{EE.17}
```

*Proof.*  Every `D\oplus b` is an admissible order-`n+1` signing, so exact
minimality and (EE.8) give

```math
\log z_{D,t,t}^0(b)
\ge F_{n+1}(t)-F_n(t)-n\log\cosh t=-\delta_n(t).
                                                               \tag{EE.18}
```

For every real `a,g`,

```math
e^{-|g|}\cosh a\le\cosh(a+g)\le e^{|g|}\cosh a.  \tag{EE.19}
```

Applying the same bounds to the normalizing constants in (EE.2) gives

```math
e^{-2|\gamma|}
\le {d\overline\mu_{D,t}^{\gamma}
          \over d\overline\mu_{D,t}^{0}}
\le e^{2|\gamma|}.                                 \tag{EE.20}
```

Equations (EE.3), (EE.14), and (EE.20) prove (EE.15).  By (EE.4) and
convexity, the denominator in (EE.13) is at least one.  Thus its density is
at most the right side of (EE.16) after exponentiation.  This proves
(EE.16)--(EE.17).  Finally, flipping one output bit changes `log z` by at
most `2u`.  Hence its range is at most `2un`; since its uniform mean is one,
`min z>=e^(-2un)`.  The same escort argument proves (EE.16a). `square`

At the contracted scale `t=beta/sqrt(N)`, with `n<=N`, (EE.12) gives

```math
\delta_n(t)\le {\beta^2n\over2N}=O_\beta(1).       \tag{EE.21}
```

Thus the unbiased canonical row has a uniform `L^infinity` density bound,
strictly stronger than the universal Renyi-two bound.  A sector row has the
same conclusion whenever the opposite child's sector bias is bounded.

The same proof also bounds the one-row inverse work itself:

```math
0\le {1\over\lambda}\log E_{U_n}
       (z_{D,t,t}^{\gamma})^{-\lambda}
\le\delta_n(t)+2|\gamma|.                           \tag{EE.21a}
```

For `gamma=0`, the deficits have an exact Bellman and telescoping
interpretation.

**Proposition EE.3 (one-vertex Bellman identity).**  At every fixed raw
temperature `t`,

```math
\boxed{
F_{n+1}(t)=\min_D\left\{
 \log\overline Z_n(D,t)+n\log\cosh t
 +\min_b\log z_{D,t,t}^0(b)\right\}.}              \tag{EE.21b}
```

If `C` is an exact order-`n+1` minimizer, then deleting any vertex gives a
pair `(D,b)` attaining the outer and inner minimum jointly; in particular,
its incident row minimizes `z_D^0` for that deletion.  Moreover,

```math
\boxed{
\sum_{k=1}^{n-1}\delta_k(t)
={n\choose2}\log\cosh t-F_n(t).}                   \tag{EE.21c}
```

*Proof.*  Every order-`n+1` signing is uniquely a pair `(D,b)` after a
vertex is distinguished, and (EE.8) gives its pressure, proving (EE.21b).
If an incident row of a minimizing `C` did not attain the inner minimum for
its deletion, replacing that row would lower the pressure.  Finally sum
(EE.11) over `k`; the pressure increments telescope and
`sum_(k=1)^(n-1)k=binom n2`, while `F_1(t)=0`. `square`

Thus the canonical row law is the entropic relaxation of the exact
one-vertex Bellman response, and `delta_n` is the adjacent contribution to
the full annealed-minus-optimized pressure budget.  This does not make
`delta_n` small, but it identifies it without a bridge or external-field
oracle.

## 4. Orientation or diffuse support

The sector biases also determine the forward orientation exactly.  Under
the two independent child Gibbs laws,

```math
E\tau_1=\tanh\gamma_A,
\qquad E\tau_2=\tanh\gamma_D,
```

so

```math
\Pi(\epsilon)
={1+\epsilon\tanh\gamma_A\tanh\gamma_D\over2}.     \tag{EE.22}
```

Its reverse entropy cost is

```math
\boxed{
D(U_\epsilon\Vert\Pi_\epsilon)
=-{1\over2}\log\{1-\tanh^2\gamma_A
                       \tanh^2\gamma_D\}.}         \tag{EE.23}
```

This is an exact nonnegative summand of the original joint compensation:

```math
D(U_\epsilon U_B\Vert\Pi_{\epsilon,B})
=D(U_\epsilon\Vert\Pi_\epsilon)
 +E_{U_\epsilon}D(U_B\Vert\Pi_{B\mid\epsilon}).   \tag{EE.23a}
```

Put `g=min(|gamma_A|,|gamma_D|)`.  Since
`1-tanh^4 g<=8e^{-2g}`,

```math
D(U_\epsilon\Vert\Pi_\epsilon)
\ge g-{3\over2}\log2.                              \tag{EE.24}
```

Combining this with EE.2 in the row direction whose opposite-child bias is
`g` proves the following actual-child alternative.

**Corollary EE.4 (orientation-or-diffuse-row dichotomy).**  For exact
minimizing children at `t=beta/sqrt(N)` and every threshold `G>=0`, one of
the following alternatives holds:

1. `g>=G`, in which case the orientation alone pays

```math
D(U_\epsilon\Vert\Pi_\epsilon)
\ge G-{3\over2}\log2;                              \tag{EE.24a}
```

2. `g<=G`, in which case one may choose a row direction for which every
canonical sector row factor satisfies

```math
H_\infty(r_{\rm row})
\ge k\log2-\lambda\{O_\beta(1)+2g\},               \tag{EE.25}
```

where `k` is the opposite shore size.  In particular, if
`D(U_epsilon||Pi_epsilon)=o(N)` along comparable splits, then `g=o(N)` and
the chosen canonical row factor has full exponential effective support,

```math
H_\infty(r_{\rm row})=k\log2-o(N).                 \tag{EE.26}
```

*Proof.*  Formula (EE.22) follows from independence and
`epsilon=tau_1tau_2`; direct evaluation gives (EE.23).  Monotonicity of
`tanh` and the displayed elementary estimate give (EE.24).  Choose the
shore whose role as the left child has bias `g`, use (EE.7), and apply
(EE.17)--(EE.21). `square`

The alternative is invariant under transposing the bridge.  It concerns
the actual induced law and uses only adjacent optimal pressures and the two
scalar sector biases.  Neither quantity contains the full external-field
or bridge landscape.

## 5. What remains

EE.2 gives a bounded-density theorem for the neutral actual extension
escort, while EE.4 shows that an extensive sector-symmetry obstruction is
already visible in the one-bit orientation channel.  Together with ES.0,
the actual law has exponentially diffuse canonical rows and uniformly
regular conditional row changes; in the sector-neutral regime the row
diffuseness is uniform rather than merely exponential-scale.

This does **not** bound the weighted path total correlation or collective
marginal retuning in ES.28.  The smallest missing lemma becomes:

> For actual minimizing children with sublinear orientation cost, use the
> one-vertex-extension structure (EE.7)--(EE.18), rather than the full
> external-field table, to prove that the weighted joint row resource in
> ES.28 is `o(N)`, or prove that a positive density of the diffuse extension
> rows carries irreducible mutual information/coherent retuning.

The extension identity makes this question strictly narrower than arbitrary
bounded-Renyi row regularity: all one-row factors now lie in a specified
all-order response family with the optimizer lower envelope (EE.14).

The independent adversarial audit is
[`../audits/actual_child_extension_escort_dichotomy_adversarial_audit.md`](../audits/actual_child_extension_escort_dichotomy_adversarial_audit.md).
