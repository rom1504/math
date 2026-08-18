# Adversarial audit of the extension-escort dichotomy

**Object audited.**
[`../drafts/actual_child_extension_escort_dichotomy.md`](../drafts/actual_child_extension_escort_dichotomy.md),
especially EE.6--EE.8, EE.14--EE.17, and EE.22--EE.26.

**Verdict: mathematical core PASS; scope needs correction.**  The sector
algebra, extension normalization, optimizer lower envelope, `D_infty`
normalization, orientation formula, and transpose choice are correct.  The
draft overstates which support conclusion is optimizer-specific.  A generic
one-bit likelihood-ratio bound already gives

```math
D_\infty(r_{\rm row}\Vert U_k)\le2\lambda u k,       \tag{EEA.1}
```

so at `u=beta/sqrt(N)` every such row has
`H_infty=k log 2-o(N)`, without minimality or an orientation dichotomy.  The
genuine optimizer-specific gain in EE.2 is the sharper `O(1)` unbiased bound
and the biased bound `O(1+|gamma|)`.

This audit does not modify the theorem draft.

## 1. Sector algebra in EE.6

Conditional on `epsilon`, CR.0 gives sector weights

```math
\pi_s^{(\epsilon)}
={Z_A^sZ_D^{\epsilon s}\over
  \sum_cZ_A^cZ_D^{\epsilon c}}
```

and sector densities

```math
{d\mu_{D,\epsilon s}\over dU_n}(y)
={e^{\epsilon stH_D(y)}\over Z_D^{\epsilon s}}.
```

Multiplication cancels the right-child sector partition exactly.  The
unnormalized mixture density is therefore

```math
Z_A^+e^{\epsilon tH_D(y)}
+Z_A^-e^{-\epsilon tH_D(y)}.                        \tag{EEA.2}
```

With `gamma_A=(1/2)log(Z_A^+/Z_A^-)`, this equals

```math
2\sqrt{Z_A^+Z_A^-}
\cosh(tH_D(y)+\epsilon\gamma_A).                    \tag{EEA.3}
```

For `epsilon=-1`, the two exponential coefficients swap exactly as required;
there is no sign error.  After normalization, (EEA.3) is EE.2 with parameter
`epsilon gamma_A`, proving EE.6.

Every sector law is invariant under `Y -> -Y`.  Hence

```math
E e^{u\langle b,Y\rangle}
=E\cosh(u\langle b,Y\rangle),
```

and the binary-channel factor

```math
\prod_j(1+\tanh(u)b_jY_j)
={e^{u\langle b,Y\rangle}\over(\cosh u)^n}
```

proves EE.7 with the advertised normalization.

The sentence before EE.8 saying that the unconditioned row prior is unbiased
is correct only when `epsilon` is averaged with its **forward-law marginal**
`Pi(epsilon)`.  It is not true under an arbitrary, uniform, or
negative-escort reweighting of orientation.  The draft should make this
conditioning convention explicit.

## 2. EE.8 and the exact-minimizer envelope

For a new spin `x_0`, direct averaging gives

```math
\begin{aligned}
\overline Z_{D\oplus b}(t)
 &=E_y\cosh(tH_D(y))\cosh(t\langle b,y\rangle)\\
 &=\overline Z_D(t)(\cosh t)^n z_{D,t,t}^0(b),       \tag{EEA.4}
\end{aligned}
```

so EE.8 is exact.  It holds for every signing `D`; neither the identity nor
the interpretation as an extension response is specific to an optimizer or
even to the actual-child selection.  Optimizer-specific content enters only
at EE.14.

If `D` attains `F_n(t)`, then every `D oplus b` is an admissible order-`n+1`
signing.  Consequently

```math
\begin{aligned}
\log z_{D,t,t}^0(b)
 &\ge F_{n+1}(t)-F_n(t)-n\log\cosh t\\
 &=-\delta_n(t),                                    \tag{EEA.5}
\end{aligned}
```

which proves EE.14 for **every** row `b`.  No assumption that `D` is the
deletion of an order-`n+1` optimizer is needed.  Conversely, this does not
say that `D` is a Bellman-optimal base for order `n+1`; the inequality may
have slack for all `b`.

The adjacent bounds `0<=delta_n<=n log cosh t` follow from the exact Bellman
recurrence: extension cannot reduce normalized augmented pressure, while a
uniformly averaged row has response `(cosh t)^n`.  Thus EE.12 and the
physical-scale estimate EE.21 are also correct.

## 3. Biased envelope and `D_infty` normalization

For all real `a,g`,

```math
e^{-|g|}\cosh a\le\cosh(a+g)\le e^{|g|}\cosh a.
```

The same comparison for the partition normalizers gives

```math
e^{-2|g|}
\le {d\bar\mu^g\over d\bar\mu^0}
\le e^{2|g|}.                                       \tag{EEA.6}
```

Since the observable `cosh(u<b,Y>)` is nonnegative, EE.15 follows from the
lower side of (EEA.6) and EE.14.

Let `c=delta_n+2|gamma|`.  Then `z^gamma>=e^{-c}` and
`E_U z^gamma=1`.  For

```math
f={dr\over dU}={z^{-\lambda}\over E_Uz^{-\lambda}},
```

convexity of `x -> x^{-lambda}` gives `E_Uz^{-lambda}>=1`, while the pointwise
lower envelope gives `z^{-lambda}<=e^{lambda c}`.  Therefore

```math
\|f\|_\infty\le e^{\lambda c},
\qquad
D_\infty(r\Vert U)\le\lambda c.                    \tag{EEA.7}
```

There is no omitted log-normalizer or additional factor of `lambda`.
Multiplying the density bound by the uniform atom mass `2^{-n}` proves the
min-entropy normalization in EE.17.  The one-row inverse-work bound EE.21a
uses the same two inequalities and is correct.

## 4. Orientation and transpose audit

Under a child's augmented Gibbs law,

```math
P(\tau=a)={Z^a\over Z^++Z^-},
\qquad E\tau=\tanh\gamma.
```

The two child variables are independent before the bridge channel, and the
channel integrates to one.  For `epsilon=tau_1tau_2`, this gives

```math
\Pi(\epsilon)
={1+\epsilon\tanh\gamma_A\tanh\gamma_D\over2},      \tag{EEA.8}
```

exactly EE.22.  If `q=tanh gamma_A tanh gamma_D`, direct binary calculation
gives

```math
D(U_\epsilon\Vert\Pi_\epsilon)
=-{1\over2}\log(1-q^2),                             \tag{EEA.9}
```

so the direction of KL in EE.23 is correct.

With `g=min(|gamma_A|,|gamma_D|)`, one has

```math
1-q^2\le1-\tanh^4g
\le2\operatorname{sech}^2g
\le8e^{-2g},                                        \tag{EEA.10}
```

which proves the constant `g-(3/2)log 2` in EE.24.

The transpose choice in EE.4 is also correct, but merits explicit wording:

- if `|gamma_A|=g`, retain rows indexed by the left child `A`; each row has
  length `|D|` and bias magnitude `g`;
- if `|gamma_D|=g`, transpose the bridge; rows are indexed by `D`, have
  length `|A|`, and bias magnitude `g`.

EE.2 is then applied to the child on which that row lives.  This requires
both children to be exact pressure minimizers at the same raw temperature,
as assumed.

## 5. Generic coordinate-oscillation ceiling

The main scope issue is that “full exponential effective support” is not new
to EE.2.  For every binary-channel row likelihood—including arbitrary
sector bias—flipping one output bit changes `log z` by at most `2u`.
Therefore the inverse escort density `f=z^{-lambda}/E z^{-lambda}` obeys

```math
|\log f(b)-\log f(b^{(j)})|\le2\lambda u.            \tag{EEA.11}
```

Any two cube points are at Hamming distance at most `k`, so
`osc(log f)<=2 lambda u k`.  Since `E_U f=1`, its maximum is at most
`exp(osc(log f))`, proving (EEA.1).  Equivalently,

```math
H_\infty(r_{\rm row})
\ge k\log2-2\lambda uk.                             \tag{EEA.12}
```

At `u=beta/sqrt(N)` and `k<=N`, the deficit is `O(sqrt(N))=o(N)` for every
row law, with no child optimality, sector-bias bound, or orientation-cost
assumption.  Hence:

1. EE.26 is true but already generic.
2. EE.4 does not newly rule out an exponentially sparse row support.
3. The statements in the status, Section 5, and the final SML should not use
   mere `H_infty=k log2-o(N)` as optimizer-specific progress.

The valid quantitative synthesis is the combined bound

```math
\boxed{
D_\infty(r_{D,t,t}^\gamma\Vert U_n)
\le\lambda\min\{\delta_n(t)+2|\gamma|,\,2tn\}.}    \tag{EEA.13}
```

The first branch is optimizer-specific; the second is generic.  At physical
scale, the unbiased case `gamma=0` improves `O(sqrt(N))` to `O_beta(1)`.
For a biased row, EE.2 improves the generic bound only when
`|gamma|=o(sqrt(N))` (or quantitatively below the generic ceiling), not merely
when orientation cost is `o(N)`.

## 6. Recommended scope corrections

The theorem draft should retain EE.1--EE.3 and the exact algebra of EE.4,
but revise the interpretation as follows:

- call EE.6--EE.8 generic channel/extension identities specialized to the
  child notation, not optimizer reductions;
- specify that “before conditioning on epsilon” means averaging with the
  forward marginal `Pi(epsilon)`;
- present EE.14 and its `O(1)` unbiased `D_infty` consequence as the genuine
  optimizer theorem;
- replace the claim that EE.4 closes sparse support by the sharper bound
  (EEA.13);
- state that sublinear orientation cost only makes the optimizer envelope
  subexponential; generic oscillation already does that, often more sharply;
- qualify “canonical row law is the Bellman relaxation” as the **unbiased**
  row law.  The conditioned canonical law has bias
  `epsilon gamma_A` and is not the neutral Bellman response unless the other
  child is sector-neutral.

After those corrections, the exact substantive advance is narrow but real:
an exact minimizing child has a uniform lower envelope on every neutral
one-vertex extension, producing a dimension-free inverse-escort density
bound.  No conclusion about collective row interaction follows, and the
support-size conclusion alone does not reset the actual-child SML.
