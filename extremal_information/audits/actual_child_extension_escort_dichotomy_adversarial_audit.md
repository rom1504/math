# Adversarial audit: actual-child extension escort dichotomy

**Object audited:**
[`../drafts/actual_child_extension_escort_dichotomy.md`](../drafts/actual_child_extension_escort_dichotomy.md).

**Verdict:** **PASS.**  The sector algebra in EE.6, extension ratio EE.8,
escort normalization and `D_infinity` constant in EE.16, and the
orientation/transpose alternative EE.23--EE.25 are correct.  No hidden use
of a parent optimizer or full bridge-response oracle occurs.

## 1. Sector mixture

Multiplying the sector density `e^(epsilon s tH_D)/Z_D^(epsilon s)` by its
weight `Z_A^s Z_D^(epsilon s)` cancels the right-child partition exactly.
Writing `Z_A^+=g e^(gamma_A)` and `Z_A^-=g e^(-gamma_A)` leaves

```math
Z_A^+e^{\epsilon tH_D}+Z_A^-e^{-\epsilon tH_D}
=2g\cosh(tH_D+\epsilon\gamma_A).
```

The same factor normalizes after averaging in `y`, so this is an equality
of probability measures, not merely proportional unnormalized densities.
This verifies EE.6 for both signs of `epsilon`.  Averaging before
conditioning on the relative orientation similarly cancels the opposite
child and gives the neutral `cosh(tH_D)` marginal.

## 2. One-vertex extension ratio

For the extension `D\oplus b`, averaging the new spin first gives exactly

```math
E_{x_0}\cosh\{tH_D(y)+tx_0\langle b,y\rangle\}
=\cosh(tH_D(y))\cosh(t\langle b,y\rangle).
```

After averaging `y`, division by `Zbar_D` produces the expectation under
`\overline\mu_(D,t)^0`; the channel likelihood contains the additional declared
normalization `(cosh t)^n`.  Thus EE.8 has neither a missing factor two nor
an inverted extension ratio.

The adjacent deficit signs are consistent: averaging extension rows gives
an increment at most `n log cosh t`, while deleting a vertex and using
`cosh>=1` gives a nonnegative optimal-pressure increment.  Hence
`0<=delta_n<=n log cosh t`, and exact order-`n+1` minimality yields the
pointwise lower bound `z^0>=exp(-delta_n)`.

## 3. Biased comparison and escort normalization

The elementary bounds on `cosh(a+gamma)/cosh(a)`, applied once to the
density numerator and once to its normalizer, give the factor
`exp(+-2|gamma|)` in EE.20.  Therefore

```math
z^\gamma(b)\ge e^{-\delta_n-2|\gamma|}.
```

If `a=delta_n+2|gamma|`, the inverse escort has density

```math
{z(b)^{-\lambda}\over E_Uz^{-\lambda}}.
```

The numerator is at most `e^(lambda a)`.  Because `E_U z=1` and
`x mapsto x^(-lambda)` is convex for positive `lambda`, its denominator is
at least one.  Hence the essential supremum of the density is at most
`e^(lambda a)`, proving EE.16 with exactly the displayed constant and KL
direction.  The min-entropy conversion EE.17 follows because
`r(b)=2^(-n)(dr/dU)(b)`.

The universal comparison EE.16a is also normalized correctly: a bit flip
changes `log z` by at most `2u`, so the total range is at most `2un`; a
positive function of uniform mean one then has minimum at least
`e^(-2un)`.

## 4. Orientation cost

For a child auxiliary sign,

```math
E\tau={Z^+-Z^-\over Z^++Z^-}=\tanh\gamma.
```

Independence before conditioning and `epsilon=tau_1 tau_2` therefore give
EE.22.  Direct substitution into reverse KL from the fair orientation law
gives

```math
D(U_\epsilon\Vert\Pi_\epsilon)
=-{1\over2}\log(1-r^2),
\qquad r=\tanh\gamma_A\tanh\gamma_D,
```

which is EE.23.  The chain rule EE.23a is in the correct direction because
the reference first law has conditional bridge law `U_B` for each
orientation.

If `g=min(|gamma_A|,|gamma_D|)`, monotonicity gives
`r^2>=tanh^4 g`.  Also

```math
1-\tanh^4g
=(1-\tanh^2g)(1+\tanh^2g)
\le8e^{-2g},
```

so EE.24 follows with `3 log(2)/2` exactly.

Finally, if the smaller bias is on child `A`, make `A` the row-index shore;
EE.7 then biases the opposite child's extension response by `gamma_A=g`.
If the smaller bias is on `D`, transpose the bridge and reverse the roles.
The row length `k` is in each case the opposite shore size, so EE.25 has the
correct dimension.  Along comparable splits, sublinear orientation cost
forces `g=o(N)` by EE.24, and the error in EE.25 is `o(N)` while
`k=Theta(N)`, proving EE.26.

## 5. Scope

The neutral `D_infinity` bound is a genuine optimizer-specific strengthening:
at physical scaling `delta_n=O_beta(1)`.  Sector rows inherit it unless the
scalar sector bias is large, in which case the orientation channel already
records that resource.  The theorem correctly stops short of controlling
the ES.28 total-correlation or collective-retuning masses; bounded density
of each canonical row does not imply row independence.
