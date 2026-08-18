# Adversarial audit: actual-child flip averaging and radial ceiling

**Object audited:**
[`drafts/actual_child_flip_averaging_ceiling.md`](../drafts/actual_child_flip_averaging_ceiling.md)

**Verdict:** **PASS, with one target-regime qualification and one SML
qualification.**  FC.7--FC.21 are exact.  The order-eight calculations and
the `t>=3` pressure-minimizer claim check out.  The witness is a genuine
ceiling for universal radial information, but it is not an asymptotic
contracted-temperature witness.  Nothing in the note controls either term
of AC.24 for actual children, and the note mostly says so.

## 1. FC.7--FC.8 and equivalence

For `y in {+-1}` and `tanh s=(1-2p)tanh t`, direct expansion gives

```math
(1-p)e^{ty}+pe^{-ty}={\cosh t\over\cosh s}e^{sy}.
```

Multiplying this coordinate identity under the child Gibbs law proves all
three expressions in FC.7, including the `P_A(r)/P_A(rho 1)` ratio.  The
inequality in FC.8 follows by averaging the vertex inequalities FC.5.

The converse is also exact.  At a corner `s_e in {+-t}`, `Z_A(s)` is the
partition function of the corresponding sign flip and the cosh normalizer
is unchanged.  Equivalently, a multiaffine real function on a box attains a
minimum at a corner.  Thus FC.5, the inequality form of the complete
Bernoulli family, FC.8, and minimization over the switching cube really are
equivalent.  This verifies the note's warning that the full inhomogeneous
box is not a strict reduction.

The tangent signs are correct: FC.11 says that
`psi(u)-K log cosh(u)` has a one-sided minimum at the right endpoint `t`, so
`psi'(t)-K tanh(t)<=0`; evenness and convexity give `psi'(t)>=0`.

## 2. Fixed-size/Krawtchouk and radial information

In the expansion of `P_A`, the spin average requires an Eulerian edge set,
and the independent `tau` average separately requires even cardinality.
Flipping a uniformly random `k`-set multiplies an `ell`-edge monomial by the
normalized Krawtchouk multiplier in FC.17.  This proves FC.18.  Counting the
numbers `(K+U)/2` and `(K-U)/2` of positive and negative `Y_e` proves the
coefficient formula FC.19.

There is a small but important normalization point in the claimed
equivalence.  Inverting the Krawtchouk transform of the **values** in FC.18
first returns `W_ell rho^ell/P_A(rho 1)`.  Its `ell=0` entry equals
`1/P_A(rho 1)` because `W_0=1`; hence the scale is recovered, and then every
`W_ell` is recovered since `rho>0`.  So the equivalence is valid.  Merely
knowing the lower bounds `>=1` would not be equivalent, and the draft does
not claim that.

FC.20 then gives `Z_A(u)` for all `u`.  Linear independence of the finitely
many functions `cosh(h u)` for distinct nonnegative integer `h` makes this
equivalent to the absolute-energy histogram.  The `L^2` identity FC.21 is
the elementary relation `cosh^2 z=(cosh 2z+1)/2`.

## 3. Exact order-eight witness

The external file `computations/results/m8_minimizer_orbits.json` records an
exhaustive enumeration of all `2^21` root-gauged signings, exactly two
signed-permutation/global-sign cap-10 classes, and the same displayed energy
histogram for both.  The verifier script reproduces the common histogram,
caps, overlap traces, and one-vertex response tables exactly.

I also independently reconstructed the two ground-state covariance matrices
over rational numbers.  Their characteristic polynomials are

```math
2^{-8}(2z-5)^2(2z-1)^6
\quad\hbox{and}\quad
2^{-8}(2z-3)^4(2z-1)^4,
```

which verifies the spectra and traces `14` and `10` in FC.24.

The `t>=3` minimization argument is valid.  Any nonminimal signing has cap at
least `12`, hence at least one of its 128 projective spin states contributes
`cosh(12t)/128`; either displayed class is bounded above by `cosh(10t)`.
Moreover

```math
{\cosh(12t)\over\cosh(10t)}\ge {e^{2t}\over2}>128
\qquad(t\ge3).
```

Thus both classes are actual pressure minimizers in the stated range.  The
field response is indeed the one-new-spin response: global spin reversal
shows that maximizing `|H_A(x)+b.x|` already includes both choices of the
new spin.

## 4. Scope of the ceiling

The witness proves that the absolute-energy histogram does not determine
overlap geometry or even a one-vertex extension response **within the class
of actual finite pressure minimizers**.  This kills any universal theorem
whose only input is FC.10--FC.21.

It does not prove the same failure along the target regime
`t=beta/sqrt(N)->0`, nor does it exclude extra rigidity of asymptotic
contracted-temperature minimizers.  Therefore “the radial route is closed”
should be read as “radial data alone cannot give a universal actual-minimizer
channel theorem,” not as a scalable small-`t` counterexample.  A theorem
using small-`t` minimizer rigidity in addition to the radial histogram
remains logically open.

## 5. Mean-field source mappings

The Augeri mapping has the correct scale.  The bridge-gradient set lies in
`lambda t` times the convex hull of rank-one sign matrices, so its Rademacher
width is at most

```math
C lambda t sqrt(mn)(sqrt m+sqrt n).
```

The diagonal-Hessian correction is at most `lambda t^2 mn`.  Both are
`O(N)` at comparable splits.  Since bit-products form a subset of
row-products, the direction
`I_row^leftarrow<=I_bit^leftarrow=O(N)` is correct and yields no `o(N)`
conclusion.

The Lacker--Mukherjee--Yeung exclusion is also sound.  Their hypotheses are
continuous strong log-concavity hypotheses.  Subtracting
`c||b||_2^2` from a smooth bounded-Hessian extension makes the continuous
extension arbitrarily strongly concave while changing every Boolean vertex
by the same constant.  The discrete Gibbs law and its reverse product gap
are unchanged, so a canonical extension/rounding theorem is indispensable.

## 6. AC.24 and the narrowed SML

FC.8 is optimizer-specific but equivalent to the complete switching
landscape.  FC.10--FC.21 retain only radial data, and the witness shows those
data cannot identify AC.24.  Thus there is no hidden proof of product-shadow
or directed-dependence control in this note.

One wording refinement matters: controlling **one** summand in AC.24 by
`o(N)` does not by itself close the pressure gain, because the other summand
may remain linear.  It is a valid intermediate narrowing only if paired with
an independent bound on the other term.  Moreover, as documented in the
companion adversarial audit, the complete AC.17 best-response oracle
reconstructs `L` on point-mass inputs.  A genuinely strict next lemma must
extract a stated nonradial statistic smaller than that oracle and use it to
control both terms, or combine control of one term with a separately proved
actual-child theorem for the other.
