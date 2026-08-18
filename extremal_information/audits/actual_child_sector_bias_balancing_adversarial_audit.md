# Adversarial audit: sector-bias balancing for actual children

## Verdict

**SB.1--SB.3 PASS, with an essential target-orientation qualification.**
The sector-mixture comparison, the universal factor `1/2`, and the
optimizer-specific `D_infinity` normalization are exact.  For any pair of
exact pressure-minimizing children there exists a choice of row direction
and relative orientation for which the canonical inverse row factor has a
dimension-free max-divergence from the fair row law.

This is an existential statement about one admissible parent presentation.
It does not control both orientations simultaneously and does not prove
that a different variational or target-reaching argument selects the same
orientation.  It also does not bound the joint interaction, row total
correlation, or reverse product projection.

The proposed SB.4 density calculation was algebraically correct but
redundant and weaker than IC.10/AC.18; it should not be counted as new
optimizer-specific progress.

## 1. Exact sector-mixture comparison

With the notation of SB.1, the neutral and field-shifted sector weights are

```math
w_a={e^{a\gamma_D}\over2\cosh\gamma_D},
\qquad
q_a={e^{a(\gamma_D+g)}\over2\cosh(\gamma_D+g)}.
```

Their pointwise ratio is

```math
{q_a\over w_a}
=e^{ag}{\cosh\gamma_D\over\cosh(\gamma_D+g)}.
```

Since every sector output likelihood is positive,

```math
\sum_aq_az_D^a
\ge\min_a{q_a\over w_a}\sum_aw_az_D^a,
```

and the minimum is exactly

```math
e^{-|g|}{\cosh\gamma_D\over\cosh(\gamma_D+g)}.
```

Thus (SB.2)--(SB.3) have the correct normalization and sign.  The penalty
is nonnegative because the `w`-average of `q_a/w_a` equals one.  No
unmentioned symmetry of the sector likelihoods is used.

## 2. Orientation and filtration choice

Suppose the children are labelled so that
`|gamma_C|<=|gamma_D|`, and use `D` as the base child whose spin word is
seen by each erased row.  Choose `epsilon` so that
`epsilon gamma_C` has sign opposite to `gamma_D`.  Writing

```math
a=|\gamma_C|,
\qquad d=|\gamma_D|,
\qquad0\le a\le d,
```

the comparison coefficient becomes

```math
\kappa=e^{-a}{\cosh d\over\cosh(d-a)}.
```

The addition formula yields

```math
\cosh d\ge\cosh(d-a)\cosh a,
```

and hence

```math
\kappa\ge e^{-a}\cosh a
={1+e^{-2a}\over2}\ge{1\over2}.
```

Therefore (SB.6) is exact and uniform in both bias magnitudes.  The constant
is sharp for this positive-mixture comparison as `a=d -> infinity`.

The two choices used here are legitimate for an upper construction:

- permuting the two vertex blocks and sending `B` to `B^T` is a bijection
  of bridge signings and permits either child to be the erased-row base;
- replacing the second child by its negative gives the other relative
  orientation, and the augmented child pressure is unchanged under
  `D -> -D`.

The first operation is a presentation change.  The second selects a
different admissible parent signing; it is not an identity between the two
orientation-conditioned bridge laws.  Consequently SB.2 proves that a
balanced presentation is **available**, not that both orientations are
balanced or that an independently optimized bridge objective necessarily
attains its best value in this presentation.  Any downstream recurrence
must analyze the same selected orientation rather than silently switch
back after using SB.2.

## 3. Optimizer envelope and inverse-escort normalization

Whichever child is selected as the base is still an exact pressure
minimizer, because both input children are assumed actual minimizers.  EE.2
therefore gives

```math
z_D^0(b;t,t)\ge e^{-\delta_n(t)}.
```

Combining this with SB.2 gives

```math
z_{C\to D}^\epsilon(b;t,t)
\ge{1\over2}e^{-\delta_n(t)}.                       \tag{A.1}
```

Every sector mixture likelihood has fair-row mean one.  If
`z=z_(C->D)^epsilon` and

```math
{dr\over dU_n}={z^{-\lambda}\over E_Uz^{-\lambda}},
```

then convexity gives `E_U z^(-lambda)>=1`, while (A.1) gives
`z^(-lambda)<=exp(lambda(delta_n+log2))`.  Hence

```math
D_\infty(r\Vert U_n)
\le\lambda\{\delta_n(t)+\log2\}.
```

The finite-Renyi and min-entropy consequences follow from this density
bound.  Finally,

```math
0\le\delta_n(t)\le n\log\cosh t\le{nt^2\over2},
```

so at `t=beta/sqrt(N)` the right side is at most

```math
\lambda\left({\beta^2\over2}+\log2\right).
```

All constants in SB.3 are therefore correct.  This is a genuine
optimizer-specific strengthening for the **canonical row factor**: its
max-divergence, not only its collision divergence, is dimension-free and
independent of both child biases.

## 4. SB.4 redundancy

The proposed SB.4 argument used

```math
f={dP\over dr},\qquad g={dr\over dU},
```

and the correct inequality

```math
E_U(fg)^2=E_r(f^2g)
\le\|g\|_\infty E_rf^2.
```

It therefore validly combined a relative `D_2(P||r)` bound with SB.3 to
obtain an absolute `D_2(P||U)` bound.  It did not, however, add a new
conclusion:

- IC.10 already proves directly, for every orientation, hybrid parameter,
  row prefix or all-other-row conditioning,

  ```math
  D_2(q_s(R_i\mid\mathcal C)\Vert U_n)
  \le\lambda^2u^2n;
  ```

- AC.18 already gives the same absolute bound for every factor of a
  globally optimal row-product shadow.

The SB.4 bound adds a positive optimizer-envelope term and a larger
coefficient, so it is weaker.  Removing it, or retaining it only as a
consistency check, is the correct evidentiary treatment.

## 5. Exact research scope

The audited implication is

```text
two exact minimizing children
  + freedom to choose block order and relative orientation
  -> one canonical iid-row presentation with O(1) component D_infinity.
```

It does **not** imply:

- max-divergence control for an orientation fixed in advance;
- simultaneous component control in the two orientation sectors;
- a tight latent representation of the full bridge escort;
- `o(N)` row total correlation or marginal retuning;
- `o(N)` canonical error `J` or reverse product gap;
- a target-reaching pressure recurrence.

Thus sector bias is no longer an obstruction to the existence of a diffuse
canonical row presentation.  The surviving obstruction is still the joint
resource in that same balanced presentation.  This is a real Level-5
sharpening, not a Level-5-to-6 closure.
