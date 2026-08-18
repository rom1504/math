# Independent audit of the projective-atom and flip-information theorem

Status: **passed**.  This note independently checks
[`actual_child_projective_atom_flip_information_ceiling.md`](actual_child_projective_atom_flip_information_ceiling.md)
and its condensation as Theorem 37.60.  It checks only the claims made there;
it does not promote the abstract flip-sign-cone witness to an actual-child
construction.

## 1. Projective normalization and the tangent bound

With uniform expectations in (PA.1),

```math
Z_A(t)=2^{-m}\sum_x\cosh(tH_A(x)).                    \tag{VPA.1}
```

The full augmented atom has mass
`exp(t tau H_A(x))/(2^(m+1)Z_A)`.  Since `H_A(-x)=H_A(x)`, summing the two
representatives of `[x]` gives exactly

```math
\bar\nu_A(\tau,[x])={e^{t\tau H_A(x)}\over2^mZ_A(t)}, \tag{VPA.2}
```

so PA.4 has neither a missing factor two nor an extra sector factor.

For one edge, the exact flip ratio is

```math
R_A(\{e\})
=\cosh(2t)-\sinh(2t)\,\mathbb E_{\bar\nu_A}Y_e\ge1.
```

Therefore

```math
\mathbb E Y_e\le {\cosh(2t)-1\over\sinh(2t)}=\tanh t. \tag{VPA.3}
```

If `tau H_A(x)>0`, at least one edge at that atom has `Y_e=1`.  If its
mass is `p`, the elementary lower bound `E Y_e>=p-(1-p)=2p-1`, combined
with (VPA.3), gives `p<=(1+tanh t)/2`.  If `tau H_A(x)<=0`, (VPA.2) and
`Z_A>=1` give `p<=2^(-m)`.  Thus PA.7 and Theorem 37.60's atom ceiling are
correct with the stated normalization.

## 2. Conditional rank-one-word normalization

Let

```math
D_\epsilon=Z_A^+Z_D^\epsilon+Z_A^-Z_D^{-\epsilon}.
```

After conditioning on `epsilon=tau_1 tau_2` and writing `s=tau_1`, the
probability of a latent triple is exactly

```math
\nu_\epsilon(s,x,y)
={e^{st(H_A(x)+\epsilon H_D(y))}\over2^{m+n}D_\epsilon}. \tag{VPA.4}
```

For a fixed signed rank-one word `Q`, each value of `s` has two preimages,
related by `(x,y)->(-x,-y)`.  If `(x,y)` represents `Q` in the positive
sector, the negative-sector representatives have the same two quadratic
energies and the opposite exponent.  Hence

```math
\Pr(Q)
={2e^{tE}+2e^{-tE}\over2^{m+n}D_\epsilon}
={2^{2-m-n}\over D_\epsilon}\cosh(tE),
\quad E=H_A(x)+\epsilon H_D(y).                       \tag{VPA.5}
```

There are `2^(m+n-1)` signed rank-one words.  Moreover

```math
D_\epsilon
=2^{2-m-n}\sum_Q\cosh(tE_Q),                         \tag{VPA.6}
```

so (VPA.5) sums to one.  This verifies PA.11, including the coefficient
`2^(2-m-n)`, and verifies PA.12 after maximizing the even function `cosh`.
The symmetry `Pr(Q)=Pr(-Q)` is also exact.  Consequently a projective word
`[Q]={Q,-Q}` has mass `2 eta_epsilon`; this is consistent with the audit's
statement that `eta_epsilon=o(1)` is equivalent to exclusion of a fixed-mass
antipodal two-word block.

## 3. The `m`-edge inversion

Fix `x_1=1`.  The star variables give

```math
x_i={Y_{1i}\over\tau a_{1i}}\qquad(2\le i\le m),       \tag{VPA.7}
```

while the triangle product gives

```math
\tau={Y_{12}Y_{13}Y_{23}\over a_{12}a_{13}a_{23}}.    \tag{VPA.8}
```

These equations recover a unique `(tau,[x])` from any sign pattern on the
`m` edges in PA.14.  Substitution into the predicted value of `Y_23`
returns the prescribed `Y_23`, so this proves surjectivity as well as
injectivity.  The map in PA.15 is therefore genuinely a bijection.

For `c=cosh(2t)` and `s=sinh(2t)`, expansion of
`exp(-2tY_e)=c-sY_e` gives PA.18.  Ordered by inclusion, the transform from
Walsh moments to flip values is subset triangular with diagonal
`(-s)^|S|`, nonzero for `t>0`.  Walsh inversion then recovers the whole law.
Thus PA.16--PA.19 and Theorem 37.60's flip-information claim are correct.

## 4. Scope of the half-atom construction

For the antipodal pair in PA.20, the contribution to every flip moment is

```math
\cosh\!\left(2t\sum_{e\in S}Y_e(z_0)\right)\ge1.
```

For the uniform component, averaging over `tau` gives the average of
`cosh(2t sum_(e in S) a_e x_i x_j)`, also at least one.  When `S` is
nonempty, the quadratic Walsh polynomial cannot vanish identically because
its distinct degree-two characters have nonzero coefficients; hence the
uniform contribution is strictly larger than one.  PA.21--PA.22 follow.

Crucially, the source audit and Theorem 37.60 both say explicitly that this
law need not have the quadratic Gibbs form and is not asserted to come from
an actual minimizing signing.  The logical conclusion is therefore only
that the **directions** of all flip inequalities cannot force diverging
min-entropy.  That abstract-versus-actual disclaimer is adequate and is not
silently lost in the theorem condensation.

## 5. Verdict

All requested checks pass:

- PA.7 uses the correctly normalized projective atom and the correct
  `tanh(t)` one-flip tangent;
- PA.11 is a normalized conditional law on signed rank-one words;
- the `m`-edge map is bijective and its transform is invertible for `t>0`;
- the half-atom witness is scoped only to the abstract inequality cone.

The source display PA.12 and the theorem condensation also use the correct
hyperbolic-cosine normalization.  No mathematical correction is needed.
